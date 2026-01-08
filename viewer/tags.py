import hashlib
import logging
from typing import Any, Generator

import numpy as np
import pandas as pd
from django.contrib.auth.models import User
from django.core.exceptions import MultipleObjectsReturned
from django.db import IntegrityError, transaction
from django.db.models import CharField, Count, Exists, F, OuterRef, Q, Subquery, Value
from django.db.models.functions import Concat

from scoring.models import SiteObservationGroup

from .models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    CompoundIdentifier,
    CompoundIdentifierType,
    Pose,
    QualityStatusType,
    QuatAssembly,
    SiteObservation,
    SiteObservationQualityStatus,
    SiteObservationTag,
    SiteObvsSiteObservationTag,
    TagCategory,
    Target,
    Xtalform,
    XtalformQuatAssembly,
    XtalformSite,
)
from .utils import alphanumerator, clean_object_id

logger = logging.getLogger(__name__)

TAG_CATEGORIES = (
    'ConformerSites',
    'CanonSites',
    'CrystalformSites',
    'Quatassemblies',
    'Crystalforms',
)
CURATED_TAG_CATEGORIES = ('Series', 'Forum', 'Other')

META_HEADER = {
    'Code': 'code',
    'Long code': 'longcode',
    'Experiment code': 'experiment__code',
    'Compound code': 'cmpd__compound_code',
    'Smiles': 'smiles',
    'Centroid res': 'canon_site_conf__canon_site__centroid_res',
    'Downloaded': 'downloaded',
}

POSE_COL = {'Pose': 'pose__display_name'}


class TagSubquery(Subquery):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, category):
        # fmt: off
        query = SiteObservationTag.objects.filter(
            pk=Subquery(
                SiteObvsSiteObservationTag.objects.filter(
                    site_observation=OuterRef(OuterRef('pk')),
                    site_obvs_tag__category=TagCategory.objects.get(
                        category=category,
                    ),
                ).values('site_obvs_tag')[:1]
            )
        ).annotate(
            combitag=Concat(
                F('tag_prefix'),
                Value(' - '),
                F('tag'),
                output_field=CharField(),
            ),
        ).values('combitag')[0:1]
        super().__init__(query)
        # fmt: on


class UploadTagSubquery(Subquery):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, category):
        # fmt: off
        query = SiteObservationTag.objects.filter(
            pk=Subquery(
                SiteObvsSiteObservationTag.objects.filter(
                    site_observation=OuterRef(OuterRef('pk')),
                    site_obvs_tag__category=TagCategory.objects.get(
                        category=category,
                    ),
                ).values('site_obvs_tag')[:1]
            )
        ).values('upload_name')[0:1]
        super().__init__(query)
        # fmt: on


class ShortTagSubquery(Subquery):
    """Annotate SiteObservation with short tag of given category"""

    def __init__(self, category):
        # fmt: off
        query = SiteObservationTag.objects.filter(
            pk=Subquery(
                SiteObvsSiteObservationTag.objects.filter(
                    site_observation=OuterRef(OuterRef('pk')),
                    site_obvs_tag__category=TagCategory.objects.get(
                        category=category,
                    ),
                ).values('site_obvs_tag')[:1]
            )
        ).values('short_tag')[0:1]
        super().__init__(query)
        # fmt: on


class CuratedTagSubquery(Exists):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, tag):
        query = SiteObvsSiteObservationTag.objects.filter(
            site_observation=OuterRef('pk'),
            site_obvs_tag=tag,
        )
        super().__init__(query)


class CustomIdentifierSubquery(Subquery):
    """Annotate SiteObservation with short tag of given category"""

    def __init__(self, identifier_type):
        # fmt: off
        query = CompoundIdentifier.objects.filter(
            type__name=identifier_type,
            compound=OuterRef('cmpd'),
        ).values('name')[0:1]
        super().__init__(query)
        # fmt: on


def get_tag_cols(columns):
    return [
        k
        for k in columns
        if any(k.startswith(f'[{cat}]') for cat in CURATED_TAG_CATEGORIES)
    ]


def strip_alias(alias: str) -> str:
    """Strip prefix from alias if present"""
    splits = alias.split('-')
    if len(splits) == 1:
        return alias
    else:
        return splits[1].strip()


def strip_catname(tag: str) -> tuple[str | None, str]:
    # tag names (tc) are coming in format like '[Other] P2 Wave 1'
    # where [Other] is category. extract the name and category
    splits = tag.split(']')
    category_name = None
    if len(splits) > 0:
        category_name = splits[0].strip('[ ')
        # if for some reason user decided to use brackets in tagname
        tagname = ']'.join(splits[1:]).strip()

    return category_name, tagname


def validate_aliases(df) -> tuple[list[tuple[str, str]], list[str]]:
    """Check upload name and alias combinations.

    Each upload name should have only one alias and vice versa.
    """
    result = []
    errors = []
    for cat in TAG_CATEGORIES:
        col_upload = f'{cat} upload name'
        col_alias = f'{cat} alias'

        tag_groups = df.loc[:, [col_upload, col_alias]].groupby(col_upload).groups
        for tag in tag_groups.keys():
            name_groups = (
                df.loc[tag_groups[tag], [col_upload, col_alias]]
                .groupby([col_upload, col_alias])
                .groups
            )
            if len(name_groups) == 1:
                # key in this name_groups is tuple(upload_name, alias)
                result.append(next(iter(name_groups.keys())))
            else:
                # this means multiple groups, multiple aliases for the same tag
                msg = f"Inconsistent aliases for tag '{tag}': " + ':'.join(
                    [
                        f"'{key[1]}' in rows {list(row_idx)}"
                        for key, row_idx in name_groups.items()
                    ]
                )

                logger.warning(msg)
                errors.append(msg)
    return result, errors


def get_metadata_fields(target: Target) -> tuple[list[str], dict[str, Any], list[str]]:
    """Compile metadata.csv header and annotation objects."""

    annotations = {}

    header: list[str] = list(META_HEADER.keys())
    values: list[str] = list(META_HEADER.values())

    # add auto-generated names first...
    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        upload_tag = get_ann_tag(f'{category.category} upload name')
        values.append(upload_tag)
        header.append(f'{category.category} upload name')
        annotations[upload_tag] = UploadTagSubquery(category.category)

    # ... then the short tags, ...
    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        short_tag = get_ann_tag(f'{category.category} short tag')
        values.append(short_tag)
        header.append(f'{category.category} short tag')
        annotations[short_tag] = ShortTagSubquery(category.category)

    # ... then aliases, ...
    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        tag = get_ann_tag(f'{category.category} alias')
        values.append(tag)
        header.append(f'{category.category} alias')
        annotations[tag] = TagSubquery(category.category)

    values.append(next(iter(POSE_COL.values())))
    header.append(next(iter(POSE_COL.keys())))

    for tag in SiteObservationTag.objects.filter(
        category__in=TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES),
        target=target,
    ):
        tagname = get_ann_tag(tag.tag)  # type: ignore[attr-defined]
        values.append(tagname)
        header.append(f'[{tag.category}] {tag.tag}')  # type: ignore[attr-defined]
        annotations[tagname] = CuratedTagSubquery(tag)

    # fmt: off

    # and finally custom identifiers
    custom_identifiers = CompoundIdentifierType.objects.filter(
        name__in=CompoundIdentifier.objects.filter(
            # NB! see comment about filter_manager in managers.py for
            # compound only fetching LHS upload compounds. not
            # convinced it's the desired behaviour here
            compound__in=Compound.filter_manager.by_target(target=target),
        ).values('type'),
    ).values_list('name', flat=True)
    header.extend(custom_identifiers)
    for identifier in custom_identifiers:
        ann_name = get_ann_tag(identifier)
        values.append(ann_name)
        annotations[ann_name] = CustomIdentifierSubquery(identifier)

    # and finally-finally quality states
    header.extend([
        'Main status',
        'GOOD count',
        'MEDIOCRE count',
        'BAD count',
    ])
    values.extend([
        'main_status',
        'count_good',
        'count_mediocre',
        'count_bad',
    ])
    annotations['main_status'] = Subquery(
        SiteObservationQualityStatus.objects.filter(
            site_observation=OuterRef('pk'),
            main_status=True,
        ).values('status')
    )
    annotations['count_good'] = Count(
        'siteobservationqualitystatus',
        filter=Q(siteobservationqualitystatus__status__status='GOOD'),
    )
    annotations['count_mediocre'] = Count(
        'siteobservationqualitystatus',
        filter=Q(siteobservationqualitystatus__status__status='MEDIOCRE'),
    )
    annotations['count_bad'] = Count(
        'siteobservationqualitystatus',
        filter=Q(siteobservationqualitystatus__status__status='BAD'),
    )

    # and finally-finally-finally refinementresolution from soakdb
    header.append('RefinementResolution')
    values.append('refinementresolution')
    annotations['refinementresolution'] = F('experiment__refinement_resolution')
    # annotations['refinementresolution'] = Subquery(
    #     SiteObservation.objects.filter(
    #         pk=OuterRef('pk'),
    #     ).values('experiment__refinement_resolution')
    # )

    return header, annotations, values


def load_tags_from_file(filename: str, target: Target, user: User | None = None) -> list[str]:  # type: ignore [return]
    # from viewer.tags import load_tags_from_file; from viewer.models import Target; target = Target.objects.get(pk=1); load_tags_from_file('metadata.csv', target)

    errors: list[str] = []

    try:
        df = pd.read_csv(filename)
    except UnicodeDecodeError:
        try:
            df = pd.read_excel(filename)
        except ValueError:
            msg = f'{filename} is not a valid CSV or XLSX file'
            errors.append(msg)
            logger.error(msg)
            return errors

    tag_cols = get_tag_cols(df.columns)
    tagnames = {strip_catname(k)[1]: strip_catname(k)[0] for k in tag_cols}

    logger.debug('tag_cols :%s', tag_cols)

    # a quick sanity check for file cols
    header, _, _ = get_metadata_fields(target)
    unknowns = [k for k in df.columns if k not in header and k not in tag_cols]
    if unknowns:
        for c in unknowns:
            msg = f"Invalid column name: '{c}'"
            logger.warning(msg)
            errors.append(msg)

    for column in tag_cols:
        df[column] = df[column].fillna(False)
        df[column], errors = sanitize_boolean_column(df[column], errors)

    tag_aliases, alias_errors = validate_aliases(df)
    errors.extend(alias_errors)

    try:
        with transaction.atomic():
            tags = SiteObservationTag.objects.filter(
                category__category__in=TAG_CATEGORIES,
                target=target,
            )

            # fmt: off
            qs = SiteObservation.filter_manager.by_target(
                target=target,
            ).filter(
                longcode__in=df['Long code'],
            )
            # fmt: on

            alias_update = []
            for upload_name, alias in tag_aliases:
                try:
                    tag = tags.get(upload_name=upload_name)
                except SiteObservationTag.DoesNotExist:
                    msg = f"Unknown tag '{upload_name}'"
                    logger.error(msg)
                    errors.append(msg)
                    continue

                tagname = strip_alias(alias)

                if tag.tag != tagname:
                    tag.tag = tagname
                    alias_update.append(tag)

            logger.debug('alias updates %s', alias_update)
            SiteObservationTag.objects.bulk_update(alias_update, ['tag'])

            poses = Pose.objects.filter(
                main_site_observation__experiment__experiment_upload__target=target,
            )

            # the spreadsheet is considered to be the new truth, make
            # the db match. procedure:
            # - delete poses that should be deleted
            # - create a new ones
            # - attach orphaned observations to poses

            # NB! display name column is not guaranteed to be unique,
            # but I have nothing else to go on from the file
            pose_current = set(poses.values_list('display_name', flat=True))
            pose_df = set(df[next(iter(POSE_COL.keys()))])

            pose_delete = pose_current.difference(pose_df)
            pose_create = pose_df.difference(pose_current)

            poses.filter(display_name__in=pose_delete).delete()
            for pose_name in pose_create:
                logger.debug('fetching pose name %s', pose_name)
                so_df = df.loc[
                    df[next(iter(POSE_COL.keys()))] == pose_name, 'Long code'
                ]
                logger.debug('so_df %s', so_df)
                so_qs = qs.filter(longcode__in=so_df)
                logger.debug('pose so queryset: %s', so_qs)

                if not so_qs.exists():
                    msg = (
                        f'No observations found for pose {pose_name}.'
                        + ' Most likely upload data is incompatible with target data.'
                    )
                    logger.error(msg)
                    errors.append(msg)
                    continue

                try:
                    main_so = so_qs.get(code=pose_name)
                except SiteObservation.DoesNotExist:
                    # if the user is using custom pose name and it's
                    # not traceable back to observation
                    main_so = so_qs.first()
                    logger.debug('first main_so: %s', main_so)

                pose = Pose(
                    canon_site=main_so.canon_site_conf.canon_site,
                    compound=main_so.cmpd,
                    main_site_observation=main_so,
                    display_name=pose_name,
                )
                logger.debug('so pose instance: %s', pose)
                pose.save()

                # attach rest of the observations
                for so in so_qs:
                    so.pose = pose
                    logger.debug('so instance: %s', so)
                    so.save()

            # refresh the poses and...
            poses = Pose.objects.filter(
                main_site_observation__experiment__experiment_upload__target=target,
            )

            # ...check if there's any observations that do not belong to correct pose
            so_update = []
            for so in qs:
                pose_name = df.loc[
                    df['Long code'] == so.longcode, next(iter(POSE_COL.keys()))
                ].to_numpy()[0]
                if so.pose is None or so.pose.display_name != pose_name:
                    try:
                        pose = poses.get(display_name=pose_name)
                    except Pose.DoesNotExist as exc:
                        msg = f'This should not have happened: {exc}'
                        logger.error(msg)
                        errors.append(msg)
                        continue
                    except MultipleObjectsReturned as exc:
                        msg = f'This should not have happened: {exc}'
                        logger.error(msg)
                        errors.append(msg)
                        continue

                    # either of these error conditions should happen, the
                    # first is taken care by creating all the missing sets
                    # above, and the other with the set operation that
                    # discards the duplicates. leavnig it in just in case
                    # something goes horribly wrong

                    # concrete pose found, continue
                    so.pose = pose
                    so_update.append(so)

            logger.debug('so bulk update %s', so_update)
            SiteObservation.objects.bulk_update(so_update, ['pose'])

            cats = TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES)
            curated_tags = SiteObservationTag.objects.filter(
                category__in=cats,
                target=target,
            )

            curated_db = set(curated_tags.values_list('tag', flat=True))
            curated_df = set(tagnames.keys())
            curated_delete = curated_db.difference(curated_df)

            # delete those missing from the uploaded file
            curated_tags.filter(tag__in=curated_delete).delete()

            # create or update new tags from file
            for tc in tag_cols:
                category_name, tagname = strip_catname(tc)
                so_group = SiteObservationGroup(target=target)
                logger.debug('so group instance: %s', so_group)
                so_group.save()

                # category_name = tagnames[tc]
                if not category_name:
                    # category not given in tagname, raise error, notify user
                    msg = f'Category name not given for tag {tc}'
                    logger.error(msg)
                    errors.append(msg)
                    continue

                try:
                    cat = cats.get(category=category_name)
                except TagCategory.DoesNotExist:
                    msg = f"Unknown category name '{category_name}'"
                    logger.error(msg)
                    errors.append(msg)
                    continue

                try:
                    so_tag = curated_tags.get(tag=tagname)
                except SiteObservationTag.DoesNotExist:
                    so_tag = SiteObservationTag(
                        tag=tagname,
                        tag_prefix='',
                        upload_name=tagname,
                        category=cat,
                        target=target,
                        mol_group=so_group,
                        short_tag=tagname,
                    )
                    logger.debug('so tag instance: %s', so_tag)
                    so_tag.save()

                so_from_db = set(so_tag.site_observations.values_list('pk', flat=True))

                site_observations = qs.filter(
                    longcode__in=df.loc[df[tc] == True]['Long code']
                )

                # compare observations from file and db, update only if different
                so_from_df = set(site_observations.values_list('pk', flat=True))
                if so_from_db != so_from_df:
                    so_group.site_observation.add(*site_observations)
                    so_tag.site_observations.add(*site_observations)

            # tags done, read quality as well
            qual_column = 'Main status'
            qual_states = []
            for status_type in QualityStatusType.objects.exclude(status='NONE'):
                codes = df.loc[df[qual_column] == status_type.status][
                    'Long code'
                ].unique()
                for code in codes:
                    try:
                        so = qs.get(longcode=code)
                    except SiteObservation.DoesNotExist:
                        msg = (
                            f'SiteObservation {code} does not exist for {target.title}'
                        )
                        logger.error(msg)
                        errors.append(msg)
                        continue

                    # only add status if current main exists and is
                    # something else
                    add_status = False
                    try:
                        main_status = so.siteobservationqualitystatus_set.get(
                            main_status=True
                        )
                        if main_status.status == status_type:
                            add_status = True
                    except SiteObservationQualityStatus.DoesNotExist:
                        add_status = True

                    # there's a constraint in the model, so I don't
                    # think multiple objects needs to be handled here

                    if add_status:
                        qual_states.append(
                            SiteObservationQualityStatus(
                                site_observation=so,
                                status=status_type,
                                user=user,
                                main_status=True,
                                comment='Loaded from metadata.csv',
                            )
                        )

            SiteObservationQualityStatus.objects.bulk_create(qual_states)

            if errors:
                # log all errors
                logger.info('Errors found processing metadata file:')
                for line in errors:
                    logger.info('err: %s', line)
                raise IntegrityError('Errors encountered when processing metadata file')

    except IntegrityError:
        # TODO: need to give user feedback what went wrong but
        # don't know which mechanim is going to be used
        return errors


def sanitize_boolean_column(column, errors):
    """
    Sanitize a DataFrame column to boolean values.
    Handles various representations of booleans such as:
    - Python booleans
    - String representations ('True', 'False', 'true', 'false', 'yes', 'no')
    - Numeric representations (1, 0)
    """
    # don't understans why it considers int and str numbers to be the same
    true_values = {  # pylint: disable=duplicate-value
        'true',
        '1',
        't',
        'y',
        'yes',
        'True',
        'TRUE',
        True,
        1,
    }
    false_values = {  # pylint: disable=duplicate-value
        'false',
        '0',
        'f',
        'n',
        'no',
        'False',
        'FALSE',
        False,
        0,
        None,
        '',
        np.nan,
    }

    def convert_to_boolean(value, column, errors):
        if value in true_values:
            return True
        elif value in false_values:
            return False
        else:
            errors.append(f"Invalid boolean value '{value}' in column '{column}'")
            return value

    return column.apply(convert_to_boolean, args=(column.name, errors)), errors


def get_ann_tag(input_str: str) -> str:
    return hashlib.md5(input_str.encode()).hexdigest()


class TagManager:
    def __init__(self, target):
        self.target = target

    def add_tags_to_canon_sites(self, canon_site_pks: list[int]):
        cat = TagCategory.objects.get(category="CanonSites")
        qs = CanonSite.objects.filter(
            pk__in=canon_site_pks,
        ).annotate(
            code_prefix=F(
                'ref_conf_site__ref_site_observation__experiment__code_prefix'
            ),
        )

        for instance in qs:
            prefix = instance.canon_site_num
            # tag = canon_name_tag_map.get(val.versioned_key, "UNDEFINED")
            so_list = SiteObservation.objects.filter(
                canon_site_conf__canon_site=instance
            )

            # tag = val.versioned_key
            tag = f"{instance.name}+{instance.version}"
            try:
                short_tag = tag.split('-')[1][1:]
                short_tag = f"{instance.code_prefix}{short_tag}"

                # memo to self: there was an elaborate shcme here to
                # catch an error if the experiment code wasn't found
                # in the metadata file. if I remember correctly, this
                # was a manifestation of a different issue altogether
                # (corrupt data?) and included here only for reporting
                # purposes

            except IndexError:
                # non-standard tag
                short_tag = tag

            self.tag_observations(
                tag,
                prefix,
                category=cat,
                site_observations=so_list,
                short_tag=short_tag,
            )

        logger.debug("canon_site objects tagged")

    def add_tags_to_conformer_sites(self, canon_site_conf_pks: list[int]):
        numerators: dict[str, Generator[str, None, None]] = {}
        cat = TagCategory.objects.get(category="ConformerSites")
        qs = CanonSiteConf.objects.filter(
            pk__in=canon_site_conf_pks,
        ).annotate(
            code_prefix=F(
                'canon_site__ref_conf_site__ref_site_observation__experiment__code_prefix'
            ),
        )

        for instance in qs:
            if instance.canon_site.canon_site_num not in numerators.keys():
                numerators[instance.canon_site.canon_site_num] = alphanumerator()

            prefix = (
                f"{instance.canon_site.canon_site_num}"
                + f"{next(numerators[instance.canon_site.canon_site_num])}"
            )

            so_list = instance.siteobservation_set.all()
            # same comment as for canon sites, there was a lot of
            # error handling for when observations were'f found, but
            # I'm pretty sure that was broken data

            tag = instance.name
            try:
                short_tag = instance.name.split('-')[1][1:]
                short_tag = f"{instance.code_prefix}{short_tag}"
            except IndexError:
                short_tag = tag

            self.tag_observations(
                tag,
                prefix,
                category=cat,
                site_observations=so_list,
                hidden=True,
                short_tag=short_tag,
            )

        logger.debug("conf_site objects tagged")

    def add_tags_to_quatassemblies(self, quatassembly_pks: list[int]):
        cat = TagCategory.objects.get(category="Quatassemblies")
        qs = QuatAssembly.objects.filter(pk__in=quatassembly_pks)

        for instance in qs:
            prefix = f"A{instance.assembly_num}"
            tag = instance.name
            so_list = SiteObservation.objects.filter(
                xtalform_site__xtalform__in=XtalformQuatAssembly.objects.filter(
                    quat_assembly=instance
                ).values("xtalform")
            )
            self.tag_observations(
                tag,
                prefix,
                category=cat,
                site_observations=so_list,
            )

        logger.debug("quat_assembly objects tagged")

    def add_tags_to_xtalforms(self, xtalform_pks: list[int]):
        cat = TagCategory.objects.get(category="Crystalforms")
        qs = Xtalform.objects.filter(pk__in=xtalform_pks)

        for instance in qs:
            prefix = f"F{instance.xtalform_num}"
            so_list = SiteObservation.objects.filter(xtalform_site__xtalform=instance)
            tag = instance.name

            self.tag_observations(
                tag,
                prefix,
                category=cat,
                site_observations=so_list,
                clean_ids=False,
            )

        logger.debug("xtalform objects tagged")

    def add_tags_to_xtalformsites(self, xtalformsite_pks: list[int]):
        cat = TagCategory.objects.get(category="CrystalformSites")
        qs = XtalformSite.objects.filter(
            pk__in=xtalformsite_pks,
        ).annotate(
            code_prefix=F(
                'canon_site__ref_conf_site__ref_site_observation__experiment__code_prefix'
            ),
        )
        for instance in qs:
            prefix = (
                f"F{instance.xtalform.xtalform_num}" + f"{instance.xtalform_site_num}"
            )

            so_list = instance.siteobservation_set.all()
            tag = f"{instance.xtalform_site_id}/{instance.version}"
            try:
                # remove protein name and 'x'
                short_tag = instance.xtalform_site_id.split('-')[1][1:]
                short_tag = f"{instance.code_prefix}{short_tag}"
            except IndexError:
                short_tag = tag

            self.tag_observations(
                tag,
                prefix,
                category=cat,
                site_observations=so_list,
                hidden=True,
                short_tag=short_tag,
            )

        logger.debug("xtalform_sites objects tagged")

    def tag_observations(
        self,
        tag: str,
        prefix: str,
        category: TagCategory,
        site_observations: list,
        hidden: bool = False,
        short_tag: str | None = None,
        clean_ids: bool = True,
    ) -> None:
        try:
            # memo to self: description is set to tag, but there's
            # no fk to tag, instead, tag has a fk to
            # group. There's no uniqueness requirement on
            # description so there's no certainty that this will
            # be unique (or remain searchable at all because user
            # is allowed to change the tag name). this feels like
            # poor design but I don't understand the principles of
            # this system to know if that's indeed the case or if
            # it is in fact a truly elegant solution
            so_group = SiteObservationGroup.objects.get(
                target=self.target, description=tag
            )
        except SiteObservationGroup.DoesNotExist:
            so_group = SiteObservationGroup(target=self.target)
            so_group.save()
        except MultipleObjectsReturned:
            SiteObservationGroup.objects.filter(
                target=self.target, description=tag
            ).delete()
            so_group = SiteObservationGroup(target=self.target)
            so_group.save()

        name = f"{prefix} - {tag}" if prefix else tag
        tag = tag if short_tag is None else short_tag
        short_name = name if short_tag is None else f"{prefix} - {short_tag}"

        if clean_ids:
            tag = clean_object_id(tag)
            name = clean_object_id(name)
            short_name = clean_object_id(short_name)

        try:
            so_tag = SiteObservationTag.objects.get(
                upload_name=name, target=self.target
            )
            # Tag already exists
            # Apart from the new mol_group and molecules, we shouldn't be
            # changing anything.
            so_tag.mol_group = so_group
        except SiteObservationTag.DoesNotExist:
            so_tag = SiteObservationTag(
                tag=tag,
                tag_prefix=prefix,
                upload_name=name,
                category=category,
                target=self.target,
                mol_group=so_group,
                hidden=hidden,
                short_tag=short_name,
            )

        so_tag.save()

        so_group.site_observation.add(*site_observations)
        so_tag.site_observations.add(*site_observations)
