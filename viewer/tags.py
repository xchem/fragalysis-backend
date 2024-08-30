import logging
import re
from typing import Any

import pandas as pd
from django.core.exceptions import MultipleObjectsReturned
from django.db import IntegrityError, transaction
from django.db.models import CharField, Exists, F, OuterRef, Subquery, Value
from django.db.models.functions import Concat

from scoring.models import SiteObservationGroup

from .models import (
    Pose,
    SiteObservation,
    SiteObservationTag,
    SiteObvsSiteObservationTag,
    TagCategory,
    Target,
)
from .utils import sanitize_boolean_column

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


class CuratedTagSubquery(Exists):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, tag):
        query = SiteObvsSiteObservationTag.objects.filter(
            site_observation=OuterRef('pk'),
            site_obvs_tag=tag,
        )
        super().__init__(query)


def get_tag_cols(columns):
    return [
        k
        for k in columns
        if any(k.startswith(f'[{cat}]') for cat in CURATED_TAG_CATEGORIES)
    ]


# def get_upload_cols(columns):
#     return  [k for k in columns if any(k == f'{cat} alias' for cat in TAG_CATEGORIES)]

# def get_alias_cols(columns):
#     return [k for k in columns if any(k == f'{cat} upload name' for cat in TAG_CATEGORIES)]


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
        logger.debug('tag_groups :%s', tag_groups.keys)
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
                msg = f"Inconsistent aliases for tag {tag}: " + ':'.join(
                    [
                        f"{key[1]} in rows {row_idx}"
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

    # add auto-generated names before...
    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        upload_tag = f'upload_tag_{category.category.lower()}'
        values.append(upload_tag)
        header.append(f'{category.category} upload name')
        annotations[upload_tag] = UploadTagSubquery(category.category)

    # ... the aliases
    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        tag = f'tag_{category.category.lower()}'
        values.append(tag)
        header.append(f'{category.category} alias')
        annotations[tag] = TagSubquery(category.category)

    values.append(next(iter(POSE_COL.values())))
    header.append(next(iter(POSE_COL.keys())))

    pattern = re.compile(r'\W+')  # non-alphanumeric characters
    for tag in SiteObservationTag.objects.filter(
        category__in=TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES),
        target=target,
    ):
        # for reasons unknown, mypy thinks tag is a string
        tagname = f'tag_{pattern.sub("_", tag.tag).strip().lower()}'  # type: ignore[attr-defined]
        values.append(tagname)
        header.append(f'[{tag.category}] {tag.tag}')  # type: ignore[attr-defined]
        annotations[tagname] = CuratedTagSubquery(tag)

    # fmt: off

    return header, annotations, values


def load_tags_from_csv(filename: str, target: Target) -> list[str]:  # type: ignore [return]
    # from viewer.tags import load_tags_from_csv; from viewer.models import Target; target = Target.objects.get(title='A71EV2A'); load_tags_from_csv('metadata_modified_ok.xlsx', target)

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

    for column in tag_cols:
        df[column] = df[column].fillna(False)
        df[column] = sanitize_boolean_column(df[column])

    # header, _, _ = get_metadata_fields(target)

    # upload_cols = get_upload_cols(df.columns)
    # alias_cols = get_alias_cols(df.columns)

    tag_aliases, alias_errors = validate_aliases(df)
    errors.extend(alias_errors)

    with transaction.atomic():
        try:
            tags = SiteObservationTag.objects.filter(
                category__category__in=TAG_CATEGORIES,
                target=target,
            )

            qs = SiteObservation.objects.filter(longcode__in=df['Long code'])

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
                pose.save()

                # attach rest of the observations
                for so in so_qs:
                    so.pose = pose
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
                ]
                if so.pose_display_name != pose_name:
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

            SiteObservation.objects.bulk_update(so_update, ['pose'])

            cats = TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES)
            curated_tags = SiteObservationTag.objects.filter(
                category__category__in=cats,
                target=target,
            )

            curated_db = set(curated_tags.values_list('tag', flat=True))
            curated_df = set(tagnames.keys())
            curated_delete = curated_db.difference(curated_df)
            # curated_create = curated_df.difference(curated_db)
            # curated_existing = curated_df.intersection(curated_db)

            # delete those missing from the uploaded file
            curated_tags.filter(tag__in=curated_delete).delete()

            # create or update new tags from file
            for tc in tag_cols:
                so_group = SiteObservationGroup(target=target)
                so_group.save()

                category_name = tagnames[tc]
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
                    so_from_db = set(
                        so_tag.site_observations.values_list('pk', flat=True)
                    )
                except SiteObservationTag.DoesNotExist:
                    so_tag = SiteObservationTag(
                        tag=tc,
                        tag_prefix='',
                        upload_name=tc,
                        category=cat,
                        target=target,
                        mol_group=so_group,
                    )
                    so_tag.save()

                site_observations = qs.filter(
                    longcode__in=df.loc[df[tc] == True]['Long code']
                )
                # compare observations from file and db, update only if different
                so_from_df = set(site_observations.values_list('pk', flat=True))
                if so_from_db != so_from_df:
                    so_group.site_observation.add(*site_observations)
                    so_tag.site_observations.add(*site_observations)

            if errors:
                # TODO: pass on collected error messages
                raise IntegrityError('smth')

        except IntegrityError:
            # TODO: need to give user feedback what went wrong but
            # don't know which mechanim is going to be used
            return errors
