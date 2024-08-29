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


def validate_aliases(df) -> tuple[list[tuple[str, str]], list[str]]:
    """Check upload name and alias combinations.

    Each upload name should have only one alias and vice versa.
    """
    result = []
    errors = []
    for cat in TAG_CATEGORIES:
        col_upload = f'{cat} upload name'
        col_alias = f'{cat} alias'

        tag_groups = df.loc[:, [col_upload, col_alias]].groupby(col_upload)
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

    for column in tag_cols:
        df[column] = sanitize_boolean_column(df[column])

    # header, _, _ = get_metadata_fields(target)

    # upload_cols = get_upload_cols(df.columns)
    # alias_cols = get_alias_cols(df.columns)

    tag_aliases, alias_errors = validate_aliases(df.columns)
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

            # poses got really complicated now. instead of just
            # checking them I need to allow create, delete and
            # reassign operations. but how?

            # if I loop over pose name and encounter one I don't have
            # then need to create one. The observation that lists it
            # will be assigned to it. also means, if this was a single
            # observation in some other pose, this needs to be deleted

            # hold on, is it really that complicated? if the
            # spreadsheet represents the new truth, then I should be
            # able to create new set of poses from it. if all cells
            # are filled (each observation belongs into a pose) and
            # there's no extra poses (each pose must have at least one
            # observation) then surely it's job done, no?

            # how to check pose equality? I think all I have is name,
            # observations belonging to pose and main observation. but
            # name is user-changeable, so.. what then? if the name is
            # changed, I consider a new pose?

            # other than that - check all poses for equality:
            # - name
            # - observations belonging to it
            # - main observation (name matches main obvs short code)

            # non-matching options:
            # - name not found: completely new pose
            # - different number of observations in pose: add or remove
            # - different main observation: change observation (if found by name)

            # procedure:

            pose_update = []
            for so in qs:
                pose_name = df.loc[df['Long code'] == so.code, POSE_COL]
                try:
                    pose = poses.get(display_name=pose_name)
                except Pose.DoesNotExist:
                    # TODO: notify user
                    # errors = True
                    continue
                except MultipleObjectsReturned:
                    # It's not actually ideal, using display_name, there's no guarantee it's unique
                    # TODO: notify user
                    # errors = True
                    continue

                # concrete pose found, continue
                if so.pose.pk != pose.pk:
                    so.pose = pose
                    pose_update.append(so)
                    # TODO: attach to update object

            SiteObservation.objects.bulk_update(pose_update, ['pose'])
            # that's nice, but can it be done at all? poses need
            # main_obvs, so any moving around would destroy this

            # the loaded file is treated as a new ground truth, all changes
            # should be propagated to database including deleting the curated
            # tags. therefore, start by deleting all previously created
            # curated tags and recreate them from file
            # update: no that won't work. deleting and recreating would kill history
            cats = TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES)
            existing_curated = SiteObservationTag.objects.filter(
                category__category__in=cats,
                target=target,
            )

            # hold on, need the correct category
            for tc in tag_cols:
                so_group = SiteObservationGroup(target=target)
                so_group.save()

                # tag names (tc) are coming in format like '[Other] P2
                # Wave 1' where [Other] is category. extract the name
                # and category
                splits = tc.split(']')
                if len(splits) > 0:
                    category_name = splits[0].strip('[ ')
                    # if for some reason user decided to use brackets in tagname
                    tagname = ']'.join(splits[1:])
                else:
                    # category not given in tagname, raise error, notify user
                    continue

                try:
                    cat = cats.get(category=category_name)
                except TagCategory.DoesNotExist:
                    # TODO: notify user
                    continue

                try:
                    so_tag = existing_curated.get(tag=tagname)
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
