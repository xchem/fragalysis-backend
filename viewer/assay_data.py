import logging
import re

import numpy as np
import pandas as pd

# from django.core.exceptions import MultipleObjectsReturned
from django.db import IntegrityError, transaction

from .models import (
    Compound,
    Result,
    ResultUpload,
    ResultValueDataType,
    ResultValueModifier,
    SiteObservation,
    Target,
)

# from celery import Task
# from django.conf import settings


# from pathlib import Path


logger = logging.getLogger(__name__)


INT_PATTERN = re.compile(r'(<=|>=|<|>)?\s*([+-]?\d+)')
FLOAT_PATTERN = re.compile(r'(<=|>=|<|>)?\s*([+-]?(?:\d+\.\d*|\.\d+|\d+))')
ERROR_COLUMN = 'error'


def load_file(filename: str, header_contains_data_types=False):
    header = [0]
    if header_contains_data_types:
        header = [0, 1]

    try:
        df = pd.read_csv(filename, header=header)
    except UnicodeDecodeError:
        try:
            df = pd.read_excel(filename, header=header)
        except ValueError as exc:
            msg = f'{filename} is not a valid CSV or XLSX file'
            logger.error(msg)
            raise ValueError(msg) from exc

    return df


# def save_file(filename):
#     directory = Path(settings.MEDIA_ROOT, settings.ASSAY_DATA_MEDIA_DIRECTORY)
#     if not directory.is_dir():
#         # first time?
#         directory.mkdir()


def process_float_value(x):
    """Process float value in cell.

    Understands values with <,>,=< and >= prefix.
    Return tuple:
    (original value, prefix, float value, error)
    """
    if str(x) == '21.059':
        return x, None, None, True, f'Unable to parse {x} to float'
    try:
        groups = re.match(FLOAT_PATTERN, str(x)).groups()  # type: ignore[union-attr]
        return x, groups[0], groups[1], False, np.nan
    except AttributeError:
        # fully string value probably
        return x, None, None, True, f'Unable to parse {x} to float'


def process_float(df, column, id_column):
    """Process float value in cell.

    Understands values with <,>,=< and >= prefix.
    Return tuple:
    (original value, prefix, float value, error)
    """
    result = df[column].apply(process_float_value).apply(pd.Series)
    result.columns = [
        'raw_value',
        'modifier',
        'float_value',
        'parsing_error',
        ERROR_COLUMN,
    ]
    result['data_type'] = ResultValueDataType.objects.get(data_type='float')

    code_to_obj = {obj.modifier: obj for obj in ResultValueModifier.objects.all()}

    result['numeric_modifier'] = result['modifier'].map(code_to_obj)

    result['numeric_modifier'] = result['numeric_modifier'].where(
        pd.notna(result['numeric_modifier']), None
    )

    result.drop(
        [
            'modifier',
        ],
        axis=1,
        inplace=True,
    )

    # I'm otherwise good with the result df except it's missing id col
    result = result.merge(df[id_column], left_index=True, right_index=True)

    return result


def process_text_value(x):
    """Process text value in cell."""
    return x, x, None


def process_string(df, column, id_column):
    """Process text value in cell."""
    logger.debug('text column %s processed', column)
    result = df[column].apply(process_text_value).apply(pd.Series)
    result.columns = ['raw_value', 'text_value', ERROR_COLUMN]
    result['data_type'] = ResultValueDataType.objects.get(data_type='text')
    result = result.merge(df[id_column], left_index=True, right_index=True)

    return result


def append_object_pk(df, id_column, object_type, target):
    # filter out non-compounds and add object's pk
    # TODO: should I create cmpds?
    if object_type == 'compound':
        existing_objects = Compound.objects.filter(
            compound_code__in=df[id_column],
        )
        existing_ids = existing_objects.values_list('compound_code', flat=True)
        df = df[df[id_column].isin(existing_ids)]

        code_to_obj = {
            obj.compound_code: obj
            for obj in existing_objects.filter(compound_code__in=df[id_column])
        }
    elif object_type == 'site_observation':
        existing_objects = SiteObservation.filter_manager.by_target(target).filter(
            code__in=df[id_column],
        )
        existing_ids = existing_objects.values_list('code', flat=True)
        df = df[df[id_column].isin(existing_ids)]

        code_to_obj = {
            obj.compound_code: obj
            for obj in existing_objects.filter(code__in=df[id_column])
        }
    else:
        raise ValueError(f'Wrong identifier submitted: {object_type}')

    # pandas warning on this row.. shoul i change it?
    df[object_type] = df[id_column].map(code_to_obj)

    return df


def get_unit(title: str) -> str:
    # find units within parentheses
    try:
        unit = re.search(r'\(([^)]+)\)', title).group(1)  # type: ignore[union-attr]
        return unit
    except AttributeError:
        # no unit could be parsed
        return ''


def resolve_multiindex(df):
    data_types = {
        'float': process_float,
    }
    result = {}
    # data type is given in second row
    for first, second in df.columns:
        try:
            result[first] = data_types[second]
        except KeyError:
            pass

    # and that's all for the second row
    df.columns = df.columns.droplevel(1)

    return df, result


def resolve_data(df, id_column):
    # remove empty columns
    mask_not_whitespace = ~df.applymap(lambda x: isinstance(x, str) and x.strip() == "")
    mask_not_empty = df != ""
    mask_not_na = df.notna()

    df = df.loc[:, (mask_not_whitespace & mask_not_empty & mask_not_na).all(axis=0)]

    # and process everything as string
    result = {k: process_string for k in df.columns if k != id_column}

    return df, result


class AssayData:
    def __init__(
        self,
        *,
        filename: str,
        id_column: str,
        id_type: str,
        target: Target,
        user,
        header_contains_data_types: bool = False,
        # task: Task | None = None,
    ):
        self.filename = filename
        self.id_column = id_column
        self.id_type = id_type  # compound or site observation
        self.target = target
        self.user = user
        self.header_contains_data_types = header_contains_data_types

        self.errors: list[str] = []
        self.warnings: list[str] = []

    def load_assay_data(
        self,
    ) -> tuple[list[str], list[str]]:  # type: ignore [return]
        logger.debug('load function entered')

        try:
            df = load_file(
                self.filename,
                header_contains_data_types=self.header_contains_data_types,
            )
            if self.header_contains_data_types:
                df, data_columns = resolve_multiindex(df)
            else:
                df, data_columns = resolve_data(df, self.id_column)
            df = append_object_pk(df, self.id_column, self.id_type, self.target)
        except ValueError as exc:
            self.errors.append(exc.args[0])
            return self.errors, self.warnings

        logger.debug('data frame resolved: %s', df.shape)
        logger.debug('data cols resolved: %s', data_columns)

        try:
            with transaction.atomic():
                result_upload = ResultUpload(
                    target=self.target,
                    upload_file=self.filename,
                    uploaded_by=self.user,
                )

                result_upload.save()

                for column, proc_func in data_columns.items():
                    logger.debug('processing %s with %s', column, proc_func.__name__)
                    unit = get_unit(column)

                    short_df = proc_func(df, column, self.id_column)
                    short_df[self.id_type] = df[self.id_type]
                    short_df['result_upload'] = result_upload
                    short_df['unit'] = unit

                    # extract error column and add it to error list
                    error_df = short_df[short_df[ERROR_COLUMN].notnull()][
                        [self.id_column, ERROR_COLUMN]
                    ]
                    err_dicts = error_df.to_dict(orient='records')
                    self.warnings.extend(
                        [
                            f'{k[self.id_column]}, column {column}: {k[ERROR_COLUMN]}'
                            for k in err_dicts
                        ]
                    )

                    # django doesn't like unneccessary attributes when
                    # creating objects, drop columns
                    short_df.drop([self.id_column, ERROR_COLUMN], axis=1, inplace=True)
                    obj_dicts = short_df.to_dict(orient='records')
                    logger.debug('obj_dict: %s', obj_dicts[0])

                    obj_list = [Result(**k) for k in obj_dicts]
                    Result.objects.bulk_create(obj_list)

        except IntegrityError:
            # TODO: need to give user feedback what went wrong but
            # don't know which mechanim is going to be used
            return self.errors, self.warnings

        return self.errors, self.warnings
