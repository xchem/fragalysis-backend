import logging
import re

import numpy as np
import pandas as pd
from django.db import IntegrityError, transaction

from .models import (
    Compound,
    Result,
    ResultProperty,
    ResultUpload,
    ResultValueDataType,
    ResultValueModifier,
    SiteObservation,
    Target,
)

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


def process_float_value(x):
    """Process float value in cell.

    Understands values with <,>,=< and >= prefix.
    Return tuple:
    (original value, prefix, float value, error)
    """
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

    return result, ResultValueDataType.objects.get(data_type='float')


def process_text_value(x):
    """Process text value in cell."""
    return x, x, None


def process_text(df, column, id_column):
    """Process text value in cell."""
    logger.debug('text column %s processed', column)
    result = df[column].apply(process_text_value).apply(pd.Series)
    result.columns = ['raw_value', 'text_value', ERROR_COLUMN]
    result = result.merge(df[id_column], left_index=True, right_index=True)

    return result, ResultValueDataType.objects.get(data_type='text')


def process_int_value(x):
    """Process float value in cell.

    Understands values with <,>,=< and >= prefix.
    Return tuple:
    (original value, prefix, float value, error)
    """
    try:
        groups = re.match(INT_PATTERN, str(x)).groups()  # type: ignore[union-attr]
        return x, groups[0], groups[1], False, np.nan
    except AttributeError:
        # fully string value probably
        return x, None, None, True, f'Unable to parse {x} to int'


def process_integer(df, column, id_column):
    """Process integer value in cell.

    Understands values with <,>,=< and >= prefix.
    Return tuple:
    (original value, prefix, int value, error)
    """
    result = df[column].apply(process_int_value).apply(pd.Series)
    result.columns = [
        'raw_value',
        'modifier',
        'int_value',
        'parsing_error',
        ERROR_COLUMN,
    ]

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

    return result, ResultValueDataType.objects.get(data_type='integer')


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
            obj.code: obj for obj in existing_objects.filter(code__in=df[id_column])
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
    result = {k: process_text for k in df.columns if k != id_column}

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

                order = 0
                for column, proc_func in data_columns.items():
                    logger.debug('processing %s with %s', column, proc_func.__name__)
                    order = order + 1
                    unit = get_unit(column)

                    short_df, data_type = proc_func(df, column, self.id_column)

                    result_property, _ = ResultProperty.objects.get_or_create(
                        result_property=column,
                        unit=unit,
                        target=self.target,
                        order=order,
                        data_type=data_type,
                    )

                    short_df[self.id_type] = df[self.id_type]
                    short_df['result_upload'] = result_upload
                    short_df['result_property'] = result_property

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


def convert(upload, property_id, new_type):
    errors = []
    warnings: list[str] = []
    result_property = ResultProperty.objects.get(pk=property_id)
    data_type = ResultValueDataType.objects.get(data_type=new_type)

    qs = Result.objects.filter(
        result_property=result_property,
        result_upload=upload,
    )

    try:
        qs, old_cols = _clear_old_value(qs, result_property.data_type)
    except TypeError as exc:
        logger.error(exc.args[0])
        errors.append(exc.args[0])
        return errors, warnings

    if new_type == 'float':
        df_proc_func = process_float
        qs_proc_func = _to_float
    elif new_type == 'text':
        df_proc_func = process_text
        qs_proc_func = _to_text
    elif new_type == 'integer':
        df_proc_func = process_integer
        qs_proc_func = _to_integer
    else:
        errors.append(f'Unknown data type: {new_type}')
        return errors, warnings

    # TODO: maybe select only necessary columns, memory
    df = pd.DataFrame.from_records(qs.values())

    proc_df, _ = df_proc_func(df, 'raw_value', 'id')

    # extract error column and add it to error list
    error_df = proc_df[proc_df[ERROR_COLUMN].notnull()][['id', ERROR_COLUMN]]
    # id is needed as value in error_df
    proc_df = proc_df.set_index('id')
    obj_dicts = proc_df.to_dict(orient='dict')
    logger.debug('obj_dicts: %s', obj_dicts.keys())

    err_dicts = error_df.to_dict(orient='records')

    warnings.extend([f'{k["id"]}: {k[ERROR_COLUMN]}' for k in err_dicts])

    qs, cols = qs_proc_func(qs, obj_dicts)

    try:
        with transaction.atomic():
            Result.objects.bulk_update(
                qs,
                cols + old_cols,
            )
            result_property.data_type = data_type
            result_property.save()
    except IntegrityError as exc:
        logger.error(exc.args[0])
        errors.append(exc.args[0])
        return errors, warnings

    return errors, warnings


def _clear_old_value(qs, old_type):
    if old_type.data_type == 'float':
        cols = ['float_value', 'numeric_modifier']
    elif old_type.data_type == 'text':
        cols = ['text_value']
    elif old_type.data_type == 'integer':
        cols = ['int_value']
    else:
        cols = []

    for obj in qs:
        for field in cols:
            setattr(obj, field, None)

    return qs, cols


def _to_float(qs, obj_dicts):
    for obj in qs:
        obj.float_value = obj_dicts['float_value'][obj.id]
        obj.numeric_modifier = obj_dicts['numeric_modifier'][obj.id]
        obj.parsing_error = obj_dicts['parsing_error'][obj.id]

    return qs, ['float_value', 'numeric_modifier', 'parsing_error']


def _to_text(qs, obj_dicts):
    for obj in qs:
        obj.text_value = obj_dicts['text_value'][obj.id]

    return qs, [
        'text_value',
    ]


def _to_integer(qs, obj_dicts):
    for obj in qs:
        obj.int_value = obj_dicts['int_value'][obj.id]
        obj.numeric_modifier = obj_dicts['numeric_modifier'][obj.id]
        obj.parsing_error = obj_dicts['parsing_error'][obj.id]

    return qs, ['float_value', 'numeric_modifier', 'parsing_error']
