#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 13:19:51 2020
@author: Warren
Script to check sdf file format for Fragalysis upload
"""
import logging

import validators

logger = logging.getLogger(__name__)

# Set .sdf format version here
version = 'ver_1.2'


def check_compound_set(description_mol, validate_dict, update=None):
    del update

    # Must have a 'generation_date'
    if not description_mol.HasProp('generation_date'):
        validate_dict = add_warning(
            molecule_name='File error',
            field='compound set',
            warning_string="Molecule has no generation_date",
            validate_dict=validate_dict,
        )
        return validate_dict
    # That's of the form "<Y>-<M>-<D>"...
    g_date = description_mol.GetProp('generation_date')
    y_m_d = g_date.split('-')
    if len(y_m_d) != 3:
        validate_dict = add_warning(
            molecule_name='File error',
            field='compound set',
            warning_string="Molecule has no generation_date is not Y-M-D (g_date)",
            validate_dict=validate_dict,
        )
        return validate_dict

    return validate_dict


def add_warning(molecule_name, field, warning_string, validate_dict):
    validate_dict['molecule_name'].append(molecule_name)
    validate_dict['field'].append(field)
    validate_dict['warning_string'].append(warning_string)

    return validate_dict


def check_sdf(sdf_file, validate_dict):
    """
    Checks if .sdf file can be read and follows naming format:
    'compound-set_<name>.sdf' with <name> replaced with
    the name you wish to give it. e.g. compound-set_fragmenstein.sdf

    :sdf_file: is the sdf in the specified format
    :return: Updates validate dictionary with pass/fail message
    """
    # Check filename
    if sdf_file.startswith("compound-set_") and sdf_file.endswith(".sdf") is False:
        validate_dict = add_warning(
            molecule_name='File error',
            field='_File_name',
            warning_string=f"Illegal filename: {str(sdf_file)} found",
            validate_dict=validate_dict,
        )

    return validate_dict


def check_ver_name(blank_mol, check_version, validate_dict):
    """
    Checks if blank mol:
    The name (title line) of this molecule should be the
    file format specification version e.g. ver_1.0 (as defined in this document)

    :blank_mol: rdkit mol of blank mol from an SD file
    :return: Updates validate dictionary with pass/fail message
    """

    ver_name = blank_mol.GetProp('_Name')
    if ver_name != check_version:
        validate_dict = add_warning(
            molecule_name=blank_mol.GetProp('_Name'),
            field='_Name',
            warning_string=f'Illegal version: {ver_name} found. Should be {check_version}',
            validate_dict=validate_dict,
        )

    return validate_dict


def check_blank_mol_props(mol, validate_dict):
    # check for compulsory fields in blank mols
    fields = [
        'ref_url',
        'submitter_name',
        'submitter_email',
        'submitter_institution',
        'generation_date',
        'method',
    ]
    for field in fields:
        validate_dict = missing_field_check(mol, field, validate_dict)

    return validate_dict


def check_blank_prop(blank_mol, validate_dict):
    """
    Checks if blank mol properties have a description

    :blank_mol: rdkit mol of blank mol from an SD file
    :return: Updates validate dictionary with pass/fail message
    """

    # Check if properties populated
    property_dict = blank_mol.GetPropsAsDict()

    # Properties to ignore
    prop_ignore_list = ['ref_mols', 'ref_pdb']

    for key, value in property_dict.items():
        if value == '' and key not in prop_ignore_list:
            validate_dict = add_warning(
                molecule_name=blank_mol.GetProp('_Name'),
                field=key,
                warning_string=f'Description for {key} missing',
                validate_dict=validate_dict,
            )
        if key == 'ref_url' and check_url(value) is False:
            validate_dict = add_warning(
                molecule_name=blank_mol.GetProp('_Name'),
                field=key,
                warning_string=f'Illegal URL {value} provided',
                validate_dict=validate_dict,
            )

    return validate_dict


def check_url(value):
    """
    Checks if url provided exists. No internet connection required.
    Checks URL using Validators package

    :value: value associated with 'ref_url' key
    :return: False if URL can not be validated
    """

    valid = validators.url(value)
    if valid is not True:
        return False


def check_name_characters(name, validate_dict):
    legal_non_alnum = ['-', '_', '.']
    for char in name:
        if not char.isalnum() and char not in legal_non_alnum:
            validate_dict = add_warning(
                molecule_name=name,
                field='_Name',
                warning_string=f'Illegal character {char} found',
                validate_dict=validate_dict,
            )

    return validate_dict


def missing_field_check(mol, field, validate_dict):
    props_dict = mol.GetPropsAsDict()
    if field not in list(props_dict.keys()):
        validate_dict = add_warning(
            molecule_name=mol.GetProp('_Name'),
            field=field,
            warning_string=f'Field {field} not found!',
            validate_dict=validate_dict,
        )

    return validate_dict


def check_mol_props(mol, validate_dict):
    # Check for (mandatory, isolated) missing fields
    fields = ['ref_mols']
    for field in fields:
        validate_dict = missing_field_check(mol, field, validate_dict)
    # More complex checks?
    # One of ref_pdb and lhs_pdb must be set
    if not (mol.HasProp('ref_pdb') or mol.HasProp('lhs_pdb')):
        validate_dict = add_warning(
            molecule_name=mol.GetProp('_Name'),
            field='ref_pdb/lhs_pdb',
            warning_string="Molecule has neither 'ref_pdb' nor 'lhs_pdb' property",
            validate_dict=validate_dict,
        )

    return validate_dict
