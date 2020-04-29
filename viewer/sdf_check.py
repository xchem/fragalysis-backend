#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 13:19:51 2020
@author: Warren
Script to check sdf file format for Fragalysis upload
"""

from rdkit import Chem
import validators
import numpy as np
import os
from viewer.models import Protein, CompoundSet
import datetime

# Set .sdf format version here
version = 'ver_1.2'

def check_compound_set(description_mol, validate_dict):
    y_m_d = description_mol.GetProp('generation_date').split('-')

    submitter_dict = {'submitter__name': description_mol.GetProp('submitter_name'),
                      'submitter__email': description_mol.GetProp('submitter_email'),
                      'submitter__institution': description_mol.GetProp('submitter_institution'),
                      'submitter__generation_date': datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2])),
                      'submitter__method': description_mol.GetProp('method')}

    query = CompoundSet.objects.filter(**submitter_dict)

    if len(query)!=0:
        validate_dict = add_warning(molecule_name='File error',
                                    field='compound set',
                                    warning_string="a compound set with the auto_generated name " + query[0].submitter.unique_name + " already exists (change method name in blank mol method field)",
                                    validate_dict=validate_dict)

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
        validate_dict = add_warning(molecule_name='File error',
                                    field='_File_name',
                                    warning_string="illegal filename: " + str(sdf_file) + " found",
                                    validate_dict=validate_dict)

    return validate_dict


def check_pdb(mol, validate_dict, target=None, zfile=None):
    """
    Checks if .pdb file can be read

    :mol: rdkit mol read from SD file
    :return: Updates validate dictionary with pass/fail message
    """

    # Check if pdb filepath given and exists
    test_fp = mol.GetProp('ref_pdb')

    # {'zip_obj': zf, 'zf_list': zip_names}

    if zfile:
        pdb_option = mol.GetProp('ref_pdb')
        # name = pdb_option.split('/')[-1]
        if zfile:
            if pdb_option not in zfile['zf_list']:
                validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                            field='ref_pdb',
                                            warning_string="path " + str(pdb_option) + " can't be found in uploaded zip file (list: " + str(zfile['zf_list'] +")"),
                                            validate_dict=validate_dict)


    # if not '_0' in test_fp or not '.pdb' in test_fp:
    #     validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
    #                                 field='ref_pdb',
    #                                 warning_string="illegal pdb assingment for " + str(test_fp),
    #                                 validate_dict=validate_dict)

    # if test_fp.endswith(".pdb"):
    #     if os.path.exists(test_fp) is False:
    #         validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
    #                                     field='ref_pdb',
    #                                     warning_string="path " + str(test_fp) + " does not exist",
    #                                     validate_dict=validate_dict)

    # else:
    if target and not test_fp.endswith(".pdb"):
        query = Protein.objects.filter(code=str(target + '-' + test_fp))
        if len(query)==0:
            validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                        field='ref_pdb',
                                        warning_string="pdb for " + str(test_fp) + " does not exist in fragalysis",
                                        validate_dict=validate_dict)

    return validate_dict


def check_SMILES(mol, validate_dict):
    """
    Checks if SMILES can be read by rdkit

    :mol: rdkit mol read from SD file
    :return: Updates validate dictionary with pass/fail message
    """
    # Check SMILES
    smi_check = mol.GetProp('original SMILES')

    m = Chem.MolFromSmiles(smi_check, sanitize=False)
    if m is None:
        validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                    field='original SMILES',
                                    warning_string="invalid SMILES %s" % (smi_check,),
                                    validate_dict=validate_dict)

    return validate_dict


def check_ver_name(blank_mol, version, validate_dict):
    """
    Checks if blank mol:
    The name (title line) of this molecule should be the
    file format specification version e.g. ver_1.0 (as defined in this document)

    :blank_mol: rdkit mol of blank mol from an SD file
    :return: Updates validate dictionary with pass/fail message
    """

    ver_name = blank_mol.GetProp('_Name')
    if ver_name != version:
        validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                    field='_Name',
                                    warning_string='illegal version: %s found. Should be %s' % (ver_name, version),
                                    validate_dict=validate_dict)

    return validate_dict


def check_blank_mol_props(mol, validate_dict):
    # check for compulsory fields in blank mols
    fields = ['ref_url', 'submitter_name', 'submitter_email', 'submitter_institution', 'generation_date', 'method']
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

    for key, value in zip(property_dict.keys(), property_dict.values()):
        if value == '' and key not in prop_ignore_list:
            validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='description for %s missing' % (key,),
                                        validate_dict=validate_dict)
        if key == 'ref_url' and check_url(value) == False:
            validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='illegal URL %s provided' % (value,),
                                        validate_dict=validate_dict)

    return validate_dict


def check_field_populated(mol, validate_dict):
    """
    Checks if all compulsory fields are populated:
        1. ref_mols - a comma separated list of the fragments
        2. ref_pdb - either (a) a filepath (relative to the sdf file)
            to an uploaded pdb file
        3. original SMILES - the original smiles of the compound
            before any computation was carried out

    :mol: rdkit mol other than blank_mol
    :return: Updates validate dictionary with pass/fail message
    """

    # Compuslory fields
    compulsory_fields = ['ref_pdb', 'ref_mols', 'original SMILES']

    property_dict = mol.GetPropsAsDict()
    for key, value in zip(property_dict.keys(), property_dict.values()):
        if value == '' and key in compulsory_fields:
            validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='value for %s missing' % (key,),
                                        validate_dict=validate_dict)

    return validate_dict


def check_url(value):
    """
    Checks if url provided exists. No internet connection required.
    Checks URL using Validators package

    :value: value associated with 'ref_url' key
    :return: False if URL can not be validated
    """

    valid = validators.url(value)
    if valid != True:
        return False


def check_name_characters(name, validate_dict):
    legal_non_alnum = ['-', '_', '.']
    for char in name:
        if not char.isalnum() and char not in legal_non_alnum:
            validate_dict = add_warning(molecule_name=name,
                                        field='_Name',
                                        warning_string='illegal character %s found' % (char,),
                                        validate_dict=validate_dict)

    return validate_dict


def missing_field_check(mol, field, validate_dict):
    props_dict = mol.GetPropsAsDict()
    if not field in props_dict.keys():
        validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                    field=field,
                                    warning_string='%s field not found!' % (field,),
                                    validate_dict=validate_dict)

    return validate_dict


def check_mol_props(mol, validate_dict):
    # check for missing fields
    fields = ['ref_pdb', 'ref_mols', 'original SMILES']
    for field in fields:
        validate_dict = missing_field_check(mol, field, validate_dict)

    return validate_dict


def validate(sdf_file, target=None, zfile=None):
    validated = True
    validate_dict = {'molecule_name': [],
                     'field': [],
                     'warning_string': []}

    # Check sdf filename & can be read
    validate_dict = check_sdf(sdf_file, validate_dict)

    suppl = Chem.SDMolSupplier(sdf_file)
    print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]
    validate_dict = check_compound_set(blank_mol, validate_dict)
    other_mols = []
    for i in range(1, len(suppl)):
        other_mols.append(suppl[i])

    # all mol checks
    # - all mols have the same properties
    all_props = []
    for mol in suppl:
        all_props.extend([key for key in mol.GetPropsAsDict().keys()])
    unique_props = list(set(all_props))
    for mol in suppl:
        props = [key for key in mol.GetPropsAsDict().keys()]
        diff_list = np.setdiff1d(props, unique_props)
        for diff in diff_list:
            add_warning(molecule_name=mol.GetProp('_Name'),
                        field='property (missing)',
                        warning_string='%s property is missing from this molecule' % (diff,),
                        validate_dict=validate_dict)

    # Check version in blank mol
    validate_dict = check_ver_name(blank_mol, version, validate_dict)

    # Check compuslory fields in blank mol props
    validate_dict = check_blank_mol_props(blank_mol, validate_dict)

    # Check properties have been described and validate url
    validate_dict = check_blank_prop(blank_mol, validate_dict)

    # main mols checks
    # - missing compulsary fields
    # - check name characters
    # - check pdb assignment and if pdb filepath exists
    # - check compulsory field populated
    # - check SMILES can be opended by rdkit
    # (check api for pdb if fragalysis)
    for m in other_mols:
        validate_dict = check_mol_props(m, validate_dict)
        validate_dict = check_name_characters(m.GetProp('_Name'), validate_dict)
        validate_dict = check_pdb(m, validate_dict, target, zfile)
        validate_dict = check_field_populated(m, validate_dict)
        validate_dict = check_SMILES(m, validate_dict)

    if len(validate_dict['molecule_name']) != 0:
        validated = False

    return validate_dict, validated