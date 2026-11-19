#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 13:19:51 2020
@author: Warren
Helper for collecting computed-set upload warnings
"""


def add_warning(molecule_name, field, warning_string, validate_dict):
    validate_dict['molecule_name'].append(molecule_name)
    validate_dict['field'].append(field)
    validate_dict['warning_string'].append(warning_string)

    return validate_dict
