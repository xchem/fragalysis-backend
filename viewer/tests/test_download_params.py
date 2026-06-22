"""Tests for download flag granularity.

Historically a single serializer flag, ``all_aligned_structures``, bundled nine
distinct aligned-structure file types, and several file groups (YAML, extra
files, PyMOL scripts, README) were always included with no flag at all.

``get_download_params`` now maps an individual serializer flag to each file type
while keeping ``all_aligned_structures`` as a backwards-compatible umbrella, and
exposes the previously-unconditional groups as toggleable flags. These tests pin
that behaviour down without building an actual archive.
"""

from viewer.download_structures import get_download_params
from viewer.serializers import DownloadStructuresSerializer

# The nine file types previously bundled under all_aligned_structures.
ALIGNED_PROTEIN_PARAMS = (
    'apo_file',
    'bound_file',
    'apo_solv_file',
    'apo_desolv_file',
    'ligand_pdb',
    'ligand_sdf',
    'ligand_smiles',
)
ALIGNED_OTHER_PARAMS = ('sdf_info', 'smiles_info')

# Groups that used to be added unconditionally; default True.
DEFAULT_GROUP_PARAMS = ('yaml_files', 'extra_files', 'pymol_scripts', 'readme')


def _params(**overrides):
    """Run the data through the serializer (applying defaults) and split it."""
    data = {'target_name': 'TestTarget', 'target_access_string': 'lb00000-1'}
    data.update(overrides)
    serializer = DownloadStructuresSerializer(data=data)
    assert serializer.is_valid(), serializer.errors
    return get_download_params(serializer.validated_data)


def test_umbrella_enables_all_aligned_structures():
    """Legacy contract: all_aligned_structures=True turns on all nine types."""
    protein_params, other_params, _ = _params(all_aligned_structures=True)
    for param in ALIGNED_PROTEIN_PARAMS:
        assert protein_params[param] is True, param
    for param in ALIGNED_OTHER_PARAMS:
        assert other_params[param] is True, param


def test_defaults_disable_aligned_structures():
    """With nothing set, no aligned-structure type is requested."""
    protein_params, other_params, _ = _params()
    for param in ALIGNED_PROTEIN_PARAMS:
        assert protein_params[param] is False, param
    for param in ALIGNED_OTHER_PARAMS:
        assert other_params[param] is False, param


def test_single_aligned_file_type_in_isolation():
    """A single individual flag yields only that file type (umbrella off)."""
    protein_params, other_params, _ = _params(ligand_sdf=True)
    assert protein_params['ligand_sdf'] is True
    for param in ALIGNED_PROTEIN_PARAMS:
        if param != 'ligand_sdf':
            assert protein_params[param] is False, param
    for param in ALIGNED_OTHER_PARAMS:
        assert other_params[param] is False, param


def test_individual_flag_ored_with_umbrella():
    """Setting both the umbrella and an individual flag still enables all."""
    protein_params, _, _ = _params(all_aligned_structures=True, apo_file=True)
    for param in ALIGNED_PROTEIN_PARAMS:
        assert protein_params[param] is True, param


def test_default_groups_included_by_default():
    """YAML/extra/scripts/README are included unless explicitly excluded."""
    _, other_params, _ = _params()
    for param in DEFAULT_GROUP_PARAMS:
        assert other_params[param] is True, param


def test_default_group_can_be_excluded():
    """Each default group can be individually switched off."""
    for param in DEFAULT_GROUP_PARAMS:
        _, other_params, _ = _params(**{param: False})
        assert other_params[param] is False, param
        # the other groups stay on
        for other in DEFAULT_GROUP_PARAMS:
            if other != param:
                assert other_params[other] is True, other
