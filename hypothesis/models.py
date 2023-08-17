from django.db import models

from hypothesis.definitions import IntTypes, VectTypes
from viewer.models import SiteObservation, Target, Compound


class TargetResidue(models.Model):
    """Model to store residue information - to curate the probes"""
    # The target it relates to
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    # The residue name
    res_name = models.CharField(max_length=20)
    # The residue number
    res_num = models.IntegerField()
    # The chain
    chain_id = models.CharField(max_length=4)

    class Meta:
        unique_together = ("target_id", "res_num", "res_name", "chain_id")


class InteractionPoint(models.Model):
    # The molecule id
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    # Set the molecule and protein identifier
    protein_atom_name = models.CharField(max_length=255)
    molecule_atom_name = models.CharField(max_length=255)
    targ_res = models.ForeignKey(TargetResidue, on_delete=models.CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "site_observation",
                    "protein_atom_name",
                    "molecule_atom_name",
                ],
                name="unique_siteobvs_protatom_molatom",
            ),
        ]        
        

class Interaction(models.Model):
    """Model to store the interaction information."""
    int_ver_choices, default_int_ver, int_type_choices, default_int_type = (
        IntTypes().define_int_types()
    )
    interaction_version = models.CharField(
        choices=int_ver_choices, max_length=2, default=default_int_ver
    )
    interaction_type = models.CharField(
        choices=int_type_choices, max_length=2, default=default_int_type
    )
    interaction_point = models.ForeignKey(InteractionPoint, on_delete=models.CASCADE)
    distance = models.FloatField(null=True)
    score = models.FloatField(null=True)
    prot_smarts = models.TextField(null=True)
    mol_smarts = models.TextField(null=True)

    class Meta:
        unique_together = (
            "interaction_type",
            "interaction_version",
            "interaction_point",
        )


class Vector(models.Model):
    # The compound it relates to
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE)
    # The smiles of the vector
    smiles = models.CharField(max_length=255)
    # Vector type
    type = models.CharField(choices=VectTypes().get_vect_types(), max_length=2)

    class Meta:
        unique_together = ("cmpd_id", "smiles", "type")


class Vector3D(models.Model):
    # The molecule it relates to
    # mol_id = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    # The vector it relates to
    vector = models.ForeignKey(Vector, on_delete=models.CASCADE)
    # The number on this
    number = models.IntegerField()
    # The start position
    start_x = models.FloatField(null=True)
    start_y = models.FloatField(null=True)
    start_z = models.FloatField(null=True)
    # The end position
    end_x = models.FloatField(null=True)
    end_y = models.FloatField(null=True)
    end_z = models.FloatField(null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["site_observation", "vector", "number"],
                name="unique_siteobvs_vector_number",
            ),
        ]
