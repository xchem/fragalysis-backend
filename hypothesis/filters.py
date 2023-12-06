from django_filters import rest_framework as filters


class Vector3DFilter(filters.FilterSet):
    site_observation = filters.CharFilter(label="SiteObservation ID")
    cmpd_id = filters.NumberFilter(label="Compound ID")
    number = filters.NumberFilter(label="Number")
    smiles = filters.CharFilter(label="Smiles")
    vector_type = filters.CharFilter(label="Type")
