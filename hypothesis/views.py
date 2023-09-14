import logging
from itertools import chain

import json

from rest_framework import viewsets
from rest_framework.views import APIView
from rest_framework.response import Response

from hypothesis.models import (
    Vector,
    Vector3D,
    Interaction,
    InteractionPoint,
    TargetResidue,
)
from hypothesis.serializers import (
    VectorSerializer,
    Vector3DSerializer,
    InteractionSerializer,
    InteractionPointSerializer,
    TargetResidueSerialzier,
)
from viewer.models import SiteObservation

from .filters import Vector3DFilter

logger = logging.getLogger(__name__)


class VectorView(viewsets.ReadOnlyModelViewSet):
    queryset = Vector.objects.filter()
    serializer_class = VectorSerializer
    filterset_fields = ("cmpd_id", "smiles", "type")


# class Vector3DView(viewsets.ReadOnlyModelViewSet):
#     queryset = Vector3D.objects.filter()
#     serializer_class = Vector3DSerializer
#     filterset_fields = ("site_observation", "vector_id", "number")

    
class Vector3DView(viewsets.ReadOnlyModelViewSet):
    queryset = SiteObservation.objects.all()
    filterset_class = Vector3DFilter
    filterset_fields = ("site_observation", "number", "smiles", "cmpd_id", "vector_type")


    def list(self, request):
        """Method to handle GET request and call discourse to list posts for a topic
        """
        logger.info('+ Vector3DView.list')
        data = [k.get_vectors(**request.query_params.dict()) for k in self.queryset]
        serializer = Vector3DSerializer(chain.from_iterable(data), many=True)
        return Response(serializer.data)        
    

class InteractionView(viewsets.ReadOnlyModelViewSet):
    queryset = Interaction.objects.filter()
    serializer_class = InteractionSerializer
    filterset_fields = (
        "interaction_point",
        "interaction_version",
        "interaction_type",
        "interaction_point__prot_res_id__targ_res_id__target_id",
        "interaction_point__mol_id",
        "interaction_point__prot_res_id__prot_id",
        "distance",
        "score",
        "prot_smarts",
        "mol_smarts",
    )


class InteractionPointView(viewsets.ReadOnlyModelViewSet):
    queryset = InteractionPoint.objects.filter()
    serializer_class = InteractionPointSerializer
    filterset_fields = ("prot_res_id", "mol_id", "protein_atom_name", "molecule_atom_name")


class TargetResidueView(viewsets.ReadOnlyModelViewSet):
    queryset = TargetResidue.objects.filter()
    serializer_class = TargetResidueSerialzier
    filterset_fields = ("target_id", "res_name", "res_num", "chain_id")
