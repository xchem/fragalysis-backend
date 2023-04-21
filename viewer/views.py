import json
import os
import zipfile
from io import StringIO
import uuid
import shlex
import shutil
from datetime import datetime
from wsgiref.util import FileWrapper
from dateutil.parser import parse
import pytz

# import the logging library
import logging
import pandas as pd

from django.db import connections
from django.http import HttpResponse, FileResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.core.mail import send_mail
from django.conf import settings
from django.http import JsonResponse
from django.views import View

from rest_framework import viewsets
from rest_framework.parsers import BaseParser
from rest_framework.exceptions import ParseError
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework import permissions

from celery.result import AsyncResult

from api.security import ISpyBSafeQuerySet

from api.utils import get_params, get_highlighted_diffs
from viewer.utils import create_squonk_job_request_url

from viewer import models
from viewer.models import (
    Molecule,
    Protein,
    Project,
    Compound,
    Target,
    ActionType,
    SessionProject,
    SessionActions,
    Snapshot,
    SnapshotActions,
    ComputedMolecule,
    ComputedSet,
    NumericalScoreValues,
    ScoreDescription,
    TagCategory,
    TextScoreValues,
    MoleculeTag,
    SessionProjectTag,
    DownloadLinks,
    JobRequest,
    JobFileTransfer,
)
from viewer import filters
from viewer.squonk2_agent import Squonk2AgentRv, Squonk2Agent, get_squonk2_agent
from viewer.squonk2_agent import AccessParams, CommonParams, SendParams, RunJobParams

from .forms import CSetForm, TSetForm
from .tasks import (
    check_services,
    erase_compound_set_job_material,
    process_compound_set,
    process_design_sets,
    process_job_file_transfer,
    process_compound_set_job_file,
    process_target_set,
    validate_compound_set,
    validate_target_set,
)
from .discourse import create_discourse_post, list_discourse_posts_for_topic, check_discourse_user
from .download_structures import (
    check_download_links,
    recreate_static_file,
    maintain_download_links
)

from .squonk_job_file_transfer import (
    check_file_transfer
)

from .squonk_job_request import (
    check_squonk_active,
    get_squonk_job_config,
    create_squonk_job,
)

from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
    MolImageSerializer,
    CmpdImageSerializer,
    ProtMapInfoSerializer,
    ProtPDBInfoSerializer,
    ProtPDBBoundInfoSerializer,
    VectorsSerializer,
    GraphSerializer,
    ActionTypeSerializer,
    SessionProjectWriteSerializer,
    SessionProjectReadSerializer,
    SessionActionsSerializer,
    SnapshotReadSerializer,
    SnapshotWriteSerializer,
    SnapshotActionsSerializer,
    ComputedSetSerializer,
    ComputedMoleculeSerializer,
    NumericalScoreSerializer,
    ScoreDescriptionSerializer,
    TextScoreSerializer,
    ComputedMolAndScoreSerializer,
    DiscoursePostWriteSerializer,
    DictToCsvSerializer,
    TagCategorySerializer,
    MoleculeTagSerializer,
    SessionProjectTagSerializer,
    TargetMoleculesSerializer,
    DownloadStructuresSerializer,
    JobFileTransferReadSerializer,
    JobFileTransferWriteSerializer,
    JobRequestReadSerializer,
    JobRequestWriteSerializer,
    JobCallBackReadSerializer,
    JobCallBackWriteSerializer,
    ProjectSerializer,
    CompoundIdentifierSerializer,
    CompoundIdentifierTypeSerializer,
)

logger = logging.getLogger(__name__)

# Fields injected in a session object to pass
# messages between views. This is used by UploadCSet
# to pass errors and other messages back to the user
# via the upload-cset.html template.
_SESSION_ERROR = 'session_error'
_SESSION_MESSAGE = 'session_message'

_SQ2A: Squonk2Agent = get_squonk2_agent()

class CompoundIdentifierTypeView(viewsets.ModelViewSet):
    queryset = models.CompoundIdentifierType.objects.all()
    serializer_class = CompoundIdentifierTypeSerializer
    permission_classes = [permissions.IsAuthenticated]


class CompoundIdentifierTypeView(viewsets.ModelViewSet):
    queryset = models.CompoundIdentifier.objects.all()
    serializer_class = CompoundIdentifierSerializer
    permission_classes = [permissions.IsAuthenticated]
    filterset_fields = ["type", "compound"]


class VectorsView(ISpyBSafeQuerySet):
    """ DjagnoRF view for vectors

    Methods
    -------
    url:
        api/vector
    queryset:
        `viewer.models.Molecule.objects.filter()`
    filter fields:
        - `viewer.models.Molecule.prot_id` - ?prot_id=<int>
        - `viewer.models.Molecule.cmpd_id` - ?cmpd_id=<int>
        - `viewer.models.Molecule.smiles` - ?smiles=<str>
        - `viewer.models.Molecule.prot_id__target_id` - ?target_id=<int>
    returns:
        vectors for a given molecule generated by `frag.network.generate.get_3d_vects_for_mol()` (JSON)

    """
    queryset = Molecule.objects.filter()
    serializer_class = VectorsSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class GraphView(ISpyBSafeQuerySet):
    """ DjagnoRF view for graph

    Methods
    -------
    url:
        api/graph
    queryset:
        `viewer.models.Molecule.objects.filter()`
    filter fields:
        - `viewer.models.Molecule.prot_id` - ?prot_id=<int>
        - `viewer.models.Molecule.cmpd_id` - ?cmpd_id=<int>
        - `viewer.models.Molecule.smiles` - ?smiles=<str>
        - `viewer.models.Molecule.prot_id__target_id` - ?target_id=<int>
        - `viewer.models.Molecule.mol_groups` - ?mol_groups=<int>,<int>
    returns:
        graph network results for given molecule from `frag.network.query.get_full_graph()` (JSON)

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 385,
                    "graph": {
                        "CC(=O)Nc1cnccc1[Xe]_1_DELETION": {
                            "vector": "CC(O)NC1CCCCC1[101Xe]",
                            "addition": [
                                {
                                    "change": "C[101Xe]",
                                    "end": "CC(=O)Nc1cccnc1",
                                    "compound_ids": [
                                        "REAL:PV-001793547821",
                                        "MOLPORT:000-165-661"
                                    ]
                                }
                            ]
                        },
                        "C[Xe].NC(=O)C[Xe]_2_LINKER": {
                            "vector": "C[101Xe].NC(O)C[100Xe]",
                            "addition": [
                                {
                                    "change": "CNC1CC([100Xe])C(O)C1[101Xe]",
                                    "end": "CN=C1SC(CC(N)=O)C(=O)N1C",
                                    "compound_ids": [
                                        "MOLPORT:000-680-640"
                                    ]
                                },
                                {
                                    "change": "[100Xe]C1CCCCC1[101Xe]",
                                    "end": "CC1CCCCN1CC(N)=O",
                                    "compound_ids": [
                                        "REAL:Z54751033",
                                        "MOLPORT:001-599-191"
                                    ]
                                }
                            ]
                        },
                        "Cc1ccnc(Cl)c1[Xe]_2_REPLACE": {
                            "vector": "CC1CCCC(Cl)C1[100Xe]",
                            "addition": [
                                {
                                    "change": "CC(O)N[100Xe]",
                                    "end": "Cc1ccnc(Cl)c1",
                                    "compound_ids": [
                                        "MOLPORT:000-140-635"
                                    ]
                                }
                            ]
                        }
                    }
                }
        ]
    """
    queryset = Molecule.objects.filter()
    serializer_class = GraphSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class MolImageView(ISpyBSafeQuerySet):
    """ DjagnoRF view for molecule images

    Methods
    -------
    url:
        api/molimg
    queryset:
        `viewer.models.Molecule.objects.filter()`
    filter fields:
        - `viewer.models.Molecule.prot_id` - ?prot_id=<int>
        - `viewer.models.Molecule.cmpd_id` - ?cmpd_id=<int>
        - `viewer.models.Molecule.smiles` - ?smiles=<str>
        - `viewer.models.Molecule.prot_id__target_id` - ?target_id=<int>
        - `viewer.models.Molecule.mol_groups` - ?mol_groups=<int>,<int>
    returns:
        SVG image text for query molecule generated by `api.utils.draw_mol()` (JSON)

    example output:

        .. code-block:: javascript

             "results": [
                    {"id": 13912,
                        "mol_image": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' nk'..."}]
    """
    queryset = Molecule.objects.filter()
    serializer_class = MolImageSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class CompoundImageView(ISpyBSafeQuerySet):
    """ DjagnoRF view for compound images

   Methods
   -------
   url:
       api/cmpdimg
   queryset:
       `viewer.models.Compound.objects.filter()`
   filter fields:
       - `viewer.models.Molecule.smiles` - ?smiles=<str>
   returns:
       SVG image text for query compound generated by `api.utils.draw_mol()` (JSON)

   example output:

       .. code-block:: javascript

        "results": [
                {"id": 13912,
                    "mol_image": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' nk'..."}]

    """
    queryset = Compound.objects.filter()
    serializer_class = CmpdImageSerializer
    filter_permissions = "project_id"
    filterset_fields = ("smiles",)


class ProteinMapInfoView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve map info (file) for a given protein

   Methods
   -------
   url:
       api/protmap
   queryset:
       `viewer.models.Protein.objects.filter()`
   filter fields:
       - `viewer.models.Protein.code` - ?code=<str>
       - `viewer.models.Protein.target_id` - ?target_id=<int>
       - `viewer.models.Protein.prot_type` - ?prot_type=<str>
   returns:
       If a map file has been uploaded for the protein `map_info.path.read()` (JSON)

   """
    queryset = Protein.objects.filter()
    serializer_class = ProtMapInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProteinPDBInfoView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve apo pdb info (file) for a given protein

   Methods
   -------
   url:
       api/protpdb
   queryset:
       `viewer.models.Protein.objects.filter()`
   filter fields:
       - `viewer.models.Protein.code` - ?code=<str>
       - `viewer.models.Protein.target_id` - ?target_id=<int>
       - `viewer.models.Protein.prot_type` - ?prot_type=<str>
   returns: JSON
       - id: id of the protein object
       - pdb_data: If a pdb file has been uploaded for the protein `bound_info.path.read()`
       - prot_type: type of protein (e.g. AP for apo - see docs for model)

   example output:

       .. code-block:: javascript

           "results": [
            {
                "id": 27387,
                "pdb_data": "REMARK warning: chains may be ommitted for alignment REMARK ...",
                "prot_type": "AP"
            },]

   """
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProteinPDBBoundInfoView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve bound pdb info (file) for a given protein

   Methods
   -------
   url:
       api/protpdbbound
   queryset:
       `viewer.models.Protein.objects.filter()`
   filter fields:
       - `viewer.models.Protein.code` - ?code=<str>
       - `viewer.models.Protein.target_id` - ?target_id=<int>
       - `viewer.models.Protein.prot_type` - ?prot_type=<str>
   returns: JSON
       - id: id of the protein object
       - pdb_data: If a pdb file has been uploaded for the protein `bound_info.path.read()`
       - prot_type: type of protein (e.g. AP for apo - see docs for model)

   example output:

       .. code-block:: javascript

           "results": [
            {
                "id": 27387,
                "pdb_data": "REMARK warning: chains may be ommitted for alignment REMARK ...",
                "prot_type": "AP"
            },]

   """
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBBoundInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProjectView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve info about projects

       Methods
       -------
       url:
           api/projects
       queryset:
           `viewer.models.Project.objects.filter()`
       returns: JSON
           - id: id of the project object
           - title: name of the target
           - init_date: The date the Project was created

       example output:

           .. code-block:: javascript

               "results": [
                {
                    "id": 2,
                    "target_access_string": "lb27156-1",
                    "init_date": "2023-01-09T15:00:00Z",
                    "authority": "DIAMOND-ISPYB",
                    "open_to_public": false
                }
            ]

       """
    queryset = Project.objects.filter()
    serializer_class = ProjectSerializer
    # Special case - Project filter permissions is blank.
    filter_permissions = ""


class TargetView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve info about targets

       Methods
       -------
       url:
           api/targets
       queryset:
           `viewer.models.Target.objects.filter()`
       filter fields:
           - `viewer.models.Target.title` - ?title=<str>
       returns: JSON
           - id: id of the target object
           - title: name of the target
           - project_id: list of the ids of the projects the target is linked to
           - protein_set: list of the ids of the protein sets the target is linked to
           - template_protein: the template protein displayed in fragalysis front-end for this target
           - metadata: link to the metadata file for the target if it was uploaded
           - zip_archive: link to the zip archive of the uploaded data
           - default_squonk_project: project identifier of project on the Squonk application for
           this target.
           - upload_status: If set, this indicates the status of the most recent reload of the
           target data. Should normally move from 'PENDING' to 'STARTED' to 'SUCCESS".

       example output:

           .. code-block:: javascript

               "results": [
                {
                    "id": 62,
                    "title": "Mpro",
                    "project_id": [
                        2
                    ],
                    "protein_set": [
                        29281,
                        29274,
                        29259,
                        29305,
                        ...,
                    ],
                    "template_protein": "/media/pdbs/Mpro-x10417_0_apo.pdb",
                    "metadata": "http://fragalysis.diamond.ac.uk/media/metadata/metadata_2FdP5OJ.csv",
                    "zip_archive": "http://fragalysis.diamond.ac.uk/media/targets/Mpro.zip"
                }
            ]

       """
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    filter_permissions = "project_id"
    filterset_fields = ("title",)


class MoleculeView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve info about molecules

   Methods
   -------
   url:
       api/molecules
   queryset:
       `viewer.models.Molecule.objects.filter()`
   filter fields:
       - `viewer.models.Molecule.prot_id` - ?prot_id=<int>
       - `viewer.models.Molecule.cmpd_id` - ?cmpd_id=<int>
       - `viewer.models.Molecule.smiles` - ?smiles=<string>
       - `viewer.models.Molecule.prot_id__target_id` - ?target_id=<int>
       - `viewer.models.Molecule.mol_type` - ?mol_type=<str>
       - `viewer.models.Molecule.mol_groups` - ?mol_groups=<int>,<int>
   returns: JSON
       - id: id of the target object
       - smiles: smiles string of the molecule
       - cmpd_id: id of the related 2D compound object
       - prot_id: id of the related protein object
       - protein_code: code of the related protein object
       - mol_type: type of molecule - see Molecule model docs
       - molecule_protein: filepath of the apo protein structure for the molecule
       - lig_id: residue label for the ligand
       - chain_id: chain in the pdb file that the ligand belongs to
       - sdf_info: 3D coordinated of the molecule in MDL file format
       - x_com: x-coordinate for molecule centre of mass
       - y_com: y-coordinate for molecule centre of mass
       - z_com: z-coordinate for molecule centre of mass
       - mw: molecular weight
       - logp: LogP
       - tpsa: Topological Polar Surface Area
       - ha: heavy atom count
       - hacc: hydrogen-bond acceptors
       - hdon: hydrogen-bond donors
       - rots: number of rotatable bonds
       - rings: number of rings
       - velec: number of valence electrons

   example output:

       .. code-block:: javascript

           "results": [
            {
                "id": 13912,
                "smiles": "CN(C)c1ccc(C(=O)Nc2ccccn2)cc1",
                "cmpd_id": 796,
                "prot_id": 13923,
                "protein_code": "NUDT7A_Crude-x2226_2",
                "mol_type": "PR",
                "molecule_protein": "/media/pdbs/NUDT7A_Crude-x2226_2_apo_x5GxiLq.pdb",
                "lig_id": "LIG",
                "chain_id": "Z",
                "sdf_info": "     RDKit          3D 18 19  0  0  0  0  0  0  0  0999...",
                "x_com": null,
                "y_com": null,
                "z_com": null,
                "mw": 241.12,
                "logp": 2.4,
                "tpsa": 45.23,
                "ha": 18,
                "hacc": 3,
                "hdon": 1,
                "rots": 3,
                "rings": 2,
                "velec": 92
            },]

   """
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = (
        "prot_id",
        "prot_id__code",
        "cmpd_id",
        "smiles",
        "prot_id__target_id",
        "prot_id__target_id__title",
        "mol_type",
        "mol_groups",
    )


class CompoundView(ISpyBSafeQuerySet):
    """ DjagnoRF view for compound info

   Methods
   -------
   url:
       api/compounds
   queryset:
       `viewer.models.Compound.objects.filter()`
   filter fields:
       - `viewer.models.Molecule.smiles` - ?smiles=<str>
   returns:
       - id: id for compound object
       - inchi: inchi key for compound
       - smiles: smiles string for compound
       - mol_log_p: LogP for compound
       - num_h_acceptors: number of hydrogen-bond acceptors
       - num_h_donors: number of hydrogen-bond donors

   example output:

       .. code-block:: javascript

        "results": [
        {
            "id": 1,
            "inchi": "InChI=1S/C9H15NOS/c1-7(11)5-10-6-9-4-3-8(2)12-9/h3-4,7,10-11H,5-6H2,1-2H3/t7-/m0/s1",
            "smiles": "Cc1ccc(CNC[C@H](C)O)s1",
            "mol_log_p": 1.52692,
            "mol_wt": 185.0874351,
            "num_h_acceptors": 3,
            "num_h_donors": 2
        },]

    """
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_permissions = "project_id"
    filterset_fields = ("smiles", "current_identifier", "inchi", "long_inchi")


class ProteinView(ISpyBSafeQuerySet):
    """ DjagnoRF view to retrieve bound pdb info (file) for a given protein

   Methods
   -------
   url:
       api/proteins
   queryset:
       `viewer.models.Protein.objects.filter()`
   filter fields:
       - `viewer.models.Protein.code` - ?code=<str>
       - `viewer.models.Protein.target_id` - ?target_id=<int>
       - `viewer.models.Protein.prot_type` - ?prot_type=<str>
   returns: JSON
       - id: id of the protein object
       - code: the code/name of the protein
       - target_id: the id of the related target object
       - prot_type: the type of protein (e.g. AP for apo)
       - pdb_info: link to the apo pdb file
       - bound_info: link to the bound pdb file
       - mtz_info: link to the mtz file
       - map_info: link to the map file
       - cif_info: link to the cif file

   example output:

       .. code-block:: javascript

           "results": [
            {
                "id": 14902,
                "code": "XX02KALRNA-x1376_1",
                "target_id": 51,
                "prot_type": "AP",
                "pdb_info": "http://fragalysis.diamond.ac.uk/media/pdbs/XX02KALRNA-x1376_1_apo_9VSCvR8.pdb",
                "bound_info": "http://fragalysis.diamond.ac.uk/media/bound/XX02KALRNA-x1376_1_bound_6xmXkUm.pdb",
                "mtz_info": null,
                "map_info": null,
                "cif_info": null
            },]

   """
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


# START HERE! - THIS IS THE FIRST API THAT THE FRONT END CALLS.
def react(request):
    """
    :param request:
    :return: viewer/react page with context
    """

    discourse_api_key = settings.DISCOURSE_API_KEY

    context = {}

    # Is the Squonk2 Agent configured?
    logger.info("Checking whether Squonk2 is configured...")
    sq2_rv = _SQ2A.configured()
    if sq2_rv.success:
        logger.info("Squonk2 is configured")
        context['squonk_available'] = 'true'
    else:
        logger.info("Squonk2 is NOT configured")
        context['squonk_available'] = 'false'

    if discourse_api_key:
        context['discourse_available'] = 'true'
    else:
        context['discourse_available'] = 'false'

    user = request.user
    if user.is_authenticated:
        context['discourse_host'] = ''
        context['user_present_on_discourse'] = 'false'
        # If user is authenticated and a discourse api key is available, then check discourse to
        # see if user is set up and set up flag in context.
        if discourse_api_key:
            context['discourse_host'] = settings.DISCOURSE_HOST
            error, error_message, user_id = check_discourse_user(user)
            if user_id:
                context['user_present_on_discourse'] = 'true'

        # If user is authenticated Squonk can be called then return the Squonk host
        # so the Frontend can navigate to it
        context['squonk_ui_url'] = ''
        if sq2_rv.success and check_squonk_active(request):
            context['squonk_ui_url'] = _SQ2A.get_ui_url()

    return render(request, "viewer/react_temp.html", context)

# Upload Compound set functions


def save_pdb_zip(pdb_file):
    zf = zipfile.ZipFile(pdb_file)
    zip_lst = zf.namelist()
    zfile = {}
    zfile_hashvals = {}
    print(zip_lst)
    for filename in zip_lst:
        # only handle pdb files
        if filename.split('.')[-1] == 'pdb':
            f = filename.split('/')[0]
            save_path = os.path.join(settings.MEDIA_ROOT, 'tmp', f)
            if default_storage.exists(f):
                rand_str = uuid.uuid4().hex
                pdb_path = default_storage.save(save_path.replace('.pdb', f'-{rand_str}.pdb'), ContentFile(zf.read(filename)))
            # Test if Protein object already exists
            # code = filename.split('/')[-1].replace('.pdb', '')
            # test_pdb_code = filename.split('/')[-1].replace('.pdb', '')
            # test_prot_objs = Protein.objects.filter(code=test_pdb_code)
            #
            # if len(test_prot_objs) > 0:
            #     # make a unique pdb code as not to overwrite existing object
            #     rand_str = uuid.uuid4().hex
            #     test_pdb_code = f'{code}#{rand_str}'
            #     zfile_hashvals[code] = rand_str
            #
            # fn = test_pdb_code + '.pdb'
            #
            # pdb_path = default_storage.save('tmp/' + fn,
            #                                 ContentFile(zf.read(filename)))
            else:
                pdb_path = default_storage.save(save_path, ContentFile(zf.read(filename)))
            test_pdb_code = pdb_path.split('/')[-1].replace('.pdb', '')
            zfile[test_pdb_code] = pdb_path

    # Close the zip file
    if zf:
        zf.close()

    return zfile, zfile_hashvals


def save_tmp_file(myfile):
    """ Save file in temporary location for validation/upload processing
    """

    name = myfile.name
    path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
    tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

    return tmp_file


class UploadCSet(View):
    """ View to render and control viewer/upload-cset.html  - a page allowing upload of computed sets. Validation and
    upload tasks are defined in `viewer.compound_set_upload`, `viewer.sdf_check` and `viewer.tasks` and the task
    response handling is done by `viewer.views.ValidateTaskView` and `viewer.views.UploadTaskView`

    Methods
    -------
    allowed requests:
        - GET: renders form
        - POST: validates, deletes or optionally uploads the computed set that the user
                provides via the template form
    url:
        viewer/upload_cset
    template:
        viewer/upload-cset.html
    request params:
        - target_name (`django.forms.CharField`): Name of the existing fragalysis target to add the computed set to
        - sdf_file (`django.forms.FileField`): SDF file of all computed molecules to upload for the computed set
        - pdb_zip (`django.forms.FileField`): zip file of apo pdb files referenced in the ref_pdb field for molecules in sdf_file (optional)
        - submit_choice (`django.forms.CharField`): validate, validate and upload, delete
    context:
        - form (`django.Forms.form`): instance of `viewer.forms.CSetForm`
        - validate_task_id (str): celery task id for validation step
        - validate_task_status (str): celery task status for validation step
        - upload_task_id (str): celery task id for upload step
        - upload_task_status (str): celery task status for upload step
    """

    def get(self, request):

        # Any messages passed to us via the session?
        # Maybe from a redirect?
        # It so take them and remove them.
        session_error = None
        if _SESSION_ERROR in request.session:
            session_error = request.session[_SESSION_ERROR]
            del request.session[_SESSION_ERROR]
        session_message = None
        if _SESSION_MESSAGE in request.session:
            session_message = request.session[_SESSION_MESSAGE]
            del request.session[_SESSION_MESSAGE]

        # Only authenticated users can upload files
        # - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context = {}
            context['error_message'] \
                = 'Only authenticated users can upload files' \
                  ' - please navigate to landing page and Login'
            return render(request, 'viewer/upload-cset.html', context)

        form = CSetForm()
        existing_sets = ComputedSet.objects.all()
        context = {'form': form,
                   'sets': existing_sets,
                   _SESSION_ERROR: session_error,
                   _SESSION_MESSAGE: session_message}
        return render(request, 'viewer/upload-cset.html', context)

    def post(self, request):

        # Only authenticated users can upload files
        # - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context = {}
            context['error_message'] \
                = 'Only authenticated users can upload files' \
                  ' - please navigate to landing page and Login'
            return render(request, 'viewer/upload-cset.html', context)

        # Celery/Redis must be running.
        # This call checks and trys to start them if they're not.
        assert check_services()

        form = CSetForm(request.POST, request.FILES)

        if form.is_valid():

            # Get all the variables needed from the form.
            # The fields we use will be based on the 'submit_choice',
            # expected to be one of V (validate), U (upload) or D delete
            choice = request.POST['submit_choice']

            # Generate run-time error if the required form fields
            # are not set based on the choice made...

            # The 'sdf_file' anf 'target_name' are only required for upload/update
            sdf_file = request.FILES.get('sdf_file')
            target = request.POST.get('target_name')
            update_set = request.POST['update_set']

            # If a set is named the ComputedSet cannot be 'Anonymous'
            # and the user has to be the owner.
            selected_set = None
            if update_set != 'None':
                computed_set_query = ComputedSet.objects.filter(unique_name=update_set)
                if computed_set_query:
                    selected_set = computed_set_query[0]
                else:
                    request.session[_SESSION_ERROR] = \
                        'The set could not be found'
                    return redirect('upload_cset')

            # If validating or uploading we need a Target and SDF file.
            # If updating or deleting we need an update set (that's not 'None')
            if choice in ['V', 'U']:
                if sdf_file is None or target is None:
                    request.session[_SESSION_ERROR] = \
                        'To Validate or Upload' \
                        ' you must provide a Target and SDF file'
                    return redirect('upload_cset')
            elif choice in ['D']:
                if update_set == 'None':
                    request.session[_SESSION_ERROR] = \
                        'To Delete you must select an existing set'
                    return redirect('upload_cset')

            # If uploading (updating) or deleting
            # the set owner cannot be anonymous
            # and the user needs to be the owner
            if choice in ['U', 'D'] and selected_set:
                if selected_set.owner_user.id == settings.ANONYMOUS_USER:
                    request.session[_SESSION_ERROR] = \
                        'You cannot Update or Delete Anonymous sets'
                elif selected_set.owner_user != user:
                    request.session[_SESSION_ERROR] = \
                        'You can only Update or Delete sets you own'
                # Something wrong?
                # If so redirect...
                if _SESSION_ERROR in request.session:
                    return redirect('upload_cset')

            # Save uploaded sdf and zip to tmp storage
            tmp_pdb_file = None
            tmp_sdf_file = None
            if 'pdb_zip' in list(request.FILES.keys()):
                pdb_file = request.FILES['pdb_zip']
                tmp_pdb_file = save_tmp_file(pdb_file)
            if sdf_file:
                tmp_sdf_file = save_tmp_file(sdf_file)

            if choice == 'V':
                # Validate
                # Start celery task
                task_params = {'user_id': user.id,
                               'sdf_file': tmp_sdf_file,
                               'target': target}
                if tmp_pdb_file:
                    task_params['zfile'] = tmp_pdb_file
                if update_set:
                    task_params['update'] = update_set
                task_validate = validate_compound_set.delay(task_params)

                # Update client side with task id and status
                context = {'validate_task_id': task_validate.id,
                           'validate_task_status': task_validate.status}
                return render(request, 'viewer/upload-cset.html', context)

            elif choice == 'U':
                # Upload
                # Start chained celery tasks. NB first function passes tuple
                # to second function - see tasks.py
                task_params = {'user_id': user.id,
                               'sdf_file': tmp_sdf_file,
                               'target': target}
                if tmp_pdb_file:
                    task_params['zfile'] = tmp_pdb_file
                if update_set:
                    task_params['update'] = update_set
                task_upload = (
                        validate_compound_set.s(task_params) |
                        process_compound_set.s()).apply_async()

                # Update client side with task id and status
                context = {'upload_task_id': task_upload.id,
                           'upload_task_status': task_upload.status}
                return render(request, 'viewer/upload-cset.html', context)

            elif choice == 'D':
                # Delete
                selected_set.delete()

                request.session[_SESSION_MESSAGE] = \
                    f'Compound set "{selected_set.unique_name}" deleted'
                return redirect('upload_cset')

        context = {'form': form}
        return render(request, 'viewer/upload-cset.html', context)
# End Upload Compound set functions


# Upload Target datasets functions
class UploadTSet(View):
    """ View to render and control viewer/upload-tset.html  - a page allowing upload of computed sets. Validation and
    upload tasks are defined in `viewer.target_set_upload`, `viewer.sdf_check` and `viewer.tasks` and the task
    response handling is done by `viewer.views.ValidateTaskView` and `viewer.views.UploadTaskView`

    Methods
    -------
    allowed requests:
        - GET: renders form
        - POST: validates and optionally uploads the computed set that the user provides via the template form
    url:
        viewer/upload_tset
    template:
        viewer/upload-tset.html
    request params:
        - target_name (`django.forms.CharField`): Name of the existing fragalysis target to add the computed set to
        - target_zip (`django.forms.FileField`): zip file of the target dataset
        - submit_choice (`django.forms.CharField`): validate, validate and upload
    context:
        - form (`django.Forms.form`): instance of `viewer.forms.TSetForm`
        - validate_task_id (str): celery task id for validation step
        - validate_task_status (str): celery task status for validation step
        - upload_task_id (str): celery task id for upload step
        - upload_task_status (str): celery task status for upload step

    """

    def get(self, request):

        # Only authenticated users can upload files - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context = {}
            context['error_message'] \
                = 'Only authenticated users can upload files - please navigate to landing page and Login'
            return render(request, 'viewer/upload-tset.html', context)

        contact_email = ''
        if user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            contact_email = user.email

        form = TSetForm(initial={'contact_email': contact_email})

        return render(request, 'viewer/upload-tset.html', {'form': form})

    def post(self, request):
        logger.info('+ UploadTSet.post')
        context = {}

        # Only authenticated users can upload files - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context['error_message'] \
                = 'Only authenticated users can upload files - please navigate to landing page and Login'
            logger.info('- UploadTSet.post - authentication error')
            return render(request, 'viewer/upload-tset.html', context)

        # Celery/Redis must be running.
        # This call checks and trys to start them if they're not.
        assert check_services()

        form = TSetForm(request.POST, request.FILES)
        if form.is_valid():
            # get all of the variables needed from the form
            target_file = request.FILES['target_zip']
            target_name = request.POST['target_name']
            choice = request.POST['submit_choice']
            proposal_ref = request.POST['proposal_ref']
            contact_email = request.POST['contact_email']

            # Create /code/media/tmp if does not exist
            media_root = settings.MEDIA_ROOT
            tmp_folder = os.path.join(media_root, 'tmp')
            if not os.path.isdir(tmp_folder):
                os.mkdir(tmp_folder)

            path = default_storage.save('tmp/' + 'NEW_DATA.zip', ContentFile(target_file.read()))
            new_data_file = str(os.path.join(settings.MEDIA_ROOT, path))

            # Settings for if validate option selected
            if choice == 'V':
                # Start celery task
                task_validate = validate_target_set.delay(new_data_file, target=target_name, proposal=proposal_ref,
                                                          email=contact_email)

                # Update client side with task id and status
                context = {'validate_task_id': task_validate.id,
                           'validate_task_status': task_validate.status}
                return render(request, 'viewer/upload-tset.html', context)

            # if it's an upload, run the validate followed by the upload target set task
            if choice == 'U':
                # Start chained celery tasks. NB first function passes tuple
                # to second function - see tasks.py
                task_upload = (validate_target_set.s(new_data_file, target=target_name, proposal=proposal_ref,
                                                     email=contact_email) | process_target_set.s()).apply_async()

                # Update client side with task id and status
                context = {'upload_task_id': task_upload.id,
                           'upload_task_status': task_upload.status}
                return render(request, 'viewer/upload-tset.html', context)

        context = {'form': form}
        return render(request, 'viewer/upload-tset.html', context)


# End Upload Target datasets functions
def email_task_completion(contact_email, message_type, target_name, target_path=None, task_id=None):
    """Notify user of upload completion
    """

    logger.info('+ email_notify_task_completion: ' + message_type + ' ' + target_name)
    email_from = settings.EMAIL_HOST_USER

    if contact_email == '' or not email_from:
        # Only send email if configured.
        return

    if message_type == 'upload-success':
        subject = 'Fragalysis: Target: '+target_name+' Uploaded'
        message = 'The upload of your target data is complete. Your target is available at: ' \
                  + str(target_path)
    elif message_type == 'validate-success':
        subject = 'Fragalysis: Target: '+target_name+' Validation'
        message = 'Your data was validated. It can now be uploaded using the upload option.'
    else:
        # Validation failure
        subject = 'Fragalysis: Target: ' + target_name + ' Validation/Upload Failed'
        message = 'The validation/upload of your target data did not complete successfully. ' \
                  'Please navigate the following link to check the errors: validate_task/' + str(task_id)

    recipient_list = [contact_email, ]
    logger.info('+ email_notify_task_completion email_from: %s', email_from )
    logger.info('+ email_notify_task_completion subject: %s',  subject )
    logger.info('+ email_notify_task_completion message: %s',  message )
    logger.info('+ email_notify_task_completion contact_email: %s', contact_email )

    # Send email - this should not prevent returning to the screen in the case of error.
    send_mail(subject, message, email_from, recipient_list, fail_silently=True)
    logger.info('- email_notify_task_completion')
    return


# Task functions common between Compound Sets and Target Set pages.
class ValidateTaskView(View):
    """ View to handle dynamic loading of validation results from `viewer.tasks.validate` - the validation of files
    uploaded to viewer/upload_cset or a target set by a user at viewer/upload_tset

    Methods
    -------
    allowed requests:
        - GET: takes a task id, checks it's status and returns the status, and result if the task is complete
    url:
        validate_task/<validate_task_id>
    template:
        viewer/upload-cset.html or viewer/upload-tset.html
    """
    def get(self, request, validate_task_id):
        """ Get method for `ValidateTaskView`. Takes a validate task id, checks it's status and returns the status,
        and result if the task is complete

        Parameters
        ----------
        request: request
            Context sent by `UploadCSet` or `UploadTset`
        validate_task_id: str
            task id provided by `UploadCSet` or `UploadTset`

        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                - if status = 'FAILURE':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - validate_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - html (str): html of task outcome - success message or html table of errors & fail message

        """
        logger.info('+ ValidateTaskView.get')
        task = AsyncResult(validate_task_id)
        response_data = {'validate_task_status': task.status,
                         'validate_task_id': task.id}

        if task.status == 'FAILURE':
            logger.info('+ ValidateTaskView.get.FAILURE')
            result = task.traceback
            response_data['validate_traceback'] = str(result)

            return JsonResponse(response_data)

        # Check if results ready
        if task.status == "SUCCESS":
            logger.info('+ ValidateTaskView.get.SUCCESS')
            results = task.get()
            # NB get tuple from validate task
            process_type = results[1]
            validate_dict = results[2]
            validated = results[3]
            if validated:
                response_data['html'] = 'Your data was validated. \n It can now be uploaded using the upload option.'
                response_data['validated'] = 'Validated'

                if process_type== 'tset':
                    target_name = results[5]
                    contact_email = results[8]
                    email_task_completion(contact_email, 'validate-success', target_name)

                return JsonResponse(response_data)

            if not validated:
                # set pandas options to display all column data
                pd.set_option('display.max_colwidth', -1)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                response_data["html"] = html_table
                response_data['validated'] = 'Not validated'
                if process_type== 'tset':
                    target_name = results[5]
                    contact_email = results[8]
                    email_task_completion(contact_email, 'validate-failure', target_name, task_id=validate_task_id)

                return JsonResponse(response_data)

        return JsonResponse(response_data)


class UpdateTaskView(View):

    def get(self, request, update_task_id):
        task = AsyncResult(update_task_id)
        response_data = {'update_task_status': task.status,
                         'update_task_id': task.id}

        result = 'Running...'

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['result'] = str(result)

        if task.status == 'SUCCESS':
            result = task.get()

        response_data['result'] = str(result)

        return JsonResponse(response_data)


class UploadTaskView(View):
    """ View to handle dynamic loading of upload results from `viewer.tasks.process_compound_set` - the upload of files
    for a computed set by a user at viewer/upload_cset or a target set by a user at viewer/upload_tset

    Methods
    -------
    allowed requests:
        - GET: takes a task id, checks it's status and returns the status, and result if the task is complete
    url:
        upload_task/<uploads_task_id>
    template:
        viewer/upload-cset.html or viewer/upload-tset.html
    """
    def get(self, request, upload_task_id):
        """ Get method for `UploadTaskView`. Takes an upload task id, checks it's status and returns the status,
        and result if the task is complete

        Parameters
        ----------
        request: request
            Context sent by `UploadCSet` or `UploadTSet`
        upload_task_id: str
            task id provided by `UploadCSet` or `UploadTSet`

        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                - if status = 'FAILURE':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - upload_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - if results are a list (data was processed - validated or uploaded):
                        if this was a validation process
                        - validated (str): 'Not validated'
                        - html (str): html table of validation errors
                        if results are a validation/upload process:
                        - validated (str): 'Validated'
                        - results (dict): results
                        For compound sets ('cset')
                        - results['cset_download_url'] (str): download url for computed set sdf file
                        - results['pset_download_url'] (str): download url for computed set pdb files (zip)
                        For target sets ('tset')
                        - results['tset_download_url'] (str): download url for processed zip file
                    - if results are not string or list:
                        - processed (str): 'None'
                        - html (str): message to tell the user their data was not processed

        """
        logger.debug('+ UploadTaskView.get')
        task = AsyncResult(upload_task_id)
        response_data = {'upload_task_status': task.status,
                         'upload_task_id': task.id}

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['upload_traceback'] = str(result)

            return JsonResponse(response_data)

        if task.status == 'SUCCESS':
            logger.debug('+ UploadTaskView.get.success')

            results = task.get()

            # Validation output for a cset or tset is a dictionary.
            if isinstance(results, list):
                if results[0] == 'validate':
                    # Get dictionary results
                    validate_dict = results[1]

                    # set pandas options to display all column data
                    pd.set_option('display.max_colwidth', -1)
                    table = pd.DataFrame.from_dict(results[2])
                    html_table = table.to_html()
                    html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                    response_data['validated'] = 'Not validated'
                    response_data['html'] = html_table

                    return JsonResponse(response_data)
                else:
                    # Upload/Update output tasks send back a tuple
                    # First element defines the source of the upload task (cset, tset)
                    response_data['validated'] = 'Validated'
                    if results[1] == 'tset':
                        target_name = results[2]
                        contact_email = results[5]
                        target_path = '/viewer/target/%s' % target_name
                        response_data['results'] = {}
                        response_data['results']['tset_download_url'] = target_path
                        logger.info('+ UploadTaskView.get.success -email: %s', contact_email)
                        email_task_completion(contact_email, 'upload-success', target_name, target_path=target_path)
                    else:
                        cset_name = results[2]
                        cset = ComputedSet.objects.get(name=cset_name)
                        submitter = cset.submitter
                        name = cset.unique_name
                        response_data['results'] = {}
                        response_data['results']['cset_download_url'] = '/viewer/compound_set/%s' % name
                        response_data['results']['pset_download_url'] = '/viewer/protein_set/%s' % name

                    return JsonResponse(response_data)

            else:
                # Error output
                html_table = '''<p> Your data was <b>not</b> processed.</p>'''
                response_data['processed'] = 'None'
                response_data['html'] = html_table
                return JsonResponse(response_data)

        return JsonResponse(response_data)
# End Task functions which hopefully can be common between Compound Sets and Target Set pages.


def img_from_smiles(request):
    """ View to generate a 2D molecule image for a given smiles string

    Methods
    -------
    allowed requests:
        - GET: generate a 2D molecule image for a given smiles string
    url:
        viewer/img_from_smiles
    request params:
        - smiles (str): smiles string to generate image for

    Returns
    -------
    HTTPResponse (str):
        - if smiles provided:
            string for SVG image of molecule
        - if smiles not provided:
            "Please insert SMILES"

    """
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        if smiles:
            return get_params(smiles, request)
        else:
            return HttpResponse("Please insert SMILES")
    else:
        return HttpResponse("Please insert SMILES")


def highlight_mol_diff(request):
    """ View to generate a 2D molecule image highlighting the difference between a reference and new molecule

    Methods
    -------
    allowed requests:
        - GET: generate a 2D molecule image highlighting the difference between a reference and new molecule
    url:
        viewer/highlight_mol_diff
    request params:
        - prb_smiles (str): smiles string to generate image for
        - ref_smiles (str): reference smiles for highlighting by MCS

    Returns
    -------
    HTTPResponse (str):
        - if smiles provided:
            string for SVG image of molecule
        - if smiles not provided:
            "Please insert smiles for reference and probe"

    """
    if 'prb_smiles' and 'ref_smiles' in request.GET:
        return HttpResponse(get_highlighted_diffs(request))
    else:
        return HttpResponse("Please insert smiles for reference and probe")


def similarity_search(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
    else:
        return HttpResponse("Please insert SMILES")
    if "db_name" in request.GET:
        db_name = request.GET["db_name"]
    else:
        return HttpResponse("Please insert db_name")
    sql_query = """SELECT sub.*
  FROM (
    SELECT rdk.id,rdk.structure,rdk.idnumber
      FROM vendordbs.enamine_real_dsi_molfps AS mfp
      JOIN vendordbs.enamine_real_dsi AS rdk ON mfp.id = rdk.id
      WHERE m @> qmol_from_smiles(%s) LIMIT 1000
  ) sub;"""
    with connections[db_name].cursor() as cursor:
        cursor.execute(sql_query, [smiles])
        return HttpResponse(json.dumps(cursor.fetchall()))


def get_open_targets(request):
    """ View to return a list of all open targets

    Methods
    -------
    allowed requests:
        - GET: return a list of all open targets
    url:
        viewer/open_targets
    request params:
        None

    Returns
    -------
    HTTPResponse (JSON/dict):
        - target_names (list): list of open targets
        - target_ids (list): list of ids for open targets in same order as target_names

    """
    targets = Target.objects.all()
    target_names = []
    target_ids = []

    for t in targets:
        for p in t.project_id.all():
            if 'OPEN' in p.title:
                target_names.append(t.title)
                target_ids.append(t.id)

    return HttpResponse(json.dumps({'target_names': target_names, 'target_ids': target_ids}))


# This is used in the URL on the process results page after uploading a compound_set
def cset_download(request, name):
    """ View to download an SDF file of a computed set by name

    Methods
    -------
    allowed requests:
        - GET: return the SDF file as a download
    url:
        viewer/compound_set/(<name>)
    request params:
        - name (str): the name of the computed set to download

    Returns
    -------
    Response (attachment; text/plain):
        - <name>.sdf: sdf file for the computed set

    """
    compound_set = ComputedSet.objects.get(unique_name=name)
    filepath = compound_set.submitted_sdf
    with open(filepath.path, 'r', encoding='utf-8') as fp:
        data = fp.read()
    filename = 'compund-set_' + name + '.sdf'
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s' % filename  # force browser to download file
    response.write(data)
    return response


def pset_download(request, name):
    """ View to download a zip file of all protein structures (apo) for a computed set

     Methods
     -------
     allowed requests:
         - GET: return the zip file as a download
     url:
         viewer/compound_set/(<name>)
     request params:
         - name (str): the name of the computed set to download

     Returns
     -------
     Response (attachment; application/zip):
         - <name>.zip: zip file for the computed set

     """
    response = HttpResponse(content_type='application/zip')
    filename = 'protein-set_' + name + '.zip'
    response['Content-Disposition'] = 'filename=%s' % filename  # force browser to download file

    compound_set = ComputedSet.objects.get(unique_name=name)
    computed = ComputedMolecule.objects.filter(computed_set=compound_set)
    pdb_filepaths = list(set([c.pdb_info.path for c in computed]))

    buff = StringIO()
    zip_obj = zipfile.ZipFile(buff, 'w')

    for fp in pdb_filepaths:
        data = open(fp, 'r', encoding='utf-8').read()
        zip_obj.writestr(fp.split('/')[-1], data)
    zip_obj.close()

    buff.flush()
    ret_zip = buff.getvalue()
    buff.close()
    response.write(ret_zip)

    return response


# This is used in the URL on the process results page after uploading a target_set
def tset_download(request, title):
    """ View to download an zip file of a target set by name

    Methods
    -------
    allowed requests:
        - GET: return the zip file as a download
    url:
        viewer/target/(<title>)
    request params:
        - title (str): the title of the target set to download

    Returns
    -------
    Response (attachment; text/plain):
        - <title>.zip: zip file for the target set

    """
    target_set = Target.objects.get(title=title)
    media_root = settings.MEDIA_ROOT
    filepath = os.path.join(media_root, target_set.zip_archive.name)
    target_zip = open(filepath, 'rb')
    filename = 'target-set_' + title + '.zip'
    response = HttpResponse(target_zip, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename="%s"' % filename  # force browser to download file
    return response


# Start of ActionType
class ActionTypeView(viewsets.ModelViewSet):
    """ Djagno view to retrieve information about action types available to users (GET)

    Methods
    -------
    url:
        api/action-types
    queryset:
        `viewer.models.ActionType.objects.filter()`
    methods:
        `get, head`
    filter fields:
        - `viewer.models.ActionType.description` - ?description=<str>
        - `viewer.models.ActionType.active` - ?active=<boolean>
        - `viewer.models.ActionType.activation_date` - ?date=<str>

    returns: JSON
        - id: id of the action type
        - description: The description of the action type
        - active: True if action type is active, False if not
        - activation_date: The datetime the action type became active

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 1,
                    "description": "Test Action TYpe",
                    "active": true,
                    "activation_date": "2020-10-06T14:42:00Z"
                }
            ]

   """
    queryset = ActionType.objects.filter()
    serializer_class = ActionTypeSerializer

    # POST method allowed for flexibility in the PoC. In the final design we may want to prevent the POST/PUT methods
    # from being used
    # for action types so that these can only be updated via the admin panel.
    #    http_method_names = ['get', 'head']

    filterset_fields = '__all__'


# Start of Session Project
class SessionProjectsView(viewsets.ModelViewSet):
    """ Djagno view to retrieve information about user projects (collection of sessions) (GET). Also used for saving
    project information (PUT, POST, PATCH)

    Methods
    -------
    url:
        api/session-projects
    queryset:
        `viewer.models.SessionProject.objects.filter()`
    filter fields:
        - `viewer.models.SessionProject.title` - ?title=<str>
        - `viewer.models.SessionProject.init_date` - ?date=<str>
        - `viewer.models.SessionProject.description` - ?description=<str>
        - `viewer.models.SessionProject.target` - ?target=<int>
        - `viewer.models.SessionProject.author` - ?author=<str>
        - `viewer.models.SessionProject.tags` - ?tags=<list>
    returns: JSON
        - id: id of the project
        - target: dict of target info:
            - id: target id
            - title: name of the target protein
            - project_id: id of the project
            - protein_set: list of protein objects associated with the target
            - template_protein: link to the file used as the reference in fragalysis frontend
            - metadata: link to the target metadata file
            - zip_archive: link to the zip archive of uploaded files
        - author: name of the author that created the session
        - title: title for the project
        - init_date: timestamp for when the project was created
        - description: author defined description of the project
        - tags: list of tags given to the project by the author

    example output:

        .. code-block:: javascript

            "results": [
            {
                "id": 122,
                "target": {
                    "id": 62,
                    "title": "Mpro",
                    "project_id": [
                        2
                    ],
                    "protein_set": [
                        29281,
                        29274,
                        29259,
                        29305,
                        29250,
                        ...,
                    ],
                    "template_protein": "/media/pdbs/Mpro-x10417_0_apo.pdb",
                    "metadata": "http://fragalysis.diamond.ac.uk/media/metadata/metadata_2FdP5OJ.csv",
                    "zip_archive": "http://fragalysis.diamond.ac.uk/media/targets/Mpro.zip"
                },
                "author": null,
                "title": "READ_ONLY",
                "init_date": "2020-07-09T19:52:10.506119Z",
                "description": "READ_ONLY",
                "tags": "[]",
                "session_project_tags": [
                {
                    "id": 3,
                    "tag": "testtag3",
                    "category_id": 1,
                    "target_id": 1,
                    "user_id": null,
                    "create_date": "2021-04-13T16:01:55.396088Z",
                    "colour": null,
                    "discourse_url": null,
                    "help_text": null,
                    "additional_info": ""
                },]
            },]

   """
    queryset = SessionProject.objects.filter()

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.SessionProjectReadSerializer`
            - if other: `viewer.serializers.SessionProjectWriteSerializer`
        """
        if self.request.method in ['GET']:
            # GET
            return SessionProjectReadSerializer
        # (POST, PUT, PATCH)
        return SessionProjectWriteSerializer

    filter_permissions = "target_id__project_id"
    filterset_fields = '__all__'


class SessionActionsView(viewsets.ModelViewSet):
    """ Djagno view to retrieve information about actions relating to sessions_project (GET). Also used for saving
    project action information (PUT, POST, PATCH)

    Methods
    -------
    url:
        api/session-actions
    queryset:
        `viewer.models.SessionActions.objects.filter()`
    filter fields:
        - `viewer.models.SessionProject.author` - ?author=<str>
        - `viewer.models.SessionProject.session_project` - ?project=<str>
        - `viewer.models.SessionProject.last_update_date` - ?date=<str>

    returns: JSON
        - id: id of the session action record
        - author: id of the user that created the session_project
        - session_project: id of the related session_project
        - last_update_date: Timestamp for when the action list was generated or updated
        - actions: JSON string containing actions related to the session_project

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 1,
                    "last_update_date": "2020-10-06T15:36:00Z",
                    "actions": {
                        "save": false,
                        "show": true,
                        "action_type": 1,
                        "object_name": "",
                        "object_type": "",
                        "action_datetime": "2020-09-30T13:44:00.000Z"
                    },
                    "author": 1,
                    "session_project": 1,
                }
            ]

   """
    queryset = SessionActions.objects.filter()
    serializer_class = SessionActionsSerializer

    #   Note: jsonField for Actions will need specific queries - can introduce if needed.
    filterset_fields = ('id', 'author', 'session_project', 'last_update_date')


class SnapshotsView(viewsets.ModelViewSet):
    """ Djagno view to retrieve information about user sessions (snapshots) (GET). Also used for saving
    session information (PUT, POST, PATCH)

    Methods
    -------
    url:
        api/snapshots
    queryset:
        `viewer.models.Snapshot.objects.filter()`
    filter class:
        `viewer.filters.SnapshotFilter`
    filter fields:
        - `viewer.models.Snapshot.id` - ?id=<int>
        - `viewer.models.Snapshot.type` - ?type=<str>
        - `viewer.models.Snapshot.author` - ?author=<str>
        - `viewer.models.Snapshot.description` ?description=<str>
        - `viewer.models.Snapshot.created` - ?created=<str>
        - `viewer.models.Snapshot.data` - ?data=<JSON str>
        - `viewer.models.Snapshot.session_project_id` - ?session_project_id=<int>
        - `viewer.models.Snapshot.parent` - ?parent=<int>
        - `viewer.models.Snapshot.children` - ?children=<list>
        - `viewer.models.Snapshot.session_project` - ?session_project=<int>
    returns: JSON
        - id: id of the snapshot
        - type: type of snapshot
        - title: title of snapshot given by author
        - description: description of the snapshot given by author
        - created: timestamp for when the snapshot was created
        - data: json string describing data needed to reproduce state of the snapshot by the front-end
        - session_project: dict describing the project that the snapshot belongs to:
            - id: project id
            - title: project title
            - init_data: timestamp for when the project was initiated
            - description: description of the project given by the author
            - tags: tags given to the project by the author
            - target: id of the target that the project is associated with
            - author: name of the author who created the project
        - parent: parent snapshot id of the current snapshot
        - children: list of children ids of the current snapshot
        - additional_info: Free format json for use by the Fragalysis frontend.

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 132,
                    "type": "INIT",
                    "title": "-- 2020-07-09 -- 16:01:35",
                    "author": null,
                    "description": "Snapshot generated by anonymous user",
                    "created": "2020-07-09T20:01:36.901552Z",
                    "data": '"{\"apiReducers\":{\"target_id_list\":[{\"id\":2,\"title\":\"NUDT7A\",\'
                    "session_project": {
                        "id": 124,
                        "title": "READ_ONLY",
                        "init_date": "2020-07-09T20:01:35.715324Z",
                        "description": "READ_ONLY",
                        "tags": "[]",
                        "target": 62,
                        "author": null
                    },
                    "parent": null,
                    "children": []
                    "additional_info": []
                },]

   """
    queryset = Snapshot.objects.filter()

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.SnapshotReadSerializer`
            - if other: `viewer.serializers.SnapshotWriteSerializer`
        """
        if self.request.method in ['GET']:
            return SnapshotReadSerializer
        return SnapshotWriteSerializer

    filter_class = filters.SnapshotFilter


class SnapshotActionsView(viewsets.ModelViewSet):
    """ Djagno view to retrieve information about actions relating to snapshots (GET). Also used for saving
    snapshot action information (PUT, POST, PATCH)

    Methods
    -------
    url:
        api/snapshot-actions
    queryset:
        `viewer.models.SnapshotActions.objects.filter()`
    filter fields:
        - `viewer.models.SnapshotActions.snapshot` - ?snapshot=<str>
        - `viewer.models.SnapshotActions.last_update_date` - ?date=<str>

    returns: JSON
        - id: id of the snapshot action record
        - author: id of the user that created the snapshot
        - session_project: id of the related session_project (if present)
        - snapshot: id of the related snapshot
        - last_update_date: Timestamp for when the action list was generated or updated
        - actions: JSON string containing actions related to the session_project

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 1,
                    "last_update_date": "2020-10-06T15:36:00Z",
                    "actions": {
                        "save": false,
                        "show": true,
                        "action_type": 1,
                        "object_name": "",
                        "object_type": "",
                        "action_datetime": "2020-09-30T13:44:00.000Z"
                    },
                    "author": 1,
                    "session_project": 1,
                    "snapshot": 1
                }
            ]

   """
    queryset = SnapshotActions.objects.filter()
    serializer_class = SnapshotActionsSerializer

    #   Note: jsonField for Actions will need specific queries - can introduce if needed.
    filterset_fields = ('id', 'author', 'session_project', 'snapshot', 'last_update_date')

# End of Session Project


# Design sets upload
# Custom parser class for a csv file
class DSetCSVParser(BaseParser):
    """
    CSV parser class specific to design set csv spec - sets media_type for DSetUploadView to text/csv
    """
    media_type = 'text/csv'


class DSetUploadView(APIView):
    """DjangoRF view to upload a design set (PUT) from a csv file

    Methods
    -------
    allowed requests:
        - PUT: takes a csv file and uploads it as a design set (of 2D compounds)
    url:
       viewer/upload_designs
    request params:
        - file (csv file): csv file containing design set information
        - type (str): design set type (options in `viewer.models.DesignSet.set_type`)
        - description (str): short description of the design set - e.g. method of creation
    csv file columns (mandatory):
        - set_name - the name of the design set to upload the molecule to
        - smiles - smiles string for the 2D molecule
        - identifier - an identifier for the molecule
        - inspirations - inspiration molecules from fragalysis used in the design of the 2D molecule

    Returns
    -------
    HTTPResponse (JSON (string))
        A message telling the user whether the upload was successful or not

    """
    parser_class = (DSetCSVParser,)

    def put(self, request, format=None):  # pylint: disable=redefined-builtin
        """Method to handle PUT request and upload a design set
        """
        # Don't need...
        del format

        f = request.FILES['file']
        set_type = request.PUT['type']
        set_description = request.PUT['description']

        # save uploaded file to temporary storage
        name = f.name
        path = default_storage.save('tmp/' + name, ContentFile(f.read()))
        tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

        df = pd.read_csv(tmp_file)
        mandatory_cols = ['set_name', 'smiles', 'identifier', 'inspirations']
        actual_cols = df.columns
        for col in mandatory_cols:
            if col not in actual_cols:
                raise ParseError("The 4 following columns are mandatory: set_name, smiles, identifier, inspirations")

        set_names, compounds = process_design_sets(df, set_type, set_description)

        string = 'Design set(s) successfully created: '

        length = len(set_names)
        string += str(length) + '; '
        for i in range(0, length):
            string += str(i + 1) + ' - ' + set_names[i] + ') number of compounds = ' + str(len(compounds[i])) + '; '

        return HttpResponse(json.dumps(string))


class ComputedSetView(viewsets.ModelViewSet):
    """DjagnoRF view to retrieve information about and delete computed sets

    Methods
    -------
    allowed requests:
        - GET: retrieve all the sets or one based on its name
        - DELETE: delete a set based on name
    url:
        api/compound-sets
    returns: JSON
        - name: name of the computed set
        - submitted_sdf: link to the uploaded sdf file
        - spec_version: version of the upload specification used to generate the computed set sdf
        - method_url: link to url describing the methodology used to create the computed set
        - unique_name: auto-generated human-readable name for the computed set
        - target: id for the associated target
        - submitter: id for the associated submitter

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "name": "100threehop2020-07-27S8vv3vx",
                    "submitted_sdf": "http://fragalysis.diamond.ac.uk/media/code/media/compound_sets/Top_100_three_hop_2020-07-27_S8vv3vx.sdf",
                    "spec_version": 1.2,
                    "method_url": "https://github.com/Waztom/xchem-xCOS",
                    "unique_name": "WT-xCOS2-ThreeHop",
                    "upload_task_id": null,
                    "upload_status": null,
                    "upload_progress": null,
                    "upload_datetime": null,
                    "target": 1,
                    "submitter": 1,
                    "owner_user": 3
                },]

    """
    queryset = ComputedSet.objects.filter()
    serializer_class = ComputedSetSerializer
    filter_permissions = "project_id"
    filterset_fields = ('target', 'target__title')

    http_method_names = ['get', 'head', 'delete']

    def destroy(self, request, pk=None):
        """User provides the name of the ComputedSet (that's its primary key).
        We simply look it up and delete it, returning a standard 204 on success.
        """
        computed_set = get_object_or_404(ComputedSet, pk=pk)
        computed_set.delete()
        return HttpResponse(status=204)


class ComputedMoleculesView(viewsets.ReadOnlyModelViewSet):
    """ DjagnoRF view to retrieve information about computed molecules - 3D info

    Methods
    -------
    url:
        api/compound-molecules
    queryset:
        `viewer.models.ComputedMolecule.objects.filter()`
    filter fields:
        - `viewer.models.ComputedMolecule.computed_set` - ?computed_set=<int>
    returns: JSON
        - id: id of the molecule
        - sdf_info: 3D coordinates of the molecule in MDL format
        - name: a name for the molecule
        - smiles: SMILES string
        - pdb_info: link to the associated pdb file (apo)
        - compound: id for the associated 2D compound
        - computed_set: name for the associated computed set
        - computed_inspirations: list of ids for the inspirations used in the design/computation of the molecule

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 1997,
                    "sdf_info": "FRA-DIA-8640f307-1    RDKit          3D 38 42  0  0  0  0  0  0  0  0999 V2000..."
                    "name": "FRA-DIA-8640f307-1",
                    "smiles": "CC(=O)NCCc1c[nH]c2c([C@H](CN(Cc3cc(C)on3)C(=O)NC3CC3)c3nnc(C)s3)cc(F)cc12",
                    "pdb_info": "http://fragalysis.diamond.ac.uk/media/pdbs/Fragmenstein_J6Sfvrs.pdb",
                    "compound": 4030,
                    "computed_set": "100XCOS2Teo2020-07-23yuZJZFY",
                    "computed_inspirations": []
                },]


    """
    queryset = ComputedMolecule.objects.filter()
    serializer_class = ComputedMoleculeSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set',)


class NumericalScoresView(viewsets.ReadOnlyModelViewSet):
    """ DjagnoRF view to retrieve information about numerical computed molecule scores

    Methods
    -------
    url:
        api/numerical-scores
    queryset:
        `viewer.models.NumericalScoreValues.objects.filter()`
    filter fields:
        - `viewer.models.NumericalScoreValues.compound` - ?compound=<int>
        - `viewer.models.NumericalScoreValues.score` - ?score=<int>
    returns: JSON
        - id: id of the score
        - score: dict of the score info:
            - id: id of the score description
            - name: name of the score
            - description: description of the score
            - computed_set: name of the computed set that the score is associated to
        - value: numerical value of the score
        - compound: id of the associated compound object

    example output:

        .. code-block:: javascript

            "results": [
                "results": [
                    {
                        "id": 8145,
                        "score": {
                            "id": 48,
                            "name": "Score_1",
                            "description": "The score is scaled by the number of bit atoms",
                            "computed_set": "100XCOS2Teo2020-07-23yuZJZFY"
                        },
                        "value": 19.8653,
                        "compound": 1997
                    },


    """

    queryset = NumericalScoreValues.objects.filter()
    serializer_class = NumericalScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('compound', 'score')


class TextScoresView(viewsets.ReadOnlyModelViewSet):
    """ DjagnoRF view to retrieve information about text computed molecule scores

    Methods
    -------
    url:
        api/text-scores
    queryset:
        `viewer.models.TextScoreValues.objects.filter()`
    filter fields:
        - `viewer.models.TextScoreValues.compound` - ?compound=<int>
        - `viewer.models.TextScoreValues.score` - ?score=<int>
    returns: JSON
        - id: id of the score
        - score: dict of the score info:
            - id: id of the score description
            - name: name of the score
            - description: description of the score
            - computed_set: name of the computed set that the score is associated to
        - value: text value of the score
        - compound: id of the associated compound object

    example output:

        .. code-block:: javascript

            "results": [
                "results": [
                    {
                        "id": 8145,
                        "score": {
                            "id": 48,
                            "name": "Score_1",
                            "description": "Desctription",
                            "computed_set": "100XCOS2Teo2020-07-23yuZJZFY"
                        },
                        "value": "Yes",
                        "compound": 1997
                    },


    """
    queryset = TextScoreValues.objects.filter()
    serializer_class = TextScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('compound', 'score')


class CompoundScoresView(viewsets.ReadOnlyModelViewSet):
    """ DjagnoRF view to retrieve descriptions of scores for a given name or computed set

    Methods
    -------
    url:
        api/compound-scores
    queryset:
        `viewer.models.ScoreDescription.objects.filter()`
    filter fields:
        - `viewer.models.ScoreDescription.computed_set` - ?computed_set=<int>
        - `viewer.models..ScoreDescription.name` - ?name=<str>
    returns: JSON
        - id: id of the score
        - score: dict of the score info:
            - id: id of the score description
            - name: name of the score
            - description: description of the score
            - computed_set: name of the computed set that the score is associated to

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 22,
                    "name": "FeatureSteinScore",
                    "description": "FeatureStein Score",
                    "computed_set": "top26062020axX6Vqt"
                },]


    """
    queryset = ScoreDescription.objects.filter()
    serializer_class = ScoreDescriptionSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set', 'name')


class ComputedMolAndScoreView(viewsets.ReadOnlyModelViewSet):
    """ DjagnoRF view to retrieve all information about molecules from a computed set, along with all of their scores

    Methods
    -------
    url:
        api/compound-mols-scores
    queryset:
        `viewer.models.ComputedMolecule.objects.filter()`
    filter fields:
        - `viewer.models.ComputedMolecule.computed_set` - ?computed_set=<int>
    returns: JSON
        - id: id of the molecule
        - sdf_info: 3D coordinates of the molecule in MDL format
        - name: a name for the molecule
        - smiles: SMILES string
        - pdb_info: link to the associated pdb file (apo)
        - compound: id for the associated 2D compound
        - computed_set: name for the associated computed set
        - computed_inspirations: list of ids for the inspirations used in the design/computation of the molecule
        - numerical_scores: dict of numerical scores, where each key is a name, and each value is the associated score
        - text_scores: dict of text scores, where each key is a score name, and each value is the associated score

    example output:

        .. code-block:: javascript

            "results": [
                {
                    "id": 1997,
                    "sdf_info": "FRA-DIA-8640f307-1     RDKit          3D 38 42  0  0  0  0  0  0  0"
                    "name": "FRA-DIA-8640f307-1",
                    "smiles": "CC(=O)NCCc1c[nH]c2c([C@H](CN(Cc3cc(C)on3)C(=O)NC3CC3)c3nnc(C)s3)cc(F)cc12",
                    "pdb_info": "http://fragalysis.diamond.ac.uk/media/pdbs/Fragmenstein_J6Sfvrs.pdb",
                    "compound": 4030,
                    "computed_set": "100XCOS2Teo2020-07-23yuZJZFY",
                    "computed_inspirations": [],
                    "numerical_scores": {
                        "Score_1": 19.8653,
                        "N_hits": 4.0
                    },
                    "text_scores": {}
                },]


    """
    queryset = ComputedMolecule.objects.filter()
    serializer_class = ComputedMolAndScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set',)


class DiscoursePostView(viewsets.ViewSet):
    """Django view to get and post to the Discourse platform

    Methods
    -------
    allowed requests:
        - POST: Takes a post and calls a function to post the details to the Discourse platform.
        - GET: (No yet fully operational - returns posts for ?post_title=<topic title>
    url:
       api/discourse_post
    params:
        - category_name: sub category_name for Discourse post (optional if post title is given)
        - parent_category_name: Discourse parent_category name - defaults to "Fragalysis targets" (Setting)
        - category_colour: Optional - defaults to '0088CC'
        - category_text_colour: Optional defaults to 'FFFFFF'
        - post_title: title of topic or post (optional if category name is given)
        - post_content: content of the post
        - post_tags: a JSON string of tags related to the Post

    Returns JSON

    example of input (GET) on local:

        http://127.0.0.1:8080/api/discourse_post/?post_title=Mpro%20First%20Project

    examples of input (POST raw data):

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "", "post_content": "", "post_tags":""}

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "New Topic Title 1",
        "post_content": "This is the first post that creates the topic - must be greater than 20 chars",
        "post_tags" :"[\\"tag1\\",\\"tag2\\"]"}

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "New Topic Title 1",
        "post_content": "This is a second post to New Topic Title 1", "post_tags" :"[]"}

    example output (POST):

        .. code-block:: javascript

          {
                "Post url": "https://discourse.xchem-dev.diamond.ac.uk/t/78/1"
          }

    example output (GET):

        .. code-block:: javascript

            {
                "Posts": {
                    "post_stream": {
                        "posts": [
                            {
                                "id": 131,
                                "name": "user",
                                "username": "user",
                                "avatar_template": "/letter_avatar_proxy/v4/letter/u/c0e974/{size}.png",
                                "created_at": "2020-12-10T14:50:56.006Z",
                                "cooked": "<p>This is a post for session-project 0005 to test if it works without parent category</p>",
                                "post_number": 1,
                                "post_type": 1,
                                "updated_at": "2020-12-10T14:50:56.006Z",
                                "reply_count": 0,
                                "reply_to_post_number": null,
                                "quote_count": 0,
                                "incoming_link_count": 1,
                                "reads": 1,
                                "readers_count": 0,
                                "score": 5.2,
                                "yours": true,
                                "topic_id": 81,
                                "topic_slug": "api-test-session-project-000005",
                                "display_username": "user",
                                "primary_group_name": null,
                                "primary_group_flair_url": null,
                                "primary_group_flair_bg_color": null,
                                "primary_group_flair_color": null,
                                "version": 1,
                                "can_edit": true,
                                "can_delete": false,
                                "can_recover": false,
                                "can_wiki": true,
                                "read": true,
                                "user_title": null,
                                "actions_summary":
                                [
                                    {
                                        "id": 3,
                                        "can_act": true
                                    },
                                    {
                                        "id": 4,
                                        "can_act": true
                                    },
                                    {
                                        "id": 8,
                                        "can_act": true
                                    },
                                    {
                                        "id": 7,
                                        "can_act": true
                                    }
                                ],
                                "moderator": false,
                                "admin": true,
                                "staff": true,
                                "user_id": 1,
                                "hidden": false,
                                "trust_level": 1,
                                "deleted_at": null,
                                "user_deleted": false,
                                "edit_reason": null,
                                "can_view_edit_history": true,
                                "wiki": false,
                                "reviewable_id": 0,
                                "reviewable_score_count": 0,
                                "reviewable_score_pending_count": 0
                            }
                        ]
                    },
                    "id": 81
                }
            }

    """

    serializer_class = DiscoursePostWriteSerializer

    def create(self, request):
        """Method to handle POST request and call discourse to create the post
        """
        logger.info('+ DiscoursePostView.post')
        data = request.data

        logger.info('+ DiscoursePostView.post %s', json.dumps(data))
        if data['category_name'] == '':
            category_details = None
        else:
            category_details = {'category_name': data['category_name'],
                                'parent_name': data['parent_category_name'],
                                'category_colour': data['category_colour'],
                                'category_text_colour': data['category_text_colour']}

        if data['post_title'] == '':
            post_details = None
        else:
            post_details = {'title': data['post_title'],
                            'content': data['post_content'],
                            'tags': json.loads(data['post_tags'])}

        error, post_url, error_message = create_discourse_post(request.user, category_details, post_details)

        logger.info('- DiscoursePostView.post')
        if error:
            return Response({"message": error_message})
        else:
            return Response({"Post url": post_url})

    def list(self, request):
        """Method to handle GET request and call discourse to list posts for a topic
        """
        logger.info('+ DiscoursePostView.get')
        query_params = request.query_params
        logger.info('+ DiscoursePostView.get %s', json.dumps(query_params))

        discourse_api_key = settings.DISCOURSE_API_KEY

        if discourse_api_key:
            post_title = request.query_params.get('post_title', None)
            error, posts = list_discourse_posts_for_topic(post_title)
        else:
            logger.info('- DiscoursePostView.no key')
            return Response({"message": "Discourse Not Available - No API key supplied"})

        logger.info('- DiscoursePostView.get')
        if error:
            return Response({"message": "No Posts Found"})
        else:
            return Response({"Posts": posts})


def create_csv_from_dict(input_dict, title=None, filename=None):
    """ Write a CSV file containing data from an input dictionary and return a URLto the file in the media
        directory.
    """

    if not filename:
        filename = 'download'

    media_root = settings.MEDIA_ROOT
    unique_dir = str(uuid.uuid4())
    # /code/media/downloads/unique_dir
    download_path = os.path.join(media_root, 'downloads', unique_dir)
    os.makedirs(download_path, exist_ok=True)

    download_file = os.path.join(download_path, filename)

    # Remove file if it already exists
    if os.path.isfile(download_file):
        os.remove(download_file)

    with open(download_file, "w", newline='', encoding='utf-8') as csvfile:
        if title:
            csvfile.write(title)
            csvfile.write("\n")

    df = pd.DataFrame.from_dict(input_dict)
    df.to_csv(download_file, mode='a', header=True, index=False)

    return download_file


class DictToCsv(viewsets.ViewSet):
    """Django view that takes a dictionary and returns a download link to a CSV file with the data.

    Methods
    -------
    allowed requests:
        - GET: Return the CSV file given the link - note that this will remove the file on the media directory.
        - POST: Return a link to a CSV file containing the input dictionary
    url:
       api/dicttocsv
    get params:
       - file_url: url returned in the post request

       Returns: CSV file when passed url.

    post params:
       - title: string to place on the first line of the CSV file.
       - input_dict: dictionary containing CSV data to place in the CSV file

       Returns: url to be passed to GET.

    example input for get

        .. code-block::

            /api/dicttocsv/?file_url=/code/media/downloads/6bc70a04-9675-4079-924e-b0ab460cb206/download

    example input for post
    ----------------------

        .. code-block:: json

            {
            "title": "https://fragalysis.xchem.diamond.ac.uk/viewer/react/landing",
            "dict": [{
                            " compound - id0 ": " CHEMSPACE - BB: CSC012451475 ",
                            " compound - id1 ": " ",
                            " smiles ": " Cc1ccncc1C(N)C(C)(C)C ",
                            " mol ": " CC( = O)Nc1cnccc1C ",
                            " vector ": " CC1CCCCC1[101Xe]",
                            " class ": " blue ",
                            " compoundClass ": " blue ",
                            " ChemPlp ": " ",
                            " MM - GBSA Nwat = 0 ": " ",
                                    " STDEV0 ": " ",
                                    " MM - GBSA Nwat = 30 ": " ",
                                    " STDEV30 ": " ",
                                    " MM - GBSA Nwat = 60 ": " "
                                },
                                {
                                    " compound - id0 ": " ",
                                    " compound - id1 ": " ",
                                    " smiles ": " CC( = O)NCCc1c[nH]c2c(C(c3ccc(Br)s3)[NH + ]3CCN(C( = O)CCl)CC3)cccc12 ",
                                    " mol ": " ",
                                    " vector ": " ",
                                    " class ": " ",
                                    " compoundClass ": " ",
                                    " ChemPlp ": -101.073,
                                    " MM - GBSA Nwat = 0 ": -38.8862,
                                    " STDEV0 ": 5.3589001,
                                    " MM - GBSA Nwat = 30 ": -77.167603,
                                    " STDEV30 ": 5.5984998,
                                    " MM - GBSA Nwat = 60 ": -84.075401
                }]
            }

    """

    serializer_class = DictToCsvSerializer

    def list(self, request):
        """Method to handle GET request
        """
        file_url = request.GET.get('file_url')

        if file_url and os.path.isfile(file_url):
            with open(file_url, encoding='utf8') as csvfile:
                # return file and tidy up.
                response = HttpResponse(csvfile, content_type='text/csv')
                response['Content-Disposition'] = 'attachment; filename=download.csv'
                shutil.rmtree(os.path.dirname(file_url), ignore_errors=True)
                return response
        else:
            return Response("Please provide file_url parameter")

    def create(self, request):
        """Method to handle POST request
        """

        logger.info('+ DictToCsv.post')
        input_dict = request.data['dict']
        input_title = request.data['title']

        if not input_dict:
            return Response({"message": "Please enter Dictionary"})
        else:
            filename_url = create_csv_from_dict(input_dict, input_title)

        return Response({"file_url": filename_url})


# Classes Relating to Tags
class TagCategoryView(viewsets.ModelViewSet):
    """ Operational Django view to set up and retrieve information about tag categories.

    Methods
    -------
    url:
        api/tag_category
    queryset:
        `viewer.models.TagCategory.objects.filter()`
    filter fields:
        - `viewer.models.TagCategory.category` - ?category=<str>
    returns: JSON

    example output:

        .. code-block:: json

            {
                "count": 1,
                "next": null,
                "previous": null,
                "results": [
                    {
                        "id": 1,
                        "category": "sites",
                        "colour": "FFFFFF",
                        "description": "site description"
                    },
                ]
            }

    """

    queryset = TagCategory.objects.filter()
    serializer_class = TagCategorySerializer
    filterset_fields = ('id', 'category')


class MoleculeTagView(viewsets.ModelViewSet):
    """ Operational Django view to set up/retrieve information about tags relating to Molecules

    Methods
    -------
    url:
        api/molecule_tag
    queryset:
        `viewer.models.MoleculeTag.objects.filter()`
    filter fields:
        - `viewer.models.MoleculeTag.tag` - ?tag=<str>
        - `viewer.models.MoleculeTag.category` - ?category=<str>
        - `viewer.models.MoleculeTag.target` - ?target=<int>
        - `viewer.models.MoleculeTag.molecules` - ?molecules=<int>
        - `viewer.models.MoleculeTag.mol_group` - ?mol_group=<int>

    returns: JSON

    example output:

        .. code-block:: json

            {
                "id": 43,
                "tag": "A9 - XChem screen - covalent hits",
                "create_date": "2021-04-20T14:16:46.850313Z",
                "colour": null,
                "discourse_url": null,
                "help_text": null,
                "additional_info": "",
                "category": 1,
                "target": 3,
                "user": null,
                "mol_group": 5468,
                "molecules": [
                    6577,
                    6578,
                    6770,
                    6771
                ]
            }

   """

    queryset = MoleculeTag.objects.filter()
    serializer_class = MoleculeTagSerializer
    filterset_fields = ('id', 'tag', 'category', 'target', 'molecules', 'mol_group')


class SessionProjectTagView(viewsets.ModelViewSet):
    """ Operational Django view to set up/retrieve information about tags relating to Session
    Projects

    Methods
    -------
    url:
        api/session_project_tag
    queryset:
        `viewer.models.SessionProjectTag.objects.filter()`
    filter fields:
        - `viewer.models.SessionProjectTag.tag` - ?tag=<str>
        - `viewer.models.SessionProjectTag.category` - ?category=<str>
        - `viewer.models.SessionProjectTag.target` - ?target=<int>
        - `viewer.models.SessionProjectTag.session_projects` - ?session_project=<int>

    returns: JSON

    example output:

        .. code-block:: json

            {
                "count": 1,
                "next": null,
                "previous": null,
                "results": [
                    {
                        "id": 3,
                        "tag": "testtag3",
                        "create_date": "2021-04-13T16:01:55.396088Z",
                        "colour": null,
                        "discourse_url": null,
                        "help_text": null,
                        "additional_info": "",
                        "category": 1,
                        "target": 1,
                        "user": null,
                        "session_projects": [
                            2
                        ]
                    }
                ]
            }

    """

    queryset = SessionProjectTag.objects.filter()
    serializer_class = SessionProjectTagSerializer
    filterset_fields = ('id', 'tag', 'category', 'target', 'session_projects')


class TargetMoleculesView(ISpyBSafeQuerySet):
    """ Django view to retrieve all Molecules and Tag information relating
    to a Target. The idea is that a single call can return all target related
    information needed by the React front end in a single call.

    Methods
    -------
    url:
        api/target_molecules/id

    returns: JSON

    example output (fragment):

        .. code-block::

            {
                "id": 4,
                "title": "nsp13",
                "project_id": [ 1 ],
                "default_squonk_project": "project-48d33e2f-6af1-42ee-a2e9-a7acf6543a1e",
                "template_protein": "/media/pdbs/nsp13-x0280_1B_apo_zOdoDll.pdb",
                "metadata": "https://127.0.0.1:8080/media/metadata/metadata_GYuEefg.csv",
                "zip_archive": "https://127.0.0.1:8080/media/targets/nsp13.zip",
                "upload_status": "SUCCESS",
                "sequences": [
                {
                    "chain": "A",
                    "sequence": ""
                }],
                "molecules": [
                {
                    "data": {
                        "id": 7012,
                        "smiles": "CS(=O)(=O)NCCc1ccccc1",
                        "cmpd_id": 185,
                        "prot_id": 6980,
                        "protein_code": "nsp13-x0176_0A",
                        "mol_type": "PR",
                        "molecule_protein": "/media/pdbs/nsp13-x0176_0A_apo.pdb",
                        "lig_id": "LIG",
                        "chain_id": "Z",
                        "sdf_info": "<SDF Block>",
                        "x_com": null,
                        "y_com": null,
                        "z_com": null,
                        "mw": 199.07,
                        "logp": 0.78,
                        "tpsa": 46.17,
                        "ha": 13,
                        "hacc": 2,
                        "hdon": 1,
                        "rots": 4,
                        "rings": 1,
                        "velec": 72
                    },
                    "tags_set": [
                        143
                    ]
                },],
                    "tags_set": [ 78 ]
                    },
                    {
                        "data": [
                            {
                    <molecule data>
                }
                ],
                "tags_info": [
                    {
                        "data": [
                            {
                                "id": 72,
                                "tag": "A - Nucleotide Site",
                                "category_id": 16,
                                "target_id": 4,
                                "user_id": null,
                                "create_date": "2021-04-22T12:11:27.315783Z",
                                "colour": null,
                                "discourse_url": null,
                                "help_text": null,
                                "additional_info": null,
                                "mol_group_id": 5498
                            }
                        ],
                        "coords": [
                            {
                                "x_com": -9.322852168872645,
                                "y_com": 3.0154678875227723,
                                "z_com": -72.34568956027785
                            }
                        ]
                    },
                    {
                        "data": [
                            {
                                "id": 73,
                    <tag data>
                }
                ],
                "tag_categories": [
                    {
                        "id": 16,
                        "category": "Sites",
                        "colour": "00CC00",
                        "description": null
                    }
                ]
            }

    """

    queryset = Target.objects.filter()
    serializer_class = TargetMoleculesSerializer
    filter_permissions = "project_id"
    filterset_fields = ("title",)
# Classes Relating to Tags - End


class DownloadStructures(ISpyBSafeQuerySet):
    """Django view that uses a selected subset of the target data
    (proteins and booleans with suggested files) and creates a Zip file
    with the contents.

    Note that old zip files are removed after one hour.

    Methods
    -------
    allowed requests:
        - GET: Return the Zip file given the link
        - POST: Return a link to a ZIP file containing the requested data.

    url:
       api/download_structures
    get params:
       - file_url: url returned in the post request

       Returns: Zip file when passed url.

    post params:
        - target_name: Selected target
        - proteins: Comma separated list of protein codes within target e.g. "Mpro-6lu7_2C,Mpro-6m0k_0A".
                    If left blank the whole target will be scanned.
        - pdb_info: True/False - include pdb file
        - bound_info: True/False - include bound file (if available)
        - cif_info: True/False - include cif file (if available)
        - mtz_info: True/False - include mtz file (if available)
        - diff_info: True/False - include diff file (if available)
        - event_info: True/False - include event file (if available)
        - sigmaa_info: True/False - include sigmaa file (if available)
        - sdf_info: True/False - include molecule sdf file (if available)
        - single_sdf_file: True/False - Also combine molecule sdf files into single file (if available)
        - trans_matrix_info: True/False - include transformation file (if available)
        - metadata_info: True/False - include metadata csv file for whole target set
        - smiles_info: True/False - include csv file containing smiles for attached molecules
        - static_link: True/False - whether zip contents will be saved as a time dependent snapshot
        - file_url: Get link to static file (this will reconstruct the zip file if it has been cleaned up.

       Returns: url to be passed to GET.

    example input for get

        .. code-block::

            /api/download_structures/?file_url=/code/media/downloads/6bc70a04-9675-4079-924e-b0ab460cb206/download

    example input for post:

        .. code-block::

            {
                "target_name": "Mpro",
                "proteins": "Mpro-6lu7_2C,Mpro-6m0k_0A",
                "pdb_info": True,
                "bound_info": True,
                "cif_info": True,
                "mtz_info": False,
                "diff_info": False,
                "event_info": False,
                "sigmaa_info": False,
                "sdf_info": False,
                "single_sdf_file": False,
                "trans_matrix_info": False,
                "metadata_info": False,
                "smiles_info": False,
                "static_link": False,
                "file_url":""
            }

    """
    queryset = Target.objects.filter()
    serializer_class = DownloadStructuresSerializer
    filter_permissions = "project_id"
    filterset_fields = ('title','id')

    def list(self, request):
        """Method to handle GET request
        """
        file_url = request.GET.get('file_url')

        if file_url:
            link = DownloadLinks.objects.filter(file_url=file_url)
            if (link and link[0].zip_file
                    and os.path.isfile(link[0].file_url)):
                logger.info('zip_file: %s', link[0].zip_file)

                # return file and tidy up.
                file_name = os.path.basename(file_url)
                wrapper = FileWrapper(open(file_url, 'rb'))
                response = FileResponse(wrapper,
                                        content_type='application/zip')
                response[
                    'Content-Disposition'] = \
                    'attachment; filename="%s"' % file_name
                response['Content-Length'] = os.path.getsize(file_url)
                return response
            elif link:
                content = {'message': 'Zip file no longer present - '
                                      'please recreate by calling '
                                      'POST/Prepare download'}
                return Response(content, status=status.HTTP_404_NOT_FOUND)

            content = {'message': 'File_url is not found'}
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        content = {'message': 'Please provide file_url parameter from '
                              'post response'}
        return Response(content, status=status.HTTP_404_NOT_FOUND)

    def create(self, request):
        """Method to handle POST request
        """
        logger.info('+ DownloadStructures.post')

        # Clear up old existing files
        maintain_download_links()

        # Static files
        # For static files, the contents of the zip file at the time of the search
        # are stored in the zip_contents field. These are used to reconstruct the
        # zip file from the time of the request.
        if request.data['file_url']:
            # This is a static link - the contents are stored in the database
            # if required.
            file_url = request.data['file_url']
            logger.info('Given file_url "%s"', file_url)
            existing_link = DownloadLinks.objects.filter(file_url=file_url)

            if existing_link and existing_link[0].static_link:
                # If the zip file is present, return it
                # Note that don't depend 100% on the zip_file flag as the
                # file might have been deleted from the media server.
                if (existing_link[0].zip_file and
                        os.path.isfile(existing_link[0].file_url)):
                    logger.info('Download is Ready!')
                    return Response({"file_url": existing_link[0].file_url},
                                    status=status.HTTP_200_OK)
                elif os.path.isfile(existing_link[0].file_url):
                    # If the file is there but zip_file is false, then it is
                    # probably being rebuilt by a parallel process.
                    logger.info('Download is under construction')
                    content = {'message': 'Zip being rebuilt - '
                                          'please try later'}
                    return Response(content,
                                    status=status.HTTP_208_ALREADY_REPORTED)
                else:
                    # Otherwise re-create the file.
                    logger.info('Recreating download...')
                    recreate_static_file (existing_link[0], request.get_host())
                    return Response({"file_url": existing_link[0].file_url},
                                    status=status.HTTP_200_OK)

            msg = 'file_url should only be provided for static files'
            logger.warning(msg)
            content = {'message': msg}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        # Dynamic files

        if 'target_name' not in request.data:
            content = {'message': 'If no file_url, a target_name (title) must be provided'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        target_name = request.data['target_name']
        target = None
        logger.info('Given target_name "%s"', target_name)

        # Check target_name is valid
        # (it should natch the title of an existing target)
        for targ in self.queryset:
            if targ.title == target_name:
                target = targ
                break

        if not target:
            msg = f'No target found with title "{target_name}"'
            logger.warning(msg)
            content = {'message': msg}
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        if request.data['proteins']:
            # Get first part of protein code
            proteins_list = [p.strip().split(":")[0]
                             for p in request.data['proteins'].split(',')]
            logger.info('Given %s proteins', len(proteins_list))
        else:
            logger.info('No proteins supplied')
            proteins_list = []

        if len(proteins_list) > 0:
            proteins = []
            # Filter by protein codes
            for code_first_part in proteins_list:
                prot = Protein.objects.filter(code__contains=code_first_part).values()
                if prot.exists():
                    proteins.append(prot.first())
        else:
            # If no protein codes supplied then return the complete list
            proteins = Protein.objects.filter(target_id=target.id).values()
        logger.info('Collected %s proteins', len(proteins))

        if len(proteins) == 0:
            content = {'message': 'Please enter list of valid protein codes '
                                  'for' + " target: {}, proteins: {} "
                .format(target.title, proteins_list) }
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        filename_url, file_exists = check_download_links(request,
                                                         target,
                                                         proteins)
        if file_exists:
            return Response({"file_url": filename_url})
        else:
            content = {'message': 'Zip being rebuilt - please try later'}
            return Response(content,
                            status=status.HTTP_208_ALREADY_REPORTED)


# Classes Relating to Squonk Jobs
class JobFileTransferView(viewsets.ModelViewSet):
    """ Operational Django view to set up/retrieve information about tags relating to Molecules

    Methods
    -------
    url:
        api/job_file_transfer
    queryset:
        `viewer.models.JobFileTransfer.objects.filter()`
    filter fields:
        - `viewer.models.JobFileTransfer.snapshot` - ?snapshot=<int>
        - `viewer.models.JobFileTransfer.target` - ?target=<int>
        - `viewer.models.JobFileTransfer.user` - ?user=<int>
        - `viewer.models.JobFileTransfer.squonk_project` - ?squonk_project=<str>
        - `viewer.models.JobFileTransfer.transfer_status` - ?transfer_status=<str>

    returns: JSON

    example input for post:

        .. code-block::

            {
                "snapshot": 2,
                "target": 5,
                "squonk_project": "project-e1ce441e-c4d1-4ad1-9057-1a11dbdccebe",
                "proteins": "CD44MMA-x0022_0A, CD44MMA-x0017_0A"
                "compounds": "PAU-WEI-b9b69149-9"
            }


    example output for post:

        .. code-block::

            {
                "id": 2,
                "transfer_status": "PENDING",
                "transfer_task_id": "d8705b7d-c065-4038-8964-c19882333247"
            }
   """

    queryset = JobFileTransfer.objects.filter()
    filter_permissions = "target__project_id"
    filterset_fields = ('id', 'snapshot', 'target', 'user',
                     'squonk_project', 'transfer_status')

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT
        or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.JobFileTransferReadSerializer`
            - if other: `viewer.serializers.JobFileTransferWriteSerializer`
        """
        if self.request.method in ['GET']:
            # GET
            return JobFileTransferReadSerializer
        # (POST, PUT, PATCH)
        return JobFileTransferWriteSerializer

    def create(self, request):
        """Method to handle POST request
        """
        logger.info('+ JobFileTransfer.post')
        # Only authenticated users can transfer files to sqonk
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can transfer files'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the Squonk2 agent is not configured
        sq2a_rv = _SQ2A.configured()
        if not sq2a_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Collect expected API parameters....
        access_id = request.data['access']  # The access ID/legacy Project record title
        target_id = request.data['target']
        snapshot_id = request.data['snapshot']
        session_project_id = request.data['session_project']

        logger.info('+ user="%s" (id=%s)', user.username, user.id)
        logger.info('+ access_id=%s', access_id)
        logger.info('+ target_id=%s', target_id)
        logger.info('+ snapshot_id=%s', snapshot_id)
        logger.info('+ session_project_id=%s', session_project_id)

        target = Target.objects.get(id=target_id)
        assert target

        # Check the user can use this Squonk2 facility.
        # To do this we need to setup a couple of API parameter objects.
        sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                        access_id=access_id,
                                                        session_id=session_project_id,
                                                        target_id=target_id)
        sq2a_send_params: SendParams = SendParams(common=sq2a_common_params,
                                                  snapshot_id=snapshot_id)
        sq2a_rv: Squonk2AgentRv = _SQ2A.can_send(sq2a_send_params)
        if not sq2a_rv.success:
            content = {f'You cannot do this ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Check the presense of the files expected to be transferred
        error, proteins, compounds = check_file_transfer(request)
        if error:
            return Response(error['message'], status=error['status'])

        # The root (in the Squonk project) where files will be written.
        # This is "understood" by the celery task (which uses this constant).
        # e.g. 'fragalysis-files'
        transfer_root = settings.SQUONK2_MEDIA_DIRECTORY
        logger.info('+ transfer_root=%s', transfer_root)

        # Create new file transfer job
        logger.info('+ Calling ensure_project() to get the Squonk2 Project...')

        # This requires a Squonk2 Project (created by the Squonk2Agent).
        # It may be an existing project, or it might be a new project.
        common_params = CommonParams(user_id=user.id,
                                        access_id=access_id,
                                        target_id=target_id,
                                        session_id=session_project_id)
        sq2_rv = _SQ2A.ensure_project(common_params)
        if not sq2_rv.success:
            msg = f'Failed to get/create a Squonk2 Project' \
                    f' for User "{user.username}", Access ID {access_id},' \
                    f' Target ID {target_id}, and SessionProject ID {session_project_id}.' \
                    f' Got "{sq2_rv.msg}".' \
                    ' Cannot continue'
            content = {'message': msg}
            logger.error(msg)
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        # The Squonk2Project record is in the response msg
        squonk2_project_uuid = sq2_rv.msg.uuid
        squonk2_unit_name = sq2_rv.msg.unit.name
        squonk2_unit_uuid = sq2_rv.msg.unit.uuid
        logger.info('+ ensure_project() returned Project uuid=%s (unit="%s" unit_uuid=%s)',
                    squonk2_project_uuid, squonk2_unit_name, squonk2_unit_uuid)

        job_transfer = JobFileTransfer()
        job_transfer.user = request.user
        job_transfer.proteins = [p['code'] for p in proteins]
        job_transfer.compounds = [c['name'] for c in compounds]
        # We should use a foreign key,
        # but to avoid migration issues with the existing code
        # we continue to use the project UUID string field.
        job_transfer.squonk_project = squonk2_project_uuid
        job_transfer.target = Target.objects.get(id=target_id)
        job_transfer.snapshot = Snapshot.objects.get(id=snapshot_id)

        # The 'transfer target' (a sub-directory of the transfer root)
        # For example the root might be 'fragalysis-files'
        # and the target may be `CD44MMA` so the targets will be written to the
        # Squonk project at fragalysis-files/CD44MMA
        assert job_transfer.target
        assert job_transfer.target.title
        transfer_target = job_transfer.target.title
        logger.info('+ transfer_target=%s', transfer_target)

        job_transfer.transfer_status = 'PENDING'
        job_transfer.transfer_datetime = None
        job_transfer.transfer_progress = None
        job_transfer.save()

        # Celery/Redis must be running.
        # This call checks and trys to start them if they're not.
        assert check_services()

        logger.info('oidc_access_token')
        logger.info(request.session['oidc_access_token'])

        logger.info('+ Starting transfer (celery) (job_transfer.id=%s)...',
                    job_transfer.id)
        job_transfer_task = process_job_file_transfer.delay(request.session['oidc_access_token'],
                                                            job_transfer.id)

        content = {'id' : job_transfer.id,
                   'transfer_root': transfer_root,
                   'transfer_target': transfer_target,
                   'transfer_status': job_transfer.transfer_status,
                   'transfer_task_id': str(job_transfer_task)}
        return Response(content,
                        status=status.HTTP_200_OK)


class JobConfigView(viewsets.ReadOnlyModelViewSet):
    """Django view that calls Squonk to get a requested job configuration

    Methods
    -------
    allowed requests:
        - GET: Get job config

    url:
       api/job_config
    get params:
       - squonk_job: name of the squonk job requested

       Returns: job details.

    example input for get

        .. code-block::

            /api/job_config/?squonk_job_name=run_smina
    """
    def list(self, request):
        """Method to handle GET request
        """
        query_params = request.query_params
        logger.info('+ JobConfigView.get: %s', json.dumps(query_params))

        # Only authenticated users can have squonk jobs
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can access squonk jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the squonk variables are not set!
        sqa_rv = _SQ2A.configured()
        if sqa_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sqa_rv.msg}'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        job_collection = request.query_params.get('job_collection', None)
        job_name = request.query_params.get('job_name', None)
        job_version = request.query_params.get('job_version', None)
        content = get_squonk_job_config(request,
                                        job_collection=job_collection,
                                        job_name=job_name,
                                        job_version=job_version)

        return Response(content)


class JobRequestView(viewsets.ModelViewSet):
    """ Operational Django view to set up/retrieve information about tags relating to Session
    Projects

    Methods
    -------
    url:
        api/job_request
    queryset:
        `viewer.models.JobRequest.objects.filter()`
    filter fields:
        - `viewer.models.JobRequest.snapshot` - ?snapshot=<int>
        - `viewer.models.JobRequest.target` - ?target=<int>
        - `viewer.models.JobRequest.user` - ?user=<int>
        - `viewer.models.JobRequest.squonk_job_name` - ?squonk_job_name=<str>
        - `viewer.models.JobRequest.squonk_project` - ?squonk_project=<str>
        - `viewer.models.JobRequest.job_status` - ?job_status=<str>

    returns: JSON

    example input for post:

        .. code-block::

            {
                "squonk_job_name": "nop",
                "snapshot": 1,
                "target": 1,
                "squonk_project": "project-e1ce441e-c4d1-4ad1-9057-1a11dbdccebe",
                "squonk_job_spec": "{\"collection\":\"im-test\",\"job\":\"nop\",\"version\":\"1.0.0\"}"
            }

    example output for post:

        .. code-block::

            {
                "id": 1,
                "squonk_url_ext": "data-manager-ui/results/instance/instance-c26fd27a-e837-4be5-af39-582b6f329f6a"
            }

    """

    queryset = JobRequest.objects.filter()
    filter_permissions = "target__project_id"
    filterset_fields = ('id', 'snapshot', 'target', 'user', 'squonk_job_name',
                     'squonk_project', 'job_status')

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT
        or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.JobRequestReadSerializer`
            - if other: `viewer.serializers.JobRequestWriteSerializer
        """
        if self.request.method in ['GET']:
            # GET
            return JobRequestReadSerializer
        # (POST, PUT, PATCH)
        return JobRequestWriteSerializer

    def create(self, request):
        """Method to handle POST request
        """
        # Celery/Redis must be running.
        # This call checks and trys to start them if they're not.
        assert check_services()

        logger.info('+ JobRequest.post')
        # Only authenticated users can create squonk job requests.
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can run jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the Squonk2 agent is not configured
        sq2a_rv = _SQ2A.configured()
        if not sq2a_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Collect expected API parameters....
        target_id = request.data['target']
        snapshot_id = request.data['snapshot']
        session_project_id = request.data['session_project']
        access_id = request.data['access']  # The access ID/legacy Project record

        logger.info('+ user="%s" (id=%s)', user.username, user.id)
        logger.info('+ access_id=%s', access_id)
        logger.info('+ target_id=%s', target_id)
        logger.info('+ snapshot_id=%s', snapshot_id)
        logger.info('+ session_project_id=%s', session_project_id)

        # Check the user can use this Squonk2 facility.
        # To do this we need to setup a couple of API parameter objects.
        # We don't (at this point) care about the Job spec or callback URL.
        sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                        access_id=access_id,
                                                        session_id=session_project_id,
                                                        target_id=target_id)
        sq2a_run_job_params: RunJobParams = RunJobParams(common=sq2a_common_params,
                                                         job_spec=None,
                                                         callback_url=None)
        sq2a_rv: Squonk2AgentRv = _SQ2A.can_run_job(sq2a_run_job_params)
        if not sq2a_rv.success:
            content = {f'You cannot do this ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        try:
            job_id, squonk_url_ext = create_squonk_job(request)
        except ValueError as error:
            logger.info('Job Request failed: %s', error)
            content = {'error': str(error)}
            return Response(content,
                            status=status.HTTP_400_BAD_REQUEST)

        logger.info('SUCCESS (job_id=%s squonk_url_ext=%s)', job_id, squonk_url_ext)

        content = {'id': job_id, 'squonk_url_ext': squonk_url_ext}
        return Response(content,
                        status=status.HTTP_200_OK)


class JobCallBackView(viewsets.ModelViewSet):
    """ View to allow the Squonk system to update the status and job information for a
    specific job identified by a UUID.

    Methods
    -------
    allowed requests:
        - GET
        - PUT - update the status or job information fields
    url:
        api/job_callback/<job_request.code>
    queryset:
        `viewer.models.JobRequest.objects.filter()`
        'lookup_value = code'

    returns: JSON

    example input:

        .. code-block::

            {
                "job_status": "SUCCESS"
            }


    """

    queryset = JobRequest.objects.all()
    lookup_field = "code"
    http_method_names = ['get', 'head', 'put']

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a PUT

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.JobCallBackWriteSerializer`
            - if other: `viewer.serializers.JobCallBackWriteSerializer`
        """
        if self.request.method in ['GET']:
            # GET
            return JobCallBackReadSerializer
        # PUT
        return JobCallBackWriteSerializer

    def update(self, request, code=None):
        """Response to a PUT on the Job-Callback.
        We're given a 'code' which we use to lookup the corresponding JobRequest
        (there'll only be one).
        """

        jr = JobRequest.objects.get(code=code)
        logger.info('+ JobCallBackView.update(code=%s) jr=%s', code, jr)

        # request.data is rendered as a dictionary
        if not request.data:
            return HttpResponse(status=204)

        j_status = request.data['job_status']
        # Get the appropriate SQUONK_STATUS...
        status_changed = False
        for squonk_status in JobRequest.SQUONK_STATUS:
            if squonk_status[0] == j_status and jr.job_status != j_status:
                jr.job_status = squonk_status[1]
                status_changed = True
                break

        if not status_changed:
            logger.info('code=%s status=%s ignoring (no status change)',
                        code, j_status)
            return HttpResponse(status=204)

        # This is now a chance to safely set the squonk_url_ext using the instance ID
        # present in the callback (if it's not already set). The instance is
        # placed at the end of the string, and is expected to be found
        # by process_compound_set_file() which loads the output file back
        # into Fragalysis
        if not jr.squonk_url_ext:
            jr.squonk_url_ext = create_squonk_job_request_url(request.data['instance_id'])
            logger.info("Setting jr.squonk_url_ext='%s'", jr.squonk_url_ext)

        # Update the state transition time,
        # assuming UTC.
        transition_time = request.data.get('state_transition_time')
        if not transition_time:
            transition_time = str(datetime.utcnow())
            logger.warning("Callback is missing state_transition_time"
                           " (using '%s')", transition_time)
        transition_time_utc = parse(transition_time).replace(tzinfo=pytz.UTC)
        jr.job_status_datetime = transition_time_utc

        logger.info('code=%s status=%s transition_time=%s (new status)',
                    code, j_status, transition_time)

        # If the Job's start-time is not set, set it.
        if not jr.job_start_datetime:
            logger.info('Setting job START datetime (%s)', transition_time)
            jr.job_start_datetime = transition_time_utc

        # Set the Job's finish time (once) if it looks lie the Job's finished.
        # We can assume the Job's finished if the status is one of a number
        # of values...
        if not jr.job_finish_datetime and j_status in ('SUCCESS', 'FAILURE', 'REVOKED'):
            logger.info('Setting job FINISH datetime (%s)', transition_time)
            jr.job_finish_datetime = transition_time_utc

        # Save the JobRequest record before going further.
        jr.save()

        if j_status != 'SUCCESS':
            # Go no further unless SUCCESS
            return HttpResponse(status=204)

        logger.info('Job finished (SUCCESS). Can we upload the results..?')

        # SUCCESS ... automatic upload?
        #
        # Only continue if the target file is 'merged.sdf'.
        # For now there must be an '--outfile' in the job info's 'command'.
        # Here we have hard-coded the expectations because the logic to identify the
        # command's outputs is not fully understood.
        # The command is a string that we split and search.
        job_output = ''
        jr_job_info_msg = jr.squonk_job_info[1]
        command = jr_job_info_msg.get('command')
        command_parts = shlex.split(command)
        outfile_index = 0
        while outfile_index < len(command_parts)\
                and command_parts[outfile_index] != '--outfile':
            outfile_index += 1
        # Found '--command'?
        if command_parts[outfile_index] == '--outfile'\
                and outfile_index < len(command_parts) - 1:
            # Yes ... the filename is the next item in the list
            job_output = command_parts[outfile_index + 1]
        job_output_path = '/' + os.path.dirname(job_output)
        job_output_filename = os.path.basename(job_output)

        logging.info('job_output_path="%s"', job_output_path)
        logging.info('job_output_filename="%s"', job_output_filename)

        # If it's not suitably named, leave
        expected_squonk_filename = 'merged.sdf'
        if job_output_filename != expected_squonk_filename:
            # Incorrectly named file - nothing to get/upload.
            logger.info('SUCCESS but not uploading.'
                        ' Expected "%s" as job_output_filename.'
                        ' Found "%s"', expected_squonk_filename, job_output_filename)
            return HttpResponse(status=204)

        if jr.upload_status != 'PENDING':
            logger.warning('SUCCESS but ignoring.'
                           ' upload_status=%s (already uploading?)', jr.upload_status)
            return HttpResponse(status=204)

        # Change of status and SUCCESS
        # - mark the job upload as 'started'
        jr.upload_status = "STARTED"
        jr.save()

        # Initiate an upload (and removal) of files from Squonk.
        # Which requires the linking of several tasks.
        # We star the process with 'process_compound_set_job_file'
        # with the path and filename already discoverd...
        task_params = {'jr_id': jr.id,
                       'transition_time': transition_time,
                       'job_output_path': job_output_path,
                       'job_output_filename': job_output_filename}
        task_upload = (
             process_compound_set_job_file.s(task_params) |
             validate_compound_set.s() |
             process_compound_set.s() |
             erase_compound_set_job_material.s(job_request_id=jr.id)
        ).apply_async()

        logger.info('Started process_job_file_upload(%s) task_upload=%s',
                    jr.id, task_upload)

        return HttpResponse(status=204)

class JobAccessView(APIView):
    """Django view that calls Squonk to allow a user (who is able to see a Job)
    the ability to access that Job in Squonk. To be successful the user
    must have access to the corresponding Fragalysis Project. This can be called by
    the Job 'owner', who always has access.

    Methods
    -------
    allowed requests:
        - GET: Get job access

    url:
       api/job_access
    get params:
       - job_request_id: The identity of the JobRequest (the Job) the user needs access to

       Returns: A structure with 'accessible' set to True on success. On failure
                an 'error' string contains a reason

    example input for get

        .. code-block::

            /api/job_access/?job_request_id=17
    """
    def get(self, request):
        """Method to handle GET request
        """
        query_params = request.query_params
        logger.info('+ JobAccessView/GET %s', json.dumps(query_params))

        err_response = {'accessible': False}
        ok_response = {'accessible': True, 'error': ''}

        # Only authenticated users can have squonk jobs
        user = self.request.user
        if not user.is_authenticated:
            content = {'accessible': False,
                       'error': 'Only authenticated users can access Squonk Jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the squonk variables are not set!
        sqa_rv = _SQ2A.configured()
        if not sqa_rv.success:
            err_response['error'] = f'Squonk is not available ({sqa_rv.msg})'
            return Response(err_response, status=status.HTTP_403_FORBIDDEN)

        # Get the JobRequest Record form the supplied record ID
        jr_id_str = request.query_params.get('job_request_id', None)
        if not jr_id_str or not jr_id_str.isdigit():
            err_response['error'] = f'The JobRequest ID ({jr_id_str}) is not valid'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)
        jr_id = int(jr_id_str)
        if jr_id < 1:
            err_response['error'] = f'The JobRequest ID ({jr_id}) cannot be less than 1'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)

        jr_list = JobRequest.objects.filter(id=jr_id)
        if len(jr_list) == 0:
            err_response['error'] = 'The JobRequest does not exist'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)
        jr = jr_list[0]

        # JobRequest must have a Squonk Project value
        if not jr.squonk_project:
            err_response['error'] = f'The JobRequest ({jr_id}) has no Squonk Project value'
            return Response(err_response, status=status.HTTP_403_FORBIDDEN)

        # User must have access to the Job's Project.
        # If the user is not the owner of the Job, and there is a Project,
        # we then check the user has access to the gievn access ID.
        #
        # If the Job has no Project (Jobs created before this chnage will not have a Project)
        # or the user is the owner of the Job we skip this check.
        if user.id != jr.user.id:
            logger.info('+ JobAccessView/GET Checking access to JobRequest %s for "%s" (%s)',
                        jr_id, user.username, jr.squonk_project)
            if not jr.project or not jr.project.title:
                logger.warning('+ JobAccessView/GET No Fragalysis Project (or title)'
                               ' for JobRequest %s - granting access',
                               jr_id)
            else:
                # The project is the Job's access ID.
                # To check access we need this and the User's ID
                access_id = jr.project.id
                access_string = jr.project.title
                sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                                access_id=access_id,
                                                                session_id=None,
                                                                target_id=None)
                sq2a_run_job_params: RunJobParams = RunJobParams(common=sq2a_common_params,
                                                                 job_spec=None,
                                                                 callback_url=None)
                sq2a_rv: Squonk2AgentRv = _SQ2A.can_run_job(sq2a_run_job_params)
                if not sq2a_rv.success:
                    err_response['error'] = f'Access to the Job for {access_id} ({access_string}) is denied. {sq2a_rv.msg}'
                    return Response(err_response, status=status.HTTP_403_FORBIDDEN)

            # All checks passed ... grant access for this user
            sq2a_access_params: AccessParams = AccessParams(username=user.username,
                                                            project_uuid=jr.squonk_project)
            sqa_rv = _SQ2A.grant_access(sq2a_access_params)
            if not sqa_rv.success:
                err_response['error'] = f'Squonk failed to grant access ({sqa_rv.msg})'
                logger.warning('+ JobAccessView/GET error=%s', content['error'])
                return Response(err_response, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        logger.info('+ JobAccessView/GET Success for %s/"%s" on %s',
                    jr_id, user.username, jr.squonk_project)
        return Response(ok_response)
