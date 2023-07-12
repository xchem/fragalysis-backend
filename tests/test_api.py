import json
import zipfile
import shutil

from deepdiff import DeepDiff
from django.contrib.auth.models import AnonymousUser, User
from django.test import RequestFactory
from django.core.management.color import no_style
from django.db import connection

from rest_framework.test import APIClient, APITestCase
from rest_framework.test import APIRequestFactory

from api.utils import draw_mol, get_token
from hotspots.models import HotspotMap
from hypothesis.models import (
    Vector3D,
    Vector,
    ProteinResidue,
    TargetResidue,
    InteractionPoint,
    Interaction,
)
from viewer.models import (
    Molecule,
    Protein,
    Target,
    Compound,
    Project,
    SessionProject,
    TagCategory,
    MoleculeTag,
    SessionProjectTag
)

# Target upload functions
from viewer.target_set_upload import (
    validate_target,
    process_target
)

# Compound set upload functions
from viewer.tasks import (
    validate_compound_set,
    process_compound_set
)

# Test all these functions


class APIUtilsTestCase(APITestCase):

    def setUp(self):
        self.factory = RequestFactory()
        self.user = User.objects.create(username="DUMMY", password="DUMMY")
        self.url_base = "/api"

    # def test_can_draw(self):
    #     output_mol = draw_mol("C1CCCCC1")
    #     svg_str = "<?xml version='1.0' encoding='iso-8859-1'?>\n<svg version='1.1' baseProfile='full'\n              xmlns='http://www.w3.org/2000/svg'\n       xmlns:rdkit='http://www.rdkit.org/xml'\n                      xmlns:xlink='http://www.w3.org/1999/xlink'\n                  xml:space='preserve'\nwidth='200px' height='200px' >\n<path d='M 190.909,100 145.455,178.73' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 190.909,100 145.455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 145.455,178.73 54.5455,178.73' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 54.5455,178.73 9.09091,100' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 9.09091,100 54.5455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 54.5455,21.2704 145.455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 190.909,100 145.455,178.73' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 190.909,100 145.455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 145.455,178.73 54.5455,178.73' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 54.5455,178.73 9.09091,100' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 9.09091,100 54.5455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path d='M 54.5455,21.2704 145.455,21.2704' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n</svg>\n"
    #     self.assertEqual(output_mol, svg_str)
    #     self.assertTrue(type(output_mol) == str)
    #     none_output_mol = draw_mol("C1CcccC1")
    #     self.assertEqual(none_output_mol, "None Mol")

    def test_can_get_token(self):
        request = self.factory.get("/viewer/react/")
        request.user = self.user
        token_one = get_token(request)
        self.assertTrue(type(token_one) == str)
        self.assertNotEqual(token_one, "")
        request.user = AnonymousUser()
        token_two = get_token(request)
        self.assertNotEqual(token_two, "")
        self.assertTrue(type(token_two) == str)


    def test_generate_csv(self):
        #Check api/dicttocsv/ is available.
        response_data_get = "Please provide file_url parameter"

        response = self.client.get(self.url_base + "/dicttocsv/")
        response.user = self.user
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.data, response_data_get)

        response = self.client.get(self.url_base + "/dicttocsv/?file_url=blah")
        response.user = self.user
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.data, response_data_get)

        # POST tests to be added
        # response = self.client.post(self.url_base + "/" + url + "/", post_data[i])
        # self.assertEqual(response.status_code, 405)
        # self.assertEqual(response.data, post_resp[i])


class APIUrlsTestCase(APITestCase):

    def setUp(self):
        self.maxDiff = None
        self.factory = APIRequestFactory()
        self.client = APIClient()
        self.user = User.objects.create(username="DUMMY", password="DUMMY")
        self.user_two = User.objects.create(username="SECURE", password="SECURE")
        self.project = Project.objects.create(id=1, title="lb00000")
        self.project_secure = Project.objects.create(id=2, title="SECURE PROJECT")
        self.project_secure.user_id.add(self.user_two)
        self.project_secure.save()
        self.client.login(username=self.user.username, password=self.user.password)
        self.target = Target.objects.create(id=1, title="DUMMY_TARGET")
        self.target.project_id.add(self.project)
        self.target_two = Target.objects.create(id=2, title="SECRET_TARGET")
        self.target_two.project_id.add(self.project_secure)
        self.target_two.save()
        self.target.save()
        self.cmpd = Compound.objects.create(
            id=1,
            inchi="DUM_INCH",
            long_inchi = "DUM_L_INCHI",
            smiles="DUM_SMI",
            current_identifier='DUMM',
            all_identifiers='DUMM',
            mol_log_p=0.1,
            mol_wt=0.2,
            tpsa=0.3,
            heavy_atom_count=1,
            heavy_atom_mol_wt=2,
            nhoh_count=3,
            no_count=4,
            num_h_acceptors=5,
            num_h_donors=6,
            num_het_atoms=7,
            num_rot_bonds=8,
            num_val_electrons=9,
            ring_count=10,
        )
        self.secret_cmpd = Compound.objects.create(
            id=2,
            inchi="SEC_INCH",
            long_inchi="SEC_L_INCHI",
            smiles="SEC_SMI",
            current_identifier='SEC_DUMM',
            all_identifiers='SEC_DUMM',
            mol_log_p=0.1,
            mol_wt=0.2,
            tpsa=0.3,
            heavy_atom_count=1,
            heavy_atom_mol_wt=2,
            nhoh_count=3,
            no_count=4,
            num_h_acceptors=5,
            num_h_donors=6,
            num_het_atoms=7,
            num_rot_bonds=8,
            num_val_electrons=9,
            ring_count=10,
        )

        self.cmpd.project_id.add(self.project)
        self.cmpd.save()
        self.secret_cmpd.project_id.add(self.project_secure)
        self.secret_cmpd.save()

        self.protein = Protein.objects.create(
            id=1, code="DUMM", target_id=self.target, pdb_info="my_pdb.pdb"
        )
        self.secret_protein = Protein.objects.create(
            id=2, code="SECC", target_id=self.target_two, pdb_info="secret_pdb.pdb"
        )
        self.mol = Molecule.objects.create(
            id=1,
            smiles="DUMMY",
            lig_id="DUM",
            chain_id="C",
            sdf_info="DUMMY_SD",
            rscc=0.1,
            occupancy=0.2,
            x_com=0.3,
            y_com=0.4,
            z_com=0.5,
            rmsd=0.6,
            prot_id=self.protein,
            cmpd_id=self.cmpd,
        )
        self.secret_mol = Molecule.objects.create(
            id=2,
            smiles="SECRET",
            lig_id="SEC",
            chain_id="C",
            sdf_info="SECRET_SD",
            rscc=0.1,
            occupancy=0.2,
            x_com=0.3,
            y_com=0.4,
            z_com=0.5,
            rmsd=0.6,
            prot_id=self.secret_protein,
            cmpd_id=self.secret_cmpd,
        )
        self.vector = Vector.objects.create(
            id=1, cmpd_id=self.cmpd, smiles="DUMMY", type="DE"
        )
        self.vector3d = Vector3D.objects.create(
            id=1, mol_id=self.mol, vector_id=self.vector, number=1
        )
        self.target_res = TargetResidue.objects.create(
            id=1, target_id=self.target, res_name="DED", res_num=1, chain_id="A"
        )
        self.prot_res = ProteinResidue.objects.create(
            id=1, prot_id=self.protein, targ_res_id=self.target_res
        )
        self.interaction_point = InteractionPoint.objects.create(
            id=1,
            mol_id=self.mol,
            prot_res_id=self.prot_res,
            protein_atom_name="A",
            molecule_atom_name="B",
        )
        self.interaction = Interaction.objects.create(
            id=1,
            interaction_version="DE",
            interaction_type="UK",
            interaction_point=self.interaction_point,
        )
        self.hotspotmap = HotspotMap.objects.create(
            id=1,
            map_type="AC",
            target_id=self.target,
            prot_id=self.protein,
            map_info="my_hotspot.map",
        )

        # Tags tests

        # TagCategory - Should have been created by the migrations.
        self.tagcategory = TagCategory.objects.get(id=1)

        # MoleculeTag
        self.moltag = MoleculeTag.objects.create(
            id = 1,
            tag = "A9 - XChem screen - covalent hits",
            create_date = "2021-04-20T14:16:46.850313Z",
            colour = "FFFFFF",
            discourse_url = "www.discoursesite.com/t/1234",
            help_text = "Some help text to display as a tooltip",
            additional_info = "{'key', 'value'}",
            category = self.tagcategory,
            target = self.target,
        )
        self.moltag.molecules.add(self.mol)

        # SessionProject created for SessionProjectTag
        self.sp = SessionProject.objects.create(
            id = 1,
            title = 'test session project',
            target = self.target
        )

        # SessionProjectTag
        self.sptag = SessionProjectTag.objects.create(
            id = 1,
            tag = "Session Project Tag",
            create_date = "2021-04-20T14:16:46.850313Z",
            colour = "FFFFFF",
            discourse_url = "www.discoursesite.com/t/1234",
            help_text = "Some help text to display as a tooltip",
            additional_info = "{'key', 'value'}",
            category = self.tagcategory,
            target = self.target,
        )
        self.sptag.session_projects.add(self.sp)

        self.url_base = "/api"

        self.get_types = ["targets"] #, "molecules"]

        self.secret_target_data = {
            "targets": {
                "count": 2,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 2,
                        "title": "SECRET_TARGET",
                        "project_id": [2],
                        "protein_set": [2],
                        "default_squonk_project": None,
                        "template_protein": "/media/secret_pdb.pdb",
                        "metadata": None,
                        "upload_status": None,
                        "zip_archive": None,
                        "sequences": [{'chain': '', 'sequence': ''}]
                    },
                    {
                        "id": 1,
                        "title": "DUMMY_TARGET",
                        "project_id": [1],
                        "protein_set": [1],
                        "default_squonk_project": None,
                        "template_protein": "/media/my_pdb.pdb",
                        "metadata": None,
                        "upload_status": None,
                        "zip_archive": None,
                        "sequences": [{'chain': '', 'sequence': ''}]
                    },
                ],
            },
            "molecules": {
                "count": 2,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "smiles": "DUMMY",
                        "cmpd_id": 1,
                        "prot_id": 1,
                        "protein_code": "DUMM",
                        "lig_id": "DUM",
                        "mol_type": "PR",
                        "molecule_protein": "/media/my_pdb.pdb",
                        "chain_id": "C",
                        "sdf_info": "DUMMY_SD",
                        "x_com": 0.3,
                        "y_com": 0.4,
                        "z_com": 0.5,
                        "mw": 0.2,
                        "logp": 0.1,
                        "tpsa": 0.3,
                        "ha": 1,
                        "hacc": 5,
                        "hdon": 6,
                        "rots": 8,
                        "rings": 10,
                        "velec": 9
                    },
                    {
                        "id": 2,
                        "smiles": "SECRET",
                        "cmpd_id": 2,
                        "prot_id": 2,
                        "protein_code": "SECC",
                        "lig_id": "SEC",
                        "mol_type": "PR",
                        "molecule_protein": "/media/secret_pdb.pdb",
                        "chain_id": "C",
                        "sdf_info": "SECRET_SD",
                        "x_com": 0.3,
                        "y_com": 0.4,
                        "z_com": 0.5,
                        "mw": 0.2,
                        "logp": 0.1,
                        "tpsa": 0.3,
                        "ha": 1,
                        "hacc": 5,
                        "hdon": 6,
                        "rots": 8,
                        "rings": 10,
                        "velec": 9
                    },
                ],
            },
            "compounds": {
                "count": 2,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "inchi": "DUM_INCH",
                        "long_inchi": "DUM_L_INCHI",
                        "smiles": "DUM_SMI",
                        "current_identifier": "DUMM",
                        'all_identifiers': "DUMM",
                        "mol_log_p": 0.1,
                        "mol_wt": 0.2,
                        "tpsa": 0.3,
                        "heavy_atom_count": 1,
                        "heavy_atom_mol_wt": 2,
                        "nhoh_count": 3,
                        "no_count": 4,
                        "num_h_acceptors": 5,
                        "num_h_donors": 6,
                        "num_het_atoms": 7,
                        'num_rot_bonds': 8,
                        'num_val_electrons': 9,
                        "ring_count": 10,
                    },
                    {
                        "id": 2,
                        "inchi": "SEC_INCH",
                        "long_inchi": "SEC_L_INCHI",
                        "smiles": "SEC_SMI",
                        "current_identifier": "SEC_DUMM",
                        'all_identifiers': "SEC_DUMM",
                        "mol_log_p": 0.1,
                        "mol_wt": 0.2,
                        "tpsa": 0.3,
                        "heavy_atom_count": 1,
                        "heavy_atom_mol_wt": 2,
                        "nhoh_count": 3,
                        "no_count": 4,
                        "num_h_acceptors": 5,
                        "num_h_donors": 6,
                        "num_het_atoms": 7,
                        'num_rot_bonds': 8,
                        'num_val_electrons': 9,
                        "ring_count": 10,
                    },
                ],
            },
            "proteins": {
                "count": 2,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "code": "DUMM",
                        "target_id": 1,
                        "apo_holo": None,
                        "prot_type": "AP",
                        "pdb_info": "http://testserver/media/my_pdb.pdb",
                        "bound_info": None,
                        "mtz_info": None,
                        "map_info": None,
                        "cif_info": None,
                        "sigmaa_info": None,
                        "diff_info": None,
                        "event_info": None,
                        "aligned": None,
                        "has_eds": None,
                        "aligned_to": None
                    },
                    {
                        "id": 2,
                        "code": "SECC",
                        "target_id": 1,
                        "apo_holo": None,
                        "prot_type": "AP",
                        "pdb_info": "http://testserver/media/secret_pdb.pdb",
                        "mtz_info": None,
                        "map_info": None,
                        "cif_info": None,
                        "sigmaa_info": None,
                        "diff_info": None,
                        "event_info": None,
                        "aligned": None,
                        "has_eds": None,
                        "aligned_to": None
                    },
                ],
            },
        }
        self.not_secret_target_data = {
            "targets": {
                "count": 1,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "title": "DUMMY_TARGET",
                        "project_id": [1],
                        "protein_set": [1],
                        'default_squonk_project': None,
                        "template_protein": "/media/my_pdb.pdb",
                        "metadata": None,
                        'upload_status': None,
                        "zip_archive": None,
                        "sequences": [{'chain': '', 'sequence': ''}]
                    }
                ],
            },
            "molecules": {
                "count": 1,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "smiles": "DUMMY",
                        "cmpd_id": 1,
                        "prot_id": 1,
                        "protein_code": "DUMM",
                        "lig_id": "DUM",
                        "mol_type": "PR",
                        "molecule_protein": "/media/my_pdb.pdb",
                        "chain_id": "C",
                        "sdf_info": "DUMMY_SD",
                        "x_com": 0.3,
                        "y_com": 0.4,
                        "z_com": 0.5,
                        "mw": 0.2,
                        "logp": 0.1,
                        "tpsa": 0.3,
                        "ha": 1,
                        "hacc": 5,
                        "hdon": 6,
                        "rots": 8,
                        "rings": 10,
                        "velec": 9
                    }
                ],
            },
            "compounds": {
                "count": 1,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "inchi": "DUM_INCH",
                        "long_inchi": "DUM_L_INCHI",
                        "smiles": "DUM_SMI",
                        "current_identifier": "DUMM",
                        'all_identifiers': "DUMM",
                        "mol_log_p": 0.1,
                        "mol_wt": 0.2,
                        "tpsa": 0.3,
                        "heavy_atom_count": 1,
                        "heavy_atom_mol_wt": 2,
                        "nhoh_count": 3,
                        "no_count": 4,
                        "num_h_acceptors": 5,
                        "num_h_donors": 6,
                        "num_het_atoms": 7,
                        'num_rot_bonds': 8,
                        'num_val_electrons': 9,
                        "ring_count": 10,
                    }
                ],
            },
            "proteins": {
                "count": 1,
                "next": None,
                "previous": None,
                "results": [
                    {
                        "id": 1,
                        "code": "DUMM",
                        "target_id": 1,
                        "apo_holo": None,
                        "prot_type": "AP",
                        "pdb_info": "http://testserver/media/my_pdb.pdb",
                        "bound_info": None,
                        "mtz_info": None,
                        "map_info": None,
                        "cif_info": None,
                        "sigmaa_info": None,
                        "diff_info": None,
                        "event_info": None,
                        "aligned": None,
                        "has_eds": None,
                        "aligned_to": None
                    }
                ],
            },
        }

    def test_API(self):
        """
        Untested but check get API works the way we want
        :return:
        """
        urls = [
            "molecules",
            "compounds",
            "targets",
            "proteins",
            "vectors",
            "vector3ds",
            "proteinres",
            "targetres",
            "interactionpoints",
            "interactions",
            "hotspots",
        ]
        response_data = [
            {
                "id": 1,
                "smiles": "DUMMY",
                "cmpd_id": 1,
                "prot_id": 1,
                "protein_code": "DUMM",
                "lig_id": "DUM",
                "mol_type": "PR",
                "molecule_protein": "/media/my_pdb.pdb",
                "chain_id": "C",
                "sdf_info": "DUMMY_SD",
                "x_com": 0.3,
                "y_com": 0.4,
                "z_com": 0.5,
                "mw": 0.2,
                "logp": 0.1,
                "tpsa": 0.3,
                "ha": 1,
                "hacc": 5,
                "hdon": 6,
                "rots": 8,
                "rings": 10,
                "velec": 9
            },
            {
                "id": 1,
                "inchi": "DUM_INCH",
                "long_inchi": "DUM_L_INCHI",
                "smiles": "DUM_SMI",
                "current_identifier": "DUMM",
                'all_identifiers': "DUMM",
                "mol_log_p": 0.1,
                "mol_wt": 0.2,
                "tpsa": 0.3,
                "heavy_atom_count": 1,
                "heavy_atom_mol_wt": 2,
                "nhoh_count": 3,
                "no_count": 4,
                "num_h_acceptors": 5,
                "num_h_donors": 6,
                "num_het_atoms": 7,
                'num_rot_bonds': 8,
                'num_val_electrons': 9,
                "ring_count": 10,
            },
            {
                "id": 1,
                "title": "DUMMY_TARGET",
                "project_id": [1],
                "protein_set": [1],
                "default_squonk_project": None,
                "template_protein": "/media/my_pdb.pdb",
                "metadata": None,
                "upload_status": None,
                "zip_archive": None,
                "sequences": [{'chain': '', 'sequence': ''}]
            },
            {
                "id": 1,
                "code": "DUMM",
                "target_id": 1,
                "apo_holo": None,
                "prot_type": "AP",
                "pdb_info": "http://testserver/media/my_pdb.pdb",
                "bound_info": None,
                "mtz_info": None,
                "map_info": None,
                "cif_info": None,
                "sigmaa_info": None,
                "diff_info": None,
                "event_info": None,
                "trans_matrix_info": None,
                "pdb_header_info": None,
                "apo_desolve_info": None,
                "aligned": None,
                "has_eds": None,
                "aligned_to": None,
                "experiment": None

            },
            {"id": 1, "cmpd_id": 1, "smiles": "DUMMY", "type": "DE"},
            {
                "id": 1,
                "mol_id": 1,
                "vector_id": 1,
                "number": 1,
                "start_x": None,
                "start_y": None,
                "start_z": None,
                "end_x": None,
                "end_y": None,
                "end_z": None,
            },
            {"id": 1, "prot_id": 1, "targ_res_id": 1},
            {"id": 1, "target_id": 1, "res_name": "DED", "res_num": 1, "chain_id": "A"},
            {
                "id": 1,
                "mol_id": 1,
                "prot_res_id": 1,
                "protein_atom_name": "A",
                "molecule_atom_name": "B",
            },
            {
                "id": 1,
                "interaction_version": "DE",
                "interaction_type": "UK",
                "interaction_point": 1,
                "distance": None,
                "score": None,
                "prot_smarts": None,
                "mol_smarts": None,
            },
            {
                "id": 1,
                "map_type": "AC",
                "target_id": 1,
                "prot_id": 1,
                "map_info": "http://testserver/media/my_hotspot.map",
                "compressed_map_info": None,
            },
        ]
        self.client.login(username=self.user.username, password=self.user.password)
        # Currently empty
        post_data = [{} for x in response_data]
        post_resp = [{"detail": 'Method "POST" not allowed.'} for x in response_data]
        for i, url in enumerate(urls):
            # GET basic request
            response = self.client.get(self.url_base + "/" + url + "/1/")
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.data, response_data[i])
            # POST shouldn't work
            response = self.client.post(self.url_base + "/" + url + "/", post_data[i])
            self.assertEqual(response.status_code, 405)
            self.assertEqual(response.data, post_resp[i])

    def do_full_scan(self, user, test_data_set):
        for get_type in self.get_types:
            if user:
                self.client.force_authenticate(user)
            response = self.client.get(self.url_base + "/" + get_type + "/")
            self.assertEqual(response.status_code, 200)
            a = json.loads(json.dumps(response.json()))
            b = json.loads(json.dumps(test_data_set[get_type]))
            self.assertFalse(DeepDiff(a, b, ignore_order=True))

    def test_secure(self):
        # Test the login user  can access secure data
        self.do_full_scan(self.user_two, self.secret_target_data)

    def test_insecure(self):
        self.do_full_scan(self.user, self.not_secret_target_data)

    def test_not_logged_in(self):
        self.do_full_scan(None, self.not_secret_target_data)

    def test_tags(self):
        """
        Check basic tag functionality works.
        The data is created in the setUp(self) method
        This could be expanded to include all the session project/snapshot stuff
        when budget/opportunity allows.
        :return:
        """
        urls = [
            "tag_category",
            "molecule_tag",
            "session_project_tag"
        ]
        response_data = [
            {
                "id": 1,
                "category": "Sites",
                "colour": "00CC00",
                "description": None
            },
            {
                "id": 1,
                "tag": "A9 - XChem screen - covalent hits",
                "user": None,
                "mol_group": None,
                "create_date": "2021-04-20T14:16:46.850313Z",
                "colour": "FFFFFF",
                "discourse_url": "www.discoursesite.com/t/1234",
                "help_text": "Some help text to display as a tooltip",
                "additional_info": "{'key', 'value'}",
                "category": 1,
                "target": 1,
                "molecules": [1]
            },
            {
                "id": 1,
                "tag": "Session Project Tag",
                "user": None,
                "create_date": "2021-04-20T14:16:46.850313Z",
                "colour": "FFFFFF",
                "discourse_url": "www.discoursesite.com/t/1234",
                "help_text": "Some help text to display as a tooltip",
                "additional_info": "{'key', 'value'}",
                "category": 1,
                "target": 1,
                "session_projects": [1]
            },
        ]
        self.client.login(username=self.user.username, password=self.user.password)
        for i, url in enumerate(urls):
            # GET same data
            response = self.client.get(self.url_base + "/" + url + "/1/")
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.data, response_data[i])

    def test_validate_target(self):

        target_zip = '/code/tests/test_data/TESTTARGET.zip'
        new_data_folder = '/code/tests/output/new_data'
        target = 'TESTTARGET'
        proposal = 'open'

        # This will create the target folder in the tmp/ location.
        with zipfile.ZipFile(target_zip, 'r') as zip_ref:
            zip_ref.extractall(new_data_folder)

        validated, validate_dict = validate_target(new_data_folder,
                                                   target,
                                                   proposal)
        self.assertEqual(validated, True)
        self.assertEqual(validate_dict, {'Location': [], 'Error': [], 'Line number': []})
        # Tidy up data if not validated
        if not validated:
            shutil.rmtree(new_data_folder)


    def test_process_target(self):

        # Reset the autoincrement keys from the tests above so the target
        # can be loaded.
        models = [Target, Project, Molecule, Protein,
                  Compound, MoleculeTag, Vector, Vector3D]
        sequence_sql = connection.ops.sequence_reset_sql(no_style(), models)

        with connection.cursor() as cursor:
            for sql in sequence_sql:
                cursor.execute(sql)

        target_zip = '/code/tests/test_data/TESTTARGET.zip'
        new_data_folder = '/code/tests/output/new_data'
        target = 'TESTTARGET'
        proposal = 'open'

        # This will create the target folder in the tmp/ location.
        with zipfile.ZipFile(target_zip, 'r') as zip_ref:
            zip_ref.extractall(new_data_folder)

        mols_loaded, mols_processed = process_target(new_data_folder,
                                                     target,
                                                     proposal)

        # Improve checks as we understand more how it works.
        self.assertEqual(mols_loaded, 8)
        self.assertEqual(mols_processed, 7)
        target =  Target.objects.filter(title='TESTTARGET').values()
        self.assertEqual(len(target), 1)
        # Check finished successfully.
        self.assertEqual(target[0]['upload_status'], 'SUCCESS')
        proteins =  Protein.objects.filter(target_id__title='TESTTARGET').values()
        self.assertEqual(len(proteins), 7)
        # Selection of files are currently saved.
        for protein in proteins:
            self.assertNotEqual(protein['pdb_info'], '')
            self.assertNotEqual(protein['bound_info'], '')
            self.assertNotEqual(protein['trans_matrix_info'], '')
            self.assertNotEqual(protein['apo_desolve_info'], '')
            self.assertEqual(protein['pdb_header_info'], '')

        # Matches the number of proteins
        molecules =  Molecule.objects.filter(prot_id__target_id__title='TESTTARGET').values()
        self.assertEqual(len(molecules), 7)

        # Matches the number of sites in Metadata.csv
        tags = MoleculeTag.objects.filter(target__title='TESTTARGET').values()
        self.assertEqual(len(tags), 3)

        # Tidy up data
        shutil.rmtree(new_data_folder)


    # def test_computed_set(self):
    #     # NOTE THIS IS COMMENTED OUT BECAUSE compund-set_test.sdf DOES NOT YET HAVE CORRECT
    #     # REF_MOLS. I TOOK IT FROM AN MPRO BASED SOURCE FILE (RATHER THAN THE CD44-BASED STUFF
    #     # IN THE TESTTARGET) AND SO SOME CHANGES NEED TO BE MADE). IT MAY ALSO BE THAT
    #     # TESTTARGET.zip HAS TO BE MODIFIED TO HAVE PROTEIN CODE CONSISTENT WITH THE TARGET NAME.
    #
    #     sdf_file = \
    #         '/code/tests/test_data/compund-set_test.sdf'
    #     target = 'TESTTARGET'
    #
    #     # Check validate step
    #     validate_output = validate_compound_set(self.user.id, sdf_file, target=target)
    #
    #     # Check if SDF validated
    #     print(validate_output)
    #     self.assertEqual(validate_output[3], True)
    #     self.assertEqual(validate_output[0], 'validate')
    #     self.assertEqual(validate_output[1], 'cset')
    #
    #     # Check process step -
    #     process_output = process_compound_set(validate_output)
    #     print(process_output)
