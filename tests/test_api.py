from rest_framework.test import APIRequestFactory
from viewer.models import Molecule,Protein,Target,Compound
from pandda.models import PanddaEvent,PanddaSite
from hypothesis.models import Vector3D, Vector, ProteinResidue, TargetResidue, InteractionPoint,Interaction
from django.contrib.auth.models import AnonymousUser, User
from django.test import TestCase, RequestFactory

from rest_framework.test import APIClient
from api.utils import draw_mol,get_token
# Test all these functions


class APIUtilsTestCase(TestCase):

    def setUp(self):
        self.factory = RequestFactory()
        self.user = User.objects.create(username="DUMMY",password="DUMMY")

    def test_can_draw(self):
        output_mol = draw_mol("C1CCCCC1")
        svg_str = '<?xml version="1.0" encoding="UTF-8"?><svg baseProfile="full" height="200px" version="1.1" width="200px" xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">\n<rect height="200" style="opacity:1.0;fill:none;stroke:none" width="200" x="0" y="0"> </rect>\n<path d="M 190.909,100 145.455,178.73" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 190.909,100 145.455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 145.455,178.73 54.5455,178.73" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 54.5455,178.73 9.09091,100" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 9.09091,100 54.5455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 54.5455,21.2704 145.455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n</svg>'
        self.assertEqual(output_mol,svg_str)
        self.assertTrue(type(output_mol)==str)
        none_output_mol = draw_mol("C1CcccC1")
        self.assertEqual(none_output_mol,"None Mol")


    def test_can_get_token(self):
        request = self.factory.get('/viewer/react/')
        request.user = self.user
        token_one = get_token(request)
        self.assertTrue(type(token_one)==unicode)
        self.assertNotEqual(token_one,"")
        request.user = AnonymousUser()
        token_two = get_token(request)
        self.assertNotEqual(token_two,"")
        self.assertTrue(type(token_two)==unicode)

class APIUrlsTestCase(TestCase):

    def setUp(self):
        self.factory = APIRequestFactory()
        self.client = APIClient()
        self.user = User.objects.create(username="DUMMY",password="DUMMY")
        self.client.login(username=self.user.username, password=self.user.password)
        self.target = Target.objects.create(id=1,title="DUMMY_TARGET")
        self.cmpd = Compound.objects.create(id=1,inchi="DUM_INCH",smiles="DUM_SMI",mol_log_p=0.1,
                                            mol_wt=0.2,tpsa=0.3,heavy_atom_count=1,heavy_atom_mol_wt=2,
                                            nhoh_count=3,no_count=4,num_h_acceptors=5,num_h_donors=6,
                                            num_het_atoms=7,num_rot_bonds=8,num_val_electrons=9,ring_count=10)
        self.protein = Protein.objects.create(id=1,code="DUMM",target_id=self.target,
                                              pdb_info="my_pdb.pdb")
        self.mol = Molecule.objects.create(id=1,smiles="DUMMY",lig_id="DUM",chain_id="C",sdf_info="DUMMY_SD",
                                rscc=0.1,occupancy=0.2,x_com=0.3,y_com=0.4,z_com=0.5,rmsd=0.6,
                                prot_id=self.protein,cmpd_id=self.cmpd)
        self.site = PanddaSite.objects.create(id=1, site_id=1, pandda_version="0.0.1-alpha",
                                              target_id=self.target,
                                              site_align_com_x=0.1, site_align_com_y=0.2, site_align_com_z=0.3)
        self.event = PanddaEvent.objects.create(id=1,xtal="DUMMY-x001",event=1,pandda_site=self.site,
                                          target_id=self.target,pdb_info="my_pdb.pdb",map_info="my_map.map",
                                                small_map_info="my_map_small.map",
                                          lig_id="LIG",event_com_x=0.1,event_com_y=0.2,event_com_z=0.3)
        self.vector = Vector.objects.create(id=1,cmpd_id=self.cmpd,smiles="DUMMY",type="DE")
        self.vector3d = Vector3D.objects.create(id=1,mol_id=self.mol,vector_id=self.vector,number=1)
        self.target_res = TargetResidue.objects.create(id=1,target_id=self.target,res_name="DED",res_num=1,chain_id="A")
        self.prot_res = ProteinResidue.objects.create(id=1,prot_id=self.protein,targ_res_id=self.target_res)
        self.interaction_point = InteractionPoint.objects.create(id=1,mol_id=self.mol,prot_res_id=self.prot_res,
                                                                 protein_atom_name="A",molecule_atom_name="B")
        self.interaction = Interaction.objects.create(id=1,interaction_version="DE",interaction_type="UK",
                                                      interaction_point=self.interaction_point)





    def testV0_1API(self):
        """
        Untested but check get API works the way we want
        :return:
        """
        url_base = "/v0.1"
        urls = ['molecules','compounds','targets','proteins','sites','events',
                'vectors','vector3ds','proteinres','targetres',
                'interactionpoints','interactions']

        # 'scorechoice','molchoice','protchoice','cmpdchoice',
        # 'viewscene',
        # 'molgroup']
        response_data = [{'id':1,'smiles':"DUMMY", 'cmpd_id':1, 'prot_id':1, 'lig_id':"DUM",
                         'chain_id':"C", 'sdf_info':"DUMMY_SD", 'x_com':0.3, 'y_com':0.4, 'z_com':0.5},
                         {'id':1, 'inchi':"DUM_INCH", 'smiles':"DUM_SMI", 'mol_log_p':0.1, 'mol_wt':0.2,
                          'num_h_acceptors':5, 'num_h_donors':6},
                         {'id':1, 'title':"DUMMY_TARGET", 'project_id': [], 'protein_set': [1]},
                         {'id':1, 'code':"DUMM", 'target_id':1, 'pdb_info':"http://testserver/media/my_pdb.pdb",
                          'mtz_info':None, 'map_info':None, 'cif_info':None},
                         {'id': 1, 'site_id': 1, 'pandda_version': "0.0.1-alpha",
                          'target_id': 1, 'site_align_com_x': None,'pandda_run':"STANDARD", 'site_align_com_x':0.1,
                          'site_align_com_y': 0.2, 'site_align_com_z':0.3, 'site_native_com_x': None,
                          'site_native_com_y': None, 'site_native_com_z': None},
                         {'id': 1, 'xtal': "DUMMY-x001", 'event': 1, 'pandda_site': 1,
                          'target_id': 1,'pdb_info': "http://testserver/media/my_pdb.pdb", 'mtz_info': None,
                          'map_info': "http://testserver/media/my_map.map",'small_map_info': "http://testserver/media/my_map_small.map",
                          'lig_id': "LIG", 'event_com_x': 0.1, 'event_com_y': 0.2, 'event_com_z': 0.3,
                          'lig_com_x': None, 'lig_com_y': None, 'lig_com_z': None,
                          'event_dist_from_site_centroid': None, 'lig_dist_from_site_centroid': None},
                         {'id':1, 'cmpd_id':1,'smiles':"DUMMY",'type':"DE"},
                         {'id':1,'mol_id':1,'vector_id':1,'number':1,'start_x':None,'start_y':None,'start_z':None,
                          'end_x': None,'end_y': None,'end_z': None},
                         {'id': 1, 'prot_id': 1, 'targ_res_id': 1},
                         {'id':1,'target_id':1,'res_name':"DED",'res_num':1,'chain_id':"A"},
                         {'id':1,'mol_id':1,'prot_res_id':1,'protein_atom_name':"A","molecule_atom_name":"B"},
                         {'id': 1, 'interaction_version':"DE", 'interaction_type':"UK",'interaction_point': 1,
                          'distance':None,'score':None,'prot_smarts': None,'mol_smarts':None}
                         ]
        # Currently empty
        post_data = [{} for x in response_data]
        post_resp = [{u'detail': u'Method "POST" not allowed.'} for x in response_data]
        for i, url in enumerate(urls):
            # GET basic request
            response = self.client.get(url_base+"/"+url+"/1/")
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.data, response_data[i])
            # POST shouldn't work
            response = self.client.post(url_base+"/"+url+"/", post_data[i])
            self.assertEqual(response.status_code, 405)
            self.assertEqual(response.data, post_resp[i])

