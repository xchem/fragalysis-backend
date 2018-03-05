from network.functions import add_keys, add_empty, order_stuctures
from django.test import TestCase
import json

class NetworkUtilsTestCase(TestCase):

    def test_can_add_key(self):
        test_d = {"D_ONE":{"TYPE_ONE":{}}}
        out_d = {}
        depth = "D_ONE"
        this_type = "TYPE_ONE"
        final_d = add_keys(out_d, depth, this_type)
        self.assertEqual(final_d,test_d)

    def test_can_add_empty(self):
        test_d = {"D_ONE":{"TYPE_ONE":{"POS_ONE": {"smiles": [],"annotation":"ANNOT"}}}}
        depth = "D_ONE"
        this_type = "TYPE_ONE"
        this_list = ["POS_ONE"]
        tot_position = []
        annotation = "ANNOT"
        out_d = {"D_ONE": {"TYPE_ONE": {}}}
        final_d = add_empty(out_d, tot_position, this_list, depth, this_type, annotation)
        self.assertEqual(final_d,test_d)

    def test_order_stuctures(self):
        test_d = {"DEPTH":{"TYPEONE":{"ONE": {"smiles": "SMILES","annotation":"ADD_DEC"}}}}
        results = {"ONE_DEPTH_TYPEONE": "SMILES"}
        decoration_list = ["ONE_ME","TWO_ME","THREE_ME"]
        out_json = order_stuctures(results, decoration_list)
        self.assertEqual(out_json,json.dumps(test_d))
