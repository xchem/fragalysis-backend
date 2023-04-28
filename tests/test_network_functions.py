from django.test import SimpleTestCase, tag

from network.functions import add_keys, add_empty, order_structures


@tag("nodb")
class NetworkFunctionsTestCase(SimpleTestCase):

    def test_can_add_key(self):
        test_d = {"D_ONE": {"TYPE_ONE": {}}}
        out_d = {}
        depth = "D_ONE"
        this_type = "TYPE_ONE"
        final_d = add_keys(out_d, depth, this_type)
        self.assertEqual(final_d, test_d)

    def test_can_add_empty(self):
        test_d = {
            "D_ONE": {"TYPE_ONE": {"POS_ONE": {"smiles": [], "annotation": "ANNOT"}}}
        }
        depth = "D_ONE"
        this_type = "TYPE_ONE"
        this_list = ["POS_ONE"]
        tot_position = []
        annotation = "ANNOT"
        out_d = {"D_ONE": {"TYPE_ONE": {}}}
        final_d = add_empty(
            out_d, tot_position, this_list, depth, this_type, annotation
        )
        self.assertEqual(final_d, test_d)

    def test_order_structures(self):
        test_json = (
            '{"DEPTH": {"TYPEONE": {"ONE": {"smiles": "SMILES", "annotation": "ADD_DEC"}}}, "-1": {"ADDITION": {"O": {"smiles": [], "annotation": "ADD_MISS"}, "N": {"smiles": [], "annotation": "ADD_MISS"}, "E": {"smiles": [], "annotation": "ADD_MISS"}, "_": {"smiles": [], "annotation": "ADD_MISS"}, "M": {"smiles": [], "annotation": "ADD_MISS"}}, "DELETION": {"T": {"smiles": [], "annotation": "DEL_MISS"}, "W": {"smiles": [], "annotation": "DEL_MISS"}, "O": {"smiles": [], "annotation": "DEL_MISS"}, "_": {"smiles": [], "annotation": "DEL_MISS"}, "M": {"smiles": [], "annotation": "DEL_MISS"}, "E": {"smiles": [], "annotation": "DEL_MISS"}}, "LINKER": {"T": {"smiles": [], "annotation": "LINK_MISS"}, "H": {"smiles": [], "annotation": "LINK_MISS"}, "R": {"smiles": [], "annotation": "LINK_MISS"}, "E": {"smiles": [], "annotation": "LINK_MISS"}, "_": {"smiles": [], "annotation": "LINK_MISS"}, "M": {"smiles": [], "annotation": "LINK_MISS"}}}}'
        )
        results = {"ONE_DEPTH_TYPEONE": "SMILES"}
        decoration_list = ["ONE_ME", "TWO_ME", "THREE_ME"]
        out_json = order_structures(results, decoration_list)
        self.maxDiff = None
        self.assertEqual(out_json, test_json,)
