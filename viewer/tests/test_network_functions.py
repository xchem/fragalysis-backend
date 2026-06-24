"""Tests for the network graph ordering (issue #984, phase 3).

``network.functions.order_structures`` re-shapes a flat ``{key: smiles}`` graph
result into a nested depth/type/position dict, decorating each position as an
ADD/DEL/LINK based on three decoration lists, and adding empty placeholders for
decorations with no matching result. The key parsing (``key.split("_")``) is
fragile, so the happy path and the "missing decoration" path are pinned here.
A pure function - no Neo4j.
"""
import json

from network.functions import order_structures


def test_order_structures_nests_by_depth_type_position():
    # key format is "position_depth_type".
    results = {"pos1_3_ADDITION": ["C"], "pos2_3_ADDITION": ["CC"]}
    decoration_list = [["pos1"], [], []]  # pos1 is an addition

    out = json.loads(order_structures(results, decoration_list))

    assert out["3"]["ADDITION"]["pos1"] == {
        "smiles": ["C"],
        "annotation": "ADD_DEC",
    }
    # pos2 is not in any decoration list -> BLANK.
    assert out["3"]["ADDITION"]["pos2"]["annotation"] == "BLANK"


def test_order_structures_annotates_each_decoration_type():
    results = {
        "a_1_T": ["C"],
        "b_1_T": ["N"],
        "c_1_T": ["O"],
    }
    decoration_list = [["a"], ["b"], ["c"]]  # add, del, link

    out = json.loads(order_structures(results, decoration_list))

    assert out["1"]["T"]["a"]["annotation"] == "ADD_DEC"
    assert out["1"]["T"]["b"]["annotation"] == "DEL_DEC"
    assert out["1"]["T"]["c"]["annotation"] == "LINK_DEC"


def test_order_structures_adds_empty_for_missing_decorations():
    """A decoration with no matching result gets an empty placeholder at depth -1."""
    results: dict = {}  # nothing came back from the graph
    decoration_list = [["addme"], ["delme"], ["linkme"]]

    out = json.loads(order_structures(results, decoration_list))

    assert out["-1"]["ADDITION"]["addme"] == {"smiles": [], "annotation": "ADD_MISS"}
    assert out["-1"]["DELETION"]["delme"] == {"smiles": [], "annotation": "DEL_MISS"}
    assert out["-1"]["LINKER"]["linkme"] == {"smiles": [], "annotation": "LINK_MISS"}
