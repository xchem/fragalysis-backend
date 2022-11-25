"""
Encoded recipes for generalised reactions
concentration (mol/L)
density (g/mL)
"""
encoded_recipes = {
    "Amidation": {
        "intramolecular": True,
        "recipes": {
            "standard": {
                "yield": 85,
                "reactionSMARTS": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6](=[#8])-[#8]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCCP1(=O)OP(=O)(OP(=O)(O1)CCC)CCC",
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H3,H2,H1]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                        "Intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6](=[#8])-[#8]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 6,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCCP1(=O)OP(=O)(OP(=O)(O1)CCC)CCC",
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 7,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 40, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O.C(=O)(O)[O-].[Na+]",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "satNaHCO3/H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 40, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 14,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 16,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O.[Na+].[Cl-]",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "Brine",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 17,
                                "content": {
                                    "platetype": "workup1",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 18,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup1",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 76, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 6,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 19,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 7,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 20,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 8,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 21,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 22,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 23,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Amide schotten - baumann": {
        "intramolecular": True,
        "recipes": {
            "standard": {
                "yield": 85,
                "reactionSMARTS": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
                "referencess": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H2,H1:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:1](=[#8:2])-[#17]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                        "intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H2,H1:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Buchwald-Hartwig amination": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#7:2].[c:3]-[F,I,Br,Cl]>>[#7:2]-[c:3]",
                "references": "https://doi.org/10.1021/acs.oprd.0c00018",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)O.C1=CC=C([C-]=C1)C2=CC=CC=C2N.[Pd]",
                                            "quantity": {
                                                "value": 0.06,
                                                "unit": "moleq",
                                            },
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:3]-[F,I,Br,Cl]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 120,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 18, "unit": "hours"},
                                },
                            },
                        ],
                    },
                                  {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 14,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # # {
                    # #     "type": "analyse",
                    # #     "driver": "robot",
                    # #     "sessionnumber": 4,
                    # #     "actions": [
                    # #         {
                    # #             "type": "add",
                    # #             "actionnumber": 16,
                    # #             "content": {
                    # #                 "plates": {
                    # #                     "fromplatetype": "spefilter",
                    # #                     "toplatetype": "lcms",
                    # #                 },
                    # #                 "material": {
                    # #                     "SMARTS": None,
                    # #                     "SMILES": None,  # Product of reaction
                    # #                     "quantity": {"value": 10, "unit": "ul"},
                    # #                     "solvent": None,
                    # #                     "density": None,
                    # #                     "concentration": None,
                    # #                 },
                    # #             },
                    # #         },
                    # #         {
                    # #             "type": "add",
                    # #             "actionnumber": 17,
                    # #             "content": {
                    # #                 "plates": {
                    # #                     "fromplatetype": "solvent",
                    # #                     "toplatetype": "lcms",
                    # #                 },
                    # #                 "material": {
                    # #                     "SMARTS": None,
                    # #                     "SMILES": "CC#N",
                    # #                     "quantity": {"value": 80, "unit": "ul"},
                    # #                     "solvent": "ACN",
                    # #                     "density": None,
                    # #                     "concentration": None,
                    # #                 },
                    # #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Buchwald-Hartwig amidation with amide-like nucleophile": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#7:2]-[#6]=[#8].[c:3]-[F,I,Br,Cl]>>[#7:2]-[c:3]",
                "references": "https://doi.org/10.1021/acs.oprd.0c00018",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)O.C1=CC=C([C-]=C1)C2=CC=CC=C2N.[Pd]",
                                            "quantity": {
                                                "value": 0.06,
                                                "unit": "moleq",
                                            },
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:3]-[F,I,Br,Cl]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7:2]-[#6]=[#8]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 120,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 18, "unit": "hours"},
                                },
                            },
                        ],
                    },
                                      {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 14,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 4,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 16,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "spefilter",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 10, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 17,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 80, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Buchwald-Hartwig (thio)etherification": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[c:1]-[#8:2].[c:3]-[F,I,Br,Cl]>>[#8:2]-[c:3]",
                "references": "https://doi.org/10.1021/acs.oprd.0c00018",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)O.C1=CC=C([C-]=C1)C2=CC=CC=C2N.[Pd]",
                                            "quantity": {
                                                "value": 0.06,
                                                "unit": "moleq",
                                            },
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:3]-[F,I,Br,Cl]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:1]-[#8:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 120,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 18, "unit": "hours"},
                                },
                            },
                        ],
                    },
                                      {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 14,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 4,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 16,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "spefilter",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 10, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 17,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 80, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Ester amidation": {
        "intramolecular": True,
        "recipes": {
            "standard": {
                "yield": 80,
                "reactionSMARTS": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
                "references": ["https://doi.org/10.3390/molecules25051040"],
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6](=[#8])-[#8]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H3,H2,H1]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "Cl[Fe](Cl)Cl",
                                            "quantity": {
                                                "value": 0.14,
                                                "unit": "moleq",
                                            },
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                        "intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6](=[#8])-[#8]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "Cl[Fe](Cl)Cl",
                                            "quantity": {
                                                "value": 0.14,
                                                "unit": "moleq",
                                            },
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Mitsunobu aryl ether synthesis": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[c:1]-[#8;H:2].[#6]-[#8:3]>>[c:1]-[#8:2]-[#6:3]",
                "references": "To do",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:1]-[#8;H:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6]-[#8:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3",
                                            "quantity": {"value": 2.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)OC(=O)N=NC(=O)OC(C)C",
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "Toluene",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 20,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 4, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 6,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 13,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "workup1",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 5, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 14,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 100, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Mitsunobu reaction with amine alcohol and thioalcohol": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[c:1]-[#7:2].[#6]-[#8:3]>>[c:1]-[#7:2]-[#6:3]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:1]-[#7:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6]-[#8:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3",
                                            "quantity": {"value": 2.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)OC(=O)N=NC(=O)OC(C)C",
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "Toluene",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 20,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 4, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 20, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 6,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 13,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "workup1",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 10, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 14,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 80, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "N-nucleophilic aromatic substitution": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:3]-[#7;H2,H1:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:1]-[F,Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 4,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 140,
                                        "unit": "degC",
                                    },  # 150, 160 and 170 test
                                    "duration": {"value": 4, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 40, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 40, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 10,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O.[Na+].[Cl-]",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "Brine",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 13,
                                "content": {
                                    "platetype": "workup1",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 14,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup1",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 25,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 76, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 6,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup1",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 5, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 16,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Nucleophilic substitution with amine": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#6:2]-[F,Cl,Br,I].[#7;H1,H2,H3:1]>>[#7:1]-[#6:2]",
                "references": "To do",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:2]-[F,Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H1,H2,H3:1]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 4,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 120,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 4, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 10,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 6,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 12,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "workup1",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 10, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 13,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 80, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Reductive amination": {
        "intramolecular": True,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H3,H2,H1:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(CC)CC",
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:2](=[#8])(-[#6:1])",
                                            "SMILES": None,
                                            "quantity": {"value": 1.1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                        "intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:2](=[#8])(-[#6:1])",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(CC)CC",
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 20, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "reaction",
                        "driver": "human",
                        "sessionnumber": 3,
                        "continuation": True,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 7,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[BH-](OC(=O)C)(OC(=O)C)OC(=O)C.[Na+]",
                                            "quantity": {"value": 1.4, "unit": "moleq"},
                                            "solvent": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                        "intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 8,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[BH-](OC(=O)C)(OC(=O)C)OC(=O)C.[Na+]",
                                            "quantity": {"value": 1.4, "unit": "moleq"},
                                            "solvent": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 9,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 20, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 6,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 12,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 7,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 14,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 15,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 16,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 8,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 17,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "workup1",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 5, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 18,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 100, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Sonogashira coupling": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[CH:1].[c:2]-[Cl,Br,I]>>[C:1]-[c:2]",
                "references": [
                    "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                    "https://pubs.acs.org/doi/pdf/10.1021/ol035632f",
                    "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                ],
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                            "quantity": {"value": 0.1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.1,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:2]-[Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[CH:1]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[Cu]I",
                                            "quantity": {"value": 0.1, "unit": "moleq"},
                                            "solvent": "ACN",
                                            "concentration": 0.1,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
                "NMP": {
                    "yield": 70,
                    "reactionSMARTS": "[CH:1].[c:2]-[Cl,Br,I]>>[C:1]-[c:2]",
                    "references": [
                        "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                        "https://pubs.acs.org/doi/pdf/10.1021/ol035632f",
                        "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                    ],
                    "actionsessions": [
                        {
                            "type": "reaction",
                            "driver": "robot",
                            "sessionnumber": 1,
                            "intermolecular": {
                                "actions": [
                                    {
                                        "type": "add",
                                        "actionnumber": 1,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "solvent",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": None,
                                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                                "quantity": {
                                                    "value": 0.1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.1,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 2,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "solvent",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": "[c:2]-[Cl,Br,I]",
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 3,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "solvent",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": "[CH:1]",
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 2,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 4,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "solvent",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": None,
                                                "SMILES": "[Cu]I",
                                                "quantity": {
                                                    "value": 0.1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "ACN",
                                                "concentration": 0.1,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 5,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "solvent",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": None,
                                                "SMILES": "CCN(C(C)C)C(C)C",
                                                "quantity": {
                                                    "value": 2,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "DMA",
                                                "concentration": 10,
                                            },
                                        },
                                    },
                                ],
                            },
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 2,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 6,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 80, "unit": "degC"},
                                        "duration": {"value": 12, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 3,
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 7,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "workup2",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C[S](C)=O",
                                            "quantity": {
                                                "value": 20,
                                                "unit": "masseq",
                                            },  # Check conc for XChem
                                            "solvent": "DMSO",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 4,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 8,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 25, "unit": "degC"},
                                        "duration": {"value": 1, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 5,
                            "actions": [
                                {
                                    "type": "extract",
                                    "actionnumber": 9,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 10, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 10,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC#N",
                                            "quantity": {"value": 80, "unit": "ul"},
                                            "solvent": "ACN",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "extract",
                                    "actionnumber": 11,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "xchem",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 150, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    ],
                },
            },
        },
    },
    "Sp2-sp2 Suzuki coupling": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 75,
                "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[c:1]-[F,Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:2]-[B]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)O.C1=CC=C([C-]=C1)C2=CC=CC=C2N.[Pd]",
                                            # Smiles for XPhosPdG3
                                            "quantity": {
                                                "value": 0.06,
                                                "unit": "moleq",
                                            },
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "starting material",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "sessionnumber": 2,
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 120, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 8,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 11,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 12,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 14,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "spefilter",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": None,
                                        "quantity": {"value": 30, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 4,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 16,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "spefilter",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 10, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 17,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 80, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    ##############################################################################################################
    ################ Reactions we will not use for, are not working on and still have to test ###########################################################
    "Sulfonamide schotten-baumann": {
        "intramolecular": True,
        "recipes": {
            "standard": {
                "yield": 75,
                "reactionSMARTS": "[#16:5](=[#8])(=[#8:7])-[#17].[#6]-[#7;H2,H1:2]>>[#16:5](=[#8])(=[#8:7])-[#7:2]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "MeOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6]-[#7;H2,H1:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "MeOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                        "intramolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "MeOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ]
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 4,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
                "DMA": {
                    "reactionSMARTS": "[#16:5](=[#8])(=[#8:7])-[#17].[#6]-[#7;H2,H1:2]>>[#16:5](=[#8])(=[#8:7])-[#7:2]",
                    "references": None,
                    "actionsessions": [
                        {
                            "type": "reaction",
                            "driver": "robot",
                            "sessionnumber": 1,
                            "intermolecular": {
                                "actions": [
                                    {
                                        "type": "add",
                                        "actionnumber": 1,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "DMA",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 2,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": "[#6]-[#7;H2,H1:2]",
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 1.2,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "DMA",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                ],
                            },
                            "intramolecular": {
                                "actions": [
                                    {
                                        "type": "add",
                                        "actionnumber": 3,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "DMA",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                ]
                            },
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 2,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 4,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 25, "unit": "degC"},
                                        "duration": {"value": 12, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 3,
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 5,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "workup2",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C[S](C)=O",
                                            "quantity": {
                                                "value": 20,
                                                "unit": "masseq",
                                            },  # Check conc for XChem
                                            "solvent": "DMSO",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 4,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 6,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 25, "unit": "degC"},
                                        "duration": {"value": 1, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 5,
                            "actions": [
                                {
                                    "type": "extract",
                                    "actionnumber": 7,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 10, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 8,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC#N",
                                            "quantity": {"value": 80, "unit": "ul"},
                                            "solvent": "ACN",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "extract",
                                    "actionnumber": 9,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "xchem",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 150, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    ],
                },
            },
        },
    },
    "Williamson ether synthesis": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[#6:1]-[#8;H:2].[#6:3]-[F,Cl,Br,I]>>[#6:1]-[#8:2]-[#6:3]",
                "references": "To do",
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:1]-[#8;H:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "reaction",
                        "driver": "human",
                        "sessionnumber": 2,
                        "continuation": True,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[H-].[Na+]",
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "continuation": True,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#6:3]-[F,Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 4,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {
                                        "value": 120,
                                        "unit": "degC",
                                    },
                                    "duration": {"value": 8, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 6,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 7,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CCOC(C)=O",
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": "EtOAc",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 10,
                                "content": {
                                    "platetype": "reaction",
                                    "repetitions": {"value": 3},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "bottomlayerquantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },
                                        "layer": "top",
                                        "SMILES": None,
                                        "quantity": {"value": 20, "unit": "masseq"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    # {
                    #     "type": "analyse",
                    #     "driver": "robot",
                    #     "sessionnumber": 8,
                    #     "actions": [
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 12,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "workup1",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": None,  # Product of reaction
                    #                     "quantity": {"value": 5, "unit": "ul"},
                    #                     "solvent": None,
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #         {
                    #             "type": "add",
                    #             "actionnumber": 13,
                    #             "content": {
                    #                 "plates": {
                    #                     "fromplatetype": "solvent",
                    #                     "toplatetype": "lcms",
                    #                 },
                    #                 "material": {
                    #                     "SMARTS": None,
                    #                     "SMILES": "CC#N",
                    #                     "quantity": {"value": 100, "unit": "ul"},
                    #                     "solvent": "ACN",
                    #                     "density": None,
                    #                     "concentration": None,
                    #                 },
                    #             },
                    #         },
                    #     ],
                    # },
                ],
            },
        },
    },
    "Boc protection": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 85,
                "reactionSMARTS": "[#7:2].[#8:3]-[#6:1](=[#8:4])-[#8:5]>>[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": [
                                                "[#8:3]-[#6:1](=[#8:4])-[#8:5]-[#6](-[#6])(-[#6])(-[#6])"
                                            ],
                                            "SMILES": None,
                                            "quantity": {"value": 3, "unit": "moleq"},
                                            "solvent": "ACN",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "n1ccc(N(C)C)cc1",
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "ACN",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": "[#7;H2:2]",
                                            "SMILES": None,
                                            "quantity": {
                                                "value": 1,
                                                "unit": "moleq",
                                            },
                                            "solvent": "ACN",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 4,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 0.5, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 6,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Boc deprotection": {  # this does not work yet, SMARTS needs rethinking
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 85,
                "reactionSMARTS": "[#6:9]-[#6:8]-[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6]).[#1]-[#17]>>[#6:9]-[#6:8]-[#7:2]",
                "references": None,
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": [
                                                "[#6:9]-[#6:8]-[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])"
                                            ],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOAc",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[#1]-[#17]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOAc",
                                            "concentration": 6,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 3,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 4,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 7,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    
    "Heck coupling": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[c:2]-[F,Cl,Br,I].[CX3;H2:1]>>[c:2]-[CX3;H1:1]",  # SMARTS for terminal alkenes only
                "references": ["Platinum Metals Rev., 1999, 43, (4), 138"],
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[CX3;H2:1]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                            "quantity": {"value": 0.1, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.1,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(CC)CC",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 0.73,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
}

####Old recipes
{
    "Buchwald-Hartwig thiolation": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 70,
                "reactionSMARTS": "[c:2]-[F,Cl,Br,I].[#6:1]-[#16;H1]>>[c:2]-[#16:3]-[#6:1]",
                "references": ["https://doi.org/10.1002/ejoc.201001393"],
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                            "quantity": {"value": 0.1, "unit": "moleq"},
                                            "solvent": "NMP",  # try NMP but references says to use 2-MeTHF
                                            "concentration": 0.1,
                                            # try same as Suzuki - similar mechanism
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[#6:1]-[#16;H1]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 120, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
            },
        },
    },
    "Buchwald-Hartwig amination": {
        "intramolecular": False,
        "recipes": {
            "standard": {
                "yield": 75,
                "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[#7;H2,H1:3]>>[c:1]-[#7:3]-[#6:2]",
                "referencess": [
                    "https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981"
                ],
                "actionsessions": [
                    {
                        "type": "reaction",
                        "driver": "robot",
                        "sessionnumber": 1,
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 1,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                            "quantity": {"value": 0.1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.1,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 2,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnmber": 3,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": ["[#6:1]-[#7;H2,H1:3]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 4,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "startingmaterial",
                                            "toplatetype": "reaction",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 10,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 2,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 5,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "add",
                                "actionnumber": 6,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "C[S](C)=O",
                                        "quantity": {
                                            "value": 20,
                                            "unit": "masseq",
                                        },  # Check conc for XChem
                                        "solvent": "DMSO",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "stir",
                        "driver": "human",
                        "sessionnumber": 4,
                        "actions": [
                            {
                                "type": "stir",
                                "actionnumber": 7,
                                "content": {
                                    "platetype": "reaction",
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 1, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    {
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 5,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 8,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "solvent",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CC#N",
                                        "quantity": {"value": 80, "unit": "ul"},
                                        "solvent": "ACN",
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 10,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup2",
                                        "toplatetype": "xchem",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 150, "unit": "ul"},
                                        "solvent": None,
                                        "density": None,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                ],
                "NMP": {
                    "yield": 75,
                    "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[#7;H2,H1:3]>>[c:1]-[#7:3]-[#6:2]",
                    "references": [
                        "https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981"
                    ],
                    "actionsessions": [
                        {
                            "type": "reaction",
                            "driver": "robot",
                            "sessionnumber": 1,
                            "intermolecular": {
                                "actions": [
                                    {
                                        "type": "add",
                                        "actionnumber": 1,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": None,
                                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                                "quantity": {
                                                    "value": 0.1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.1,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 2,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 1,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnmber": 3,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": ["[#6:1]-[#7;H2,H1:3]"],
                                                "SMILES": None,
                                                "quantity": {
                                                    "value": 2,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "NMP",
                                                "concentration": 0.5,
                                            },
                                        },
                                    },
                                    {
                                        "type": "add",
                                        "actionnumber": 4,
                                        "content": {
                                            "plates": {
                                                "fromplatetype": "startingmaterial",
                                                "toplatetype": "reaction",
                                            },
                                            "material": {
                                                "SMARTS": None,
                                                "SMILES": "C1CCN2CCCN=C2CC1",
                                                "quantity": {
                                                    "value": 2,
                                                    "unit": "moleq",
                                                },
                                                "solvent": "DMA",
                                                "concentration": 10,
                                            },
                                        },
                                    },
                                ],
                            },
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 2,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 5,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 80, "unit": "degC"},
                                        "duration": {"value": 12, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 3,
                            "actions": [
                                {
                                    "type": "add",
                                    "actionnumber": 6,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "workup2",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C[S](C)=O",
                                            "quantity": {
                                                "value": 20,
                                                "unit": "masseq",
                                            },  # Check conc for XChem
                                            "solvent": "DMSO",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                        {
                            "type": "stir",
                            "driver": "human",
                            "sessionnumber": 4,
                            "actions": [
                                {
                                    "type": "stir",
                                    "actionnumber": 7,
                                    "content": {
                                        "platetype": "reaction",
                                        "temperature": {"value": 25, "unit": "degC"},
                                        "duration": {"value": 1, "unit": "hours"},
                                    },
                                },
                            ],
                        },
                        {
                            "type": "analyse",
                            "driver": "robot",
                            "sessionnumber": 5,
                            "actions": [
                                {
                                    "type": "extract",
                                    "actionnumber": 8,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 10, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "actionnumber": 9,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "solvent",
                                            "toplatetype": "lcms",
                                        },
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC#N",
                                            "quantity": {"value": 80, "unit": "ul"},
                                            "solvent": "ACN",
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                                {
                                    "type": "extract",
                                    "actionnumber": 10,
                                    "content": {
                                        "plates": {
                                            "fromplatetype": "workup2",
                                            "toplatetype": "xchem",
                                        },
                                        "material": {
                                            "layer": "bottom",
                                            "SMILES": None,  # Product of reaction
                                            "quantity": {"value": 150, "unit": "ul"},
                                            "solvent": None,
                                            "density": None,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    ],
                },
            },
        },
    },
}
