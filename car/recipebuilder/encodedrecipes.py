"""
Encoded recipes for generalised reactions
concentration (mol/L)
density (g/mL)
"""
encoded_recipes = {
    "Amidation": {
        "reactionSMARTS": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
        "intramolecular": True,
        "recipes": {
            "standard": {
                "reference": None,
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
                                            "solvent": None,
                                            "density": 0.74,
                                            "concentration": None,
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
                                            "solvent": None,
                                            "density": 0.74,
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
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": 1.00,
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
                                        "SMILES": "CO",
                                        "quantity": {
                                            "value": 80,
                                            "unit": "ul",
                                        },  # Check conc for XChem
                                        "solvent": "MeOH",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 4,
                        "actions": [
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
                                        "SMILES": "C(Cl)Cl",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "DCM",
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
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 13,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,
                                        "quantity": {"value": 23, "unit": "masseq"},
                                        "solvent": None,
                                        "density": 1.00,
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
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O.[Na+].[Cl-]",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "H2O/NaCl",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 15,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup1",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,
                                        "quantity": {"value": 23, "unit": "masseq"},
                                        "solvent": None,
                                        "density": 1.00,
                                        "concentration": None,
                                    },
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
                                "type": "add",
                                "actionnumber": 16,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 17,
                                "content": {
                                    "platetype": "workup2",
                                    "repetitions": {"value": 3},
                                    "quantity": {"value": 20, "unit": "masseq"},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 18,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 19,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 20,
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
                                        "density": 1.00,
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
        "reactionSMARTS": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
        "intramolecular": True,
        "recipes": {
            "standard": {
                "referencess": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 0.74,
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
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 3.5, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 0.74,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 4,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 5,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Buchwald-Hartwig amination": {
        "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[#7;H2,H1:3]>>[c:1]-[#7:3]-[#6:2]",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "referencess": [
                    "https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981"
                ],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
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
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 1.02,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
            "NMP": {
                "references": ["https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981"],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": ["[#6:1]-[#7;H2,H1:3]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 1.02,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
            "EtOH": {
                "references": ["https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981"],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
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
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 1.02,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Ester amidation": {
        "reactionSMARTS": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
        "intramolecular": True,
        "recipes": {
            "standard": {
                "references": ["https://doi.org/10.3390/molecules25051040"],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
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
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 4,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 5,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "N-nucleophilic aromatic substitution": {
        "reactionSMARTS": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
        "intramolecular": False,
        "recipes": {
            "standard": {
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
                                            "solvent": None,
                                            "density": 0.74,
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
                        "type": "analyse",
                        "driver": "robot",
                        "sessionnumber": 3,
                        "actions": [
                            {
                                "type": "extract",
                                "actionnumber": 5,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,  # Product of reaction
                                        "quantity": {"value": 10, "unit": "ul"},
                                        "solvent": None,
                                        "density": 1.00,
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
                                        "toplatetype": "lcms",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "CO",
                                        "quantity": {
                                            "value": 80,
                                            "unit": "ul",
                                        },  # Check conc for XChem
                                        "solvent": "MeOH",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                        ],
                    },
                    {
                        "type": "workup",
                        "driver": "robot",
                        "sessionnumber": 4,
                        "actions": [
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
                                        "SMILES": "C(Cl)Cl",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "DCM",
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
                                        "toplatetype": "reaction",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "H2O",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 9,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "reaction",
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,
                                        "quantity": {"value": 23, "unit": "masseq"},
                                        "solvent": None,
                                        "density": 1.00,
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
                                        "toplatetype": "workup1",
                                    },
                                    "material": {
                                        "SMARTS": None,
                                        "SMILES": "O.[Na+].[Cl-]",
                                        "quantity": {"value": 25, "unit": "masseq"},
                                        "solvent": "H2O/NaCl",
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 11,
                                "content": {
                                    "plates": {
                                        "fromplatetype": "workup1",
                                        "toplatetype": "workup2",
                                    },
                                    "material": {
                                        "layer": "bottom",
                                        "SMILES": None,
                                        "quantity": {"value": 23, "unit": "masseq"},
                                        "solvent": None,
                                        "density": 1.00,
                                        "concentration": None,
                                    },
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
                                "type": "add",
                                "actionnumber": 12,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "mix",
                                "actionnumber": 13,
                                "content": {
                                    "platetype": "workup2",
                                    "repetitions": {"value": 3},
                                    "quantity": {"value": 20, "unit": "masseq"},
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 14,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "add",
                                "actionnumber": 15,
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
                                        "density": 1.00,
                                        "concentration": None,
                                    },
                                },
                            },
                            {
                                "type": "extract",
                                "actionnumber": 16,
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
                                        "density": 1.00,
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
    "Reductive amination": {
        "reactionSMARTS": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
        "intramolecular": True,
        "recipes": {
            "standard": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": "[#7;H2,H1:3]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[Na+].CC(=O)O[BH-](OC(C)=O)OC(C)=O",
                                            "quantity": {"value": 1.4, "unit": "moleq"},
                                            "solvent": "MeOH",
                                            "concentration": 0.25,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(=O)O",
                                            "quantity": {"value": 1.0, "unit": "moleq"},
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
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "[Na+].CC(=O)O[BH-](OC(C)=O)OC(C)=O",
                                            "quantity": {"value": 1.4, "unit": "moleq"},
                                            "solvent": "MeOH",
                                            "concentration": 0.25,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(=O)O",
                                            "quantity": {"value": 1.0, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Sonogashira coupling": {
        "reactionSMARTS": "[CH:1].[c:2]-[Cl,Br,I]>>[C:1]-[c:2]",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": [
                    "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                    "https://pubs.acs.org/doi/pdf/10.1021/ol035632f",
                    "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                ],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": ["[c:2]-[Cl,Br,I]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": ["[CH:1]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 4,
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
                                    "content": {
                                        "number": 5,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 0.74,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 6,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 7,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
            "NMP": {
                "references": [
                    "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                    "https://pubs.acs.org/doi/pdf/10.1021/ol035632f",
                    "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                ],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": ["[c:2]-[Cl,Br,I]"],
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": ["[CH:1]"],
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "NMP",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 4,
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
                                    "content": {
                                        "number": 5,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CCN(C(C)C)C(C)C",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 0.74,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 6,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 7,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Sp2-sp2 Suzuki coupling": {
        "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
                                        "material": {
                                            "SMARTS": "[c:1]-[F,Cl,Br,I]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": "[#6:2]-[B]",
                                            "SMILES": None,
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)[O-].C1=CC=C(C=C1)C2=CC=CC=C2N.[Pd]",
                                            # Smiles for XPhosPdG3
                                            "quantity": {
                                                "value": 0.1,
                                                "unit": "moleq",
                                            },  # 10mol% catalyst
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": "EtOH",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 100, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    ##############################################################################################################
    ################ Reactions we will not use for, are not working on and still have to test ###########################################################
    "Sulfonamide schotten-baumann": {
        "reactionSMARTS": "[#16:5](=[#8])(=[#8:7])-[#17].[#6]-[#7;H2,H1:2]>>[#16:5](=[#8])(=[#8:7])-[#7:2]",
        "intramolecular": True,
        "recipes": {
            "standard": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 1,
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
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 3,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 4,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
            "DMA": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
                                        "material": {
                                            "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 2,
                                        "material": {
                                            "SMARTS": "[#6]-[#7;H2,H1:2]",
                                            "SMILES": None,
                                            "quantity": {"value": 1.2, "unit": "moleq"},
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
                                    "content": {
                                        "number": 1,
                                        "material": {
                                            "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                            "SMILES": None,
                                            "quantity": {"value": 1, "unit": "moleq"},
                                            "solvent": "DMA",
                                            "concentration": 0.5,
                                        },
                                    },
                                },
                            ]
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 3,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 4,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Boc protection": {
        "reactionSMARTS": "[#7:2].[#8:3]-[#6:1](=[#8:4])-[#8:5]>>[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
                                        "material": {
                                            "SMARTS": ["[#7;H2:2]"],
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
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 4,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 0.5, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 5,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Boc deprotection": {  # this does not work yet, SMARTS needs rethinking
        "reactionSMARTS": "[#6:9]-[#6:8]-[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6]).[#1]-[#17]>>[#6:9]-[#6:8]-[#7:2]",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": None,
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 3,
                                    "temperature": {"value": 25, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 4,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Buchwald-Hartwig thiolation": {
        "reactionSMARTS": "[c:2]-[F,Cl,Br,I].[#6:1]-[#16;H1]>>[c:2]-[#16:3]-[#6:1]",
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": ["https://doi.org/10.1002/ejoc.201001393"],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
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
                                    "content": {
                                        "number": 4,
                                        "material": {
                                            "SMARTS": None,
                                            "SMILES": "C1CCN2CCCN=C2CC1",
                                            "quantity": {"value": 2, "unit": "moleq"},
                                            "solvent": None,
                                            "density": 1.02,
                                            "concentration": None,
                                        },
                                    },
                                },
                            ],
                        },
                    },
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
    "Heck coupling": {
        "reactionSMARTS": "[c:2]-[F,Cl,Br,I].[CX3;H2:1]>>[c:2]-[CX3;H1:1]",  # SMARTS for terminal alkenes only
        "intramolecular": False,
        "recipes": {
            "standard": {
                "references": ["Platinum Metals Rev., 1999, 43, (4), 138"],
                "actionsessions": {
                    "reaction": {
                        "driver": "robot",
                        "intermolecular": {
                            "actions": [
                                {
                                    "type": "add",
                                    "content": {
                                        "number": 1,
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
                                    "content": {
                                        "number": 2,
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
                                    "content": {
                                        "number": 3,
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
                                    "content": {
                                        "number": 4,
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
                    "stir": {
                        "driver": "human",
                        "actions": [
                            {
                                "type": "stir",
                                "content": {
                                    "number": 5,
                                    "temperature": {"value": 80, "unit": "degC"},
                                    "duration": {"value": 12, "unit": "hours"},
                                },
                            },
                        ],
                    },
                    "analyse": {
                        "driver": "robot",
                        "actions": [
                            {
                                "type": "analyse",
                                "content": {
                                    "number": 6,
                                    "method": "LCMS",
                                    "samplevolume": {"value": 10, "unit": "ul"},
                                    "solvent": "formic acid/acetonitrile",
                                    "solventvolume": {"value": 80, "unit": "ul"},
                                },
                            },
                        ],
                    },
                },
            },
        },
    },
}
