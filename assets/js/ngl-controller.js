// modified for local use - add rcsb url
const rcsbUrl = 'http://www.rcsb.org';

const NglController = function(params) {
    var signals = {
        taskCountChanged: new NGL.Signal(),
        fullscreenChanged: new NGL.Signal(),
        structureChanged: new NGL.Signal(),
        symmetryDataLoaded: new NGL.Signal(),
        validationDataLoaded: new NGL.Signal(),
        colorSchemeChanged: new NGL.Signal(),
        modelChanged: new NGL.Signal(),
        hydrogenVisibilityChanged: new NGL.Signal(),
        ionVisibilityChanged: new NGL.Signal(),
        waterVisibilityChanged: new NGL.Signal(),
        clashVisibilityChanged: new NGL.Signal(),
        qualityChanged: new NGL.Signal(),
        assemblyChanged: new NGL.Signal(),
        symmetryChanged: new NGL.Signal(),
        interactionChanged: new NGL.Signal(),
        styleChanged: new NGL.Signal(),
        ligandStyleChanged: new NGL.Signal(),
        focusChanged: new NGL.Signal(),
        clicked: new NGL.Signal(),
        hovered: new NGL.Signal()
    };
    var pdbid;
    var reduced;
    var structureComponent;
    var symmetryBuffer;
    var symmetryData = {};
    var spatialHash;
    var clashBuffer;
    var validationData;
    var validationDataLoading;
    var atomCount;
    var instanceCount;
    var axesRepr;
    var isBackboneOnly;

    // assign properties from params
    var p = Object.assign({}, params);
    var colorScheme = p.colorScheme || "chainname";
    var assembly = p.assembly || "BU1";
    var style = p.style !== undefined ? p.style : "cartoon";
    var ligandStyle = p.ligandStyle !== undefined ? p.ligandStyle : "ball+stick";
    var model = p.model || 0;
    var symmetry = p.symmetry || 0;
    var interaction = p.interaction || false;
    var hydrogenVisibility = p.hydrogenVisibility === undefined ? true : p.hydrogenVisibility;
    var ionVisibility = p.ionVisibility === undefined ? false : p.ionVisibility;
    var waterVisibility = p.waterVisibility === undefined ? false : p.waterVisibility;
    var clashVisibility = p.clashVisibility === undefined ? false : p.clashVisibility;
    var quality = p.quality === undefined ? "auto" : p.quality;
    var spin = p.spin === undefined ? false : p.spin;
    var focus = p.focus === undefined ? 0 : p.focus;
    var hasStructureFactors = p.hasStructureFactors === undefined ? false : p.hasStructureFactors;

    // rl - stage already created in ngl-ui.js
    //var stage = new NGL.Stage(id, {
    //    backgroundColor: "white",
    //    hoverTimeout: 500
    //});
    //this.stage = stage;
    var stage;
    var tasks;
    // rl - setStage (and tasks)  - called from ngl-ui.js after stage is created
    function setStage(_stage) {
        stage = _stage;
        stage.signals.fullscreenChanged.add(function(fullscreen) {
            signals.fullscreenChanged.dispatch(fullscreen);
        });
        stage.signals.clicked.add(function(pickingData) {
            signals.clicked.dispatch(pickingData);
        });
        stage.signals.hovered.add(function(pickingData) {
            signals.hovered.dispatch(pickingData);
        });
        window.addEventListener("resize", function() {
            stage.handleResize();
        }, false);
        tasks = new NGL.Counter();
        tasks.listen(stage.tasks);
        tasks.signals.countChanged.add(function(delta, count) {
            signals.taskCountChanged.dispatch(delta, count);
        });
        setSpin(spin);

        // rl - added from.jsp
        signals.clicked.add( function( pd ){
            console.log('pd=' + pd );
            var msg;
            if( pd.atom ){
                msg = "Atom: " + pd.atom.qualifiedName();
            } else if( pd.bond ) {
                msg = "Bond: " +
                    pd.bond.atom1.qualifiedName() + " - " + pd.bond.atom2.qualifiedName();
            } else {
                msg = "";
            }
            document.getElementById('pickingInfo').innerHTML = msg;
        } );
    }

    var polymerReprDict = {};
    var polymerReprDefs = {
        "unitcell": {
            disableImpostor: true,
            radiusSegments: 16
        },
        "cartoon": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            aspectRatio: 5,
            scale: 0.7,
            quality: "custom",
            subdiv: function() {
                if (quality === "auto") {
                    if (atomCount < 15000) {
                        return 12;
                    } else if (atomCount < 70000) {
                        return 6;
                    } else {
                        return 3;
                    }
                } else {
                    if (quality === "high") {
                        return 12;
                    } else if (quality === "medium") {
                        return 6;
                    } else {
                        return 3;
                    }
                }
            },
            radialSegments: function() {
                if (quality === "auto") {
                    if (atomCount < 15000) {
                        return 20;
                    } else if (atomCount < 70000) {
                        return 10;
                    } else {
                        return 6;
                    }
                } else {
                    if (quality === "high") {
                        return 20;
                    } else if (quality === "medium") {
                        return 10;
                    } else {
                        return 6;
                    }
                }
            },
            sele: function() {
                var sele = "";
                if (model !== "all") {
                    sele += "/" + model;
                }
                return sele;
            }
        },
        "base": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            quality: "custom",
            sphereDetail: function() {
                if (quality === "auto") {
                    return atomCount < 15000 ? 1 : 0;
                } else {
                    if (quality === "high") {
                        return 1;
                    } else if (quality === "medium") {
                        return 1;
                    } else {
                        return 0;
                    }
                }
            },
            radialSegments: function() {
                if (quality === "auto") {
                    if (atomCount < 15000) {
                        return 20;
                    } else if (atomCount < 70000) {
                        return 10;
                    } else {
                        return 5;
                    }
                } else {
                    if (quality === "high") {
                        return 20;
                    } else if (quality === "medium") {
                        return 10;
                    } else {
                        return 5;
                    }
                }
            },
            sele: function() {
                var sele = "polymer and nucleic";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                return sele;
            }
        },
        "backbone": {
            lineOnly: function() {
                if (quality === "auto") {
                    return atomCount > 250000;
                } else {
                    return quality === "low";
                }
            },
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            scale: 2.0,
            sele: function() {
                var sele = "";
                if (model !== "all") {
                    sele += "/" + model;
                }
                return sele;
            }
        },
        "surface": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            surfaceType: "sas",
            probeRadius: 1.4,
            useWorker: true,
            scaleFactor: function() {
                var sf;
                if (quality === "low") {
                    sf = 0.1;
                } else if (quality === "medium") {
                    sf = 0.7;
                } else if (quality === "high") {
                    sf = 1.7;
                } else {
                    sf = Math.min(1.5, Math.max(0.1, 50000 / atomCount));
                }
                return sf;
            },
            sele: function() {
                var sele = "polymer and ( protein or nucleic )";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                return sele;
            }
        },
        "spacefill": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            quality: "custom",
            sphereDetail: function() {
                if (quality === "auto") {
                    return atomCount < 15000 ? 1 : 0;
                } else {
                    if (quality === "high") {
                        return 1;
                    } else if (quality === "medium") {
                        return 1;
                    } else {
                        return 0;
                    }
                }
            },
            sele: function() {
                var sele = "polymer and ( protein or nucleic )";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                return sele;
            }
        },
        "licorice": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            multipleBond: "symmetric",
            quality: "custom",
            sphereDetail: function() {
                if (quality === "auto") {
                    return atomCount < 15000 ? 1 : 0;
                } else {
                    if (quality === "high") {
                        return 1;
                    } else if (quality === "medium") {
                        return 1;
                    } else {
                        return 0;
                    }
                }
            },
            radialSegments: function() {
                if (quality === "auto") {
                    if (atomCount < 15000) {
                        return 20;
                    } else if (atomCount < 70000) {
                        return 10;
                    } else {
                        return 5;
                    }
                } else {
                    if (quality === "high") {
                        return 20;
                    } else if (quality === "medium") {
                        return 10;
                    } else {
                        return 5;
                    }
                }
            },
            sele: function() {
                var sele = "polymer and ( protein or nucleic )";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                return sele;
            }
        },
        "ball+stick": {
            colorScheme: getColorScheme,
            colorScale: getColorScale,
            multipleBond: "symmetric",
            quality: "custom",
            sphereDetail: function() {
                if (quality === "auto") {
                    return atomCount < 15000 ? 1 : 0;
                } else {
                    if (quality === "high") {
                        return 1;
                    } else if (quality === "medium") {
                        return 1;
                    } else {
                        return 0;
                    }
                }
            },
            radialSegments: function() {
                if (quality === "auto") {
                    if (atomCount < 15000) {
                        return 20;
                    } else if (atomCount < 70000) {
                        return 10;
                    } else {
                        return 5;
                    }
                } else {
                    if (quality === "high") {
                        return 20;
                    } else if (quality === "medium") {
                        return 10;
                    } else {
                        return 5;
                    }
                }
            },
            sele: function() {
                var sele;
                if (!reduced && validationData) {
                    sele = "( " + validationData.clashSele + " )";
                    if (model !== "all") {
                        sele += " and /" + model;
                    }
                    if (hydrogenVisibility === false) {
                        sele += " and not hydrogen";
                    }
                } else {
                    sele = "NONE";
                }
                return sele;
            }
        }
    };
    var ligandReprDict = {};
    var ligandReprDefs = {
        "spacefill": {
            colorScheme: getLigandColorScheme,
            quality: function() {
                return "medium";
            },
            sele: function() {
                var sele = "( not polymer or not ( protein or nucleic ) )";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                if (ionVisibility === false) {
                    sele += " and not ion";
                }
                if (waterVisibility === false) {
                    sele += " and not water";
                }
                return sele;
            }
        },
        "ball+stick": {
            multipleBond: "symmetric",
            colorScheme: getLigandColorScheme,
            quality: function() {
                return "medium";
            },
            scale: 2.5,
            aspectRatio: 1.2,
            sele: function() {
                var sele = "( not polymer or not ( protein or nucleic ) )";
                if (model !== "all") {
                    sele += " and /" + model;
                }
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                if (ionVisibility === false) {
                    sele += " and not ion";
                }
                if (waterVisibility === false) {
                    sele += " and not water";
                }
                return sele;
            }
        }
    };
    var interactionReprDict = {};
    var interactionReprDefs = {
        "licorice": {
            multipleBond: "symmetric",
            colorScheme: getLigandColorScheme,
            colorValue: "lightsteelblue",
            quality: function() {
                return "medium";
            },
            scale: 1.5,
            sele: function() {
                if (interaction && structureComponent) {
                    var s = structureComponent.structure;
                    var withinSele = new NGL.Selection(interaction + (model !== "all" ? " and /" + model : ""));
                    var as = s.getAtomSetWithinSelection(withinSele, 5);
                    if (model !== "all") {
                        as.intersection(s.getAtomSet(new NGL.Selection("/" + model)));
                    }
                    var asGroup = s.getAtomSetWithinGroup(as);
                    var sele = asGroup.toSeleString();
                    if (hydrogenVisibility === false) {
                        sele += " and not hydrogen";
                    }
                    if (ionVisibility === false) {
                        sele += " and not ion";
                    }
                    if (waterVisibility === false) {
                        sele += " and not water";
                    }
                    if (model !== "all") {
                        sele += " and /" + model;
                    }
                    return sele;
                } else {
                    return "NONE";
                }
            }
        },
        "ball+stick": {
            multipleBond: "symmetric",
            colorScheme: getLigandColorScheme,
            quality: function() {
                return "medium";
            },
            scale: 2.5,
            aspectRatio: 1.2,
            sele: function() {
                var sele = interaction ? interaction : "NONE";
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                if (model !== "all") {
                    sele += " and /" + model;
                }
                return sele;
            }
        },
        "label": {
            color: "#333333",
            zOffset: 2.0,
            attachment: "middle-center",
            showBackground: true,
            backgroundColor: "white",
            backgroundOpacity: 0.5,
            scale: 0.6,
            sele: function() {
                if (interaction && structureComponent) {
                    var s = structureComponent.structure;
                    var withinSele = new NGL.Selection(interaction + (model !== "all" ? " and /" + model : ""));
                    var as = s.getAtomSetWithinSelection(withinSele, 5);
                    if (model !== "all") {
                        as.intersection(s.getAtomSet(new NGL.Selection("/" + model)));
                    }
                    var asGroup = s.getAtomSetWithinGroup(as);
                    var sele = asGroup.toSeleString();
                    if (hydrogenVisibility === false) {
                        sele += " and not hydrogen";
                    }
                    if (ionVisibility === false) {
                        sele += " and not ion";
                    }
                    if (waterVisibility === false) {
                        sele += " and not water";
                    }
                    if (model !== "all") {
                        sele += " and /" + model;
                    }
                    sele += " and ( ( protein and .CA ) or ( nucleic and .P ) )";
                    return sele;
                } else {
                    return "NONE";
                }
            }
        },
        "surface": {
            colorScheme: getLigandColorScheme,
            opacity: 0.5,
            side: "front",
            surfaceType: "sas",
            probeRadius: 1.4,
            useWorker: false,
            scaleFactor: 2.5,
            sele: function() {
                var sele = interaction ? interaction : "NONE";
                if (hydrogenVisibility === false) {
                    sele += " and not hydrogen";
                }
                if (model !== "all") {
                    sele += " and /" + model;
                }
                return sele;
            }
        },
    };

    function evalParam(paramValue) {
        if (typeof paramValue === "function") {
            return paramValue();
        } else {
            return paramValue;
        }
    }

    function getPolymerParam(reprName, paramName) {
        return evalParam(polymerReprDefs[reprName][paramName]);
    }

    function getLigandParam(reprName, paramName) {
        return evalParam(ligandReprDefs[reprName][paramName]);
    }

    function getInteractionParam(reprName, paramName) {
        return evalParam(interactionReprDefs[reprName][paramName]);
    }

    function makeRepresentations(reprDefs, reprDict) {
        if (!structureComponent) return;
        for (var reprName in reprDefs) {
            var reprParams = reprDefs[reprName];
            var params = {
                lazy: true,
                visible: false,
                assembly: assembly
            };
            for (var paramName in reprParams) {
                params[paramName] = evalParam(reprParams[paramName]);
            }
            if (reprDict[reprName]) {
                reprDict[reprName].dispose();
            }
            reprDict[reprName] = structureComponent.addRepresentation(reprName, params);
        }
    }

    function makeCounts() {
        var structure = structureComponent.structure;
        var _assembly = structure.biomolDict[assembly];
        if (_assembly) {
            atomCount = _assembly.getAtomCount(structure)
            if (model !== "all") {
                atomCount /= structure.modelStore.count;
            }
            instanceCount = _assembly.getInstanceCount();
        } else {
            if (model === "all") {
                atomCount = structure.atomStore.count;
            } else {
                atomCount = structure.getModelProxy(0).atomCount;
            }
            instanceCount = 1;
        }
        if (typeof window.orientation !== 'undefined') {
            atomCount *= 4;
        }
        isBackboneOnly = structure.atomStore.count / structure.residueStore.count < 2;
        if (isBackboneOnly) {
            atomCount *= 10;
        }
    }

    function setStructure(comp) {
        console.log('setStructure');
        structureComponent = comp;
        initStructure();
        signals.structureChanged.dispatch(structureComponent);
    }

    function initStructure() {
        console.log('initStructure');
        if (!structureComponent) return;
        var s = structureComponent.structure;
        spatialHash = new NGL.SpatialHash(s.atomStore, s.boundingBox);
        makeCounts();
        if (getStyleNames(true)[style] === undefined) {
            style = getDefaultStyle();
        }
        if (getSymmetryNames()[symmetry] === undefined) {
            symmetry = 0;
        }
        if (getSymmetryInfo() === undefined) {
            loadSymmetryData(assembly).then(updateSymmetry);
        }
        setSymmetry(symmetry);
        makeRepresentations(polymerReprDefs, polymerReprDict);
        makeRepresentations(ligandReprDefs, ligandReprDict);
        makeRepresentations(interactionReprDefs, interactionReprDict);
        structureComponent.setDefaultAssembly(assembly);
        axesRepr = structureComponent.addRepresentation("axes", {
            visible: false
        });
        setStyle(style);
        setLigandStyle(ligandStyle);
        setColorScheme(colorScheme);
        if (assembly === "UNITCELL") {
            polymerReprDict["unitcell"].setVisibility(true);
        } else if (assembly === "SUPERCELL") {
            setLigandStyle("");
        } else {
            axesRepr.repr.align();
        }
        setClashVisibility();
        centerView();
        setInteraction(interaction);
    }
    

    // rl - added setPdbid because loadPdbid is not called
    function setPdbid(_pdbid, _assembly, _reduced) {
        console.log('setPdbid');
        pdbid = _pdbid;
        reduced = _reduced;
        symmetryData = {};
        if (symmetryBuffer) {
            symmetryBuffer.dispose();
            symmetryBuffer = undefined;
        }
        validationData = undefined;
        if (clashBuffer) {
            clashBuffer.dispose();
            clashBuffer = undefined;
        }
        if (_assembly !== undefined) {
            assembly = _assembly;
        }
    }


    function loadSymmetryData(_assembly) {
        tasks.increment();
        return retrieveSymmetryData(pdbid, _assembly).then(function(data) {
            tasks.decrement();
            if (!data || !data.symmetries || !data.symmetries.length) {
                symmetry = -1;
            } else {
                data.symmetries = data.symmetries.filter(function(sym) {
                    return sym.pointGroup !== "C1";
                });
                if (!data.symmetries.length) {
                    symmetry = -1;
                }
            }
            symmetryData[_assembly] = data;
            signals.symmetryDataLoaded.dispatch(data);
        }).catch(function(e) {
            tasks.decrement();
            console.error(e);
        });
    }

    function loadValidationData() {
        validationData = undefined;
        validationDataLoading = true;
        tasks.increment();
        return retrieveValidationData(pdbid).then(function(data) {
            tasks.decrement();
            validationData = data;
            validationDataLoading = false;
            signals.validationDataLoaded.dispatch(data);
        });
    }

    function updateSelections() {
        var name;
        for (name in polymerReprDict) {
            polymerReprDict[name].setSelection(getPolymerParam(name, "sele"));
        }
        for (name in ligandReprDict) {
            ligandReprDict[name].setSelection(getLigandParam(name, "sele"));
        }
        for (name in interactionReprDict) {
            interactionReprDict[name].setSelection(getInteractionParam(name, "sele"));
        }
    }

    function setStyle(value) {
        style = value;
        for (var name in polymerReprDict) {
            if (name === "unitcell" && assembly === "UNITCELL") {
                polymerReprDict["unitcell"].setVisibility(true);
            } else if (name === "ball+stick" && clashVisibility) {
                polymerReprDict["ball+stick"].setVisibility(true);
            } else if (name === "base" && style === "cartoon") {
                polymerReprDict["base"].setVisibility(true);
            } else {
                polymerReprDict[name].setVisibility(name === style);
            }
        }
        signals.styleChanged.dispatch(style);
    }

    function setLigandStyle(value) {
        ligandStyle = value;
        for (var name in ligandReprDict) {
            ligandReprDict[name].setVisibility(name === ligandStyle);
        }
        signals.ligandStyleChanged.dispatch(ligandStyle);
    }

    function setColorScheme(value) {
        if (value !== undefined) colorScheme = value;
        if (["fit", "geo"].includes(colorScheme) && validationData === undefined) {
            if (validationDataLoading) {
                signals.validationDataLoaded.addOnce(setColorScheme);
            } else {
                loadValidationData().then(setColorScheme);
            }
        } else {
            for (var name in polymerReprDict) {
                polymerReprDict[name].setParameters({
                    colorScheme: getPolymerParam(name, "colorScheme"),
                    colorScale: getPolymerParam(name, "colorScale")
                });
            }
            for (var name in ligandReprDict) {
                ligandReprDict[name].setParameters({
                    colorScheme: getLigandParam(name, "colorScheme"),
                    colorScale: getLigandParam(name, "colorScale")
                });
            }
            for (var name in interactionReprDict) {
                interactionReprDict[name].setParameters({
                    colorScheme: getInteractionParam(name, "colorScheme"),
                    colorScale: getInteractionParam(name, "colorScale")
                });
            }
        }
        signals.colorSchemeChanged.dispatch(colorScheme);
    }

    function setModel(value) {
        if (value !== "all") value = parseInt(value);
        if (value !== model && (model === "all" || value === "all")) {
            model = value;
            initStructure();
        } else {
            model = value;
            updateSelections();
        }
        if (clashBuffer) {
            clashBuffer.dispose();
            clashBuffer = undefined;
        }
        setClashVisibility();
        signals.modelChanged.dispatch(model);
    }

    function setHydrogenVisibility(value) {
        hydrogenVisibility = value;
        updateSelections();
        signals.hydrogenVisibilityChanged.dispatch(hydrogenVisibility);
    }

    function setIonVisibility(value) {
        ionVisibility = value;
        updateSelections();
        signals.ionVisibilityChanged.dispatch(ionVisibility);
    }

    function setWaterVisibility(value) {
        waterVisibility = value;
        updateSelections();
        signals.waterVisibilityChanged.dispatch(waterVisibility);
    }

    function setClashVisibility(value) {
        if (value !== undefined) clashVisibility = value;
        if (clashVisibility) {
            if (validationData === undefined) {
                if (validationDataLoading) {
                    signals.validationDataLoaded.addOnce(setClashVisibility);
                } else {
                    loadValidationData().then(setClashVisibility);
                }
            } else {
                if (clashBuffer) {
                    clashBuffer.setVisibility(true);
                } else {
                    var sele = (model === "all" ? undefined : "/" + model);
                    clashBuffer = new ClashBuffer(validationData.clashDict, structureComponent.structure, {
                        sele: sele
                    });
                    clashBuffer.attach(structureComponent);
                }
                updateSelections();
                polymerReprDict["ball+stick"].setVisibility(true);
            }
        } else {
            if (clashBuffer) {
                clashBuffer.setVisibility(false);
            }
            polymerReprDict["ball+stick"].setVisibility(false);
        }
        signals.clashVisibilityChanged.dispatch(clashVisibility);
    }

    function setQuality(value) {
        quality = value;
        polymerReprDict["surface"].setParameters({
            scaleFactor: getPolymerParam("surface", "scaleFactor")
        });
        polymerReprDict["cartoon"].setParameters({
            subdiv: getPolymerParam("cartoon", "subdiv"),
            radialSegments: getPolymerParam("cartoon", "radialSegments")
        });
        polymerReprDict["backbone"].setParameters({
            lineOnly: getPolymerParam("backbone", "lineOnly")
        });
        polymerReprDict["spacefill"].setParameters({
            sphereDetail: getPolymerParam("spacefill", "sphereDetail"),
        });
        signals.qualityChanged.dispatch(quality);
    }

    function getDefaultStyle() {
        if (atomCount < 200000) {
            return "cartoon";
        } else {
            return "surface";
        }
    }

    function getColorScale() {
        if (colorScheme === "hydrophobicity") {
            return "RdYlGn";
        } else if (colorScheme === "bfactor") {
            return "OrRd";
        } else {
            return "RdYlBu";
        }
    }

    function getColorScheme() {
        if (colorScheme === "fit") {
            return validationData ? validationData.fitScheme : "chainname";
        } else if (colorScheme === "geo") {
            return validationData ? validationData.geoScheme : "chainname";
        } else {
            return colorScheme;
        }
    }

    function getLigandColorScheme() {
        if (colorScheme === "bfactor") {
            return "bfactor";
        } else if (colorScheme === "fit") {
            return validationData ? validationData.fitScheme : "element";
        } else if (colorScheme === "geo") {
            return validationData ? validationData.geoScheme : "chainname";
        } else {
            return "element";
        }
    }

    function setAssembly(value) {
        assembly = value;
        initStructure();
        signals.assemblyChanged.dispatch(assembly);
    }

    function setInteraction(value) {
        interaction = value;
        if (structureComponent) {
            var s = structureComponent.structure;
            var available = Object.keys(getInteractionNames());
            if (interaction && !available.includes(interaction)) {
                var residues = [];
                s.eachResidue(function(rp) {
                    var sele = "";
                    if (rp.resno !== undefined) sele += rp.resno;
                    if (rp.inscode) sele += "^" + rp.inscode;
                    if (rp.chain) sele += ":" + rp.chainname;
                    residues.push(sele);
                }, new NGL.Selection(interaction + " and ( not polymer or not ( protein or nucleic ) )"));
                if (available.includes(residues[0])) {
                    interaction = residues[0]
                } else {
                    interaction = "";
                }
            }
            for (var name in interactionReprDict) {
                var sele = getInteractionParam(name, "sele");
                if (interaction) {
                    interactionReprDict[name].setSelection(sele);
                    interactionReprDict[name].setVisibility(true);
                } else {
                    interactionReprDict[name].setVisibility(false);
                    interactionReprDict[name].setSelection(sele);
                }
            }
            if (interaction) {
                structureComponent.centerView(true, interaction);
                stage.viewer.zoom(5);
                var sceneSize = stage.viewer.boundingBox.size().length();
                var interactionSize = sceneSize;
                if (structureComponent && interaction) {
                    var bb = s.getBoundingBox(new NGL.Selection(interaction));
                    interactionSize = Math.max(10, Math.min(sceneSize, bb.size().length() - 5));
                }
                if (ligandStyle) {
                    setLigandStyle("");
                }
                setFocus(Math.min(95, 100 - ((interactionSize / sceneSize) * 100)));
            } else {
                setFocus(0);
                if (!style) {
                    setStyle(getDefaultStyle());
                }
                if (!ligandStyle) {
                    setLigandStyle("ball+stick");
                }
                centerView();
            }
        }
        signals.interactionChanged.dispatch(interaction);
    }

    function setSpin(value) {
        spin = value;
        if (spin === true) {
            stage.setSpin([0, 1, 0], 0.005);
        } else if (value === false) {
            stage.setSpin(null, null);
        }
    }

    // Add in Stage as Argument
    function setSpinStage(stage, value) {
        spin = value;
        if (spin === true) {
            stage.setSpin([0, 1, 0], 0.005);
        } else if (value === false) {
            stage.setSpin(null, null);
        }
    }

    function setFocus(value) {
        focus = parseInt(value);
        stage.setParameters({
            clipNear: focus / 2,
            fogFar: 100 - (focus / 2),
            fogNear: 100 - (focus / 2) - (focus / 20)
        });
        signals.focusChanged.dispatch(focus);
    }

    function getSymmetryInfo() {
        var data = symmetryData[assembly];
        if (data && data.nrSymmetries) {
            return data.symmetries[symmetry];
        } else {
            return undefined;
        }
    }

    function updateSymmetry() {
        if (symmetryBuffer) {
            symmetryBuffer.dispose();
            symmetryBuffer = undefined;
        }
        var data = getSymmetryInfo();
        if (data && data.symmetryAxes) {
            symmetryBuffer = new SymmetryBuffer(data.symmetryAxes, {});
            symmetryBuffer.attach(structureComponent);
        }
        if (data && data.rotation && data.center) {
            var r = data.rotation;
            var m3 = new NGL.Matrix3().set(parseFloat(r.m00), parseFloat(r.m01), parseFloat(r.m02), parseFloat(r.m10), parseFloat(r.m11), parseFloat(r.m12), parseFloat(r.m20), parseFloat(r.m21), parseFloat(r.m22));
            var c = new NGL.Vector3().copy(data.center);
            var v1 = new NGL.Vector3(parseFloat(r.m00), parseFloat(r.m01), parseFloat(r.m02));
            var v2 = new NGL.Vector3(parseFloat(r.m10), parseFloat(r.m11), parseFloat(r.m12));
            var v3 = new NGL.Vector3(parseFloat(r.m20), parseFloat(r.m21), parseFloat(r.m22));
            stage.viewer.alignView(v3, v1, c, false);
            stage.viewer.centerView(true);
        }
    }

    function setSymmetry(value) {
        symmetry = parseInt(value) || 0;
        var data = getSymmetryInfo();
        if (symmetry === -1) {
            if (symmetryBuffer) {
                symmetryBuffer.dispose();
                symmetryBuffer = undefined;
            }
        } else if (data === undefined) {
            if (symmetryBuffer) {
                symmetryBuffer.dispose();
                symmetryBuffer = undefined;
            }
            loadSymmetryData(assembly).then(updateSymmetry);
        } else {
            updateSymmetry();
            setFocus(0);
            if (!style) {
                setStyle(getDefaultStyle());
            }
        }
        signals.symmetryChanged.dispatch(symmetry);
    }

    function centerView() {
        stage.centerView();
    }

    function downloadScreenshot() {
        stage.makeImage({
            factor: 2,
            antialias: true,
            trim: false,
            transparent: false
        }).then(function(blob) {
            NGL.download(blob, pdbid + "_screenshot.png");
        });
    }

    function toggleFullscreen(element) {
        stage.toggleFullscreen(element);
    }

    function getStyleNames(recommended) {
        var styleDict = {
            "": "None",
            backbone: "Backbone",
            surface: "Surface",
        };
        if (recommended) {
            if (atomCount < 100000) {
                styleDict["cartoon"] = "Cartoon";
            }
            if (atomCount < 80000) {
                styleDict["spacefill"] = "Spacefill";
            }
            if (atomCount < 80000) {
                styleDict["licorice"] = "Licorice";
            }
        } else {
            styleDict["cartoon"] = "Cartoon";
            if (!isBackboneOnly) {
                styleDict["spacefill"] = "Spacefill";
                styleDict["licorice"] = "Licorice";
            }
        }
        return styleDict;
    }

    function getLigandStyleNames() {
        return {
            "": "None",
            "ball+stick": "Ball & Stick",
            spacefill: "Spacefill"
        };
    }

    function getModelNames() {
        var modelDict = {};
        if (structureComponent) {
            var modelStore = structureComponent.structure.modelStore;
            if (modelStore.count > 1) {
                modelDict["all"] = "All Models";
            }
            for (var i = 0; i < modelStore.count; ++i) {
                modelDict[i] = "Model " + (i + 1);
            }
        }
        return modelDict;
    }

    function getAssemblyNames() {
        var assemblyDict = {};
        if (structureComponent) {
            var structure = structureComponent.structure;
            var biomolDict = structure.biomolDict;
            if (!structure.unitcell && Object.keys(biomolDict).length === 1 && biomolDict["BU1"] && biomolDict["BU1"].isIdentity(structure)) {
                assemblyDict["BU1"] = "Full Structure";
            } else {
                assemblyDict["__AU"] = (structure.unitcell ? "Asymmetric Unit" : "Full Structure");
                for (var name in biomolDict) {
                    if (name === "UNITCELL") {
                        assemblyDict[name] = "Unitcell";
                    } else if (name === "SUPERCELL") {
                        assemblyDict[name] = "Supercell";
                    } else if (name.substr(0, 2) === "BU") {
                        assemblyDict[name] = "Bioassembly " + name.substr(2);
                    } else {
                        assemblyDict[name] = name;
                    }
                }
            }
        }
        return assemblyDict;
    }

    function getColorSchemeNames() {
        var schemeDict = {
            chainname: "By Chain",
            residueindex: "Rainbow",
            element: "By Element",
            bfactor: "By B-factor",
            sstruc: "By Secondary Structure",
            hydrophobicity: "By Hydrophobicity",
            fit: "By Density Fit",
            geo: "By Geometry Quality"
        };
        if (structureComponent) {
            var methods = structureComponent.structure.header.experimentalMethods;
            if (methods && !methods.includes("X-RAY DIFFRACTION") && !methods.includes("ELECTRON CRYSTALLOGRAPHY") && !methods.includes("NEUTRON DIFFRACTION")) {
                delete schemeDict.bfactor;
                delete schemeDict.fit;
            }
        }
        if (!hasStructureFactors) {
            delete schemeDict.fit;
        }
        return schemeDict;
    }

    function getSymmetryNames() {
        var symmetryDict = {
            "-1": "None"
        };
        var data = symmetryData[assembly];
        if (data && data.symmetries) {
            for (var i = 0; i < data.symmetries.length; ++i) {
                var sym = data.symmetries[i];
                var type = sym.local ? "local" : "global";
                if (sym.pseudoSymmetric) {
                    type += ", pseudo";
                }
                symmetryDict[i] = sym.pointGroup + " (" + type + ")";
            }
        }
        return symmetryDict;
    }

    function getInteractionNames() {
        var interactionDict = {
            "": "None"
        };
        if (structureComponent) {
            var s = structureComponent.structure;
            var ligandSele = "( not polymer or not ( protein or nucleic ) ) and not ( water or ACE or NH2 )";
            var _assembly = s.biomolDict[assembly];
            if (_assembly) {
                ligandSele += " and (" + _assembly.getSelection().string + ")";
            }
            var ligandSelection = new NGL.Selection();
            s.eachResidue(function(rp) {
                if (rp.isWater()) return;
                var sele = "";
                if (rp.resno !== undefined) sele += rp.resno;
                if (rp.inscode) sele += "^" + rp.inscode;
                if (rp.chain) sele += ":" + rp.chainname;
                var name = (rp.resname ? "[" + rp.resname + "]" : "") + sele;
                interactionDict[sele] = name;
            }, new NGL.Selection(ligandSele));
        }
        return interactionDict;
    }
    // functions
    this.signals = signals;
    //this.loadPdbid = loadPdbid; // rl - not used
    this.centerView = centerView;
    this.downloadScreenshot = downloadScreenshot;
    this.toggleFullscreen = toggleFullscreen;
    this.setStyle = setStyle;
    this.setLigandStyle = setLigandStyle;
    this.setModel = setModel;
    this.setHydrogenVisibility = setHydrogenVisibility;
    this.setIonVisibility = setIonVisibility;
    this.setWaterVisibility = setWaterVisibility;
    this.setClashVisibility = setClashVisibility;
    this.setQuality = setQuality;
    this.setAssembly = setAssembly;
    this.setColorScheme = setColorScheme;
    this.setSpin = setSpin;
    this.setSpinStage = setSpinStage;
    this.setSymmetry = setSymmetry;
    this.setInteraction = setInteraction;
    this.setFocus = setFocus;
    this.getStyleNames = getStyleNames;
    this.getLigandStyleNames = getLigandStyleNames;
    this.getModelNames = getModelNames;
    this.getAssemblyNames = getAssemblyNames;
    this.getColorSchemeNames = getColorSchemeNames;
    this.getSymmetryNames = getSymmetryNames;
    this.getInteractionNames = getInteractionNames;
    // rl - added these
    this.setStage = setStage;
    this.setStructure = setStructure;
    this.setPdbid = setPdbid;
    // vars
    this.getStyle = function() {
        return style;
    };
    this.getLigandStyle = function() {
        return ligandStyle;
    };
    this.getModel = function() {
        return model;
    };
    this.getHydrogenVisibility = function() {
        return hydrogenVisibility;
    };
    this.getIonVisibility = function() {
        return ionVisibility;
    };
    this.getWaterVisibility = function() {
        return waterVisibility;
    };
    this.getClashVisibility = function() {
        return clashVisibility;
    };
    this.getQuality = function() {
        return quality;
    };
    this.getAssembly = function() {
        return assembly;
    };
    this.getSymmetry = function() {
        return symmetry;
    };
    this.getInteraction = function() {
        return interaction;
    };
    this.getColorScheme = function() {
        return colorScheme;
    };
    this.getSpin = function() {
        return spin;
    };
    this.getFocus = function() {
        return focus;
    };
    // rl added
    this.getPdbid = function() { return pdbid; };
    this.getAtomCount = function() { return atomCount; };
};

// END NglController ===================================================================================================

// helper functions fro NglController ==================================================================================
function xhrPromise(url, responseType) {
    return new Promise(function(resolve, reject) {
        var xhr = new XMLHttpRequest();
        xhr.addEventListener("load", function() {
            if (xhr.status === 200 || xhr.status === 304) {
                resolve(xhr.response);
            } else {
                reject(xhr.status);
            }
        }, true);
        xhr.addEventListener("error", function(event) {
            reject(event);
        }, true);
        xhr.responseType = responseType;
        xhr.open("GET", url);
        xhr.send();
    });
}

function retrieveSymmetryData(pdbid, bioassembly) {
    var bionumber = bioassembly2bionumber(bioassembly);
    var basePath = rcsbUrl + "/pdb/json/symmetryOrientation";
    var url = basePath + "?pdbID=" + pdbid + "&bioassembly=" + bionumber;
    return xhrPromise(url, "json");
}

function retrieveValidationData(pdbid) {
    pdbid = pdbid.toLowerCase();
    var baseUrl = "http://ftp.wwpdb.org/pub/pdb/validation_reports/";
    var filename = pdbid + "_validation.xml.gz";
    var url = baseUrl + pdbid.substr(1, 2) + "/" + pdbid + "/" + filename;
    return NGL.autoLoad(url, {
        ext: "xml",
        useDomParser: true,
        compressed: false
    }).then(function(xml) {
        return new ValidationData(xml);
    }).catch(function(e) {
        var xml = document.implementation.createDocument("http://wwpdb.org/validation/schema/wwpdb_validation_v002.xsd", "wwPDB-validation-information");
        return new ValidationData({
            data: xml
        });
    });
}

function bionumber2bioassembly(bionumber) {
    if (!bionumber || bionumber === "asym") {
        return "__AU";
    } else {
        return "BU" + bionumber;
    }
}

function bioassembly2bionumber(bioassembly) {
    if (!bioassembly || bioassembly === "__AU") {
        return "asym";
    } else {
        return bioassembly.substr(2);
    }
}

var SymmetryBuffer = function(axes, params) {
    var p = Object.assign({}, params);
    var c = new NGL.Color(p.color || "lime");
    var radius = p.radius || 0.5;
    var shape = new NGL.Shape("symmetry", {
        disableImpostor: false,
        openEnded: true
    });
    axes.forEach(function(ax) {
        shape.addSphere(ax.start, c, radius);
        shape.addSphere(ax.end, c, radius);
        shape.addCylinder(ax.start, ax.end, c, radius);
    });
    this.attach = function(component) {
        shapeRepr = component.addBufferRepresentation(shape.getBufferList());
    };
    this.dispose = function() {
        if (shapeRepr) shapeRepr.dispose();
    };
};

var ValidationData = function(xml) {
    var rsrzDict = {};
    var rsccDict = {};
    var clashDict = {};
    var geoDict = {};
    var geoAtomDict = {};

    function getSele(a, atomname, useAltcode) {
        var icode = a.icode.value;
        var chain = a.chain.value;
        var altcode = a.altcode.value;
        var sele = a.resnum.value;
        if (icode.trim()) sele += "^" + icode;
        if (chain.trim()) sele += ":" + chain;
        if (atomname) sele += "." + atomname;
        if (useAltcode && altcode.trim()) sele += "%" + altcode;
        sele += "/" + (parseInt(a.model.value) - 1);
        return sele;
    }

    function setBitDict(dict, key, bit) {
        if (dict[key] === undefined) {
            dict[key] = bit;
        } else {
            dict[key] |= bit;
        }
    }

    function countSetBits(i) {
        i = i - ((i >> 1) & 0x55555555);
        i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
        return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    }
    var groups = xml.data.getElementsByTagName("ModelledSubgroup");
    for (var i = 0, il = groups.length; i < il; ++i) {
        var g = groups[i];
        var ga = g.attributes;
        var sele = getSele(ga);
        var clashes = g.getElementsByTagName("clash");
        for (var j = 0, jl = clashes.length; j < jl; ++j) {
            var ca = clashes[j].attributes;
            var cid = parseInt(ca.cid.value);
            if (clashDict[cid] === undefined) {
                clashDict[cid] = {
                    mag: parseFloat(ca.clashmag.value),
                    dist: parseFloat(ca.dist.value),
                    sele1: getSele(ga, ca.atom.value, true),
                    atom1: ca.atom.value,
                    res1: sele
                };
            } else {
                clashDict[cid].sele2 = getSele(ga, ca.atom.value, true);
                clashDict[cid].atom2 = ca.atom.value;
                clashDict[cid].res2 = sele;
            }
        }
    }
    for (var i = 0, il = groups.length; i < il; ++i) {
        var g = groups[i];
        var ga = g.attributes;
        var sele = getSele(ga);
        if (ga.rsrz !== undefined) {
            rsrzDict[sele] = parseFloat(ga.rsrz.value);
        }
        if (ga.rscc !== undefined) {
            rsccDict[sele] = parseFloat(ga.rscc.value);
        }
        var isPolymer = ga.seq.value !== ".";
        var clashAtoms = [];
        var geoProblemCount = 0;
        var clashes = g.getElementsByTagName("clash");
        for (var j = 0, jl = clashes.length; j < jl; ++j) {
            var ca = clashes[j].attributes;
            var cid = parseInt(ca.cid.value);
            if (clashDict[cid] !== undefined) {
                var c = clashDict[cid];
                if (c.res1 === c.res2 || c.atom1 === undefined || c.atom2 === undefined || NGL.guessElement(c.atom1) === "H" || NGL.guessElement(c.atom2) === "H") {
                    delete clashDict[cid];
                } else {
                    clashAtoms.push(ca.atom.value);
                }
            }
        }
        if (isPolymer) {
            if (clashAtoms.length > 0) {
                geoProblemCount += 1;
            }
            var angleOutliers = g.getElementsByTagName("angle-outlier");
            if (angleOutliers.length > 0) {
                geoProblemCount += 1;
            }
            var bondOutliers = g.getElementsByTagName("bond-outlier");
            if (bondOutliers.length > 0) {
                geoProblemCount += 1;
            }
            var planeOutliers = g.getElementsByTagName("plane-outlier");
            if (planeOutliers.length > 0) {
                geoProblemCount += 1;
            }
            if (ga.rota !== undefined && ga.rota.value === "OUTLIER") {
                geoProblemCount += 1;
            }
            if (ga.rama !== undefined && ga.rama.value === "OUTLIER") {
                geoProblemCount += 1;
            }
            if (ga.RNApucker !== undefined && ga.RNApucker.value === "outlier") {
                geoProblemCount += 1;
            }
            if (geoProblemCount > 0) {
                geoDict[sele] = geoProblemCount;
            }
        } else {
            var mogBondOutliers = g.getElementsByTagName("mog-bond-outlier");
            var mogAngleOutliers = g.getElementsByTagName("mog-angle-outlier");
            if (mogBondOutliers.length > 0 || mogAngleOutliers.length > 0 || clashes.length > 0) {
                var atomDict = {};
                geoAtomDict[sele] = atomDict;
                for (var j = 0, jl = clashAtoms.length; j < jl; ++j) {
                    setBitDict(atomDict, clashAtoms[j], 1);
                }
                for (var j = 0, jl = mogBondOutliers.length; j < jl; ++j) {
                    var mbo = mogBondOutliers[j].attributes;
                    mbo.atoms.value.split(",").forEach(function(atomname) {
                        setBitDict(atomDict, atomname, 2);
                    });
                }
                for (var j = 0, jl = mogAngleOutliers.length; j < jl; ++j) {
                    var mao = mogAngleOutliers[j].attributes;
                    mao.atoms.value.split(",").forEach(function(atomname) {
                        setBitDict(atomDict, atomname, 4);
                    });
                }
            }
        }
    }
    var clashList = [];
    for (var k in clashDict) {
        var c = clashDict[k];
        clashList.push(c.res1, c.res2);
    }
    var clashSele = clashList.length ? clashList.join(" OR ") : "NONE";
    var fitScheme = NGL.ColorMakerRegistry.addScheme(function(params) {
        this.scale = "RdYlBu";
        this.domain = [2, 0];
        var rsrzScale = this.getScale();
        this.domain = [0.678, 1.0];
        var rsccScale = this.getScale();
        this.atomColor = function(atom) {
            var sele = atom.resno;
            if (atom.inscode) sele += "^" + atom.inscode;
            if (atom.chainname) sele += ":" + atom.chainname;
            sele += "/" + atom.modelIndex;
            var rsrz = rsrzDict[sele];
            if (rsrz !== undefined) {
                return rsrzScale(rsrz);
            }
            var rscc = rsccDict[sele];
            if (rscc !== undefined) {
                return rsccScale(rscc);
            }
            return 0x909090;
        };
    });
    var geoScheme = NGL.ColorMakerRegistry.addScheme(function(params) {
        this.atomColor = function(atom) {
            var geoProblemCount;
            var sele = atom.resno;
            if (atom.inscode) sele += "^" + atom.inscode;
            if (atom.chainname) sele += ":" + atom.chainname;
            sele += "/" + atom.modelIndex;
            geoAtom = geoAtomDict[sele];
            if (geoAtom !== undefined) {
                var atomProblems = geoAtom[atom.atomname] || 0;
                geoProblemCount = countSetBits(atomProblems);
            } else {
                geoProblemCount = geoDict[sele] || 0;
            }
            if (geoProblemCount === 0) {
                return 0x2166ac;
            } else if (geoProblemCount === 1) {
                return 0xfee08b;
            } else if (geoProblemCount === 2) {
                return 0xf46d43;
            } else if (geoProblemCount >= 3) {
                return 0xa50026;
            }
            return 0x909090;
        };
    });
    this.clashDict = clashDict;
    this.clashSele = clashSele;
    this.fitScheme = fitScheme;
    this.geoScheme = geoScheme;
};

var ClashBuffer = function(clashDict, structure, params) {
    var p = Object.assign({}, params);
    var color = new NGL.Color(p.color || "#f0027f");
    var sele = p.sele ? " and (" + p.sele + ")" : undefined;
    var shape = new NGL.Shape("clashes", {
        openEnded: false,
        disableImpostor: true
    });
    var s = structure;
    var ap1 = s.getAtomProxy();
    var ap2 = s.getAtomProxy();
    var vDir = new NGL.Vector3();
    var vPos1 = new NGL.Vector3();
    var vPos2 = new NGL.Vector3();
    for (var k in clashDict) {
        var c = clashDict[k];
        var sele1 = c.sele1;
        var sele2 = c.sele2;
        if (sele !== undefined) {
            sele1 += sele;
            sele2 += sele;
        }
        ap1.index = s.getAtomIndices(new NGL.Selection(sele1))[0];
        ap2.index = s.getAtomIndices(new NGL.Selection(sele2))[0];
        if (ap1.index === undefined || ap2.index === undefined) continue;
        vDir.subVectors(ap2, ap1).setLength(ap1.vdw);
        vPos1.copy(ap1).add(vDir);
        vDir.subVectors(ap1, ap2).setLength(ap2.vdw);
        vPos2.copy(ap2).add(vDir);
        var dHalf = ap1.distanceTo(ap2) / 2;
        var r1 = Math.sqrt(ap1.vdw * ap1.vdw - dHalf * dHalf);
        var r2 = Math.sqrt(ap2.vdw * ap2.vdw - dHalf * dHalf);
        shape.addCylinder(vPos1, vPos2, color, (r1 + r2) / 2);
    }
    this.attach = function(component) {
        shapeRepr = component.addBufferRepresentation(shape.getBufferList());
    };
    this.setVisibility = function(value) {
        if (shapeRepr) shapeRepr.setVisibility(value);
    };
    this.dispose = function() {
        if (shapeRepr) shapeRepr.dispose();
    };
};

export { NglController, bionumber2bioassembly }
