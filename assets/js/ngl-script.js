var structureComponent;
var cartoonRepr, sphereRepr, surfaceRepr;
var ligandRepr, ligandSphereRepr, ligandSurfaceRepr;
var currentColorScheme = 'rainbow';
var currentStyle = 'cartoon';
var currentSpin = false;
var atomCount, instanceCount;
var ligandSele = "";

function initRepr(_structureComponent, groupNames) {
    structureComponent = _structureComponent;
    if (structureComponent.structure.biomolDict["BU1"]) {
        var assembly = structureComponent.structure.biomolDict["BU1"];
        atomCount = assembly.getAtomCount(structureComponent.structure);
        instanceCount = assembly.getInstanceCount();
    } else {
        atomCount = structureComponent.structure.getModelProxy(0).atomCount;
        instanceCount = 1;
    }

    // feature mostly only available on phones/tablets
    if (typeof window.orientation !== 'undefined') {
        atomCount *= 4;  // bump atomCount which is used for deciding quality level
    }
    //
    ligandSele = 'none';
    if (groupNames) {
        groupNames = groupNames.split(/,\s*/);
        if (groupNames.length) {
            ligandSele = '[' + groupNames.join('] OR [') + ']';
        }
    }
    //
    initLigand();
    initSpheres();
    initCartoon();
    initSurface();

    // Set Style and Spin
    setStyle(currentStyle);
    setSpin(currentSpin);
    colorRainbow();
}

// init representations
function initLigand() {
    ligandRepr = structureComponent.addRepresentation('ball+stick', {
        color: 'element',
        scale: 3.0,
        aspectRatio: 1.3,
        sele: ligandSele,
        visible: false,
        lazy: true
    });
}
function initSpheres() {
    var quality = atomCount < 15000 ? "medium" : "low";
    sphereRepr = structureComponent.addRepresentation('spacefill', {
        colorScheme: currentColorScheme,
        colorScale: 'RdYlBu',
        quality: quality,
        sele: 'polymer',
        visible: false,
        lazy: true
    });
    //
    ligandSphereRepr = structureComponent.addRepresentation('spacefill', {
        color: 'element',
        sele: ligandSele,
        visible: false,
        lazy: true
    });
}
function initCartoon() {
    var quality = "low";
    if (atomCount < 15000) {
        quality = "high";
    } else if (atomCount < 70000) {
        quality = "medium";
    }
    cartoonRepr = structureComponent.addRepresentation('cartoon', {
        colorScheme: currentColorScheme,
        colorScale: 'RdYlBu',
        aspectRatio: 5,
        scale: 0.7,
        quality: quality,
        visible: false,
        lazy: true
    });
}
function initSurface() {
    surfaceRepr = structureComponent.addRepresentation('surface', {
        colorScheme: currentColorScheme,
        colorScale: 'RdYlBu',
        surfaceType: "sas",
        probeRadius: 1.4,
        scaleFactor: Math.min(1.5, Math.max(0.1, 20000 / atomCount)),
        sele: 'polymer',
        visible: false,
        lazy: true
    });
    ligandSurfaceRepr = structureComponent.addRepresentation('surface', {
        colorScheme: 'element',
        surfaceType: "sas",
        probeRadius: 1.4,
        scaleFactor: Math.min(1.5, Math.max(0.1, 20000 / atomCount)),
        sele: ligandSele,
        visible: false,
        lazy: true
    });
}

function styleSpheres() {
    sphereRepr.setVisibility(true);
    cartoonRepr.setVisibility(false);
    surfaceRepr.setVisibility(false);
    ligandRepr.setVisibility(false);
    ligandSurfaceRepr.setVisibility(false);
    ligandSphereRepr.setVisibility(true);
}
function styleCartoon() {
    sphereRepr.setVisibility(false);
    cartoonRepr.setVisibility(true);
    surfaceRepr.setVisibility(false);
    ligandRepr.setVisibility(true);
    ligandSurfaceRepr.setVisibility(false);
    ligandSphereRepr.setVisibility(false);
}
function styleSurface() {
    sphereRepr.setVisibility(false);
    cartoonRepr.setVisibility(false);
    surfaceRepr.setVisibility(true);
    ligandRepr.setVisibility(false);
    ligandSurfaceRepr.setVisibility(true);
    ligandSphereRepr.setVisibility(false);
}

// set color schemes
function colorRainbow() {
    var p = {colorScheme: 'residueindex', colorScale: 'RdYlBu'};
    cartoonRepr.setParameters(p);
    sphereRepr.setParameters(p);
    surfaceRepr.setParameters(p);
    currentColorScheme = 'residueindex';
    console.log('colorRainbow(): currentColorScheme=' + currentColorScheme);
}
function colorSecondaryStructure() {
    var p = {colorScheme: 'sstruc'};
    cartoonRepr.setParameters(p);
    sphereRepr.setParameters(p);
    surfaceRepr.setParameters(p);
    currentColorScheme = 'sstruc';
    console.log('colorSecondaryStructure(): currentColorScheme=' + currentColorScheme);
}
function colorChain() {
    var p = {colorScheme: 'chainindex'};
    cartoonRepr.setParameters(p);
    sphereRepr.setParameters(p);
    surfaceRepr.setParameters(p);
    currentColorScheme = 'chainindex';
    console.log('colorChain(): currentColorScheme=' + currentColorScheme);
}


// added for ngl-ui
function setColor(color) {
    console.log('setColor: color=' + color);
    switch (color) {
        case 'rainbow':
            colorRainbow();
            break;
        case 'secondaryStructure':
            colorSecondaryStructure();
            break;
        case 'chain':
            colorChain();
            break;
        default:
            colorRainbow();
    }
}

// Spin - NGL
function setSpin(spin) {
    if (spin === true) {
        console.log("Spin Set True: stage.setSpin([0, 1, 0], 0.005)");
    } else if (spin === false) {
        console.log("Spin Set False: stage.setSpin(null, null)");
    }
}

// Hydrogen
function setHydrogen(hydrogen) {
    if (hydrogen === true) {
        console.log("Hydrogen Set True: getHydrogenVisibility");
    } else if (hydrogen === false) {
        console.log("Hydrogen Set False");
    }
}

// Water
function setWater(water) {
    if (water === true) {
        console.log("Water Set True: getWaterVisibility");
    } else if (water === false) {
        console.log("Water Set False");
    }
}


// set styles
function setStyle(style) {
    console.log('setStyle: style=' + style);
    currentStyle = style;
    switch (style) {
        case 'cartoon':
            styleCartoon();
            break;
        case 'spheres':
            styleSpheres();
            break;
        case 'surface':
            styleSurface();
            break;
        default:
            colorCartoon();
    }
}

export {initRepr, setColor, setStyle, setSpin, setHydrogen, setWater}
