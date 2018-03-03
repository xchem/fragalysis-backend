console.log("ACTIONS file passed thru");

import {TOGGLE_COMPLEX, TOGGLE_SPHERE } from '../actonTypes';

export const toggleComplex = function (protein_id,mol_id,interactions,load) {
    return {
        type: TOGGLE_COMPLEX,
        load: load,
        protein: protein_id,
        mol: mol_id,
        interactions: interactions
    };
}

export const toggleSphere = function (coordinates,groupId,load) {
    return {
        type: TOGGLE_SPHERE,
        coordinates: coordinates,
        groupId: groupId,
        load: load,
    };
}

