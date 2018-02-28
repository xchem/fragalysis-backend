import * as actions from '../actions/actions.js';

function app(state, action) {
    console.log('REDUCERS FIRED OFF. OLD STATE');
    console.log(state);
    console.log('action.type=' + action.type);

    switch (action.type) {
        // Defined in initialState - but may be needed if we want to load a different structure
        case actions.SET_STRUCTURE:
            return Object.assign({}, state, {
                pdbId: action.pdbId
            });

        case actions.SET_ASSEMBLY:
            return Object.assign({}, state, {
                assembly: action.assembly
            });
        case actions.SET_ASSEMBLY_OPTIONS:
            return Object.assign({}, state, {
                assemblyOptions: action.assemblyOptions
            });

        case actions.SET_COLOR:
            return Object.assign({}, state, {
                color: action.color
            });

        case actions.SET_STYLE:
            return Object.assign({}, state, {
                style: action.style
            });
        case actions.SET_SPIN:
            return Object.assign({}, state, {
                spin: action.spin
            });
        case actions.SET_WATER:
            return Object.assign({}, state, {
                water: action.water
            });
        case actions.SET_HYDROGEN:
            return Object.assign({}, state, {
                hydrogen: action.hydrogen
            });
        // Cases like: @@redux/INIT
        default:
            return state;
    }
}

export { app };
