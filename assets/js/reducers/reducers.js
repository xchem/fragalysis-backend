import * as actions from '../actions/actions.js';

function app(state, action) {
    console.log('REDUCERS FIRED OFF. OLD STATE');
    console.log(state);
    console.log('action.type=' + action.type);

    switch (action.type) {
        // Defined in initialState - but may be needed if we want to load a different structure
            
        case actions.TOGGLE_COMPLEX:
            return Object.assign({}, state, {
                visible: action.visible,
                protein: action.protein,
                mol: action.mol,
                interactions: action.interactions
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
