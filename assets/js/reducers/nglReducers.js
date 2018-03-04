import * as actions from '../actions/actonTypes'

const INITIALSTATE = {
      // Lists storing the information of what is in the viewer
      objectsToLoad: {},
      objectsToDelete: {},
      objectsInView: {},
      // Set the basic things about NGL
      visible: true,
      interactions: true,
      color: "blue",
      style: "xstick",
      spin: false,
      water: true,
      hydrogen: true,
}


export default function nglReducers(state = INITIALSTATE, action) {
    console.log('REDUCERS FIRED OFF. OLD STATE');
    console.log(state);
    console.log('action.type=' + action.type);

    switch (action.type) {
        // Defined in initialState - but may be needed if we want to load a different structure
            
        case actions.LOAD_OBJECT:
            // Append the input to objectsToLoad list
            var objectsToLoad = state.objectsToLoad
            objectsToLoad[action.loadObj.name]=action.loadObj
            return Object.assign({}, state, {
                objectsToLoad: objectsToLoad
            });
        // NEED TEST FOR EVERYTHING BELOW
        case actions.LOAD_OBJECT_SUCCESS:
            // Remove from objectsToLoad List
            var objectsToLoad = state.objectsToLoad
            delete objectsToLoad[action.loadObj.name]
            // Add to Objects in view list
            var objectsInView = state.objectsInView
            objectsInView[action.loadObj.name]=action.loadObj
            return Object.assign({}, state, {
                objectsInView: objectsInView,
                objectsToLoad: objectsToLoad
            });
        case actions.LOAD_OBJECT_FAILURE:
            // Don't change state - but can be used later to count number of failures
            console.log("Failed to load - "+action.loadObj.name)
            return Object.assign({}, state, {
            });
        case actions.DELETE_OBJECT:
            // Append the input to objectsToDelete list
            var objectsToDelete = state.objectsToDelete
            objectsToDelete[action.loadObj.name]=action.loadObj
            return Object.assign({}, state, {
            objectsToDelete: objectsToDelete
        });

        case actions.DELETE_OBJECT_SUCCESS:
            // Remove from objectsToDelete list
            var objectsToDelete = state.objectsToDelete
            delete objectsToDelete[action.loadObj.name]
            // Remove from ObjecsIn view list
            var objectsInView = state.objectsInView
            delete objectsInView[action.loadObj.name]
            return Object.assign({}, state, {
                objectsToDelete:objectsToDelete,
                objectsInView: objectsInView
            });
        case actions.DELETE_OBJECT_FAILURE:
            // Don't change state - but can be used later to count number of failures
            console.log("Failed to delete - "+action.loadObj.name)
            return Object.assign({}, state, {
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
