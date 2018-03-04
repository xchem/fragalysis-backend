/**
 * Created by abradley on 03/03/2018.
 */
import nglReducers,{INITIALSTATE} from  './nglReducers'
import * as types from '../actions/actonTypes'
 
const TEST_OBJ_ONE = {name: "TESTOBJ", loadMe: "STRING"}
const TEST_OBJ_TWO = {name: "TESTOBJ2", loadMe: "STRING_TWO"}

function getInitialState(){
    return {
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
}


// All the possible states
const TEST_LOAD_ONE_STATE = Object.assign({}, getInitialState(), {
    objectsToLoad: {"TESTOBJ": TEST_OBJ_ONE}
});
const TEST_LOADED_ONE_STATE = Object.assign({}, getInitialState(), {
    objectsInView: {"TESTOBJ": TEST_OBJ_ONE}
});
const TEST_DEL_ONE_STATE = Object.assign({}, getInitialState(), {
    objectsToDelete: {"TESTOBJ": TEST_OBJ_ONE},
    objectsInView: {"TESTOBJ": TEST_OBJ_ONE},
});



describe('NGL reducer', () => {
  it('should return the initial state', () => {
    expect(nglReducers(undefined, {})).toEqual(
        getInitialState()
    )
  })
 
  it('should handle LOAD_OBJECT', () => {
    expect(nglReducers(undefined, {
        type: types.LOAD_OBJECT,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        TEST_LOAD_ONE_STATE
    )
  })
    it('should handle LOAD_OBJECT_SUCCESS', () => {
    expect(nglReducers(TEST_LOAD_ONE_STATE, {
        type: types.LOAD_OBJECT_SUCCESS,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        TEST_LOADED_ONE_STATE
    )
  })
    it('should handle LOAD_OBJECT_FAILURE', () => {
    expect(nglReducers(TEST_LOAD_ONE_STATE, {
        type: types.LOAD_OBJECT_FAILURE,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        TEST_LOAD_ONE_STATE
    )
  })
    it('should handle DELETE_OBJECT', () => {
    expect(nglReducers(TEST_LOADED_ONE_STATE, {
        type: types.DELETE_OBJECT,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        TEST_DEL_ONE_STATE
    )
    })
    it('should handle DELETE_OBJECT_SUCCESS', () => {
    expect(nglReducers(TEST_DEL_ONE_STATE, {
        type: types.DELETE_OBJECT_SUCCESS,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        getInitialState()
    )
    })
    it('should handle DELETE_OBJECT_FAILURE', () => {
    expect(nglReducers(TEST_DEL_ONE_STATE, {
        type: types.DELETE_OBJECT_SUCCESS,
        loadObj:  TEST_OBJ_ONE
      })
    ).toEqual(
        TEST_DEL_ONE_STATE
    )
    })
})