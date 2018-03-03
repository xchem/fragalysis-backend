/**
 * Created by abradley on 03/03/2018.
 */
import nglReducers from  './nglReducers'
import * as types from '../actions/actonTypes'
 

describe('NGL reducer', () => {
  it('should return the initial state', () => {
    expect(nglReducers(undefined, {})).toEqual(
      {
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
    )
  })
 
  it('should handle LOAD_OBJECT', () => {
    expect(nglReducers(undefined, {
        type: types.LOAD_OBJECT,
        loadObj: {name: "TESTOBJ", loadMe: "STRING"}
      })
    ).toEqual(
      {
          // Lists storing the information of what is in the viewer
          objectsToLoad: {"TESTOBJ": {name: "TESTOBJ", loadMe: "STRING"}},
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
    )
 
    expect(
      nglReducers(
          {
              // Lists storing the information of what is in the viewer
              objectsToLoad: {"TESTOBJ": {name: "TESTOBJ", loadMe: "STRING"}},
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
          },
        {
        type: types.LOAD_OBJECT,
        loadObj: {name: "TESTOBJ_TWO", loadMe: "STRING"}
        }
      )
    ).toEqual({
          // Lists storing the information of what is in the viewer
          objectsToLoad: {"TESTOBJ": {name: "TESTOBJ", loadMe: "STRING"},
              "TESTOBJ_TWO": {name: "TESTOBJ_TWO", loadMe: "STRING"}},
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
    )
  })
})