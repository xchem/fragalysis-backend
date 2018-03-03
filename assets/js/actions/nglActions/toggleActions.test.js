/**
 * Created by abradley on 03/03/2018.
 */
import * as actions from './toggleActions'
import * as types from '../actonTypes'
 
describe('actions', () => {
  it('Should create an action to load a new complex', () => {
      const protein = 1
      const mol = 12
      const interactions = true
      const load = false
      // const load = true
      const expectedAction = {
          type: types.TOGGLE_COMPLEX,
          protein,
          mol,
          interactions,
          load
    }
    expect(actions.toggleComplex(protein,mol,interactions,load)).toEqual(expectedAction)
  })
})


describe('actions', () => {
  it('Should create an action to load a new complex', () => {
      const coordinates = [1.0,2.0,3.0]
      const groupId = 12
      const load = true
      // const load = true
      const expectedAction = {
          type: types.TOGGLE_SPHERE,
          coordinates,
          groupId,
          load
    }
    expect(actions.toggleSphere(coordinates, groupId, load)).toEqual(expectedAction)
  })
})

