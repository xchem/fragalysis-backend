import { createStore } from 'redux';
import { h, render } from 'preact';
import NGL from 'ngl';

// Subscribe > Action > Reducer > Component
import { Header, DisplayOptions } from './components/components.js';
import { app } from './reducers/reducers.js';
import { setAssemblyOptions } from "./actions/actions.js";
import { initRepr, setColor, setStyle, setSpin, setHydrogen, setWater } from "./ngl-script.js";
import { NglController } from './ngl-controller';
import { getQueryStringParameterByName, inspect } from "./util.js";

// Kick start the initialization process


export function init() {
    // Get the pdbId from the parameterString or set a default pdbId
    const pdbId = (getQueryStringParameterByName('pdbId') !== '') ? getQueryStringParameterByName('pdbId') : '4cup';

    // Set the NGL Stage, argument is the id of HTML element
    const stage = new NGL.Stage('viewport');
    console.log('********** NGL stage created **********');

    // NGL remove all component, load file
    stage.removeAllComponents();

    // Promise resolves to 'StructureComponent'
    stage.loadFile( "rcsb://" + pdbId, {
            defaultRepresentation: true
        } ).then( function( _structureComponent ) {
            // Utilizing the return promise
            const structureComponent = _structureComponent;
            console.log(structureComponent);

            // Function initRepr() - returns console logs
            initRepr(structureComponent, null);

            // Set Initial State for Redux STORE
            // TODO this could be combined with setAssemblyOptions
            const initialState = {
                pdbId: pdbId,
                structureTitle: structureComponent.structure.title,
                spin: false,
                hydrogen: true,
                water: false,
                assembly: structureComponent.defaultAssembly,
                color: 'rainbow',
                colorOptions: [
                    { value: "rainbow", label: "Rainbow" },
                    { value: "secondaryStructure", label: "SecondaryStructure" },
                    { value: "chain", label: "Chain" }
                ],
                style: 'cartoon',
                styleOptions: [
                    { value: "cartoon", label: "Cartoon" },
                    { value: "spheres", label: "Spheres" },
                    { value: "surface", label: "Surface" }
                ]
            };


            const nglController = new NglController( {
                spin: false,
            } );
            nglController.setStage(stage);
            console.log(nglController.getSpin());
            // nglController.setSpin(true);

            // Creates a REDUX store that holds the complete state tree of app
            // createStore(reducer, [preloadedState])

            // Go to reducers.js
            console.log('********* Store created with the initialState being 2nd variable **********');
            const store = createStore(app, initialState);
            console.log('********* Store creation DONE **********');


            // Dispatching the action: setAssemblyOptions
            store.dispatch(
                setAssemblyOptions(structureComponent.structure)
            );

            console.log('The assembly options for the specific PDB ID');
            console.log(store.getState().assemblyOptions);

            // Any change in STATE will trigger call to update NGL stage
            store.subscribe(() => {
                updateStageFromReduxStore(structureComponent, store, nglController);
            });


            // Render the UI in index.html
            initUi(store);
    });
}

// Initialize the UI, second argument is the id of the HTML structure
function initUi (store) {

    // Rendering React components to root DOM node
    const header = document.getElementById('header');
    const displayOptions = document.getElementById('displayOptions');

    // id, title stored into this.props in Header component
    render(
        <Header store={store}/>, header
    );
    // Using the store
    render(
        <DisplayOptions store={store} />, displayOptions
    );
}

// Function to update the NGL STAGE whenever REDUX state changes
function updateStageFromReduxStore(structureComponent, store, nglController) {

    console.log('UPDATING NGL STAGE WITH REDUX STORE. USING NEW STATE');
    console.log(store.getState());

    // NGL structure component, using functions in ngl-script.js
    if ( structureComponent && structureComponent.type==="structure" ) {

        // Set Assembly, Color, Style
        structureComponent.setDefaultAssembly( store.getState().assembly );
        setColor(store.getState().color);
        setStyle(store.getState().style);
        setSpin(store.getState().spin);
        console.log("+++++ from NGL controller");
        
        // NOT WORKING
        nglController.setSpinStage(true);
        console.log("+++++ end from NGL controller");
        setHydrogen(store.getState().hydrogen);
        setWater(store.getState().water);
    }
}
