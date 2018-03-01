import React from 'react';
import ReactDOM from 'react-dom';
// TODO - fix this
import NGL from 'ngl';
import '../css/index.css';
import FillMe from 'generic_components';


function showPick (pickingProxy) {
    if (pickingProxy) {
        if (pickingProxy.object.name) {
            if (pickingProxy.object.name.startsWith("arrow") == true || pickingProxy.object.name.startsWith("cylinder") == true) {
                var smiles_str = pickingProxy.object.name.split(" ")[2].slice(1).slice(0, -1);

                for (var key in glob_data) {
                    if (key.startsWith(smiles_str)) {
                        var mols = glob_data[key];
                        var type = key.split("_")[2];
                        // Fill the relevant compound div
                        // We need compound divs defined
                        // TODO Define compound div in here
                    }
                }

            }
        }
    }
}


class NGLView extends React.Component {

    constructor(props) {
         // Create NGL Stage object
        var stage = new NGL.Stage(props.div_id);
        // Handle window resizing
        window.addEventListener("resize", function (event) {
            stage.handleResize();
        }, false);
        stage.mouseControls.add("clickPick-left",showPick);
    }

    render() {
    return (
        <div> </div>
    );
  }
}

class NGLSimpleFooter extends React.Component {
    render () {
        // The simple footer we have in most of the apps
        
    }
}



// The links between data and what is rendered
ReactDOM.render(<NGLView url="/path/to/single/view" />, document.getElementById('viewport'));