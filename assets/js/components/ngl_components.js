/**
 * Created by abradley on 01/03/2018.
 */
import { Stage } from 'ngl';
import React from 'react';
import { ShowMolLigs} from '../views/ngl_views';


export class NGLView extends React.Component {


    constructor(props) {
        super(props);
        // Create NGL Stage object
        this.div_id = "viewport";
        this.height = "600px";
        this.interval = 100;
        this.old_dict = {}
        this.showMol = this.showMol.bind(this);
    }


    componentDidMount(){
        this.stage = new Stage(this.div_id);
        // Handle window resizing
        window.addEventListener("resize", function (event) {
            this.stage.handleResize();
        }, false);
        this.showMol();
        setInterval(this.showMol,
            this.interval)
    }


    /**
     * Function to deal with the logic of showing molecules
     */
    showMol() {
        if (this.props.mol_dict && this.props.mol_dict != this.old_dict) {
            var showMolLigs = new ShowMolLigs(this.stage, this.props.mol_dict);
            showMolLigs.run();
        }
        // Now update the dicts
        this.old_dict = this.props.mol_dict;
        if(this.props.clear_all==true){
            this.stage.removeAllComponents();
            this.props.onMolCleared();
        }
    }

    render(){
        return <div style={{height: this.height}} id={this.div_id}></div>
    }
}