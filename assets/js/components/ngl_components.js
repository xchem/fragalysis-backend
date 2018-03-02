/**
 * Created by abradley on 01/03/2018.
 */
import { Stage, concatStructures, Selection } from 'ngl';
import React from 'react';


export class NGLView extends React.Component {


    constructor(props) {
        super(props);
        // Create NGL Stage object
        this.div_id = "viewport";
        this.focus_var = 95;
        this.height = "600px";
        this.interval = 100;
        this.old_dict = {}
        this.mol_url = "/viewer/mol_from_pk/"
        this.prot_url = "/viewer/prot_from_pk/"
        this.showMol = this.showMol.bind(this);
        this.removeComponentByName = this.removeComponentByName.bind(this);
        this.molLigInteraction = this.molLigInteraction.bind(this);
        this.genMolLigInt = this.genMolLigInt.bind(this);
        this.getInputDict = this.getInputDict.bind(this);
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


    removeComponentByName(inputDict){
        var local_stage = this.stage
        var compName = inputDict["object_name"]
        local_stage.eachComponent(
            function(comp){
                if (comp.name==compName){
                    local_stage.removeComponent(comp)
                }
            })
    }

    molLigInteraction(input_dict){
        Promise.all([
                this.stage.loadFile(input_dict["prot_url"], {ext: "pdb"}),
                this.stage.loadFile(input_dict["mol_url"], {ext: "sdf"}),
                this.stage, this.focus_var, input_dict["object_name"]]
            ).then(genMolLigInt(ol));
    }

    genMolLigInt(ol){
            var cs = concatStructures(
                ol[4],
                ol[0].structure.getView(new Selection("not ligand")),
                ol[1].structure.getView(new Selection(""))
            )
            var stage = ol[2];
            var focus_var = ol[3];
            // Set the object name
            var comp = stage.addComponentFromObject(cs)
            comp.addRepresentation("cartoon")
            comp.addRepresentation("contact", {
                masterModelIndex: 0,
                weakHydrogenBond: true,
                maxHbondDonPlaneAngle: 35,
                sele: "/0 or /1"
            })
            comp.addRepresentation("licorice", {
                sele: "ligand and /1",
                multipleBond: "offset"
            })
            comp.addRepresentation("line", {
                sele: "/0"
            })
            comp.autoView("ligand");
            stage.setFocus(focus_var);
            NProgress.done();
    }

    getInputDict(mol_dict){
        var prot_id = this.props.mol_dict["prot_id"]
        var mol_id = this.props.mol_dict["mol_id"]
        var toggle = this.props.mol_dict["toggle"]
        var prot_url = this.prot_url + prot_id.toString() + "/"
        var mol_url = this.mol_url + mol_id.toString() + "/"
        var object_name = mol_id.toString()+"_mol"
        var inputDict ={"mol_url": mol_url, "prot_url": prot_url, "object_name":object_name, "toggle": toggle}
        return inputDict
    }

    showMol() {
        if (this.props.mol_dict && this.props.mol_dict != this.old_dict) {
            var inputDict = getInputDict(this.props.mol_dict);
            NProgress.start();
            if(inputDict["toggle"]==true) {
                this.molLigInteraction(inputDict)
            }
            else{
                this.removeComponentByName(inputDict);
                NProgress.done();
            }
        }
        // Now update the dicts
        this.old_dict = this.props.mol_dict;
    }

    render(){
        return <div style={{height: this.height}} id={this.div_id}></div>
    }
}