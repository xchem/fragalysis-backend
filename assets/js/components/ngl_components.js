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
        this.show_mol = this.show_mol.bind(this);
        this.removeComponentByName = this.removeComponentByName.bind(this);
        this.molLigInteraction = this.molLigInteraction.bind(this);
    }

    removeComponentByName(compName){
        var local_stage = this.stage
        local_stage.eachComponent(
            function(comp){
                if (comp.name==object_name){
                    local_stage.removeComponent(comp)
                }
            })
    }

    componentDidMount(){
        this.stage = new Stage(this.div_id);
        // Handle window resizing
        window.addEventListener("resize", function (event) {
            this.stage.handleResize();
        }, false);
        this.show_mol();
        setInterval(this.show_mol,
            this.interval)
    }
    molLigInteraction(prot_url,mol_url,object_name){
        Promise.all([
                this.stage.loadFile(prot_url, {ext: "pdb"}),
                this.stage.loadFile(mol_url, {ext: "sdf"}),
                this.stage, this.focus_var, object_name]
            ).then(function (ol) {
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
        });
    }


    show_mol() {
        if (this.props.mol_dict && this.props.mol_dict != this.old_dict) {
            const prot_id = this.props.mol_dict["prot_id"]
            const mol_id = this.props.mol_dict["mol_id"]
            const toggle = this.props.mol_dict["toggle"]
            const prot_url = this.prot_url + prot_id.toString() + "/"
            const mol_url = this.mol_url + mol_id.toString() + "/"
            const object_name = mol_id.toString()+"_mo"
            NProgress.start();
            if(toggle==true) {
                this.molLigInteraction(prot_url,mol_url,object_name)
            }
            else{
                this.remove_component_by_name(object_name);
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