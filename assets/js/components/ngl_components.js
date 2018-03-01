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
        this.show_mol = this.show_mol.bind(this);
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

    show_mol() {
        const MOL_URL = "/viewer/mol_from_pk/"
        const PROT_URL = "/viewer/prot_from_pk/"
        if (this.props.mol_dict && this.props.mol_dict != this.old_dict) {
            const prot_id = this.props.mol_dict["prot_id"]
            const mol_id = this.props.mol_dict["mol_id"]
            const toggle = this.props.mol_dict["toggle"]
            const prot_url = PROT_URL + prot_id.toString() + "/"
            const mol_url = MOL_URL + mol_id.toString() + "/"
            const object_name = mol_id.toString()
            NProgress.start();
            if(toggle==true) {
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
                    this.stage.setFocus(focus_var);
                });
            }
            else{
                this.stage.eachComponent(
                  function(comp){
                  if (comp.name==object_name){
                      this.stage.removeComponent(comp)
                  }
                })

            }
        }
        NProgress.done();
        this.old_dict = this.props.mol_dict;
    }
    render(){
        return <div style={{height: this.height}} id={this.div_id}></div>

    }
}