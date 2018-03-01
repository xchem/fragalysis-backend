/**
 * Created by abradley on 01/03/2018.
 */
import { Stage, concatStructures, Selection } from 'ngl';
import React from 'react';


class NGLView extends React.Component {


    constructor(props) {
        super(props);
        // Create NGL Stage object
        this.div_id = "viewport";
        this.focus_var = 95;
        this.height = "600px";
        // this.stage.mouseControls.add("clickPick-left",showPick);
    }

    componentDidMount(){
        this.stage = new Stage(this.div_id);
        // Handle window resizing
        window.addEventListener("resize", function (event) {
            this.stage.handleResize();
        }, false);

    }

    show_mol() {
        NProgress.start();
        Promise.all([
            this.stage.loadFile(prot_url, {ext: "pdb"}),
            this.stage.loadFile(lig_data, {ext: "sdf"})]
        ).then(function (ol) {
            var cs = concatStructures(
                "concat",
                ol[0].structure.getView(new Selection("not ligand")),
                ol[1].structure.getView(new Selection(""))
            )
            cs.path = ol[0].structure.path + ":::" + ol[1].structure.path
            var comp = this.stage.addComponentFromObject(cs)
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
            NProgress.done();
            this.stage.setFocus(this.focus_var);
        })
    }
    render(){

        return <div style={{height: this.height}} id={this.div_id}></div>

    }
}