/**
 * Created by abradley on 02/03/2018.
 */
class BaseNGLDisplay {

    constructor(stage, mol_dict) {
        this.stage = stage
        this.mol_dict = mol_dict
        this.focus_var = 95;
        this.mol_url = "/viewer/mol_from_pk/"
        this.prot_url = "/viewer/prot_from_pk/"
    }

    run(){
        var inputDict = this.getInputDict();
        NProgress.start();
        if(inputDict["toggle"]==true) {
            this.generatePromise(inputDict)
        }
        else{
            this.removeComponentByName(inputDict);
            NProgress.done();
        }
    }

    removeComponentByName(inputDict){
        var local_stage = this.stage
        var compName = inputDict["object_name"]
        local_stage.eachComponent(
            function(comp){
                if (comp.name==compName){
                    local_stage.removeComponent(comp)
                }
            }
        )
    }

}


class ShowMolLigs extends BaseNGLDisplay {

    generatePromise(input_dict){
        Promise.all([
                this.stage.loadFile(input_dict["prot_url"], {ext: "pdb"}),
                this.stage.loadFile(input_dict["mol_url"], {ext: "sdf"}),
                this.stage, this.focus_var, input_dict["object_name"]]
            ).then(ol => this.renderDisplay(ol));
    }

    renderDisplay(ol){
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

    getInputDict(){
        var prot_id = this.mol_dict["prot_id"]
        var mol_id = this.mol_dict["mol_id"]
        var toggle = this.mol_dict["toggle"]
        var prot_url = this.prot_url + prot_id.toString() + "/"
        var mol_url = this.mol_url + mol_id.toString() + "/"
        var object_name = mol_id.toString()+"_mol"
        return {"mol_url": mol_url, "prot_url": prot_url, "object_name":object_name, "toggle": toggle}
    }
}
export { ShowMolLigs }