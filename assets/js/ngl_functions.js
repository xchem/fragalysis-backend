
class NGLFunctions {
    // Class of NGL Helper functions

    constructor(stage){
        this.stage = stage;
    }

    getJsonView() {
        // Get component (path) -> representation (parameters)
        var out_data = new Array()
        this.stage.eachComponent(function (component) {
            // TODO Sort this out for objects
            if (component.structure) {
                var component_data = {
                    "file_path": component.structure.path,
                    "repr": new Array(),
                };
                component.eachRepresentation(function (repr) {
                    var my_params = {"data": repr.getParameters()};
                    my_params["name"] = repr.name;
                    component_data["repr"].push(my_params);
                })
                out_data.push(component_data)
            }
        });
        return JSON.stringify({
            "components": out_data,
            "orientation": this.stage.viewerControls.getOrientation()
        })
    }

    setJsonView(input_json) {
        var tot_obj = JSON.parse(input_json);
        var obj = tot_obj["components"];
        var ori = tot_obj["orientation"];
        // Loop through dict
        for (var index = 0; index < obj.length; ++index) {
            // Generate the ey value pairs
            var propt = obj[index];
            var ext = "";
            if (propt["file_path"].includes("mol_from_pk")) {
                ext = "sdf";
            }
            else {
                ext = "pdb";
            }

            if (propt["file_path"].includes(":::")) {
                var f_p = propt["file_path"];
                file_path_1 = f_p.split(":::")[0]
                file_path_2 = f_p.split(":::")[1]
                var this_repr = propt["repr"];
                Promise.all([
                    this.stage.loadFile(file_path_1, {ext: "pdb"}),
                    this.stage.loadFile(file_path_2, {ext: "sdf"})
                ]).then(
                    function (ol) {
                        var cs = NGL.concatStructures("concat",
                            ol[0].structure.getView(new NGL.Selection("not ligand")),
                            ol[1].structure.getView(new NGL.Selection(""))
                        )
                        var comp = this.stage.addComponentFromObject(cs)
                        for (j = 0; j < this.length; ++j) {
                            comp.addRepresentation(this[j]["name"], this[j]["data"]);
                        }
                    }.bind(this_repr));
            }
            else {
                var this_repr = propt["repr"];
                this.stage.loadFile(propt["file_path"], {ext: ext})
                    .then(function (comp) {
                        for (j = 0; j < this.length; ++j) {
                            comp.addRepresentation(this[j]["name"], this[j]["data"]);
                        }
                    }.bind(this_repr));
            }
        }
        // Set the orientation
        var curr_orient = this.stage.viewerControls.getOrientation();
        for (i = 0; i < curr_orient["elements"].length; i++) {
            curr_orient["elements"][i] = ori["elements"][i];
        }
        this.stage.viewerControls.orient(curr_orient);
        this.stage.setFocus(focus_var);
    }

    post_view(title) {
        NProgress.start();
        var scene = getJsonView();
        $.ajax({
            type: "POST",
            url: "/viewer/post_view/",
            data: {
                "title": title,
                "scene": scene
            },
            success: function (response) {
                $("#clip_target").val(response);
                NProgress.done();
            }
        }).fail(function () {
            alert("error");
        });
    }

    get_view(pk) {
        $.ajax({
            type: "POST",
            url: "/viewer/get_view/" + pk.toString() + "/",
            success: function (response) {
                display_view(response)
            }
        }).fail(function () {
            alert("error");
        });
    }

    display_view(response) {
        var res = JSON.parse(response);
        var title = res["title"];
        setJsonView(res["scene"]);
    }

    display_arrow(start, end, colour, key) {
        var shape = new NGL.Shape(key.split("_")[0]);
        shape.addArrow(start, end, colour, 0.5);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer");
    }

    display_cylinder(start, end, colour, key) {
        var shape = new NGL.Shape(key.split("_")[0]);
        shape.addCylinder(start, end, colour, 0.8);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer");
    }
}