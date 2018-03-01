import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import $ from 'jquery';
import SVGInline from "react-svg-inline"
import { ListGroup, ListGroupItem, Tooltip, OverlayTrigger, Col, Row} from 'react-bootstrap';

import {toggleComplex} from './actions/actions'

import { createStore } from 'redux'
import { app }from './reducers/reducers'

let store = createStore(app);

// Basic config of the API
const BASE_API = "/v0.1/"
const TARGET_URL = BASE_API+"targets/"
const COMPOUNDS_URL = BASE_API+"compounds/"
const MOLECULE_URL = BASE_API+"molecules/"
const PROTEIN_URL = BASE_API+"proteins/"
const SVG_CMPD = '/viewer/img_from_cmpd_pk/'
const SVG_MOL = '/viewer/img_from_mol_pk/'



// Actions
/*
 * action types
 */
export const SELECT_MOL = 'SELECT_MOL'
export const TOGGLE_MOL = 'TOGGLE_MOL'
export const SET_VISIBILITY_FILTER = 'SET_VISIBILITY_FILTER'
/*
 * other constants
 */
export const VisibilityFilters = {
  SHOW_ALL: 'SHOW_ALL',
  SHOW_COMPLETED: 'SHOW_SELECTED',
  SHOW_ACTIVE: 'SHOW_ACTIVE'
}
export function selectMol(index) {
  return { type: SELECT_MOL, index }
}
 
export function toggleMol(index) {
  return { type: TOGGLE_MOL, index }
}
 
export function setVisibilityFilter(filter) {
  return { type: SET_VISIBILITY_FILTER, filter }
}

// Generic Classes
export class GenericList extends React.Component {

    constructor(props) {
    super(props);
        this.url = BASE_API
        this.interval = 1000
        this.loadFromServer = this.loadFromServer.bind(this);
        this.state = { data: [] };
  }

  loadFromServer() {
        $.ajax({
            url: this.url,
            datatype: 'json',
            cache: false,
            success: function (data) {
                this.setState({data: data})
            }.bind(this)
        })
    }

    componentDidMount() {
        this.loadFromServer();
        setInterval(this.loadFromServer,
            this.interval)
    }


    render() {
        if (this.state.data) {
            console.log(this.props.message)
            //
            return <ListGroup>
                {
                this.state.data.map((data, index) => (this.render_method(data,index)))
                 }
            </ListGroup>;
        }
        else {
            return (<FillMe />)
        }
    }

}

export class GenericView extends React.Component{


    constructor(props) {
    super(props);
        this.loadFromServer = this.loadFromServer.bind(this);
        this.handleClick = this.handleClick.bind(this);
        this.state = {isToggleOn: false,
            data: '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="110px" height="110px"><g>' +
        '<circle cx="50" cy="0" r="5" transform="translate(5 5)"/>' +
        '<circle cx="75" cy="6.6987298" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="93.3012702" cy="25" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="100" cy="50" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="93.3012702" cy="75" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="75" cy="93.3012702" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="50" cy="100" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="25" cy="93.3012702" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="6.6987298" cy="75" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="0" cy="50" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="6.6987298" cy="25" r="5" transform="translate(5 5)"/> ' +
        '<circle cx="25" cy="6.6987298" r="5" transform="translate(5 5)"/> ' +
        '<animateTransform attributeType="xml" attributeName="transform" type="rotate" from="0 55 55" to="360 55 55" dur="3s" repeatCount="indefinite" /> </g> ' +
        '</svg>'};
  }

    loadFromServer() {
        $.ajax({
            url: this.url,
            datatype: 'json',
            cache: false,
            success: function (data) {
                this.setState({data: data})
            }.bind(this)
        })
    }


    componentDidMount() {
        this.loadFromServer();
    }

    handleClick() {
        this.setState(prevState => ({
          isToggleOn: !prevState.isToggleOn
        }));
      }

    render() {
        if (this.state.data) {

            const svg_image = <SVGInline svg={this.state.data}/>;
            console.log(this.props.message)
            if (this.state.isToggleOn){
                return <div onClick={this.handleClick} style={{border: "1px solid black"}}>{svg_image}</div>
            }
            else{
                return <div onClick={this.handleClick}>{svg_image}</div>
            }

        }
        else {
            return (<FillMe />)
        }
    }

}


function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

// Specific Lists

class CompoundList extends GenericList {
        constructor(props) {
            super(props);
            this.url = COMPOUNDS_URL
            this.interval = 10000
            this.render_method = function (data, index) {
                return <CompoundView key={data.id} my_id={data.id} />
            }
        }
};

class TargetList extends GenericList {
        constructor(props) {
            super(props);
            this.url = TARGET_URL
            this.interval = 1000
            this.render_method = function (data, index) {
                return <ListGroupItem key={index}>{data.title}</ListGroupItem>
            }
        }
};

class MoleculeList extends GenericList {
        constructor(props) {
            super(props);
            this.url = MOLECULE_URL
            this.payload = {"target_id": null, "prot_id": null}
            this.interval = 1000
            this.render_method = function (data, index) {
                return <MoleculeView key={data.id} my_id={data.id} />
            }
        }
};

class ProteinList extends GenericList {

    constructor(props) {
            super(props);
            this.url = PROTEIN_URL
            this.interval = 1000
            this.render_method = function (data, index) {
                return <ListGroupItem key={data.id}>{data.code}</ListGroupItem>
            }
        }
}


// Specific Views

class CompoundView extends GenericView {

    constructor(props) {
        super(props);
        this.url = SVG_CMPD + props.my_id + '/'
    }
}

class MoleculeView extends GenericView {

    constructor(props) {
        super(props);
        this.url = SVG_MOL + props.my_id + '/'
    }
}


class NGLView extends React.Component {


    constructor(props) {
        super(props);
        // Create NGL Stage object
        this.stage = new NGL.Stage("viewport");
        // Handle window resizing
        window.addEventListener("resize", function (event) {
            this.stage.handleResize();
        }, false);
        this.focus_var = 95;
        // this.stage.mouseControls.add("clickPick-left",showPick);
    }



    show_mol() {
        NProgress.start();
        Promise.all([
            this.stage.loadFile(prot_url, {ext: "pdb"}),
            this.stage.loadFile(lig_data, {ext: "sdf"})]
        ).then(function (ol) {
            var cs = NGL.concatStructures(
                "concat",
                ol[0].structure.getView(new NGL.Selection("not ligand")),
                ol[1].structure.getView(new NGL.Selection(""))
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
}

class TotalView extends React.Component {


    constructor(props) {
        super(props);
        this.div_id = props.div_id;
    }

    render() {
        return <Row>
            <Col xs={2} md={2}>
                <TargetList />
            </Col>
            <Col xs={4} md={4}>
                <MoleculeList />
            </Col>
            <Col xs={4} md={4}>
                <NGLView />
            </Col>
        </Row>
    }
}


// The links between data and what is rendered
// ReactDOM.render(target_div, document.getElementById('targets'))
ReactDOM.render(<TotalView key="main_app"/>, document.getElementById('app'))

