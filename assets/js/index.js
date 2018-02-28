import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import $ from 'jquery';
import SVGInline from "react-svg-inline"

const BASE_API = "/v0.1/"
const TARGET_URL = BASE_API+"targets/"
const COMPOUNDS_URL = BASE_API+"compounds/"
const MOLECULE_URL = BASE_API+"molecules/"
const PROTEIN_URL = BASE_API+"proteins/"
const SVG_CMPD = '/viewer/img_from_cmpd_pk/'
const SVG_MOL = '/viewer/img_from_mol_pk/'


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
            return this.state.data.map((data, index) => (
                this.render_method(data,index)
            ));
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
        this.state = {data: '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="110px" height="110px"><g>' +
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

    render() {
        if (this.state.data) {
            console.log(this.props.message)
            return <SVGInline svg={this.state.data}/>
        }
        else {
            return (<FillMe />)
        }
    }

}


function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

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
                return <h3 key={index}>{data.title}</h3>
            }
        }
};

class MoleculeList extends GenericList {
        constructor(props) {
            super(props);
            this.url = MOLECULE_URL
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
                return <h3 key={data.id} >{data.code}</h3>
            }
        }
}

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


// The different DIV elements
const compound_div = <CompoundList />;
const target_div = <TargetList />;
const protein_div = <ProteinList />;
const molecule_div = <MoleculeList />;
// The links between data and what is rendered
ReactDOM.render(element, document.getElementById('react'));
ReactDOM.render(target_div, document.getElementById('targets'))
ReactDOM.render(molecule_div, document.getElementById('molecules'))
ReactDOM.render(protein_div, document.getElementById('proteins'))
ReactDOM.render(compound_div, document.getElementById('compounds'))