/**
 * Created by abradley on 01/03/2018.
 */
import { ListGroupItem } from 'react-bootstrap';
import { GenericView, GenericList } from './general_components';
import React from 'react';

import {toggleComplex} from '../actions/actions'

// Basic config of the API
const BASE_API = "/v0.1/"
const TARGET_URL = BASE_API+"targets/"
const COMPOUNDS_URL = BASE_API+"compounds/"
const MOLECULE_URL = BASE_API+"molecules/"
const PROTEIN_URL = BASE_API+"proteins/"
const SVG_CMPD = '/viewer/img_from_cmpd_pk/'
const SVG_MOL = '/viewer/img_from_mol_pk/'
const GENERIC_INTERVAL = 100000
    
export class CompoundList extends GenericList {
        constructor(props) {
            super(props);
            this.url = COMPOUNDS_URL
            this.interval = GENERIC_INTERVAL
            this.render_method = function (data, index) {
                return <CompoundView key={data.id} my_id={data.id} />
            }
        }
};

export class TargetList extends GenericList {


    constructor(props) {
            super(props);
            this.url = TARGET_URL
            this.interval = GENERIC_INTERVAL
            this.render_method = function (data, index) {
                return <ListGroupItem key={index} >
                    <label>
                        <input type="radio" value={data.title} checked={this.state.target === data.title} onChange={this.handleOptionChange}/>
                        {data.title}
                    </label>
                </ListGroupItem>
            }
        }
    onComponentMount(){

        this.setState(prevState => ({
                targetOn: ""
            }))
    }

};

export class MoleculeList extends GenericList {
        constructor(props) {
            super(props);
            this.url = MOLECULE_URL
            this.interval = 1000
            this.render_method = function (data, index) {
                return <MoleculeView key={data.id} my_id={data.id} />
            }
        }
};

export class ProteinList extends GenericList {

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

