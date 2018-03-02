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
const GROUP_URL = BASE_API+"molgroup/"
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
        this.handleOptionChange = this.handleOptionChange.bind(this);
            this.render_method = function (data, index) {
                return <ListGroupItem key={index} >
                    <label>
                        <input type="radio" value={data.id} checked={this.state.targetOn == data.id} onChange={this.handleOptionChange}/>
                        {data.title}
                    </label>
                </ListGroupItem>
            }
        }

    handleOptionChange(changeEvent) {
        const new_value = changeEvent.target.value;
        this.setState(prevState => ({
          targetOn: new_value
        }));
        this.props.communicateChecked(new_value);
      }

};


export class MolGroupList extends GenericList {

    constructor(props) {
        super(props);
        this.url = GROUP_URL
        this.interval = 1000
        this.handleOptionChange = this.handleOptionChange.bind(this);
        this.render_method = function (data, index) {
            return <ListGroupItem key={index} >
                <label>
                    <input type="radio" value={data.id} checked={this.state.groupOn == data.id} onChange={this.handleOptionChange}/>
                    SITE {data.index}
                </label>
            </ListGroupItem>
        }
    }

    handleOptionChange(changeEvent) {
        const new_value = changeEvent.target.value;
        this.setState(prevState => ({
          groupOn: new_value
        }));
        this.props.communicateChecked(new_value);
      }

}

export class MoleculeList extends GenericList {
        constructor(props) {
            super(props);
            this.url = MOLECULE_URL
            this.get_params={"width": 100, "height": 100}
            this.communicateChecked = props.communicateChecked;
            this.interval = 1000
            this.render_method = function (data, index) {
                return <MoleculeView get_params={this.get_params} communicateChecked={this.communicateChecked} key={data.id} prot_id={data.prot_id} my_id={data.id} />
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

