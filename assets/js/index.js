import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import GenericList from 'generic_components.js';
import GenericView from 'generic_components.js';


function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

class CompoundList extends GenericList {
        constructor(props) {
            super(props);
            this.url = '/v0.1/compounds/'
            this.interval = 10000
            this.render_method = function (data, index) {
                return <CompoundView key={data.id} my_id={data.id} />
            }
        }
};

class TargetList extends GenericList {
        constructor(props) {
            super(props);
            this.url = '/v0.1/targets/'
            this.interval = 1000
            this.render_method = function (data, index) {
                return <h3 key={index}>{data.title}</h3>
            }
        }
};


class CompoundView extends GenericView {

    constructor(props) {
        super(props);
        this.url = '/viewer/img_from_cmpd_pk/' + props.my_id + '/'
    }

}

class TargetListLink extends TargetList {
        constructor(props) {
            super(props);
            this.render_method = function (data, index) {
                return <li key={index}><a>{data.title}</a></li>
            }
        }
};


function Welcome(props) {
  return <h1>Hello there, {props.name}</h1>;
}

const element = <Welcome name="anthony" />;
const target_div = <CompoundList />;

// The links between data and what is rendered
ReactDOM.render(element, document.getElementById('react'));
ReactDOM.render(target_div, document.getElementById('targets'))