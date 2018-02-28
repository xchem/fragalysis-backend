import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import $ from 'jquery';
import SVGInline from "react-svg-inline"


function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

class GenericList extends React.Component {

    constructor(props) {
    super(props);
        this.url = '/v0.1/'
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

class GenericView extends React.Component{


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
            return <div><SVGInline svg={this.state.data}/></div>
        }
        else {
            return (<FillMe />)
        }
    }

}

class CompoundList extends GenericList {
        constructor(props) {
            super(props);
            this.url = '/v0.1/compounds/'
            this.interval = 1000
            this.render_method = function (data, index) {
                return <CompoundView key={data.id} my_id={data.id} />
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