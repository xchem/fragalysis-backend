import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import $ from 'jquery';

function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

class GenericList extends React.Component {

    constructor(props) {
    super(props);
    this.loadFromServer = this.loadFromServer.bind(this);
    this.state = { data: [] };
    this.render_method = <h3 key={props.data.id}>{props.data.title}</h3>
  }

  loadFromServer() {
        $.ajax({
            url: this.props.url,
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
            this.pollInterval)
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
            render_method = render_title
        }
}
class TargetListTwo extends GenericList {
        constructor(props) {
            super(props);
            render_method = <RenderTitleTwo data={data}/>
        }
}

function render_title(props){
    return <h3 key={props.data.id}>{props.data.title}</h3>
}


function Welcome(props) {
  return <h1>Hello there, {props.name}</h1>;
}

const element = <Welcome name="anthony" />;
const target_div = <TargetList url='/v0.1/targets/' pollInterval={5000} />;

// The links between data and what is rendered
ReactDOM.render(element, document.getElementById('react'));
ReactDOM.render(target_div, document.getElementById('targets'))