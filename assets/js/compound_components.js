import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import 'generic_components';

class GenericList extends React.Component {

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

    getInitialState() {
        return {data: []}
    }

    componentDidMount() {
        this.loadFromServer();
        setInterval(this.loadFromServer,
            this.props.pollInterval)
    }
}


class ProjectList extends GenericList {

    render() {
        if (this.state.data) {
            console.log("Refreshing compound load")
            //
            return this.props.data.map(data => (
                <ProjectView data={data}/>
            ));
        }
        else {
            return (<FillMe />)
        }
    }
}

class ProjectView extends React.Component {

  render() {
    return (this.data);
  }
}

function App(props) {
    return <h1>Hello Anthony!!!</h1>;
}


// The links between data and what is rendered
ReactDOM.render(<ProjectList url='/projects/' pollInterval={1000} />, document.getElementById('projects'))