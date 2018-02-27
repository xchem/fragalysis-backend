import React from 'react';
import ReactDOM from 'react-dom';
import '../css/index.css';
import 'generic_components';

class CompoundList extends React.Component {
    loadCompoundsFromServer() {
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
        this.loadCompoundsFromServer();
        setInterval(this.loadCompoundsFromServer,
            this.props.pollInterval)
    }

    render() {
        if (this.state.data) {
            console.log("Refreshing compound load")
            //
            return this.props.data.map(data => (
                <CompoundView data={data}/>
            ));
        }
        else {
            return (<FillMe />)
        }
    }
}

class CompoundView extends React.Component {

    loadImageFromServer() {
        $.ajax({
            url: "/viewer/img_from_smiles",
            // POST THE SMILES
            datatype: 'json',
            cache: false,
            success: function (data) {
                this.setState({data: data})
            }.bind(this)
        })
    }

  render() {
    return (this.data);
  }
}


// The links between data and what is rendered
ReactDOM.render(<CompoundList url='/insert/path/here/' pollInterval={1000} />, document.getElementById('container'))