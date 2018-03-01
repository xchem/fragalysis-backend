/**
 * Created by abradley on 01/03/2018.
 */
/**
 * Created by abradley on 01/03/2018.
 */
import $ from 'jquery';
import SVGInline from "react-svg-inline"
import React from 'react';
import { ListGroup, Col } from 'react-bootstrap';

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

function FillMe(props) {
    return <h1>FILL ME UP PLEASE</h1>;
}

// Generic Classes
export class GenericList extends React.Component {

    constructor(props) {
    super(props);
        this.url = ''
        this.interval = 10000
        this.state = {data:[]}

        this.loadFromServer = this.loadFromServer.bind(this);
        this.handleOptionChange = this.handleOptionChange.bind(this);
  }

    handleOptionChange(changeEvent) {
        const new_value = changeEvent.target.value;
        this.setState(prevState => ({
          targetOn: new_value
        }));
        this.props.communicateChecked(new_value);
      }

  loadFromServer() {
        $.ajax({
            url: this.url,
            datatype: 'json',
            data: this.props.get_params,
            cache: false,
            success: function (data) {
                this.setState({data: data["results"]})
            }.bind(this)
        })
    }

    componentDidMount() {
        this.setState(prevState => (
        {data: []}
        ));
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
        this.selected_style = {border: "1px solid black"}
        this.not_selected_style = {}
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
        this.props.communicateChecked(this.props.my_id, this.props.prot_id, isToggleOn);
    }

    render() {
        if (this.state.data) {

            const svg_image = <SVGInline svg={this.state.data}/>;
            console.log(this.props.message)
            if (this.state.isToggleOn){
                return <Col xs={6} onClick={this.handleClick} style={this.selected_style}>{svg_image}</Col>
            }
            else{
                return <Col xs={6} onClick={this.handleClick} style={this.not_selected_style}>{svg_image}</Col>
            }
        }
        else {
            return (<FillMe />)
        }
    }
}