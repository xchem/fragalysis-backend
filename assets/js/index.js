import React from 'react'
import ReactDOM from 'react-dom'


function Welcome(props) {
  return <h1>Hello there, {props.name}</h1>;
}

const element = <Welcome name="anthony" />;
ReactDOM.render(
  element,
  document.getElementById('react')
);