import React from 'react';
import ReactDOM from 'react-dom';
import { createStore } from 'redux'
import { Col, Row} from 'react-bootstrap'
import { NGLView } from './components/ngl_components'
import { TargetList, MoleculeList } from './components/api_components'
import { app } from './reducers/reducers'

let store = createStore(app);

class TotalView extends React.Component {

    constructor(props) {
        super(props);
        this.div_id = props.div_id;
    }

    render() {
        return <Row>
            <Col xs={2} md={2}>
                <TargetList />
            </Col>
            <Col xs={4} md={4}>
                <MoleculeList />
            </Col>
            <Col xs={6} md={6}>
                <NGLView />
            </Col>
        </Row>
    }
}


// The links between data and what is rendered
ReactDOM.render(<TotalView key="main_app"/>, document.getElementById('app'))