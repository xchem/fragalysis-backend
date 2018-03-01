import React from 'react';
import ReactDOM from 'react-dom';
import { connect, createStore } from 'redux'
import { Col, Row} from 'react-bootstrap'
import { NGLView } from './components/ngl_components'
import { TargetList, MoleculeList } from './components/api_components'
import { app } from './reducers/reducers'

let store = createStore(app);

class TotalView extends React.Component {

    constructor(props) {
        super(props);
        this.div_id = props.div_id;
        this.onTargetChecked = this.onTargetChecked.bind(this)
        this.onMolChecked = this.onMolChecked.bind(this)
    }
    
    onTargetChecked(target){
        alert(target);
    }

    onMolChecked(mol){
        alert(mol);
    }


    render() {
        return <Row>
            <Col xs={2} md={2}>
                <TargetList communicateChecked={this.onTargetChecked}/>
            </Col>
            <Col xs={4} md={4}>
                <MoleculeList communicateChecked={this.onMolChecked}/>
            </Col>
            <Col xs={6} md={6} >
                <NGLView />
            </Col>
        </Row>
    }
}


// The links between data and what is rendered
ReactDOM.render(<TotalView key="main_app" />, document.getElementById('app'))