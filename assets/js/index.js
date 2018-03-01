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
        this.mol_list = new Array();
        this.state = {mol_params: {"target_id": 1}
        }
    }
    
    onTargetChecked(target){
        // Now pass this to the molecule div
        this.setState(prevState => ({
          mol_params:  {"target_id": target}
        }));
    }

    onMolChecked(mol){
        // Now add or remove this to the list
    }


    render() {
        return <Row>
            <Col xs={2} md={2}>
                <TargetList communicateChecked={this.onTargetChecked}/>
            </Col>
            <Col xs={4} md={4}>
                <MoleculeList get_params={this.state.mol_params} communicateChecked={this.onMolChecked}/>
            </Col>
            <Col xs={6} md={6} >
                <NGLView mol_id={this.mol_list}/>
            </Col>
        </Row>
    }
}


// The links between data and what is rendered
ReactDOM.render(<TotalView key="main_app" />, document.getElementById('app'))