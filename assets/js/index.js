import React from 'react';
import ReactDOM from 'react-dom';
import { connect, createStore } from 'redux'
import { Col, Row} from 'react-bootstrap'
import { NGLView } from './components/ngl_components'
import { TargetList, MoleculeList } from './components/api_components'

class TotalView extends React.Component {

    constructor(props) {
        super(props);
        this.div_id = props.div_id;
        this.onTargetChecked = this.onTargetChecked.bind(this)
        this.onMolChecked = this.onMolChecked.bind(this)
        this.state = {mol_params: {"prot_id__target_id": 1},
            mol_dict: false
        }
    }
    
    onTargetChecked(target){
        // Now pass this to the molecule div
        this.setState(prevState => ({
          mol_params:  {"prot_id__target_id": target}
        }));
    }

    onMolChecked(mol,prot_id,isToggleOn){
        // Now add this to NGL
        this.setState(prevState => ({
            mol_dict: {"mol_id": mol, "prot_id": prot_id, "toggle": isToggleOn}
        }));
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
                    <NGLView mol_dict={this.state.mol_dict}/>
                </Col>
        </Row>
    }
}


// The links between data and what is rendered
ReactDOM.render(<TotalView key="main_app" />, document.getElementById('app'))