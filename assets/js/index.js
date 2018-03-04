import React from 'react';
import ReactDOM from 'react-dom';
import { connect, createStore } from 'redux'
import { Col, Row} from 'react-bootstrap'
import { NGLView } from './components/nglComponents'
import { TargetList, MoleculeList, MolGroupList } from './components/apiComponents'
import { createStore } from 'redux'
import { Provider } from 'react-redux'
import { nglReducers } from './reducers/nglReducers'
 
let store = createStore(nglReducers)

class TotalView extends React.Component {

    constructor(props) {
        super(props);
        this.div_id = props.div_id;
        this.onTargetChecked = this.onTargetChecked.bind(this)
        this.onMolChecked = this.onMolChecked.bind(this)
        this.onGroupChecked = this.onGroupChecked.bind(this)
        this.dummy_target = 1
        this.state = {mol_params: {"prot_id__target_id": this.dummy_target},
            group_params: {"target_id": this.dummy_target},
            mol_dict: false, clear_all: false
        }
    }
    
    onTargetChecked(target){
        // Now pass this to the molecule div
        this.setState(prevState => ({
          mol_params:  {"prot_id__target_id": target},
            group_params: {"target_id": target,"group_type": "MC"}
        }));
    }

    onMolChecked(mol,prot_id,isToggleOn){
        // Now add this to NGL
        this.setState(prevState => ({
            mol_dict: {"mol_id": mol, "prot_id": prot_id, "toggle": isToggleOn}
        }));
    }

    onGroupChecked(mol_list){
        this.setState(prevState => ({
          mol_params:  {"mol_groups": mol_list}
        }));
    }

    onMolCleared(){
        this.setState(prevState => ({
            clear_all: false
        }));
    }


    render() {
        return <Provider store={store}>
            <Row>
                <Col xs={1} >
                    <TargetList communicateChecked={this.onTargetChecked}/>
                </Col>
                <Col xs={1}>
                    <MolGroupList get_params={this.state.group_params} communicateChecked={this.onGroupChecked}/>
                </Col>
                <Col xs={4}>
                    <MoleculeList get_params={this.state.mol_params} communicateChecked={this.onMolChecked}/>
                </Col>
                <Col xs={6} md={6} >
                    <NGLView mol_dict={this.state.mol_dict} communicateCleared={this.onMolCleared}/>
                </Col>
            </Row>
        </Provider>
    }
}


// The links between data and what is rendered
ReactDOM.render(<TotalView key="main_app" />, document.getElementById('app'))