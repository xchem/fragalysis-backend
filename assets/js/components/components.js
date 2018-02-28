import { h, Component } from 'preact';
import * as actions from '../actions/actions.js';
import SelectGroup from './SelectOptions';
import Checkbox from './Checkbox';

// TODO: Break up into different files
class Header extends Component {
    render() {
        const store = this.props.store;
        const storeState = store.getState();
        console.log("------------ HEADER SET -------------");
        console.log(storeState);
        return (
            <div>
                <h1>{storeState.pdbId}</h1>
                <h4>{storeState.structureTitle}</h4>
            </div>
        );
    }
}

// Display Options Panel
class DisplayOptions extends Component {
    render() {
        const store = this.props.store;
        const storeState = store.getState();
        console.log("------------ DISPLAY OPTIONS -------------");
        console.log(storeState);

        // Uses the SelectGroup class to generate dropdown selection
        return (
            <div className='form-horizontal'>
                <SelectGroup
                    label='Assembly'
                    name='assembly'
                    options={storeState.assemblyOptions}
                    selected={storeState.assembly}
                    onChange={store.dispatch}
                    action={actions.setAssembly}/>
                <SelectGroup
                    label='Color'
                    name='color'
                    options={storeState.colorOptions}
                    selected={storeState.color}
                    onChange={store.dispatch}
                    action={actions.setColor}/>
                <SelectGroup
                    label='Style'
                    name='style'
                    options={storeState.styleOptions}
                    selected={storeState.style}
                    onChange={store.dispatch}
                    action={actions.setStyle}/>
                <Checkbox
                    label='Spin'
                    isChecked={storeState.spin}
                    id="spinCheckbox"
                    onChange={store.dispatch}
                    action={actions.setSpin}/>
                <Checkbox
                    label='Water'
                    isChecked={storeState.water}
                    id="waterVisibilityCheckbox"
                    onChange={store.dispatch}
                    action={actions.setWater}/>
                <Checkbox
                    label='Hydrogen'
                    isChecked={storeState.hydrogen}
                    id="hydrogenVisibilityCheckbox"
                    onChange={store.dispatch}
                    action={actions.setHydrogen}/>
            </div>
        );
    }
}

export { Header, DisplayOptions };

