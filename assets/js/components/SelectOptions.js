import { h, Component } from 'preact';

// TODO: Programatically subscribe to updates to NGL
class SelectGroup extends Component {
    // Action Creator
    handleChange(e) {
        const { onChange, action } = this.props;
        onChange(action(e.target.value));
    }
    render() {
        // Arguments supplied by SelectGroup component
        const { label, name, options, selected } = this.props;
        // Push options into an array
        const optionsList = [];

        console.log('***** ' + label + ' selected: ' + selected + ' ******');
        console.log(options);

        // Loop thru each Options
        options.forEach((option) => {
            optionsList.push(
                <Option
                    option={option}
                    selected={selected}
                />
            );
        });

        return (
            <div className='form-group'>
                <label className='col-sm-5 control-label' for={name + 'Select'}>
                    {label}
                </label>
                <div class='col-sm-7'>
                    <select onChange={this.handleChange.bind(this)} className='form-control input-sm'>
                        {optionsList}
                    </select>
                </div>
            </div>
        );
    }
}

// List of Options from Select Group
class Option extends Component {
    render() {
        const { option, selected } = this.props;
        console.log('++++++ Selected: ' + selected + ' ++++++');
        // Ternary operator
        const isSelected = (option.value === selected) ? 'selected' : null;
        return (
            <option selected={isSelected} value={option.value}>
                {option.label}
            </option>
        );
    }
}

export default SelectGroup;