import { h, Component } from 'preact';

class Checkbox extends Component {
    handleChange(e) {
        const { onChange, action } = this.props;
        // console.log(onChange);
        // console.log(action);
        // console.log(e.target.checked);
        onChange(action(e.target.checked));
    }
    render() {
        const { label, isChecked, id } = this.props;
        return (
            <div className="checkbox">
                <label>
                    <input
                        type="checkbox"
                        value={label}
                        id={id}
                        checked={isChecked}
                        onChange={this.handleChange.bind(this)}
                    />
                    {label}
                </label>
            </div>
        )
    }
}

export default Checkbox;