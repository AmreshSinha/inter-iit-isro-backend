import React from "react";
import "./Switch.css";

const ToggleSwitch = () => {
return (
	<div className="container">
	<div className="toggle-switch">
		<input type="checkbox" className="checkbox"
			name='GraphChange' id='GraphChange' />
		<label className="label" htmlFor='GraphChange'>
		<span className="inner" />
		<span className="switch" />
		</label>
	</div>
	</div>
);
};

export default ToggleSwitch;
