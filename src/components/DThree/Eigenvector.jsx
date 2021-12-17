import Table from 'react-bootstrap/Table';
import {Component} from 'react';

class Eigenvector extends Component {

	render() {
		//console.log(this.props.vector);
		let values = {};
		if (this.props.vector) { 
		  values = {
            pair1: this.props.vector.slice(0, 6),
            pho1:  this.props.vector.slice(6, 12),
            step:  this.props.vector.slice(12, 18),
            pho2:  this.props.vector.slice(18, 24),
            pair2: this.props.vector.slice(24, 30),
          };
        }
		return (
		<><b>Eigenvector (linear combination, norm = 1)</b>
		<Table striped bordered><tbody>
		<tr>
		<td>Pair 1</td>
		{values.pair1.map((element) => (<td key={element}>{element.toFixed(7)}</td>))}
		</tr><tr>
		<td>Pho 1</td>
		{values.pho1.map((element) => (<td key={element}>{element.toFixed(7)}</td>))}
		</tr><tr>
        	<td>PairStep</td>
        	{values.step.map((element) => (<td key={element}>{element.toFixed(7)}</td>))}
		</tr><tr>
        	<td>Pho 2</td>
        	{values.pho2.map((element) => (<td key={element}>{element.toFixed(7)}</td>))}
		</tr><tr>
        	<td>Pair 2</td>
        	{values.pair2.map((element) => (<td key={element}>{element.toFixed(7)}</td>))}
		</tr>
		</tbody></Table></>
		);

	}

}

//{"Norm = " + Math.sqrt(this.props.vector.reduce((acc, value) => acc+value**2,0))}
export default Eigenvector;

