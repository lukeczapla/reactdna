import numeric from 'numeric';
import * as ref from '../References/References.jsx';

let data = [];

function persistenceLength() {
	let cov = data.cov;
	let F = numeric.inv(cov);
	let eigen = ref.jeigen(F);
	let factor = [[]];
	factor[0].length = 30;
	for (let i = 0; i < 30; i++) {
	  factor[0][i] = 1/Math.sqrt(eigen.eigenvalues[i]);
	}
	let step0 = [data.mean];
	let A = numeric.identity(4);
	for (let i = 0; i < 4; i++) A[i][i] = 0;
	let v = [[]];
	v[0].length = 30;
	let coord = [];
	coord.length = 30;

	for (let i = 0; i < parseInt(data.numSamples); i++) {
	  for (let i = 0; i < 30; i++) v[0][i] = factor[0][i]*ref.randg();
	  let step = numeric.add(numeric.transpose(numeric.dot(eigen.eigenvectors, numeric.transpose(v))), step0);
	  //console.log(step);
	  for (let i = 0; i < 30; i++) coord[i] = step[0][i];
	  ref.get30Coordinates(coord, data.tetramer);
	  let bpstep = ref.getStep();
	  A = numeric.add(A, bpstep);
	}
	A = numeric.mul(1/parseFloat(data.numSamples), A);
	for (let i = 0; i < 2000; i++) {
	  A = numeric.dot(A, A);
	}
	console.log(A[0][3] + " " + A[1][3]);
	return A[2][3];

}

onmessage = (e) => {
	data = e.data;
	//console.log(data);
	//console.log("calculating");
	let P = persistenceLength();
	postMessage(P);
}

