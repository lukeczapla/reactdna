import numeric from 'numeric';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';

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
	  ref.get30Coordinates(coord, data.tetramer, true);
	  // start new code
	  let PhoW = ref.getAtomSets().atoms[1][1];
	  let PhoC = ref.getAtomSets().atoms[4][1];
      ref3.getBasePlanes(ref.getAtomSets());
      let stepParameters = ref3.getParameters();
	  let midframe = ref.getMidFrame();
	  let px = PhoW.x - midframe[0][3];
	  let py = PhoW.y - midframe[1][3];
	  let pz = PhoW.z - midframe[2][3];
	  let pw = [
		px*midframe[0][0]+py*midframe[1][0]+pz*midframe[2][0],
		px*midframe[0][1]+py*midframe[1][1]+pz*midframe[2][1],
		px*midframe[0][2]+py*midframe[1][2]+pz*midframe[2][2]
	  ];
	  px = PhoC.x - midframe[0][3];
	  py = PhoC.y - midframe[1][3];
	  pz = PhoC.z - midframe[2][3];	
	  let pc = [
		px*midframe[0][0]+py*midframe[1][0]+pz*midframe[2][0],
		px*midframe[0][1]+py*midframe[1][1]+pz*midframe[2][1],
		px*midframe[0][2]+py*midframe[1][2]+pz*midframe[2][2]
	  ];
	  // end new code  
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

