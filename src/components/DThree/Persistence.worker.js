import numeric from 'numeric';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';

let data = [];

function binZ(pz, bin) {
  if (pz < -2.5 || pz > 2.5) console.log(pz);
  pz += 5.0;
  pz *= 10;
  bin[parseInt(pz)]++
}

function createBinZ(bin) {
  let result = [];
  for (let i = 0; i < 100; i++) {
    result.push([(i-50)/10.0 + 0.05, bin[i]]);
  }
  return result;
}

function binX(px, bin, crick = false) {
  if (px > -1.5) console.log("significant watson px (Å) = " + px);
  if (!crick) {
    px += 6.0;
    px *= 10;
    bin[parseInt(px)]++;
  } else {
    px *= 10;
    bin[parseInt(px)]++;
  }
}

function createBinX(bin, crick = false) {
  let result = [];
  for (let i = 0; i < 60; i++) {
    if (!crick) result.push([(i-60)/10.0 + 0.05, bin[i]]);
    else result.push([i/10.0 + 0.05, bin[i]]);
  }
  return result;
}

function binY(py, bin, crick = false) {
  if (!crick) {
    py -= 6;
    py *= 10;
    bin[parseInt(py)]++;
  } else {
    py += 12;
    py *= 10;
    bin[parseInt(py)]++;
  }
}

function createBinY(bin, crick = false) {
  let result = [];
  for (let i = 0; i < 60; i++) {
    if (!crick) result.push([i/10.0 + 6 + 0.05, bin[i]]);
    else result.push([(i/10) - 12 + 0.05, bin[i]]);
  }
  return result;
}

function binTwist(twist, bin) {
    twist -= 15;
    twist *= 2.5;
    bin[parseInt(twist)]++;
}

function createBinTwist(bin) {
  let result = [];
  for (let i = 0; i < 100; i++) {
    result.push([i/2.5 + 15 + 0.200, bin[i]]);
  }
  return result;
}

function binRoll(roll, bin) {
    roll += 25;
    roll *= 2.0;
    bin[parseInt(roll)]++;
}

function createBinRoll(bin) {
  let result = [];
  for (let i = 0; i < 100; i++) {
    result.push([(i-50)/2.0+0.250, bin[i]])
  }
  return result;
}

function binSlide(slide, bin) {
    if (slide < -2.0) console.log("SLIDE: " + slide);
    slide += 4;
    slide *= 10;
    bin[parseInt(slide)]++;
}

function createBinSlide(bin) {
  let result = [];
  for (let i = 0; i < 80; i++) {
    result.push([(i-40)/10.0+0.05, bin[i]]);
  }
  return result;
}

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
	//***start new code (histograms)
    //let watson = [];
    //let crick = [];
    let binTilt1 = [];
    let binRoll1 = [];
    let binTwist1 = [];
    let binSlide1 = [];
    let binWZ = [];
    let binCZ = [];
    let binWX = [];
    let binWY = [];
    binTilt1.length = 100;
    binRoll1.length = 100;
    binTwist1.length = 100;
    binWZ.length = 100;
    binCZ.length = 100;
    binSlide1.length = 80;
    binWX.length = 60;
    binWY.length = 60;
    for (let i = 0; i < 100; i++) {
      binCZ[i] = 0;
      binWZ[i] = 0;
      binTilt1[i] = 0;
      binRoll1[i] = 0;
      binTwist1[i] = 0;
      if (i < 60) binWY[i] = 0;
      if (i < 60) binWX[i] = 0;
      if (i < 80) binSlide1[i] = 0;
    }
	//***end new code
	for (let i = 0; i < parseInt(data.numSamples); i++) {
	  for (let i = 0; i < 30; i++) v[0][i] = factor[0][i]*ref.randg();
	  let step = numeric.add(numeric.transpose(numeric.dot(eigen.eigenvectors, numeric.transpose(v))), step0);
	  //console.log(step);
	  for (let i = 0; i < 30; i++) coord[i] = step[0][i];
	  ref.get30Coordinates(coord, data.tetramer, true);
	  //***start new code
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
	  binZ(pw[2], binWZ);
	  binX(pw[0], binWX);
	  binY(pw[1], binWY);
	  binZ(pc[2], binCZ);
	  binRoll(stepParameters[1][0], binTilt1);
	  binRoll(stepParameters[1][1], binRoll1);
	  binTwist(stepParameters[1][2], binTwist1);
	  binSlide(stepParameters[1][4], binSlide1);
      //watson.push(pw);
	  //crick.push(pc);
	  //***end new code
	  let bpstep = ref.getStep();
	  A = numeric.add(A, bpstep);
	}
	A = numeric.mul(1/parseFloat(data.numSamples), A);
	for (let i = 0; i < 2000; i++) {
	  A = numeric.dot(A, A);
	}
	let watsonBinZ = createBinZ(binWZ);
	let crickBinZ = createBinZ(binCZ);
	let watsonBinX = createBinX(binWX);
	let watsonBinY = createBinY(binWY);
	let BinTilt = createBinRoll(binTilt1);
	let BinRoll = createBinRoll(binRoll1);
	let BinTwist = createBinTwist(binTwist1);
	let BinSlide = createBinSlide(binSlide1);
	console.log(JSON.stringify(createBinZ(binWZ)));
	console.log(JSON.stringify(createBinZ(binCZ)));
	console.log("remainder: " + A[0][3] + " " + A[1][3]);
	return [A[2][3], watsonBinZ, crickBinZ, BinRoll, BinSlide, watsonBinX, watsonBinY, BinTwist, BinTilt];
}

onmessage = (e) => {
	data = e.data;
	let P = persistenceLength();
	postMessage(P);
}

