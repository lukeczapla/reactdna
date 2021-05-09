import './DThree.css';
import React, {Component} from 'react';
import PersistenceWorker from './Persistence.worker';
import * as trois from 'three';
import numeric from 'numeric';
import Plot from 'react-plotly.js';
import Eigenvector from './Eigenvector.jsx';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';

class DThree extends Component {

  constructor(props) {
    super(props);
    this.state = {
      modeNum : '0',
      stepBy: 0,
      showIC: false,
      isDragging: false,
      useScale: true,
      numSamples: 50000,
      parameters3: [],
      phoW: [0,0,0],
      phoC: [0,0,0],
      startX: 0,
      startY: 0,
      translateX: 0,
      translateY: 0,
      lastTranslateX: 0,
      lastTranslateY: 0,
      evalue: '',
      midVal: '',
      mean: props.mean,
      cov: props.cov,
      tetramer: props.tetramer,
      eigenvalues: [],
      eigenvectors: [],
      eigenlist: [],
      processing: false,
      showPlot: false,
      pdbText: ''
    };
    this.scene = null;
    this.camera = null;
    this.renderer = null;
  }

  onMouseMove = ({clientX, clientY}) => {
    if (!this.state.isDragging) return;
    this.setState({
      translateX: clientX - this.state.startX + this.state.lastTranslateX,
      translateY: clientY - this.state.startY + this.state.lastTranslateY
    });
  };

  makePlots = (P) => {
        let watson = P[1];
        let crick = P[2];
        let roll = P[3];
        let slide = P[4];
        let watsonX = P[5];
        let watsonY = P[6];
        let twist = P[7];
        let tilt = P[8];
        let data1 = []; let data2 = []; let data3 = []; let data4 = []; let data5 = []; let data6 = []; let data7 = []; let data8 = [];
        let xvalues = []; let yvalues = [];
        let xvalues2 = []; let yvalues2 = [];
        let xvalues3 = []; let yvalues3 = [];
        let xvalues4 = []; let yvalues4 = [];
        let xvalues5 = []; let yvalues5 = [];
        let xvalues6 = []; let yvalues6 = [];
        let xvalues7 = []; let yvalues7 = [];
        let xvalues8 = []; let yvalues8 = [];        
        let max = 0; let max2 = 0; let max3 = 0; let max4 = 0; let max5 = 0; let max6 = 0; let max7 = 0; let max8 = 0;
        for (let i = 0; i < P[1].length; i++) {
            xvalues.push(watson[i][0]);
            yvalues.push(watson[i][1]);
            if (watson[i][1] > max) max = watson[i][1];
            xvalues2.push(crick[i][0]);
            yvalues2.push(crick[i][1]);
            if (crick[i][1] > max2) max2 = crick[i][1];
        }
        for (let i = 0; i < P[3].length; i++) {
            xvalues3.push(roll[i][0]);
            yvalues3.push(roll[i][1]);
            if (roll[i][1] > max3) max3 = roll[i][1];
        }
        for (let i = 0; i < P[4].length; i++) {
            xvalues4.push(slide[i][0]);
            yvalues4.push(slide[i][1]);
            if (slide[i][1] > max4) max4 = slide[i][1];
        }
        for (let i = 0; i < P[5].length; i++) {
            xvalues5.push(watsonX[i][0]);
            yvalues5.push(watsonX[i][1]);
            if (watsonX[i][1] > max5) max5 = watsonX[i][1];
        }
        for (let i = 0; i < P[6].length; i++) {
            xvalues6.push(watsonY[i][0]);
            yvalues6.push(watsonY[i][1]);
            if (watsonY[i][1] > max6) max6 = watsonY[i][1];
        }
        for (let i = 0; i < P[7].length; i++) {
            xvalues7.push(twist[i][0]);
            yvalues7.push(twist[i][1]);
            if (twist[i][1] > max7) max7 = twist[i][1];
        }
        for (let i = 0; i < P[8].length; i++) {
        	xvalues8.push(tilt[i][0]);
        	yvalues8.push(tilt[i][1]);
			if (tilt[i][1] > max8) max8 = tilt[i][1];
        }
        data1.push({
            x: xvalues,
            y: yvalues,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data2.push({
            x: xvalues2,
            y: yvalues2,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data3.push({
            x: xvalues3,
            y: yvalues3,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data4.push({
            x: xvalues4,
            y: yvalues4,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data5.push({
            x: xvalues5,
            y: yvalues5,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data6.push({
            x: xvalues6,
            y: yvalues6,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data7.push({
            x: xvalues7,
            y: yvalues7,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        data8.push({
            x: xvalues8,
            y: yvalues8,
            name: this.state.tetramer,
            type: 'bar',
            mode: 'markers',
            marker: {
                symbol: 'square',
                size: 5,
                line: {
                    width: 2
                }
            }
        });
        let dtick1 = 10; let dtick2 = 10; let dtick3 = 10; let dtick4 = 10; let dtick5 = 10; let dtick6 = 10; let dtick7 = 10; let dtick8 = 10;
        if (max > 200) dtick1 = 50; if (max > 500) dtick1 = 100; if (max > 2000) dtick1 = 500;
        if (max > 5000) dtick1 = 1000; if (max > 20000) dtick1 = 5000; if (max > 50000) dtick1 = 10000;
        if (max2 > 200) dtick2 = 50; if (max2 > 500) dtick2 = 100; if (max2 > 2000) dtick2 = 500;
        if (max2 > 5000) dtick2 = 1000; if (max2 > 20000) dtick2 = 5000; if (max2 > 50000) dtick2 = 10000;
        if (max3 > 200) dtick3 = 50; if (max3 > 500) dtick3 = 100; if (max3 > 2000) dtick3 = 500;
        if (max3 > 5000) dtick3 = 1000; if (max3 > 20000) dtick3 = 5000; if (max3 > 50000) dtick3 = 10000;
        if (max4 > 200) dtick4 = 50; if (max4 > 500) dtick4 = 100; if (max4 > 2000) dtick4 = 500;
        if (max4 > 5000) dtick4 = 1000; if (max4 > 20000) dtick4 = 5000; if (max4 > 50000) dtick4 = 10000;
        if (max5 > 200) dtick5 = 50; if (max5 > 500) dtick5 = 100; if (max5 > 2000) dtick5 = 500;
        if (max5 > 5000) dtick5 = 1000; if (max5 > 20000) dtick5 = 5000; if (max5 > 50000) dtick5 = 10000;
        if (max6 > 200) dtick6 = 50; if (max6 > 500) dtick6 = 100; if (max6 > 2000) dtick6 = 500;
        if (max6 > 5000) dtick6 = 1000; if (max6 > 20000) dtick6 = 5000; if (max6 > 50000) dtick6 = 10000;
        if (max7 > 200) dtick7 = 50; if (max7 > 500) dtick7 = 100; if (max7 > 2000) dtick7 = 500;
        if (max7 > 5000) dtick7 = 1000; if (max7 > 20000) dtick7 = 5000; if (max7 > 50000) dtick7 = 10000;
        if (max8 > 200) dtick8 = 50; if (max8 > 500) dtick8 = 100; if (max8 > 2000) dtick8 = 500;
        if (max8 > 5000) dtick8 = 1000; if (max8 > 20000) dtick8 = 5000; if (max8 > 50000) dtick8 = 10000;
      	this.setState({
      	        showPlot: true,
      	        data1: data1,
                layout1: {
                  title: "zp values of Watson Phosphate",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
    	          xaxis: {
                    range: [-5, 5],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: 'zp value (Å)',
                    automargin: true,
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max+10],
                    showgrid: true,
                    showline: true,
                    dtick: dtick1,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
    			},
    			data2: data2,
                layout2: {
                  title: "zp values of Crick Phosphate",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [-5, 5],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    //title: 'number of days before today',
                    automargin: true,
                    title: "zp value (Å)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max2+10],
                    showgrid: true,
                    showline: true,
                    dtick: dtick2,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data3: data3,
                layout3: {
                  title: "Base-pair Step Roll",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [-25, 25],
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "roll value (degrees)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max3+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick3,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data4: data4,
                layout4: {
                  title: "Base-pair Step Slide",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [-4, 4],
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "slide value (Å)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max4+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick4,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data5: data5,
                layout5: {
                  title: "xp values of Watson phosphate",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [-6, 0],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "xp value (Å)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max5+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick5,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data6: data6,
                layout6: {
                  title: "yp values of Watson phosphate",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [6, 12],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "yp value (Å)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max6+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick6,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data7: data7,
                layout7: {
                  title: "Base-pair Step Twist",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [10, 60],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "twist value (degrees)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max7+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick7,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                },
                data8: data8,
                layout8: {
                  title: "Base-pair Step Tilt",
                  width: 1080,
                  height: 800,
                  margin: {
                    b: 150
                  },
                  font: {
                    size: 18
                  },
                  xaxis: {
                    range: [-25, 25],
                    //autorange: 'reversed',
                    showgrid: true,
                    showline: true,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    automargin: true,
                    title: "tilt value (degrees)",
                    linewidth: 1
                  },
                  yaxis: {
                    range: [0, max8+5],
                    showgrid: true,
                    showline: true,
                    dtick: dtick8,
                    gridwidth: 2,
                    gridcolor: '#777777',
                    title: "number of samples",
                    linewidth: 1
                  }
                }

    	});

  }

  persistenceLength = () => {
    this.setState({processing: true});
    let worker = new PersistenceWorker();
    worker.postMessage(this.state);
    worker.addEventListener("message", (e) => {
    	let P = e.data;
    	console.log(P);
    	alert("Persistence length in Angstroms: " + P[0]);
    	this.makePlots(P);
    	this.setState({processing: false});
    });
/*  // blocking (non-webworker, non-threaded) version:
	let cov = this.state.cov;
	//console.log(cov);
	let F = numeric.inv(cov);
	let eigen = ref.jeigen(F);
	let factor = [[]];
	factor[0].length = 30;
	for (let i = 0; i < 30; i++) {
	  factor[0][i] = 1/Math.sqrt(eigen.eigenvalues[i]);
	}
	let step0 = [this.state.mean];
	let A = numeric.identity(4);
	for (let i = 0; i < 4; i++) A[i][i] = 0;
	let v = [[]];
	v[0].length = 30;
	let coord = [];
	coord.length = 30;

	for (let i = 0; i < parseInt(this.state.numSamples); i++) {
	  for (let i = 0; i < 30; i++) v[0][i] = factor[0][i]*ref.randg();
	  let step = numeric.add(numeric.transpose(numeric.dot(eigen.eigenvectors, numeric.transpose(v))), step0);
	  //console.log(step);
	  for (let i = 0; i < 30; i++) coord[i] = step[0][i];
	  ref.get30Coordinates(coord, this.state.tetramer);
	  let bpstep = ref.getStep();
	  A = numeric.add(A, bpstep);
	}
	A = numeric.mul(1/parseFloat(this.state.numSamples), A);
	for (let i = 0; i < 2000; i++) {
	  A = numeric.dot(A, A);
	}
	alert("Persistence length in Angstroms: " + A[2][3]);*/
  }

  onMouseUp = ({clientX, clientY}) => {
    this.mount.removeEventListener('mousemove', this.onMouseMove);
    this.mount.removeEventListener('mouseup', this.onMouseUp);
    this.setState({
      lastTranslateX: this.state.translateX,
      lastTranslateY: this.state.translateY,
      startX: 0,
      startY: 0,
      isDragging: false
    });
  };

  onMouseDown = ({clientX,  clientY}) => {
    this.mount.addEventListener('mousemove', this.onMouseMove);
    this.mount.addEventListener('mouseup', this.onMouseUp);
    this.setState({
      startX: clientX,
      startY: clientY,
      isDragging: true
    });
  };
 
  moveEigen() {
    let cov = this.state.cov;
    let F = numeric.inv(cov);
    let eigen = ref.jeigen(F);
    let demo = this.state.mean;
    demo = numeric.add(demo, numeric.mul(parseFloat(this.state.stepBy) * (this.state.useScale ? 1/Math.sqrt(eigen.eigenvalues[parseInt(this.state.modeNum)]) : 1.0), eigen.eigenvectors[parseInt(this.state.modeNum)]));
    let points = ref.get30Coordinates(demo, this.state.tetramer, true);
    let PhoW = ref.getAtomSets().atoms[1][1];
    let PhoC = ref.getAtomSets().atoms[4][1];
    ref3.getBasePlanes(ref.getAtomSets());
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
	this.setState({
		phoW: pw, 
		phoC: pc, 
		parameters3: ref3.getParameters(),
		pdbText: ref.writePDB()
	});
  }

  componentWillUnmount() {
    this.mount.removeEventListener('mousemove', this.onMouseMove);
    this.mount.removeEventListener('mouseup', this.onMouseUp);
    this.mount.removeEventListener('mousedown', this.onMouseDdown);
  }

  componentDidMount() {
    let cov = this.state.cov;
    let F = numeric.inv(cov);
    let eigen = ref.jeigen(F);
    let eigenlist = [];
    eigen.eigenvalues.forEach((item, index) => {
	  let val = {index: index, value: item};
	  eigenlist.push(val);
    });
    eigenlist.sort((a, b) => (a.value > b.value ? 1 : -1));
    let modeN = eigenlist[0].index + '';
    this.setState({evalue: eigen.eigenvalues[parseInt(modeN)],
    eigenvalues: eigen.eigenvalues, eigenvectors: eigen.eigenvectors, eigenlist: eigenlist, modeNum: modeN});
    let demo = this.state.mean;
    let points = ref.get30Coordinates(demo, this.state.tetramer, true);
    let PhoW = ref.getAtomSets().atoms[1][1];
    let PhoC = ref.getAtomSets().atoms[4][1];
    ref3.getBasePlanes(ref.getAtomSets());
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
	this.setState({
		phoW: pw, 
		phoC: pc, 
		parameters3: ref3.getParameters(),
		pdbText: ref.writePDB()
	});
    document.body.style.backgroundColor = "white";
    this.scene = new trois.Scene();
    this.scene.background = new trois.Color( 0xffffff );
    this.camera = new trois.OrthographicCamera(-20, 20, 20, -20, -40, 500);
    this.scene.add(this.camera);
    this.renderer = new trois.WebGLRenderer();
    this.renderer.setSize(600, 600);
//    this.renderer.setSize( window.innerWidth, window.innerHeight);
    this.mount.appendChild(this.renderer.domElement);
    this.mount.addEventListener('mousedown', this.onMouseDown);

    this.camera.zoom = 1.0;

    let avgx=0.0, avgy=0.0, avgz=0.0;
    for (let i = 0; i < points.length; i++) {
      avgx += points[i].x;
      avgy += points[i].y;
      avgz += points[i].z;
    }
    avgx /= points.length; avgy /= points.length; avgz /= points.length;
    let spheres = [];
    let val = '';
    points.forEach((atom) => {
      //  console.log(atom);
        let geometry = new trois.SphereGeometry(1.4, 10, 10, 10);
        let color = 0xAA00AA;
        if (atom.name.charAt(0) === 'O') color = 0xCC0000;
        if (atom.name.charAt(0) === 'C') color = 0x777777;
        if (atom.name.charAt(0) === 'N') color = 0x0000CC;

	if (atom.name === 'P') {
		let px = atom.x - ref.getMidBasis()[0][3];
		let py = atom.y - ref.getMidBasis()[1][3];
		let pz = atom.z - ref.getMidBasis()[2][3];
		let pxx = px*ref.getMidBasis()[0][0] + py*ref.getMidBasis()[1][0] + pz*ref.getMidBasis()[2][0];
                let pyy = px*ref.getMidBasis()[0][1] + py*ref.getMidBasis()[1][1] + pz*ref.getMidBasis()[2][1];
                let pzz = px*ref.getMidBasis()[0][2] + py*ref.getMidBasis()[1][2] + pz*ref.getMidBasis()[2][2];
		val += pxx + ", " + pyy + ", " + pzz + " ... ";
	}

        let material = new trois.MeshLambertMaterial({color: color});
        let sphere = new trois.Mesh(geometry, material);
        sphere.position.set(atom.x-avgx, atom.y-avgy, atom.z-avgz);
        this.scene.add(sphere);
        spheres.push(sphere);
    });
    this.setState({ midVal: val });
    const light1 = new trois.PointLight( 0xffffff, 5, 200 );
    light1.position.set(0, 0, 50);
    this.scene.add(light1);
    const light2 = new trois.PointLight( 0xffffff, 5, 200 );
    light2.position.set(0, 0, -50);
    this.scene.add(light2);

    let angle = 0;
    //his.camera.position.z = 0.0;
    this.camera.position.x = 20.0*Math.cos(angle);
    this.camera.position.z = 20.0*Math.sin(angle);

    this.camera.lookAt(0.0, 0.0, 0.0);

    let scale = 1/Math.sqrt(eigen.eigenvalues[parseInt(modeN)]);
    let x = 0;
    let dir = 1;
    const animate = () => {
      if (this.state.useScale) scale = 1.0/Math.sqrt(eigen.eigenvalues[parseInt(this.state.modeNum)]);
      else scale = 1.0;
      let icnew = numeric.add(demo, numeric.mul(x/50*scale, eigen.eigenvectors[parseInt(this.state.modeNum)]));
      let points = ref.get30Coordinates(icnew, this.state.tetramer, false);
      let i = 0;
      points.forEach((atom) => {
        spheres[i].position.set(atom.x-avgx, atom.y-avgy, atom.z-avgz);
        i++;
      })
      let theta = this.state.translateX / 100.0; let phi = -this.state.translateY / 100.0;
      this.camera.position.x = 30.0*Math.cos(theta)*Math.cos(phi);
      this.camera.position.z = 30.0*Math.sin(theta)*Math.cos(phi);
      this.camera.position.y = 30.0*Math.sin(phi);
      this.camera.lookAt(0.0, 0.0, 0.0);
      requestAnimationFrame(animate);
      this.renderer.render(this.scene, this.camera);
      x += dir;
      if (x === 100) {
        dir = -1;
      }
      if (x === -100) {
        dir = 1;
      }
    }
    console.log("animating.....");
    animate();
  }

  inputChanged = (event) => {
      const target = event.target;
      const value = target.type === 'checkbox' ? target.checked : target.value;
      const name = target.name;
      this.setState({
        [name] : value
      }, () => {
        if (name === 'stepBy' && !isNaN(parseFloat(value))) this.moveEigen();
        
        if ((name === 'useScale' || name === 'modeNum') && parseFloat(this.state.stepBy) !== 0 && !isNaN(parseFloat(this.state.stepBy))) this.moveEigen();
      });
      if (name === 'modeNum')
        this.setState({evalue: this.state.eigenvalues[parseInt(value)]});

  }


  render() {
  
    const valtitles = ["basepair 1", "phosphateW", "bp step", "phosphateC", "basepair 2"];
    const val3titles = ["basepair 1", "bp step", "basepair 2"];
    let meanvals = [this.state.mean.slice(0,6), this.state.mean.slice(6, 12), this.state.mean.slice(12, 18), this.state.mean.slice(18, 24), this.state.mean.slice(24, 30)];
    
    
    return (<>
    	<div className="data-window">
    Number of Gaussian samples<input type="number" step="1" value={this.state.numSamples} name="numSamples" onChange={this.inputChanged}/><button disabled={this.state.processing} onClick={this.persistenceLength}>Calculate Statistics</button><br/>
    {this.state.showPlot ? <><Plot layout={this.state.layout1} data={this.state.data1} /><br/><Plot layout={this.state.layout5} data={this.state.data5} /><br/><Plot layout={this.state.layout6} data={this.state.data6} /><br/><Plot layout={this.state.layout2} data={this.state.data2} /><br/><Plot layout={this.state.layout8} data={this.state.data8} /><br/><Plot layout={this.state.layout3} data={this.state.data3} /><br/><Plot layout={this.state.layout7} data={this.state.data7} /><br/><Plot layout={this.state.layout4} data={this.state.data4} /><br/></> : null}
	<b>Mean EPFL/Curves+ state:</b><input type="checkbox" name="showIC" checked={this.state.showIC} onChange={this.inputChanged}/>Show original internal (rad/5) form<table><tbody><tr><td><b>PARAM</b></td><td>rot1 (°)</td><td>rot2 (°)</td><td>rot3 (°)</td><td>trans1 (Å)</td><td>trans2 (Å)</td><td>trans3 (Å)</td></tr>{meanvals.map((value,index) => (<tr><td>{valtitles[index]}</td><td>{(value[0]*(this.state.showIC ? 1 : 11.4591559*ref.scale(value[0], value[1], value[2]))).toFixed(5)}</td><td>{(value[1]*(this.state.showIC ? 1 : 11.4591559*ref.scale(value[0], value[1], value[2]))).toFixed(5)}</td><td>{(value[2]*(this.state.showIC ? 1 : 11.4591559*ref.scale(value[1], value[2], value[3]))).toFixed(5)}</td><td>{value[3].toFixed(5)}</td><td>{(value[4]).toFixed(5)}</td><td>{(value[5]).toFixed(5)}</td></tr>))}</tbody></table>
	<b>Phosphates:</b><table><thead style={{border: "0 none"}}><tr><td></td><td>xp (Å)</td><td>yp (Å)</td><td>zp (Å)</td></tr></thead><tbody><tr><td>watson</td><td>{this.state.phoW[0].toFixed(5)}</td><td>{this.state.phoW[1].toFixed(5)}</td><td>{this.state.phoW[2].toFixed(5)}</td></tr><tr><td>crick</td><td>{this.state.phoC[0].toFixed(5)}</td><td>{this.state.phoC[1].toFixed(5)}</td><td>{this.state.phoC[2].toFixed(5)}</td></tr></tbody></table>
	{this.state.eigenvalues.length > 0 ? <><select value={this.state.modeNum} name="modeNum" onChange={this.inputChanged}>
        {this.state.eigenlist.map((value) => (
            <option key={value.value} value={value.index}>{value.index + ": eigenvalue: " + value.value}</option>
        ))}</select> <input type="checkbox" name="useScale" onChange={this.inputChanged} checked={this.state.useScale}/> Scale by standard deviation<br/>Step Along Eigenvector<input type="number" step="0.1" onChange={this.inputChanged} name="stepBy" value={this.state.stepBy} /></>
        : null}
	<div style={{border: "solid", margin: "auto", width:"600px", height: "600px", justifyContent: "center", textAlign: "center"}} ref={ref => (this.mount = ref)}></div>
	  {this.state.eigenvectors.length > 0 ? <Eigenvector vector={this.state.eigenvectors[parseInt(this.state.modeNum)]}/> : null}<br/><br/>
	  	<b>3DNA state:</b><table><tbody><tr><td><b>PARAM</b></td><td>rot1 (°)</td><td>rot2 (°)</td><td>rot3 (°)</td><td>trans1 (Å)</td><td>trans2 (Å)</td><td>trans3 (Å)</td></tr>{this.state.parameters3.map((value,index) => (<tr><td>{val3titles[index]}</td><td>{value[0].toFixed(5)}</td><td>{value[1].toFixed(5)}</td><td>{value[2].toFixed(5)}</td><td>{value[3].toFixed(5)}</td><td>{value[4].toFixed(5)}</td><td>{value[5].toFixed(5)}</td></tr>))}</tbody></table><br/><br/>
	  <u><b>PDB text file</b></u>
	  <pre>{this.state.pdbText}</pre><br/>
	  </div>
	  </>);
  }

}

export default DThree;
