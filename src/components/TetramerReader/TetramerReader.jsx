import React from 'react';
import DThree from '../DThree/DThree.jsx';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';
import Plot from 'react-plotly.js';
import "./TetramerReader.css";
import numeric from 'numeric';

class TetramerReader extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      loaded: "",
      selected: "1",
      tetramer: "",
      mean: {},
      covariance: {},
      items: [],
      check1: false, check2: false, check3: false, check4: false,
      plotItem: "twist",
      analysis: [[], [], [], []],
      analysisC: [[], [], [], []],
      engaged: false,
      layout: null,
      data: null
    };
  }

  engage = () => {
    if (this.state.tetramer.length > 3) {
      this.setState({
        engaged: false
      }, () => {
        this.setState({
          engaged: true
        });
      });
    }
  }
  
  analyzeData() {
  	let source = [
  		["mean_TX3_3S_C1_cg_136.txt", "cov_TX3_3S_C1_cg_136.txt"],
  		["mean_cgDNA_136.txt", "cov_cgDNA_136.txt"],
  		["MD_mean_bstj_cgf_136.txt", "MD_cov_bstj_cgf_136.txt"],
  		["mean_TX2_3S_C1_cg_136.txt", "cov_TX2_3S_C1_cg_136.txt"]
  	];
  	let currentAnalysis = [[], [], [], []];
  	let currentAnalysisC = [[], [], [], []];
  	let count = 0;
    let count2 = 0;
  	source.forEach((set, index) => {
  	let dataSet = [];
  	let dataSet2 = [];
	fetch(set[0]).then((r) => r.text()).then(text => {
      let lines = text.split("\n");
      let sumTwist = 0;
      let sumTwist2 = 0;
      lines.forEach(line => {
        let litems = line.split(" ");
        if (litems.length < 31) return;

        let steps = [];
        for (let i = 0; i < 30; i++) {
          steps.push(parseFloat(litems[i+1]));
        }
        
        if (ref.complement(litems[0]) === litems[0]) sumTwist += steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        else sumTwist += 2*steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        
        let points = ref.get30Coordinates(steps, litems[0], true);
	    let PhoW = ref.getAtomSets().atoms[1][1];
	    let PhoC = ref.getAtomSets().atoms[4][1];
		ref3.getBasePlanes(ref.getAtomSets());
        let stepParameters = ref3.getParameters();

        if (ref.complement(litems[0]) === litems[0]) sumTwist2 += stepParameters[1][2];
        else sumTwist2 += 2*stepParameters[1][2];

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
		let dataItem = [litems[0], stepParameters[1][2], pw[0], pw[1], pw[2], pc[0], pc[1], pc[2], steps[14]*11.4591559*ref.scale(steps[12],steps[13],steps[14]), steps[12]*11.4591559*ref.scale(steps[12],steps[13],steps[14]), steps[13]*11.4591559*ref.scale(steps[12],steps[13],steps[14])];
        dataSet.push(dataItem);
      });
      console.log(set[0] + " twist = " + sumTwist/256.0);
      console.log(set[0] + " twist 3DNA = " + sumTwist2/256.0);
      //console.log(JSON.stringify(dataSet));
      currentAnalysis[index] = dataSet;
	  count++;
	  if (count === 4) {
		this.setState({analysis: currentAnalysis});	  	
	  }
    });

    let lastitem = "";
	fetch(set[1]).then((r) => r.text()).then(text => {
      let lines = text.split("\n");

      let covT = [];
      lines.forEach(line => {
        let litems = line.split(" ");
        if (litems.length < 31) return;

        if (litems[0] !== lastitem) {
          if (lastitem !== "") {
          	let eig = ref.jeigen(covT);
          	let prod = 1;
          	for (let i = 0; i < 30; i++) prod *= Math.sqrt(eig.eigenvalues[i]);
          	//console.log(Math.sqrt(numeric.det(covT)) + " " + prod);

          	let dataItem = [lastitem, prod];
          	dataSet2.push(dataItem);
          	covT = [];
          }
          lastitem = litems[0];
        }
        let row = [];
        for (let i = 0; i < 30; i++) {
          row.push(parseFloat(litems[i+1]));
        }
        covT.push(row);
      });
      //console.log(covT);
  	  let eig = ref.jeigen(covT);
  	  let prod = 1;
  	  for (let i = 0; i < 30; i++) prod *= Math.sqrt(eig.eigenvalues[i]);
  	  let dataItem = [lastitem, prod];
  	  dataSet2.push(dataItem);
  	  currentAnalysisC[index] = dataSet2;
      count2++;
      //console.log(count2);
      if (count2 === 4) {
      	this.setState({
        	analysisC: currentAnalysisC
        });
      }
    });
    
  	});
  	//console.log(currentAnalysis);
  }


  clicked = () => {
  	this.analyzeData();
  }


  inputChanged = (e) => {
    const value = e.target.type === "checkbox" ? e.target.checked : e.target.value;
    this.setState({
      [e.target.name]: value
    });
    if (e.target.name === "tetramer") this.engage();
  };

  loadSet = () => {
    let f1 = "", f2 = "";
    let itemsT = [];
    let covarianceT = {};
    let meanT = {};
    let field = 'unknown';
    if (parseInt(this.state.selected) === 1) {
      f1 = process.env.PUBLIC_URL+"mean_TX3_3S_C1_cg_136.txt";
      f2 = process.env.PUBLIC_URL+"cov_TX3_3S_C1_cg_136.txt";
      field = "TX3-S3-C1 X-ray analysis of protein-DNA structures"
    }
    if (parseInt(this.state.selected) === 2) {
      f1 = process.env.PUBLIC_URL+"mean_cgDNA_136.txt";
      f2 = process.env.PUBLIC_URL+"cov_cgDNA_136.txt";
      field = "cgDNA+ model trained on MD palindrome library with AMBER parmbsc1"
    }
    if (parseInt(this.state.selected) === 3) {
      f1 = process.env.PUBLIC_URL+"MD_mean_bstj_cgf_136.txt";
      f2 = process.env.PUBLIC_URL+"MD_cov_bstj_cgf_136.txt";
      field = "MD data directly from palindrome library with AMBER parmbsc1"
    }
    if (parseInt(this.state.selected) === 4) {
      f1 = process.env.PUBLIC_URL+"mean_TX2_3S_C1_cg_136.txt";
      f2 = process.env.PUBLIC_URL+"cov_TX2_3S_C1_cg_136.txt";
	  field = "TX2-S3-C1 X-ray analysis of protein-DNA structures";
    }
    fetch(f1).then((r) => r.text()).then(text => {
      let lines = text.split("\n");
      lines.forEach(line => {
        let litems = line.split(" ");
        if (litems.length < 31) return;
        itemsT.push(litems[0]);
        let steps = [];
        for (let i = 0; i < 30; i++) {
          steps.push(parseFloat(litems[i+1]));
        }
        
        meanT[litems[0]] = steps;
      });
    });
    fetch(f2).then((r) => r.text()).then(text => {
      let lines = text.split("\n");
      let lastitem = "";
      lines.forEach(line => {
        let litems = line.split(" ");
        if (litems.length < 31) return;
        if (!itemsT.includes(litems[0])) {
          alert("Step does not exist! " + litems[0]);
        }
        if (litems[0] !== lastitem) {
          covarianceT[litems[0]] = [];
          lastitem = litems[0];
        }
        let row = [];
        for (let i = 0; i < 30; i++) {
          row.push(parseFloat(litems[i+1]));
        }
        covarianceT[litems[0]].push(row);
      });
      this.setState({
        tetramer: itemsT[0],
        items: itemsT,
        covariance: covarianceT,
        mean: meanT,
        loaded: field
      }, () => this.forceUpdate());
      console.log(meanT);
      console.log(covarianceT);
      console.log(itemsT);
    });

  }
  
  makePlot = () => {
  	let data = [];
  	let names = ["TX3", "cgDNA+", "MD", "TX2"];
  	let useData = [false, false, false, false];
  	let labels = {
  		"twist": "Twist in degrees",
  		"S": "log(configurational volume) (Å^(15/2)*(rad/5)^(15/2))",
  		"roll": "Roll in degrees",
  		"tilt": "Tilt in degrees",
  		"px": "xp of Watson phosphate (Å)",
  		"px2": "xp of Crick phosphate (Å)",
  		"py": "yp of Watson phosphate (Å)",
  		"py2": "yp of Crick phosphate (Å)",
  		"pz": "zp of Watson phosphate (Å)",
  		"pz2": "zp of Crick phosphate (Å)"
  	}
  	if (this.state.check1) useData[0] = true;
  	if (this.state.check2) useData[1] = true;
  	if (this.state.check3) useData[2] = true;
  	if (this.state.check4) useData[3] = true;
  	if (!this.state.check1 && !this.state.check2 && !this.state.check3 && !this.state.check4) return;
  	let min = 1000, max = -1000;

  	for (let i = 0; i < 4; i++) {

  		if (useData[i]) {
	  	  	let xvalues = [];
 	 		let yvalues = [];
  			let max_diff = 0;
  			if (this.state.plotItem === "S") {
  				//console.log(this.state.analysisC[i]);
  				this.state.analysisC[i].forEach((value, index) => {
  					xvalues.push(value[0]);
  					yvalues.push(Math.log(value[1])/Math.log(10));
  					if (Math.log(value[1])/Math.log(10) > max) max = Math.log(value[1])/Math.log(10);
  					if (Math.log(value[1])/Math.log(10) < min) min = Math.log(value[1])/Math.log(10);
  				});
  			} else this.state.analysis[i].forEach((value, index) => {

  				xvalues.push(value[0]);
  				if (this.state.plotItem === "tilt") {
  					yvalues.push(value[9]);
  					if (value[9] > max) max = value[9];
  					if (value[9] < min) min = value[9];
  				}
  				else if (this.state.plotItem === "roll") {
  					yvalues.push(value[10]);
  					if (value[10] > max) max = value[10];
  					if (value[10] < min) min = value[10];
  				}
  				else if (this.state.plotItem === "twist") {
  					yvalues.push(value[1]);
  					if (value[1] > max) max = value[1];
  					if (value[1] < min) min = value[1];
  					xvalues.push(value[0]);
  					yvalues.push(value[8]);
  					if (Math.abs(value[8]-value[1]) > max_diff) max_diff = Math.abs(value[8] - value[1]);
  				}
  				else if (this.state.plotItem === "px") {
  					yvalues.push(value[2]);
  					if (value[2] > max) max = value[2];
  					if (value[2] < min) min = value[2];
  				}
  				else if (this.state.plotItem === "py") {
  					yvalues.push(value[3]);
  					if (value[3] > max) max = value[3];
  					if (value[3] < min) min = value[3];
  				}
  				else if (this.state.plotItem === "pz") {
  					yvalues.push(value[4]);
  					if (value[4] > max) max = value[4];
  					if (value[4] < min) min = value[4];
  				}
  				else if (this.state.plotItem === "px2") {
  					yvalues.push(value[5]);
  					if (value[5] > max) max = value[5];
  					if (value[5] < min) min = value[5];
  				}
  				else if (this.state.plotItem === "py2") {
  					yvalues.push(value[6]);
  					if (value[6] > max) max = value[6];
  					if (value[6] < min) min = value[6];
  				}
  				else if (this.state.plotItem === "pz2") {
  					yvalues.push(value[7]);
  					if (value[7] > max) max = value[7];
  					if (value[7] < min) min = value[7];
  				}
  			});
  			data.push({
            	x: xvalues,
            	y: yvalues,
            	name: names[i],
            	type: 'scatter',
            	mode: 'markers+lines',
            	marker: {
            	    size: 5,
            	    line: {
                  		width: 2
                	}
            	}
            });
            //console.log(max_diff);
  		}
  	}

  	let dtick = 0.1;
  	if (max-min > 2) dtick = 0.25;
  	if (max-min > 4) dtick = 0.5;
  	if (max-min > 10) dtick = 1;
  	this.setState({data: data,
            layout: {
              title: "Comparison of Data Sets",
              width: 1080,
              height: 800,
              margin: {
                b: 150
              },
              font: {
                size: 18
              },
	          xaxis: {
                range: [0, 136],
                //autorange: 'reversed',
                showgrid: true,
                showline: true,
                gridwidth: 2,
                gridcolor: '#777777',
                //title: 'number of days before today',
                automargin: true,
                linewidth: 1
              },
              yaxis: {
                range: [min-(max-min)/20.0, max+(max-min)/20.0],
                showgrid: true,
                showline: true,
                dtick: dtick,
                gridwidth: 2,
                gridcolor: '#777777',
                title: labels[this.state.plotItem],
                linewidth: 1
              }
			}
	});
  }

  clearData = () => {
  	this.setState({analysis: [[],[],[],[]], layout: null});
  }

  render() {
    return (<div>
    	  <div className="analyze-data">
    	  <button onClick={this.clicked}>Collect Data</button>{this.state.analysis[0].length > 0 ? <><input type="checkbox" name="check1" checked={this.state.check1} onChange={this.inputChanged}/>TX-3 <input type="checkbox" name="check2" checked={this.state.check2} onChange={this.inputChanged}/>cgDNA+ <input type="checkbox" name="check3" checked={this.state.check3} onChange={this.inputChanged}/>MD <input type="checkbox" name="check4" checked={this.state.check4} onChange={this.inputChanged}/>TX-2 <select name="plotItem" value={this.state.plotItem} onChange={this.inputChanged}>
    	  	<option id="twist" value="twist" key="1">Plot Twist</option>
    	  	<option id="px" value="px" key="2">Phosphate xp</option>
    	  	<option id="py" value="py" key="3">Phosphate yp</option>
    	  	<option id="pz" value="pz" key="4">Phosphate zp</option>
    	  	<option id="px2" value="px2" key="5">Phosphate xp (Crick)</option>
    	  	<option id="py2" value="py2" key="6">Phosphate yp (Crick)</option>
    	  	<option id="pz2" value="pz2" key="7">Phosphate zp (Crick)</option>
    	  	<option id="tilt" value="tilt" key="8">Tilt</option>
    	  	<option id="roll" value="roll" key="9">Roll</option>
    	  	<option id="S" value="S" key="10">S (config. volume)</option>
    	  </select><button onClick={this.makePlot}>Create Plot!</button></> : null}
    	  {this.state.layout ? <Plot layout={this.state.layout} data={this.state.data} /> : null} <button onClick={this.clearData}>Clear Data and Plots</button>
    	  </div>
          <br/>{this.state.loaded}<br/><br/>
          <select name="selected" value={this.state.selected} onChange={this.inputChanged}>
            <option id="1" value="1" key="1">TX3 X-ray</option>
            <option id="2" value="2" key="2">cgDNA+ 136</option>
            <option id="3" value="3" key="3">MD</option>
            <option id="4" value="4" key="4">TX2 X-ray</option>
          </select>
          <button onClick={this.loadSet}>Load Dataset</button>
          <select name="tetramer" value={this.state.tetramer} onChange={this.inputChanged}>
            {this.state.items.length > 0 ? this.state.items.map((name) => (
              <option key={name} value={name}>{name}</option>
            )) : null}
          </select>
          <button onClick={this.engage}>Run</button>
          {this.state.engaged ? <DThree mean={this.state.mean[this.state.tetramer]} cov={this.state.covariance[this.state.tetramer]} tetramer={this.state.tetramer}/> : null}
          </div>
        )
  }
}

export default TetramerReader;

