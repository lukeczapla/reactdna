import React from 'react';
import DThree from '../DThree/DThree.jsx';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';
import Plot from 'react-plotly.js';
import "./TetramerReader.css";
import numeric from 'numeric';
import PersistenceAllWorker from './PersistenceAll.worker.js';

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
      processing: false,
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
        
        if (ref.complementTetramer(litems[0]) === litems[0]) sumTwist += steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        else sumTwist += 2*steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        
        let points = ref.get30Coordinates(steps, litems[0], true);
	    let PhoW = ref.getAtomSets().atoms[1][1];
	    let PhoC = ref.getAtomSets().atoms[4][1];
		ref3.getBasePlanes(ref.getAtomSets());
        let stepParameters = ref3.getParameters();

        if (ref.complementTetramer(litems[0]) === litems[0]) sumTwist2 += stepParameters[1][2];
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
  	//console.log(currentAnalysis); console.log(currentAnalyisC);
  }

  persistenceAll = () => {
    if (Object.keys(this.state.mean).length < 130) {
    	alert("No model is currently loaded");
    	return;
    }
    this.setState({processing: true});
    let worker = new PersistenceAllWorker();
    worker.postMessage(this.state);
    worker.addEventListener("message", (e) => {
    	let r = e.data;
    	console.log(r);
    	alert("Global mean persistence length = " + r.P); 
		document.getElementById("resultPersistence").innerHTML = JSON.stringify(r.tetramerData);
    	 
    	this.setState({processing: false});
    });
  }

  meanTwist = () => {
  	let seq = this.state.sequence.toUpperCase();
	let total = 0;
	if (seq.length < 4) {
		alert("Bad sequence input");
		return;
	}
  	for (let i = 1; i < seq.length-2; i++) {
  	  let tetramer = seq.substring(i-1, i+3);
  	  console.log(tetramer);
  	  console.log(this.state.mean);
  	  if (this.state.mean[tetramer] === undefined) {
  	    tetramer = ref.complementTetramer(tetramer);
  	    if (this.state.mean[tetramer] === undefined) {
  	    	alert("Bad sequence input");
  	    	return;
  	    }
  	  }
  	  console.log(tetramer);
  	  let state = this.state.mean[tetramer];
  	  let twist = state[14]*11.4591559*ref.scale2(state[14]);
  	  total += twist;
  	}
  	alert("Mean twist in degrees " + total/(seq.length-3));
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
  		"pz2": "zp of Crick phosphate (Å)",
  		"persistence": "Persistence length of repeat (Å)"
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

  				if (this.state.plotItem !== "persistence") xvalues.push(value[0]);
  				if (this.state.plotItem === "persistence") {
  					xvalues.push(persistenceData[i][index][0]);
  					yvalues.push(persistenceData[i][index][1]);
  					if (persistenceData[i][index][1] > max) max = persistenceData[i][index][1];
  					if (persistenceData[i][index][1] < min) min = persistenceData[i][index][1];
  				}
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
  	if (max-min > 100) dtick = 50;
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
    	  	<option id="persistence" value="persistence" key="persistence">Persistence Length</option>
    	  </select><button onClick={this.makePlot}>Create Plot!</button></> : null}
    	  {this.state.layout ? <Plot layout={this.state.layout} data={this.state.data} /> : null} <button onClick={this.clearData}>Clear Data and Plots</button>
    	  <br/>
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
          {this.state.engaged ? <><br/>Sequence: <input type="text" name="sequence" value={this.state.sequence} onChange={this.inputChanged}/>
    	  <button onClick={this.meanTwist}>Calculate Twist</button>
    	  <br/><button disabled={this.state.processing} onClick={this.persistenceAll}>All Persistence Lengths</button>
    	  <div id="resultPersistence"></div></> : null}
          </div>
        )
  }
}

let persistenceData = [[["ATAT",416.7287745403619],["TTAA",320.601004387896],["CTAG",332.1513061028878],["GTAC",244.47104511431863],["ATAC",551.080734575081],["GTAG",233.68351211585534],["TTAG",489.27197297442797],["CTAT",387.8126577206288],["TTAC",395.72987047177315],["TTAT",451.7536107818927],["GCGC",449.8893730326439],["CCGG",584.4018552032469],["TCGA",244.6704452059255],["ACGT",377.65058548676376],["CCGC",436.18061417448297],["CCGT",598.3960955645102],["ACGC",494.28567753029546],["TCGT",579.8289486281102],["TCGC",665.4382040298821],["TCGG",682.6275081904797],["GCAA",551.4750191015536],["ACAA",528.7515201551973],["GCAG",402.09084770831765],["CCAA",187.73455896068737],["ACAG",495.46535796442663],["TCAA",350.74258424361744],["GCAC",231.38057548929666],["ACAT",526.8421919677572],["CCAG",469.63780007478266],["GCAT",1002.8613053760585],["CCAC",569.364278391309],["ACAC",565.0820105473265],["TCAG",799.3374031952591],["CCAT",622.2158757666708],["TCAC",404.95276981953555],["TCAT",679.2092362414846],["GAAA",810.8107575430605],["AAAA",765.5705697678369],["CAAA",625.5807146267397],["GAAG",765.2473848079637],["GAAC",890.27206507643],["TAAA",509.299103127015],["AAAG",744.0426430231677],["AAAC",661.5587685673415],["AAAT",1036.9908362580325],["TAAC",693.3433297426034],["GAAT",682.339027775706],["TAAG",533.2617730646128],["CAAG",739.4700348740901],["CAAC",586.1480305486532],["TAAT",450.5505159802673],["CAAT",959.3586799380713],["AAGA",973.3649817456575],["GAGA",644.8106046306476],["CAGA",986.5786493223319],["GAGG",895.1665408453484],["AAGG",661.6496972143235],["TAGA",493.8361966940738],["AAGC",775.9387178031968],["CAGC",1004.9807674574662],["AAGT",730.6717522904349],["GAGC",488.4460346919492],["TAGG",700.4640101007274],["GAGT",720.7963472170151],["CAGG",736.2475123662331],["TAGT",646.8324276022374],["TAGC",501.4204675702566],["CAGT",911.2844152693207],["GGAA",559.9444948385385],["AGAA",870.0060419889872],["AGAT",539.8568984946007],["TGAA",591.9182482186252],["AGAC",669.5756926699519],["AGAG",686.0580287193301],["GGAC",1190.3619098803845],["CGAA",963.3624603161996],["GGAG",551.3394434575715],["CGAG",524.8928324581894],["GGAT",764.0948346742987],["TGAC",835.8442982656694],["TGAT",938.5536629376852],["TGAG",568.7703941926347],["CGAT",610.5847936419043],["CGAC",1045.0607204057767],["AGGA",1006.2722721461321],["AGGG",634.3662438751642],["GGGA",583.2215110356864],["AGGT",617.9815537836457],["AGGC",603.8734715385459],["CGGG",487.37277004252684],["CGGA",672.5088731945871],["GGGT",599.4960894091419],["GGGG",406.9001044964782],["GGGC",629.1879083783908],["TGGA",631.296762339112],["TGGG",799.3069389158634],["CGGT",438.5492364076849],["TGGC",810.6509163007644],["CGGC",580.5334829495255],["TGGT",660.5234766642524],["AGCT",620.8216542251791],["TGCA",532.4136360020335],["CGCG",652.056279307327],["GGCC",767.7017406175797],["CGCC",606.1403228657902],["TGCG",1558.589907004925],["TGCC",622.560121137302],["GGCT",987.0001109442076],["CGCT",980.2749971897907],["TGCT",550.9330254694288],["CATG",1085.2449122149887],["TATA",615.2290729971702],["GATC",1111.9269854158993],["AATT",1246.7630237457142],["CATA",898.2103634073467],["GATT",1051.2701378918407],["CATC",1460.2956047608293],["TATC",882.2622010620862],["TATT",1038.116472371246],["CATT",1028.6786240167487],["AACC",1076.3246242284804],["GACA",938.0470275610215],["AACG",732.8892113327843],["GACG",1317.5380827832162],["AACA",788.6331849739886],["CACG",943.6919288833128],["CACA",803.8383153176028],["TACA",690.0256471698434],["GACC",1102.1621535982817],["TACC",662.0232658743964],["CACC",796.9406197249397],["TACG",583.7541137539739],["GACT",1073.8842613408492],["AACT",821.9882102612412],["TACT",784.7447543592313],["CACT",1125.511768861923]],
[["ATAT",319.01807502103156],["TTAA",258.5958865482089],["CTAG",262.42626864447476],["GTAC",278.2300340691773],["ATAC",296.5388050936105],["GTAG",273.16687191800986],["TTAG",261.9582246895267],["CTAT",291.4353404972914],["TTAC",270.1113653823751],["TTAT",292.6063162410017],["GCGC",337.76347004028594],["CCGG",320.0705654163649],["TCGA",304.6442727464671],["ACGT",374.4615674498627],["CCGC",331.56447709451464],["CCGT",349.9757233597836],["ACGC",355.7356354990418],["TCGT",347.3141201343603],["TCGC",328.9674320157178],["TCGG",313.2023869625751],["GCAA",293.4460511528658],["ACAA",319.52532548915036],["GCAG",305.74857286641225],["CCAA",296.85943492898565],["ACAG",324.8818188630538],["TCAA",280.91699877355],["GCAC",312.4343496740689],["ACAT",355.48473230943137],["CCAG",308.61692840131275],["GCAT",338.7566271396223],["CCAC",311.55204988600036],["ACAC",330.16367514693934],["TCAG",294.00744113026576],["CCAT",339.8307626933394],["TCAC",302.1505640949217],["TCAT",332.2756803439738],["GAAA",575.5619657907635],["AAAA",585.8984220780127],["CAAA",565.7303024577756],["GAAG",568.9512747510088],["GAAC",574.926680766056],["TAAA",560.995611619632],["AAAG",581.6967679094639],["AAAC",594.7334399283153],["AAAT",618.5766397633167],["TAAC",566.3585383364838],["GAAT",612.9407253066057],["TAAG",561.2965805995406],["CAAG",565.5775911041126],["CAAC",583.7568133483078],["TAAT",607.3497837549295],["CAAT",610.389341394848],["AAGA",562.4307134134765],["GAGA",555.0923050949444],["CAGA",545.2448230564959],["GAGG",555.1700194823644],["AAGG",564.2429505874886],["TAGA",528.1614676906479],["AAGC",581.0595439063605],["CAGC",564.0126946289845],["AAGT",594.5532728694432],["GAGC",565.3906980544634],["TAGG",534.3185613650663],["GAGT",586.0966666680689],["CAGG",552.1336803854049],["TAGT",566.6598867370199],["TAGC",551.9678272938861],["CAGT",579.7815181060735],["GGAA",543.215339298579],["AGAA",582.3256489456774],["AGAT",600.649487274424],["TGAA",515.5236416459887],["AGAC",576.7111414979742],["AGAG",583.9952559996626],["GGAC",539.8330017978265],["CGAA",511.5475277148904],["GGAG",545.6756783014397],["CGAG",504.84946380713956],["GGAT",563.3608570472974],["TGAC",517.5090086572372],["TGAT",528.5255228027331],["TGAG",508.0871937962104],["CGAT",524.4596975327979],["CGAC",500.5419524586712],["AGGA",531.2128256153396],["AGGG",545.0012818599229],["GGGA",510.7937617173989],["AGGT",537.4418382806194],["AGGC",553.5657886989119],["CGGG",495.21138177733894],["CGGA",489.604947793191],["GGGT",525.7216685660823],["GGGG",515.5676144693069],["GGGC",528.6991687007843],["TGGA",499.1156458659142],["TGGG",503.27249592275297],["CGGT",496.74927812391553],["TGGC",513.5126894795753],["CGGC",499.1135616143504],["TGGT",510.42860317245186],["AGCT",681.7790475263874],["TGCA",601.1591785589178],["CGCG",594.74676165827],["GGCC",618.3327879520142],["CGCC",611.8821334612957],["TGCG",602.4985366769291],["TGCC",609.5420524594351],["GGCT",652.7535232461132],["CGCT",639.8040219967252],["TGCT",645.9566460222836],["CATG",720.3417157347606],["TATA",678.9145718647762],["GATC",757.0680226670268],["AATT",795.7926747341292],["CATA",710.2604332665387],["GATT",775.2391094053665],["CATC",731.4972005486474],["TATC",716.7993543056108],["TATT",755.0576115096196],["CATT",753.2347481972793],["AACC",656.3079742020261],["GACA",591.6320264755162],["AACG",615.1288440808344],["GACG",594.4011957056633],["AACA",620.4878556790381],["CACG",584.3237004833879],["CACA",588.9512564250203],["TACA",569.4250915467006],["GACC",647.031156422624],["TACC",612.3307218206193],["CACC",607.4386381232266],["TACG",589.0592865560096],["GACT",661.6493577173081],["AACT",697.6485659677437],["TACT",649.198065721569],["CACT",658.6604642528158]],
[["ATAT",285.74660478265065],["TTAA",228.7286414145127],["CTAG",244.48384641097059],["GTAC",262.2361148036629],["ATAC",262.28047003433653],["GTAG",250.19208129033717],["TTAG",224.88435203457132],["CTAT",253.91500530567023],["TTAC",231.02622532765886],["TTAT",248.267255381557],["GCGC",312.75992718209096],["CCGG",291.96156466053856],["TCGA",273.04348490344086],["ACGT",283.34550034850406],["CCGC",295.9501675953176],["CCGT",294.0977091257742],["ACGC",296.3668540509002],["TCGT",279.0566503869202],["TCGC",307.5249534910934],["TCGG",305.0028899261057],["GCAA",286.679707399371],["ACAA",261.8646795482499],["GCAG",265.9869308292012],["CCAA",252.74489794400733],["ACAG",273.59745532579944],["TCAA",248.31983891621394],["GCAC",279.25813486437517],["ACAT",318.782678980545],["CCAG",277.8466143200518],["GCAT",273.3272643213285],["CCAC",278.6642909576198],["ACAC",276.8335696127179],["TCAG",287.94176989923017],["CCAT",292.3185858910762],["TCAC",286.0864815505253],["TCAT",271.01537581364653],["GAAA",521.8529064699128],["AAAA",544.8658198674274],["CAAA",541.8140557920942],["GAAG",571.0689909016535],["GAAC",576.8588618596398],["TAAA",505.5924283638281],["AAAG",505.99542216554727],["AAAC",543.5603647625461],["AAAT",625.7021977549485],["TAAC",499.95727830399954],["GAAT",601.9992112758417],["TAAG",483.2187125356188],["CAAG",505.52202183369263],["CAAC",529.1226073025804],["TAAT",520.0112122582742],["CAAT",563.7960416762467],["AAGA",478.5074643666814],["GAGA",467.45654147445805],["CAGA",470.129925315031],["GAGG",485.5901285830671],["AAGG",491.60183292180506],["TAGA",432.0137187244915],["AAGC",530.1178433173247],["CAGC",566.8858683949205],["AAGT",563.2266441901625],["GAGC",542.1614111086941],["TAGG",513.320673222814],["GAGT",524.7830348375738],["CAGG",485.4144585269071],["TAGT",517.0186314193295],["TAGC",477.33875343826963],["CAGT",533.9771487511213],["GGAA",485.3035948237757],["AGAA",452.2103780307637],["AGAT",504.36709951582515],["TGAA",473.8816285691705],["AGAC",481.5515379424191],["AGAG",461.928419611233],["GGAC",494.63890401873726],["CGAA",487.0343365385525],["GGAG",467.6667307032311],["CGAG",482.93925013119775],["GGAT",481.0443455099914],["TGAC",517.3077407237997],["TGAT",525.9170886993361],["TGAG",485.0822929858733],["CGAT",473.009153293215],["CGAC",497.3145826009677],["AGGA",412.59462612757847],["AGGG",528.5834381578433],["GGGA",432.37454766445626],["AGGT",519.061904030197],["AGGC",513.8845544187095],["CGGG",444.78435080357553],["CGGA",396.069550363611],["GGGT",534.8305577156473],["GGGG",449.3302504455544],["GGGC",492.7122185085492],["TGGA",425.3677129099936],["TGGG",417.22022532452525],["CGGT",410.65065946796346],["TGGC",440.3256103014547],["CGGC",457.4583126493936],["TGGT",456.877742986078],["AGCT",765.8822255812353],["TGCA",544.5401332634416],["CGCG",524.8070048711364],["GGCC",653.0053146583001],["CGCC",586.0231313587641],["TGCG",522.9924664465469],["TGCC",552.8785044951691],["GGCT",695.5805459379733],["CGCT",551.4011521607433],["TGCT",577.2511889901037],["CATG",616.2677414218713],["TATA",540.2694700205044],["GATC",681.451507274799],["AATT",882.336260506322],["CATA",559.6569551438288],["GATT",778.9457686019559],["CATC",631.5792963764388],["TATC",696.9041429985645],["TATT",680.5483768439095],["CATT",751.6026171455189],["AACC",632.1354177629571],["GACA",605.6072105689013],["AACG",557.968906680198],["GACG",575.6013476975946],["AACA",567.4434499137556],["CACG",485.3884007892327],["CACA",503.9012679807989],["TACA",493.76116057132487],["GACC",631.8244357846878],["TACC",547.428861523219],["CACC",535.4332911774025],["TACG",512.8344744744508],["GACT",650.3504995010642],["AACT",674.5499020830605],["TACT",535.706794928708],["CACT",550.7181038118384]],
[["ATAT",538.2814835126306],["TTAA",438.20021647072235],["CTAG",334.3622410576122],["GTAC",230.9841001360835],["ATAC",595.6700275729867],["GTAG",291.41912123801023],["TTAG",442.09986509659814],["CTAT",357.6109185924768],["TTAC",352.86275139473196],["TTAT",480.7067044983743],["GCGC",454.48319539856976],["CCGG",604.622924378564],["TCGA",236.44088893487228],["ACGT",464.30608166923895],["CCGC",582.6750628550093],["CCGT",641.0598900637225],["ACGC",674.3849939314354],["TCGT",647.917412116562],["TCGC",797.0888261815262],["TCGG",754.3550139464058],["GCAA",565.0843851797512],["ACAA",567.5991328277993],["GCAG",466.1057680411427],["CCAA",222.71790241353938],["ACAG",622.1681182305862],["TCAA",402.7498895020523],["GCAC",447.1126989501901],["ACAT",553.9082657570913],["CCAG",524.5236588226609],["GCAT",1343.9505298574838],["CCAC",581.1845364794075],["ACAC",698.3635826625464],["TCAG",915.7374462357736],["CCAT",709.1805803103173],["TCAC",443.39174192379437],["TCAT",739.5691964360647],["GAAA",884.3423600452617],["AAAA",901.5531223305215],["CAAA",676.6975738138582],["GAAG",859.3039610197187],["GAAC",903.1085100091143],["TAAA",388.0897607643798],["AAAG",789.2449409385388],["AAAC",704.5651287406755],["AAAT",1149.4859650789347],["TAAC",799.6608239066338],["GAAT",845.9757916090725],["TAAG",738.7412909566469],["CAAG",982.9215986771495],["CAAC",661.7891384881921],["TAAT",577.9768727259824],["CAAT",1082.291355013455],["AAGA",1151.7281405189963],["GAGA",1059.4423175775319],["CAGA",1119.7447100915413],["GAGG",973.1039037052999],["AAGG",989.3944749667387],["TAGA",572.3097331003154],["AAGC",877.2178621999734],["CAGC",1390.5203557804898],["AAGT",794.083505956009],["GAGC",510.94449965942874],["TAGG",783.4050178793673],["GAGT",740.0345899070284],["CAGG",845.8332830258147],["TAGT",701.194648451981],["TAGC",637.5384116555512],["CAGT",1215.1535739583032],["GGAA",595.6393240714175],["AGAA",1069.4160446837793],["AGAT",495.7688308514843],["TGAA",640.4342034119588],["AGAC",837.608837135415],["AGAG",724.1501897002079],["GGAC",1312.3613072225432],["CGAA",1126.0499065650597],["GGAG",652.2238238226007],["CGAG",628.4988695779516],["GGAT",885.4933114742408],["TGAC",957.0689239120956],["TGAT",1022.6453983660837],["TGAG",675.3023627378543],["CGAT",627.2897800892255],["CGAC",1271.4523907479697],["AGGA",1073.9449273412981],["AGGG",780.1577378889814],["GGGA",824.1545595784463],["AGGT",604.033053715964],["AGGC",640.2804794041816],["CGGG",724.157971522406],["CGGA",791.5864032904893],["GGGT",653.9763037884629],["GGGG",749.8746795379108],["GGGC",627.047894209208],["TGGA",778.0084791544706],["TGGG",838.5146796762775],["CGGT",482.50924911801116],["TGGC",848.9902701275769],["CGGC",705.5883256579581],["TGGT",691.3864178323963],["AGCT",667.743543368885],["TGCA",879.01439304592],["CGCG",688.7749620796347],["GGCC",735.3276045529273],["CGCC",646.4217437423315],["TGCG",2219.9158556870107],["TGCC",763.140254515059],["GGCT",1155.4896737146544],["CGCT",1199.9229607508662],["TGCT",762.8179024469657],["CATG",1365.5149634077482],["TATA",504.5798844974184],["GATC",1129.3288605693363],["AATT",1454.0284431905536],["CATA",1003.3362873458394],["GATT",1103.5702207732063],["CATC",1688.5800705951488],["TATC",923.5172341033989],["TATT",1000.636098871376],["CATT",1142.90623941199],["AACC",1128.2831973173963],["GACA",1096.433966212035],["AACG",990.7090210846092],["GACG",1521.1242397291132],["AACA",818.6197620936623],["CACG",1060.5342763716956],["CACA",934.1554030256639],["TACA",937.8710827361855],["GACC",1205.877839358905],["TACC",662.6787944308423],["CACC",969.6725703492782],["TACG",841.0295439306683],["GACT",1181.1638206632333],["AACT",931.9308391869702],["TACT",1007.4585674546461],["CACT",1180.0709098542739]]];


export default TetramerReader;

