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
  	  let twist = state[14]*11.4591559*ref.scale(state[12], state[13], state[14]);
  	  total += twist;
  	}
  	alert("Mean twist in degrees " + total/(seq.length-2));
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

let persistenceData = [[["ATAT",425.93898115005675],["TTAA",326.92917059457034],["CTAG",322.65636325429057],["GTAC",238.5412365280469],["ATAC",529.1344533217449],["GTAG",230.85198471813507],["TTAG",483.5967379912384],["CTAT",404.78368232723443],["TTAC",403.51630679894294],["TTAT",448.08873332050723],["GCGC",461.4170368491386],["CCGG",575.8364024750817],["TCGA",243.32840469354733],["ACGT",362.90500317684075],["CCGC",446.52920307734365],["CCGT",601.6816293212323],["ACGC",483.8551515796741],["TCGT",582.6806790378151],["TCGC",686.1841100227429],["TCGG",656.8357825737598],["GCAA",581.048966111188],["ACAA",521.8948922361886],["GCAG",408.2210914746432],["CCAA",186.64006776920687],["ACAG",518.0027340899206],["TCAA",349.7017192751979],["GCAC",234.5869932644138],["ACAT",525.3865248923676],["CCAG",470.44799338733276],["GCAT",985.5666753688057],["CCAC",566.5674937247613],["ACAC",568.5028347824799],["TCAG",817.4506458750325],["CCAT",640.2571634244347],["TCAC",397.70154716987895],["TCAT",668.6339226690906],["GAAA",788.8997645619855],["AAAA",743.0016903917902],["CAAA",636.4015699013413],["GAAG",780.6168114022666],["GAAC",873.9156134306882],["TAAA",543.393463075244],["AAAG",771.7759502050017],["AAAC",668.5192990213443],["AAAT",1086.9491752574033],["TAAC",687.3968758473237],["GAAT",675.0687981675072],["TAAG",534.3876987857864],["CAAG",761.459744589082],["CAAC",580.5157063311533],["TAAT",442.86316607161024],["CAAT",950.859538199493],["AAGA",1027.5989609799424],["GAGA",659.6190514914109],["CAGA",979.0994844873939],["GAGG",888.2117612410802],["AAGG",658.707132148642],["TAGA",495.5659940206004],["AAGC",779.8092464543433],["CAGC",992.1270099686078],["AAGT",760.3869004643511],["GAGC",484.74948924580184],["TAGG",685.4262642463077],["GAGT",716.5776826557067],["CAGG",762.9172480595151],["TAGT",625.496123383388],["TAGC",508.2609442011775],["CAGT",905.7505885809304],["GGAA",543.2903086395506],["AGAA",829.2622471002754],["AGAT",560.9798977432645],["TGAA",595.6758962339688],["AGAC",676.6304596003424],["AGAG",681.4963531554599],["GGAC",1230.123467443683],["CGAA",951.6961285243837],["GGAG",543.6707001930193],["CGAG",536.0224868673077],["GGAT",768.5152577764495],["TGAC",869.6867311287718],["TGAT",910.8155518774063],["TGAG",541.1843537194519],["CGAT",615.0140206204477],["CGAC",1060.2099433079452],["AGGA",1000.9145256909925],["AGGG",609.3786101584908],["GGGA",606.7150023445789],["AGGT",634.9811999003107],["AGGC",589.6411179005497],["CGGG",501.3413801728817],["CGGA",654.620706886379],["GGGT",599.812046761991],["GGGG",407.3889417458898],["GGGC",644.28556338577],["TGGA",635.2927309843265],["TGGG",806.0348474160946],["CGGT",444.9187859982743],["TGGC",822.5199024552654],["CGGC",581.0812767373926],["TGGT",648.170670979672],["AGCT",606.857409172752],["TGCA",537.4640311926432],["CGCG",658.3649893548009],["GGCC",779.0439082901316],["CGCC",598.7501314491084],["TGCG",1540.5061720733556],["TGCC",628.4343173943978],["GGCT",1006.9472502554221],["CGCT",976.796821835213],["TGCT",548.9128978062785],["CATG",1129.0697003057326],["TATA",624.8714385478396],["GATC",1114.6208190883751],["AATT",1249.6234376448876],["CATA",889.9018549935611],["GATT",1027.7684320811725],["CATC",1516.311877692751],["TATC",890.2533412731832],["TATT",1053.2485549872722],["CATT",1042.9804181077154],["AACC",1064.5063722555092],["GACA",930.9780783853836],["AACG",743.4282005018255],["GACG",1315.8180734018695],["AACA",784.5238248722588],["CACG",942.4920642129209],["CACA",805.2114044999586],["TACA",681.4241938150626],["GACC",1127.4326805595492],["TACC",674.3646574634587],["CACC",816.6789906076485],["TACG",574.4843851574681],["GACT",1116.5187245619636],["AACT",864.1787646820358],["TACT",788.2462985663572],["CACT",1129.0194607796507]],
[["ATAT",326.4867479730945],["TTAA",253.36832697962333],["CTAG",268.7870255118676],["GTAC",285.1666818827286],["ATAC",293.2919233027393],["GTAG",268.5373488446511],["TTAG",269.38522172063034],["CTAT",295.41041822442673],["TTAC",262.1281368771223],["TTAT",310.9540668146072],["GCGC",339.8470252734544],["CCGG",317.02111821541666],["TCGA",297.20763188938616],["ACGT",368.9224778586616],["CCGC",325.4099603107123],["CCGT",343.7090843139967],["ACGC",355.2266967002073],["TCGT",336.2071728047369],["TCGC",322.5245811617884],["TCGG",316.7784797514125],["GCAA",294.82091090257376],["ACAA",319.3678154731538],["GCAG",308.94146640432723],["CCAA",294.8187079997643],["ACAG",335.76058391664674],["TCAA",282.4607910106764],["GCAC",302.20086320120595],["ACAT",342.9765471193588],["CCAG",307.241297132686],["GCAT",330.8210809333155],["CCAC",320.51129749186833],["ACAC",321.53213019376216],["TCAG",295.8730828576098],["CCAT",345.5637357681693],["TCAC",312.5639361140569],["TCAT",318.8076520044289],["GAAA",569.2615462641772],["AAAA",575.7966168436393],["CAAA",578.2063632102163],["GAAG",577.334051730207],["GAAC",582.6019753877141],["TAAA",580.8737587185585],["AAAG",582.7522708495077],["AAAC",569.1432402191911],["AAAT",618.0501729894856],["TAAC",555.5882341105261],["GAAT",602.6239595467542],["TAAG",560.9928422698114],["CAAG",556.8913325968286],["CAAC",581.8446870532136],["TAAT",611.8858308851213],["CAAT",587.4470761500509],["AAGA",559.0532318674636],["GAGA",546.0454355747664],["CAGA",540.7158454266287],["GAGG",549.8540315862763],["AAGG",572.1226761602459],["TAGA",534.6642319501984],["AAGC",578.8956566597601],["CAGC",555.7989330099911],["AAGT",606.426652709143],["GAGC",572.6459159428384],["TAGG",523.8736065100426],["GAGT",584.9800728133463],["CAGG",588.5926833227311],["TAGT",587.1441611609778],["TAGC",556.3086829943538],["CAGT",584.3162059096511],["GGAA",559.1474666424288],["AGAA",588.0279496123472],["AGAT",573.429092962544],["TGAA",514.9418354943755],["AGAC",576.0102757793557],["AGAG",576.6544353259709],["GGAC",555.3121400132874],["CGAA",521.0341574015893],["GGAG",527.7904363583067],["CGAG",502.85762942567504],["GGAT",575.5414934919031],["TGAC",506.95732712909665],["TGAT",520.4136980303704],["TGAG",521.8723606285404],["CGAT",527.5096144325335],["CGAC",490.9712096303611],["AGGA",534.7242959461367],["AGGG",538.2024600052999],["GGGA",495.428285786664],["AGGT",564.4144935707337],["AGGC",534.1380141411405],["CGGG",481.8852176049788],["CGGA",468.7108711878456],["GGGT",537.8517706109654],["GGGG",536.9129175236595],["GGGC",535.4593086778119],["TGGA",482.0330149421424],["TGGG",517.6483615154416],["CGGT",505.0468235604081],["TGGC",500.4115039126043],["CGGC",505.53382326894155],["TGGT",498.256737965092],["AGCT",696.4288384015608],["TGCA",586.3650131590653],["CGCG",581.0781530088236],["GGCC",593.2511532245425],["CGCC",625.1488213245414],["TGCG",613.5370012038055],["TGCC",611.4294579631563],["GGCT",654.6830318135139],["CGCT",641.4329382541613],["TGCT",630.451710031909],["CATG",731.3560263398417],["TATA",712.1627693523399],["GATC",757.4867256260643],["AATT",794.5292948194681],["CATA",704.241558683763],["GATT",750.3625285362886],["CATC",748.8650556931598],["TATC",740.3782090460303],["TATT",752.1825438126048],["CATT",799.0100807282668],["AACC",665.4436764405859],["GACA",594.2426218024542],["AACG",622.5506799916452],["GACG",597.3605103258493],["AACA",596.191398998363],["CACG",599.222484918893],["CACA",574.7053520471042],["TACA",580.7822759395063],["GACC",636.2452234532392],["TACC",586.703006034546],["CACC",620.2953588288462],["TACG",606.0483702287462],["GACT",689.9657600565192],["AACT",710.4403675847647],["TACT",625.1841390559822],["CACT",635.895396451646]],
[["ATAT",268.9449212281795],["TTAA",228.96210456678835],["CTAG",244.55069940114424],["GTAC",281.70824735760516],["ATAC",274.5091227368108],["GTAG",253.0063397220033],["TTAG",220.02538421785138],["CTAT",257.04789364153135],["TTAC",223.3428943030924],["TTAT",248.11751495977876],["GCGC",307.56744004815096],["CCGG",283.89768298129144],["TCGA",276.8125272553022],["ACGT",283.7222932112805],["CCGC",301.7124850829594],["CCGT",291.86580357324505],["ACGC",299.5859357631433],["TCGT",285.18401769141906],["TCGC",307.6682522270348],["TCGG",310.2263362890348],["GCAA",287.13286105333304],["ACAA",245.43033729001067],["GCAG",266.51902403539793],["CCAA",242.93028085750825],["ACAG",275.683951836052],["TCAA",249.65980046899102],["GCAC",276.91510585396213],["ACAT",323.4929589412155],["CCAG",273.4954865965792],["GCAT",278.16801673068966],["CCAC",294.65956060131936],["ACAC",280.1297822650488],["TCAG",288.629788643687],["CCAT",295.9059730634788],["TCAC",280.6063646503778],["TCAT",275.2781543524274],["GAAA",520.9552773040484],["AAAA",537.8873792523922],["CAAA",540.2994063391162],["GAAG",564.776885708879],["GAAC",581.1711128075256],["TAAA",513.0323245335585],["AAAG",507.87803493019936],["AAAC",543.4848119670518],["AAAT",630.056472794869],["TAAC",539.983199942082],["GAAT",606.3684617628114],["TAAG",487.8345700337187],["CAAG",497.30024185151984],["CAAC",524.3505946606663],["TAAT",495.1051824294994],["CAAT",542.6661286757859],["AAGA",498.48685283702696],["GAGA",474.19116600829705],["CAGA",476.6958715004074],["GAGG",478.33503267650025],["AAGG",496.85159604786014],["TAGA",426.982938470066],["AAGC",532.9807172288339],["CAGC",579.3338303230585],["AAGT",542.2135829745944],["GAGC",539.1119349303901],["TAGG",524.1003527758535],["GAGT",520.3650103292996],["CAGG",490.40398639009743],["TAGT",516.5047772606646],["TAGC",479.3666078155934],["CAGT",498.6853616238388],["GGAA",482.1570600385695],["AGAA",441.7850266007751],["AGAT",520.0022012728793],["TGAA",472.7051349389742],["AGAC",485.89885472493603],["AGAG",453.1994152223279],["GGAC",512.6997281948537],["CGAA",474.2612011726804],["GGAG",463.7563694645009],["CGAG",480.5150532050123],["GGAT",505.4974594491417],["TGAC",534.3308402220944],["TGAT",518.8785601273489],["TGAG",486.99207269123343],["CGAT",462.5020598597132],["CGAC",495.6953773266352],["AGGA",425.33340138466065],["AGGG",521.3710526091614],["GGGA",433.0325172772776],["AGGT",517.1647059142049],["AGGC",535.0162133221123],["CGGG",421.99483812458004],["CGGA",397.138904085087],["GGGT",545.1711764278335],["GGGG",429.05331847987827],["GGGC",495.8739462926691],["TGGA",422.6600561567077],["TGGG",409.62860503966004],["CGGT",418.33929163641386],["TGGC",455.3920725797652],["CGGC",456.62829247623046],["TGGT",433.08430167871416],["AGCT",765.41464705954],["TGCA",567.0426841426578],["CGCG",513.4039115571998],["GGCC",648.4034319810311],["CGCC",580.646905135433],["TGCG",523.8595550526506],["TGCC",545.3663092310819],["GGCT",704.4974561973696],["CGCT",543.1483018903308],["TGCT",603.8427960176073],["CATG",616.2133124930212],["TATA",579.2542086902055],["GATC",702.4057346106475],["AATT",884.267826737832],["CATA",580.9380482349126],["GATT",792.9845754076421],["CATC",637.7932404176391],["TATC",694.6654247625612],["TATT",676.4992545306911],["CATT",736.7564811583072],["AACC",645.4954211148353],["GACA",586.4612493358064],["AACG",574.9773800349764],["GACG",566.0807357782122],["AACA",558.1675663468383],["CACG",489.5254071484027],["CACA",507.065784236008],["TACA",496.39776674634317],["GACC",631.1857594657258],["TACC",555.8838035936827],["CACC",540.2067429909484],["TACG",500.0825977806356],["GACT",644.8784077912243],["AACT",682.3132747009471],["TACT",540.1422110954521],["CACT",546.2291612112095]],
[["ATAT",537.4279921406265],["TTAA",427.40690829272444],["CTAG",347.46382520372396],["GTAC",219.80342912393007],["ATAC",604.2052265724153],["GTAG",309.10676930675726],["TTAG",438.6718367362833],["CTAT",373.26007124779215],["TTAC",353.0216953350303],["TTAT",485.8190873174683],["GCGC",461.2143861747551],["CCGG",608.2839060030636],["TCGA",245.9769191669172],["ACGT",458.36299255350997],["CCGC",591.7619032465163],["CCGT",637.1454035321711],["ACGC",668.6545888933351],["TCGT",662.053927539966],["TCGC",788.2535769937382],["TCGG",741.6550356793792],["GCAA",580.2366780420499],["ACAA",543.0548198957998],["GCAG",461.41825600213684],["CCAA",220.35573587312624],["ACAG",607.5631303035304],["TCAA",422.26988782976747],["GCAC",442.8909959702859],["ACAT",554.0798258353037],["CCAG",532.935900339684],["GCAT",1328.6642893619096],["CCAC",573.1861796936201],["ACAC",719.1034251315627],["TCAG",938.5737566395549],["CCAT",721.5419617507066],["TCAC",452.3669661461534],["TCAT",756.1234398712822],["GAAA",912.4293956336155],["AAAA",886.6442350787878],["CAAA",689.4861531333607],["GAAG",853.2076302592542],["GAAC",920.1969066192469],["TAAA",393.5040920033749],["AAAG",798.9380159584334],["AAAC",706.0027103315113],["AAAT",1097.1343024013934],["TAAC",801.7684281336741],["GAAT",842.7202090430396],["TAAG",770.5323174044947],["CAAG",1001.3584057625258],["CAAC",690.9550447937069],["TAAT",590.6565978117422],["CAAT",1133.4726610315142],["AAGA",1158.5394216906666],["GAGA",1088.285300031415],["CAGA",1147.8607097917823],["GAGG",949.1115918272972],["AAGG",1000.5408324633054],["TAGA",592.8549093408823],["AAGC",872.9143807459446],["CAGC",1346.3148118872587],["AAGT",809.4073134248143],["GAGC",514.9646262477745],["TAGG",741.6097992140703],["GAGT",753.4837350141839],["CAGG",855.0721517238804],["TAGT",693.6557186288694],["TAGC",633.7409753050382],["CAGT",1213.687506881904],["GGAA",583.5856798517055],["AGAA",1055.096264405084],["AGAT",523.4782641647554],["TGAA",627.2223893146678],["AGAC",850.5999730205076],["AGAG",711.455887188288],["GGAC",1350.4050267195162],["CGAA",1127.2769080068445],["GGAG",664.418108521027],["CGAG",634.056951411199],["GGAT",914.9265933906252],["TGAC",976.1652139353001],["TGAT",1010.7637574501999],["TGAG",674.3936685160421],["CGAT",632.7767623742087],["CGAC",1253.0399678898946],["AGGA",1051.02595067049],["AGGG",784.5486248392415],["GGGA",826.4267921095242],["AGGT",582.4763358429835],["AGGC",609.4173668488418],["CGGG",733.1969026081872],["CGGA",772.194181755357],["GGGT",637.242367664885],["GGGG",812.509370893757],["GGGC",630.7503478975971],["TGGA",779.8802569411283],["TGGG",853.0333600174819],["CGGT",502.6470628892272],["TGGC",887.3453495612],["CGGC",708.2636741016911],["TGGT",700.6723741849283],["AGCT",687.6680474288708],["TGCA",880.8086497691205],["CGCG",641.5799598651137],["GGCC",744.2367770322154],["CGCC",653.248234546969],["TGCG",2262.3242676018335],["TGCC",760.6819141363856],["GGCT",1146.0215647561515],["CGCT",1191.1743864216796],["TGCT",794.9140635005359],["CATG",1391.0807648953164],["TATA",517.4050609515792],["GATC",1190.78325304783],["AATT",1481.1568284297691],["CATA",1002.0350794399493],["GATT",1112.79712974776],["CATC",1673.9585867691635],["TATC",920.1174144817708],["TATT",1069.8251833703412],["CATT",1192.8267167359647],["AACC",1154.0896327768062],["GACA",1128.2634268234701],["AACG",1014.5100565252662],["GACG",1551.176520650037],["AACA",814.8147817735136],["CACG",1035.6454331213902],["CACA",923.4652527735141],["TACA",929.167091917696],["GACC",1135.0392140241092],["TACC",684.6392180725885],["CACC",945.8734080643328],["TACG",841.095555526342],["GACT",1127.3150728271019],["AACT",965.4167198576826],["TACT",996.5955590052674],["CACT",1172.3899537341695]]];


export default TetramerReader;

