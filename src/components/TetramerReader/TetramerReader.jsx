import React from 'react';
import DThree from '../DThree/DThree.jsx';
import * as ref from '../References/References.jsx';
import * as ref3 from '../References/References3DNA.jsx';
import Plot from 'react-plotly.js';
import Form from 'react-bootstrap/Form'
import InputGroup from 'react-bootstrap/InputGroup';
import Table from 'react-bootstrap/Table';
import Button from 'react-bootstrap/Button';
import SearchBar from '../SearchBar.js';
import 'bootstrap/dist/css/bootstrap.min.css';
import styles from './SearchBar.module.css';

import "./TetramerReader.css";
import numeric from 'numeric';
import PersistenceAllWorker from './PersistenceAll.worker.js';

class TetramerReader extends React.Component {


  constructor(props) {
    super(props);
    this.state = {
      loaded: "",
      analyzed: false,
      tetramerText: "",
      selected: "1",
      tetramer: "",
      mean: {},
      covariance: {},
      items: [],
      check1: false, check2: false, check3: false, check4: false,// check5: false,
      plotItem: "twist",
      analysis: [[], [], [], []],
      analysisC: [[], [], [], []],
      engaged: false,
      layout: null,
      processing: false,
      data: null,
      rahul: null,
      resultField: "",
      resultField2: "",
      suggestions: []
    };
    this.tetramerList = ["ATAT","TTAA","CTAG","GTAC","ATAC","GTAG","TTAG","CTAT","TTAC","TTAT","GCGC","CCGG","TCGA","ACGT","CCGC","CCGT","ACGC","TCGT","TCGC","TCGG","GCAA","ACAA","GCAG","CCAA","ACAG","TCAA","GCAC","ACAT","CCAG","GCAT","CCAC","ACAC","TCAG","CCAT","TCAC","TCAT","GAAA","AAAA","CAAA","GAAG","GAAC","TAAA","AAAG","AAAC","AAAT","TAAC","GAAT","TAAG","CAAG","CAAC","TAAT","CAAT","AAGA","GAGA","CAGA","GAGG","AAGG","TAGA","AAGC","CAGC","AAGT","GAGC","TAGG","GAGT","CAGG","TAGT","TAGC","CAGT","GGAA","AGAA","AGAT","TGAA","AGAC","AGAG","GGAC","CGAA","GGAG","CGAG","GGAT","TGAC","TGAT","TGAG","CGAT","CGAC","AGGA","AGGG","GGGA","AGGT","AGGC","CGGG","CGGA","GGGT","GGGG","GGGC","TGGA","TGGG","CGGT","TGGC","CGGC","TGGT","AGCT","TGCA","CGCG","GGCC","CGCC","TGCG","TGCC","GGCT","CGCT","TGCT","CATG","TATA","GATC","AATT","CATA","GATT","CATC","TATC","TATT","CATT","AACC","GACA","AACG","GACG","AACA","CACG","CACA","TACA","GACC","TACC","CACC","TACG","GACT","AACT","TACT","CACT"].sort();
  }


  engage = () => {
    if (this.state.tetramer.length > 3 && this.tetramerList.includes(this.state.tetramer)) {
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
  		["means_internal.txt", "cov_internal.txt"],
        ["means-XNE-EV2.txt", "cov-XNE-EV2.txt"]
  	];
  	let currentAnalysis = [[], [], [], []];
  	let currentAnalysisC = [[], [], [], []];
  	let count = 0;
    let count2 = 0;
    let svalues = [];
    this.setState({analyzed: true});
    /*fetch("1TGH.CURVE.coord").then(r => r.text()).then(text => {
        let seq = "ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT";
        let lines = text.split("\n");
        lines.forEach(line => {
            svalues.push(parseFloat(line));
        });
        let v30 = svalues.slice(0, 30);
        let pt = ref.get30Coordinates(v30, "A" + seq.substring(0, 2) + "T", true);
        ref3.getBasePlanes(ref.getAtomSets());
        let stepParameters2 = ref3.getParameters();
        this.state.resultField += stepParameters2[0][0] + " " + stepParameters2[0][1] + " " + stepParameters2[0][2] + " " + stepParameters2[0][3] + " " + stepParameters2[0][4] + " " + stepParameters2[0][5] + " 0 0 0 0 0 0<br/>";
        this.state.resultField += stepParameters2[2][0] + " " + stepParameters2[2][1] + " " + stepParameters2[2][2] + " " + stepParameters2[2][3] + " " + stepParameters2[2][4] + " " + stepParameters2[2][5] + " " + stepParameters2[1][0] + " " + stepParameters2[1][1] + " " + stepParameters2[1][2] + " " + stepParameters2[1][3] + " " + stepParameters2[1][4] + " " + stepParameters2[1][5] + "<br/>";
        this.state.resultField2 += (v30[0]*11.4591559*ref.scale(v30[0], v30[1], v30[2])) + " " + (v30[1]*11.4591559*ref.scale(v30[0], v30[1], v30[2])) + " " + (v30[2]*11.4591559*ref.scale(v30[0], v30[1], v30[2])) + " " + v30[3] + " " + v30[4] + " " + v30[5] + " " + "0.0 0.0 0.0 0.0 0.0 0.0" + "<br/>";
        this.state.resultField2 += (v30[24]*11.4591559*ref.scale(v30[24], v30[25], v30[26])) + " " + (v30[25]*11.4591559*ref.scale(v30[24], v30[25], v30[26])) + " " + (v30[26]*11.4591559*ref.scale(v30[25], v30[26], v30[7])) + " " + v30[27] + " " + v30[28] + " " + v30[29] + " " + (v30[12]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + (v30[13]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + (v30[14]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + v30[15] + " " + v30[16] + " " + v30[17] + "<br/>";

        let chain = ref.writePDB(true, true, true);
        for (let i = 1; i < 12; i++) {
            v30 = svalues.slice(24 * i, 24 * i + 30);

            pt = ref.get30Coordinates(v30, "A" + seq.substring(0, 2) + "T", true);
            ref3.getBasePlanes(ref.getAtomSets());
            stepParameters2 = ref3.getParameters();
            this.state.resultField += stepParameters2[2][0] + " " + stepParameters2[2][1] + " " + stepParameters2[2][2] + " " + stepParameters2[2][3] + " " + stepParameters2[2][4] + " " + stepParameters2[2][5] + " " + stepParameters2[1][0] + " " + stepParameters2[1][1] + " " + stepParameters2[1][2] + " " + stepParameters2[1][3] + " " + stepParameters2[1][4] + " " + stepParameters2[1][5] + "<br/>";
            this.state.resultField2 += (v30[24]*11.4591559*ref.scale(v30[24], v30[25], v30[26])) + " " + (v30[25]*11.4591559*ref.scale(v30[24], v30[25], v30[26])) + " " + (v30[26]*11.4591559*ref.scale(v30[25], v30[26], v30[7])) + " " + v30[27] + " " + v30[28] + " " + v30[29] + " " + (v30[12]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + (v30[13]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + (v30[14]*11.4591559*ref.scale(v30[12], v30[13], v30[14])) + " " + v30[15] + " " + v30[16] + " " + v30[17] +  "<br/>";

            //let pt = ref.get30Coordinates(v30, "A" + seq.substring(i, i+2) + "T", true, ref.getA());
            //if (i < 144) chain += ref.writePDB(true, true, true);
            //else chain += ref.writePDB(true, false, true);
            //ref3.getBasePlanes(ref.getAtomSets());
            //let stepParameters2 = ref3.getParameters();
            //this.state.resultField += stepParameters2[2][3] + " " + stepParameters2[2][4] + " " + stepParameters2[2][5] + " " + stepParameters2[2][0] + " " + stepParameters2[2][1] + " " + stepParameters2[2][2] + " " + stepParameters2[1][3] + " " + stepParameters2[1][4] + " " + stepParameters2[1][5] + " " + stepParameters2[1][0] + " " + stepParameters2[1][1] + " " + stepParameters2[1][2] + "<br/>";
        }
        //this.state.resultField2 = resultField2;
        ref.setLastAtom(1);
        ref.setLastRes(1);
        let Acopy = JSON.parse(JSON.stringify(ref.getA()));
        Acopy[0][2] = -Acopy[0][2];
        Acopy[1][2] = -Acopy[1][2];
        Acopy[2][2] = -Acopy[2][2];
        Acopy[0][1] = -Acopy[0][1];
        Acopy[1][1] = -Acopy[1][1];
        Acopy[2][1] = -Acopy[2][1];
        ref.setA(Acopy);
        let seq2 = ref.complementLong(seq);
        chain = "";
        for (let i = 144; i >= 0; i--) {
            let v30 = ref.reverse30(svalues.slice(24 * i, 24 * i + 30));
            let pt = ref.get30Coordinates(v30, "A" + seq2.substring(i, i+2) + "T", true, ref.getA());
            if (i > 0) chain += ref.writePDB(true, true, true, 'B');
            else chain += ref.writePDB(true, false, true, 'B');
        }
        this.state.resultField2 += chain;


//        console.log("1KX5 3DNA-prime " + stepParameters2[1][3] + " " + stepParameters2[1][4] + " " + stepParameters2[1][5] + " " + stepParameters2[1][0] + " " + stepParameters2[1][1] + " " + stepParameters2[1][2]);
//        console.log("1KX5 3DNA-prime " + stepParameters2[0][3] + " " + stepParameters2[0][4] + " " + stepParameters2[0][5] + stepParameters2[0][0] + " " + stepParameters2[0][1] + " " + stepParameters2[0][2]);
    });*/
  	source.forEach((set, index) => {
  	let dataSet = [];
  	let dataSet2 = [];
	fetch(set[0]).then((r) => r.text()).then(text => {
      let lines = text.split("\n");
      let sumTwist = 0;
      let sumTwist2 = 0;
      let sumTilt = 0, sumTilt2 = 0, sumRoll = 0, sumRoll2 = 0;
      lines.forEach(line => {
        let litems = line.split(" ");
        if (litems.length < 31) return;

        let steps = [];
        for (let i = 0; i < 30; i++) {
          steps.push(parseFloat(litems[i+1]));
        }
        
	    //let v = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 3.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        //let pt = ref.get30Coordinates(v, litems[0], true);
        //ref3.getBasePlanes(ref.getAtomSets());
        //let stepParameters2 = ref3.getParameters();
        //console.log("tilt = " + stepParameters2[1][0] + ", roll = " + stepParameters2[1][1] + ", twist = " + stepParameters2[1][2]);

        if (ref.complementTetramer(litems[0]) === litems[0]) sumTwist += steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        else sumTwist += 2*steps[14]*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        if (ref.complementTetramer(litems[0]) === litems[0]) sumTilt += Math.abs(steps[12])*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        else sumTilt += 2*Math.abs(steps[12])*11.4591559*ref.scale(steps[12], steps[13], steps[14]);        
        if (ref.complementTetramer(litems[0]) === litems[0]) sumRoll += Math.abs(steps[13])*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
        else sumRoll += 2*Math.abs(steps[13])*11.4591559*ref.scale(steps[12], steps[13], steps[14]);
	    let points = ref.get30Coordinates(steps, litems[0], true);
	    let PhoW = ref.getAtomSets().atoms[1][1];
	    let PhoC = ref.getAtomSets().atoms[4][1];
        ref3.getBasePlanes(ref.getAtomSets());
        let stepParameters = ref3.getParameters();

        if (ref.complementTetramer(litems[0]) === litems[0]) sumTwist2 += stepParameters[1][2];
        else sumTwist2 += 2*stepParameters[1][2];
        if (ref.complementTetramer(litems[0]) === litems[0]) sumTilt2 += Math.abs(stepParameters[1][0]);
        else sumTilt2 += 2*Math.abs(stepParameters[1][0]);
        if (ref.complementTetramer(litems[0]) === litems[0]) sumRoll2 += Math.abs(stepParameters[1][1]);
        else sumRoll2 += 2*Math.abs(stepParameters[1][1]);

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
		let dataItem = [litems[0], stepParameters[1][2], pw[0], pw[1], pw[2], pc[0], pc[1], pc[2], steps[14]*11.4591559*ref.scale(steps[12],steps[13],steps[14]), steps[12]*11.4591559*ref.scale(steps[12],steps[13],steps[14]), steps[13]*11.4591559*ref.scale(steps[12],steps[13],steps[14]), stepParameters[1][0], stepParameters[1][1], stepParameters[1][3], stepParameters[1][4], stepParameters[1][5]];
	        dataSet.push(dataItem);
      });
      console.log(set[0] + " twist = " + sumTwist/256.0);
      console.log(set[0] + " twist 3DNA = " + sumTwist2/256.0);
      console.log(set[0] + " average |tilt|: " + sumTilt/256 + " and 3DNA " + sumTilt2/256);
      console.log(set[0] + " average |roll| " + sumRoll/256 + " and 3DNA " + sumRoll2/256);
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


  handleClear = () => {
      this.setState({suggestions: [], tetramerText: ""});
  }

  handleSearch = (value) => {
      if (!this.state.engaged) return;
      if (this.tetramerList.includes(value)) {
          this.setState({
              tetramer: value
          });
          this.engage();
      } else if (this.tetramerList.includes(ref.complementTetramer(value))) {
          this.setState({
              tetramer: ref.complementTetramer(value)
          });
          this.engage();
      }
  }


  inputChanged = (e) => {
    let value = e.target.type === "checkbox" ? e.target.checked : e.target.value;
    if (e.target.name === "tetramerText") {
      if (value.length > 4) value = value.slice(0,4);
	  value = value.toUpperCase();
	  if (this.state.mean[value] !== undefined) {
	  	this.setState({
	  	  tetramer: value
	  	});
	  	this.engage();
	  } else if (this.state.mean[ref.complementTetramer(value)] !== undefined) {
	  	value = ref.complementTetramer(value);
	  	this.setState({
	  	  tetramer: value
	  	});
	  	this.engage();	  	
	  }
    }
    this.setState({
      [e.target.name]: value
    });
    if (e.target.name === "tetramerText") {
    	
    }
    if (e.target.name === "tetramer") this.engage();
    if (e.target.name === "selected") {
        this.loadSet(value);
    }
  };

  loadSet = (value = null) => {
    let f1 = "", f2 = "";
    let itemsT = [];
    let covarianceT = {};
    let meanT = {};
    let field = 'unknown';
    if (value == null) value = this.state.selected;
    if (this.state.tetramer === "") this.setState({tetramer:"ATAT"});
    if (parseInt(value) === 1) {
      f1 = process.env.PUBLIC_URL+"mean_TX3_3S_C1_cg_136.txt";
      f2 = process.env.PUBLIC_URL+"cov_TX3_3S_C1_cg_136.txt";
      field = "TX3-S3-C1 X-ray analysis of protein-DNA structures"
    }
    if (parseInt(value) === 2) {
      f1 = process.env.PUBLIC_URL+"mean_cgDNA_136.txt";
      f2 = process.env.PUBLIC_URL+"cov_cgDNA_136.txt";
      field = "cgDNA+ model trained on MD palindrome library with AMBER parmbsc1"
    }
    if (parseInt(value) === 3) {
      f1 = process.env.PUBLIC_URL+"MD_mean_bstj_cgf_136.txt";
      f2 = process.env.PUBLIC_URL+"MD_cov_bstj_cgf_136.txt";
      field = "MD data directly from palindrome library with AMBER parmbsc1"
    }
    if (parseInt(value) === 4) {
      f1 = process.env.PUBLIC_URL+"means_internal.txt";
      f2 = process.env.PUBLIC_URL+"cov_internal.txt";
	  field = "X-ray/NMR/CryoEM analysis of protein-DNA structures";
    }
    if (parseInt(value) === 5) {
        f1 = process.env.PUBLIC_URL+"means-XNE-EV2.txt";
        f2 = process.env.PUBLIC_URL+"cov-XNE-EV2.txt";
        field = "X-ray/NMR/CryoEM EV/EV analysis of protein-DNA structures";
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
        tetramerText: "",
        items: itemsT,
        covariance: covarianceT,
        mean: meanT,
        loaded: field
      }, () => this.forceUpdate());
      this.engage();
      console.log(meanT);
      console.log(covarianceT);
      console.log(itemsT);
    });

  }
  
  makePlot = () => {
  	let data = [];
  	let names = ["TX3", "cgDNA+", "MD", "X/N/EM", "X/N/EM EV2"];
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
  		"persistence": "Persistence length of repeat (Å)",
  		"slide": "Slide (Å)",
  		"shift": "Shift (Å)",
  		"rise": "Rise (Å)"
  	}
  	if (this.state.check1) useData[0] = true;
  	if (this.state.check2) useData[1] = true;
  	if (this.state.check3) useData[2] = true;
  	if (this.state.check4) useData[3] = true;
    //if (this.state.check5) useData[4] = true;
  	if (!this.state.check1 && !this.state.check2 && !this.state.check3 && !this.state.check4 /*&& !this.state.check5*/) return;
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
					xvalues.push(value[0]);
					yvalues.push(value[11]);
  				}
  				else if (this.state.plotItem === "roll") {
  					yvalues.push(value[10]);
  					if (value[10] > max) max = value[10];
  					if (value[10] < min) min = value[10];
  					xvalues.push(value[0]);
					yvalues.push(value[12]);
				}
  				else if (this.state.plotItem === "twist") {
  					yvalues.push(value[1]);
  					if (value[1] > max) max = value[1];
  					if (value[1] < min) min = value[1];
  					xvalues.push(value[0]);
  					yvalues.push(value[8]);
  					if (Math.abs(value[8]-value[1]) > max_diff) max_diff = Math.abs(value[8] - value[1]);
  				}
                                else if (this.state.plotItem === "shift") {
                                        yvalues.push(value[13]);
                                        if (value[13] > max) max = value[13];
                                        if (value[13] < min) min = value[13];
                                }
                                else if (this.state.plotItem === "slide") {
                                        yvalues.push(value[14]);
                                        if (value[14] > max) max = value[14];
                                        if (value[14] < min) min = value[14];
                                }
                                else if (this.state.plotItem === "rise") {
                                        yvalues.push(value[15]);
                                        if (value[15] > max) max = value[15];
                                        if (value[15] < min) min = value[15];
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
  	this.setState({analysis: [[],[],[],[]], layout: null, analyzed: false});
  }


  binarySearch = (x, val, n= 1, eq= false) => {
        let results = [];
        let start = 0;
        let end = x.length;
        let mid = Math.trunc(x.length/2);
        n = Math.min(x[mid].length, val.length);
        if (!eq) {
            while (x[mid].substring(0,n) !== val.substring(0,n) && end-start > 1) {
                if (val.substring(0,n) < x[mid].substring(0,n)) end = mid;
                if (val.substring(0,n) > x[mid].substring(0,n)) start = mid;
                mid = Math.trunc((start + end)/2);
                n = Math.min(x[mid].length, val.length);
            }
        } else {
            while (x[mid] !== val && end-start > 1) {
                if (val < x[mid]) end = mid;
                if (val > x[mid]) start = mid;
                mid = Math.trunc((start + end)/2);
                n = Math.min(x[mid].length, val.length);
            }
        }

        if (!eq) {
            if (x[mid].substring(0,n) === val.substring(0,n)) results.push(x[mid]);
        } else {
            if (x[mid] === val) results.push(x[mid]);
        }
        let mid_o = mid;
        if (!eq) {
            while (mid > 0 && x[mid-1].substring(0,n) === val.substring(0,n)) {
                mid--;
                results.push(x[mid]);
            }
        } else {
            while (mid > 0 && x[mid-1] === val) {
                mid--;
                results.push(x[mid]);
            }
        }
        mid = mid_o;
        if (!eq) {
            while (mid < x.length-1 && x[mid+1].substring(0,n) === val.substring(0,n)) {
                mid++;
                results.push(x[mid]);
            }
        } else {
            while (mid < x.length-1 && x[mid+1] === val) {
                mid++;
                results.push(x[mid]);
            }
        }
        return {values: results, mid: mid_o};
    }

    handleClear = () => {
      this.setState({
          suggestions: []
      });
    }

    handleChange = (input) => {
        let og = input.toUpperCase();
        if (og.length === 0) return;
        this.setState({tetramerText: og});
        if (og.length < 1) return;
        let values1 = this.binarySearch(this.tetramerList, og).values;
        //console.log(og + " " + values1);
        let items = [];
        let exact = false;
        if (values1.length > 0)
            values1.forEach((x) => { items.push(x); if (og === x) exact = true; });
        //console.log(items);
        this.setState({suggestions: items});
        if (og.length === 4) {
            if (this.state.mean[og] !== undefined) {
                this.setState({
                    tetramer: og
                });
                this.engage();
            } else if (this.state.mean[ref.complementTetramer(og)] !== undefined) {
                let value = ref.complementTetramer(og);
                this.setState({
                    tetramer: value
                });
                this.engage();
            }
        }
        //if (input.indexOf(" ") === -1) input = "\"" + input + "\"";
    }

  render() {
    return (<div>
    	  <div className="analyze-data">
    	  <Button hidden={this.state.analyzed} onClick={this.clicked}>Collect Data</Button>{this.state.analysis[0].length > 0 ?
              <><Table striped bordered><thead></thead><tbody><tr>
                  <td><input type="checkbox" name="check1" checked={this.state.check1} onChange={this.inputChanged}/>TX3 X-ray protein-DNA</td>
                  <td><input type="checkbox" name="check2" checked={this.state.check2} onChange={this.inputChanged}/>cgDNA+ 136</td>
                  <td><input type="checkbox" name="check3" checked={this.state.check3} onChange={this.inputChanged}/>MD (AMBER parmbsc1) </td>
                  <td><input type="checkbox" name="check4" checked={this.state.check4} onChange={this.inputChanged}/>X+NMR+Cryo-EM</td>
                  {/*<td><input type="checkbox" name="check5" checked={this.state.check5} onChange={this.inputChanged}/>X+NMR+Cryo-EM EV2</td>*/}
              </tr></tbody></Table>
          <Form.Select name="plotItem" value={this.state.plotItem} onChange={this.inputChanged}>
    	  	<option id="twist" value="twist" key="1">Plot Twist</option>
    	  	<option id="px" value="px" key="2">Phosphate xp</option>
    	  	<option id="py" value="py" key="3">Phosphate yp</option>
    	  	<option id="pz" value="pz" key="4">Phosphate zp</option>
    	  	<option id="px2" value="px2" key="5">Phosphate xp (Crick)</option>
    	  	<option id="py2" value="py2" key="6">Phosphate yp (Crick)</option>
    	  	<option id="pz2" value="pz2" key="7">Phosphate zp (Crick)</option>
    	  	<option id="tilt" value="tilt" key="8">Tilt</option>
    	  	<option id="roll" value="roll" key="9">Roll</option>
            <option id="shift" value="shift" key="10">Shift</option>
		    <option id="slide" value="slide" key="11">Slide</option>
		    <option id="rise" value="rise" key="12">Rise</option>
		    <option id="S" value="S" key="13">S (config. volume)</option>
    	  	<option id="persistence" value="persistence" key="persistence">Persistence Length</option>
    	  </Form.Select><Button onClick={this.makePlot}>Create Plot!</Button></> : null}
    	  {this.state.layout ? <Plot layout={this.state.layout} data={this.state.data} /> : null} <Button onClick={this.clearData}>Clear Data and Plots</Button>
    	  <br/>
    	  </div>
          <br/>{this.state.loaded}<br/><br/>
          <Form.Select name="selected" value={this.state.selected} onChange={this.inputChanged}>
            <option id="1" value={1} key="1">TX3 X-ray protein-DNA analysis</option>
            <option id="2" value={2} key="2">cgDNA+ 136</option>
            <option id="3" value={3} key="3">MD AMBER parmbsc1</option>
            <option id="4" value={4} key="4">X-ray/NMR/CryoEM protein-DNA</option>
            {/*<option id="5" value={5} key="5">X-ray/NMR/CryoEM EV^2 protein-DNA</option>*/}
          </Form.Select>
          {!this.state.engaged && <Button onClick={this.loadSet}>Load Dataset</Button>}
          {this.state.engaged && <Form.Select name="tetramer" value={this.state.tetramer} onChange={this.inputChanged}>
            {this.state.items.length > 0 ? this.state.items.map((name) => (
              <option key={name} value={name}>{name}</option>
            )) : null}
          </Form.Select>}
        {/*<pre><div dangerouslySetInnerHTML={{__html: this.state.resultField}}></div></pre>
          <pre><div dangerouslySetInnerHTML={{__html: this.state.resultField2}}></div></pre>*/}
          {this.state.engaged && <SearchBar autoFocus count={this.clearSearch} shouldRenderClearButton={true} shouldRenderSearchButton={false} placeholder="type tetramer" suggestions={this.state.suggestions} onClear={this.handleClear} onSearch={this.handleSearch} onChange={this.handleChange} styles={styles} />}

        {/*<input type="text" value={this.state.tetramerText} onChange={this.inputChanged} name="tetramerText" size="6"/>*/}
          {this.state.engaged ? <DThree mean={this.state.mean[this.state.tetramer]} cov={this.state.covariance[this.state.tetramer]} tetramer={this.state.tetramer}/> : null}
          {this.state.engaged ? <><br/>Sequence: <input type="text" name="sequence" value={this.state.sequence} onChange={this.inputChanged}/>
    	  <Button onClick={this.meanTwist}>Calculate Twist of DNA chain</Button>
    	  <br/><Button disabled={this.state.processing} onClick={this.persistenceAll}>All Persistence Lengths</Button>
    	  <div id="resultPersistence"></div></> : null}
          </div>
        )
  }
}

let persistenceData = [[["ATAT",416.7287745403619],["TTAA",320.601004387896],["CTAG",332.1513061028878],["GTAC",244.47104511431863],["ATAC",551.080734575081],["GTAG",233.68351211585534],["TTAG",489.27197297442797],["CTAT",387.8126577206288],["TTAC",395.72987047177315],["TTAT",451.7536107818927],["GCGC",449.8893730326439],["CCGG",584.4018552032469],["TCGA",244.6704452059255],["ACGT",377.65058548676376],["CCGC",436.18061417448297],["CCGT",598.3960955645102],["ACGC",494.28567753029546],["TCGT",579.8289486281102],["TCGC",665.4382040298821],["TCGG",682.6275081904797],["GCAA",551.4750191015536],["ACAA",528.7515201551973],["GCAG",402.09084770831765],["CCAA",187.73455896068737],["ACAG",495.46535796442663],["TCAA",350.74258424361744],["GCAC",231.38057548929666],["ACAT",526.8421919677572],["CCAG",469.63780007478266],["GCAT",1002.8613053760585],["CCAC",569.364278391309],["ACAC",565.0820105473265],["TCAG",799.3374031952591],["CCAT",622.2158757666708],["TCAC",404.95276981953555],["TCAT",679.2092362414846],["GAAA",810.8107575430605],["AAAA",765.5705697678369],["CAAA",625.5807146267397],["GAAG",765.2473848079637],["GAAC",890.27206507643],["TAAA",509.299103127015],["AAAG",744.0426430231677],["AAAC",661.5587685673415],["AAAT",1036.9908362580325],["TAAC",693.3433297426034],["GAAT",682.339027775706],["TAAG",533.2617730646128],["CAAG",739.4700348740901],["CAAC",586.1480305486532],["TAAT",450.5505159802673],["CAAT",959.3586799380713],["AAGA",973.3649817456575],["GAGA",644.8106046306476],["CAGA",986.5786493223319],["GAGG",895.1665408453484],["AAGG",661.6496972143235],["TAGA",493.8361966940738],["AAGC",775.9387178031968],["CAGC",1004.9807674574662],["AAGT",730.6717522904349],["GAGC",488.4460346919492],["TAGG",700.4640101007274],["GAGT",720.7963472170151],["CAGG",736.2475123662331],["TAGT",646.8324276022374],["TAGC",501.4204675702566],["CAGT",911.2844152693207],["GGAA",559.9444948385385],["AGAA",870.0060419889872],["AGAT",539.8568984946007],["TGAA",591.9182482186252],["AGAC",669.5756926699519],["AGAG",686.0580287193301],["GGAC",1190.3619098803845],["CGAA",963.3624603161996],["GGAG",551.3394434575715],["CGAG",524.8928324581894],["GGAT",764.0948346742987],["TGAC",835.8442982656694],["TGAT",938.5536629376852],["TGAG",568.7703941926347],["CGAT",610.5847936419043],["CGAC",1045.0607204057767],["AGGA",1006.2722721461321],["AGGG",634.3662438751642],["GGGA",583.2215110356864],["AGGT",617.9815537836457],["AGGC",603.8734715385459],["CGGG",487.37277004252684],["CGGA",672.5088731945871],["GGGT",599.4960894091419],["GGGG",406.9001044964782],["GGGC",629.1879083783908],["TGGA",631.296762339112],["TGGG",799.3069389158634],["CGGT",438.5492364076849],["TGGC",810.6509163007644],["CGGC",580.5334829495255],["TGGT",660.5234766642524],["AGCT",620.8216542251791],["TGCA",532.4136360020335],["CGCG",652.056279307327],["GGCC",767.7017406175797],["CGCC",606.1403228657902],["TGCG",1558.589907004925],["TGCC",622.560121137302],["GGCT",987.0001109442076],["CGCT",980.2749971897907],["TGCT",550.9330254694288],["CATG",1085.2449122149887],["TATA",615.2290729971702],["GATC",1111.9269854158993],["AATT",1246.7630237457142],["CATA",898.2103634073467],["GATT",1051.2701378918407],["CATC",1460.2956047608293],["TATC",882.2622010620862],["TATT",1038.116472371246],["CATT",1028.6786240167487],["AACC",1076.3246242284804],["GACA",938.0470275610215],["AACG",732.8892113327843],["GACG",1317.5380827832162],["AACA",788.6331849739886],["CACG",943.6919288833128],["CACA",803.8383153176028],["TACA",690.0256471698434],["GACC",1102.1621535982817],["TACC",662.0232658743964],["CACC",796.9406197249397],["TACG",583.7541137539739],["GACT",1073.8842613408492],["AACT",821.9882102612412],["TACT",784.7447543592313],["CACT",1125.511768861923]],
[["ATAT",319.01807502103156],["TTAA",258.5958865482089],["CTAG",262.42626864447476],["GTAC",278.2300340691773],["ATAC",296.5388050936105],["GTAG",273.16687191800986],["TTAG",261.9582246895267],["CTAT",291.4353404972914],["TTAC",270.1113653823751],["TTAT",292.6063162410017],["GCGC",337.76347004028594],["CCGG",320.0705654163649],["TCGA",304.6442727464671],["ACGT",374.4615674498627],["CCGC",331.56447709451464],["CCGT",349.9757233597836],["ACGC",355.7356354990418],["TCGT",347.3141201343603],["TCGC",328.9674320157178],["TCGG",313.2023869625751],["GCAA",293.4460511528658],["ACAA",319.52532548915036],["GCAG",305.74857286641225],["CCAA",296.85943492898565],["ACAG",324.8818188630538],["TCAA",280.91699877355],["GCAC",312.4343496740689],["ACAT",355.48473230943137],["CCAG",308.61692840131275],["GCAT",338.7566271396223],["CCAC",311.55204988600036],["ACAC",330.16367514693934],["TCAG",294.00744113026576],["CCAT",339.8307626933394],["TCAC",302.1505640949217],["TCAT",332.2756803439738],["GAAA",575.5619657907635],["AAAA",585.8984220780127],["CAAA",565.7303024577756],["GAAG",568.9512747510088],["GAAC",574.926680766056],["TAAA",560.995611619632],["AAAG",581.6967679094639],["AAAC",594.7334399283153],["AAAT",618.5766397633167],["TAAC",566.3585383364838],["GAAT",612.9407253066057],["TAAG",561.2965805995406],["CAAG",565.5775911041126],["CAAC",583.7568133483078],["TAAT",607.3497837549295],["CAAT",610.389341394848],["AAGA",562.4307134134765],["GAGA",555.0923050949444],["CAGA",545.2448230564959],["GAGG",555.1700194823644],["AAGG",564.2429505874886],["TAGA",528.1614676906479],["AAGC",581.0595439063605],["CAGC",564.0126946289845],["AAGT",594.5532728694432],["GAGC",565.3906980544634],["TAGG",534.3185613650663],["GAGT",586.0966666680689],["CAGG",552.1336803854049],["TAGT",566.6598867370199],["TAGC",551.9678272938861],["CAGT",579.7815181060735],["GGAA",543.215339298579],["AGAA",582.3256489456774],["AGAT",600.649487274424],["TGAA",515.5236416459887],["AGAC",576.7111414979742],["AGAG",583.9952559996626],["GGAC",539.8330017978265],["CGAA",511.5475277148904],["GGAG",545.6756783014397],["CGAG",504.84946380713956],["GGAT",563.3608570472974],["TGAC",517.5090086572372],["TGAT",528.5255228027331],["TGAG",508.0871937962104],["CGAT",524.4596975327979],["CGAC",500.5419524586712],["AGGA",531.2128256153396],["AGGG",545.0012818599229],["GGGA",510.7937617173989],["AGGT",537.4418382806194],["AGGC",553.5657886989119],["CGGG",495.21138177733894],["CGGA",489.604947793191],["GGGT",525.7216685660823],["GGGG",515.5676144693069],["GGGC",528.6991687007843],["TGGA",499.1156458659142],["TGGG",503.27249592275297],["CGGT",496.74927812391553],["TGGC",513.5126894795753],["CGGC",499.1135616143504],["TGGT",510.42860317245186],["AGCT",681.7790475263874],["TGCA",601.1591785589178],["CGCG",594.74676165827],["GGCC",618.3327879520142],["CGCC",611.8821334612957],["TGCG",602.4985366769291],["TGCC",609.5420524594351],["GGCT",652.7535232461132],["CGCT",639.8040219967252],["TGCT",645.9566460222836],["CATG",720.3417157347606],["TATA",678.9145718647762],["GATC",757.0680226670268],["AATT",795.7926747341292],["CATA",710.2604332665387],["GATT",775.2391094053665],["CATC",731.4972005486474],["TATC",716.7993543056108],["TATT",755.0576115096196],["CATT",753.2347481972793],["AACC",656.3079742020261],["GACA",591.6320264755162],["AACG",615.1288440808344],["GACG",594.4011957056633],["AACA",620.4878556790381],["CACG",584.3237004833879],["CACA",588.9512564250203],["TACA",569.4250915467006],["GACC",647.031156422624],["TACC",612.3307218206193],["CACC",607.4386381232266],["TACG",589.0592865560096],["GACT",661.6493577173081],["AACT",697.6485659677437],["TACT",649.198065721569],["CACT",658.6604642528158]],
[["ATAT",285.74660478265065],["TTAA",228.7286414145127],["CTAG",244.48384641097059],["GTAC",262.2361148036629],["ATAC",262.28047003433653],["GTAG",250.19208129033717],["TTAG",224.88435203457132],["CTAT",253.91500530567023],["TTAC",231.02622532765886],["TTAT",248.267255381557],["GCGC",312.75992718209096],["CCGG",291.96156466053856],["TCGA",273.04348490344086],["ACGT",283.34550034850406],["CCGC",295.9501675953176],["CCGT",294.0977091257742],["ACGC",296.3668540509002],["TCGT",279.0566503869202],["TCGC",307.5249534910934],["TCGG",305.0028899261057],["GCAA",286.679707399371],["ACAA",261.8646795482499],["GCAG",265.9869308292012],["CCAA",252.74489794400733],["ACAG",273.59745532579944],["TCAA",248.31983891621394],["GCAC",279.25813486437517],["ACAT",318.782678980545],["CCAG",277.8466143200518],["GCAT",273.3272643213285],["CCAC",278.6642909576198],["ACAC",276.8335696127179],["TCAG",287.94176989923017],["CCAT",292.3185858910762],["TCAC",286.0864815505253],["TCAT",271.01537581364653],["GAAA",521.8529064699128],["AAAA",544.8658198674274],["CAAA",541.8140557920942],["GAAG",571.0689909016535],["GAAC",576.8588618596398],["TAAA",505.5924283638281],["AAAG",505.99542216554727],["AAAC",543.5603647625461],["AAAT",625.7021977549485],["TAAC",499.95727830399954],["GAAT",601.9992112758417],["TAAG",483.2187125356188],["CAAG",505.52202183369263],["CAAC",529.1226073025804],["TAAT",520.0112122582742],["CAAT",563.7960416762467],["AAGA",478.5074643666814],["GAGA",467.45654147445805],["CAGA",470.129925315031],["GAGG",485.5901285830671],["AAGG",491.60183292180506],["TAGA",432.0137187244915],["AAGC",530.1178433173247],["CAGC",566.8858683949205],["AAGT",563.2266441901625],["GAGC",542.1614111086941],["TAGG",513.320673222814],["GAGT",524.7830348375738],["CAGG",485.4144585269071],["TAGT",517.0186314193295],["TAGC",477.33875343826963],["CAGT",533.9771487511213],["GGAA",485.3035948237757],["AGAA",452.2103780307637],["AGAT",504.36709951582515],["TGAA",473.8816285691705],["AGAC",481.5515379424191],["AGAG",461.928419611233],["GGAC",494.63890401873726],["CGAA",487.0343365385525],["GGAG",467.6667307032311],["CGAG",482.93925013119775],["GGAT",481.0443455099914],["TGAC",517.3077407237997],["TGAT",525.9170886993361],["TGAG",485.0822929858733],["CGAT",473.009153293215],["CGAC",497.3145826009677],["AGGA",412.59462612757847],["AGGG",528.5834381578433],["GGGA",432.37454766445626],["AGGT",519.061904030197],["AGGC",513.8845544187095],["CGGG",444.78435080357553],["CGGA",396.069550363611],["GGGT",534.8305577156473],["GGGG",449.3302504455544],["GGGC",492.7122185085492],["TGGA",425.3677129099936],["TGGG",417.22022532452525],["CGGT",410.65065946796346],["TGGC",440.3256103014547],["CGGC",457.4583126493936],["TGGT",456.877742986078],["AGCT",765.8822255812353],["TGCA",544.5401332634416],["CGCG",524.8070048711364],["GGCC",653.0053146583001],["CGCC",586.0231313587641],["TGCG",522.9924664465469],["TGCC",552.8785044951691],["GGCT",695.5805459379733],["CGCT",551.4011521607433],["TGCT",577.2511889901037],["CATG",616.2677414218713],["TATA",540.2694700205044],["GATC",681.451507274799],["AATT",882.336260506322],["CATA",559.6569551438288],["GATT",778.9457686019559],["CATC",631.5792963764388],["TATC",696.9041429985645],["TATT",680.5483768439095],["CATT",751.6026171455189],["AACC",632.1354177629571],["GACA",605.6072105689013],["AACG",557.968906680198],["GACG",575.6013476975946],["AACA",567.4434499137556],["CACG",485.3884007892327],["CACA",503.9012679807989],["TACA",493.76116057132487],["GACC",631.8244357846878],["TACC",547.428861523219],["CACC",535.4332911774025],["TACG",512.8344744744508],["GACT",650.3504995010642],["AACT",674.5499020830605],["TACT",535.706794928708],["CACT",550.7181038118384]],
[["ATAT",139.3558220965114],["TTAA",249.6856045540212],["CTAG",227.7691678316397],["GTAC",322.2312393682881],["ATAC",441.8162297560926],["GTAG",173.6569054919804],["TTAG",445.1158553989212],["CTAT",81.04412263186433],["TTAC",242.3916610402668],["TTAT",273.4195928110714],["GCGC",384.1822046153868],["CCGG",471.07597460732865],["TCGA",234.08876099038983],["ACGT",310.3814167795879],["CCGC",308.97000883780476],["CCGT",253.94138803633396],["ACGC",437.4733498474364],["TCGT",404.6896041975404],["TCGC",555.9083753826739],["TCGG",408.39541131108473],["GCAA",538.1761263943125],["ACAA",393.3964499991858],["GCAG",410.457658802886],["CCAA",162.223564897129],["ACAG",413.1214643786801],["TCAA",241.95221251251473],["GCAC",286.0142173371032],["ACAT",497.6780735877335],["CCAG",389.6940040252438],["GCAT",667.171075783049],["CCAC",399.6377696820286],["ACAC",370.22320152266593],["TCAG",557.370822692077],["CCAT",580.6051077065713],["TCAC",216.33737125913711],["TCAT",463.0958602879626],["GAAA",669.1195275591707],["AAAA",652.9911381144914],["CAAA",563.545045293609],["GAAG",782.877054112511],["GAAC",1009.9531275407845],["TAAA",347.4124047797656],["AAAG",480.92799860943023],["AAAC",748.1191771244712],["AAAT",755.3179671730117],["TAAC",481.01278638608164],["GAAT",598.4885967749293],["TAAG",450.72733047395093],["CAAG",526.506879041268],["CAAC",429.96523847238035],["TAAT",454.26435187277394],["CAAT",621.4928108223351],["AAGA",671.7907772118486],["GAGA",630.4066743661783],["CAGA",847.6082943570628],["GAGG",733.0976913395597],["AAGG",540.9956052623006],["TAGA",406.3783736686873],["AAGC",886.9589596663114],["CAGC",942.2222475884523],["AAGT",554.5068678394989],["GAGC",699.8768065377091],["TAGG",494.5294769039792],["GAGT",647.1432680575427],["CAGG",344.87925592902434],["TAGT",556.8483744703406],["TAGC",655.4150310239563],["CAGT",660.4138049730828],["GGAA",545.8483035220535],["AGAA",659.7083542384672],["AGAT",772.4023559744576],["TGAA",534.472297054069],["AGAC",520.6189415209344],["AGAG",639.6222281132743],["GGAC",1062.1587865148808],["CGAA",895.5065829543599],["GGAG",642.0145610464724],["CGAG",470.61720247145917],["GGAT",605.5127040245048],["TGAC",621.9801730163599],["TGAT",740.3092840529732],["TGAG",293.9907014125445],["CGAT",648.1275082438137],["CGAC",663.3679268725781],["AGGA",788.2485293718836],["AGGG",428.5762256728407],["GGGA",353.63617362062837],["AGGT",582.8677452130856],["AGGC",464.7591048416629],["CGGG",488.3200191386441],["CGGA",526.854963285258],["GGGT",469.17567330302046],["GGGG",382.22110471527105],["GGGC",587.1597098191047],["TGGA",586.3450795605385],["TGGG",654.2313941071625],["CGGT",471.82433975881713],["TGGC",616.5378595727291],["CGGC",417.68370859396254],["TGGT",579.2332296182782],["AGCT",832.4808582671873],["TGCA",508.62187203425844],["CGCG",746.1698265365601],["GGCC",952.0269062110722],["CGCC",656.1768773759802],["TGCG",1197.2036052733665],["TGCC",633.8458780996334],["GGCT",805.763214485882],["CGCT",819.2630499622653],["TGCT",491.98286649093797],["CATG",815.0530625845583],["TATA",281.5490239187541],["GATC",946.0409584238058],["AATT",889.5524065353035],["CATA",925.2206387211438],["GATT",1020.6778264048671],["CATC",1279.961186859604],["TATC",609.6683135807915],["TATT",675.3388448597151],["CATT",732.6844774149299],["AACC",705.6809599840628],["GACA",656.8869852898539],["AACG",742.277843590248],["GACG",922.8960334737271],["AACA",688.8883908228605],["CACG",970.0619876084336],["CACA",644.81978828756],["TACA",651.2273080125776],["GACC",1052.5424050425083],["TACC",462.35290310284813],["CACC",915.505316558305],["TACG",483.19487872873145],["GACT",927.4703169508509],["AACT",722.6183735435474],["TACT",665.5075148668983],["CACT",586.3448202101646]]];

export default TetramerReader;

