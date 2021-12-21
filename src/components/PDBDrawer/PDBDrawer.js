import React, {useState, useEffect, createRef} from 'react';
import Button from 'react-bootstrap/Button';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import * as $3Dmol from '3dmol';
import Form from 'react-bootstrap/Form';
import * as ref from '../References/References.jsx';

const PDBDrawer = (props) => {

    const [viewer, setViewer] = useState(null);
    const [graphItem, setGraphItem] = useState(0);
    const [active, setActive] = useState(false);
    const [coordinates, setCoordinates] = useState([]);
    const [sequences, setSequences] = useState([]);
    const [proteins, setProteins] = useState([]);
    const [PDBText, setPDBText] = useState("");
    const [menuItems, setMenuItems] = useState([]);
    const [includeProtein, setIncludeProtein] = useState(false);
    const child1 = createRef();

    useEffect(() => {
        let sources = [["1KX5.CURVE.coord", "ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGA", "1KX5 - X-Ray Structure of the Nucleosome Core Particle, NCP147, at 1.9 A Resolution", "1kx5p.pdb"],
            ["1TGH.CURVE.coord", "CGTATATATACG", "TATA Binding Protein (TBP)/DNA Complex", "1tghp.pdb"],
            ["6ON0.CURVE.coord", "TTTATAGCTAGCTATAA", "Structure of N15 Cro Complexed with Consensus Operator DNA", "6on0p.pdb"],
            ["5WV7.CURVE.coord", "CCGGGGTACCCCGG", "Crystal Structure of d(CCGGGGTACCCCGG)2 at 1.4 Å resolution"],
            ["1TC3.CURVE.coord", "AGTTCTATAGGACCCCCC", "Transposase TC3A1-65 from C. elegans", "1tc3p.pdb"],
            ["1JGR.CURVE.coord", "CGCGAATTCGCG", "Crystal Structure Analysis of the B-DNA Dodecamer CGCGAATTCGCG with Thallium Ions"]];
        let items = [];
        items.length = sources.length;
        let values = [];
        let coord = [];
        let seq = [];
        let protein = [];
        protein.length = sources.length;
        coord.length = sources.length;
        seq.length = sources.length;
        sources.forEach((set, index) => {
            items[index] = set[2];
            seq[index] = set[1];
            let values = [];
            fetch(set[0]).then((r) => r.text()).then(text => {
                //console.log(text);
                let lines = text.split("\n");
                values.length = lines.length;
                let i = 0;
                lines.forEach(line => {
                    values[i++] = parseFloat(line);
                });
                coord[index] = values;
            });
            if (set[3] !== undefined) fetch(set[3]).then((r) => r.text()).then(text => {
                console.log(text);
                protein[index] = [text];
            });
        });
        //console.log(coord, items, seq);
        setMenuItems(items);
        setCoordinates(coord);
        setSequences(seq);
        setProteins(protein);
    }, []);

    const inputChanged = (e) => {
        if (e.target.name === "includeProtein") {
            setIncludeProtein(e.target.checked);
        }
        if (e.target.name === "graphItem") {
            setGraphItem(e.target.value);
        }
        if (e.target.name === "graphItem" || e.target.name === "loadButton") {
            setActive(true);
            let val = e.target.value;
            if (e.target.name === "loadButton") val = graphItem;
            //console.log(coordinates[e.target.value]);
            let data = ref.writePDBFull(coordinates[val], sequences[val], sequences[val].length);
            if (proteins[val] === undefined || !includeProtein) setPDBText(data);
            if (proteins[val] !== undefined) data += proteins[val][0];
            if (includeProtein) setPDBText(data);

            let v;
            let element = child1.current;
            let config = {backgroundColor: "white"};
            if (viewer === null) {
                v = $3Dmol.createViewer(element, config);
                setViewer(v);
            } else {
                v = viewer;
                v.clear();
            }
            v.addModel(data, "pdb");
            //v.setStyle({}, {stick: {colorscheme: 'default'}});
            v.setStyle({}, {cartoon: {color: '#FF00FF'}});
            v.setStyle({chain: 'X'}, {stick: {colorscheme: 'default'}});
            v.setStyle({chain: 'Y'}, {stick: {colorscheme: 'default'}});
            v.zoomTo();
            v.render();
            v.zoom(1.2, 1000);
        }
    }

    return (
      <>
          <Container>
          {menuItems.length > 0 && <Form.Select name="graphItem" value={graphItem} onChange={inputChanged}>
              {menuItems.map((name, i) => <option id={i} value={i} key={i}>{name}</option>)}
          </Form.Select>}
              <Row><Button name="loadButton" onClick={inputChanged}>Load/Build Structure</Button><Form.Check inline type="checkbox" name="includeProtein" label="Include protein in PDB output" onChange={inputChanged} checked={includeProtein}/></Row>
          <center><div id="molbox" style={{width: "1000px", height: "800px", position: "relative"}} ref={child1}></div></center>
          {active && <pre><div dangerouslySetInnerHTML={{__html: PDBText}}></div></pre>}
          </Container>
      </>
    );
}


export default PDBDrawer;