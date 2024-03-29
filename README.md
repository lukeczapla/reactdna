
# reactdna - Movie animation visualizer for DNA and Monte Carlo simulator
Animates modes of motion of DNA bases, phosphates and base-pair step regions (136 x 30, or 4080 movies per model), from data extracted from the RCSB Protein Databank and from 10 microsecond simulations of a library of DNA sequences with AMBER parmbsc1.  Reports simulated properties of central junctions for three models (X-ray from a survey of protein-DNA structures, cgDNA+ 136 parameterized from AMBER parmbsc1, and the direct data from AMBER parmbsc1 MD simulations).

Using WebGL, it animates the motions of DNA as expressed by eigenvectors of the inverse covariance matrix for a subset of 30 coordinates representing the state of the central junction of a DNA tetramer duplex (e.g., ATAG representing 5'-ATAG-3' paired with 5'-CTAT-3'), for a set of 136 non-redundant/unique tetramers (out of 256), and displays the means, eigenvectors (normalized), and eigenvalues (effective rigidity) of the inverse covariance matrices, and it plots variations over the sets.  Computes persistence length, statistical properties of simulated central junctions of tetramers (using the inverse covariance matrix and Gaussian random number generators) and also allows perturbing the 3D structure along eigenvectors (outputting PDB files compatible with pymol, UCSF Chimera, VMD and other 3D visualization software programs), and reports parameters in Curves+/cgDNA+ and 3DNA coordinate systems.

Requires React, three, numeric, plotly.js, and react-plotly.js installed with npm for node, and the input files (30x30 covariance matrices and 1x30 row vector of mean values for dimeric or tetrameric sequences), and the browser should be WebGL enabled for the 3D rendering of movies live in the browser.

The Monte Carlo simulation worker thread (and the thread for running all persistence lengths for the 136 tetramers) requires two packages installed:
```bash
npm install react-app-rewired worker-loader --save-dev
```

The build folder contains the entire package (with the start page index.html) for deployment on a web server or locally with "Load File", there is no compilation required.

![Screenshot](https://github.com/lukeczapla/reactdna/blob/master/Website.png?raw=true)


