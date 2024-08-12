![Logo LDD Suite](./images/bistoniheader.png)

<h1 align='center'>
London Dispersion Density (LDD) Suite
</h1>

[![Github](https://img.shields.io/badge/GitHub-bistonigroup-181717.svg?style=social&logo=github)](https://github.com/bistonigroup)
[![X (formerly Twitter) Follow](https://img.shields.io/twitter/follow/easyaspython?label=Follow&style=social)](https://x.com/BistoniGroup)
[![Website](https://img.shields.io/badge/Giovanni%20Bistoni%20Group-blue?style=social&logo=wordpress&labelColor=grey&label=WordPress)](https://giovannibistoni.wordpress.com/)
[![Python 3.x](https://img.shields.io/badge/3.x-blue?style=flat-square&logo=Python&logoColor=blue&label=Python&labelColor=grey)](https://www.python.org/download/releases/3.0/)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg?style=flat-square)](https://www.gnu.org/licenses/lgpl-3.0)

The LDD suite is designed to compute atomic contributions to the London dispersion energy, the London dispersion density function, and the London dispersion density difference function. It includes the following tools: _lddensityd4.py_, _lddensityd3.py_, and _lddensitydifference.py_, each utilizing different dispersion corrections or computational approaches. The suite includes a Visual Molecular Dynamics (VMD) script, _PubQualityVMD.tcl_, which is used to generate images suitable for scientific publications and to visually analyze atomic contributions.
## Table of contents
- [ Tools Description ](#tools-description)
- [ Prerequisites ](#prerequisites)
- [ Usage ](#usage)
- [ Arguments Description ](#arguments-description)
- [ Example Images ](#example-images)
- [ Credits ](#credits)
- [ License ](#license)
- [ Contact ](#contact)

## Tools Description

### _lddensityd3.py_ and _lddensityd4.py_
These Python scripts are designed to compute atomic contributions to the London dispersion energy and the London dispersion density function, relying on D3 and D4 corrections, respectively.

- **Input**: 
  - `{basename}.xyz` : a classic _xyz_ file containing atomic coordinates of the system in ångström. 
- **Output**:
  - `.d3atomwise.txt` or `.d4atomwise.txt` : a file containing atomic coordinates of the system with an additional column indicating each atom contribution to the London dispersion energy.
  - `.d3out.txt` or `.d4out.txt` : the original output generated by Grimme's DFT-D3 or DFT-D4 executable, respectively.
  - `.d3omega.cube` or `.d4omega.cube` : a _cube_ file that stores volumetric data of the LDD function, which can be used to easily visualize atomic contributions to London dispersion.

### _lddensitydifference.py_
This script computes the LDD difference function, which can be used to easily visualize the difference between atomic LD contributions for two different molecular structures (e.g., a pair of structures along a reaction profile).

- **Input**:
  - `{basename}.atomwise.txt` : a file containing atomic coordinates in the _xyz_ format and in the 5th column the differences of the atomic London dispersion energies.
- **Output**:
  - `{basename}.omega.cube`: the LDD difference function in _cube_ format.

### _PubQualityVMD.tcl_
This script is designed mainly to visualize LD density (difference) function from a _cube_ file, but it also allows you to load an _xyz_ file to view a molecule structure. The script adjusts VMD settings to generate publication-quality images.
Additionally, it includes an improved pick event feature that displays atomic contributions to dispersion energy in the terminal, extracted from the _*atomwise.txt_ file.

## Prerequisites
- For _lddensityd4.py_ and _lddensityd3.py_, ensure that _xyz_ file and the respective DFT-D3 or DFT-D4 executable are in the current path. The executable can be downloaded from the [official website](https://www.chemie.uni-bonn.de/grimme/de/software).
- For _lddensitydifference.py_, ensure `{basename}.atomwise.txt` is in the current path.
- To use _PubQualityVMD.tcl_ you must have VMD installed on your computer. VMD can be downloaded from the [official website](https://www.ks.uiuc.edu/Research/vmd/).
- Python 3.x with standard libraries.

## Usage

### _lddensityd3.py_ and _lddensityd4.py_

1. Ensure `{basename}.xyz` file is at the current path.
2. Run the script:
- For _lddensityd3.py_:

    <pre><code style="font-size: 13px;">python3 lddensityd3.py &lt;basename&gt; [--npoints NP] [--func FUNC] [--damp DAMP] [--nprocs NPROCS]</code></pre>
    
- For _lddensityd4.py_:

  <pre><code style="font-size: 13px;">python3 lddensityd4.py &lt;basename&gt; [--npoints NP] [--func FUNC] [--charge CHARGE] [--s9 S9] [--nprocs NPROCS]</code></pre>

### _lddensitydifference.py_

1. Ensure `{basename}.atomwise.txt` file is at the current path.
2. Run the script:

    <pre><code style="font-size: 13px;">python3 lddensitydifference.py &lt;basename&gt; [--npoints NP] [--nprocs NPROCS]</code></pre>

### _PubQualityVMD.tcl_
1. Ensure the input file (_cube_ or _xyz_) you wish to visualize is correctly named and located within an accessible directory (it's recommended to work within a single directory).
2. Modify the script to include the path name in `set file_name` command.
3. Run the script with the following command:

    <pre><code style="font-size: 13px;">vmd -e PubQualityVMD.tcl</code></pre>

When a _cube_ file is loaaded, if the _*atomwise.txt_ file is present in the same folder, you can use the _pick_ feature in VMD (_Mouse_ > _Pick_) to visualize the dispersion energy of individual atoms upon selection.


> [!NOTE] 
> Some functions within the script are commented out by default for safety. However, if you are interested in exploring or utilizing these functions, you can easily uncomment them. \
> Auto color scaling is enabled by default. \
> Rendering is disabled by default. Set _render_mode_ to 1 to enable it ([ImageMagick](https://imagemagick.org/script/download.php) is required for converting to _png_)

> [!TIP]
> Change the rotation angles (_rot_x_, _rot_y_, _rot_z_) and the moleculare size (_scale_val_) to display the molecule according to your preference.
 
## Arguments Description

Below is a detailed table of the arguments that can be used with our script. Each entry provides you with the name of the argument, a brief description, information on whether the argument is optional, and the default value it takes if not specified by the user. This table is designed to help you quickly understand how to configure the script to meet your specific requirements.

| Argument  | Description | Optional | Default Value |
|-----------|-------------|----------|---------------|
| ***basename*** | The base name for the input _.xyz_ file. | No | *None* |
| ***npoints*** | Specifies the number of grid points for each dimension in the density calculation. | Yes | 80 |
| ***func***    | Defines the functional to be used in the calculation. Refer to [DFT-D3](https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/man.pdf)/[DFT-D4](https://github.com/dftd4/dftd4/tree/main) manual for all available functionals. | Yes | `b3-lyp` |
| ***damp***    | Specifies the damping function to be used in the D3 calculation, refer to [DT-D3 manual](https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/man.pdf) for all available damping functions. | Yes | `bj` |
| ***charge***  | Sets the overall charge of the molecule being analyzed in the D4 calculation. | Yes | 0 |
| ***s9***      | A coefficient that scale the ATM (Axilrod-Teller-Muto) term for three-body dispersion in the D4 calculation. Allows users to adjust the contribution of three-body interactions in the dispersion energy. Set to 0 to eliminate this contribution. | Yes | 0 |
| ***nprocs***  | Determines the number of processors used for parallel computation. | Yes | 1 |

## Example Images
The images displayed below are different visual representations of a Benzene-Lithium complex generated using the _PubQualityVMD.tcl_ script. On the left, you can see the London dispersion density (_cube_ file), while on the right the simple molecular structure (_xyz_ file). Regarding LDD, red indicates high dispersion contribution, while blue indicates low contribution.

<p align="center">
  <img src="./images/BenzLi_40_b3-lyp_0_abc1.0_d4omega.cube.png" alt="London dispersion density for Benzene-Li" width="400" style="margin-right: 100px;"/>
  <img src="./images/BenzLi.xyz.png" alt="Water" width="400"/> 
</p>

## Credits
If you use the LDD suite in your research or any publication, please cite the author(s).

## License
Distributed under GNU Lesser General Public License. See _LICENSE_ for more information.

## Contact
For general inquiries, feedback, or assistance with using the LDD suite, please contact us at:

Gianluca Regni - [_gianluca.regni@studenti.unipg.it_](mailto:gianluca.regni@studenti.unipg.it)  
Lorenzo Baldinelli - [_lorenzo.baldinelli@studenti.unipg.it_](mailto:lorenzo.baldinelli@studenti.unipg.it)  
Giovanni Bistoni - [_giovanni.bistoni@unipg.it_](mailto:giovanni.bistoni@unipg.it) 
