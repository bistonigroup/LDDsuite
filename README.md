![Header](./images/header.png)

<h1 align='center'>
London Dispersion Density (LDD) Suite
</h1>

[![X (formerly Twitter) Follow](https://img.shields.io/twitter/follow/easyaspython?label=Follow&style=social)](https://x.com/BistoniGroup)
[![Website](https://img.shields.io/badge/Bistoni%20Group-Official%20Website-blue?style=flat-square&logo=world&logoColor=white)](https://bistonigroup.org/)
[![Python 3.x](https://img.shields.io/badge/3.x-blue?style=flat-square&logo=Python&logoColor=blue&label=Python&labelColor=grey)](https://www.python.org/download/releases/3.0/)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg?style=flat-square)](https://www.gnu.org/licenses/lgpl-3.0)

The **LDD suite**  computes atomic contributions to the London dispersion energy and the dispersion density. It is based on the [**Atomic Decomposition of London Dispersion energy (ADLD)**](https://pubs.acs.org/doi/full/10.1021/acscentsci.5c00356) method and includes the following tools:

- `lddensity.py` – computes atomic dispersion contributions using the selected Grimme's dispersion correction and generates a `.cube` file with the London dispersion density function.
- `vmd.tcl` – VMD script to automatically set visualization parameters for publication-quality images.

# Table of contents
- [ Tools Description ](#tools-description)
- [ Installation ](#installation)
- [ Basic Usage ](#basic-usage)
- [ Requirements ](#requirements)
- [ Arguments Description ](#arguments-description)
- [ Generate Only the Cube File](#generate-only-the-cube-file)
- [ Dispersion Density Difference ](#dispersion-density-difference)
- [ Renderings ](#renderings)
- [ Credits ](#credits)
- [ License ](#license)
- [ Contacts ](#contacts)

# Tools Description

## lddensity.py
Python script to compute atomic contributions to the London dispersion energy and the London dispersion density function.

- **Input**: 
  - `{basename}.xyz`: a classic XYZ file containing atomic coordinates of the system in ångström. 
- **Output**:
  - `{basename}.atomwise.txt`: a file containing atomic coordinates of the system with an additional column indicating each atom contribution to the London dispersion energy.
  - `{basename}.omega.cube`: a `.cube` file that stores volumetric data of the LDD function, which can be used to easily visualize atomic contributions to London dispersion.

> [!NOTE]
> If atomic contributions are computed or if you want to plot the dispersion density difference between two molecular systems, you can use the script solely to generate the `.cube` file.  
> In the latter case, first compute the atomic contributions for both systems, then build a new `{basename}.atomwise.txt` file containing the difference of the atomic contributions between the two structures.  
> See [Generate Only the Cube File](#generate-only-the-cube-file) for more information.

## vmd.tcl
Script designed to visualize LD density function from a `.cube` file, adjusting VMD settings to generate publication-quality images.

# Installation

Navigate to the directory you want to run the code and clone the repository:
```bash
cd /path/to/target/directory
git clone git@github.com:bistonigroup/LDDsuite.git .
```

> [!IMPORTANT]  
> To ensure the code runs correctly and to avoid conflicts with your global Python packages, it is highly recommended to use a **virtual environment**.

Create the environment:

```bash
python -m venv venv
```

Load the environment:
```bash
source venv/bin/activate
```
> [!NOTE]
> Now you should see `(venv)` at the beginning of your terminal prompt.

Ensure `pip` is up to date and install the required packages:
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

Once you are done, you can exit the virtual environment by simply typing:
```bash
deactivate
```

You don't need to reinstall python packages every time you want to run the code. Before run the code be sure to reload the environment with:
```bash
source venv/bin/activate
```

# Basic Usage

## 1. Generate Density Data
First, run the Python script to perform ADLD and generate the volumetric data (e.g., a .cube file). Use the following command:
```bash
python lddensity.py {basename}                
```

To see all available options and advanced settings, use the help flag:
```bash
python lddensity.py -h                   
```
## 2. Visualization with VMD
Once the `.cube` file is generated, you can visualize the results with publication quality using the provided VMD script.  
First, open `vmd.tcl` in a text editor and set the `file_name` variable to match your generated `.cube` file

```tcl
set file_name "your_output_file.cube"
```
Finally, run the following command to load `.cube` file and apply the visualization settings:
```bash
vmd -e vmd.tcl
```

> [!NOTE] 
> When a `.cube` file is loaded, if the `{basename}.atomwise.txt` file is present in the same folder, you can use the _pick_ feature in VMD (_Mouse_ > _Pick_) to visualize the dispersion energy of individual atoms upon selection.
> Auto color scaling is enabled by default.

> [!TIP]
> For better visualization, set `autoscale` to `1` and choose symmetric values for `colorscale_min` and `colorscale_max`.

# Requirements
- Refer to `requirements.txt` for Python script.
- VMD must be installed for running `vmd.tcl` ([download](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)).
 
# Arguments Description

Below is a detailed table of the arguments that can be used with our script. Each entry provides you with the name of the argument, a brief description, information on whether the argument is optional, and the default value it takes if not specified by the user. This table is designed to help you quickly understand how to configure the script to meet your specific requirements.

| Argument  | Description | Optional | Default Value |
|-----------|-------------|----------|---------------|
| **basename** | The base name for the input `.xyz` file. | No | *None* |
| **d** | Specifies the D correction (D3 or D4) | Yes | 4 |
| **npoints** | Specifies the number of grid points for each dimension in the density calculation. | Yes | 80 |
| **func**   | Defines the functional to be used in the calculation. Refer to [DFT-D3](https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/man.pdf)/[DFT-D4](https://github.com/dftd4/dftd4/tree/main) manual for all available functionals. | Yes | `b3-lyp` |
| **damp**    | Specifies the damping function to be used in the D3 calculation, refer to [DT-D3 manual](https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/man.pdf) for all available damping functions. | Yes | `bj` |
| **charge**  | Sets the overall charge of the molecule being analyzed in the D4 calculation. | Yes | 0 |
| **abc**     | Three-body contributions. | Yes | `False` |
| **nprocs**  | Determines the number of processors used for parallel computation. | Yes | 1 |
| **onlycube**  | Generates only the `.cube` file. See [Generate Only the Cube File](#generate-only-the-cube-file) section | Yes | False |

# Generate Only the Cube File
If you have already computed the atomic contributions (for example using another software), or if you want to plot the **dispersion density difference function**, you can use the `--onlycube` flag in `lddensity.py`.  

In this case, the input file is **NOT** a standard `.xyz`, but rather an `.xyz` file with an extra column containing the atomic dispersion contributions in kcal/mol.  
The file must be named `{basename}.atomwise.txt`, which corresponds to the file produced by `lddensity.py` when the `--onlycube` flag is **NOT** used.

An example of `.atomwise.txt` file for the $C_6H_6–Li$ system is:

```
#
atom   x(Å)         y(Å)        z(Å)        Edisp(kcal/mol)
C      -2.378672     1.005050   -0.079453   -2.128400
C      -2.318663    -0.385217   -0.085888   -2.128400
C      -1.085955    -1.028059   -0.147588   -2.128400
C       0.086589    -0.280636   -0.203051   -2.128300
C       0.026549     1.109641   -0.196866   -2.127300
C      -1.206135     1.752463   -0.134993   -2.127300
H      -3.337288     1.505059   -0.038209   -0.486950
H      -3.230736    -0.966137   -0.050005   -0.487450
H       0.937922     1.690612   -0.246784   -0.486450
H      -1.253275     2.833395   -0.136867   -0.486550
H      -1.039495    -2.108952   -0.159689   -0.487350
H       1.044500    -0.780664   -0.258085   -0.486850
Li     -1.018560     0.350445    2.545363   -3.707400
```

# Dispersion Density Difference
To generate the **dispersion density difference** between two molecular systems:  

1. Compute the atomic dispersion contributions for both systems separately.  
2. Be sure there's an atom-to-atom mapping between the two structures. 
3. Build a new `{basename}.atomwise.txt` file where the atomic dispersion contribution for each atom corresponds to the **difference** between its values in the two systems.
4. Use this difference file as input to generate the `.cube` file with `--onlycube` (refer to [Generate Only the Cube File](#generate-only-the-cube-file) section).  

This procedure allows you to visualize the difference between atomic LD contributions for two different molecular structures (e.g., a pair of structures along a reaction profile).

An example command could be:

```bash
python3 lddensity.py system_diff --onlycube --npoints 80 --nprocs 2
```

# Renderings
The image below is a rendering generated with `vmd.tcl` for a $C_6H_6–Li$ system. Red/blu color indicates high/low atomic dispersion contribution.

<p align="center">
  <img src="./images/benzli-d4.80.omega.cube.png" alt="London dispersion density for Benzene-Lithium" width="400" style="margin-right: 100px;"/>
</p>

# Credits
If you use the LDD suite in your research or any publication, please cite:
- [https://pubs.acs.org/doi/full/10.1021/acscentsci.5c00356](https://pubs.acs.org/doi/full/10.1021/acscentsci.5c00356)

- [https://pubs.acs.org/doi/abs/10.1021/acs.jctc.3c00977](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.3c00977)

- [https://github.com/bistonigroup/LDDsuite](https://github.com/bistonigroup/LDDsuite)

# License
Distributed under GNU Lesser General Public License. See _LICENSE_ for more information.

# Contacts
For general inquiries, feedback, or assistance with using the LDD suite, please contact us at:
 
Gianluca Regni - [_gianluca.regni@dottorandi.unipg.it_](mailto:gianluca.regni@dottorandi.unipg.it) \
Lorenzo Baldinelli - [_lorenzo.baldinelli@dottorandi.unipg.it_](mailto:lorenzo.baldinelli@dottorandi.unipg.it)   
Giovanni Bistoni - [_giovanni.bistoni@unipg.it_](mailto:giovanni.bistoni@unipg.it)
