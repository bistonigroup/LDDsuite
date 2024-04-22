"""
(c) 2024, Gianluca Regni, Lorenzo Baldinelli, Filippo De Angelis, Giovanni Bistoni
License: MIT License
Credits: If you use this code, please cite the work ###

Description:
This script computes the atomic contributions to the London dispersion energy as well as the London
dispersion density function based on the D3 dispersion correction, as described in Ref. ####

Prerequisites:
- Ensure the DFT-D3 executable is available at at the specified path (currently, PATHD3 = "./dftd3").
  Download from https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3.
- Python 3.x with standard libraries.
- Please cite the related papers of Prof. Grimme and coworkers

Input/Output:
- Input  : a "{basename}.xyz" file as input is required (coordinates specified in Ångström).
- Output : at the end of the calculation, three output files will be generated:
    - ".d3atomwise.txt" : a file where coordinates are accompanied by a column indicating the contribution to the London dispersion energy for each atom.
    - ".d3out.txt"      : the original DFT-D3 output (J. Chem. Phys. 132, 154104 (2010); DOI: 10.1063/1.3382344).
    - ".d3omega.cube"   : a function of spatial coordinates, facilitating visualization and analysis of dispersion energy contributions.

Usage   : python3 lddensityd3.py <basename> [--npoints NPOINTS] [--func FUNC] [--damp DAMP] [--nprocs NPROCS] 
Example : python3 lddensityd3.py water --npoints 80 --func b3-lyp --damp bj --nprocs 2

Arguments Explanation:
- basename : Base name for the input .xyz file. It's the only required argument without a default value.
- npoints  : The number of grid points for each dimension. Optional, defaults to 80.
- func     : Specifies the functional parameters used in the D4 dispersion calculation. Optional, defaults to 'b3-lyp'.
- damp     : Specifies the damping scheme used in the D3 dispersion calculation. Optional, defaults to 'bj'.
- nprocs   : The number of processors to use for parallel calculations. Optional, defaults to 1.
"""

def read_xyz(xyz):
    """
    Reads .xyz file and extract atomic coordinates
    
    Args:
        xyz (str): The path of the .xyz file to read.
    
    Returns:
        atoms (list): A list containing the atomic symbols.
        x (list): A list containing the x-coordinates of the atoms.
        y (list): A list containing the y-coordinates of the atoms.
        z (list): A list containing the z-coordinates of the atoms.
        nat (int): The number of atoms.
    """
    atoms = []
    x, y, z = [], [], []

    with open(xyz, 'r', encoding="utf-8") as fp:
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atoms.append(data[0])
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            nat = len(atoms)

    return atoms, x, y, z, nat

def box_gen(xmin, ymin, zmin, xstep, ystep, zstep, npoints):
    """
    Generates a 3D grid of points for computing the London dispersion density function
    
    Args:
        xmin (float): Minimum x-coordinate of the grid.
        ymin (float): Minimum y-coordinate of the grid.
        zmin (float): Minimum z-coordinate of the grid.
        xstep (float): Step size in the x-direction.
        ystep (float): Step size in the y-direction.
        zstep (float): Step size in the z-direction.
        npoints (int): Number of points along each axis.
        
    Returns:
        xyzbox (list): A list of 3D points in the grid.
    """
    xyzbox = []
    for p in range(npoints):
        xpt = xmin + p * xstep
        for q in range(npoints):
            ypt = ymin + q * ystep
            for r in range(npoints):
                zpt = zmin + r * zstep
                xyzbox.append((xpt,ypt,zpt))
    return xyzbox

def omega_comp_wrapper(args):
    """
    Computes the London dispersion density at a single grid point

    Args:
        args (tuple): A tuple containing the following:
        - boxpt (list): XYZ coordinates of the point where omega is to be computed.
        - xyzcoord_au (list): XYZ coordinates of atoms in atomic unit.
        - atwdisp (list): List of atomic-wise dispersion.
        
    Returns:
        omegaval (float): The computed London dispersion density value.
    """
    boxpt, xyzcoord_au, atwdisp = args
    omegaval, ith = 0.0, 0
    a = 0.5
    norm = (1.0 / (math.pi / a)) ** 1.5
    for syspt in xyzcoord_au:
        rbox_rith = (boxpt[0] - syspt[0]) ** 2 + \
                    (boxpt[1] - syspt[1]) ** 2 + \
                    (boxpt[2] - syspt[2]) ** 2
        exp_norm = norm * math.exp(- a * rbox_rith)
        ld_rbox_rith = atwdisp[ith] * exp_norm
        omegaval += ld_rbox_rith
        ith += 1
    return omegaval

def omega_comp(nproc, ANG_AU, xyzbox, x, y, z, atwdisp):
    """
    Compute the London dispersion density values across the entire grid in parallel

    Args:
        nproc (int): The number of processors to use for parallel calculations.
        ANG_AU (float): The conversion factor from Angstroem to atomic units.
        xyzbox (list): A list of 3D grid points for computing the omega function.
        x (list): A list containing the x-coordinates of atoms in Angstroem.
        y (list): A list containing the y-coordinates of atoms in Angstroem.
        z (list): A list containing the z-coordinates of atoms in Angstroem.
        atwdisp (list): List of atomic-wise dispersion.

    Returns:
        omega (list): A list containing the computed London dispersion density values for each grid point.
    """
    xyzcoord_au = [[x[i] * ANG_AU, y[i] * ANG_AU, z[i] * ANG_AU]  for i in range(len(x))]

    with Pool(nproc) as pool:
        args = [(boxpt, xyzcoord_au, atwdisp) for boxpt in xyzbox]
        results = pool.map(omega_comp_wrapper, args)
    omega = results
    return omega

def output(basename, func, damp, Esyskcal, omegaintegral, atwdisptot):
    """
    Prints DFT-D3 atom-wise analysis information for the given molecule.

    Args:
        basename (str): The base name of the molecule file.
        func (str): The functional used in the calculation.
        damp (str): The damping function used in the calculation.
        Esyskcal (float): The total dispersion energy in kcal/mol.
        omegaintegral (float): The integral of the London dispersion density function.
        atwdisptot (float): The total atom-wise contribution to dispersion energy.

    Returns:
        None
    """
    separator = '=' * 60
    sub_separator = '-' * 60
    
    print(f"\n{separator}")
    print(f"DFT-D3 atom-wise analysis info for {basename}.xyz")
    print(f"{separator}")
    print("Functional:".ljust(20) + f"{func}")
    print("Damping:".ljust(20) + f"{damp}\n")
    print(f"{sub_separator}")
    print("Total Dispersion Energy:".ljust(40) + f"{Esyskcal:10.3f} kcal/mol")
    print(f"{sub_separator}\n")
    print("Total atom-wise contribution:".ljust(40) + f"{atwdisptot:10.3f} kcal/mol")
    print("Integral of Omega Function:".ljust(40) + f"{omegaintegral:10.3f}")
    print(f"\n{separator}\n")

if __name__ == "__main__":
    import os
    import subprocess
    from multiprocessing import Pool
    import sys
    import time
    import math
    import argparse

    start_time = time.time()

    # Define constants for unit conversion
    ANG_AU = 1.0/0.5291772083
    AU_ANG = 0.5291772083

    elements = [None,
         "H", "He",
         "Li", "Be",
         "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg",
         "Al", "Si", "P", "S", "Cl", "Ar",
         "K", "Ca",
         "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
         "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr",
         "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
         "In", "Sn", "Sb", "Te", "I", "Xe",
         "Cs", "Ba",
         "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
         "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
         "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra",
         "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
         "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub"]
    
    # Initialize command-line argument parser and define expected arguments
    parser = argparse.ArgumentParser(description="Compute atomic contributions to the London dispersion energy using DFT-D3")
    parser.add_argument('basename', type=str, help='Base name for the input .xyz file')
    parser.add_argument('--npoints', type=int, default=80, help='Number of grid points for each dimension (default: 80)')
    parser.add_argument('--func', type=str, default='b3-lyp', help='Functional parameters used in D4 dispersion calculation (default: b3-lyp)')
    parser.add_argument('--damp', type=str, default='bj', help='Damping scheme used in the D3 dispersion calculation (default: bj)')
    parser.add_argument('--nprocs', type=int, default=1, help='Number of processors for parallel calculations (default: 1)')
    
    args = parser.parse_args()

    # Extract arguments into variables
    basename = args.basename
    npoints = args.npoints
    func = args.func
    damp = args.damp
    nproc = args.nprocs

    # Construct file names
    xyz = f"{basename}.xyz"
    atw = f"{basename}.{npoints}.{func}.{damp}.d3atomwise.txt"
    d3out = f"{basename}.{npoints}.{func}.{damp}.d3out.txt"
    omegaout = f"{basename}.{npoints}.{func}.{damp}.d3omega.cube"

    if not os.path.isfile(xyz):
        sys.exit("Could not find the .xyz file")

    atoms, x, y, z, natoms = read_xyz(xyz)

    EXTENT = 7.0
    xmin = min(x) * ANG_AU - EXTENT
    xmax = max(x) * ANG_AU + EXTENT
    ymin = min(y) * ANG_AU - EXTENT
    ymax = max(y) * ANG_AU + EXTENT
    zmin = min(z) * ANG_AU - EXTENT
    zmax = max(z) * ANG_AU + EXTENT

    # Set the path to the DFT-D3 executable
    PATHD3 = "./dftd3"

    # Execute DFT-D3 calculation with specified parameters
    subprocess.run([f"{PATHD3} "
                    f"{xyz} -func {func} -{damp} -anal > "
                    f"{d3out}"], shell=True, check=True)

    with open(f"{d3out}", "r", encoding="utf-8") as d3realsys:
        d3rs = d3realsys.readlines()
        # generate list to store pairwise terms:
        atwdisp = [0.000000] * natoms
        # calculate the total number of possible pairs
        COMBO = int(natoms * (natoms - 1) / 2)

        for j, line in enumerate(d3rs):
            if 'Edisp /kcal,au:' in line:
                erealsys = line.strip().split()
                Esyskcal = float(erealsys[2])
                Esysau = float(erealsys[3])

            if 'analysis of pair-wise terms (in kcal/mol)' in line:
                for l in range(j + 2, j + 2 + COMBO):
                    epair = d3rs[l].strip().split()
                    if len(epair) >= 9:
                        iat, jat, eij = int(epair[0]), int(epair[1]), float(epair[-1])
                        atwdisp[iat - 1] += eij / 2.00000
                        atwdisp[jat - 1] += eij / 2.00000
                    else:
                        break

    # Generate the atomwise txt file where the atomic contributions are listed
    atwdisptot = 0
    with open(f"{atw}", "w", encoding="utf-8") as pwout:
        pwout.write("analysis of atom-wise contributions (in kcal/mol)\n")
        pwout.write("          X            Y           Z         Edisp\n")
        for el, atwdisp_value in enumerate(atwdisp):
            pwout.write(f"{atoms[el]} {x[el]:12.6f} {y[el]:12.6f}"
                        f"{z[el]:12.6f}{atwdisp_value:12.6f}\n")
            atwdisptot += atwdisp_value

    with open(f"{omegaout}", "w", encoding="utf-8") as fp:

        fp.write(f"LD density ({npoints} grid points)\n")
        fp.write(f"input file: {basename}.xyz ({func} {damp})\n")
        fp.write(f"{len(atoms):5d}{xmin:12.6f}{ymin:12.6f}{zmin:12.6f}\n")

        # calculate step sizes for each dimension
        xstep = (xmax - xmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{xstep:12.6f}{0:12.6f}{0:12.6f}\n")
        ystep = (ymax - ymin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{ystep:12.6f}{0:12.6f}\n")
        zstep = (zmax - zmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{0:12.6f}{zstep:12.6f}\n")

        # write the atomic number of each atom and its corresponding coordinates
        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            xi, yi, zi = x[i] * ANG_AU, y[i] * ANG_AU, z[i] * ANG_AU
            fp.write(f"{index:5d}{index:12.6f}{xi:12.6f}{yi:12.6f}{zi:12.6f}\n")

        xyzbox = box_gen(xmin, ymin, zmin, xstep, ystep, zstep, npoints)

        omega = omega_comp(nproc, ANG_AU, xyzbox, x, y, z, atwdisp)
        
        m = 0
        n = 0
        num = 0
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    fp.write(f"{omega[num]:14.5e}")
                    m += 1
                    n += 1
                    num += 1
                    if n > 5:
                        fp.write("\n")
                        n = 0
                if n != 0:
                    fp.write("\n")
                    n = 0
        fp.write(" ")

        # calculate integral of LD density function
        omegaintegral = sum(omega) * xstep * ystep * zstep
  
        output(basename, func, damp, Esyskcal, omegaintegral, atwdisptot)

exetime = time.time() - start_time
print(f"execution time: {exetime:10.5f} seconds\n")
