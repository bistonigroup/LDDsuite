"""
(c) 2025, Gianluca Regni, Lorenzo Baldinelli, Giovanni Bistoni
License: LGPL-3.0 with restriction (see LICENSE)
Credits: If you use this code, please cite the work(s)
DOI: 10.1021/acscentsci.5c00356; DOI: 10.1021/acs.jctc.3c00977

Description:
The method is based on the Atomic Decomposition of London Dispersion energy (ADLD).
The script computes the atomic contributions to the London dispersion energy as well as the London
dispersion density function based on the D3/D4 dispersion correction, as described in Ref.
ACS Cent. Sci. 2025, 11, 6, 890–898;
J. Chem. Theory Comput. 2024, 20, 5, 1923–1931
Then generate a .cube file to visualize these contributions.

Prerequisites:
- Ensure the DFT-D3/DFT-D4 executable is available at at the specified path (in the current folder by default).
  Download from
  https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3 or
  https://www.chemie.uni-bonn.de/grimme/de/software/dft-d4
- Python 3.x with standard libraries.
- Please cite the related papers of Prof. Grimme and coworkers

Input/Output:
- Input  : a "{basename}.xyz" file as input is required (coordinates specified in Ångström).
- Output : at the end of the calculation, three output files will be generated:
    - ".d{n}atomwise.txt" : a file where coordinates are accompanied by a column indicating the contribution to the London dispersion energy for each atom.
    - ".d{n}out.txt"      : the original DFT-D3/DFT-D4 output (J. Chem. Phys. 132, 154104 (2010) DOI: 10.1063/1.3382344; J. Chem. Phys. 147, 034112 (2017) DOI: 10.1063/1.4993215).
    - ".d{n}omega.cube"   : a file containing volumetric data to easily visualize the atomic dispersion energy contributions.

Usage   : python3 lddensity.py <basename> [--d {3,4}] [--npoints NPOINTS] [--func FUNC] [--damp DAMP] [--charge CHARGE] [--s9 S9] [--nprocs NPROCS] 
Example : python3 lddensity.py water --d 4 --npoints 80 --func b3-lyp --damp bj --charge 0 --s9 1.0 --nprocs 2
          python3 lddensity.py water --d 3 --npoints 80 --func b3-lyp --damp bj --nprocs 2
          python3 lddensity.py water --onlycube --npoints 80 --nprocs 2

Arguments Explanation:
- basename : Base name for the input .xyz file. It's the only required argument without a default value.
- npoints  : The number of grid points for each dimension. Optional, defaults to 80.
- func     : Specifies the functional to be used in the calculation. Optional, defaults to 'b3-lyp'.
- damp     : Specifies the damping function to be used in the calculation calculation. Optional, defaults to 'bj'.
- charge   : The overall charge of the molecule (only D4). Optional, defaults to 0.
- s9       : Coefficient for ATM (Axilrod-Teller-Muto) three-body dispersion. To eliminate the three-body contribution, set this to 0. Optional, defaults to 1.0. 
- nprocs   : The number of processors to use for parallel calculations. Optional, defaults to 1.
- onlycube : Flag to generate only the .cube file and skip other outputs. Note that {basename}.atomwise.txt file is required.
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
            atom = str(data[0]).capitalize()
            atoms.append(atom)
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            nat = len(atoms)

    return atoms, x, y, z, nat

def read_xyzdisp(xyzdisp):
    
    """
    Reads geometric coordinates and the atomic contribution/dispersion
    
    Args:
        xyz (str): The path of the .atomwise.txt file to read.
    
    Returns:
        atoms (list): A list containing the atomic symbols.
        x (list): A list containing the x-coordinates of the atoms.
        y (list): A list containing the y-coordinates of the atoms.
        z (list): A list containing the z-coordinates of the atoms.
        nat (int): The number of atoms.
        atwdisp (list): A list containing atomwise contributon/dispersion.
    """
    atoms = []
    x, y, z = [], [], []
    atwdisp = []
    
    print(f"INFO: Be sure atomic contributions in {xyzdisp} are in kcal/mol.")

    with open(xyzdisp, 'r', encoding='utf-8') as fp:
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atom = str(data[0]).capitalize()
            atoms.append(atom)
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            atwdisp.append(float(data[4]))
            nat = len(atoms)

    return atoms, x, y, z, nat, atwdisp

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
        ANG_AU (float): The conversion factor from angstroems to atomic units.
        xyzbox (list): A list of 3D grid points for computing the omega function.
        x (list): A list containing the x-coordinates of atoms in angstroems.
        y (list): A list containing the y-coordinates of atoms in angstroems.
        z (list): A list containing the z-coordinates of atoms in angstroems.
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

def output(omegaintegral, atwdisptot, Esyskcal=None):
    """
    Print results of the calculation.

    Args:
        Esyskcal (float): The total dispersion energy in kcal/mol.
        omegaintegral (float): The integral of the London dispersion density function.
        atwdisptot (float): The total atom-wise contribution to dispersion energy.

    Returns:
        None
    """
    
    if not args.onlycube:
        print("Total Dispersion Energy:".ljust(40) + f"{Esyskcal:10.1f} kcal/mol")
    print("Total atom-wise contribution:".ljust(40) + f"{atwdisptot:10.1f} kcal/mol")
    print("Integral of Omega Function:".ljust(40) + f"{omegaintegral:10.1f}")
    
def format_time(seconds):
    """
    Format a time duration in seconds into HH:MM:SS string.

    Args:
    seconds (int or float): total number of seconds to format

    Returns:
    time_str (str) : string representing the time in "HH:MM:SS" format
    """
    
    hours   = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs    = seconds % 60

    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

if __name__ == "__main__":
    import os
    import subprocess
    from multiprocessing import Pool
    import sys
    import time
    import math
    import argparse
    import numpy as np

    start_time = time.time()

    # Define constants for unit conversion
    ANG_AU = 1.0/0.5291772083
    AU_ANG = 0.5291772083
    AU_KCALMOL = 627.5096080305927

    elements = [None,
         "H", "He",
         "Li", "Be", "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
         "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
         "Zn","Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
         "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
         "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
         "Dy", "Ho", "Er", "Tm", "Yb",
         "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
         "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
         "Cf", "Es", "Fm", "Md", "No",
         "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
         "Fl", "Mc", "Lv", "Ts", "Og"]
    
    # Creation of command-line arguments
    parser = argparse.ArgumentParser(description="Generate dispersion density using D3/D4 correction")
    parser.add_argument('basename',  type=str,                     help='Base name for the input .xyz file')
    parser.add_argument('--d',       type=int,   default=4,        help='D3/D4 correction (default: 4)', choices=[3, 4])
    parser.add_argument('--npoints', type=int,   default=80,       help='Number of grid points for each dimension (default: 80)')
    parser.add_argument('--func',    type=str,   default='b3-lyp', help='Functional to be used in the calculation (default: b3-lyp)')
    parser.add_argument('--damp',    type=str,   default='bj',     help='Damping function to be used in calculation (default: bj)')
    parser.add_argument('--charge',  type=int,   default=0,        help='Molecule charge for D4 (default: 0)')
    parser.add_argument('--s9',      type=float, default=0,        help='Coefficient for ATM three-body dispersion for D4 (default: 0). Set 1 or 0 to include or not the three-body contribution in D3.')
    parser.add_argument('--nprocs',  type=int,   default=1,        help='Number of processors for parallel calculations (default: 1)')
    
    parser.add_argument('--onlycube', action='store_true', help='If set, generate only the .cube file. Note that {basename}.atomwise.txt is the input file in this case.')
    
    args = parser.parse_args()

    # Extract arguments into variables
    basename = args.basename
    d        = args.d
    npoints  = args.npoints
    func     = args.func
    damp     = args.damp
    nproc    = args.nprocs
    charge   = args.charge
    s9       = args.s9
    
    # set three-body calculation
    if s9 == 0:
        abc = False
    else:
        abc = True
    # if d == 3 and (0.0 < s9 < 1.0):
    #     print(f"WARNING: Requested D3 with a fractional s9. Three-body calculation will not be performed.")
    #     abc = False
    
    if d == 3 and s9 != 0:
        print(f"WARNING: Three-body calculation is not yet supported for D3.")
        abc = False
    
    if d == 4 and damp != "bj":
        print("WARNING: damping function for D4 must be Becke-Johnson. Setting damp as 'bj'")
        damp = 'bj'

    if not args.onlycube:
               
        # Construct file names
        xyz      = f"{basename}.xyz"
        atw      = f"{basename}.{npoints}.{func}.{damp}.d{d}atomwise.txt"
        dout     = f"{basename}.{npoints}.{func}.{damp}.d{d}out.txt"
        omegaout = f"{basename}.{npoints}.{func}.{damp}.d{d}omega.cube"

        # Input check
        if not os.path.isfile(xyz):
            sys.exit("ERROR: could not find the .xyz file")
        
        # Settings summary  
        print("===SETTINGS SUMMARY===")
        print(f"Input File: {xyz}")
        print(f"D:          {d}")
        print(f"N Points:   {npoints}")
        print(f"Functional: {func}")
        print(f"Damping:    {damp}")
        print(f"Three-Body: {abc}")
        if d == 4:
            print(f"Charge:     {charge}")
        print(f"CPU:        {nproc}")
        print("======================")

        # read input
        atoms, x, y, z, natoms = read_xyz(xyz)

        # Set path
        if   d == 3:
            PATHD = "./dftd3"
        elif d == 4:
            PATHD = "./dftd4"
        
        # Check if path exists
        if not os.path.exists(PATHD):
            sys.exit(f"ERROR: directory not valid. Please make sure DFT-D binary is placed in '{PATHD}'.")
            
        # Compute dispersion
        if   d == 3:
            if charge != 0:
                    print(f"WARNING: A molecular charge of {charge} was set, but D3 is not sensitive to molecular charge")
            # Grimme's code give weird errors with -abc flag
            # if s9 == 1.0:
            #     print("NOTE: Three-body term will be included in the D3 calculation")
            #     subprocess.run([f"{PATHD} "
            #             f"{xyz} -func {func} -{damp} -anal -abc > "
            #             f"{dout}"], shell=True, check=True)
            # else:
            subprocess.run([f"{PATHD} "
                    f"{xyz} -func {func} -{damp} -anal> "
                    f"{dout}"], shell=True, check=True)
        elif d == 4:
            subprocess.run([f"{PATHD} "
                            f"-f {func} -c {charge} --pair-resolved "
                            f"--property --noedisp --mbdscale {s9} {xyz} "
                            f"> {dout}"], shell=True, check=True)


        with open(f"{dout}", "r", encoding="utf-8") as drealsys:
            drs = drealsys.readlines()
            # generate list to store pairwise terms:
            atwdisp = [0.000000] * natoms
            if d == 3:
                # calculate the total number of possible pairs
                COMBO = int(natoms * (natoms - 1) / 2)

                for j, line in enumerate(drs):
                    if 'Edisp /kcal,au:' in line:
                        erealsys = line.strip().split()
                        Esyskcal = float(erealsys[2])
                        Esysau = float(erealsys[3])

                    if 'analysis of pair-wise terms (in kcal/mol)' in line:
                        for l in range(j + 2, j + 2 + COMBO):
                            epair = drs[l].strip().split()
                            if len(epair) >= 9:
                                iat, jat, eij = int(epair[0]), int(epair[1]), float(epair[-1])
                                atwdisp[iat - 1] += eij / 2.00000
                                atwdisp[jat - 1] += eij / 2.00000
                            else:
                                break
            elif d == 4:
                # calculate the total number of possible pairs
                COMBO = int(natoms * natoms)
                
                for j, line in enumerate(drs):
                    if 'Dispersion energy:' in line:
                        erealsys = line.strip().split()
                        Esyseh = float(erealsys[-2])
                        Esyskcal = Esyseh * AU_KCALMOL

                    if 'Pairwise representation of dispersion (in kcal/mol):' in line:
                        for l in range(j + 4, (j + 4 + COMBO)):
                            epair = drs[l].strip().split()
                            if len(epair) >= 11:
                                iat, jat, eij = int(epair[0]), int(epair[3]), float(epair[-1])
                                atwdisp[iat - 1] += eij / 2.00000
                                atwdisp[jat - 1] += eij / 2.00000
                            else:
                                break


        # Generate the .d{n}atomwise.txt file where the atomic contributions are listed
        atwdisptot = 0
        with open(f"{atw}", "w", encoding="utf-8") as pwout:
            pwout.write("analysis of atom-wise contributions (in kcal/mol)\n")
            pwout.write("          X            Y           Z         Edisp\n")
            for el, atwdisp_value in enumerate(atwdisp):
                pwout.write(f"{atoms[el]} {x[el]:12.6f} {y[el]:12.6f}"
                            f"{z[el]:12.6f}{atwdisp_value:12.6f}\n")
                atwdisptot += atwdisp_value

    else: # Only cube file required
        # Construct file names
        xyzdisp  = f"{basename}.atomwise.txt"
        omegaout = f"{basename}.{npoints}.omega.cube"
        
        # Input check
        if not os.path.isfile(xyzdisp):
            sys.exit("ERROR: could not find {basename}.atomwise.txt file")
    
        # Settings summary  
        print("===SETTINGS SUMMARY===")
        print(f"Input File: {xyzdisp}.atomwise.txt")
        print(f"N Points:   {args.npoints}")
        print(f"CPU:        {args.nprocs}")
        print("======================")
        
        # read input
        atoms, x, y, z, natoms, atwdisp = read_xyzdisp(xyzdisp)
        
        atwdisptot = np.sum(atwdisp)
        Esyskcal = None
        
    EXTENT = 7.0
    xmin = min(x) * ANG_AU - EXTENT
    xmax = max(x) * ANG_AU + EXTENT
    ymin = min(y) * ANG_AU - EXTENT
    ymax = max(y) * ANG_AU + EXTENT
    zmin = min(z) * ANG_AU - EXTENT
    zmax = max(z) * ANG_AU + EXTENT
    
    # Generate the .d{n}omega.cube file to visualize the London dispersion density
    with open(f"{omegaout}", "w", encoding="utf-8") as fp:

        fp.write(f"LD density ({npoints} grid points)\n")
        if not args.onlycube:
            if d == 3:
                fp.write(f"input file: {basename}.xyz ({func} {damp}, E(3)={abc})\n")
            elif d == 4:
                fp.write(f"input file: {basename}.xyz ({func} bj, charge={charge}, s9={s9})\n")
        else:
            fp.write(f"input file: {basename}.atomwise.txt\n")
        fp.write(f"{len(atoms):5d}{xmin:12.6f}{ymin:12.6f}{zmin:12.6f}\n")

        # calculate step sizes for each dimension
        xstep = (xmax - xmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{xstep:12.6f}{0:12.6f}{0:12.6f}\n")
        ystep = (ymax - ymin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{ystep:12.6f}{0:12.6f}\n")
        zstep = (zmax - zmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{0:12.6f}{zstep:12.6f}\n")

        # write the atomic number, nuclear charge and the corresponding coordinates for each atom
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
    print("")
    output(omegaintegral, atwdisptot, Esyskcal)

exetime = int(round(time.time() - start_time))
print(f"\nexecution time: {format_time(exetime)} (hh:mm:ss)\n")
