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
- Refer to requirements.txt

Input/Output:
- Input  : a "{basename}.xyz" file as input is required (coordinates specified in Ångström).
- Output : at the end of the calculation, three output files will be generated:
    - ".d{n}atomwise.txt" : a file where coordinates are accompanied by a column indicating the contribution to the London dispersion energy for each atom.
    - ".d{n}omega.cube"   : a file containing volumetric data to easily visualize the atomic dispersion energy contributions.

Usage   : python3 lddensity.py <basename> [--d {3,4}] [--npoints NPOINTS] [--func FUNC] [--damp DAMP] [--charge CHARGE] [--abc] [--nprocs NPROCS] 
Example : python3 lddensity.py water --d 4 --npoints 80 --func b3-lyp --damp bj --charge 0 --abc --nprocs 2
          python3 lddensity.py water --d 3 --npoints 80 --func b3-lyp --damp bj --nprocs 2
          python3 lddensity.py water --onlycube --npoints 80 --nprocs 2

Arguments Explanation:
- basename : Base name for the input .xyz file. It's the only required argument without a default value.
- npoints  : The number of grid points for each dimension. Optional, defaults to 80.
- func     : Specifies the functional to be used in the calculation. Optional, defaults to 'b3-lyp'.
- damp     : Specifies the damping function to be used in the calculation calculation. Optional, defaults to 'bj'.
- charge   : The overall charge of the molecule (only D4). Optional, defaults to 0.
- abc      : Include three-body contribution. Optional, defaults to False. 
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
        coords (np.ndarray float64): (N,3) array of atomic coordinates (in angstroem)
        nat (int): The number of atoms.
    """
    if not os.path.exists(xyz):
        raise FileNotFoundError(f"XYZ file '{xyz}' not found")

    with open(xyz, 'r') as fp:
        lines = fp.readlines()
    nat = int(lines[0])
    atoms, coords = [], []
    for line in lines[2:2+nat]:
        el, x, y, z = line.split()
        atoms.append(el.capitalize())
        coords.append([float(x), float(y), float(z)])
    return atoms, np.array(coords), nat

def read_xyzdisp(xyzdisp):
    
    """
    Reads geometric coordinates and the atomic contribution/dispersion
    
    Args:
        xyzdisp (str): The path of the .atomwise.txt file to read.
    
    Returns:
        atoms (list): A list containing the atomic symbols.
        coords (np.ndarray float64): (N,3) array of atomic coordinates (in angstroem)
        nat (int): The number of atoms.
        atwdisp (np.ndarray float64): (N) array containing atomwise contributon/dispersion.
    """
    if not os.path.exists(xyzdisp):
        raise FileNotFoundError(f"file '{xyzdisp}' not found")
    
    print(f"INFO: Be sure atomic contributions in {xyzdisp} are in kcal/mol.")
    
    atoms   = []
    coords  = []
    atwdisp = []
        
    with open(xyzdisp, 'r') as fp:
            for line in fp:
                parts = line.split()
                if len(parts) == 5:
                    try:
                        el, x, y, z, atomdisp = parts
                        atoms.append(el.capitalize())
                        coords.append([float(x), float(y), float(z)])
                        atwdisp.append(float(atomdisp))
                    except ValueError:
                        continue
                    
    nat = len(atoms)
    
    if nat == 0:
        print(f"WARNING: No atomic data found in {xyzdisp}")
        
    return atoms, np.array(coords), nat, np.array(atwdisp)

from numba import njit, prange
@njit(parallel=True)
def omega_comp(xmin, ymin, zmin, xstep, ystep, zstep, npoints, xyzcoord_au, atwdisp):
    """
    Computes the London dispersion density across a 3D grid in parallel.
    To maximize efficiency, grid points are generated on-the-fly, avoiding 
    the creation of large coordinate arrays and optimizing memory usage.

    Args:
    xmin            [float]      : minimum x-coordinate of the grid
    ymin            [float]      : minimum y-coordinate of the grid
    zmin            [float]      : minimum z-coordinate of the grid
    xstep           [float]      : step size along the x-axis
    ystep           [float]      : step size along the y-axis
    zstep           [float]      : step size along the z-axis
    npoints         [int]        : number of points along each axis
    xyzcoord_au     [np.ndarray] : (N,3) array of atomic coordinates (in atomic units)
    atwdisp         [np.ndarray] : length N array of atomic dispersion contributions

    Returns:
    omega           [np.ndarray] : 1D array of shape (npoints**3,) containing computed density values for cube file
    """
    # initilization
    omega = np.zeros(npoints**3)
    
    # constants
    a = 0.5
    norm = (a / np.pi) ** 1.5
    
    # disp density calculation while we create the cube box
    
    for ix in prange(npoints):  # prange is for numba
        vx = xmin + ix * xstep
        for iy in range(npoints):
            vy = ymin + iy * ystep
            for iz in range(npoints):
                vz = zmin + iz * zstep
                
                # unique index for each point
                idx = ix * (npoints**2) + iy * npoints + iz
                
                # dispersion density calculation
                val = 0.0
                for j in range(xyzcoord_au.shape[0]):
                    dx = vx - xyzcoord_au[j, 0]
                    dy = vy - xyzcoord_au[j, 1]
                    dz = vz - xyzcoord_au[j, 2]
                    r2 = dx*dx + dy*dy + dz*dz
                    val += atwdisp[j] * np.exp(-a * r2)
                
                omega[idx] = val * norm
                
    return omega
    
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
    
    import time
    import argparse

    start_time = time.time()

    # Define constants for unit conversion
    ANG_AU     = 1.889725989
    AU_ANG     = 0.5291772083
    AU_KCALMOL = 627.509474

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
    parser.add_argument('--abc',     action='store_true',          help='Compute three-body interactions')
    parser.add_argument('--nprocs',  type=int,   default=1,        help='Number of processors for parallel calculations (default: 1)')
    parser.add_argument('--onlycube', action='store_true', help='If set, generate only the .cube file. Note that {basename}.atomwise.txt is the input file in this case.')
    
    args = parser.parse_args()
    
    import os
    import numpy as np

    # Extract arguments into variables
    basename = args.basename
    d        = args.d
    npoints  = args.npoints
    func     = args.func
    damp     = args.damp
    nproc    = args.nprocs
    charge   = args.charge
    abc      = args.abc
       
    if d == 4 and damp != "bj":
        print("WARNING: damping function for D4 must be Becke-Johnson. Setting damp as 'bj'")
        damp = 'bj'
        
    if (args.d == 3):
        import dftd3.interface as d3
        if (func == "b3-lyp"):
            func = "b3lyp"
    elif (args.d == 4):
        import dftd4.interface as d4

    if not args.onlycube:
               
        # Construct file names
        xyz      = f"{basename}.xyz"
        atw      = f"{basename}.atomwise.txt"
        omegaout = f"{basename}.omega.cube"
       
        # Settings summary  
        print("\n--- Settings ---\n")
        print(f"Input File: {xyz}")
        print(f"D:          {d}")
        print(f"N Points:   {npoints}")
        print(f"Functional: {func}")
        print(f"Damping:    {damp}")
        print(f"Three-Body: {abc}")
        if d == 4:
            print(f"Charge:     {charge}")
        print(f"CPU:        {nproc}")

        # read input
        atoms, coords, natoms = read_xyz(xyz)
            
        # Compute dispersion
        print("\ncomputing dispersion...")
        numbers = np.array([elements.index(a) for a in atoms], dtype=int)
        
        if (args.d == 3):
            if   damp == "bj":
                param = d3.RationalDampingParam(method=func, atm=abc)
            elif damp == "zero":
                param = d3.ZeroDampingParam(method=func, atm=abc)
            elif damp == "bjm":
                param = d3.ModifiedRationalDampingParam(method=func, atm=abc)
            elif damp == "zerom":
                param = d3.ModifiedZeroDampingParam(method=func, atm=abc)
            else:
                raise ValueError(f"Unknown damping method: {damp}")
            model = d3.DispersionModel(numbers, coords/AU_ANG)
        elif (args.d == 4):
                param = d4.DampingParam(method=func, atm=abc)
                model = d4.DispersionModel(numbers, coords/AU_ANG, charge=charge)
            
            
        pair = model.get_pairwise_dispersion(param)
        E2   = pair["additive pairwise energy"]
        
        atwdisp  = np.zeros(natoms)
        atwdisp += np.sum(E2, axis=1)
        
        if (abc):
            E3       = pair["non-additive pairwise energy"]
            atwdisp += np.sum(E3, axis=1)
            
        atwdisp *= AU_KCALMOL
        Esyskcal = np.sum(atwdisp)

        # E2 has the following shape for water molecule:
        # E2: 
        #            O                H              H
        # O  [[ 0.00000000e+00 -2.98587192e-07 -2.98570820e-07]
        # H  [-2.98587192e-07   0.00000000e+00 -5.04702069e-06]
        # H  [-2.98570820e-07  -5.04702069e-06  0.00000000e+00]]
        
        # it is the same for three-body contributons (E3)

        # Generate the .d{n}atomwise.txt file where the atomic contributions are listed
        atwdisptot = 0
        with open(f"{atw}", "w") as pwout:
            pwout.write("analysis of atom-wise contributions (in kcal/mol)\n")
            pwout.write("          X            Y           Z         Edisp\n")
            for el, atwdisp_value in enumerate(atwdisp):
                pwout.write(f"{atoms[el]} {coords[el,0]:12.6f} {coords[el,1]:12.6f}"
                            f"{coords[el,2]:12.6f}{atwdisp_value:12.6f}\n")
                atwdisptot += atwdisp_value

    else: # Only cube file required
        # Construct file names
        xyzdisp  = f"{basename}.atomwise.txt"
        omegaout = f"{basename}.{npoints}.omega.cube"
    
        # Settings summary  
        print("--- Settings ---\n")
        print(f"Input File: {xyzdisp}.atomwise.txt")
        print(f"N Points:   {args.npoints}")
        print(f"CPU:        {args.nprocs}")
        
        # read input
        atoms, coords, natoms, atwdisp = read_xyzdisp(xyzdisp)
        
        atwdisptot = np.sum(atwdisp)
        Esyskcal = None
        
    print("computing dispersion density...")   
    # we need coordinates in bohr
    coords = coords * ANG_AU
    
    # define the grid boundaries. An extra padding (EXTENT) is added
    EXTENT = 7.0
    xmin, ymin, zmin = np.min(coords, axis=0) - EXTENT
    xmax, ymax, zmax = np.max(coords, axis=0) + EXTENT
    
    # compute step sizes along each axis
    xstep = (xmax - xmin) / float(npoints - 1)
    ystep = (ymax - ymin) / float(npoints - 1)
    zstep = (zmax - zmin) / float(npoints - 1)
    
    with open(omegaout, "w") as fp:
        # header
        fp.write(f"LD density ({npoints} grid points)\n")
        fp.write(f"input file: {basename}.xyz ({func} {damp}, E(3)={abc})\n")
        fp.write(f"{len(atoms):5d}{xmin:12.6f}{ymin:12.6f}{zmin:12.6f}\n")
        fp.write(f"{npoints:5d}{xstep:12.6f}{0:12.6f}{0:12.6f}\n")
        fp.write(f"{npoints:5d}{0:12.6f}{ystep:12.6f}{0:12.6f}\n")
        fp.write(f"{npoints:5d}{0:12.6f}{0:12.6f}{zstep:12.6f}\n")

        # list of atoms in the header
        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            xi, yi, zi = coords[i]
            fp.write(f"{index:5d}{index:12.6f}{xi:12.6f}{yi:12.6f}{zi:12.6f}\n")

        omega = omega_comp(
            xmin, ymin, zmin, xstep, ystep, zstep, 
            npoints, coords, atwdisp
        )
                    
        # body
        text = []
        n = 0
        num = 0
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    text.append(f"{omega[num]:14.5e}")
                    n += 1
                    num += 1
                    if n > 5:
                        text.append("\n")
                        n = 0
                if n != 0:
                    text.append("\n")
                    n = 0
            fp.write("".join(text))
            text = []
        fp.write(" ")

    # calculate integral of LD density function
    omegaintegral = sum(omega) * xstep * ystep * zstep
    
    print("\n--- Results ---\n")
    if not args.onlycube or Esyskcal!=None:
        print("Total Dispersion Energy:".ljust(40)  + f"{Esyskcal:10.1f} kcal/mol")
    print("Total atom-wise contribution:".ljust(40) + f"{atwdisptot:10.1f} kcal/mol")
    print("Integral of Omega Function:".ljust(40)   + f"{omegaintegral:10.1f} kcal/mol")
    
if not args.onlycube:
    print("\n--- References for this run ---")
    print("""
ACS Cent. Sci. 2025, 11, 6, 890–898;
J. Chem. Theory Comput. 2024, 20, 5, 1923–1931;""")
    if d == 4:
        print("""J. Chem. Phys., 2017, 147, 034112;
J. Chem. Phys., 2019, 150, 154122;"""
        )
    else:
        print("""J. Open Source Softw., 2024, 9(103), 7169;
J. Chem. Phys., 2010, 132, 154104;
""")
    

exetime = int(round(time.time() - start_time))
print(f"\nexecution time: {format_time(exetime)} (hh:mm:ss)\n")
