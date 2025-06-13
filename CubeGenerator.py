"""
(c) 2024, Gianluca Regni
License: LGPL-3.0
Credits: If you use this code, please cite the the author(s)

Description:
This script is developed to generate the .cube file to visualize the atomic London dispersion contributions from the ADLD method.

Prerequisites:
- Python 3.x with standard libraries.

Input/Output:
- Input  : an xyz file with an additional column with the respective atomic contribution in hartree, named "{basename}.atomwise.txt" (coordinates in angstroems).
- Output : the .cube file

Usage   : python3 CubeGenerator.py <basename> [--npoints NP]  [--nprocs NPROCS]
Example : python3 CubeGenerator.py water --npoints 40 --nprocs 2

Arguments Explanation:
- basename : Base name for the input file. It's the only required argument without a default value.
- npoints  : The number of grid points for each dimension. Optional, defaults to 80.
- nprocs   : The number of processors to use for parallel calculations. Optional, defaults to 1.
"""

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

    with open(xyzdisp, 'r', encoding='utf-8') as fp:
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atom = data[0]
            atoms.append(atom)
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            atwdisp.append(float(data[4])*627.5)
            nat = len(atoms)
    atwdisptot = sum(atwdisp)
    print("Total atom-wise contribution:".ljust(40) + f"{atwdisptot:10.3f} kcal/mol")

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

if __name__ == "__main__":
    import os
    from multiprocessing import Pool
    import sys
    import math
    import argparse
    
    # Define constants for unit conversion
    ANG_AU = 1.0/0.5291772083
    AU_ANG = 0.5291772083
    AU_KCALMOL = 627.5096080305927

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
    parser = argparse.ArgumentParser(description="Generator of .cube file, provide the xyz file with the dispersion associated with the atoms next to it")
    parser.add_argument('basename', type=str, help='Base name for the input atomwise.txt file')
    parser.add_argument('--npoints', type=int, default=80, help='Number of grid points for each dimension (default: 80)')
    parser.add_argument('--nprocs', type=int, default=1, help='Number of processors for parallel calculations (default: 1)')
    
    args = parser.parse_args()

    # Extract arguments into variables
    basename = args.basename
    npoints = args.npoints
    nproc = args.nprocs
    
    # Construct file names
    xyzdisp = f"{basename}.atomwise.txt"
    omegaout = f"{basename}.{npoints}.omega.cube"
    
    if not os.path.isfile(xyzdisp):
        sys.exit("Error: could not find the .atomwise.txt file")
        
    atoms, x, y, z, natoms, atwdisp = read_xyzdisp(xyzdisp)
    
    EXTENT = 7.0
    xmin = min(x) * ANG_AU - EXTENT
    xmax = max(x) * ANG_AU + EXTENT
    ymin = min(y) * ANG_AU - EXTENT
    ymax = max(y) * ANG_AU + EXTENT
    zmin = min(z) * ANG_AU - EXTENT
    zmax = max(z) * ANG_AU + EXTENT
    
    with open(f"{omegaout}", "w", encoding="utf-8") as fp:

        fp.write(f"LD density ({npoints} grid points)\n")
        fp.write(f"input file: {basename}.atomwise.txt\n")
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

        # compute integral of LD density difference function
        omegaintegral = sum(omega) * xstep * ystep * zstep
        print("Integral of Omega Function:".ljust(40) + f"{omegaintegral:10.3f}\n")
