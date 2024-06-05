"""
(c) 2024, Gianluca Regni, Lorenzo Baldinelli, Filippo De Angelis, Giovanni Bistoni
License: MIT License
Credits: If you use this code, please cite the work ###

Description:
This script computes the London dispersion density difference function, which is described in Ref. ###

Prerequisites:
- Python 3.x with standard libraries.

Input/Output:
- Input  : A "{basename}.atomwise.txt" file as input is required.
           This contains the spatial coordinates of a chemical system in XYZ format and,
           as a 5th column, the difference between the atomic LD energy of each atom minus
           that obtained for a different molecular structure.  
           For example, the atomic LD energies can be obtained using the lddensity.py tool.
           However, this tool can be used in principle in conjuction with any arbitrary atom-wise
           decomposition scheme.
- Output : At the end of the calculation, one output file will be generated:
    - ".omega.cube": The London dispersion density difference function in cube format.

Usage   : python3 lddensitydifference.py <basename> [--npoints NP] [--nprocs NPROCS]
Example : python3 lddensitydifference.py water --npoints 80 --nprocs 2

Arguments Explanation:
- basename : Base name for the input ".atomwise.txt" file. It's the only required argument without a default value.
- npoints  : The number of grid points for each dimension. Optional, defaults to 80. 
- nprocs   : The number of processors to use for parallel calculations. Optional, defaults to 1.
"""

def box_gen(xmin, ymin, zmin, xstep, ystep, zstep, npoints):
    """
    Generates a 3D grid of points for computing the LD density difference function
    
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
    Computes the LD density difference function at a single grid point

    Args:
        args (tuple): A tuple containing the following:
        - boxpt (list): XYZ coordinates of the point where omega is to be computed.
        - xyzcoord_au (list): XYZ coordinates of atoms in atomic unit.
        - atwdisp (list): List of atomic-wise dispersion.
        
    Returns:
        omegaval (float): The computed value of LD density difference function.
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

# Compute the omega values across the entire grid in parallel
def omega_comp(nproc, ANG_AU, xyzbox, x, y, z, atwdisp):
    """
    Compute the values of LD density different function across the entire grid in parallel

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

def output(basename, omegaintegral, atwdisptot):
    """
    Prints analysis information

    Args:
        basename (str): The base name of the molecule file.
        omegaintegral (float): The integral of the London dispersion density function.
        atwdisptot (float): The total atom-wise contribution to dispersion energy.

    Returns:
        None
    """
    separator = '=' * 60
    
    print(f"\n{separator}")
    print(f"Dispersion difference analysis info for {basename}.atomwise.txt")
    print(f"{separator}\n")
    print("Total atom-wise contribution:".ljust(40) + f"{atwdisptot:10.3f} kcal/mol")
    print("Integral of Omega Function:".ljust(40) + f"{omegaintegral:10.3f}")
    print(f"\n{separator}\n")

if __name__ == "__main__":
    import os
    from multiprocessing import Pool
    import sys
    import time
    import math
    import argparse

    start_time = time.time()

    # Define constants for unit conversion
    ANG_AU = 1.0 / 0.5291772083
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
    parser = argparse.ArgumentParser(description="Compute the values of London dispersion density difference function")
    parser.add_argument('basename', type=str, help='Base name for the input .atomwise.txt file')
    parser.add_argument('--npoints', type=int, default=80, help='Number of grid points for each dimension (default: 80)')
    parser.add_argument('--nprocs', type=int, default=1, help='Number of processors for parallel calculations (default: 1)')
    
    args = parser.parse_args()

    # Extract arguments into variables
    basename = args.basename
    npoints = args.npoints
    nproc = args.nprocs

    # Construct file names
    atwinp = f"{basename}.atomwise.txt"
    omegaout = f"{basename}.omega.cube"

    if not os.path.isfile(atwinp):
        sys.exit("Could not find the *.atomwise.txt file")

    # Parses atomic coordinates and dispersion data from the input file
    with open(atwinp, "r", encoding="utf-8") as atwread:
        atw = atwread.readlines()
        nlines = len(atw)
        at, xat, yat, zat, atwdisp = [], [], [], [], []
        for line in range(2, nlines):
            atom, x, y, z, atwd3 = atw[line].strip().split()
            at.append(str(atom))
            xat.append(float(x))
            yat.append(float(y))
            zat.append(float(z))
            atwdisp.append(float(atwd3))
    atwdisptot = sum(atwdisp)

    EXTENT = 7.0
    xmin = min(xat) * ANG_AU - EXTENT
    xmax = max(xat) * ANG_AU + EXTENT
    ymin = min(yat) * ANG_AU - EXTENT
    ymax = max(yat) * ANG_AU + EXTENT
    zmin = min(zat) * ANG_AU - EXTENT
    zmax = max(zat) * ANG_AU + EXTENT


    with open(f"{omegaout}", "w", encoding="utf-8") as fp:
        
        fp.write(f"LD difference density ({npoints} grid points)\n")
        fp.write(f"input file: {basename}.atomwise.txt\n")
        fp.write(f"{len(at):5d}{xmin:12.6f}{ymin:12.6f}{zmin:12.6f}\n")

        # calculate step sizes for each dimension
        xstep = (xmax - xmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{xstep:12.6f}{0:12.6f}{0:12.6f}\n")
        ystep = (ymax - ymin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{ystep:12.6f}{0:12.6f}\n")
        zstep = (zmax - zmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{0:12.6f}{zstep:12.6f}\n")

        # write the atomic number of each atom and its corresponding coordinates
        for i, atom in enumerate(at):
            index = elements.index(atom)
            xi, yi, zi = xat[i] * ANG_AU, yat[i] * ANG_AU, zat[i] * ANG_AU
            fp.write(f"{index:5d}{index:12.6f}{xi:12.6f}{yi:12.6f}{zi:12.6f}\n")

        xyzbox = box_gen(xmin, ymin, zmin, xstep, ystep, zstep, npoints)

        omega = omega_comp(nproc, ANG_AU, xyzbox, xat, yat, zat, atwdisp)
        
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
        
        output(basename, omegaintegral, atwdisptot)

exetime = time.time() - start_time
print(f"execution time: {exetime:10.5f} seconds\n")
