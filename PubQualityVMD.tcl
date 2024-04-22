# VMD Script for Publication-Quality Images
# Author: Gianluca Regni
# Copyright (c) 2024 Gianluca Regni
# License: MIT (https://opensource.org/licenses/MIT)
# Credits: Please cite the author and the References if you use or reference this program.

# References:
# VMD manual: https://www.ks.uiuc.edu/Research/vmd/current/ug.pdf
# Shiv material settings adapted from: https://shivupa.github.io/blog/making-publication-quality-images-with-vmd/

# Description:
# This script is designed to visualize LD density (difference) function and generate publication-quality
# images using VMD (Visual Molecular Dynamics). 
# Additionally, it includes a re-designed pick event feature that also displays atomic contributions to dispersion
# energy in the terminal, extracted from the atomwise.txt file.
# It is designed to be used in combination with other programs of LD Density suite ####.

# Usage:
# vmd -e PubQualityVMD.tcl 

# ------------------------------------
# File Name Declarations
# ------------------------------------

puts "File provided:"
set file_name "water.xyz"
puts "File extension:"
set file_extension [string range [file extension $file_name] 1 end]

if {$file_extension eq "cube"} {
    puts "Base name:"
    set base_name [string range [file rootname $file_name] 0 [expr {[string length [file rootname $file_name]] - 6}]]
    puts $base_name

    # Check if the atomwise file is provided
    puts "Atomwise file name:"
    set atomwise_file "${base_name}atomwise.txt"
    puts "${base_name}atomwise.txt"
    set has_atomwise_file [file exists $atomwise_file]

    if {!$has_atomwise_file} {
        puts "Warning: atomwise.txt file not found. Atomic contributions to dispersion energy will not be displayed."
    }
}

# ------------------------------------
# Procedure Declarations
# ------------------------------------

proc handle_atom_pick {} {
    global atom_info
    global vmd_pick_atom
    
    foreach info $atom_info {
        lassign $info atom_ID atomName x y z edisp
        
        if {$atom_ID == $vmd_pick_atom} {
            puts "You have selected atom $atomName (ID: $atom_ID)"
            puts "Coordinates:  x = $x,  y = $y,  z = $z"
            puts "Dispersion energy (kcal/mol): $edisp"
            break
        }
    }
}

# ------------------------------------
# Main Script
# ------------------------------------

# Load Molecule
mol new "${file_name}" type ${file_extension}

# Display Settings
display rendermode GLSL
display projection Orthographic
display resetview
display resize 800 600
display ambientocclusion on
display shadows on
display height 5
display resize 1600 1200
display reposition 800 1600

# Set the background color to white
color Display Background white

# Material Customization
material add Shiv
material change ambient Shiv 0.00
material change diffuse Shiv 0.85
material change specular Shiv 0.00
material change shininess Shiv 0.00
material change mirror Shiv 0.00
material change opacity Shiv 1
#material change outline Shiv 0.00
#material change outlinewidth Shiv 0.00

# Rotation and Scaling
#rotate y by 135
#rotate z by 45
#scale by 1.3

# Molecular Representation
color Element C gray

mol delrep 0 top
mol addrep 0
mol modcolor 0 0 Element
mol modstyle 0 0 CPK 0.8 0.3 72 72
if {$file_extension eq "cube"} {
    mol modcolor 0 0 Volume 0   
    #set colorscale_min -0.308149
    #set colorscale_max 0.00
    #mol scaleminmax 0 0 -0.308149 0.00
    #puts "Color scale set to $colorscale_min $colorscale_max"
}

# Dispersion Section
# Read the file containing atomic contributions
if {$file_extension eq "cube" && $has_atomwise_file} {
    set file_id [open $atomwise_file r]

    if {[catch {

        gets $file_id
        gets $file_id

        set lines [split [read $file_id] "\n"]

        close $file_id

        set atom_info {}
        set atom_id 0  ;

        foreach line $lines {
            if {[regexp {^\s*(\w+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)} $line - atom_name x y z edisp]} {
                lappend atom_info [list $atom_id $atom_name $x $y $z $edisp]
                incr atom_id
            }
        }

    } error]} {
        puts "An error occurred while reading the file: $error"
    }

    # Define the Tcl procedure to handle atom pick event
    proc write_edisp {args} { 

        global atom_info
        global vmd_pick_atom

        foreach info $atom_info {

            lassign $info atom_ID atomName x y z edisp

            if {$atom_ID == $vmd_pick_atom} {
                puts "You have selected atom $atomName (ID: $atom_ID)"
                puts "Coordinates:  x = $x,  y = $y,  z = $z"
                puts "Dispersion energy (kcal/mol): $edisp"
                break
            }
        }
    }

    # Associate the Tcl procedure with atom pick event
    trace add variable ::vmd_pick_event write write_edisp
}

# Rendering
display update ui
#render TachyonInternal ${file_name}.tga

# Convert from tga to png (ImageMagick® required)
#exec convert "${file_name}.tga" "${file_name}.png"
#exec rm "${file_name}.tga"

# Cleanup
# material delete Shiv
