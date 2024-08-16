# VMD Script for Publication-Quality Images
# Author: Gianluca Regni
# Copyright (c) 2024 Gianluca Regni
# License: LGPL-3.0
# Credits: Please cite the author(s) and the references if you use this program.

# References:
# VMD manual: https://www.ks.uiuc.edu/Research/vmd/current/ug.pdf
# Shiv material settings adapted from https://shivupa.github.io/blog/making-publication-quality-images-with-vmd/
# Color scale bar adapted from https://github.com/thatchristoph/vmd-cvs-github/blob/master/plugins/colorscalebar/colorscalebar.tcl

# Description:
# This script is designed to visualize LD density (difference) function.
# It adjusts VMD settings to generate publication-quality images.
# Additionally,  it includes an improved pick event feature that displays atomic contributions to dispersion
# energy in the terminal, extracted from the _*atomwise.txt_ file.
# It is designed to be used in combination with other programs of LDD suite (https://github.com/bistonigroup/LDDsuite).

# Usage:
# vmd -e PubQualityVMD.tcl 

# ------------------------------------
# Variables Declarations
# ----------------------------------

set file_name "water.xyz"; # file name
set autoscale 0; # use 1 to disable auto color scale
set colorscale_min 0.00; # min value for manual color scale
set colorscale_max 0.00; # max value for manual color scale
set scale_val 1.3; # increase molecular size 
set rot_y 0;  # rotation on y axe
set rot_x 0;  # rotation on x axe
set rot_z 0;  # rotation on z axe
set render_mode 0; # use 0 to disable rendering

# ------------------------------------
# File Names Assignation
# ------------------------------------

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


# Color Scale bar

namespace eval ::ColorScaleBar:: {
  namespace export color_scale_bar delete_color_scale_bar
  variable version      1.4; # version number of this package.
  variable w
  variable molid        "0"; # id of the molecule
  variable moltxt        ""; # title the molecule
  variable repid        "0"; # id of the representation
  variable reptxt        ""; # title/name of the rep
  variable bar_mol       -1
  variable bar_text "Color Scale:"; # legend text for the colorbar
  variable text_show    "0"; # whether to show the legend.
  # variable bar_orient   "0"; # orientation of the colorbar (vertical=0, horizontal=1)
  # variable bar_text_pos "s"; # position of the colorbar (n,s,w,e)
  # variable text_orient  "0"; # orientation of the text label to the bar.
  variable lengthsetting  0.8
  variable widthsetting   0.05
  variable autoscale      0
  variable fixedsetting   1
  # variable minvalue     0
  # variable maxvalue     100
  variable axislabels     5
  variable textcolor      white
  variable fpformat       1      # 0=decimal, 1=scientific
}

package provide colorscalebar $ColorScaleBar::version

proc ::ColorScaleBar::color_scale_bar {{length 0.3} {width 0.05} {auto_scale $autoscale} {fixed 1} {label_num 1} {text 16} {fp_format 0} {x_pos 1.3} {y_pos 0.0} {replacebar 1} {molid top} {repid 0} {showlegend 0} {legend "Dispersion Color Scale:"}} {
  variable bar_mol

  set m $molid
  if {[string equal "$molid" "top"]} {set m [molinfo top]}

  # if there's already a color scale bar molecule, delete and replace it
  #if {$replacebar == 1} {
  #  delete_color_scale_bar
  #}

  # So that the draw cmds will work right, must save top mol and set top
  # to our new created mol, then set it back later.
  set old_top [molinfo top]
  if { $old_top == -1 } {
    vmdcon -error "Color Scale Bar Plugin: No molecules loaded"
    return -1;
  }

  array set viewpoints {}
  foreach mol [molinfo list] {
  # save orientation and zoom parameters for each molecule.
    set viewpoints($mol) [molinfo $mol get { 
    center_matrix rotate_matrix scale_matrix global_matrix}]
  } 
 
  # don't update the display while we do this since otherwise there
  # will be thousands of individual draw commands going on 
  display update off
  display resetview

  # XXX: the previous heuristics for setting the min/max values only worked, if there
  # was exactly one molecule using a color scale and the colorized rep was the first
  # representation. so quite often it could not work. color scales are per rep (and mol),
  # so we have to pick _one_representation.
  if {$auto_scale == 1} {
    set min  999999
    set max -999999
    set minmax {0.0 0.0}
    if {[catch {mol scaleminmax $m $repid} minmax]} {
      #XXX: print error message, if in text mode.
      set min 0
      set max 0
    } else {
      set min [lindex $minmax 0]
      set max [lindex $minmax 1]
    }
  } else {
    catch {mol scaleminmax $m $repid} minmax
    # Parse the minmax list to get the min and max coordinates
    set min [lindex $minmax 0]
    set max [lindex $minmax 1]
    # Output the min and max coordinates
    puts "Min: $min"
    puts "Max: $max"
  }

  # Convert Hartree to kcal/mol
  #set min [expr $min * 627.503]
  #set max [expr $max * 627.503]
  #puts "... Min and Max of color scale converted from Eh to kcal/mol"

  # check for situation where all mols were skipped by the catch statement
  #if { $min > $max } {
  #  set min 0
  #  set max 0
  #}

  # Create a seperate molid to draw in, so it's possible for the user to 
  # delete the bar.
  set bar_mol [mol new]
  mol top $bar_mol
  mol rename top "Color Scale Bar"

  # If a fixed bar was requested...
  if {$fixed == 1} {
    mol fix $bar_mol
  }

  # set position relative to top molecule 
  # We want to draw relative to the location of the top mol so that the bar 
  # will always show up nicely.
  #set center [molinfo $old_top get center]
  #set center [regsub -all {[{}]} $center ""]
  #set center [split $center]
  #set start_y [expr [lindex $center 1] - (0.5 * $length)]
  #set use_x [expr 1+[lindex $center 0]]
  #set use_z [lindex $center 2]

  # set in absolute screen position
  set start_y [expr (-0.5 * $length) + $y_pos]
  set use_x $x_pos
  set use_z 0

  # draw background border behind bar, same color as text
  draw color $text

  # disable material properties for the color scale bar background
  # so that it looks truly black (no specular) when it's set black
  set bw [expr $width * 0.05]
  set lx [expr $use_x             - $bw]
  set rx [expr $use_x   + $width  + $bw] 
  set ly [expr $start_y           - $bw]
  set uy [expr $start_y + $length + $bw]
  set bz [expr $use_z - 0.00001]
  
  draw line "$lx $ly $bz" "$lx $uy $bz" width 2
  draw line "$lx $uy $bz" "$rx $uy $bz" width 2
  draw line "$rx $uy $bz" "$rx $ly $bz" width 2
  draw line "$rx $ly $bz" "$lx $ly $bz" width 2

  # draw the color bar
  set mincolorid [colorinfo num] 
  set maxcolorid [expr [colorinfo max] - 1]
  set numscaleids [expr $maxcolorid - $mincolorid]
  set step [expr $length / double($numscaleids)]
  for {set colorid $mincolorid } { $colorid <= $maxcolorid } {incr colorid 1 } {
    draw color $colorid
    set cur_y [ expr $start_y + ($colorid - $mincolorid) * $step ]
    draw line "$use_x $cur_y $use_z"  "[expr $use_x+$width] $cur_y $use_z"
  }

  # draw the labels
  set coord_x [expr (1.2*$width)+$use_x];
  set step_size [expr $length / $label_num]
  set color_step [expr double($numscaleids)/$label_num]
  set value_step [expr ($max - $min ) / double ($label_num)]
  for {set i 0} {$i <= $label_num } { incr i 1} {
    draw color $text
    set coord_y [expr $start_y+$i * $step_size ]
    set cur_text [expr $min + $i * $value_step ]

    set labeltxt ""
    if { $fp_format == 0 } {
      # format the string in decimal notation
      # we save a spot for a leading '-' sign
      set labeltxt [format "% 6.2f"  $cur_text]
    } else {
      # format the string in scientific notation
      # we save a spot for a leading '-' sign
      # since there are only 1024 distinct colors, there's no point in 
      # displaying more than 3 decimal places after the decimal point
      set labeltxt [format "% #.3e"  $cur_text]
    }
    draw text  "$coord_x $coord_y $use_z" "$labeltxt"
    draw line "[expr $use_x+$width] $coord_y $use_z" "[expr $use_x+(1.45*$width)] $coord_y $use_z" width 2
  }

  if {$showlegend == 1} {
      # set in absolute screen position
      draw color 6
      set use_y [expr {$start_y + $length + 0.15}]
      draw text  "$x_pos $use_y $use_z" "$legend"
  }
  # re-set top
  if { $old_top != -1 } {
    mol top $old_top
  }

  foreach mol [molinfo list] {
    if {$mol == $bar_mol} continue
    # restore orientation and zoom
    molinfo $mol set {center_matrix rotate_matrix scale_matrix
    global_matrix} $viewpoints($mol)
  }

  display update on

  return 0
}

# ------------------------------------
# Main Script
# ------------------------------------

# Load Molecule
mol new "${file_name}" type ${file_extension}

# Display Settings
display rendermode GLSL
display projection Orthographic
display ambientocclusion on
display shadows on
if {$file_extension eq "cube"} {
    display depthcue off
}
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

# Molecular Representation
color Element C gray

mol delrep 0 top
mol addrep 0
mol modcolor 0 0 Element
mol modstyle 0 0 CPK 0.8 0.3 72 72
if {$file_extension eq "cube"} {
    mol modcolor 0 0 Volume 0
    if {$autoscale == 1} {
      mol scaleminmax 0 0 $colorscale_min $colorscale_max
      puts "Color scale set to $colorscale_min $colorscale_max"
    }
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

    # Associate the Tcl procedure with atom pick event
    trace add variable ::vmd_pick_event write write_edisp
}

# Draw color scale bar
if {$file_extension eq "cube"} {
namespace import ::ColorScaleBar::*
color_scale_bar
}

# Rotation and Scaling
rotate y by $rot_y
rotate x by $rot_x
rotate z by $rot_z
scale by $scale_val

# Cleanup
# material delete Shiv

# Rendering
display update ui

if {$render_mode} {
  render TachyonInternal ${file_name}.tga
  # Convert from tga to png (ImageMagick required)
  exec convert "${file_name}.tga" "${file_name}.png"
  exec rm "${file_name}.tga"
}
