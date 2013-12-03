# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
# 
source "tests_common.tcl"

puts "------------------------------------------------"
puts "- Testcase thermostat_z_only.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

puts [code_info]

require_feature THERMOSTAT_RESCALING
require_feature THERMOSTAT_Z_ONLY
require_feature LEES_EDWARDS


set tcl_precision 5
set n_part      100

proc print_ke_by_axis { n_part } {

        

    if { [regexp "ROTATION" [code_info]] } { 
        set deg_free 6
    } else {
        set deg_free 3
    }


        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set vel [ part $i print v ] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }
        set KEx [expr 0.5 * $vx2 / ( $n_part )]
        set KEy [expr 0.5 * $vy2 / ( $n_part )]
        set KEz [expr 0.5 * $vz2 / ( $n_part )]

        set toteng [analyze energy total]
        set cureng [analyze energy kin] 

        ##equipartition: mean energy = 0.5 * temperature * deg_free
        set curtemp [expr (2.0*$cureng/[setmd n_part]) / $deg_free] 
        puts "$toteng $cureng $curtemp $KEx $KEy $KEz"
}

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos v omega} 
    close $f
}


    if { [regexp "ROTATION" [code_info]] } {
        puts "##rotation found, 6 degrees of freedom"
        puts "##...the rotational degrees of freedom are not treated directly."
        set deg_free 6
    } else {
        puts "##rotation not active, therefore 3 degrees of freedom"
        set deg_free 3
    }

    set filename "thermostat.data"

    read_data $filename

    # seed the RNG
    set cmd "t_random seed"
    for {set i 1} {$i <= [setmd n_nodes]} { incr i } {
#	lappend cmd [expr [pid] + $i] 
        lappend cmd $i 
    }
    eval $cmd


    # generate some particles
    set box_l 10
    setmd box_l $box_l $box_l $box_l
    for {set i 0} {$i < $n_part} {incr i} {
    	part $i pos [expr rand()*$box_l] [expr rand()*$box_l] [expr rand()*$box_l] 
    	part $i v [expr 1] [expr 1] [expr 1] 
    }


    ##Turn on the thermostat
    thermostat langevin 1.0 1.0
    setmd time_step 0.01
    setmd skin 0.5

    puts "##Initial heating: Simple Langevin"
    puts "##Expected temperature is 1.0"
    for {set t 0} { $t < 10} { incr t 1 } {

        print_ke_by_axis $n_part
        integrate 1

    }
   
    # Try with some Lees-Edwards shear
    puts "##Shearing in the x-y plane"
    for {set i 0} { $i < 10} { incr i } {

        set offset [expr $i * 0.05]
        lees_edwards_offset $offset
 
        print_ke_by_axis $n_part
        puts "#LE offset: $offset"

	integrate 1
    }

    thermostat z_only on
    puts "##Langevin Z only"
    for {set t 0} { $t < 10} { incr t 1 } {

        set offset [expr $i * 0.1]
        lees_edwards_offset $offset

        print_ke_by_axis $n_part
        puts "#LE offset: $offset"

        integrate 1
    }


    # switch on the rescaling thermostat (turns off Langevin)
    thermostat velocity_rescaling 1.0 0.05
    thermostat z_only on
    puts "##Rescaling thermostat: Z only"
    for {set i 0} { $i < 10} { incr i } {

        set offset [expr $i * 0.05]
        lees_edwards_offset $offset

        print_ke_by_axis $n_part
        puts "#LE offset: $offset"
 
	integrate 1
    }

exit 0
