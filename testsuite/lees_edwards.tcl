
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
#
# To generate the test system, use gen_lees_edwards.tcl.
#
source "tests_common.tcl"

require_feature "LEES_EDWARDS" on

puts "------------------------------------------"
puts "- Testcase lees_edwards.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------"


# system setup
set L                 10     ;#box size
set temperature        1.0   ;#temperature in reduced units
set T                500     ;#temperature in Kelvin

set k_temp [expr $T * 0.00198721]
set k_b    [expr $k_temp / $temperature]
set mass  1.67e-27        ; # to convert mass from AMU to kg
set fac   6.94e-21        ; # To convert energy from Kcal/mol to Joules

set time_step              0.001
set num_steps_equil 10000000
set skin                   0.2

set length_polymer  50
set N_chains         5
set n_part [expr $length_polymer * $N_chains]


#############################################################
# Energy Parameters                                         #
#############################################################

# Epsilon in KCal/mol divided by KbT in KCal/mol
set lj_epsilon [expr  0.112/$k_temp]    

#############################################################
#  Lennard Jones Parameters                                 #
#############################################################
set sigma 4.01
set lj1_sig 1.0
set lj1_cut [expr 2.5*$lj1_sig]
#set lj1_cut [expr 15.0-$skin]
set cut [expr 1/$lj1_cut]
set shift [expr 2* pow($cut, 6) - pow($cut, 12)]
set lj1_shift [expr $shift]

set harm_k [expr  18.0 / $k_temp]  
set harm_r0 3.0
set sigma   1.0

# set up global parameters
setmd box_l $L $L $L
setmd time_step $time_step
setmd skin $skin

# Seed the Random Number Generator
set ran_seed 1902
set cmd "t_random seed"
for {set i 0} {$i < [setmd n_nodes]} { incr i } { 
   lappend cmd [expr $ran_seed + $i] 
}
eval $cmd
expr srand($ran_seed)


#############################################################
#  Interaction Setup                                        #
#############################################################

inter 0 0 lennard-jones $lj_epsilon $lj1_sig $lj1_cut auto 0  

# Place particles randomly
for { set i 0 } { $i < $n_part } { incr i } {
    part $i pos [expr rand()*$L] [expr rand()*$L] [expr rand()*$L]
    part $i type 0
}

##################################################
# simulation

##start by shearing the system, and see what happens
set max_step_shear 1000000
set shear_equil       5000
set mean_from         1000
set shear_per            1
set write_per          100
set shear_rate           0.003
set offset               0.0

set f [open "thermalising.vtf" w]
writevsf $f

thermostat langevin $temperature 5.0
inter forcecap  10
puts "#Shear displacement   |  <vx^2>    <vy^2>   <vz^2>"
set mean_x2 0.0
set mean_y2 0.0
set mean_z2 0.0
set count   0.0
for { set step 0 } { $step < $shear_equil } { incr step $shear_per } {
    if { [expr $step % $write_per] == 0 } then {
        writevcf $f 
#folded
    }

    lees_edwards_offset  $offset
    integrate            $shear_per

    if { [expr $step % $write_per] == 0 } then {
        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set pos [ part $i print pos ] 
            set vel [ part $i print v ] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KE_horses_mouth [analyze energy kinetic]
        set KE_horses_mouth [expr 2.0 * $KE_horses_mouth / (3 * $n_part )]
        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]

        if { $step >= $mean_from } {
            set mean_x2 [expr $mean_x2 + $KEx ] 
            set mean_y2 [expr $mean_y2 + $KEy ] 
            set mean_z2 [expr $mean_z2 + $KEz ] 
            set count   [expr $count + 1.0 ] 
        }

        puts "$offset     $KEx $KEy $KEz"
        
    }


    set offset [expr $offset + $shear_rate * $shear_per]
}
close $f

set mean_x2 [expr $mean_x2 / $count ] 
set mean_y2 [expr $mean_y2 / $count ] 
set mean_z2 [expr $mean_z2 / $count ] 
puts "#means: $mean_x2 $mean_y2 $mean_z2"

##Below is the loop for an extended test, to
##see if the code is stable without a forcecap for a long time.

set f [open "noCap.vtf" w]
writevsf $f
inter forcecap  0
##thermostat langevin $temperature 0.01
for { set step 0 } { $step < $max_step_shear } { incr step $shear_per } {
    if { [expr $step % $write_per] == 0 } then {
        writevcf $f folded
    }
    lees_edwards_offset  $offset
    integrate            $shear_per

    if { [expr $step % $write_per] == 0 } then {
        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set pos [ part $i print pos ] 
            set vel [ part $i print v ] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KE_horses_mouth [analyze energy kinetic]
        set KE_horses_mouth [expr 2.0 * $KE_horses_mouth / (3 * $n_part )]
        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]

        puts "$offset     $KEx $KEy $KEz"
        
    }
    set offset [expr $offset + $shear_rate * $shear_per]
}





exit

