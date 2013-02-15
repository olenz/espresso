


####reduced units are:
## energy:   kcal/mol
## distance: Angstrom
## time:     

##################################################
# system setup
set file "harm_system.data.gz"

set L          30

###kB in kcal/mol/K is 0.00198721
set temperature 1.0   ;#temperature in reduced units
set T         500     ;#temperature in Kelvin
set k_temp [expr $T * 0.00198721]
set k_b    [expr $k_temp / $temperature]
set mass  1.67e-27        ; # to convert mass from AMU to kg
set fac   6.94e-21        ; # To convert energy from Kcal/mol to Joules

set time_step              0.001
set num_steps_equil 10000000
set skin                   0.2

set length_polymer  50
set N_chains        28
set n_part [expr $length_polymer * $N_chains]


#############################################################
# Energy Parameters                                         #
#############################################################

# Epsilon in KCal/mol divided by KbT in KCal/mol
set lj_epsilon [expr  0.112/$k_temp]    

#############################################################
#  Lennard Jones Parameter                                  #
#############################################################
set sigma 4.01
set lj1_sig 1.0
set lj1_cut [expr 2.5*$lj1_sig]
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

set therm_1  1000
set therm_2  7000
set therm_3 10000


#############################################################
# Units conversion                                          #
#############################################################
set fac 6.94e-21        ; # To convert energy from Kcal/mol to Joules
set length_conversion   [expr $sigma * 10e-10]
set time_conversion     [expr $length_conversion * $time_step * sqrt($mass/($T*$k_b*$fac))] 
set velocity_conversion [expr sqrt($fac * 2 * $mass)/$mass ]

puts "#time     conversion factor to   s: $time_conversion"
puts "#velocity conversion factor to m/s: $velocity_conversion"


#############################################################
#  Interaction Setup                                        #
#############################################################
# set up harm interaction
#inter 0 harm $harm_k $harm_r0

inter 0 0 lennard-jones $lj_epsilon $lj1_sig $lj1_cut auto 0  
#inter 0 1 lennard-jones $lj2_epsilon $lj1_sig $lj1_cut auto 0  
#inter 1 1 lennard-jones $lj3_epsilon $lj1_sig $lj1_cut auto 0 
  
# set up polymer
#polymer $N_chains $length_polymer $harm_r0 mode RW bond 0

# Particle types (0,1)
for { set i 0 } { $i < $n_part } { incr i } {
    part $i pos [expr rand()*$L] [expr rand()*$L] [expr rand()*$L]
    part $i type 0
}

##################################################
# simulation

# run a number of integration steps
puts "#thermalising:"

set stepCount 0
set step    100
puts "##step || 2x(mean KE) || T?"

#lees_edwards_offset  0.0
#puts [lees_edwards_offset  7.5]

##start by shearing the system, and see what happens
set max_step_shear 1000000
set shear_equil      10000
set shear_per            1
set write_per          100
set therm_per           10
set shear_rate           0.01
set offset               0.0

set f [open "thermalising.vtf" w]
writevsf $f


thermostat langevin $temperature 5.0
inter forcecap  10
integrate 100

for { set step 0 } { $step < $shear_equil } { incr step $shear_per } {
    if { [expr $step % $write_per] == 0 } then {
        writevcf $f folded
    }

    lees_edwards_offset  $offset
    integrate            $shear_per

    if { [expr $step % $write_per] == 0 } then {
        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        set mX  0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set pos [ part $i print pos ] 
            set vel [ part $i print v ] 
            set mX  [expr $mX + [expr [lindex $vel 0]*([lindex $pos 1]-0.5*$L)]] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KE_horses_mouth [analyze energy kinetic]
        set KE_horses_mouth [expr 2.0 * $KE_horses_mouth / (3 * $n_part )]
        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]
        set mX              [expr $mX / $n_part ]

        puts "offset: $offset temp: $KEx $KEy $KEz $mX"
        
    }
    set offset [expr $offset + $shear_rate * $shear_per]
    invalidate_system
}
close $f
set f [open "noCap.vtf" w]
writevsf $f
inter ljforcecap  0
set write_per    10
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
        set mX  0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set pos [ part $i print pos ] 
            set vel [ part $i print v ] 
            set mX  [expr $mX + [expr [lindex $vel 0]*([lindex $pos 1]-0.5*$L)]] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KE_horses_mouth [analyze energy kinetic]
        set KE_horses_mouth [expr 2.0 * $KE_horses_mouth / (3 * $n_part )]
        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]
        set mX              [expr $mX / $n_part ]

        puts "offset: $offset temp: $KEx $KEy $KEz $mX"
        
    }
    set offset [expr $offset + $shear_rate * $shear_per]
    invalidate_system
}





exit

