#!/bin/bash

############################################################################################
##                                  evaporate.sh                                          ##
############################################################################################
#                                                                                          #
# Citation:                                                                                #
#  - R. Alessandri, J. J. Uusitalo, A. H. de Vries, R. W. A. Havenith, S. J. Marrink,      #
#    "Bulk Heterojunction Morphologies with Atomistic Resolution from Coarse-Grain         #
#     Solvent Evaporation Simulations", JACS 2017, 139, 3697 [DOI: 10.1021/jacs.6b11717]   #
#                                                                                          #
############################################################################################

# USE:      ./evaporation.sh      <number_mol1> <number_mol2> [start]
# USE:      sbatch evaporation.sh <number_mol1> <number_mol2> [start]  (SLURM)
#  
# > <number_mol1> = number of mol1 molecules you want in the system
# > <number_mol2> = number of mol2 molecules you want in the system
# > [start]       = the *optional* third argument "start" tells the script whether you are
#                   starting an evaporation; if no third argument is given to the script,
#                   the evaporation will be a restart of a unfinished evaporation 

# List of files that need to be in this folder:
# COORDINATES    :  files indicated as 'mol1_geom_file', 'mol2_geom_file', and 'solv_box'
#                   (e.g., P3HT_n48.gro, PCBM.gro, box_CB.gro)
# TOPOLOGY       :  SYSTEM_empty.top and itp files needed
#                   (e.g., PCBM.itp, P3HT_n48.itp, martini_v2.2_CNP.itp, martini_v2.0_solvents.itp) 
# RUN PARAMETERS :  martini_v2.x_new_minimize.mdp, martini_v2.x_new_eq_NVT.mdp, 
#                   martini_v2.x_new_eq_NPT_semiiso.mdp, martini_v2.x_new_run_semiiso.mdp
 
############################################################################################


#################################
##        Load GROMACS 5.x     ##
#################################
## module load gromacs/5.x
## source gromacs5

# Set flags for mdrun (-cpi is set by the script when needed, no need to include it here)
FLAGS='-dlb yes -rdd 1.4'
# Set working directory
echo "The job is started here:"
pwd


#################################
## System-dependent parameters ##
#################################
# The following four lines will have to be adapted depending on the molecules in the system; names are self-explanatory
beads_per_polym_molecule=288
beads_per_fulle_molecule=21
polym_geom_file=P3HT_n48.gro
fulle_geom_file=PCBM.gro
# Simulation box dimensions
box_X=30
box_Y=30
box_Z=88


#################
## Check INPUT ##
#################
# Check if the number of polymer and fullerene molecules has been given
polym_mol="$1"
fulle_mol="$2"
polym_beads=$(expr $polym_mol \* $beads_per_polym_molecule)
fulle_beads=$(expr $fulle_mol \* $beads_per_fulle_molecule)
size_1=${#polym_mol}
size_2=${#polym_mol}
if [ $size_1 -eq 0 -o $size_2 -eq 0 ] ; then
   echo ""
   echo "Missing number of polymer and/or fullerene molecules. Check that"
   exit
fi



##########################################
## START a new evaporation OR RESTART ? ##
##########################################
is_a_start="$3"
size_3=${#is_a_start}

if [ "$is_a_start" == "start" ] ; then
####################
## It's a start!  ##
####################
   mkdir step000
   echo "Moving to folder step000"
   cd    step000
   # Generate a box of P3HT:PCBM:CB
   echo "Generate the initial box"
   srun -n  1 gmx_mpi insert-molecules -ci ../$polym_geom_file -nmol $polym_mol -box $box_X $box_Y $box_Z  -radius 0.27 # Martini beads radius
   wait
   srun -n  1 gmx_mpi insert-molecules -ci ../$fulle_geom_file -f out.gro -nmol $fulle_mol                 -radius 0.27 # Martini beads radius
   wait
   srun -n  1 gmx_mpi solvate          -cp out.gro -cs ../box_CB.gro                                       -radius 0.18 # leads to a good starting density for CB
   wait
   # Count molecules of CB that are in the box and add this to the TOP file
   cb_mol=$(grep "CB      Cl" out.gro | wc -l)
   cp ../SYSTEM_empty.top SYSTEM_step000.top
   echo "P3HT           $polym_mol" >> SYSTEM_step000.top
   echo "PCBM           $fulle_mol" >> SYSTEM_step000.top
   echo "CB             $cb_mol" >> SYSTEM_step000.top
   # Minimisation, equilibration and first run
   echo "Minimisation (x4), eq (NVT, NPT), and first run"
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_minimize.mdp          -p SYSTEM_step000.top -c out.gro          -o min_step000.tpr  -po min_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -v -deffnm min_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_minimize.mdp          -p SYSTEM_step000.top -c min_step000.gro  -o min2_step000.tpr -po min2_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -v -deffnm min2_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_minimize.mdp          -p SYSTEM_step000.top -c min2_step000.gro -o min3_step000.tpr -po min3_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -v -deffnm min3_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_minimize.mdp          -p SYSTEM_step000.top -c min3_step000.gro -o min4_step000.tpr -po min4_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -v -deffnm min4_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_eq_NVT.mdp           -p SYSTEM_step000.top -c min4_step000.gro -o NVT_step000.tpr  -po NVT_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -deffnm NVT_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_eq_NPT_semiiso.mdp   -p SYSTEM_step000.top -c NVT_step000.gro  -o NPT_step000.tpr  -po NPT_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -deffnm NPT_step000 >> mdrun.log 2>&1
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_run_semiiso.mdp       -p SYSTEM_step000.top -c NPT_step000.gro  -o run_step000.tpr  -po run_step000.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -deffnm run_step000 >> mdrun.log 2>&1
   #
   touch finished
   # Define where to start the evaporation loop from (i.e., step000, since it's a start)
   starting_i=000
   checkpoint_restart="no"

elif [ $size_3 = 0 ] ; then
######################
## It's a restart!  ##
######################
   echo -e "\n That's an evaporation restart"
   # Where do we have to start from? Look for folders and the file 'finished'
   for i in {000..700}  # MAX NUMBER OF EVAPORATION STEPS - change if needed
   do
   find ./step$i -name 'finished' 1> /dev/null
   findres=$?
   if [ $findres != 0 ] ; then
   # No folder has been found for this step. The previous step had completed.
      echo "Folder not found; we have to restart from" $i
      # Define where to start the evaporation loop from
      previous=$(expr $i - 1)               # previous step variable
      previous=$(printf "%03d" $previous)   # makes sure also "previous" is a 3 digits long variable 
      echo "Let's move to the previous folder (number" $previous ")"
      cd step$previous
      starting_i=$previous
      checkpoint_restart="no"
      break
   fi
   found_a_finished=$(find ./step$i -name 'finished' | wc -l)
   if [ $found_a_finished = 0 ] ; then
   # No file 'finished' has been found in the current folder --> start from CHECKPOINT (state_prev.cpt)
      echo "we have to restart from step" $i
      # Define where to start the evaporation loop from
      echo "Let's move to the folder (number" $i ")"
      cd step$i
      # Let's treat "i" as the previous step, for compatibility with EVAPORATION
      previous=$(expr $i - 1)               # previous step variable
      previous=$(printf "%03d" $previous)   # makes sure also "previous" is a 3 digits long variable 
      starting_i=$previous
      checkpoint_restart="yes"
      break
   elif [ $found_a_finished = 1 ] ; then
      echo "step" $i "had finished"
   else
      echo "Something is wrong? More than one 'finished' file has been found in folder 'step"$i"'"
      exit
   fi
   done

else
###############
## Error !   ##
###############
   echo "If you want to start a new evaporation process pass "start" (as third argument) to the script."
   echo "If you want to restart a unfinished evaporation process pass no arguments after the number of polymer and fullerene molecules."
   exit
fi



##########################################
##            EVAPORATION               ##
##########################################
i=$starting_i
nmax=700
#----------------------------------------------------------------------#
while [ "$i" -le "$nmax" ]
do

j=$(expr $i + 1)        # "j" = 'next step'
j=$(printf "%03d" $j)   # makes sure also "j" is a 3 digits long variable                     


if [ "$checkpoint_restart" != "yes" ] ; then
######################################################
## 1a. Setting up a new step: checks                ##
######################################################
   cb_mol=$(grep "CB      Cl" run_step$i.gro | wc -l) # how many solvent molecules
   # Check grep: if grep not succeeded, exit
   grep -q "CB      Cl" run_step$i.gro
   grepres=$? # if=0 -> OK, if=1 -> nothing matched "CB      Cl", if=2 -> file not found
   if [ $grepres != 0 ] ; then
      echo "Something went wrong."
      echo "Either the file has not been found or there are no CB molecules in the previous-step GRO file."
      exit
   fi
   # Print the number of lines in the current GRO and the last line number (for checking)
   n_lines=$(( 2 + $polym_beads + $fulle_beads + 3 * $cb_mol))
   echo "number of lines in the current GRO" $n_lines "+1"
   n_last_line=$(expr $n_lines + 1)
   echo "the last line is line number" $n_last_line 

#####################################################
## 1b. Setting up a new step: prepare input files  ##
#####################################################
   mkdir ../step$j

# Compute how many solvent molecules have to be removed and generate next GRO without those molecules
   evaporated=$(expr $cb_mol / 80) # let's remove 1.25% of the solvent
   if [ $evaporated -lt 10 ] ; then
      evaporated=10
   fi
   echo $evaporated "CB molecules have evaporated"
   n_lines=$(( $n_lines - 3 * $evaporated))
   echo "number of lines in the new GRO" $n_lines "+1"
   echo "run_step$i.gro"
   sed -n '1,'$n_lines' p' run_step$i.gro > ../step$j/start_step$j.gro
   # Has the evaporation finished?
   cb_mol=$(expr $cb_mol - $evaporated)
   if [ $cb_mol -lt 1 ] ; then
      echo "evaporation has finished (step number" $i ")"
      exit
   fi
   # Change the number of atoms in the (second line of the) new GRO file accordingly
   n_atoms=$(( $polym_beads + $fulle_beads + 3 * $cb_mol))
   sed -i '2 c\  '$n_atoms' ' ../step$j/start_step$j.gro
   # Add the box size to the new GRO file
   sed -n ''$n_last_line' p' run_step$i.gro >> ../step$j/start_step$j.gro

   # Move to the new step folder
   cd ../step$j
   # Generate the new TOP file
   cp ../SYSTEM_empty.top SYSTEM_step$j.top
   echo "P3HT           $polym_mol" >> SYSTEM_step$j.top
   echo "PCBM           $fulle_mol" >> SYSTEM_step$j.top
   echo "CB             $cb_mol"    >> SYSTEM_step$j.top

   ### NPT equilibration before run ###
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_eq_NPT_semiiso.mdp    -p SYSTEM_step$j.top -c start_step$j.gro -o NPT_step$j.tpr  -po NPT_step$j.mdp -maxwarn 10
   wait
   srun gmx_mpi mdrun $FLAGS -deffnm NPT_step$j >> mdrun.log 2>&1
   ### 
   srun -n 1 gmx_mpi grompp -f ../martini_v2.x_new_run_semiiso.mdp       -p SYSTEM_step$j.top -c NPT_step$j.gro -o run_step$j.tpr -po run_step$j.mdp -maxwarn 10
   wait


######################################
## 2. Run the new evaporation step  ##
######################################
   echo "Run step number" $j

else
   echo "Run step number" $j "from checkpoint"
fi

srun gmx_mpi mdrun $FLAGS -cpi run_step${j}_prev.cpt -deffnm run_step$j >> mdrun.log 2>&1

# remove trajectory files if step != step000 (optional)
#if [ $j != 0 ]; then
#   rm *_step$j.xtc  *_step$j.trr
#fi
#rm \#*

touch finished
checkpoint_restart="no"


i=$(expr $i + 1)
i=$(printf "%03d" $i)   # makes sure the new "i" is a 3 digits long variable 
#----------------------------------------------------------------------#
done


echo "Evaporation not finished yet! Increase number of steps."
