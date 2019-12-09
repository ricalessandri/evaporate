#!/bin/bash

############################################################################################
##                                  evaporate.sh                                          ##
############################################################################################
                                                                                            
# CITATION:                                                                                 
#  - R. Alessandri, J. J. Uusitalo, A. H. de Vries, R. W. A. Havenith, S. J. Marrink,       
#    "Bulk Heterojunction Morphologies with Atomistic Resolution from Coarse-Grain          
#     Solvent Evaporation Simulations", JACS 2017, 139, 3697 [DOI: 10.1021/jacs.6b11717]    
#                                                                                           

# USE:      ./evaporation.sh      <number_mol1> <number_mol2> [start]
# USE:      sbatch evaporation.sh <number_mol1> <number_mol2> [start]  (SLURM)
#  
# > <number_mol1> = number of mol1 molecules you want in the system
# > <number_mol2> = number of mol2 molecules you want in the system
# > [start]       = the *optional* third argument "start" tells the script whether you are
#                   starting an evaporation; if no third argument is given to the script,
#                   the evaporation will be a restart of a unfinished evaporation 
#
# List of FILES that need to be in this folder:
# COORDINATES    :  files indicated as 'mol1_geom_file', 'mol2_geom_file', and 'solv_box'
#                   (e.g., P3HT_n48.gro, PCBM.gro, box_CB.gro)
# TOPOLOGY       :  SYSTEM_empty.top and itp files needed
#                   (e.g., PCBM.itp, P3HT_n48.itp, martini_v2.2_CNP.itp, martini_v2.0_solvents.itp) 
# RUN PARAMETERS :  martini_v2.x_new_minimize.mdp, martini_v2.x_new_eq_NVT.mdp, 
#                   martini_v2.x_new_eq_NPT_semiiso.mdp, martini_v2.x_new_run_semiiso.mdp
 


############################################################################################
#                             GROMACS: Load and set GROMACS                                # 
############################################################################################

# QUEUEING system:
#SBATCH ...
#SBATCH ... 
 
# module load gromacs/5.x (or later)
# source gromacs5 (or later)
 
# SET GMX COMMANDS aliases
INSMOL="srun -n  1 gmx_mpi insert-molecules" # Cartesius
SOLVATE="srun -n  1 gmx_mpi solvate"         # Cartesius
GROMPP="srun -n 1 gmx_mpi grompp -normvsbds" # Cartesius ("-normvsbds" needed because of https://www.mail-archive.com/gromacs.org_gmx-users@maillist.sys.kth.se/msg35762.html) 
MDRUN="srun gmx_mpi mdrun"                   # Cartesius



############################################################################################
#                              System-dependent parameters                                 # 
############################################################################################

# The following lines will have to be adapted depending on the molecules in the system. 
# Names are mostly self-explanatory.
 
# Solute molecules
beads_per_mol1_molecule=288 # P3HT (48-mer)
beads_per_mol2_molecule=21  # PCBM
mol1_geom_file=P3HT_n48.gro
mol2_geom_file=PCBM.gro 
mol1_name=P3HT
mol2_name=PCBM
 
# Solvent
## solv_box=box_CLF.gro
## solv_name=CLF
## no_solv_beads=1
## solv_pattern="CLF     CX"  # grepping CLF molecules
## solv_density_param=27      # good starting density for CLF (?)
solv_box=box_CB.gro
solv_name=CB
no_solv_beads=3
solv_pattern="CB      Cl"  # grepping CLBZ molecules
solv_density_param=18      # good starting density for CB (Martini 2.2)
 
# (Initial) simulation box dimensions
box_X=30
box_Y=30
box_Z=88 # to get an initial box dimension of ~80 with CB (Martini 2.2) 
#box_Z=80 # CLF
 
# Set flags for mdrun 
FLAGS='-dlb yes -rdd 1.4'
NMAX=700

############################################################################################
############################################################################################


# TODO:
# - remove hardcoded 700 (must be = to NMAX)
# - control printing of XTC files via a flag (now, default is not to produce them;
#   could be just controlled via the mdps (using a SED command or so, for example))




#################
## Check INPUT ##
#################
# Check if the number of molecule_1 and molecule_2 molecules has been given
number_mol1="$1"
number_mol2="$2"
mol1_beads=$(expr $number_mol1 \* $beads_per_mol1_molecule)
mol2_beads=$(expr $number_mol2 \* $beads_per_mol2_molecule)
size_1=${#number_mol1}
size_2=${#number_mol2}
if [ $size_1 -eq 0 -o $size_2 -eq 0 ] ; then
   echo ""
   echo "Missing number of molecule_1 and/or molecule_2 molecules. Check that"
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
   # Generate a box of mol1:mol2:solvent
   echo "Generate the initial box"
   $INSMOL -ci ../$mol1_geom_file -nmol $number_mol1 -box $box_X $box_Y $box_Z  -radius 0.27 # Martini beads radius
   wait
   $INSMOL -ci ../$mol2_geom_file -f out.gro -nmol $number_mol2                 -radius 0.27 # Martini beads radius
   wait
   $SOLVATE          -cp out.gro -cs ../$solv_box                            -radius 0.${solv_density_param}
   wait
   # Count molecules of solvent that are in the box and add this to the TOP file
   no_solv_mol=$(grep "${solv_pattern}" out.gro | wc -l)
   cp ../SYSTEM_empty.top SYSTEM_step000.top
   echo "${mol1_name}           $number_mol1" >> SYSTEM_step000.top
   echo "${mol2_name}           $number_mol2" >> SYSTEM_step000.top
   echo "${solv_name}           $no_solv_mol" >> SYSTEM_step000.top
   # Minimisation, equilibration and first run
   echo "Minimisation (x3), eq (NVT, NPT), and first run"
   $GROMPP -f ../martini_v2.x_new_minimize.mdp  -p SYSTEM_step000.top -c out.gro          -o min_step000.tpr  -po min_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -v -deffnm min_step000 >> mdrun.log 2>&1
   wait
   $GROMPP -f ../martini_v2.x_new_minimize.mdp  -p SYSTEM_step000.top -c min_step000.gro  -o min2_step000.tpr -po min2_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -v -deffnm min2_step000 >> mdrun.log 2>&1
   wait
   $GROMPP -f ../martini_v2.x_new_minimize.mdp  -p SYSTEM_step000.top -c min2_step000.gro -o min3_step000.tpr -po min3_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -v -deffnm min3_step000 >> mdrun.log 2>&1
   wait
   $GROMPP -f ../martini_v2.x_new_eq_NVT.mdp    -p SYSTEM_step000.top -c min3_step000.gro -o NVT_step000.tpr  -po NVT_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -deffnm NVT_step000 >> mdrun.log 2>&1
   $GROMPP -f ../martini_v2.x_new_eq_NPT_semiiso.mdp  -p SYSTEM_step000.top -c NVT_step000.gro  -o NPT_step000.tpr  -po NPT_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -deffnm NPT_step000 >> mdrun.log 2>&1
   $GROMPP -f ../martini_v2.x_new_run_semiiso.mdp  -p SYSTEM_step000.top -c NPT_step000.gro  -o run_step000.tpr  -po run_step000.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -deffnm run_step000 >> mdrun.log 2>&1
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
   echo "If you want to restart a unfinished evaporation process pass no arguments after the number of molecule_1 and molecule_2 molecules."
   exit
fi



#####################
##   EVAPORATION   ##
#####################
i=$starting_i
#----------------------------------------------------------------------#
while [ "$i" -le "$NMAX" ]
do

j=$(expr $i + 1)        # "j" = 'next step'
j=$(printf "%03d" $j)   # makes sure also "j" is a 3 digits long variable                     


if [ "$checkpoint_restart" != "yes" ] ; then
######################################################
## 1a. Setting up a new step: checks                ##
######################################################
   no_solv_mol=$(grep "${solv_pattern}" run_step$i.gro | wc -l) # how many solvent molecules
   # Check grep: if grep not succeeded, exit
   grep -q "${solv_pattern}" run_step$i.gro
   grepres=$? # if=0 -> OK, if=1 -> nothing matched "${solv_pattern}", if=2 -> file not found
   if [ $grepres != 0 ] ; then
      echo "Something went wrong."
      echo "Either the file has not been found or there are no solvent molecules in the previous-step GRO file."
      exit
   fi
   # Print the number of lines in the current GRO and the last line number (for checking)
   n_lines=$(( 2 + $mol1_beads + $mol2_beads + $no_solv_beads * $no_solv_mol))
   echo "number of lines in the current GRO" $n_lines "+1"
   n_last_line=$(expr $n_lines + 1)
   echo "the last line is line number" $n_last_line 

#####################################################
## 1b. Setting up a new step: prepare input files  ##
#####################################################
   # Compute how many solvent molecules have to be removed and generate next GRO without those molecules
   evaporated=$(expr $no_solv_mol / 80) # let's remove 1.25% of the solvent
   if [ $evaporated -lt 10 ] ; then
      evaporated=10
      if [ $no_solv_mol -lt 10 ] ; then
         evaporated=${no_solv_mol}
         # Hold on, has the evaporation finished?
         if [ $no_solv_mol -lt 1 ] ; then
            echo "evaporation has finished (step number" $i ")"
            touch evaporation_is_DONE
            exit
         fi
      fi
   fi

   # Set up the next step
   mkdir ../step$j
   no_solv_mol=$(expr $no_solv_mol - $evaporated)
   echo $evaporated " solvent molecules have evaporated"
   n_lines=$(( $n_lines - $no_solv_beads * $evaporated))
   echo "number of lines in the new GRO" $n_lines "+1"
   echo "run_step$i.gro"
   sed -n '1,'$n_lines' p' run_step$i.gro > ../step$j/start_step$j.gro
   # Change the number of atoms in the (second line of the) new GRO file accordingly
   n_atoms=$(( $mol1_beads + $mol2_beads +$no_solv_beads* $no_solv_mol))
   sed -i '2 c\  '$n_atoms' ' ../step$j/start_step$j.gro
   # Add the box size to the new GRO file
   sed -n ''$n_last_line' p' run_step$i.gro >> ../step$j/start_step$j.gro

   # Move to the new step folder
   cd ../step$j
   # Generate the new TOP file
   cp ../SYSTEM_empty.top SYSTEM_step$j.top
   echo "${mol1_name}           $number_mol1" >> SYSTEM_step$j.top
   echo "${mol2_name}           $number_mol2" >> SYSTEM_step$j.top
   echo "${solv_name}           $no_solv_mol" >> SYSTEM_step$j.top

   ### NPT equilibration before run ###
   $GROMPP -f ../martini_v2.x_new_eq_NPT_semiiso.mdp -p SYSTEM_step$j.top -c start_step$j.gro -o NPT_step$j.tpr -po NPT_step$j.mdp -maxwarn 10
   wait
   $MDRUN $FLAGS -deffnm NPT_step$j >> mdrun.log 2>&1
   ### 
   $GROMPP -f ../martini_v2.x_new_run_semiiso.mdp    -p SYSTEM_step$j.top -c NPT_step$j.gro -o run_step$j.tpr -po run_step$j.mdp -maxwarn 10
   wait


######################################
## 2. Run the new evaporation step  ##
######################################
   echo "Run step number" $j

else
   echo "Run step number" $j "from checkpoint"
fi

$MDRUN $FLAGS -cpi run_step${j}_prev.cpt -deffnm run_step$j >> mdrun.log 2>&1

# remove trajectory files 
rm *_step$j.xtc  *_step$j.trr \#*

touch finished
checkpoint_restart="no"


i=$(expr $i + 1)
i=$(printf "%03d" $i)   # makes sure the new "i" is a 3 digits long variable 
#----------------------------------------------------------------------#
done

echo "Evaporation not finished yet! Increase number of steps. NMAX currently = " $NMAX

