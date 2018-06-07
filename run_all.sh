#!/bin/bash

#molecList="C D E F H I K L M N Q R S T V W Y" 
molecList="Q"

for mut in $molecList ; do 
    molec=RasRalC18CNC_Q61$mut
    if [ ! -f StartingStructures/$molec.pdb ] ; then 
        echo "$molec.pdb not found!" 
        continue 
        fi 
    
    printf "\n\n\t\t*** $molec ***\n\n" 


    echo "#!/bin/bash" > submit_Q61$mut
    echo >> submit_Q61$mut
    echo "#SBATCH -J Q61$mut " >> submit_Q61$mut
    echo "#SBATCH -o Q61$mut.o%j" >> submit_Q61$mut
    echo "#SBATCH -N 1" >> submit_Q61$mut
    echo "#SBATCH -n 16 " >> submit_Q61$mut
    echo "#SBATCH -p skx-normal " >> submit_Q61$mut
    echo "#SBATCH -t 48:00:00" >> submit_Q61$mut
    echo "#SBATCH -A Ras" >> submit_Q61$mut
    echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_Q61$mut
    echo "#SBATCH --mail-type=all" >> submit_Q61$mut
    
    echo >> submit_Q61$mut
    echo "module load gromacs " >> submit_Q61$mut 
   
    echo >> submit_Q61$mut
    echo "bash run_RasRal_Umbrella.sh StartingStructures/${molec}.pdb" >> submit_Q61$mut
   
    sbatch submit_$molec
    #bash run_RasRal_Umbrella.sh StartingStructures/$molec.pdb


    done 
