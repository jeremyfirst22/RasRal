#!/bin/bash
angBinDist=30

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} >"  
    exit 
}

if [ -z $1 ] ; then 
    usage 
    fi 

fileName=$1 
if [ ! -f $fileName ] ; then 
    echo "ERROR: $fileName not found " 
    exit 
    fi 
if [[ $fileName == *.pdb ]] ; then 
    MOLEC=$(basename $fileName)
    MOLEC=${MOLEC%.*}
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -d $MOLEC ] ; then mkdir $MOLEC ; fi 
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/. ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log 
errFile=$TOP/$MOLEC/$MOLEC.err 
FF=$TOP/GMXFF
forceField=amber03

if [ ! -d $FF/$forceField.ff ] ; then 
    echo ; echo "ERROR: FF not found" 
    exit
    fi  

check(){
    for arg in $@ ; do  
         if [ ! -s $arg ] ; then 
             echo ; echo "ERROR: $arg missing. Exitting" 
             exit 
             fi  
         done 
}

clean(){
    if [ -d $forceField.ff ] ; then rm -r $forceField.ff *.dat ; fi  
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi  

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi  
                                                        
    if [ ! -d $dirName/$forceField.ff ] ; then 
        if [ -d $FF/$forceField.ff ] ; then 
            cp -r $FF/$forceField.ff $dirName
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found" 
            exit 
            fi  
        fi  
}

prep(){
    create_dir Prep 
    cd Prep ; clean 

    protein_steep
    solvate
    solvent_steep
    solvent_nvt
    solvent_npt

    cd ../
}

prep_windows(){
    printf "\t\tPrepping windows: \n" 
    create_dir Prep_windows 
    cd Prep_windows ; clean 
    for angle in `seq 0 $angBinDist 359` ; do 
        prep_angle $angle 
    done 
    cd ../
} 

production_windows(){
    printf "\t\tProduction runs (2 ns): \n"
    create_dir Production
    cd Production ; clean
    for angle in `seq 0 $angBinDist 359` ; do 
        production_run $angle 
    done 
    cd ../
}

protein_steep(){
    printf "\t\tProtein steep............................." 
    if [ ! -f Protein_steep/protein_steep.gro ] ; then 
        create_dir Protein_steep
        
        cp ../$MOLEC.pdb Protein_steep/.
        cd Protein_steep

        gmx pdb2gmx -f $MOLEC.pdb \
            -p $MOLEC.top \
            -ff $forceField \
            -water tip3p \
            -merge all \
            -o $MOLEC.gro >> $logFile 2>> $errFile 
        check $MOLEC.gro 

        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -d 1.5 \
            -bt dodecahedron \
            -o boxed.gro >> $logFile 2>> $errFile
        check boxed.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c boxed.gro \
            -p $MOLEC.top \
            -o protein_steep.tpr >> $logFile 2>> $errFile 
        check protein_steep.tpr 

        gmx mdrun -deffnm protein_steep >> $logFile 2>> $errFile 
        check protein_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi
} 

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        cp Protein_steep/*.itp Solvate/. 
        cd Solvate

        gmx solvate -cp protein_steep.gro \
            -p $MOLEC.top \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        gmx grompp -f $MDP/vac_md.mdp \
            -p $MOLEC.top \
            -c solvated.gro \
            -o genion.tpr >> $logFile 2>> $errFile 
        check genion.tpr
        
        cp $MOLEC.top neutral.top
        echo 'SOL' | gmx genion -s genion.tpr \
            -neutral \
            -nname 'CL' \
            -pname 'NA' \
            -p neutral.top \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

solvent_steep(){
    printf "\t\tSolvent steep............................." 
    if [ ! -f Solvent_steep/solvent_steep.gro ] ; then 
        create_dir Solvent_steep
        
        cp Solvate/neutral.gro Solvent_steep/. 
        cp Solvate/neutral.top Solvent_steep/. 
        #cp Solvate/neutral*.itp Solvent_steep/. 
        cp Solvate/posre*.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep.tpr >> $logFile 2>> $errFile 
        check solvent_steep.tpr 

        gmx mdrun -deffnm solvent_steep >> $logFile 2>> $errFile 
        check solvent_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt(){
    printf "\t\tSolvent NVT relaxation...................." 
    if [ ! -f Solvent_nvt/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt/. 
        cp Solvent_steep/neutral.top Solvent_nvt/. 
        cp Solvent_steep/*.itp Solvent_nvt/. 
        cd Solvent_nvt

        gmx grompp -f $MDP/solvent_nvt_relax.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt(){
    printf "\t\tSolvent NPT relaxation...................." 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_npt_relax.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

prep_angle(){
    if [ -z $1 ] ; then 
        echo "ERROR: Function requires argument." 
        echo "Usage: $0 < center angle of window >" 
        exit 
    fi 
    window=$1 

    printf "\t\t\t%3i..............................." $window
    if [ ! -f $window/$MOLEC.$window.gro ] ; then 
        create_dir $window  
        cp ../Prep/Solvent_npt/solvent_npt.gro $window/. 
        cp ../Prep/Solvent_npt/neutral.top $window/. 
        cd $window 

        create_restraint solvent_npt.gro $window 500 0 

        awk '/Include Position restraint file/{
            print "\n"; 
            print "; Include Chi 1 dihedral restraint" 
            print "#include \"dihrestraint.itp\" \n";
            print;
            next
        }1' neutral.top > $MOLEC.$window.top

        gmx grompp -f $MDP/force_probe_nvt.mdp \
            -p $MOLEC.$window.top \
            -c solvent_npt.gro \
            -o $MOLEC.$window.tpr >> $logFile 2>> $errFile 
        check $MOLEC.$window.tpr 

        gmx mdrun -deffnm $MOLEC.$window >> $logFile 2>> $errFile 
        check $MOLEC.$window.gro 

        mkdir Verify ; cd Verify 
        gmx angle -f ../$MOLEC.$window.xtc \
            -n ../dihedral.ndx \
            -type dihedral \
            -ov >> $logFile 2>> $errFile 
        check angaver.xvg 
        cd ../

        clean
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
    fi
}

create_restraint(){
    if [[ -z $1 || -z $2 || -z $3 || -z $4 ]] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < Structure file > < window (degrees) > < force constant (kJ/mol/rad^2) >  < flat potential range (deg) " 
        exit 
    fi 

    Q61=165

    sFile=$1
    r1=$2
    kfac=$3
    dPhi=$4

    if [[ $Q61 -lt 100 || $Q61 -gt 999 ]] ; then 
        echo "ERROR: Indexing for Q61 is going to fail. We are assuming it is three digits here." 
        exit 
    fi 

    ai=`grep "^  $Q61..." $sFile | grep -w N  | awk '{print $3}'`
    aj=`grep "^  $Q61..." $sFile | grep -w CA | awk '{print $3}'`
    ak=`grep "^  $Q61..." $sFile | grep -w CB | awk '{print $3}'`  #These three atoms are consistent for all A.A.s

    letterMut="${MOLEC: -1}"

    if   [ "$letterMut" = "C" ] ; then atom4='SG' 
    elif [ "$letterMut" = "I" ] ; then atom4='CG1'
    elif [ "$letterMut" = "S" ] ; then atom4='OG'
    elif [ "$letterMut" = "T" ] ; then atom4='OG1' 
    elif [ "$letterMut" = "V" ] ; then atom4='CG1' 
    else atom4='CG' ; fi 

    if ! grep -sq "^  $Q61..." $sFile | grep -w $atom4 | awk '{print $3}' ; then 
        echo "ERROR: Atom 4 of Chi 1 not found" 
        exit 
    else 
        al=`grep "^  $Q61..." $sFile | grep -w $atom4 | awk '{print $3}'`
    fi 
    
    echo "[ dihedral_restraints ]" > dihrestraint.itp 
    printf ";%6s%6s%6s%6s%8s%8s%8s%12s\n" ai aj ak al func phi dphi kfac >> dihrestraint.itp
    printf " %6i%6i%6i%6i%8i%8i%8i%12i\n" $ai $aj $ak $al 1 $r1 $dPhi $kfac >> dihrestraint.itp
    
    echo "[ X1 ]" > dihedral.ndx 
    echo "$ai $aj $ak $al" >> dihedral.ndx
} 

production_run(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%3i..............................." $window
    if [ ! -f $window/$MOLEC.$window.run_$run.gro ] ; then
        create_dir $angle 
        cp ../Prep_windows/$window/$MOLEC.$window.gro $window/equilibrated_starting.$window.gro
        cp ../Prep_windows/$window/$MOLEC.$window.top $window/. 
        cd $window
        
        ##Each window is flat within 15 degrees, and has a 70 kJ/mol/rad^2 harmondic restraint outside
        create_restraint equilibrated_starting.$window.gro $window 70 15
        check dihrestraint.itp dihedral.ndx 

        if [ ! -f $MOLEC.$window.gro ] ; then 
            if [ ! -f $MOLEC.$window.tpr ] ; then 
                gmx grompp -f $MDP/production.mdp \
                    -p $MOLEC.$window.top \
                    -c equilibrated_starting.$window.gro \
                    -o $MOLEC.$window.tpr >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.$window.tpr

            if [ -f $MOLEC.$window.cpt ] ; then 
                ibrun mdrun_mpi -deffnm $MOLEC.$window \
                    -cpi $MOLEC.$window.cpt >> $logFile 2>> $errFile  
            else 
                ibrun mdrun_mpi -deffnm $MOLEC.$window >> $logFile 2>> $errFile 
                fi 
            fi 
        #check $MOLEC.$window.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
prep
prep_windows
production_windows
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
