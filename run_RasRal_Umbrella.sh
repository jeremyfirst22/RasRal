#!/bin/bash
angBinDist=30
Q61=165
FORCE_TOOLS=/Users/jfirst/force_calc_tools


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

analysis_windows(){
    printf "\t\tAnalysis: \n" 
    create_dir Analysis
    cd Analysis ; clean 
    for angle in `seq 0 $angBinDist 359` ; do 
    printf "\t\t\t%3i :\n" $angle
        sasa $angle
        polar_sasa $angle 
        chi1 $angle 
        nopbc $angle
    done
    wham
    boltzmann_weight
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

        create_restraint solvent_npt.gro $window 100 0 

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
    if [ ! -f $window/$MOLEC.$window.gro ] ; then
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
        check $MOLEC.$window.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

sasa(){ 
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "SASA"
    if [[ ! -f sasa/area.$window.xvg || ! -f sasa/sidechain.$window.xvg ]]  ; then 
        create_dir sasa
        cd sasa

        if [ ! -f area.$window.xvg ] ; then 
            gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -surface 'Protein' \
                -output "resindex $Q61" \
                -ndots 240 \
                -o area.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check area.$window.xvg 

        if [ ! -f sidechain.$window.xvg ] ; then 
            gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -surface 'Protein' \
                -output "group Sidechain and resindex $Q61" \
                -ndots 240 \
                -o sidechain.$window.xvg >> $logFile 2>> $errFile 
        fi
        check sidechain.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

polar_sasa(){ 
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    sideChainPolar=true
    case "${MOLEC: -1}" in 
        D) 
            polarAtoms="resindex $Q61 and name OD1 OD2" 
          ;; 
        E) 
            polarAtoms="resindex $Q61 and name OE1 OE2" 
          ;; 
        H)  ##Q61H was assigned rtp HIE entry, so we need those atoms
            polarAtoms="resindex $Q61 and name ND1 NE2 HE2"
          ;; 
        K) 
            polarAtoms="resindex $Q61 and name NZ HZ1 HZ2 HZ3"
          ;; 
        N) 
            polarAtoms="resindex $Q61 and name OD1 ND2 HD21 HD22"
          ;; 
        Q) 
            polarAtoms="resindex $Q61 and name OE1 NE2 HE21 HE22"
          ;; 
        R) 
            polarAtoms="resindex $Q61 and name NE NH1 NH2 HH11 HH12 HH21 HH22 HE"
          ;;
        S) 
            polarAtoms="resindex $Q61 and name OG HG" 
          ;; 
        T) 
            polarAtoms="resindex $Q61 and name OG1 HG1"
          ;; 
        W) 
            polarAtoms="resindex $Q61 and name NE1 HE1" 
          ;; 
        Y) 
            polarAtoms="resindex $Q61 and name OH HH" 
          ;; 
        *) 
            sideChainPolar=false
            polarAtoms="resindex $Q61 and name"
          ;; 
    esac 

    printf "\t\t\t%10s........................" "Polar SASA"
    if [ ! -f polar_sasa/polar.$window.xvg ]  ; then 
        create_dir polar_sasa
        cd polar_sasa

        if [ $sideChainPolar = true ] ; then 
            if [ ! -f sc_polar.$window.xvg ] ; then 
                gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                    -s ../../Production/$window/$MOLEC.$window.tpr \
                    -surface 'Protein' \
                    -output "$polarAtoms" \
                    -ndots 240 \
                    -o sc_polar.$window.xvg >> $logFile 2>> $errFile 
            fi 
            check sc_polar.$window.xvg 

            if [ ! -f davids.$window.xvg ] ; then 
                gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                    -select "$polarAtoms" \
                    -on selection.$window.ndx >> $logFile 2>> $errFile 
                check selection.$window.ndx 

                echo '0' | gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                    -s ../../Production/$window/$MOLEC.$window.tpr \
                    -n selection.$window.ndx \
                    -o davids.$window.xvg >> $logFile 2>> $errFile 
                check davids.$window.xvg 
            fi 
        fi 

        polarAtoms=$polarAtoms" N H O"   ##Append three backbone polar atoms
        if [ ! -f polar.$window.xvg ] ; then 
            gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -surface 'Protein' \
                -output "$polarAtoms" \
                -ndots 240 \
                -o polar.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check polar.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

chi1(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "Chi 1" 
    if [ ! -f chi1/angaver.$window.xvg ]  ; then 
        create_dir chi1 
        cp ../Production/$window/dihedral.ndx chi1/. 
        cd chi1

        ##XVG none required for WHAM analysis (easier than cleaning each file) 
        gmx angle -f ../../Production/$window/$MOLEC.$window.xtc \
            -n dihedral.ndx \
            -type dihedral \
            -xvg none \
            -od angdist.$window.xvg \
            -ov angaver.$window.xvg >> $logFile 2>> $errFile 
        check angaver.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

nopbc(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "No pbc" 
    if [ ! -f nopbc/nopbc.$window.xtc ]  ; then 
        create_dir nopbc
        cd nopbc

        echo 'Protein System' | gmx trjconv -f ../../Production/$window/equilibrated_starting.$window.gro \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -center \
            -pbc mol \
            -ur compact \
            -o nopbc.$window.gro >> $logFile 2>> $errFile 
        check nopbc.$window.gro 

        echo 'Protein System' | gmx trjconv -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -center \
            -pbc mol \
            -ur compact \
            -o nopbc.$window.xtc >> $logFile 2>> $errFile 
        check nopbc.$window.xtc 
            
        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

wham(){
    printf "\t\t\t%-4s.............................." "WHAM" 
    if [[ ! -f wham/$MOLEC.output.prob && -f ~/wham/wham-1.0/WHAM ]]  ; then 
        create_dir wham
        cd wham ; clean 

        i='0' 
        if [ -f $MOLEC.wham.inp ] ; then rm $MOLEC.wham.inp ; fi 
        touch $MOLEC.wham.inp 

        for window in `seq 0 $angBinDist 359` ; do 
            ##    index      fileName       center dphi kfac Temp 
            echo "$i ../chi1/angaver.$window.xvg $window 15 70 300" >> $MOLEC.wham.inp  
            ((i++)) 
        done 

        ~/wham/wham-1.0/WHAM --f $MOLEC.wham.inp --o $MOLEC.output --b 1 >> $logFile 2>> $errFile 

        check $MOLEC.output.prob 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

boltzmann_weight(){
    printf "\t\t\t%-16s.................." "Boltzmann weight"
    if [[ ! -f boltzmann/angaver.weighted.out \
        || ! -f boltzmann/area.weighted.out \
        || ! -f boltzmann/sidechain.weighted.out ]] \
        || [[ ! -f boltzmann/sc_polar.weighted.out && -f polar_sasa/sc_polar.330.xvg ]] \
        || [[ ! -f boltzmann/davids.weighted.out && -f polar_sasa/davids.330.xvg ]] \
        || [[ ! -f boltzmann/polar.weighted.out && -f polar_sasa/polar.330.xvg ]] \
        && [ -f $FORCE_TOOLS/boltzmann_weight ]  ; then 
        create_dir boltzmann
        cd boltzmann ; clean 

        if [ ! -f angaver.weighted.out ] ; then 
            if [ -f angaver.boltzmann.inp ] ; then rm angaver.boltzmann.inp ; fi 

            i=0
            touch angaver.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                echo "../WHAM/$MOLEC.output.$i.bin   ../chi1/angaver.$window.xvg" >> angaver.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l angaver.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o angaver.weighted.out >> $logFile 2>> $errFile 
        fi 
        check angaver.weighted.out 

        if [ ! -f area.weighted.out ] ; then 
            if [ -f area.boltzmann.inp ] ; then rm area.boltzmann.inp ; fi 

            i=0
            touch area.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../sasa/area.$window.xvg area.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   area.$window.xvg" >> area.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l area.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o area.weighted.out >> $logFile 2>> $errFile 
        fi 
        check area.weighted.out 

        if [ ! -f sidechain.weighted.out ] ; then 
            if [ -f sidechain.boltzmann.inp ] ; then rm sidechain.boltzmann.inp ; fi 

            i=0
            touch sidechain.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../sasa/sidechain.$window.xvg sidechain.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   sidechain.$window.xvg" >> sidechain.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l sidechain.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o sidechain.weighted.out >> $logFile 2>> $errFile 
        fi 
        check sidechain.weighted.out 

        if [[ ! -f sc_polar.weighted.out && -f ../polar_sasa/sc_polar.$window.xvg ]] ; then 
            if [ -f sc_polar.boltzmann.inp ] ; then rm sc_polar.boltzmann.inp ; fi 

            i=0
            touch sc_polar.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../polar_sasa/sc_polar.$window.xvg sc_polar.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   sc_polar.$window.xvg" >> sc_polar.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l sc_polar.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o sc_polar.weighted.out >> $logFile 2>> $errFile 
            check sc_polar.weighted.out 
        fi 

        if [[ ! -f davids.weighted.out && -f ../polar_sasa/davids.$window.xvg ]] ; then 
            if [ -f davids.boltzmann.inp ] ; then rm davids.boltzmann.inp ; fi 

            i=0
            touch davids.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../polar_sasa/davids.$window.xvg davids.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   davids.$window.xvg" >> davids.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l davids.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o davids.weighted.out >> $logFile 2>> $errFile 
            check davids.weighted.out 
        fi 

        if [[ ! -f polar.weighted.out && -f ../polar_sasa/polar.$window.xvg ]] ; then 
            if [ -f polar.boltzmann.inp ] ; then rm polar.boltzmann.inp ; fi 

            i=0
            touch polar.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../polar_sasa/polar.$window.xvg polar.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   polar.$window.xvg" >> polar.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l polar.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o polar.weighted.out >> $logFile 2>> $errFile 
            check polar.weighted.out 
        fi 

        #check endFile.$window.xvg
        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

clean_xvg(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: clean_xvg < target file > < new cleaned file> "
    fi 
    cat $1 | grep -v "^@" | grep -v "^#" > $2 
}

template(){
    ##Call from analysis_windows()
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "Template" 
    if [ ! -f template/endFile.$window.xvg ]  ; then 
        create_dir template
        cd template

        #analysis function that produces endFile

        check endFile.$window.xvg
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
analysis_windows
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
