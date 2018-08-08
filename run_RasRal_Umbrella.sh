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
        nearest_water $angle 
        dist_to_Q69 $angle 
        dist_to_Y32 $angle 
        rmsd_gtp $angle 
        rmsd_l4loop $angle 
        gmx_hbond $angle 
        volume_Q61 $angle
        hbond_Q61 $angle 
        hbond_Y32 $angle 
        chi1_Y32 $angle
        scount $angle 
        scount_O2G $angle 
        scount_vary_O1G $angle 
        scount_vary_SCD $angle 
        scount_both $angle 
        scount_polar $angle 
        scount_2 $angle 
        force_nitrile $angle 
        g12_os3_dist $angle 
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
    if [[ ! -f sasa/area.$window.xvg || ! -f sasa/sidechain.$window.xvg || ! -f sasa/isolated.$window.xvg ]]  ; then 
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

        if [ ! -f isolated.$window.xvg ] ; then 
            if [ ! -f residue.ndx ] ; then 
                gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                    -select "group \"SideChain\" and resindex $Q61" \
                    -on residue.ndx >> $logFile 2>> $errFile 
            fi 
            echo '0' | gmx sasa -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -n residue.ndx \
                -ndots 240 \
                -o isolated.$window.xvg >> $logFile 2>> $errFile 
        fi
        check isolated.$window.xvg 

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
        X)  ##rtp HIP; charged Histidine  
            polarAtoms="resindex $Q61 and name ND1 NE2 HE2 HD1"
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
            -pbc whole \
            -ur compact \
            -o nopbc.$window.gro >> $logFile 2>> $errFile 
        check nopbc.$window.gro 

        echo 'Protein System' | gmx trjconv -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -center \
            -pbc whole \
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
        || ! -f boltzmann/isolated_sc.weighted.out \
        || ! -f boltzmann/mindist.weighted.out \
        || ! -f boltzmann/distQ69.weighted.out \
        || ! -f boltzmann/distY32.weighted.out \
        || ! -f boltzmann/rmsd.weighted.out \
        || ! -f boltzmann/l4loop.weighted.out \
        || ! -f boltzmann/active.weighted.out \
        || ! -f boltzmann/hbnum.weighted.out \
        || ! -f boltzmann/volume.weighted.out \
        || ! -f boltzmann/hbnum_Q61.weighted.out \
        || ! -f boltzmann/hbnum_Y32.weighted.out \
        || ! -f boltzmann/angaver_Y32.weighted.out \
        || ! -f boltzmann/size_245.weighted.out \
        || ! -f boltzmann/size_3.weighted.out \
        || ! -f boltzmann/size_4.weighted.out \
        || ! -f boltzmann/size_5.weighted.out \
        || ! -f boltzmann/both_size_245.weighted.out \
        || ! -f boltzmann/both_size_3.weighted.out \
        || ! -f boltzmann/both_size_4.weighted.out \
        || ! -f boltzmann/both_size_5.weighted.out \
        || ! -f boltzmann/size2.weighted.out \
        || ! -f boltzmann/external_field.weighted.out \
        || ! -f boltzmann/sidechain.weighted.out ]] \
        || [[ ! -f boltzmann/sc_polar.weighted.out && -f polar_sasa/sc_polar.330.xvg ]] \
        || [[ ! -f boltzmann/size_polar_245.weighted.out && -f scount_polar/size.330.xvg ]] \
        || [[ ! -f boltzmann/size_polar_3.weighted.out && -f scount_polar/size.330.xvg ]] \
        || [[ ! -f boltzmann/size_polar_4.weighted.out && -f scount_polar/size.330.xvg ]] \
        || [[ ! -f boltzmann/size_polar_5.weighted.out && -f scount_polar/size.330.xvg ]] \
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

        if [ ! -f isolated_sc.weighted.out ] ; then 
            if [ -f isolated_sc.boltzmann.inp ] ; then rm isolated_sc.boltzmann.inp ; fi 

            i=0
            touch isolated_sc.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../sasa/isolated.$window.xvg isolated_sc.$window.xvg 
                echo "../wham/$MOLEC.output.$i.bin   isolated_sc.$window.xvg" >> isolated_sc.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l isolated_sc.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o isolated_sc.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check isolated_sc.weighted.out 

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

        if [ ! -f mindist.weighted.out ] ; then 
            if [ -f mindist.boltzmann.inp ] ; then rm mindist.boltzmann.inp ; fi 

            i=0
            touch mindist.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../mindist/mindist.$window.xvg mindist.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   mindist.$window.xvg" >> mindist.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l mindist.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o mindist.weighted.out >> $logFile 2>> $errFile 
        fi 
        check mindist.weighted.out 

        if [ ! -f distQ69.weighted.out ] ; then 
            if [ -f distQ69.boltzmann.inp ] ; then rm distQ69.boltzmann.inp ; fi 

            i=0
            touch distQ69.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../distQ69/distQ69.$window.xvg distQ69.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   distQ69.$window.xvg" >> distQ69.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l distQ69.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o distQ69.weighted.out >> $logFile 2>> $errFile 
        fi 
        check distQ69.weighted.out 

        if [ ! -f distY32.weighted.out ] ; then 
            if [ -f distY32.boltzmann.inp ] ; then rm distY32.boltzmann.inp ; fi 

            i=0
            touch distY32.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../distY32/distY32.$window.xvg distY32.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   distY32.$window.xvg" >> distY32.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l distY32.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o distY32.weighted.out >> $logFile 2>> $errFile 
        fi 
        check distY32.weighted.out 

        if [ ! -f rmsd.weighted.out ] ; then 
            if [ -f rmsd.boltzmann.inp ] ; then rm rmsd.boltzmann.inp ; fi 

            i=0
            touch rmsd.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../rmsd_gtp/rmsd.$window.xvg rmsd.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   rmsd.$window.xvg" >> rmsd.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l rmsd.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o rmsd.weighted.out >> $logFile 2>> $errFile 
        fi 
        check rmsd.weighted.out 

        if [ ! -f l4loop.weighted.out ] ; then 
            if [ -f l4loop.boltzmann.inp ] ; then rm l4loop.boltzmann.inp ; fi 

            i=0
            touch l4loop.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../rmsd_l4loop/rmsd.$window.xvg l4loop.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   l4loop.$window.xvg" >> l4loop.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l l4loop.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o l4loop.weighted.out >> $logFile 2>> $errFile 
        fi 
        check l4loop.weighted.out 

        if [ ! -f active.weighted.out ] ; then 
            if [ -f active.boltzmann.inp ] ; then rm active.boltzmann.inp ; fi 

            i=0
            touch active.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../rmsd_gtp/active.$window.xvg active.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   active.$window.xvg" >> active.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l active.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o active.weighted.out >> $logFile 2>> $errFile 
        fi 
        check active.weighted.out 

        if [ ! -f hbnum.weighted.out ] ; then 
            if [ -f hbnum.boltzmann.inp ] ; then rm hbnum.boltzmann.inp ; fi 

            i=0
            touch hbnum.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../gmx_hbond/hbnum.$window.xvg hbnum.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   hbnum.$window.xvg" >> hbnum.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l hbnum.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o hbnum.weighted.out >> $logFile 2>> $errFile 
        fi 
        check hbnum.weighted.out 

        if [ ! -f volume.weighted.out ] ; then 
            if [ -f volume.boltzmann.inp ] ; then rm volume.boltzmann.inp ; fi 

            i=0
            touch volume.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../volume_Q61/volume.$window.xvg volume.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   volume.$window.xvg" >> volume.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l volume.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o volume.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check volume.weighted.out 

        if [ ! -f hbnum_Q61.weighted.out ] ; then 
            if [ -f hbnum_Q61.boltzmann.inp ] ; then rm hbnum_Q61.boltzmann.inp ; fi 

            i=0
            touch hbnum_Q61.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../hbond_Q61/hbnum_Q61.$window.xvg hbnum_Q61.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   hbnum_Q61.$window.xvg" >> hbnum_Q61.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l hbnum_Q61.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o hbnum_Q61.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check hbnum_Q61.weighted.out 

        if [ ! -f hbnum_Y32.weighted.out ] ; then 
            if [ -f hbnum_Y32.boltzmann.inp ] ; then rm hbnum_Y32.boltzmann.inp ; fi 

            i=0
            touch hbnum_Y32.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../hbond_Y32/hbnum_Y32.$window.xvg hbnum_Y32.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   hbnum_Y32.$window.xvg" >> hbnum_Y32.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l hbnum_Y32.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o hbnum_Y32.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check hbnum_Y32.weighted.out 

        if [ ! -f angaver_Y32.weighted.out ] ; then 
            if [ -f angaver_Y32.boltzmann.inp ] ; then rm angaver_Y32.boltzmann.inp ; fi 

            i=0
            touch angaver_Y32.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../chi1_Y32/angaver_Y32.$window.xvg angaver_Y32.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   angaver_Y32.$window.xvg" >> angaver_Y32.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l angaver_Y32.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o angaver_Y32.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check angaver_Y32.weighted.out 

        if [ ! -f size_245.weighted.out ] ; then 
            if [ -f size_245.boltzmann.inp ] ; then rm size_245.boltzmann.inp ; fi 

            i=0
            touch size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_245.$window.xvg size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_245.$window.xvg" >> size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check size_245.weighted.out 

        if [ ! -f size_3.weighted.out ] ; then 
            if [ -f size_3.boltzmann.inp ] ; then rm size_3.boltzmann.inp ; fi 

            i=0
            touch size_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_3.$window.xvg size_3.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_3.$window.xvg" >> size_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check size_3.weighted.out 

        if [ ! -f size_4.weighted.out ] ; then 
            if [ -f size_4.boltzmann.inp ] ; then rm size_4.boltzmann.inp ; fi 

            i=0
            touch size_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_4.$window.xvg size_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_4.$window.xvg" >> size_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check size_4.weighted.out 

        if [ ! -f size_5.weighted.out ] ; then 
            if [ -f size_5.boltzmann.inp ] ; then rm size_5.boltzmann.inp ; fi 

            i=0
            touch size_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_5.$window.xvg size_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_5.$window.xvg" >> size_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check size_5.weighted.out 

        if [ ! -f size_6.weighted.out ] ; then 
            if [ -f size_6.boltzmann.inp ] ; then rm size_6.boltzmann.inp ; fi 

            i=0
            touch size_6.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_6.$window.xvg size_6.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_6.$window.xvg" >> size_6.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_6.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_6.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check size_6.weighted.out 

        if [ ! -f O1G_size_245.weighted.out ] ; then 
            if [ -f O1G_size_245.boltzmann.inp ] ; then rm O1G_size_245.boltzmann.inp ; fi 

            i=0
            touch O1G_size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O1G/size_245.$window.xvg O1G_size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O1G_size_245.$window.xvg" >> O1G_size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O1G_size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O1G_size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O1G_size_245.weighted.out 

        if [ ! -f O1G_size_3.weighted.out ] ; then 
            if [ -f O1G_size_3.boltzmann.inp ] ; then rm O1G_size_3.boltzmann.inp ; fi 

            i=0
            touch O1G_size_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O1G/size_3.$window.xvg O1G_size_3.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O1G_size_3.$window.xvg" >> O1G_size_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O1G_size_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O1G_size_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O1G_size_3.weighted.out 

        if [ ! -f O1G_size_4.weighted.out ] ; then 
            if [ -f O1G_size_4.boltzmann.inp ] ; then rm O1G_size_4.boltzmann.inp ; fi 

            i=0
            touch O1G_size_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O1G/size_4.$window.xvg O1G_size_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O1G_size_4.$window.xvg" >> O1G_size_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O1G_size_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O1G_size_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O1G_size_4.weighted.out 

        if [ ! -f O1G_size_5.weighted.out ] ; then 
            if [ -f O1G_size_5.boltzmann.inp ] ; then rm O1G_size_5.boltzmann.inp ; fi 

            i=0
            touch O1G_size_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O1G/size_5.$window.xvg O1G_size_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O1G_size_5.$window.xvg" >> O1G_size_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O1G_size_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O1G_size_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O1G_size_5.weighted.out 

        if [ ! -f O1G_size_6.weighted.out ] ; then 
            if [ -f O1G_size_6.boltzmann.inp ] ; then rm O1G_size_6.boltzmann.inp ; fi 

            i=0
            touch O1G_size_6.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O1G/size_6.$window.xvg O1G_size_6.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O1G_size_6.$window.xvg" >> O1G_size_6.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O1G_size_6.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O1G_size_6.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O1G_size_6.weighted.out 

        if [ ! -f SCD_size_245.weighted.out ] ; then 
            if [ -f SCD_size_245.boltzmann.inp ] ; then rm SCD_size_245.boltzmann.inp ; fi 

            i=0
            touch SCD_size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_SCD/size_245.$window.xvg SCD_size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   SCD_size_245.$window.xvg" >> SCD_size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l SCD_size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o SCD_size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check SCD_size_245.weighted.out 

        if [ ! -f O2G_size_245.weighted.out ] ; then 
            if [ -f O2G_size_245.boltzmann.inp ] ; then rm O2G_size_245.boltzmann.inp ; fi 

            i=0
            touch O2G_size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_245.$window.xvg O2G_size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_245.$window.xvg" >> O2G_size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_245.weighted.out 
        if [ ! -f O2G_size_3.weighted.out ] ; then 
            if [ -f O2G_size_3.boltzmann.inp ] ; then rm O2G_size_3.boltzmann.inp ; fi 

            i=0
            touch O2G_size_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_3.$window.xvg O2G_size_3.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_3.$window.xvg" >> O2G_size_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_3.weighted.out 
        if [ ! -f O2G_size_4.weighted.out ] ; then 
            if [ -f O2G_size_4.boltzmann.inp ] ; then rm O2G_size_4.boltzmann.inp ; fi 

            i=0
            touch O2G_size_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_4.$window.xvg O2G_size_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_4.$window.xvg" >> O2G_size_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_4.weighted.out 
        if [ ! -f O2G_size_5.weighted.out ] ; then 
            if [ -f O2G_size_5.boltzmann.inp ] ; then rm O2G_size_5.boltzmann.inp ; fi 

            i=0
            touch O2G_size_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_5.$window.xvg O2G_size_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_5.$window.xvg" >> O2G_size_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_5.weighted.out 
        if [ ! -f O2G_size_6.weighted.out ] ; then 
            if [ -f O2G_size_6.boltzmann.inp ] ; then rm O2G_size_6.boltzmann.inp ; fi 

            i=0
            touch O2G_size_6.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_6.$window.xvg O2G_size_6.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_6.$window.xvg" >> O2G_size_6.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_6.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_6.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_6.weighted.out 


        if [ ! -f SCD_size_3.weighted.out ] ; then 
            if [ -f SCD_size_3.boltzmann.inp ] ; then rm SCD_size_3.boltzmann.inp ; fi 

            i=0
            touch SCD_size_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_SCD/size_3.$window.xvg SCD_size_3.$window.xvg 
        if [ ! -f O2G_size_245.weighted.out ] ; then 
            if [ -f O2G_size_245.boltzmann.inp ] ; then rm O2G_size_245.boltzmann.inp ; fi 

            i=0
            touch O2G_size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_O2G/size_245.$window.xvg O2G_size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   O2G_size_245.$window.xvg" >> O2G_size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l O2G_size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o O2G_size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check O2G_size_245.weighted.out 
                echo "../WHAM/$MOLEC.output.$i.bin   SCD_size_3.$window.xvg" >> SCD_size_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l SCD_size_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o SCD_size_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check SCD_size_3.weighted.out 

        if [ ! -f SCD_size_4.weighted.out ] ; then 
            if [ -f SCD_size_4.boltzmann.inp ] ; then rm SCD_size_4.boltzmann.inp ; fi 

            i=0
            touch SCD_size_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_SCD/size_4.$window.xvg SCD_size_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   SCD_size_4.$window.xvg" >> SCD_size_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l SCD_size_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o SCD_size_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check SCD_size_4.weighted.out 

        if [ ! -f SCD_size_5.weighted.out ] ; then 
            if [ -f SCD_size_5.boltzmann.inp ] ; then rm SCD_size_5.boltzmann.inp ; fi 

            i=0
            touch SCD_size_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_SCD/size_5.$window.xvg SCD_size_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   SCD_size_5.$window.xvg" >> SCD_size_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l SCD_size_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o SCD_size_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check SCD_size_5.weighted.out 

        if [ ! -f SCD_size_6.weighted.out ] ; then 
            if [ -f SCD_size_6.boltzmann.inp ] ; then rm SCD_size_6.boltzmann.inp ; fi 

            i=0
            touch SCD_size_6.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_SCD/size_6.$window.xvg SCD_size_6.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   SCD_size_6.$window.xvg" >> SCD_size_6.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l SCD_size_6.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o SCD_size_6.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check SCD_size_6.weighted.out 

        if [ ! -f both_size_245.weighted.out ] ; then 
            if [ -f both_size_245.boltzmann.inp ] ; then rm both_size_245.boltzmann.inp ; fi 

            i=0
            touch both_size_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_both/size_245.$window.xvg both_size_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   both_size_245.$window.xvg" >> both_size_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l both_size_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o both_size_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check both_size_245.weighted.out 

        if [ ! -f both_size_3.weighted.out ] ; then 
            if [ -f both_size_3.boltzmann.inp ] ; then rm both_size_3.boltzmann.inp ; fi 

            i=0
            touch both_size_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_both/size_3.$window.xvg both_size_3.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   both_size_3.$window.xvg" >> both_size_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l both_size_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o both_size_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check both_size_3.weighted.out 

        if [ ! -f both_size_4.weighted.out ] ; then 
            if [ -f both_size_4.boltzmann.inp ] ; then rm both_size_4.boltzmann.inp ; fi 

            i=0
            touch both_size_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount/size_4.$window.xvg both_size_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   both_size_4.$window.xvg" >> both_size_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l both_size_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o both_size_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check both_size_4.weighted.out 

        if [ ! -f both_size_5.weighted.out ] ; then 
            if [ -f both_size_5.boltzmann.inp ] ; then rm both_size_5.boltzmann.inp ; fi 

            i=0
            touch both_size_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_both/size_5.$window.xvg both_size_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   both_size_5.$window.xvg" >> both_size_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l both_size_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o both_size_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check both_size_5.weighted.out 

        if [[ ! -f size_polar_245.weighted.out && -f ../scount_polar/size_245.$window.xvg ]] ; then 
            if [ -f size_polar_245.boltzmann.inp ] ; then rm size_polar_245.boltzmann.inp ; fi 

            i=0
            touch size_polar_245.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_polar/size_245.$window.xvg size_polar_245.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_polar_245.$window.xvg" >> size_polar_245.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_polar_245.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_polar_245.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
            check size_polar_245.weighted.out
        fi 

        if [[ ! -f size_polar_3.weighted.out && -f ../scount_polar/size_3.$window.xvg ]] ; then 
            if [ -f size_polar_3.boltzmann.inp ] ; then rm size_polar_3.boltzmann.inp ; fi 

            i=0
            touch size_polar_3.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_polar/size_3.$window.xvg size_polar_3.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_polar_3.$window.xvg" >> size_polar_3.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_polar_3.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_polar_3.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
            check size_polar_3.weighted.out
        fi 

        if [[ ! -f size_polar_4.weighted.out && -f ../scount_polar/size_4.$window.xvg ]] ; then 
            if [ -f size_polar_4.boltzmann.inp ] ; then rm size_polar_4.boltzmann.inp ; fi 

            i=0
            touch size_polar_4.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_polar/size_4.$window.xvg size_polar_4.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_polar_4.$window.xvg" >> size_polar_4.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_polar_4.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_polar_4.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
            check size_polar_4.weighted.out
        fi 

        if [[ ! -f size_polar_5.weighted.out && -f ../scount_polar/size_5.$window.xvg ]] ; then 
            if [ -f size_polar_5.boltzmann.inp ] ; then rm size_polar_5.boltzmann.inp ; fi 

            i=0
            touch size_polar_5.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_polar/size_5.$window.xvg size_polar_5.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size_polar_5.$window.xvg" >> size_polar_5.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size_polar_5.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size_polar_5.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
            check size_polar_5.weighted.out
        fi 

        if [ ! -f size2.weighted.out ] ; then 
            if [ -f size2.boltzmann.inp ] ; then rm size2.boltzmann.inp ; fi 

            i=0
            touch size2.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../scount_2/size.$window.xvg size2.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   size2.$window.xvg" >> size2.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l size2.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o size2.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 

        if [ ! -f external_field_nitrile.weighted.out ] ; then 
            if [ -f external_field_nitrile.boltzmann.inp ] ; then rm external_field_nitrile.boltzmann.inp ; fi 

            i=0
            touch external_field_nitrile.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                #clean_xvg ../force_nitrile/$MOLEC.$window.external_field_nitrile.projected.xvg external_field_nitrile.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   ../force_nitrile/$MOLEC.$window.external_field.projected.xvg" >> external_field_nitrile.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l external_field_nitrile.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o external_field_nitrile.weighted.out >> $logFile 2>> $errFile 
            #rm *.xvg 
        fi 
        check external_field_nitrile.weighted.out 

        if [ ! -f g12_hbnum.weighted.out ] ; then 
            if [ -f g12_hbnum.boltzmann.inp ] ; then rm g12_hbnum.boltzmann.inp ; fi 

            i=0
            touch g12_hbnum.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../g12_os3_dist/hbnum.$window.xvg g12_hbnum.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   g12_hbnum.$window.xvg" >> g12_hbnum.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l g12_hbnum.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o g12_hbnum.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check g12_hbnum.weighted.out 

        if [ ! -f g12_distave.weighted.out ] ; then 
            if [ -f g12_distave.boltzmann.inp ] ; then rm g12_distave.boltzmann.inp ; fi 

            i=0
            touch g12_distave.boltzmann.inp 
            for window in `seq 0 $angBinDist 359` ; do 
                clean_xvg ../g12_os3_dist/distave.$window.xvg g12_distave.$window.xvg 
                echo "../WHAM/$MOLEC.output.$i.bin   g12_distave.$window.xvg" >> g12_distave.boltzmann.inp 
                ((i++)) 
            done 

            $FORCE_TOOLS/boltzmann_weight -l g12_distave.boltzmann.inp \
                -p ../wham/$MOLEC.output.prob \
                -o g12_distave.weighted.out >> $logFile 2>> $errFile 
            rm *.xvg 
        fi 
        check g12_distave.weighted.out 

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

nearest_water(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "mindist" 
    if [ ! -f mindist/mindist.$window.xvg ]  ; then 
        create_dir mindist
        cd mindist

        if [ ! -f gtp.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "resname GTP and name PG" \
                -on gtp.ndx >> $logFile 2>> $errFile 
        fi 
        check gtp.ndx 

        if [ ! -f water.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "Water" \
                -on water.ndx >> $logFile 2>> $errFile 
        fi
        check water.ndx 

        cat gtp.ndx > index.ndx 
        cat water.ndx >> index.ndx 

        echo '0 1' | gmx mindist -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -n index.ndx \
            -od mindist.$window.xvg >> $logFile 2>> $errFile 
        check mindist.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

dist_to_Q69(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "dist Q69" 
    if [ ! -f distQ69/distQ69.$window.xvg ]  ; then 
        create_dir distQ69
        cd distQ69

        if [ ! -f gtp.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "resname GTP and name PG" \
                -on gtp.ndx >> $logFile 2>> $errFile 
        fi 
        check gtp.ndx 

        if [ ! -f q69.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "resindex $Q61 and not name \"H*\"" \
                -on q69.ndx >> $logFile 2>> $errFile 
        fi
        check q69.ndx 

        cat gtp.ndx > index.ndx 
        cat q69.ndx >> index.ndx 

        echo '0 1' | gmx mindist -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -n index.ndx \
            -od distQ69.$window.xvg >> $logFile 2>> $errFile 
        check distQ69.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

dist_to_Y32(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "dist Y32" 
    if [ ! -f distY32/distY32.$window.xvg ]  ; then 
        create_dir distY32
        cd distY32

        if [ ! -f Q61.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"SideChain-H\" and resindex $Q61" \
                -on Q61.ndx >> $logFile 2>> $errFile 
        fi 
        check Q61.ndx 

        if [ ! -f Y32.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "resindex 136 and name OH" \
                -on Y32.ndx >> $logFile 2>> $errFile 
        fi
        check Y32.ndx 

        cat Q61.ndx > index.ndx 
        cat Y32.ndx >> index.ndx 

        echo '0 1' | gmx mindist -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -n index.ndx \
            -od distY32.$window.xvg >> $logFile 2>> $errFile 
        check distY32.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

rmsd_gtp(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "rmsd_gtp" 
    if [[ ! -f rmsd_gtp/rmsd.$window.xvg || ! -f rmsd_gtp/active.$window.xvg ]]  ; then 
        create_dir rmsd_gtp
        cd rmsd_gtp

        if [ ! -f gtp.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'resname GTP and not name "H*"' \
                -on gtp.ndx >> $logFile 2>> $errFile 
        fi 
        check gtp.ndx 

        if [ ! -f rmsd.$window.xvg ] ; then
            gmx rms -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -n gtp.ndx \
                -o rmsd.$window.xvg >> $logFile 2>> $errFile 
        fi
        check rmsd.$window.xvg 

        if [ ! -f active.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select '(resname GTP and name PG O1G O2G O3G OS3 PB O1B O2B OS3 PA O1A O2A OS1) or resname MG' \
                -on active.ndx >> $logFile 2>> $errFile 
        fi 
        check active.ndx 

        if [ ! -f active.$window.xvg ] ; then 
            gmx rms -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -n active.ndx \
                -o active.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check active.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

rmsd_l4loop(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s......................." "rmsd_l4loop" 
    if [[ ! -f rmsd_l4loop/rmsd.$window.xvg ]]  ; then 
        create_dir rmsd_l4loop
        cd rmsd_l4loop

        if [ ! -f l4loop.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"Backbone\" and residue 163 to 169" \
                -on l4loop.ndx >> $logFile 2>> $errFile 
        fi 
        check l4loop.ndx 

        if [ ! -f rmsd.$window.xvg ] ; then
            gmx rms -f ../../Production/$window/$MOLEC.$window.xtc \
                -s ../../Production/$window/$MOLEC.$window.tpr \
                -n l4loop.ndx \
                -o rmsd.$window.xvg >> $logFile 2>> $errFile 
        fi
        check rmsd.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

gmx_hbond(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "gmx_hbond"
    if [ ! -f gmx_hbond/hbnum.$window.xvg ]  ; then 
        create_dir gmx_hbond
        cd gmx_hbond

        if [ ! -f acceptors.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'resname GTP and name O1G O2G O3G' \
                -on acceptors.ndx >> $logFile 2>> $errFile 
        fi 
        check acceptors.ndx 

        if [ ! -f water.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'Water' \
                -on water.ndx >> $logFile 2>> $errFile 
        fi 
        check water.ndx 

        if [ ! -f center.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'resname GTP and name PB' \
                -on center.ndx >> $logFile 2>> $errFile 
        fi 
        check center.ndx 

        cat acceptors.ndx > index.ndx 
        cat water.ndx >> index.ndx 
        cat center.ndx >> index.ndx 

        echo '0 1 2' | gmx hbond -f ../../Production/$window/$MOLEC.$window.xtc \
            -s ../../Production/$window/$MOLEC.$window.tpr \
            -n index.ndx \
            -shell 0.7 \
            -num hbnum.$window.xvg >> $logFile 2>> $errFile 
        check hbnum.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

volume_Q61(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "Volume Q61" 
    if [ ! -f volume_Q61/volume.$window.xvg ]  ; then 
        create_dir volume_Q61
        cd volume_Q61

        #analysis function that produces endFile
        if [ ! -f q61.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"SideChain\" and residue $Q61" \
                -on q61.ndx >> $logFile 2>> $errFile 
        fi 

        echo '0' | gmx sasa -s ../../Production/$window/$MOLEC.$window.tpr \
            -f ../../Production/$window/$MOLEC.$window.xtc \
            -n q61.ndx \
            -ndots 240 \
            -probe 0 \
            -tv volume.$window.xvg >> $logFile 2>> $errFile 

        check volume.$window.xvg
        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

hbond_Q61(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "hbond_Q61" 
    if [ ! -f hbond_Q61/hbnum_Q61.$window.xvg ]  ; then 
        create_dir hbond_Q61
        cd hbond_Q61

        #analysis function that produces endFile    
        if [ ! -f index.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'Water' \
                -on water.ndx >> $logFile 2>> $errFile 
            check water.ndx 

            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"SideChain\" and residue $Q61" \
                -on q61.ndx >> $logFile 2>> $errFile 
            check q61.ndx 

            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"C-alpha\" and residue $Q61" \
                -on center.ndx >> $logFile 2>> $errFile 
            check center.ndx 

            cat water.ndx > index.ndx 
            cat q61.ndx >> index.ndx 
            cat center.ndx >> index.ndx 
        fi 

        echo '0 1 2' | gmx hbond -s ../../Production/$window/$MOLEC.$window.tpr \
            -f ../../Production/$window/$MOLEC.$window.xtc \
            -n index.ndx \
            -shell 3 \
            -num hbnum_Q61.$window.xvg >> $logFile 2>> $errFile 
        check hbnum_Q61.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

hbond_Y32(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "hbond_Y32" 
    if [ ! -f hbond_Y32/hbnum_Y32.$window.xvg ]  ; then 
        create_dir hbond_Y32
        cd hbond_Y32

        #analysis function that produces endFile    
        if [ ! -f index.ndx ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select 'Water' \
                -on water.ndx >> $logFile 2>> $errFile 
            check water.ndx 

            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "group \"SideChain\" and residue 136" \
                -on y32.ndx >> $logFile 2>> $errFile 
            check y32.ndx 

            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -select "residue 136 and name OH" \
                -on center.ndx >> $logFile 2>> $errFile 
            check center.ndx 

            cat water.ndx > index.ndx 
            cat y32.ndx >> index.ndx 
            cat center.ndx >> index.ndx 
        fi 

        echo '0 1 2' | gmx hbond -s ../../Production/$window/$MOLEC.$window.tpr \
            -f ../../Production/$window/$MOLEC.$window.xtc \
            -n index.ndx \
            -shell 3 \
            -num hbnum_Y32.$window.xvg >> $logFile 2>> $errFile 
        check hbnum_Y32.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

chi1_Y32(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "Chi1 Y32" 
    if [ ! -f chi1_Y32/angaver_Y32.$window.xvg ]  ; then 
        create_dir chi1_Y32
        cd chi1_Y32

        #analysis function that produces endFile
        if [ ! -f index.ndx ] ; then 
            N=`grep "136TYR" ../../Production/$window/$MOLEC.$window.gro | grep " N " | awk '{print $3}'`
            CA=`grep "136TYR" ../../Production/$window/$MOLEC.$window.gro | grep " CA " | awk '{print $3}'`
            CB=`grep "136TYR" ../../Production/$window/$MOLEC.$window.gro | grep " CB " | awk '{print $3}'`
            CG=`grep "136TYR" ../../Production/$window/$MOLEC.$window.gro | grep " CG " | awk '{print $3}'`
            echo "[ chi1_Y32 ]" > index.ndx 
            echo "$N $CA $CB $CG" >> index.ndx 
        fi 

        gmx angle -f ../../Production/$window/$MOLEC.$window.xtc \
            -n index.ndx \
            -type dihedral \
            -od angdist_Y32.$window.xvg \
            -ov angaver_Y32.$window.xvg >> $logFile 2>> $errFile  
        check angaver_Y32.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount" 
    if [[ ! -f scount/size_245.$window.xvg || ! -f scount/size_3.$window.xvg || ! -f scount/size_4.$window.xvg || ! -f scount/size_5.$window.xvg || scount/size_6.$window.xvg ]]  ; then 
        create_dir scount
        cd scount

        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.245 of group \"SideChain\" and resindex $Q61) and (within 0.245 of resname GTP and name O1G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.3 of group \"SideChain\" and resindex $Q61) and (within 0.3 of resname GTP and name O1G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.4 of group \"SideChain\" and resindex $Q61) and (within 0.4 of resname GTP and name O1G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        if [ ! -f size_6.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.6 of group \"SideChain\" and resindex $Q61) and (within 0.6 of resname GTP and name O1G))" \
                -os size_6.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_6.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount_O2G(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount O2G" 
    if [[ ! -f scount_O2G/size_245.$window.xvg || ! -f scount_O2G/size_3.$window.xvg || ! -f scount_O2G/size_4.$window.xvg || ! -f scount_O2G/size_5.$window.xvg || scount_O2G/size_6.$window.xvg ]]  ; then 
        create_dir scount_O2G
        cd scount_O2G

        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.245 of group \"SideChain\" and resindex $Q61) and (within 0.245 of resname GTP and name O2G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.3 of group \"SideChain\" and resindex $Q61) and (within 0.3 of resname GTP and name O2G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.4 of group \"SideChain\" and resindex $Q61) and (within 0.4 of resname GTP and name O2G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O2G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        if [ ! -f size_6.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.6 of group \"SideChain\" and resindex $Q61) and (within 0.6 of resname GTP and name O2G))" \
                -os size_6.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_6.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount_vary_O1G(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount O1G" 
    if [[ ! -f scount_O1G/size_245.$window.xvg || ! -f scount_O1G/size_3.$window.xvg || ! -f scount_O1G/size_4.$window.xvg || ! -f scount_O1G/size_5.$window.xvg || scount_O1G/size_6.$window.xvg ]]  ; then 
        create_dir scount_O1G
        cd scount_O1G

        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.245 of resname GTP and name O1G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.3 of resname GTP and name O1G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.4 of resname GTP and name O1G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        if [ ! -f size_6.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.6 of resname GTP and name O1G))" \
                -os size_6.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_6.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount_both(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount_both" 
    if [[ ! -f scount_both/size_245.$window.xvg || ! -f scount_both/size_3.$window.xvg || ! -f scount_both/size_4.$window.xvg || ! -f scount_both/size_5.$window.xvg ]]  ; then 
        create_dir scount_both
        cd scount_both

        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.245 of group \"SideChain\" and resindex $Q61) and (within 0.245 of resname GTP and name O1G O2G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.3 of group \"SideChain\" and resindex $Q61) and (within 0.3 of resname GTP and name O1G O2G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.4 of group \"SideChain\" and resindex $Q61) and (within 0.4 of resname GTP and name O1G O2G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G O2G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount_vary_SCD(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount SCD" 
    if [[ ! -f scount_SCD/size_245.$window.xvg || ! -f scount_SCD/size_3.$window.xvg || ! -f scount_SCD/size_4.$window.xvg || ! -f scount_SCD/size_5.$window.xvg || scount_SCD/size_6.$window.xvg ]]  ; then 
        create_dir scount_SCD
        cd scount_SCD

        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.245 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.3 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.4 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        if [ ! -f size_6.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.6 of group \"SideChain\" and resindex $Q61) and (within 0.5 of resname GTP and name O1G))" \
                -os size_6.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_6.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

scount_polar(){
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
        X)  ##rtp HIP; charged Histidine  
            polarAtoms="resindex $Q61 and name ND1 NE2 HE2 HD1"
          ;; 
        Y) 
            polarAtoms="resindex $Q61 and name OH HH" 
          ;; 
        *) 
            sideChainPolar=false
            polarAtoms="resindex $Q61 and name"
          ;; 
    esac 


    if [ "$sideChainPolar" != false ] ; then 
    printf "\t\t\t%10s........................" "scount_polar" 
        if [[ ! -f scount_polar/size_245.$window.xvg || ! -f scount_polar/size_3.$window.xvg || ! -f scount_polar/size_4.$window.xvg || ! -f scount_polar/size_5.$window.xvg ]]  ; then 
        create_dir scount_polar
        cd scount_polar

        #analysis function that produces endFile
        if [ ! -f size_245.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.245 of $polarAtoms) and (within 0.245 of resname GTP and name O1G))" \
                -os size_245.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_245.$window.xvg

        if [ ! -f size_3.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.3 of $polarAtoms) and (within 0.3 of resname GTP and name O1G))" \
                -os size_3.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_3.$window.xvg

        if [ ! -f size_4.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.4 of $polarAtoms) and (within 0.4 of resname GTP and name O1G))" \
                -os size_4.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_4.$window.xvg

        if [ ! -f size_5.$window.xvg ] ; then 
            gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -select "group \"Water\" and same residue as ((within 0.5 of $polarAtoms) and (within 0.5 of resname GTP and name O1G))" \
                -os size_5.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check size_5.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
    fi 
}

scount_2(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "scount_2" 
    if [ ! -f scount_2/size.$window.xvg ]  ; then 
        create_dir scount_2
        cd scount_2

        #analysis function that produces endFile
        gmx select -s ../../Production/$window/$MOLEC.$window.tpr \
            -f ../../Production/$window/$MOLEC.$window.xtc \
            -select "group \"Water\" and same residue as ((within 0.4 of group \"Backbone\" and resindex 164) and (within 0.3 of resname GTP and name O1G))" \
            -os size.$window.xvg >> $logFile 2>> $errFile 
        check size.$window.xvg

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
}

force_nitrile(){
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\n\t\t\tCalculating force:\n"
    if [[ ! -f force_nitrile/$MOLEC.$window.solvent_rxn_field.projected.xvg || ! -f force_nitrile/$MOLEC.$window.external_field.projected.xvg || ! -f force_nitrile/$MOLEC.$window.total_field.projected.xvg ]] ; then
        create_dir force_nitrile
        cp ../Production/$window/*.itp force_nitrile/. 
        cp ../Prep/Solvate/neutral.gro force_nitrile/. 
        cp ../Production/$window/$MOLEC.$window.top force_nitrile/. 

        cd force_nitrile

        #analysis function that produces endFile
        ## We use version 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
        ## We allow for warnings, since we are generated .tpr from a gromacs 2016 mdp file. We are only inserting
        ## atoms this should not matter.
        grompp -f $MDP/vac_md.mdp \
            -c neutral.gro \
            -p $MOLEC.$window.top \
            -maxwarn 3 \
            -o v4.$window.tpr >> $logFile 2>> $errFile
        check v4.$window.tpr

        CD=`grep CNC neutral.gro | grep CD | awk '{print $3}'`
        NE=`grep CNC neutral.gro | grep NE | awk '{print $3}'`

        printf "\t\t\t\tInserting dummy atoms....."
        if [ ! -s $MOLEC.$window.with_dummy.top ] ; then 
            if [ ! -s $MOLEC.$window.with_dummy.xtc ] ; then
                source /usr/local/gromacs/bin/GMXRC
                $FORCE_TOOLS/g_insert_dummy_atom -f ../nopbc/nopbc.$window.xtc \
                    -s v4.$window.tpr \
                    -a1 $CD \
                    -a2 $NE \
                    -o $MOLEC.$window.with_dummy.xtc >> $logFile 2>> $errFile
               fi 
            check $MOLEC.$window.with_dummy.xtc

            ## We use the initial configuration so that titration states are conserved (ie, at the end of the production run, pdb2gmx might assign a different titration state to a histidine, which causes it to fail.
            if [ ! -s $MOLEC.$window.with_dummy.pdb ] ; then
                echo 'Protein System' | gmx trjconv -f neutral.gro \
                    -s v4.$window.tpr \
                    -center \
                    -ur compact \
                    -pbc mol \
                    -o $MOLEC.$window.nopbc.gro >> $logFile 2>> $errFile
                check $MOLEC.$window.nopbc.gro

                source /usr/local/gromacs/bin/GMXRC
                $FORCE_TOOLS/g_insert_dummy_atom -f $MOLEC.$window.nopbc.gro \
                    -s v4.$window.tpr \
                    -a1 $CD \
                    -a2 $NE \
                    -o $MOLEC.$window.with_dummy.gro >> $logFile 2>> $errFile
            
                gmx editconf -f $MOLEC.$window.with_dummy.gro \
                    -o $MOLEC.$window.with_dummy.pdb >> $logFile 2>> $errFile


                fi
            check $MOLEC.$window.with_dummy.pdb 
            ##Annoyingly, editconf puts the Na atoms behind the SOL molecules, which causes probles for pdb2gmx -merge all
            ##   We have to put the Na atoms directly after the GTP molecule. 
            cat $MOLEC.$window.with_dummy.pdb | grep -v NA > no_na.pdb 
            headlines=`cat -n $MOLEC.$window.with_dummy.pdb | grep OC2 | awk '{print $1}' | tail -n1`
            totalLines=`wc -l no_na.pdb | awk '{print $1}'`
            taillines=`echo "$totalLines - $headlines" | bc -l`

            head -n $headlines no_na.pdb > temp.pdb 
            cat $MOLEC.$window.with_dummy.pdb | grep NA >> temp.pdb 
            tail -n $taillines no_na.pdb >> temp.pdb 
            mv temp.pdb $MOLEC.$window.with_dummy.pdb 

            ## We have to use a PDB with TER breaks so that the termini get assigned correctly.
            awk '/OC2/{print;print "TER";next}1' $MOLEC.$window.with_dummy.pdb > temp.pdb 
            mv temp.pdb $MOLEC.$window.with_dummy.pdb

            ##All His residues are rtp HIE in these simulations (luckily) 
            if [ ${MOLEC: -1} == "H" ] ; then 
                echo '1 1 1 1 1 1' | gmx pdb2gmx -f $MOLEC.$window.with_dummy.pdb \
                    -his \
                    -water tip3p \
                    -ff amber03 \
                    -merge all \
                    -p $MOLEC.$window.with_dummy.top \
                    -o $MOLEC.$window.with_dummy.gro >> $logFile 2>> $errFile
            else 
                echo '1 1 1 1 1' | gmx pdb2gmx -f $MOLEC.$window.with_dummy.pdb \
                    -his \
                    -water tip3p \
                    -ff amber03 \
                    -merge all \
                    -p $MOLEC.$window.with_dummy.top \
                    -o $MOLEC.$window.with_dummy.gro >> $logFile 2>> $errFile
            fi 
            check $MOLEC.$window.with_dummy.top
            printf "Done\n"
        else 
            printf "Skipped\n" 
            fi 

        ##Find new atom numbers
        CD=`grep CNC $MOLEC.$window.with_dummy.gro | grep CD | awk '{print $3}'`
        NE=`grep CNC $MOLEC.$window.with_dummy.gro | grep NE | awk '{print $3}'`

        echo "[ probe ]" > probe.ndx
        echo "$CD $NE" >> probe.ndx

        echo "[ protein ]" > protein.ndx
        grep -v TCHG $MOLEC.$window.with_dummy.gro | grep -v SOL | grep -v HOH | grep -v NA | grep -v CL | grep -v MG | grep -v GTP | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx

        cp $MOLEC.$window.with_dummy.top $MOLEC.$window.total_field.top

        if [ ! -s $MOLEC.$window.solvent_rxn_field.top ] ; then
            $FORCE_TOOLS/zero_charges.py $MOLEC.$window.with_dummy.top protein.ndx $MOLEC.$window.solvent_rxn_field.top >> $logFile 2>> $errFile
            fi

        if [ ! -s $MOLEC.$window.external_field.top ] ; then
            $FORCE_TOOLS/zero_charges.py $MOLEC.$window.with_dummy.top probe.ndx $MOLEC.$window.external_field.top >> $logFile 2>> $errFile
            fi
        check $MOLEC.$window.total_field.top $MOLEC.$window.external_field.top $MOLEC.$window.solvent_rxn_field.top

        for field in total_field external_field solvent_rxn_field ; do 
            printf "\t\t%20s..." $field 

            ##Extract forces 
            if [ ! -s $MOLEC.$window.$field.projected.xvg ] ; then 
                printf "forces..." 
                if [ ! -s $MOLEC.$window.$field.xvg ] ; then 
                    if [ ! -s $MOLEC.$window.$field.tpr ] ; then 
                        gmx grompp -f $MDP/rerun.mdp \
                            -p $MOLEC.$window.$field.top \
                            -c $MOLEC.$window.with_dummy.gro \
                            -o $MOLEC.$window.$field.tpr  >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$window.$field.tpr 
 
                    if [ ! -s $MOLEC.$window.$field.trr ] ; then 
                        gmx mdrun -s $MOLEC.$window.$field.tpr \
                            -nt 1 \
                            -rerun $MOLEC.$window.with_dummy.xtc \
                            -deffnm $MOLEC.$window.$field >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$window.$field.trr 

                    echo 2 | gmx traj -f $MOLEC.$window.$field.trr \
                        -s $MOLEC.$window.$field.tpr \
                        -xvg none \
                        -of $MOLEC.$window.$field.xvg >> $logFile 2>> $errFile 
                    rm $MOLEC.$window.$field.trr 
                fi 
                check $MOLEC.$window.$field.xvg 

                ##extract postions for bond vector
                printf "positions..." 
                if [ ! -s $MOLEC.$window.positions.xvg ] ; then 
                    gmx traj -f $MOLEC.$window.with_dummy.xtc \
                        -s $MOLEC.$window.$field.tpr \
                        -n probe.ndx \
                        -xvg none \
                        -ox $MOLEC.$window.positions.xvg >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.$window.positions.xvg 

                ##project force along bond vector 
                printf "Projecting..." 
                $FORCE_TOOLS/get_force.py $MOLEC.$window.positions.xvg $MOLEC.$window.$field.xvg $MOLEC.$window.$field.projected.xvg 
                check $MOLEC.$window.$field.projected.xvg 
                printf "Done\n"  
            else 
                printf "...................Skipped\n" 
                fi 
            done 

        check $MOLEC.$window.total_field.projected.xvg $MOLEC.$window.external_field.projected.xvg $MOLEC.$window.solvent_rxn_field.projected.xvg
        clean 
        cd ../
    else 
        printf "\t\t\t\t  ........................Skipped\n"
        fi 
}

g12_os3_dist(){
    ##Call from analysis_windows()
    if [ -z $1 ] ; then 
        echo "ERROR: Argument missing." 
        echo "Usage: $0 < window (degrees) > "
    fi 
    window=$1

    printf "\t\t\t%10s........................" "G12 OS3 Dist" 
    if [[ ! -f g12_os3_dist/hbnum.$window.xvg || ! -f distave.$window.xvg ]]  ; then 
        create_dir g12_os3_dist
        cd g12_os3_dist

        #analysis function that produces endFile

        if [ ! -f hbnum.$window.xvg ] ; then 
            echo "[ GTP ]" > index.ndx 
            cat ../../Production/$window/$MOLEC.$window.gro | grep "GTP" | grep "OS3" | awk '{print $3}' >> index.ndx 
            echo "[ G12 ]" >> index.ndx 
            cat ../../Production/$window/$MOLEC.$window.gro | grep "117GLY" | grep " N " | awk '{print $3}' >> index.ndx 
            cat ../../Production/$window/$MOLEC.$window.gro | grep "117GLY" | grep " H " | awk '{print $3}' >> index.ndx 

            echo '1 0' | gmx hbond -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -n index.ndx \
                -num hbnum.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check hbnum.$window.xvg 

        if [ ! -f distave.$window.xvg ] ; then 
            echo "[ OS3_G12 ]" > index.ndx 
            cat ../../Production/$window/$MOLEC.$window.gro | grep "GTP" | grep "OS3" | awk '{print $3}' >> index.ndx 
            cat ../../Production/$window/$MOLEC.$window.gro | grep "117GLY" | grep " N " | awk '{print $3}' >> index.ndx 

            echo '0' | gmx distance -s ../../Production/$window/$MOLEC.$window.tpr \
                -f ../../Production/$window/$MOLEC.$window.xtc \
                -n index.ndx \
                -oav distave.$window.xvg >> $logFile 2>> $errFile 
        fi 
        check distave.$window.xvg 

        clean 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n"
        fi 
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
