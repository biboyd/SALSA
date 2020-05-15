#!/bin/bash
refine=$1
ds=$2

# run through the silicons
for i in 1 2 3;
do

	for dist in 100 500; 
	do 
		sbatch -J "rot${dist} Si${i} $ds" runRotProj.sb $refine $ds Si_p${i}_number_density 312 30 $dist
		sbatch -J "rot${dist} Si${i} $ds" runRotProj_cgm.sb $refine $ds Si_p${i}_number_density 312 30 $dist
	done
done

# run throught the carbon
for i in 1 3;
do

    for dist in 100 500; 
    do 
        sbatch -J "rot${dist} C${i} $ds" runRotProj.sb $refine $ds C_p${i}_number_density 312 30 $dist
        sbatch -J "rot${dist} C${i} $ds" runRotProj_cgm.sb $refine $ds C_p${i}_number_density 312 30 $dist
    done
done

#run through everything else
for fld in density temperature metallicity H_p0_number_density O_p5_number_density;
do
    for dist in 100 500; 
    do 
        sbatch -J "rot${dist} $fld $ds" runRotProj.sb $refine $ds $fld 312 30 $dist
        sbatch -J "rot${dist} $fld $ds" runRotProj_cgm.sb $refine $ds $fld 312 30 $dist
    done
done	
