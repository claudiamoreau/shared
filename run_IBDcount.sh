#! /bin/sh

#Oct 2, 2015

param1=150000
#for i in {1..4}
#do
#	infile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/chr"${i}".IBD.match
#	mapfile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/2_PLINK/cMcorrected/chr"${i}".cM.map
#	outfile=/Users/cmoreau/Documents/Kyoto/results/IBD/IBDsharing_positions_chr"${i}"_allIBD_bp.txt
#	segnum=$(wc -l < "${infile}")

#	python /Users/cmoreau/Documents/Kyoto/scripts/IBDsharing_byPosition.py ${infile} ${mapfile} ${outfile}

#	Rscript /Users/cmoreau/Documents/Kyoto/scripts/countIBDseg.r ${i} ${segnum} ${outfile} ${param1}
#done

#Oct 5, 2015

#for i in {5..22}
#do
#i=3
#	infile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/"${param1}"_0.1_3.0/chr"${i}".bp2cM.IBD_filtered.match
#	mapfile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/2_PLINK/cMcorrected/chr"${i}".cM.map
#	outfile=/Users/cmoreau/Documents/Kyoto/results/IBD/IBDsharing_positions_chr"${i}"_"${param1}".txt
#	segnum=$(wc -l < "${infile}")

#	python  /Users/cmoreau/Documents/Kyoto/scripts/IBDsharing_byPosition.py ${infile} ${mapfile} ${outfile}

#	Rscript /Users/cmoreau/Documents/Kyoto/scripts/countIBDseg.r ${i} ${segnum} ${outfile} ${param1}

#done

#Run Rscript with only 1 argument to put all graphs in the same file
#and to plot total IBD length distribution
Rscript /Users/cmoreau/Documents/Kyoto/scripts/countIBDseg.r ${param1}
