! /bin/sh

#Oct 7, 2015


infile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/150000_0.1_3.0/

mapfile=/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/2_PLINK/cMcorrected/
samplefile=/Volumes/project/gravel/moreau_projects/Kyoto/data/unphased/Kyoto_unrelated_notinbred_nohh_nohet_250Ksnps.fam
outfile=/Users/cmoreau/Documents/Kyoto/results/IBD/clauses_byPosition_allchr_150000.txt
numind=10

python /Users/cmoreau/Documents/Kyoto/scripts/buildFile4panel/IBDsharing_byPosition_4panel.py ${infile} ${mapfile} ${samplefile} ${outfile} ${numind}


