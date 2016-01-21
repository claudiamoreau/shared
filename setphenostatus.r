################################
#
#Claudia Moreau
#
#14 avril 2015
#
#Called with 3 arguments: Plot IBD segments at each position for each chr
#Caleed with 1 argument: Plot IBD segments for all chr and 
#plot total IBD sharing length between pairs of individudals
#
################################

#This script can be called with /home/claudia/Documents/Achilee/Users/cmoreau/Documents/Kyoto/scripts/run_IBDcount.sh

rm(list = ls());
setwd("/Volumes/project/gravel/moreau_projects/Kyoto/data/")
getwd()

#Take fam file
filename<-"unphased/Kyoto_unrelated_notinbred_nohh_nohet_250Ksnps.fam"
fam<-read.table(filename,sep=" ",header=F,colClasses=c(rep("character",2),rep("numeric",4) ))
fam[1:10,]
length(fam[,1])

#Set as 2-affected one st of data at a time...
filename<-"rawdata/Kyoto_genotype/extractSNPs/set1_610.fam"
fam1<-read.table(filename,sep=" ",header=F,colClasses=c(rep("character",2),rep("numeric",4) ))
fam1[1:10,]
length(fam1[,1])

fam[,6] <- rep(1,length(fam[,1]))
fam[which(is.element(fam[,1],fam1[,1])),6] <- 2

print (length(fam[fam[,6]==2,1]))
print (length(fam1[,1]))

myfile<-"unphased/Kyoto_affected_set1_610.fam"
write.table(fam,myfile,quote=F,row.names=F,col.names=F)

#set1_660.fam, set1_2500e.fam,set1_2500k.fam, set1_2500o.fam, set1_5000e.fam
#set2_2500o.fam, set2_610.fam,
#set3_610.fam, set3_2500k.fam,set3_2500o.fam,set3_5000e.fam


q()
