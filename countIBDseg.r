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
setwd("/Users/cmoreau/Documents/Kyoto/results/IBD/")
getwd()

args <- commandArgs(trailingOnly = TRUE)
print (paste(length(args), " arguments"))
#If called with 4 arguments:
#Do 1 plot per page
if (length(args)>1){
	ichr <- as.numeric(args[1]) #chr
	segnum <- as.numeric(args[2]) #num of IBD segments
	#this file is created in python -see /home/claudia/Documents/Achilee/Users/cmoreau/Documents/Kyoto/scripts/run_IBDcount.sh
	filename<- (args[3]) #input file name - 
	sufix<- as.numeric(args[4]) #filtered...or not...

	print (ichr)
	print (segnum)
	print (filename)
	print (sufix)


	#Take file of number of pair of inds who share at each position
	sharing<-read.table(filename,sep="\t",header=F,colClasses=c("character","numeric" ))
	#sharing[1:10,]
	myfile<-paste("chr",ichr,"_",sufix,".jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	plot(x=1:length(sharing[,1]),y=sharing[,2],type="l",main=paste("Chr",ichr,"\n",segnum," IBD segments",sep=""),ylab="Count",xlab="Markers")
	dev.off()

#If called with 1 argument:
#Do all plots in the same 
}else{

	sufix<- (args[1]) #filtering parameter
	print (sufix)

	print("Making figure of IBD sharing by position for all chr...")

	#Takes file of number of pair of inds who share at each position
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBDsharing_positions_allchr_",sufix,".jpg",sep="")
	jpeg(filename=myfile,width=14,height=11,units="in",res=600)
	par(mfrow=c(4,6))

	for (i in 1:22){
		
		filename<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBDsharing_positions_chr",i,"_",sufix,".txt",sep="")
		if (!file.exists (filename)) {
			filename<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBDsharing_positions_chr",i,"_",sufix,"_bp.txt",sep="")
		}

		sharing<-read.table(filename,sep="\t",header=F,colClasses=c("character","numeric" ))
		plot(x=1:length(sharing[,1]),y=sharing[,2],type="l",main=paste("Chr",i,sep=""),ylab="Count",xlab="Markers")
	}
	
	dev.off()



	## Calculate total length of IBD for each pair of inds

	print ("Calculating total IBD sharing length and segment num distribution...")

	#Takes file of IBD segment sharing
		
	filename<-paste("/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/",sufix,"_0.1_3.0/chr1.bp2cM.IBD_filtered.match",sep="")
	IBD<-read.table(filename,sep="\t",header=F,colClasses=c(rep("character",4),rep("numeric",4),rep("character",2),rep("numeric",2),"character",rep("numeric",3) ))
	dimnames(IBD)[[2]]<-c("ID1","ID1chr","ID2","ID2chr","chr","startcM","stopcM","idontknow","startsnp","stopsnp","snpnum","lengthcM","unit","toto1","toto2","toto3")

	#filename<-paste("/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/chr1.bp2cM.IBD.match",sep="")
	#IBD<-read.table(filename,sep="\t",header=F,colClasses=c(rep("character",4),rep("numeric",3),rep("character",2),rep("numeric",2),"character",rep("numeric",3) ))
	#dimnames(IBD)[[2]]<-c("ID1","ID1chr","ID2","ID2chr","chr","startcM","stopcM","startsnp","stopsnp","snpnum","lengthcM","unit","toto1","toto2","toto3")
	
	IBD$ID12<-paste(IBD$ID1,IBD$ID2,sep="-")

	print (1)
	print (length(IBD[,1]))
	print (IBD[which(IBD$lengthcM==max(IBD$lengthcM)),])

	## Total length
	lengthdf <- as.data.frame(xtabs(lengthcM ~ ID12, data = IBD))

	## Segment no
	segnumdf <- as.data.frame(table(IBD$ID12))
	
	for (i in 2:22){
		print (i)
		filename<-paste("/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/",sufix,"_0.1_3.0/chr",i,".bp2cM.IBD_filtered.match",sep="")
		IBD<-read.table(filename,sep="\t",header=F,colClasses=c(rep("character",4),rep("numeric",4),rep("character",2),rep("numeric",2),"character",rep("numeric",3) ))
		dimnames(IBD)[[2]]<-c("ID1","ID1chr","ID2","ID2chr","chr","startcM","stopcM","idontknow","startsnp","stopsnp","snpnum","lengthcM","unit","toto1","toto2","toto3")
		#filename<-paste("/Volumes/project/gravel/moreau_projects/Kyoto/data/phased/3_GERMLINE/chr",i,".IBD.match",sep="")
		#IBD<-read.table(filename,sep="\t",header=F,colClasses=c(rep("character",4),rep("numeric",3),rep("character",2),rep("numeric",2),"character",rep("numeric",3) ))
		#dimnames(IBD)[[2]]<-c("ID1","ID1chr","ID2","ID2chr","chr","startcM","stopcM","startsnp","stopsnp","snpnum","lengthcM","unit","toto1","toto2","toto3")

		print (length(IBD[,1]))
		print (IBD[which(IBD$lengthcM==max(IBD$lengthcM)),])

		IBD$ID12<-paste(IBD$ID1,IBD$ID2,sep="-")

		## Total length
		tmplength <- as.data.frame(xtabs(lengthcM ~ ID12, data = IBD))
		lengthdf <- merge(lengthdf,tmplength,by=1)

		## Segment no
		tmpsegnum <- as.data.frame(table(IBD$ID12))
		segnumdf <- merge(segnumdf,tmpsegnum,by=1)
		
	}

	## Total length
	lengthdf$total <- apply(data=lengthdf[,-1],1,sum)

	## Segment no
	segnumdf$total <- apply(data=segnumdf[,-1],1,sum)
	
	tmpcordf <- merge(lengthdf[,c(1,24)],segnumdf[,c(1,24)],by=1)
	print (tmpcordf[1:10,])



	print ("Making plots...")
	#Plot correlation between segment no and total legnth
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_total_length_seg_num_correlation_",sufix,".jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(tmpcordf[,1])
	#Adding 0s for pairs that don't share
	plot(tmpcordf[,2],tmpcordf[,3],main=paste("Correlation between total IBD length and no. of IBD segments shared\n",pairnum," individuals' pairs (",sufix,")\nPearson r=",format(cor(tmpcordf[,2],tmpcordf[,3]),digits=2),sep=""),xlab="Length (cM)",ylab="Segment no.")
	dev.off()



	#Plot distribution for all segments
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_total_length_distribution_",sufix,".jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data - take in the plinf.fam file
	filename<-"/Volumes/project/gravel/moreau_projects/Kyoto/data/unphased/Kyoto_unrelated_notinbred_nohh_250Ksnps.fam"
	numindsdf<-read.table(filename,sep="",header=F,colClasses=c(rep("character",2),rep("numeric",4) ))
	numinds <- length(numindsdf[,1])
	print(paste(numinds," inds",sep="") )
	
	pairnum<-(numinds*(numinds-1)/2)
	#Adding 0s for pairs that don't share
	hist(c(lengthdf[,2],rep(0,pairnum-length(lengthdf[,1]))),breaks=100,main=paste("Total IBD length for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Length (cM)",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 3 cM
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_total_length_distribution_",sufix,"_3cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(lengthdf[,2])
	hist(lengthdf[,2],breaks=100,main=paste("Total IBD length (>=3cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Length (cM)",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 8 cM
	lengthdf <- as.data.frame(xtabs(lengthcM ~ ID12, data = IBD[IBD$lengthcM>=8,]))

	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_total_length_distribution_",sufix,"_8cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(lengthdf[,2])
	hist(lengthdf[,2],breaks=100,main=paste("Total IBD length (>=8cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Length (cM)",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 100 cM
	lengthdf <- as.data.frame(xtabs(lengthcM ~ ID12, data = IBD[IBD$lengthcM>=18,]))

	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_total_length_distribution_",sufix,"_18cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(lengthdf[,2])
	hist(lengthdf[,2],breaks=100,main=paste("Total IBD length (>=18cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Length (cM)",ylab="Counts")
	dev.off()




	## Calculate no of IBD segments shared for each pair of inds

	#Plot distribution for all segments
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_segment_num_distribution_",sufix,".jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-(numinds*(numinds-1)/2)
	#Adding 0s for pairs that don't share
	hist(c(segnumdf[,2],rep(0,pairnum-length(segnumdf[,1]))),breaks=100,main=paste("IBD segment no. for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Segment no.",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 3 cM
	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_segment_num_distribution_",sufix,"_3cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(segnumdf[,2])
	hist(segnumdf[,2],breaks=100,main=paste("IBD segment no. (>=3cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Segment no.",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 8 cM
	segnumdf <- as.data.frame(table(IBD$ID12[IBD$lengthcM>=8]))

	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_segment_num_distribution_",sufix,"_8cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(segnumdf[,2])
	hist(segnumdf[,2],breaks=100,main=paste("IBD segment no. (>=8cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Segment no.",ylab="Counts")
	dev.off()

	#Plot distribution for segments > 18 cM
	segnumdf <- as.data.frame(table(IBD$ID12[IBD$lengthcM>=18]))

	myfile<-paste("/Users/cmoreau/Documents/Kyoto/results/IBD/IBD_segment_num_distribution_",sufix,"_18cM.jpg",sep="")
	jpeg(filename=myfile,width=8.5,height=8.5,units="in",res=600)
	#Number of pairs in the whole data
	pairnum<-length(segnumdf[,2])
	hist(segnumdf[,2],breaks=100,main=paste("IBD segment no. (>=18cM) for ",pairnum," individuals' pairs\n",sufix,sep=""),xlab="Segment no.",ylab="Counts")
	dev.off()


}









#Count IBD segments per position - this doesn't work and is now done in python...
#Try with a loop...
#i<-400
#countIBDvec<-c()
#print ("IBD length: ")
#print (length(IBD[,1]))
#for (i in 1:length(IBD[,1])){
#	if (i %% 1000 == 0) {print (i)}
#	countIBDvec <- c(countIBDvec,map[map[,3]>=IBD[i,6] & map[,3]<=IBD[i,7],3])
#}
#length(countIBDvec)

#Look if IBD segs are spread with no significant pics in...
#IBDdf <- as.data.frame(table(countIBDvec ))
#length(IBDdf [,1])

#myfile<-"/Users/cmoreau/Documents/Kyoto/data/chr21_allIBD.jpg"
#jpeg(filename=myfile,width=8.5,height=11,units="in",res=600)
#plot(x=c(1:length(IBDdf [,1])),y=IBDdf[,2],type="l",main=paste("Chr21\n",length(IBD[,1])," IBD segments",sep=""),ylab="Count",xlab="Markers")
#dev.off()

q()
