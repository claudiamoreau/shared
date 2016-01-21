#Created by Claudia
#Oct 7, 2015

import sys
from sys import argv
import numpy as np

#Files path must be passed in argument
datafile = sys.argv[1] 
print "Data is in : " + datafile
mapfile = sys.argv[2] 
print "Map is in : " + mapfile
samplefilename = sys.argv[3] 
print "Sample file : " + samplefilename
outputfilename = sys.argv[4] 
print "Output file : " + outputfilename
indseq = int(sys.argv[5]) 
print "No of individuals to output in the panel : " + str(indseq)

#loads IBD file
#Receives IBD file name (GERMLINE formatted)
#Returns list of 4 lists: ID1,ID2,SNP1,SNP2
def loadfile (datafile):
    alldata = np.genfromtxt(fname = datafile,dtype=str)

    #indices = range(len(alldata[:,0]))
    ID1 = map(str, alldata[:,0].tolist())
    #print alldata[:,0]
    #ID1chr = map(str, alldata[:,1].tolist())
    ID2 = map(str, alldata[:,2].tolist())
    #ID2chr = map(str, alldata[:,3].tolist())
    #chr = map(int, alldata[:,4].tolist())
    #cM1 = map(float, alldata[:,5].tolist())
    #cM2 = map(float, alldata[:,6].tolist())
    ##ATTENTION - filtered data have an empty column in 7 so here we should shift everything
    SNP1 = map(str, alldata[:,7].tolist())
    SNP2 = map(str, alldata[:,8].tolist())
    #notsure = map(int, alldata[:,9].tolist())
    #length = map(float, alldata[:,10].tolist())

    datalist = []
    datalist.extend([ID1, ID2, SNP1, SNP2])

    return datalist

#loads map file
#Receives map file name (PLINK cM)
#Returns SNP list
def loadmapfile (mapfile):
    alldata = np.genfromtxt(fname = mapfile, delimiter='\t',dtype=str)

    #indices = range(len(alldata[:,0]))
    #chr = map(int, alldata[:,0].tolist())
    SNP = map(str, alldata[:,1].tolist())
    #cM = map(float, alldata[:,2].tolist())
    #bp = map(int, alldata[:,3].tolist())

    return SNP

#loads sample file
#Receives sample file name (PLINK.fam)
#Returns sample list
def loadsamplefile (samplefile):
    alldata = np.genfromtxt(fname = samplefile,dtype=str)

    indices = range(len(alldata[:,0]))
    sample = map(str, alldata[:,0].tolist())

    return sample

#Writes file
#Receives out file name ,dict of samples sharing and the list of most representative in order
#Writes output file of ind with coverage
def writefile (outputfilename,ind_dict,final_big_list):
    outputfile = open(outputfilename, "w") #To append to file foreach chr
    for i in final_big_list:
        outputfile.write(str(i) + "\t" + str(len(ind_dict[i])) + "\n")

    outputfile.close()


#Looks at the most representative guy...
#Receives dict of ind-segment with value=list of cluases for all IBD segments
#Returns the ind who is the most representative
def findbiggest (ind_dict,final_big_list):
    biggestlength=0
    chosensample=''

    for sample in ind_dict:
        if sample not in final_big_list:
            tmplength=len(ind_dict[sample])
            if tmplength>=biggestlength:
                biggestlength=tmplength
                chosensample=sample

    return (chosensample)

#Removes the clauses from the list of samples
#Receives dict of ind-segment with value=list of cluases for all IBD segments
#Receives the most representative guy
#Returns dict without clauses
def removeclauses (ind_dict,big,final_big_list):
    big_list = ind_dict[big]

    for sample in big_list:
        ## I don't want to remove clauses for samples already chosen as most representative
        if sample not in final_big_list:
            #print sample
            #print (len(ind_dict[sample]))
            #print "Removing clauses..."
            for clause in big_list:
                ## I only want to remove one occurence of a clause in case
                ## there are more than one as indivuals are considered as haploid chromosomes
                if clause in ind_dict[sample]:
                    ind_dict[sample].remove(clause)
            if big in ind_dict[sample]:
                ind_dict[sample].remove(big)
            #print (len(ind_dict[sample]))

    return (ind_dict)



########## Now start doing something! ##########
print 'loading sample file...'
samples = loadsamplefile(samplefilename)
print str(len(samples)) + " samples"

## Build dict of inds with value = list of inds who share with this guy
ind_dict = dict()

#For all chr...
for j in range(1,22):

    #Loading files
    print 'loading data file...'
    data = loadfile(datafile + "chr" + str(j) + ".bp2cM.IBD_filtered.match")
    print str(len(data[0])) + " segments"
    print 'loading map file...'
    positions = loadmapfile(mapfile + "chr" + str(j) + ".cM.map")
    print str(len(positions)) + " positions"



    ## fill dict with all samples for all segments in IBD sharing data
    for i in range(len(data[0])):
    #for i in range(10):

        #print range(positions.index(data[2][i]),(positions.index(data[3][i])+1),1)
        ## OK finally we decided that is it better to take the whole segment instead of each position
        #for j in range(positions.index(data[2][i]),(positions.index(data[3][i])+1),1):
        #print data[2][i]
        if data[0][i] not in ind_dict:
            ind_dict[data[0][i]] = [data[1][i]]
        else:
            ind_dict[data[0][i]].append(data[1][i])
        if data[1][i] not in ind_dict:
            ind_dict[data[1][i]] = [data[0][i]]
        else:
            ind_dict[data[1][i]].append(data[0][i])


## List of most representative samples in order
final_big_list = list()
print 'Finding most representative samples...'
## Take the one with the biggest no of relateds 
## Take the second biggest and so on...

for i in range(indseq):
    if i%10==0:
        print str(i) + ' done on ' + str(indseq)
        
    biggest = findbiggest(ind_dict,final_big_list)
    if biggest == '':
        break
    final_big_list.append(biggest)
    #print final_big_list

    ## Remove from the list all those that were already accounted for
    ind_dict = removeclauses(ind_dict,biggest,final_big_list)


print 'Writing output file...'
writefile (outputfilename,ind_dict,final_big_list)

print "Finished!!"
