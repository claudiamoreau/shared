#Created by Claudia
#Jan 28 2015

import sys
from sys import argv
import numpy as np

#Filename and path in argument and output filename
datafilename = sys.argv[1] #0 is python script!
print "Data is in : " + datafilename
mapfilename = sys.argv[2] #0 is python script!
print "Map is in : " + mapfilename
outputfilename = sys.argv[3] #0 is python script!
print "Output file : " + outputfilename


#load file
def loadfile (datafile):
    alldata = np.genfromtxt(fname = datafile,dtype=str)

    #indices = range(len(alldata[:,0]))
    #ID1 = map(str, alldata[:,0].tolist())
    #print alldata[:,0]
    #ID1chr = map(int, alldata[:,1].tolist())
    #ID2 = map(str, alldata[:,2].tolist())
    #ID2chr = map(int, alldata[:,3].tolist())
    #chr = map(int, alldata[:,4].tolist())
    #cM1 = map(float, alldata[:,5].tolist())
    #cM2 = map(float, alldata[:,6].tolist())
    ##ATTENTION - filtered data have an empty column in 7 so here we should shift everything
    SNP1 = map(str, alldata[:,7].tolist())
    SNP2 = map(str, alldata[:,8].tolist())
    #notsure = map(int, alldata[:,9].tolist())
    #length = map(float, alldata[:,10].tolist())

    #ind_dict = dict(zip(alldata[:,0], indices))

    datalist = []
    datalist.extend([SNP1, SNP2])
    print SNP1[1:10]
    print SNP2[1:10]

    return datalist

#load map file
def loadmapfile (mapfile):
    alldata = np.genfromtxt(fname = mapfile, delimiter='\t',dtype=str)

    #indices = range(len(alldata[:,0]))
    #chr = map(int, alldata[:,0].tolist())
    SNP = map(str, alldata[:,1].tolist())
    #cM = map(float, alldata[:,2].tolist())
    #bp = map(int, alldata[:,3].tolist())

    #ind_dict = dict(zip(alldata[:,0], indices))
    #print SNP[1:10]

    return SNP

#Write file
def writefile (outputfilename,dict,pos):
    outputfile = open(outputfilename, "w")
    for i in range(len(pos)):
        outputfile.write(str(pos[i]) + "\t" + str(dict[pos[i]]) + "\n")

    outputfile.close()

	
#Create the bincounts dict
#def bincount (data):
#	y = np.bincount(data)
#	print y
#	y2 = np.nonzero(y)[0]
#	print y2
#
#	bin_dict = zip(y2,y[y2])
#	return bin_dict

	
#let's go!
print 'loading data file...'
data = loadfile(datafilename)
print len(data[0])
print 'loading map file...'
positions = loadmapfile(mapfilename)
print len(positions)
#print mymap

#Build dict of positions
pos_dict = dict(zip(positions, [0]*len(positions)))
#print pos_dict[positions[0]]

#fill dict with all positions found in IBD sharing data


for i in range(len(data[0])):

    for j in range(positions.index(data[0][i]),(positions.index(data[1][i])+1),1):
        pos_dict[positions[j]]+=1
        #print str(positions[368]) + ":" + str( pos_dict[positions[368]])

print 'Writing output file...'
writefile (outputfilename,pos_dict,positions)

print "Finished!!"
