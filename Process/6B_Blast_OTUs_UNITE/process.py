import sys
#import uc
#import die
#import fasta

#!FileName = sys.argv[1] #reads filename from the console
FileName= "/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Blast_Unite_OTUS.txt"
FileOut= "/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Processed_Blast_Unite_OTUS.txt"
#
#
# PARSE UNITE OUTPUT TO NEW FILE
# RETAIN OTU LIST
#
OTUNames=["query","size","ref","kingdom", "division", "classe", "order", "family", "genus", "species","1","2","3","4","5","6","7","8","9","10"]
OTUTable=["0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"]
#line='OTU_1;size=5787,JN906464|k__Fungi;p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__uncultured,93.67,221,7,6,1,220,24,238,1e-87,  324'
with open(FileOut,'a') as f:
    with open(FileName) as infile:
        for line in infile:
            splitline0=line.split(",")              # Parse most fields
            OTUTable[10:20]=splitline0[2:12]
            splitline1=splitline0[0].split(";")     # Parse OTU name from repetitions FIRST TO KEEP
            OTUTable[0]=splitline1[0]
            OTUTable[1]=splitline1[1][5:]
            splitline1=splitline0[1].split("|k__")     # Separate sequence number from Taxonomy
            OTUTable[2]=splitline1[0]
            splitline1=splitline1[1].split(";p__")     # Split taxonomy kingdom
            OTUTable[3]=splitline1[0]
            if "c__" in splitline1[1]:
                splitline1=splitline1[1].split(";c__")     # Split taxonomy
                OTUTable[4]=splitline1[0]
            else:
                OTUTable[4]="NA"
            
            if "o__" in splitline1[1]:
                splitline1=splitline1[1].split(";o__")     # Split taxonomy
                OTUTable[5]=splitline1[0]
            else:
                OTUTable[5]="NA"
            
            if "f__" in splitline1[1]:
                splitline1=splitline1[1].split(";f__")     # Split taxonomy
                OTUTable[6]=splitline1[0]
            else:
                OTUTable[6]="NA"
            
            if "g__" in splitline1[1]:
                splitline1=splitline1[1].split(";g__")     # Split taxonomy
                OTUTable[7]=splitline1[0]
            else:
                OTUTable[7]="NA"
            
            if "s__" in splitline1[1]:
                splitline1=splitline1[1].split(";s__")     # Split taxonomy
                OTUTable[8]=splitline1[0]
                OTUTable[9]=splitline1[1]
            else:
                OTUTable[8]="NA"
                OTUTable[9]="NA"
            f.write(';'.join(OTUTable))
            




#Define functions
# def GetSampleId(Label):
# 	Fields = Label.split(";")
# 	for Field in Fields:
# 		if Field.startswith("barcodelabel="):
# 			return Field[13:]
# 	die.Die("barcodelabel= not found in read label '%s'" % Label)
# 
# def OnRec():
# 	global OTUs, Samples, OTUTable
# 	if uc.Type != 'H':
# 		return
# 
# 	OTUId = uc.TargetLabel
# 	if OTUId not in OTUIds:
# 		OTUIds.append(OTUId)
# 		OTUTable[OTUId] = {}
# 
# 	SampleId = GetSampleId(uc.QueryLabel)
# 	if SampleId not in SampleIds:
# 		SampleIds.append(SampleId)
# 
# 	N = fasta.GetSizeFromLabel(uc.QueryLabel, 1)
# 	try:
# 		OTUTable[OTUId][SampleId] += N
# 	except:
# 		OTUTable[OTUId][SampleId] = N
# #Define elements
# OTUIds = []
# SampleIds = []
# OTUTable = {}
# 
# >>> if x < 0:
# ...     x = 0
# ...     print 'Negative changed to zero'
# ... elif x == 0:
# ...     print 'Zero'
# ... elif x == 1:
# ...     print 'Single'
# ... else:
# ...     print 'More'
# ...
# More
# 
# 
# 
# 
# 
# uc.ReadRecs(FileName, OnRec)
# 
# s = "OTUId"
# for SampleId in SampleIds:
# 	s += "\t" + SampleId
# print s
# 
# for OTUId in OTUIds:
# 	s = OTUId
# 	for SampleId in SampleIds:
# 		try:
# 			n = OTUTable[OTUId][SampleId]
# 		except:
# 			n = 0
# 		s += "\t" + str(n)
# 	print s
# 