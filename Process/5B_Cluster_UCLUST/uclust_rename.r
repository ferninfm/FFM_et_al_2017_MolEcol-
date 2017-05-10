#!/usr/bin/env Rscript
library(ape)
file<-read.dna("./All.fa",format="fasta")
names<-strsplit(names(file)," ")
l<-length(names)
dum<-mat.or.vec(l,1)
for (i in 1:l)
{
	if(length(names[[i]])==1)
		{
			names[[i]]<-paste (names[[i]][1],";size=1;",sep="")
			dum[i]<-1
		}
	if(length(names[[i]])>1)
		{
			size<-substr(names[[i]][3],start=12,stop=17)
			names[[i]]<-paste (names[[i]][1],";size=",(size+1),";",names[[i]][4], sep="")
			dum[i]<-size
		}
}
dum<-cbind(seq(1,l),dum)
names(file)<-names
file<-file[as.numeric(dum[order(as.numeric(dum[,2]),as.numeric(dum[,1]),decreasing=TRUE),1])]
write.dna(file,"./All_renamed.fa",format="fasta")