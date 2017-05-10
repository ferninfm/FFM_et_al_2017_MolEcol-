load("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Start_Workspace")
seq.dataset[,5:8]<-0
library(ape)
seqs.its<-read.dna("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/4C_ITSx/All_its.fa.ITS1.fasta",format="fasta")
its.names<-strsplit(names(seqs.its),";")
for (i in 1:length(its.names))
{
	its.names[i]<-unlist(its.names[i])[1]
}
its.names<-unlist(its.names)
names(seqs.its)<-its.names
include.its<-as.numeric(names(seqs)%in%its.names)
foo<-readLines("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/5C_SWARM/its1.otus.fas")
for (i in 1:dim(seq.dataset)[1])
	{
	if(include.its[i]==1)
		{
		seq.dataset[i,5]<-grep(seq.dataset[i,2],foo)
		}
	}
colnames(seq.dataset)<-c("sample", "reference", "size", "megan.taxonomy", "OTU", "exclude.itsx", "exclude.wirt", "exclude.lichen")
seq.dataset[,6]<-as.numeric(!include.its)
seq.dataset[,8]<-exclude
seqs.otus<-read.dna("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/5C_SWARM/ITS1_OTUs",format="fasta")
otus.names<-unlist(strsplit(names(seqs.otus),"_"))
dim(otus.names)<-c(2,length(otus.names)/2)
otus.names<-t(otus.names)
cbind(seq.dataset,as.numeric(!(seq.dataset[,5]%in%which(otus.names[,2]>1))))
seq.dataset<-cbind(seq.dataset,as.numeric(!(seq.dataset[,5]%in%which(otus.names[,2]>1))))
colnames(seq.dataset)[9]<-"exclude.smallotus"
save.image("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Second_Workspace")
#
rownames(unite_min_new)<-unite_min_new[,1]
foobar<-unite_min_new[otus.names[,1],]
rownames(foobar)<-otus.names[,1]
seq.dataset<-cbind(seq.dataset,mat.or.vec(dim(seq.dataset)[1],dim(unite_min_new)[2]))
for (i in 1:dim(seq.dataset)[1])
	{
	j<-as.numeric(seq.dataset[i,5])
	if (j!=0)
		{
		seq.dataset[i,10:26]<-as.matrix(foobar[j,])
		}
	}
save.image("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Corrected_Workspace")

# I edited the dataset to reinclude the chimeras from uclust that I had mistakenly taken out
counts.per.otu<-table(seq.dataset[,5])
for (i in names(counts.per.otu))
counts.per.otu[i]<-sum(as.numeric(seq.dataset[seq.dataset[,5]==i,3]))
seq.subset<-seq.dataset[seq.dataset[,6]==0&seq.dataset[,8]==0,]
save.image("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Corrected_Final_Workspace")
table.otus<-table(seq.subset[,1],as.numeric(seq.subset[,5]))
for(i in rownames(table.otus))
{
	for (j in colnames(table.otus))
	{
		if (table.otus[i,j]!=0)
		{
		table.otus[i,as.character(j)]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
table.clases<-table(seq.subset[,1],seq.subset[,19])
for(i in rownames(table.clases))
{
	for (j in colnames(table.clases))
	{
		if (table.clases[i,j]!=0)
		{
		table.clases[i,as.character(j)]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,19]==j,3]),na.rm=T)
		}
	}
}
# Reinterpret Fibulobasidium
seq.dataset[seq.dataset[,4]==" Fibulobasidium",13]<-"Fungi"
seq.dataset[seq.dataset[,4]==" Fibulobasidium",15]<-"Basidiomycota"
seq.dataset[seq.dataset[,4]==" Fibulobasidium",17]<-"Tremellomycetes"
seq.dataset[seq.dataset[,4]==" Fibulobasidium",19]<-"Tremellales"

seq.subset<-seq.dataset[seq.dataset[,6]==0&seq.dataset[,8]==0,]

seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",13]<-"Fungi"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",15]<-"Ascomycota"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",17]<-"Eurotiomycetes"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",19]<-"Verrucariales"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",21]<-"Verrucariaceae"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",23]<-"Endococcus?"
seq.dataset[seq.dataset[,5]%in%c(17,120,1810,1817)|seq.dataset[,4]==" Bagliettoa marmorea",25]<-"Endococcus macrosporus"

seq.subset<-seq.dataset[seq.dataset[,6]==0&seq.dataset[,8]==0,]

#Export datasets
colnames(seq.dataset)<-c("sample", "reference", "size", "megan.taxonomy", "OTU", "exclude.itsx", "exclude.lc", "exclude.lichen", "exclude.smallotus", "ref.otu", "size.otu", "hit.otu", "kingdom", "evk", "division", "evd", "class", "evc", "order", "evo", "family", "evf", "genus", "evg", "species", "evs")
colnames(seq.subset)<-c("sample", "reference", "size", "megan.taxonomy", "OTU", "exclude.itsx", "exclude.lc", "exclude.lichen", "exclude.smallotus", "ref.otu", "size.otu", "hit.otu", "kingdom", "evk", "division", "evd", "class", "evc", "order", "evo", "family", "evf", "genus", "evg", "species", "evs")
write.table(cbind(seq.dataset,megan.taxonomy),"/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Full_dataset1.txt")
write.table(seq.subset,"/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Reduced_dataset1.txt")
write.table(table.otus,"/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/table_otus.txt")
write.table(table.clases,"/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/table_clases.txt")

#pie(table.clases[10,table.clases[10,]!=0])
library(vegan)
library(vegetarian)

tabulate.and.print<-function(data=seq.subset,samples=1,groups=19,counts=3,category="Class")
{
	table.clases<-table(data[,samples],data[,groups])
	for(i in rownames(table.clases))
	{
		for (j in colnames(table.clases))
			{
				if (table.clases[i,j]!=0)
				{
					table.clases[i,as.character(j)]<-sum(as.numeric(data[data[,samples]==i&data[,groups]==j,counts]),na.rm=T)
				}
			}
	}
require(vegetarian)
table.clases[,]<-normalize.rows(table.clases)
dc<-data.frame(table.clases)
dc<-dc[dc$Freq!=0,]
require(ggplot2)
p1<-ggplot(dc, aes(x = factor(Var1), y=Freq,fill=factor(Var2))) + geom_bar(stat="identity",position=("stack"))+ theme(axis.text.x=element_text(angle=90), legend.position="right", legend.text=element_text(size=6)) + guides(fill=guide_legend(ncol=3))+ labs(x="" , y="", title=paste("Ascription of seqences to",category, sep=" "))+ scale_fill_discrete()
print(p1)
return(list(table.clases,p1))
}

opla<-tabulate.and.print(groups=19)


# I redid the unite_min step
load(".RData")
load("~/Desktop/ANTONIA/2_Not_collapsed/6C_BLAST_UNITE_ITSx_SWARM/Workspace_full.Rdata")
rownames(unite_min_new)<-as.character(unite_min_new[,1])
unite_min_new<-as.matrix(unite_min_new)
seq.dataset<-as.matrix(seq.dataset)
seq.dataset[is.na(seq.dataset[,10]),10]<-"NO_OTU"
for(i in rownames(unite_min_new))
{
	for (j in 10:26)
	{
		seq.dataset[seq.dataset[,10]==i,j]<-unite_min_new[i,(j-9)]
	}
}
save.image("/Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/7C_Alltogether/Corrected_Final_Workspace")
#Again parse MEGAN
parse.megan<-function(data=megan.taxonomy,end.val=length(data))
	{
	megan.taxonomy.parsed<-mat.or.vec(length(data),7)
	colnames(megan.taxonomy.parsed)<-c("Kingdom","Division","Class","Order","Family","Genus","Species")
	for (i in 1:end.val)
		{
			if(length(grep(" Not assigned; 100;",data[i]))==1)
			{
			megan.taxonomy.parsed[i,]<-"Unknown"
			} else {
			if(length(grep(" No hits; 100;",data[i]))==1)
			{
			megan.taxonomy.parsed[i,]<-"Unknown"
			} else {
			if(length(grep("Bacteria",data[i]))==1)
			{
			megan.taxonomy.parsed[i,1]<-"Bacteria"
			} else {
			if(length(grep("Animalia",data[i]))==1)
			{
			megan.taxonomy.parsed[i,1]<-"Animalia"
			} else {
			if(length(grep("Fungi",data[i]))==1)
				{
				megan.taxonomy.parsed[i,1]<-"Fungi"
				control1<-strsplit(data[i],";")
				megan.taxonomy.parsed[i,2]<-substr(grep("mycota",control1[[1]],value=TRUE)[1],2,100)
				megan.taxonomy.parsed[i,3]<-substr(grep("mycetes",control1[[1]],value=TRUE)[1],2,100)
				megan.taxonomy.parsed[i,4]<-substr(grep("ales",control1[[1]],value=TRUE)[1],2,100)
				megan.taxonomy.parsed[i,5]<-substr(grep("aceae",control1[[1]],value=TRUE)[1],2,100)
				megan.taxonomy.parsed[i,7]<-substr(control1[[1]][(length(control1[[1]])-1)],2,100)
				megan.taxonomy.parsed[i,6]<-substr(control1[[1]][(length(control1[[1]])-3)],2,100)
				if(length(grep("aceae",megan.taxonomy.parsed[i,6]))==1)
				{
					megan.taxonomy.parsed[i,6]<-megan.taxonomy.parsed[i,7]
				}
					
			}
			}
			}
		}
	}
		}
		megan.taxonomy.parsed[megan.taxonomy.parsed=="0"]<-"Unknown"
		megan.taxonomy.parsed[is.na(megan.taxonomy.parsed)]<-"Unknown"
		return(megan.taxonomy.parsed)
}

taxonomia<-parse.megan(megan.taxonomy)