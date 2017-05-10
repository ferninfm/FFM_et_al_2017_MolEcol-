#------------------
# 1. Batch of code intended to do phylogenetic evaluation of OTU assignments. Not used in the final version
#------------------------
#
#library(Biostrings)
#seqs<-readLines("/Users/ferninfm/Desktop/ANTONIA/All_Reads_fer.fasta")
#sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
#rownames(sigma)<-c("a","c","g","t")
#colnames(sigma)<-c("a","c","g","t")
#pairwiseAlignment(seqs[2],seqs[88], type="global",substitutionMatrix=sigma, fuzzyMatrix=NULL, gapOpening=-10, gapExtension=-4, scoreOnly=TRUE)
#out<-mat.or.vec(length(seqs)/2,length(seqs)/2)
#Calculate distance matrix Watch out, takes 5 days
#for(i in seq(2,length(seqs),2))
#{
#	for (j in seq(2,length(seqs),2))
#	{
#	out[i/2,j/2]<-pairwiseAlignment(seqs[i],seqs[j], type="global",substitutionMatrix=sigma, fuzzyMatrix=NULL, gapOpening=-1, gapExtension=-1, scoreOnly=TRUE)
#	}
#}
#Normalize distance matrix in % of the shortest sequence NOPE
#out.norm<-mat.or.vec(dim(out)[1],dim(out)[2])
#for(i in 1:dim(out)[1])
#{
#	for (j in i:dim(out)[1])
#	{
#	normalizing.size<-min(c(out[i,i],out[j,j]))
#	out.norm[i,j]<-(2*(out[i,j]+normalizing.size)/(3*normalizing.size))
#	out.norm[j,i]<-out.norm[i,j]
#	}
#}
seqs<-read.dna("/Users/ferninfm/Desktop/ANTONIA/All_Reads_fer.fasta",format="fasta")
#locs<-unlist(strsplit(names(seqs),"_"))
#dim(locs)<-c(2,dim(out)[1])
#locs<-locs[1,]
#locs<-table(names(seqs),locs)
#colnames(out)<-names(seqs)
#rownames(out)<-names(seqs)
#upgma.tree<-hclust(dist(out),method="centroid") #watch out changed it
#locs<-locs[upgma.tree$labels[upgma.tree$order],]
#Plot things
#pdf("~/Desktop/ANTONIA/FINAL_heatmap.pdf",paper="a4")
#heatmap(out)
#dev.off()
#pdf("~/Desktop/ANTONIA/first_hclust.pdf",width=150,height=40)
#layout(matrix(c(1,2), 2, 1, byrow = TRUE), heights=c(5,1))
#par(mai=c(0.5,0.5,0.5,0.5))
#plot(upgma.tree)
#par(mai=c(0.5,0.5,0.5,0.5))
#plot(c(1,2,3,4),xlim=c(1,dim(locs)[1]),ylim=c(1,dim(locs)[2]),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
#for (i in 1:dim(locs)[1])
#	{for (j in 1:dim(locs)[2])
#		{points(i,j,type="p",pch=21,bg=locs[i,j])
#		}
#	}
#dev.off()
#
#
# Input a distance based alignment
# NOT USED
#
#pairwise.identity.alignment<-read.table(file = "/Users/ferninfm/Desktop/ANTONIA/alignment_pairwise_similarity.csv", header = TRUE, row.names=1,sep=",",)
#
#
#Prepare stuff for Megan
for (i in 1:dim(locs)[2])
{
write.dna(seqs[rownames(locs)[locs[,i]==1]],paste("/Users/ferninfm/Desktop/ANTONIA/Splidata_",colnames(locs)[i],".fas",sep=""),format="fasta",nbcol=-1,colsep="")
}

#
# 2. INPUT MEGAN TAXONOMY
#
# 2.1. In the Shell I concatenated all text files using 'cat Taxonomy*.txt>All_Tax_FINAL.txt'
#
taxonomy.megan<-readLines("/Users/ferninfm/Desktop/ANTONIA/1_Collapsed/2_Meganv2/All_Tax_FINAL.txt")
taxonomy.megan<-unlist(strsplit(taxonomy.megan,"; ; "))
dim(taxonomy.megan)<-c(2,length(taxonomy.megan)/2)
taxonomy.megan<-t(taxonomy.megan)
rownames(taxonomy.megan)<-taxonomy.megan[,1]
taxonomy.megan<-taxonomy.megan[,2]
#taxonomy.megan<-gsub("Not assigned; 100;","Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100; Not assigned; 100;",taxonomy.megan)
#
# 2.2. Plot Annotated tree
#
pdf("~/Desktop/ANTONIA/Megan_hclust.pdf",width=30,height=180)
foo<-strsplit(taxonomy.megan,";")
foo.labels<-mat.or.vec(1,length(foo))
for (i in 1:length(foo))
{
	if (length(foo[[i]])==2)
	{
		foo.labels[i]<-foo[[i]][1]
	}
	if (length(foo[[i]])>3&length(foo[[i]])<15)
	{
		foo.labels[i]<-foo[[i]][length(foo[[i]])-1]
	}
	else if (length(foo[[i]])>15)
	{
		foo.labels[i]<-paste(foo[[i]][length(foo[[i]])-5],foo[[i]][length(foo[[i]])-3],foo[[i]][length(foo[[i]])-1],sep="_")
	}
}
rm(foo)
names(foo.labels)<-names(taxonomy.megan)
upgma.tree.2<-upgma.tree
upgma.tree.2$labels<-paste(upgma.tree.2$labels,foo.labels[upgma.tree.2$labels],sep="_")
par(mai=c(0.5,0.5,0.5,0.5))
library(ape)
plot(as.phylo(upgma.tree.2))
write.tree(as.phylo(upgma.tree.2),"~/Desktop/ANTONIA/Megan_hclust.tre")
dev.off()
#
# 3. INPUT UNITE TAXONOMY
# I prepared a file in the terminal
#cat *.txt>All_UNITE_FINAL.txt
#cat All_UNITE_FINAL.txt| sed 's/##/E/g'|sed 's/|/;/g'|sed 's/s__;/s__/g'>All_UNITE_FILT.txt
#
blast.unite<-read.table("/Users/ferninfm/Desktop/ANTONIA/1_Collapsed/3_UNITE_Blast/All_UNITE_FILT.txt",sep=",")
#par(mfcol=c(2,2))
pdf("/Users/ferninfm/Desktop/ANTONIA/1_Collapsed/3_UNITE_Blast/boxplots.pdf")
#
# Function category.exclude
#
category.exclude<-function(category)
{
cat.list<-c("c__unidentified","p__unidentified","c__unidentified","o__unidentified","f__unidentified","g__unidentified","s__unidentified")
return(c(cat.list[category],"0"))
}	
#
# DEFINE OUTPUT MATRIX
#
unite.classes<-mat.or.vec(length(names(taxonomy.megan)),7) #PREPARE DATASET
rownames(unite.classes)<-names(taxonomy.megan)
#
# PARSE TAXONOMY
#
for (i in names(taxonomy.megan))
		{
			foo<-blast.unite[blast.unite[,1]==i,]
			if (dim(foo)[1]>0)
				{
				foo2<-strsplit(as.character(foo[,2]),";")
				fool<-length(foo2)
				for (j in 1:fool)
				{
					if (length(foo2[[j]])<8)
					{
						foo2[[j]]<-c(foo2[[j]],mat.or.vec(8-length(foo2[[j]]),1))
					}
					if (length(foo2[[j]])>8)
					{
						foo2[[j]]<-c(foo2[[j]][c(1:7)],paste(foo2[[j]][8],foo2[[j]][9]))
					}
				}
				foo2<-unlist(foo2)
				dim(foo2)<-c(8,fool)
				foo2<-t(foo2)
				#
				# Down to here is the same for every Taxonomic category
				#
				#
				for (rank in 2:7)
					{
						boxplot(log(foo[,11])~foo2[,rank+1],main=paste("Sequence",i))#FIRST USE OF RANK
						groups.class<-foo2[!duplicated(foo2[,rank+1]),rank+1]
						if (dim(foo)[1]==0)
							{
							print (paste(i,"has no matches or was not blasted"))
							}
						if (length(groups.class)==1)
							{
								unite.classes[i,rank]<-groups.class
							}
						if (length(groups.class)>1)
							{
								a<-mat.or.vec(1,length(groups.class))
								for (k in 1:length(a))
									{
										a[k]<-mean(log(foo[foo2[,rank+1]==groups.class[k],11]))
									}
								c.e<-category.exclude(rank)
								if (length (groups.class[a==min(a[groups.class%in%c.e==FALSE])])==1)
								{
									unite.classes[i,rank]<-groups.class[a==min(a[groups.class%in%c.e==FALSE])]
								}
								if (length (groups.class[a==min(a[groups.class%in%c.e==FALSE])])>1)
								{
									unite.classes[i,rank]<-paste(groups.class[a==min(a[groups.class%in%c.e==FALSE])],collapse="_")
								}
							}	
					}
				}
			}
dev.off()	






















#unite.classes<-mat.or.vec(length(names(taxonomy.megan)),7)
#names(unite.classes)<-names(taxonomy.megan)
#for (rank in 1:7)
#for (i in names(taxonomy.megan))
#	{
#		foo<-blast.unite[blast.unite[,1]==i,]
#		if (dim(foo)[1]>0)
#			{
#			foo2<-strsplit(as.character(foo[,2]),";")
#			fool<-length(foo2)
#			for (j in 1:fool)
#			{
#				if (length(foo2[[j]])<8)
#				{
#					foo2[[j]]<-c(foo2[[j]],mat.or.vec(8-length(foo2[[j]]),1))
#				}
#				if (length(foo2[[j]])>8)
#				{
#					foo2[[j]]<-c(foo2[[j]][c(1:7)],paste(foo2[[j]][8],foo2[[j]][9]))
#				}
#			}
#			foo2<-unlist(foo2)
#			dim(foo2)<-c(8,fool)
#			foo2<-t(foo2)
#			boxplot(log(foo[,11])~foo2[,4],main=paste("Sequence",i))
#			groups.class<-foo2[!duplicated(foo2[,4]),4]
#			if (length(groups.class)==1)
#			{
#				unite.classes[i]<-groups.class
#				}
#			if (length(groups.class)>1)
#				{
#					a<-mat.or.vec(1,length(groups.class))
#					for (k in 1:length(a))
#						{
#							a[k]<-mean(log(foo[foo2[,4]==groups.class[k],11]))
#						}
#					unite.classes[i]<-(groups.class[a==min(a[groups.class!="c__unidentified"])])
#					}	
#			}
#			if (dim(foo)[1]==0)
#			{
#				print (paste(i,"has no matches or was not blasted"))
#			}
#	}
#	
