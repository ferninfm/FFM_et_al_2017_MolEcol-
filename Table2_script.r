#___________________________#
#                           #
# Table 2, Manuscript       #
#                           #
#___________________________#

#
# 1. All sequences
#
table_1foo<-tabulate.and.print(data=seq.dataset,groups=5,normalize=FALSE,print=FALSE)
table_1foo<-table_1foo[[1]]
table_1foo<-table_1foo[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26),]
table_1<-rowSums(table_1foo)
#
# 2. Of which unique
#
table_1<-cbind(table_1,table(seq.dataset[,1])[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26)])
#
# 3. ITs 1
#
table_1<-cbind(table_1,rowSums(table_1foo[,-1]))
#
# 4. Not Host
#
table_1foo<-tabulate.and.print(data=seq.dataset[exclude[,8]%in%c("000","001"),],groups=5,normalize=FALSE,print=FALSE)
table_1foo<-table_1foo[[1]]
table_1foo<-table_1foo[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26),]
table_1<-cbind(table_1,rowSums(table_1foo))
#
# 5. Lichen
#
table_1foo<-tabulate.and.print(data=seq.subset,groups=5,normalize=FALSE,print=FALSE)
table_1foo<-table_1foo[[1]]
table_1foo<-table_1foo[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26),]
table_1<-cbind(table_1,rowSums(table_1foo))
#
##
# Threshold
##
#
#
#
#
#
#
table_1<-cbind(table_1,table(seq.dataset[exclude[,8]%in%c("000","001","011"),1])[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26)])
table_1<-cbind(table_1,table(seq.dataset[exclude[,8]%in%c("000","001"),1])[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26)])
table_1<-cbind(table_1,table(seq.dataset[exclude[,8]%in%c("000"),1])[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26)])
#
foo.len<-dim(table_1)[2]
table_1foo<-mat.or.vec(3,foo.len)
table_1foo[2,]<-colMeans(table_1)
for (i in 1:foo.len)
	{
	table_1foo[1,i]<-max(table_1[,i])
	table_1foo[3,i]<-min(table_1[,i])
	}
table_1<-rbind(table_1,table_1foo)
rm(table_1foo,foo.len)


#
#
#
#table_3
#
shared.reads<-table(seq.subset[,1],seq.subset[,5])
for (i in rownames(shared.reads))
{
	for (j in colnames(shared.reads))
	{
		if (shared.reads[i,j]!=0)
		{
		shared.reads[i,j]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
shared.otus<-shared.reads
shared.otus[,]!=0<-1
shared.otus.control<-colSums(shared.otus)
taxon.control<-table(seq.subset[,5],seq.subset[,19])[,c("Tremellales","Capnodiales","Chaetothyriales","Botryosphaeriales")]
taxon.control[rowSums(taxon.control)!=0,]
table_3<-cbind(
	rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,1]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,1]!=0]]>1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,1]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,1]!=0]]==1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,1]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,1]!=0]]>1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,1]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,1]!=0]]==1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,2]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,2]!=0]]>1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,2]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,2]!=0]]==1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,2]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,2]!=0]]>1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,2]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,2]!=0]]==1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,3]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,3]!=0]]>1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,3]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,3]!=0]]==1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,3]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,3]!=0]]>1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,3]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,3]!=0]]==1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,4]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,4]!=0]]>1])
	,rowSums(shared.otus[,row.names(taxon.control)[taxon.control[,4]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,4]!=0]]==1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,4]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,4]!=0]]>1])
	,rowSums(shared.reads[,row.names(taxon.control)[taxon.control[,4]!=0]][,shared.otus.control[row.names(taxon.control)[taxon.control[,4]!=0]]==1])
	,rowSums(shared.otus[,shared.otus.control>1])
	,rowSums(shared.otus[,shared.otus.control==1])
	,rowSums(shared.reads[,shared.otus.control>1])
	,rowSums(shared.reads[,shared.otus.control==1])
	)
	table_3<-table_3[c(1,2,4,6,7,12,13,15,16,17,20,22,25,3,5,8,9,10,11,14,18,19,21,23,24,26),]