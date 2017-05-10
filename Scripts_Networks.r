#----------------------------------------------------------------------------#
#                                                                            #
#                      Script to draw Circos plots                           #
#                                 FFM  8.4.2017                              #
#                                                                            #
#----------------------------------------------------------------------------#
library(circlize)
setwd("/Volumes/Dropbox/Dropbox/11_Github_repositories/FFM_et_al_2017/FFM_et_al_2017/")
#load("/Volumes/Dropbox/Dropbox/11_Github_repositories/FFM_et_al_2017/FFM_et_al_2017/Workspace_corrected_V5")
pdf("./Fig4_THRESHOLD_mean.pdf")

#------------------
# 1. GEt total reads
#------------------

shared.reads<-table(seq.dataset[,1],seq.dataset[,5])
for (i in rownames(shared.reads))
{
	for (j in colnames(shared.reads))
	{
		if (shared.reads[i,j]!=0)
		{
		shared.reads[i,j]<-sum(as.numeric(seq.dataset[seq.dataset[,1]==i&seq.dataset[,5]==j,3]))
		}
	}
}
total.reads<-rowSums(shared.reads)

#------------------
#2. Get all.shared reads in the SUBSET
#------------------

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
reads.foo<-cbind(rowSums(shared.reads[,colSums(shared.reads !=0)>1]), rowSums(shared.reads))
reads.foo<-cbind(reads.foo[,1]/reads.foo[,2],(reads.foo[,2]-reads.foo[,1])/reads.foo[,2])
#
# OLD MANUAL FIX FOR HAVING NO READS IN R405
#
#shared.reads<-rbind(shared.reads[1:18,],0,shared.reads[19:25,])
#rownames(shared.reads)[19]<-names(total.reads)[19]

#------------------
# 3.Normalize
#------------------

subset.reads<-rowSums(shared.reads)
total.reads.threshold<-subset.reads
for (i in names(total.reads.threshold)) total.reads.threshold[i]<-sum(as.numeric((seq.dataset[seq.dataset[,5]!=0&seq.dataset[,1]==i,3])))
#FIRST IDEA WILL NOT WORK
#for (i in rownames(shared.reads)) shared.reads[i,]<-shared.reads[i,]/total.reads[i]
# SECOND WILL
for (i in rownames(shared.reads)) shared.reads[i,]<-shared.reads[i,]/subset.reads[i]
total.normalized<-total.reads/sum(total.reads)

#------------------
# 4. Reorder for plot
#
# Note that sample Rhi_geo+Mue_pyg+A405 is from here on excluded from the graphs
#
#------------------

reorder<- c(1,2,3,4,6,5,7,12,8,9,10,11,13,15,14,16,17,20,18,22,21,23,24,25,26)
total.normalized<-total.normalized[reorder]
total.reads<-total.reads[reorder]
shared.reads<-shared.reads[reorder,]
subset.reads<-subset.reads[reorder]


#----------------------------------------------
# 5. Get all.shared TREMELLALES reads in the SUBSET
#----------------------------------------------

shared.trem.reads<-table(seq.subset[seq.subset[,17]=="Tremellomycetes",1],seq.subset[seq.subset[,17]=="Tremellomycetes",5])
top.trem<-names(colSums(shared.trem.reads)[order(colSums(shared.trem.reads),decreasing=T)][1:3])
for (i in rownames(shared.trem.reads))
{
	for (j in colnames(shared.trem.reads))
	{
		if (shared.trem.reads[i,j]!=0)
		{
		shared.trem.reads[i,j]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
reads.trem<-cbind(rowSums(shared.trem.reads[,colSums(shared.trem.reads !=0)>1]), rowSums(shared.trem.reads))
reads.trem<-cbind(reads.trem[,1]/reads.trem[,2],(reads.trem[,2]-reads.trem[,1])/reads.trem[,2])
for (i in rownames(shared.trem.reads)) shared.trem.reads[i,]<-shared.trem.reads[i,]/subset.reads[i]
#------------------
# 6. Get Botryosphaeriales
#------------------
shared.bot.reads<-table(seq.subset[seq.subset[,19]=="Botryosphaeriales",1],seq.subset[seq.subset[,19]=="Botryosphaeriales",5])
top.bot<-names(colSums(shared.bot.reads)[order(colSums(shared.bot.reads),decreasing=T)][1:3])
for (i in rownames(shared.bot.reads))
{
	for (j in colnames(shared.bot.reads))
	{
		if (shared.bot.reads[i,j]!=0)
		{
		shared.bot.reads[i,j]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
reads.bot<-cbind(rowSums(shared.bot.reads[,colSums(shared.bot.reads !=0)>1]), rowSums(shared.bot.reads))
reads.bot<-cbind(reads.bot[,1]/reads.bot[,2],(reads.bot[,2]-reads.bot[,1])/reads.bot[,2])
for (i in rownames(shared.bot.reads)) shared.bot.reads[i,]<-shared.bot.reads[i,]/subset.reads[i]
#------------------
# 7. Get Chaetothyriales
#------------------
shared.chae.reads<-table(seq.subset[seq.subset[,19]=="Chaetothyriales",1],seq.subset[seq.subset[,19]=="Chaetothyriales",5])
top.chae<-names(colSums(shared.chae.reads)[order(colSums(shared.chae.reads),decreasing=T)][1:3])
for (i in rownames(shared.chae.reads))
{
	for (j in colnames(shared.chae.reads))
	{
		if (shared.chae.reads[i,j]!=0)
		{
		shared.chae.reads[i,j]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
reads.chae<-cbind(rowSums(shared.chae.reads[,colSums(shared.chae.reads !=0)>1]), rowSums(shared.chae.reads))
reads.chae<-cbind(reads.chae[,1]/reads.chae[,2],(reads.chae[,2]-reads.chae[,1])/reads.chae[,2])
for (i in rownames(shared.chae.reads)) shared.chae.reads[i,]<-shared.chae.reads[i,]/subset.reads[i]
#------------------
# 8. Get Capnodiales
#------------------
shared.cap.reads<-table(seq.subset[seq.subset[,19]=="Capnodiales",1],seq.subset[seq.subset[,19]=="Capnodiales",5])
top.cap<-names(colSums(shared.cap.reads)[order(colSums(shared.cap.reads),decreasing=T)][1:3])
for (i in rownames(shared.cap.reads))
{
	for (j in colnames(shared.cap.reads))
	{
		if (shared.cap.reads[i,j]!=0)
		{
		shared.cap.reads[i,j]<-sum(as.numeric(seq.subset[seq.subset[,1]==i&seq.subset[,5]==j,3]))
		}
	}
}
reads.cap<-cbind(rowSums(shared.cap.reads[,colSums(shared.cap.reads !=0)>1]), rowSums(shared.cap.reads))
reads.cap<-cbind(reads.cap[,1]/reads.cap[,2],(reads.cap[,2]-reads.cap[,1])/reads.cap[,2])
for (i in rownames(shared.cap.reads)) shared.cap.reads[i,]<-shared.cap.reads[i,]/subset.reads[i]

#------------------
# 9. Reorder for plot
#------------------
shared.trem.reads<-shared.trem.reads[names(total.normalized)[names(total.normalized)%in%row.names(shared.trem.reads)],]
reads.trem<-reads.trem[names(total.normalized)[names(total.normalized)%in%row.names(reads.trem)],]
rownames(reads.trem)<-c(1:25)[names(total.normalized)%in%row.names(reads.trem)]
shared.bot.reads<-shared.bot.reads[names(total.normalized)[names(total.normalized)%in%row.names(shared.bot.reads)],]
reads.bot<-reads.bot[names(total.normalized)[names(total.normalized)%in%row.names(reads.bot)],]
rownames(reads.bot)<-c(1:25)[names(total.normalized)%in%row.names(reads.bot)]
shared.chae.reads<-shared.chae.reads[names(total.normalized)[names(total.normalized)%in%row.names(shared.chae.reads)],]
reads.chae<-reads.chae[names(total.normalized)[names(total.normalized)%in%row.names(reads.chae)],]
rownames(reads.chae)<-c(1:25)[names(total.normalized)%in%row.names(reads.chae)]
shared.cap.reads<-shared.cap.reads[names(total.normalized)[names(total.normalized)%in%row.names(shared.cap.reads)],]
reads.cap<-reads.cap[names(total.normalized)[names(total.normalized)%in%row.names(reads.cap)],]
rownames(reads.cap)<-c(1:25)[names(total.normalized)%in%row.names(reads.cap)]
reads.foo<-reads.foo[names(total.normalized)[names(total.normalized)%in%row.names(reads.foo)],]
rownames(reads.foo)<-c(1:25)[names(total.normalized)%in%row.names(reads.foo)]
#
starting.points<-shared.reads
rownames(starting.points)<-c(1:25)
library(reshape2)
starting.points<-melt(starting.points)
starting.points<-starting.points[starting.points[,3]!=0,]
starting.points<-starting.points[order(starting.points[,1],starting.points[,2]),]
starting.points<-cbind(starting.points,0,0)
for (i in 1:25)
{
	lilo<-cumsum(starting.points[starting.points[,1]==i,3])
	if(length(lilo)!=0)
	{
	starting.points[starting.points[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
	starting.points[starting.points[,1]==i,5]<-lilo
	}
}
colnames(starting.points)<-c("X1","X2","value","out","in")
#-----------------------------------------------------
#
# Threshold criteria for networks
#
#-----------------------------------------------------
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
exclude.lichens<-seq.dataset[liq.ds==1,]
exclude.lichens<-data.frame(exclude.lichens)
exclude.lichens$size<-as.numeric(as.matrix(exclude.lichens$size))
exclude.lichens<-aggregate(size~sample+OTU,exclude.lichens,sum)
names(exclude.lichens)[3]<-"size"
# Optional substract OTU 0
exclude.lichens<-exclude.lichens[exclude.lichens[,2]!=0,]
#
exclude.lichens<-aggregate(size~sample,exclude.lichens,mean)
#OTHER OPTIONS
#exclude.lichens<-aggregate(size~sample,exclude.lichens,sum)
#exclude.lichens<-aggregate(size~sample,exclude.lichens,max)
rownames(exclude.lichens)<-exclude.lichens[,1]
#for (i in rownames(exclude.lichens)) exclude.lichens[i,2]<-exclude.lichens[i,2]/total.reads[i]
for (i in rownames(exclude.lichens)) exclude.lichens[i,2]<-exclude.lichens[i,2]/subset.reads[i]
foo.exclude<-subset.reads
for (i in rownames(exclude.lichens)) foo.exclude[i]<-exclude.lichens[i,2]
foo.exclude[foo.exclude>1]<-0
foo.exclude[foo.exclude<0.001]<-0.001
#------------------------------------
#
# 10. Start Circos Plots
#
#------------------------------------

#------------------
# 10.1 Declare variables
#------------------
colorinos<-c("#a6cee3","#1f78b4","#1f78b4","#b2df8a","#33a02c","#33a02c","#21ac6a","#fb9a99","#fb9a99","#fb9a99","#fb9a99","#fb9a99","#e31a1c","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00","#cab2d6","#cab2d6","#6a3d9a","#6a3d9a","#6a3d9a","#6a3d9a","#266cdd","#4636e1")
colorinos.trans<-c("#a6cee344","#1f78b444","#1f78b444","#b2df8a44","#33a02c44","#33a02c44","#21ac6a44","#fb9a9944","#fb9a9944","#fb9a9944","#fb9a9944","#fb9a9944","#e31a1c44","#fdbf6f44","#fdbf6f44","#ff7f0044","#ff7f0044","#cab2d644","#cab2d644","#6a3d9a44","#6a3d9a44","#6a3d9a44","#6a3d9a44","#266cdd44","#4636e144")
infected<-c(0,0,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,0,1)
foo<-unlist(strsplit(names(total.normalized),"\\+"))
dim(foo)<-c(3,25)
foo<-t(foo)
foo<-gsub("_"," ",foo)
foo[,1]<-c("A. fuscata","A. myrinii","A. myrinii","C. vitellina","L. biccinta","L. biccinta","L. intricata","L. polytropa","L. polytropa","L. polytropa","L. polytropa","L. polytropa","L. swartzii","L. lapicida","L. lapicida","P. conglomerata","P. conglomerata","R. geographicum","R. geographicum","T. atra","T. atra","T. atra","T. atra", "U. cylindrica","V. lactea")
text.pos<-c(0.5,1,0.5,0.5,1,0.5,0.5,2.5,0.5,0.5,0.5,0.5,0.5,1,0.5,1,0.5,1.0,0.5,2.5,0.5,0.5,0.5,0.5,0.5)#
#----------
# 10.2 Make function to simplify
#----------
base.for.circos<-function()
{
	#------------------
	# Initialize
	#------------------
	circos.initialize(1:25,xlim=c(0,1))
	#------------------
	# Start text
	#------------------
	circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
	j<-0
	for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
	{
	circos.text(text.pos[i],0.6,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
	if (j!=0)
		{
		circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
		}
		j<-i
		if (j==25)
			{
			circos.lines(c(0.1,1),c(0.4,0.4),j)
			}
	}
	for (i in c(1:25)[infected==1])
	{
	circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
	}
	#------------------
	# Samples
	#------------------
	circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
	for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
}
#------------------------------------
#
# 11. Start Circos Tremellomycetes
#
#------------------------------------
base.for.circos()
#------
# Percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.trem))
{	circos.polygon(c(0,reads.trem[i,1],reads.trem[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.trem[i,1],reads.trem[i,1]+reads.trem[i,2],reads.trem[i,1]+reads.trem[i,2],reads.trem[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
#
# Generate Connector dataset
#
index<-c(1:26)[rownames(shared.reads)%in%rownames(shared.trem.reads)]
starting.trem<-starting.points[starting.points[,1]%in%index,]
starting.trem<-starting.trem[starting.trem[,2]%in%colnames(shared.trem.reads),]
for (i in index)
{
	lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
	if(length(lilo)>1)
	{
	starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
	if(length(lilo)==1)
	{
	starting.trem[starting.trem[,1]==i,4]<-0
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
}
#
# Plot Connectors
#
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}

text(-0.8,1,"Tremellales",cex=1.3,font=4)
#
#
#
#
base.for.circos()
#------
# Percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.trem))
{	circos.polygon(c(0,reads.trem[i,1],reads.trem[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.trem[i,1],reads.trem[i,1]+reads.trem[i,2],reads.trem[i,1]+reads.trem[i,2],reads.trem[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
#
# Generate Connector dataset
#
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
#
# Plot Connectors
#
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}

text(-0.8,1,"Tremellales",cex=1.3,font=4)

#------------------------------------
#
# 12. Start Circos BOT
#
#------------------------------------
#
base.for.circos()
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.bot))
	{	circos.polygon(c(0,reads.bot[i,1],reads.bot[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
		circos.polygon(c(reads.bot[i,1],reads.bot[i,1]+reads.bot[i,2],reads.bot[i,1]+reads.bot[i,2],reads.bot[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
	}
# Botryospheriales PLOT
index<-c(1:25)[rownames(shared.reads)%in%rownames(shared.bot.reads)]
starting.trem<-starting.points[starting.points[,1]%in%index,]
starting.trem<-starting.trem[starting.trem[,2]%in%colnames(shared.bot.reads),]
for (i in index)
{
	lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
	if(length(lilo)>1)
	{
	starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
	if(length(lilo)==1)
	{
	starting.trem[starting.trem[,1]==i,4]<-0
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
}
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Botryosphaeriales",cex=1.3,font=4)
#
#END CIRCOS
#
#
#
#
#
#
base.for.circos()
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.bot))
	{	circos.polygon(c(0,reads.bot[i,1],reads.bot[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
		circos.polygon(c(reads.bot[i,1],reads.bot[i,1]+reads.bot[i,2],reads.bot[i,1]+reads.bot[i,2],reads.bot[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
	}
# Botryospheriales PLOT
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Botryosphaeriales",cex=1.3,font=4)
#
#END CIRCOS
#--------------------------------------
#
#------------------------------------
#
# 13. Start Circos Chaetothyriales
#
#------------------------------------
#
base.for.circos()
	circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.chae))
{	circos.polygon(c(0,reads.chae[i,1],reads.chae[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.chae[i,1],reads.chae[i,1]+reads.chae[i,2],reads.chae[i,1]+reads.chae[i,2],reads.chae[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
# Chaetothyriales PLOT
index<-c(1:25)[rownames(shared.reads)%in%rownames(shared.chae.reads)]
starting.trem<-starting.points[starting.points[,1]%in%index,]
starting.trem<-starting.trem[starting.trem[,2]%in%colnames(shared.chae.reads),]
for (i in index)
{
	lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
	if(length(lilo)>1)
	{
	starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
	if(length(lilo)==1)
	{
	starting.trem[starting.trem[,1]==i,4]<-0
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
}
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Chaetothyriales",cex=1.3,font=4)
#
#
#
#
#
#
base.for.circos()
	circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.chae))
{	circos.polygon(c(0,reads.chae[i,1],reads.chae[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.chae[i,1],reads.chae[i,1]+reads.chae[i,2],reads.chae[i,1]+reads.chae[i,2],reads.chae[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
# Chaetothyriales PLOT
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
#
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Chaetothyriales",cex=1.3,font=4)
#END CIRCOS
#--------------------------------------
#
#------------------------------------
#
# 14. Start Circos Capnodiales
#
#------------------------------------
#
base.for.circos()
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.cap))
{	circos.polygon(c(0,reads.cap[i,1],reads.cap[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.cap[i,1],reads.cap[i,1]+reads.cap[i,2],reads.cap[i,1]+reads.cap[i,2],reads.cap[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}

#
# Capnodiales PLOT
index<-c(1:25)[rownames(shared.reads)%in%rownames(shared.cap.reads)]
starting.trem<-starting.points[starting.points[,1]%in%index,]
starting.trem<-starting.trem[starting.trem[,2]%in%colnames(shared.cap.reads),]
for (i in index)
{
	lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
	if(length(lilo)>1)
	{
	starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
	if(length(lilo)==1)
	{
	starting.trem[starting.trem[,1]==i,4]<-0
	starting.trem[starting.trem[,1]==i,5]<-lilo
	}
}
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Capnodiales",cex=1.3,font=4)
#
#
#
#
#
#
#
base.for.circos()
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.cap))
{	circos.polygon(c(0,reads.cap[i,1],reads.cap[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.cap[i,1],reads.cap[i,1]+reads.cap[i,2],reads.cap[i,1]+reads.cap[i,2],reads.cap[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}

#
# Capnodiales PLOT
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
#
for (i in 1:(length(index)-1))
{
	for (j in (i+1):length(index))
	{
	uno<-starting.trem[starting.trem[,1]==index[i],]
	dos<-starting.trem[starting.trem[,1]==index[j],]
	tres<-uno[,2][uno[,2]%in%dos[,2]]
	for (k in tres)
		{
			if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
			if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
			circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
		}
	rm(uno,dos,tres)
	}
}
text(-0.8,1,"Capnodiales",cex=1.3,font=4)
####
#END CIRCOS
#--------------------------------------

#------------------------------------
#
#15. Start Circos All
#
#------------------------------------
#
circos.initialize(1:25,xlim=c(0,1))
# Start text
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
#
j<-0
for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
{
circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
if (j!=0)
	{
	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
	}
	j<-i
	if (j==25)
		{
		circos.lines(c(0.1,1),c(0.4,0.4),j)
		}
}
for (i in c(1:25)[infected==1])
{
circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
}

#
#
# Start Plot connectors
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#------
#NEW Sector percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.foo))
{	circos.polygon(c(0,reads.foo[i,1],reads.foo[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.foo[i,1],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
starting.trem<-starting.points
index<-c(1:25)
	for (i in 1:25)
	{
		lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
		if(length(lilo)>1)
		{
		starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
		starting.trem[starting.trem[,1]==i,5]<-lilo
		}
		if(length(lilo)==1)
		{
		starting.trem[starting.trem[,1]==i,4]<-0
		starting.trem[starting.trem[,1]==i,5]<-lilo
		}
	}
	for (i in 1:(length(index)-1))
	{
		for (j in (i+1):length(index))
		{
		uno<-starting.trem[starting.trem[,1]==index[i],]
		dos<-starting.trem[starting.trem[,1]==index[j],]
		tres<-uno[,2][uno[,2]%in%dos[,2]]
		for (k in tres)
			{
				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
			}
		rm(uno,dos,tres)
		}
	}

text(-0.8,1,"Non_lichen_OTUs",cex=1.3,font=4)
#

#------------------------------------
#
#15.2 Second plot All
#
#------------------------------------
#
circos.initialize(1:25,xlim=c(0,1))
# Start text
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
#
j<-0
for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
{
circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
if (j!=0)
	{
	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
	}
	j<-i
	if (j==25)
		{
		circos.lines(c(0.1,1),c(0.4,0.4),j)
		}
}
for (i in c(1:25)[infected==1])
{
circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
}

#
#
# Start Plot connectors
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#------
#NEW Sector percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.foo))
{	circos.polygon(c(0,reads.foo[i,1],reads.foo[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.foo[i,1],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
	for (i in 1:(length(index)-1))
	{
		for (j in (i+1):length(index))
		{
		uno<-starting.trem[starting.trem[,1]==index[i],]
		dos<-starting.trem[starting.trem[,1]==index[j],]
		tres<-uno[,2][uno[,2]%in%dos[,2]]
		for (k in tres)
			{
				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
			}
		rm(uno,dos,tres)
		}
	}

text(-0.8,1,"Non_lichen_OTUs",cex=1.3,font=4)
#

#------------------------------------
#
#16. Start Circos Lichens
#
#------------------------------------
#seq.lichen<-seq.dataset[seq.dataset[,8]==1,]
#shared.lichen<-table(seq.lichen[,1],seq.lichen[,5])
#for (i in rownames(shared.lichen))
#{
#	for (j in colnames(shared.lichen))
#	{
#		if (shared.lichen[i,j]!=0)
#		{
#		shared.lichen[i,j]<-sum(as.numeric(seq.lichen[seq.lichen[,1]==i&seq.lichen[,5]==j,3]))
#		}
#	}
#}
##reorder and normalize
#shared.lichen<-shared.lichen[reorder,]
#for (i in rownames(shared.lichen)) shared.lichen[i,]<-shared.lichen[i,]/total.reads[i]
#rownames(shared.lichen)<-c(1:25)
#shared.lichen<-melt(shared.lichen)
#shared.lichen<-shared.lichen[shared.lichen[,3]!=0,]
#shared.lichen<-shared.lichen[order(shared.lichen[,1],shared.lichen[,2]),]
#shared.lichen<-cbind(shared.lichen,0,0)
#for (i in 1:25)
#{
#	lilo<-cumsum(shared.lichen[shared.lichen[,1]==i,3])
#	if(length(lilo)!=0)
#	{
#	shared.lichen[shared.lichen[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
#	shared.lichen[shared.lichen[,1]==i,5]<-lilo
#	}
#}
#
#
##
## START PLOT
##-------
##
#circos.initialize(1:25,xlim=c(0,1))
## Start text
##
#circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
##
#j<-0
#for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
#{
#circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
#if (j!=0)
#	{
#	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
#	}
#	j<-i
#	if (j==25)
#		{
#		circos.lines(c(0.1,1),c(0.4,0.4),j)
#		}
#}
#for (i in c(1:25)[infected==1])
#{
#circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
#}
###
#circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
##
##
##
## NEW
##
##
##
##
#
#for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#start<-mat.or.vec(1,25)
#starting.trem<-shared.lichen
#index<-c(1:25)
#	for (i in 1:25)
#	{
#		lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
#		if(length(lilo)>1)
#		{
#		starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
#		starting.trem[starting.trem[,1]==i,5]<-lilo
#		}
#		if(length(lilo)==1)
#		{
#		starting.trem[starting.trem[,1]==i,4]<-0
#		starting.trem[starting.trem[,1]==i,5]<-lilo
#		}
#	}
#	for (i in 1:(length(index)-1))
#	{
#		for (j in (i+1):length(index))
#		{
#		uno<-starting.trem[starting.trem[,1]==index[i],]
#		dos<-starting.trem[starting.trem[,1]==index[j],]
#		tres<-uno[,2][uno[,2]%in%dos[,2]]
#		for (k in tres)
#			{
#				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
#				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
#				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
#			}
#		rm(uno,dos,tres)
#		}
#	}
#
#text(-0.8,1,"Lichen_OTUs",cex=1.3,font=4)
#
##------------------------------------
##
##16.2 Second Circos Lichens
##
##------------------------------------
#seq.lichen<-seq.dataset[seq.dataset[,8]==1,]
#shared.lichen<-table(seq.lichen[,1],seq.lichen[,5])
#for (i in rownames(shared.lichen))
#{
#	for (j in colnames(shared.lichen))
#	{
#		if (shared.lichen[i,j]!=0)
#		{
#		shared.lichen[i,j]<-sum(as.numeric(seq.lichen[seq.lichen[,1]==i&seq.lichen[,5]==j,3]))
#		}
#	}
#}
##reorder and normalize
#shared.lichen<-shared.lichen[reorder,]
#for (i in rownames(shared.lichen)) shared.lichen[i,]<-shared.lichen[i,]/total.reads[i]
#rownames(shared.lichen)<-c(1:25)
#shared.lichen<-melt(shared.lichen)
#shared.lichen<-shared.lichen[shared.lichen[,3]!=0,]
#shared.lichen<-shared.lichen[order(shared.lichen[,1],shared.lichen[,2]),]
#shared.lichen<-cbind(shared.lichen,0,0)
#for (i in 1:25)
#{
#	lilo<-cumsum(shared.lichen[shared.lichen[,1]==i,3])
#	if(length(lilo)!=0)
#	{
#	shared.lichen[shared.lichen[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
#	shared.lichen[shared.lichen[,1]==i,5]<-lilo
#	}
#}
#
#
##
## START PLOT
##-------
##
#circos.initialize(1:25,xlim=c(0,1))
## Start text
##
#circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
##
#j<-0
#for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
#{
#circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
#if (j!=0)
#	{
#	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
#	}
#	j<-i
#	if (j==25)
#		{
#		circos.lines(c(0.1,1),c(0.4,0.4),j)
#		}
#}
#for (i in c(1:25)[infected==1])
#{
#circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
#}
###
#circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
##
##
##
## NEW
##
##
##
##
#
#for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#
#	for (i in 1:(length(index)-1))
#	{
#		for (j in (i+1):length(index))
#		{
#		uno<-starting.trem[starting.trem[,1]==index[i],]
#		dos<-starting.trem[starting.trem[,1]==index[j],]
#		tres<-uno[,2][uno[,2]%in%dos[,2]]
#		for (k in tres)
#			{
#				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
#				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
#				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
#			}
#		rm(uno,dos,tres)
#		}
#	}
#
#text(-0.8,1,"Lichen_OTUs",cex=1.3,font=4)
#
#
#
#
#END CIRCOS
#--------------------------------------
#------------------------------------
#
#17. Start Circos LICHENICOLOUS
#
#------------------------------------
#
circos.initialize(1:25,xlim=c(0,1))
# Start text
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
#
j<-0
for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
{
circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
if (j!=0)
	{
	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
	}
	j<-i
	if (j==25)
		{
		circos.lines(c(0.1,1),c(0.4,0.4),j)
		}
}
for (i in c(1:25)[infected==1])
{
circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
}

#
#
# Start Plot connectors
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#------
#NEW Sector percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.foo))
{	circos.polygon(c(0,reads.foo[i,1],reads.foo[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.foo[i,1],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
starting.trem<-starting.points[starting.points[,2]%in%c(21,34,18,32,17,46,33,10,12),]
index<-c(1:25)
	for (i in 1:25)
	{
		lilo<-cumsum(starting.trem[starting.trem[,1]==i,3])
		if(length(lilo)>1)
		{
		starting.trem[starting.trem[,1]==i,4]<-c(0,lilo[1:(length(lilo)-1)])
		starting.trem[starting.trem[,1]==i,5]<-lilo
		}
		if(length(lilo)==1)
		{
		starting.trem[starting.trem[,1]==i,4]<-0
		starting.trem[starting.trem[,1]==i,5]<-lilo
		}
	}
	for (i in 1:(length(index)-1))
	{
		for (j in (i+1):length(index))
		{
		uno<-starting.trem[starting.trem[,1]==index[i],]
		dos<-starting.trem[starting.trem[,1]==index[j],]
		tres<-uno[,2][uno[,2]%in%dos[,2]]
		for (k in tres)
			{
				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
			}
		rm(uno,dos,tres)
		}
	}

text(-0.8,1,"Lichenicolous fungi1",cex=1.3,font=4)
#

#------------------------------------
#
#17.2 Second plot LICHENICOLOUS FUNGI
#
#------------------------------------
#
circos.initialize(1:25,xlim=c(0,1))
# Start text
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.1,bg.col=NA,bg.border=NA)
#
j<-0
for (i in c(1,2,4,5,7,8,13,14,16,18,20,24,25))
{
circos.text(text.pos[i],0.8,sector=i, facing="bending.inside",cex=0.7, foo[i,1],niceFacing=F,font=4)
if (j!=0)
	{
	circos.lines(c(0.1,((i-j+(i-j)*0.1))),c(0.4,0.4),j)
	}
	j<-i
	if (j==25)
		{
		circos.lines(c(0.1,1),c(0.4,0.4),j)
		}
}
for (i in c(1:25)[infected==1])
{
circos.text(0.5,-0.1,sector=i, facing="bending.inside",cex=0.6, foo[i,2],niceFacing=F,font=3)
}

#
#
# Start Plot connectors
#
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.col=colorinos,bg.border=infected)
for (i in 1:25)	circos.text(0.5,0.5,sector=i, cex=0.8, facing="bending.inside", foo[i,3])
#------
#NEW Sector percentage of shared reads
#-------
circos.trackPlotRegion(ylim=c(0,1), track.height=0.05,bg.border=F)
for (i in rownames(reads.foo))
{	circos.polygon(c(0,reads.foo[i,1],reads.foo[i,1],0), c(0,0,1,1), sector.index = i,col=colorinos.trans[as.numeric(i)],border=F)
	circos.polygon(c(reads.foo[i,1],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]+reads.foo[i,2],reads.foo[i,1]), c(0,0,1,1), sector.index = i,col=colorinos[as.numeric(i)],border=F)
}
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
starting.trem<-starting.trem[starting.trem$value>foo.exclude[as.numeric(starting.trem$X1)],]
	for (i in 1:(length(index)-1))
	{
		for (j in (i+1):length(index))
		{
		uno<-starting.trem[starting.trem[,1]==index[i],]
		dos<-starting.trem[starting.trem[,1]==index[j],]
		tres<-uno[,2][uno[,2]%in%dos[,2]]
		for (k in tres)
			{
				if (uno[uno[,2]==k,3]>=dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[i]]}
				if (uno[uno[,2]==k,3]<dos[dos[,2]==k,3]) {colorin<-colorinos.trans[index[j]]}
				circos.link(index[i],as.numeric(as.vector(uno[uno[,2]==k,c(4,5)])),index[j],as.numeric(as.vector(dos[dos[,2]==k,c(4,5)])),col=colorin)		
			}
		rm(uno,dos,tres)
		}
	}

text(-0.8,1,"Lichenicolous fungi",cex=1.3,font=4)
#


dev.off()

#----------------------------------------------------------------------------#
#                                                                            #
#                      Script to analyse bipartite networks                  #
#                                 FFM  6.4.2017                              #
#                                                                            #
#----------------------------------------------------------------------------#

library(bipartite)

#----------------------------------------------------------------------------
# 1.Complete network
#----------------------------------------------------------------------------

eraseme<-starting.points
eraseme<-eraseme[,1:4]
eraseme[,4]<-1
eraseme$value<-as.numeric(eraseme$value)
#eraseme<-frame2webs(eraseme,varnames=c("Var1","Var2","0","value"))#eraseme<-frame2webs(eraseme,varnames=c("Var1","Var2","0","value"))
#eraseme<-frame2webs(eraseme,varnames=c("X1","X2","0","value"))
eraseme<-frame2webs(eraseme,varnames=c("X1","X2","out","value"))
# NOT USE eraseme.shared<-eraseme[[1]][,colSums(eraseme[[1]]>0)>1]
eraseme<-eraseme[[1]]
#----------------------------------------------------------------------------
# Create category names
#----------------------------------------------------------------------------
# The important
#----------------------------------------------------------------------------
names.foo.eraseme<-seq.subset[seq.subset[,5]%in%colnames(eraseme)&!duplicated(seq.subset[,5]),19]
subset.network.list<-cbind(seq.subset[seq.subset[,5]%in%colnames(eraseme)&!duplicated(seq.subset[,5]),5],names.foo.eraseme)
names.foo.eraseme<-gsub("Tremellales",colorinos[13],names.foo.eraseme)
names.foo.eraseme<-gsub("Botryosphaeriales",colorinos[3],names.foo.eraseme)
names.foo.eraseme<-gsub("Capnodiales",colorinos[17],names.foo.eraseme)
names.foo.eraseme<-gsub("Chaetothyriales",colorinos[6],names.foo.eraseme)
#----------------------------------------------------------------------------
# The rest
#----------------------------------------------------------------------------
names.foo.eraseme<-gsub("Verrucariales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Lecanorales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("unidentified", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Ostropales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Myriangiales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Pleosporales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Hypocreales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Peltigerales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Candelariales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Russulales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Rhizocarpales", colorinos[19],names.foo.eraseme)#
names.foo.eraseme<-gsub("Polyporales", colorinos[19],names.foo.eraseme)
#----------------------------------------------------------------------------
# The other rest
#----------------------------------------------------------------------------
names.foo.eraseme<-gsub("Helotiales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Cystofilobasidiales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Cantharellales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Taphrinales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Auriculariales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Diaporthales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Lecideales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Hymenochaetales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Agaricales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Umbilicariales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Pertusariales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Sebacinales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Eurotiales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Corticiales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Sporidiobolales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Saccharomycetales", colorinos[19],names.foo.eraseme)
names.foo.eraseme<-gsub("Blastocladiales", colorinos[19],names.foo.eraseme)
subset.network.list<-cbind(subset.network.list,names.foo.eraseme)
#----------------------------------------------------------------------------
# Plot bipartite network
#----------------------------------------------------------------------------
pdf("./Bipartite_networks.pdf",paper="a4r")
#----------------------------------------------------------------------------
# Calculate network metrics Complete network
#----------------------------------------------------------------------------
op<-networklevel(eraseme, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE)
#cop<-degreedistr(eraseme)
#bop<-betweenness_w (eraseme)
plotweb(eraseme,method="normal",ybig=1,col.low=colorinos,col.high=names.foo.eraseme,col.interaction=names.foo.eraseme,bor.col.interaction=0,bor.col.high=0)
title(main = "All non lichen OTUs",cex=1.3,font=4,line=-2)
visweb(eraseme,type="diagonal",circles=TRUE,circle.max=1,circle.min=0.1,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3])
title(main = "All non lichen OTUs",cex=1.3,font=4,line=-2)
#----------------------------------------------------------------------------
# Calculate network metrics Complete network with threshold
#----------------------------------------------------------------------------
eraseme.trem<-eraseme
for (i in 1:25) eraseme.trem[i,eraseme.trem[i,]<foo.exclude[i]]<-0
eraseme.trem<-eraseme.trem[,colSums(eraseme.trem)>0]
op<-cbind(op,op<-networklevel(eraseme.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(eraseme.trem))
#bop<-c(bop,betweenness_w (eraseme.trem))
#bipartite:::nestedcontribution(eraseme.trem, nsimul = 499)
plotweb(eraseme.trem,method="normal",ybig=1,col.low=colorinos[rowSums(eraseme.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(eraseme.trem),3],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(eraseme.trem),3],bor.col.interaction=0)
title(main = "All non lichen OTUs with Treshold",cex=1.3,font=4,line=-2)
visweb(eraseme.trem,type="diagonal",circles=TRUE,circle.max=1,circle.min=0.1,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3])
title(main = "All non lichen OTUs with Treshold",cex=1.3,font=4,line=-2)
#----------------------------------------------------------------------------
# Calculate network metrics for Tremellales only
#----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Tremellales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Tremellales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Tremellales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
for (i in 1:25) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Tremellales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Tremellales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Capnodiales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Capnodiales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Capnodiales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Capnodiales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:25) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Capnodiales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Capnodiales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Botryosphaeriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Botryosphaeriales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Botryosphaeriales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Botryosphaeriales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:25) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Botryosphaeriales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Botryosphaeriales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Chaetothyriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Chaetothyriales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Chaetothyriales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Chaetothyriales ",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:25) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op<-cbind(op,op<-networklevel(network.trem, index=c("number of species","connectance","links per species","number of compartments","nestedness","weighted nestedness","weighted NODF","mean number of shared partners","cluster coefficient", "weighted cluster coefficient","C score"),weighted=TRUE))
#cop<-c(cop,degreedistr(network.trem))
#bop<-c(bop,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Chaetothyriales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Chaetothyriales with Threshold",cex=1.3,font=4,line=-2)
dev.off()
write.table(op, "./OP2_Sample_network1s.txt")
#----------------------------------------------------------------------------#
#                                                                            #
#                      Script to analyse UNIPARTITE networks                 #
#                                 FFM  6.4.2017                              #
#                                                                            #
#----------------------------------------------------------------------------#

library(bipartite)
#-----------------------------------------------------
#
# Threshold criteria for networks
#
#-----------------------------------------------------
# Option 1 Averaged Threshold
# starting.trem<-starting.trem[starting.trem[,3]>0.0022,]
# Option 2 Specific Threshold per sample
exclude.lichens<-seq.dataset[liq.ds==1,]
exclude.lichens<-data.frame(exclude.lichens)
exclude.lichens$size<-as.numeric(as.matrix(exclude.lichens$size))
exclude.lichens<-aggregate(size~sample+OTU,exclude.lichens,sum)
names(exclude.lichens)[3]<-"size"
# Optional substract OTU 0
exclude.lichens<-exclude.lichens[exclude.lichens[,2]!=0,]
#
exclude.lichens<-aggregate(size~sample,exclude.lichens,mean)
#OTHER OPTIONS
#exclude.lichens<-aggregate(size~sample,exclude.lichens,sum)
#exclude.lichens<-aggregate(size~sample,exclude.lichens,max)
rownames(exclude.lichens)<-exclude.lichens[,1]
#for (i in rownames(exclude.lichens)) exclude.lichens[i,2]<-exclude.lichens[i,2]/total.reads[i]
for (i in rownames(exclude.lichens)) exclude.lichens[i,2]<-exclude.lichens[i,2]/subset.reads[i]
foo.exclude<-subset.reads
for (i in rownames(exclude.lichens)) foo.exclude[i]<-exclude.lichens[i,2]
foo.exclude[foo.exclude>1]<-0
foo.exclude[foo.exclude<0.001]<-0.001
#----------------------------------------------------------------------------
# 1.Complete network
#----------------------------------------------------------------------------

unipartite.all<-mat.or.vec(25,25)
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.points[starting.points$X1==i,][starting.points[starting.points$X1==i,2]%in%starting.points[starting.points$X1==j,2],3])
		}
}

for (i in 1:25)
{
		unipartite.all[i,i]<-0
		for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}

pipi<-cbind(PDI(unipartite.all),0)
unip<-networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE)
#----------------------------------------------------------------------------
# 1.Trimmed network
#----------------------------------------------------------------------------
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.points[starting.points$value>foo.exclude[starting.points$X1],]
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))
pipi<-rbind(pipi,cbind(PDI(unipartite.all),1))
boxplot(pipi[-c(26,27),1]~pipi[-c(26,27),2])
dev.off()
#----------------------------------------------------------------------------
# Calculate network metrics for Tremellales only
#----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Tremellales",]
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))
##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem[starting.trem$value>foo.exclude[starting.trem$X1],]
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))

##----------------------------------------------------------------------------
##Capnodiales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Capnodiales",]
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))

##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem[starting.trem$value>foo.exclude[starting.trem$X1],]
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))

##----------------------------------------------------------------------------
##Botryosphaeriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Botryosphaeriales",]
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))

##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem[starting.trem$value>foo.exclude[starting.trem$X1],]
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))

##----------------------------------------------------------------------------
##Chaetothyriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Chaetothyriales",]
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))
##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
starting.trem<-starting.points[starting.points[,2]%in%subset.network.trem[,1],]
unipartite.all<-mat.or.vec(25,25)
starting.trim<-starting.trem[starting.trem$value>foo.exclude[starting.trem$X1],]
for (i in 1:25)
{
	for (j in 1:25)
		{
			unipartite.all[i,j]<-sum(starting.trim[starting.trim$X1==i,][starting.trim[starting.trim$X1==i,2]%in%starting.trim[starting.trim$X1==j,2],3])
		}
}
for (i in 1:25)
{
	unipartite.all[i,i]<-0
	for (j in 1:25)
		{
			unipartite.all[i,j]<-unipartite.all[j,i]<-unipartite.all[i,j]*unipartite.all[j,i]
		}
}
unip<-cbind(unip,networklevel(unipartite.all, index=c("number of species","connectance","links per species","number of compartments","mean number of shared partners"),weighted=TRUE))
write.table(unip, "./OP2_Sample_Unipartite_networks.txt")
#
##----------------------------------------------------------------------------#
##                                                                            #
##                      Script to analyse bipartite networks at species level #
##                                 FFM  6.4.2017                              #
##                                                                            #
##----------------------------------------------------------------------------#
eraseme<-starting.points
eraseme<-eraseme[,1:4]
eraseme[,4]<-1
eraseme$value<-as.numeric(eraseme$value)
eraseme$X1<-gsub(10,"L_polytropa",eraseme$X1)
eraseme$X1<-gsub(11,"L_polytropa",eraseme$X1)
eraseme$X1<-gsub(12,"L_polytropa",eraseme$X1)
eraseme$X1<-gsub(13,"L_swartzii",eraseme$X1)
eraseme$X1<-gsub(14,"L_lapicida",eraseme$X1)
eraseme$X1<-gsub(15,"L_lapicida",eraseme$X1)
eraseme$X1<-gsub(16,"P_conglomerata",eraseme$X1)
eraseme$X1<-gsub(17,"P_conglomerata",eraseme$X1)
eraseme$X1<-gsub(18,"R_geographicum",eraseme$X1)
eraseme$X1<-gsub(19,"R_geographicum",eraseme$X1)
eraseme$X1<-gsub(20,"T_atra",eraseme$X1)
eraseme$X1<-gsub(21,"T_atra",eraseme$X1)
eraseme$X1<-gsub(22,"T_atra",eraseme$X1)
eraseme$X1<-gsub(23,"T_atra",eraseme$X1)
eraseme$X1<-gsub(24,"U_cylindrica",eraseme$X1)
eraseme$X1<-gsub(25,"V_lactea",eraseme$X1)
eraseme$X1<-gsub(1,"A_fuscata",eraseme$X1)
eraseme$X1<-gsub(2,"A_myrinii",eraseme$X1)
eraseme$X1<-gsub(3,"A_myrinii",eraseme$X1)
eraseme$X1<-gsub(4,"Ca_vitellina",eraseme$X1)
eraseme$X1<-gsub(5,"L_bicincta",eraseme$X1)
eraseme$X1<-gsub(6,"L_bicincta",eraseme$X1)
eraseme$X1<-gsub(7,"L_intricata",eraseme$X1)
eraseme$X1<-gsub(8,"L_polytropa",eraseme$X1)
eraseme$X1<-gsub(9,"L_polytropa",eraseme$X1)

eraseme<-frame2webs(eraseme,varnames=c("X1","X2","out","value"))
#eraseme<-frame2webs(eraseme,varnames=c("X1","X2","0","value"))

## NOT USE eraseme.shared<-eraseme[[1]][,colSums(eraseme[[1]]>0)>1]
eraseme<-eraseme[[1]]
#
#
##----------------------------------------------------------------------------
## Plot bipartite networks
##----------------------------------------------------------------------------
pdf("./Bipartite_networks_species.pdf",paper="a4r")
foo.exclude<-as.numeric(c(foo.exclude[1],mean(as.numeric(foo.exclude[2:3])),foo.exclude[4],mean(as.numeric(foo.exclude[5:6])),foo.exclude[7],mean(as.numeric(foo.exclude[8:12])),foo.exclude[13],mean(as.numeric(foo.exclude[14:15])),mean(as.numeric(foo.exclude[16:17])),mean(as.numeric(foo.exclude[18:19])),mean(as.numeric(foo.exclude[20:23])),foo.exclude[24],foo.exclude[25]))
##----------------------------------------------------------------------------
## Calculate network metrics Complete network
##----------------------------------------------------------------------------
d.complete.1<-dfun(eraseme)
op2<-networklevel(eraseme)
#cop2<-degreedistr(eraseme)
#bop2<-betweenness_w (eraseme)
plotweb(eraseme,method="cca",ybig=1,col.low=colorinos,col.high=names.foo.eraseme,col.interaction=names.foo.eraseme,bor.col.interaction=0,bor.col.high=0)
title(main = "All non lichen OTUs",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
## Calculate network metrics Complete network with threshold
##----------------------------------------------------------------------------
eraseme.trem<-eraseme
for (i in 1:13) eraseme.trem[i,eraseme.trem[i,]<foo.exclude[i]]<-0
eraseme.trem<-eraseme.trem[,colSums(eraseme.trem)>0]
op2<-cbind(op2,networklevel(eraseme.trem))
d.complete.2<-dfun(eraseme.trem)
#cop2<-c(cop2,degreedistr(eraseme.trem))
#bop2<-c(bop2,betweenness_w (eraseme.trem))
#bipartite:::nestedcontribution(eraseme.trem, nsimul = 99)
plotweb(eraseme.trem,method="cca",ybig=1,col.low=colorinos[rowSums(eraseme.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(eraseme.trem),3],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(eraseme.trem),3],bor.col.interaction=0)
title(main = "All non lichen OTUs with Treshold",cex=1.3,font=4,line=-2)
visweb(eraseme.trem,type="diagonal",circles=TRUE,circle.max=1,circle.min=0.1,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3])
title(main = "All non lichen OTUs with Treshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
## Calculate network metrics for Tremellales only
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Tremellales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Tremellales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Tremellales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
## Threshold
##----------------------------------------------------------------------------
for (i in 1:13) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Tremellales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Tremellales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Capnodiales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Capnodiales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Capnodiales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Capnodiales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:13) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Capnodiales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Capnodiales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Botryosphaeriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Botryosphaeriales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Botryosphaeriales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Botryosphaeriales",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:13) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Botryosphaeriales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Botryosphaeriales with Threshold",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Chaetothyriales
##----------------------------------------------------------------------------
subset.network.trem<-subset.network.list[subset.network.list[,2]=="Chaetothyriales",]
network.trem<-eraseme[,colnames(eraseme)%in%subset.network.trem[,1]]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Chaetothyriales",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Chaetothyriales ",cex=1.3,font=4,line=-2)
##----------------------------------------------------------------------------
##Threshold
##----------------------------------------------------------------------------
for (i in 1:13) network.trem[i,network.trem[i,]<foo.exclude[i]]<-0
network.trem<-network.trem[,colSums(network.trem)>0]
op2<-cbind(op2,networklevel(network.trem))
#cop2<-c(cop2,degreedistr(network.trem))
#bop2<-c(bop2,betweenness_w (network.trem))
#bipartite:::nestedcontribution(network.trem, nsimul = 499)
plotweb(network.trem,method="normal",ybig=1,col.low=colorinos[rowSums(network.trem)!=0],col.high=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],col.interaction=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1],bor.col.interaction=0)
title(main = "Chaetothyriales with Threshold",cex=1.3,font=4,line=-2)
visweb(network.trem,type="diagonal",circles=TRUE,circle.max=2,circle.min=0.2,circle.col=subset.network.list[subset.network.list[,1]%in%colnames(network.trem),3][1])
title(main = "Chaetothyriales with Threshold",cex=1.3,font=4,line=-2)
dev.off()
write.table(op2, "./OP2_Species_network1s.txt")
##------------------------------------
##
## d for lichenicoles 
##
##------------------------------------
network.lichenicoles<-eraseme[,c("21","34","18","32","17","46","33","10","12")]
d.complete.3<-dfun(network.lichenicoles)
d.complete.5<-dfun(t(network.lichenicoles))
for (i in 1:13)
	{
	network.lichenicoles[i,network.lichenicoles[i,]<foo.exclude[i]]<-0
	}
d.complete.4<-dfun(network.lichenicoles)
d.complete.6<-dfun(t(network.lichenicoles))

#
# All together
#
cbind(
d.complete.1$dprime,d.complete.1$d,d.complete.1$dmin,d.complete.1$dmax,
d.complete.2$dprime,d.complete.2$d,d.complete.2$dmin,d.complete.2$dmax,
d.complete.3$dprime[names(d.complete.1[[1]])],d.complete.3$d[names(d.complete.1[[1]])],d.complete.3$dmin[names(d.complete.1[[1]])],d.complete.3$dmax[names(d.complete.1[[1]])],
d.complete.4$dprime[names(d.complete.1[[1]])],d.complete.4$d[names(d.complete.1[[1]])],d.complete.4$dmin[names(d.complete.1[[1]])],d.complete.4$dmax[names(d.complete.1[[1]])])