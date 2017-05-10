#!/usr/bin/env Rscript
library(ggplot2)
#pdf("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/boxplots1_fin_2248.pdf", onefile = TRUE)
#values_mean<-mat.or.vec(24960,14)
#values_min<-mat.or.vec(24960,14)
#otu_names<-mat.or.vec(24960,2)
#for (file in c(1:2248))
#{
#	if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size==0)
#	{
#		cat (paste("The file OTU_",file,".txt no tiene datos",sep=""))
#	}
#	if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size>0) #Check file size>0
#	{
#		unite_otus<-read.table(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""),sep=";")
#		category.names<-c("Query","Size","Ref","Kingdoms","Divisions","Classes","Orders","Families","Genera","Species")
#		plots<-list()
#		otu_names[file,1]<-unite_otus[1,1]
#		otu_names[file,2]<-unite_otus[1,2]
#		for (category in c(4:10))
#		{
#		p1<-ggplot(unite_otus,aes(x=factor(unite_otus[,category]),y=log(unite_otus[,19]),colour=unite_otus[,category])) + geom_boxplot(outlier.size = 0) + geom_jitter() + theme(axis.text.x=element_text(angle=90)) + guides(colour=FALSE)+ labs(x="" , y="log E value", title=paste ("Adscription of OTU",file,"to",category.names[category],sep=" "))
#		print(p1)
#		foo_mean_E<-mat.or.vec(length(levels(unite_otus[,category])),1)
#		names(foo_mean_E)<-levels(unite_otus[,category])
#		foo_min_E<-foo_mean_E
#		for (i in names(foo_mean_E))
#			{
#				if (sum(unite_otus[,category]==i,na.rm=T)==1)
#				{
#					foo_mean_E[i]<-log(unite_otus[unite_otus[,category]==i,19])[is.na(log(unite_otus[unite_otus[,category]==i,19]))==F]
#					foo_min_E[i]<-foo_mean_E[i]
#				}
#				if (sum(unite_otus[,category]==i,na.rm=T)!=1)
#				{
#					foo_mean_E[i]<-mean(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)
#					foo_min_E[i]<-min(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)				
#				}
#			}
#			if (length(foo_mean_E)!=1)
#				{
#					if (length(foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))])!=0)
#							{
#								foo_mean_E<-foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))]
#								foo_min_E<-foo_min_E[-which(names(foo_min_E)%in%c("uncultured","NA","unidentified","unknown"))]
#							}
#				}
#			foo_bar<-c(names(foo_mean_E)[foo_mean_E==min(foo_mean_E)],min(foo_mean_E))	
#			if (length(foo_bar)==2)
#				{
#					values_mean[file,((category-3)*2)-1]<-foo_bar[1]
#					values_mean[file,(category-3)*2]<-foo_bar[2]
#				}
#			if (length(foo_bar)>2)
#				{
#					values_mean[file,(category-3)*2]<-foo_bar[length(foo_bar)]
#					values_mean[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
#				}
#			foo_bar<-c(names(foo_min_E)[foo_min_E==min(foo_min_E)],min(foo_min_E))	
#			if (length(foo_bar)==2)
#				{
#					values_min[file,((category-3)*2)-1]<-foo_bar[1]
#					values_min[file,(category-3)*2]<-foo_bar[2]
#					
#				}
#			if (length(foo_bar)>2)
#				{
#					values_min[file,(category-3)*2]<-foo_bar[length(foo_bar)]
#					values_min[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
#				}
#	}
#	}
#}
#dev.off()
#save(values_mean,values_min,otu_names,file="~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Workspace_singletons.Rdata")
load ("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Workspace_singletons.Rdata")
pdf("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/boxplots_the_rest1.pdf", onefile = TRUE)
for (file in c(2249:10000))
{
	if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size==0|file.exists(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))==F)
	{
		cat (paste("The file OTU_",file,".txt no tiene datos o no existe",sep=""))
	}
	if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size>0) #Check file size>0
	{
		unite_otus<-read.table(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""),sep=";")
		category.names<-c("Query","Size","Ref","Kingdoms","Divisions","Classes","Orders","Families","Genera","Species")
		plots<-list()
		otu_names[file,1]<-unite_otus[1,1]
		otu_names[file,2]<-unite_otus[1,2]
		for (category in c(4:10))
		{
		p1<-ggplot(unite_otus,aes(x=factor(unite_otus[,category]),y=log(unite_otus[,19]),colour=unite_otus[,category])) + geom_boxplot(outlier.size = 0) + geom_jitter() + theme(axis.text.x=element_text(angle=90)) + guides(colour=FALSE)+ labs(x="" , y="log E value", title=paste ("Adscription of OTU",file,"to",category.names[category],sep=" "))
		print(p1)
		foo_mean_E<-mat.or.vec(length(levels(unite_otus[,category])),1)
		names(foo_mean_E)<-levels(unite_otus[,category])
		foo_min_E<-foo_mean_E
		for (i in names(foo_mean_E))
			{
				if (sum(unite_otus[,category]==i,na.rm=T)==1)
				{
					foo_mean_E[i]<-log(unite_otus[unite_otus[,category]==i,19])[is.na(log(unite_otus[unite_otus[,category]==i,19]))==F]
					foo_min_E[i]<-foo_mean_E[i]
				}
				if (sum(unite_otus[,category]==i,na.rm=T)!=1)
				{
					foo_mean_E[i]<-mean(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)
					foo_min_E[i]<-min(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)				
				}
			}
			if (length(foo_mean_E)!=1)
				{
					if (length(foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))])!=0)
							{
								foo_mean_E<-foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))]
								foo_min_E<-foo_min_E[-which(names(foo_min_E)%in%c("uncultured","NA","unidentified","unknown"))]
							}
				}
			foo_bar<-c(names(foo_mean_E)[foo_mean_E==min(foo_mean_E)],min(foo_mean_E))	
			if (length(foo_bar)==2)
				{
					values_mean[file,((category-3)*2)-1]<-foo_bar[1]
					values_mean[file,(category-3)*2]<-foo_bar[2]
				}
			if (length(foo_bar)>2)
				{
					values_mean[file,(category-3)*2]<-foo_bar[length(foo_bar)]
					values_mean[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
				}
			foo_bar<-c(names(foo_min_E)[foo_min_E==min(foo_min_E)],min(foo_min_E))	
			if (length(foo_bar)==2)
				{
					values_min[file,((category-3)*2)-1]<-foo_bar[1]
					values_min[file,(category-3)*2]<-foo_bar[2]
					
				}
			if (length(foo_bar)>2)
				{
					values_min[file,(category-3)*2]<-foo_bar[length(foo_bar)]
					values_min[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
				}
	}
	}
}
save(values_mean,values_min,otu_names,file="~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Workspace_full1.Rdata")
dev.off()
#pdf("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/boxplots_the_rest2.pdf", onefile = TRUE)
for (file in c(10001:24960))
{
	if (file.exists(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))==F)
		{
			print (paste("The file OTU_",file,".txt does not exist ",sep=""))
		}
	if (file.exists(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))==T)
		{
		if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size==0)
		{
			print (paste("The file OTU_",file,".txt has no data",sep=""))
		}
		if (file.info(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""))$size>0) #Check file size>0
		{
			unite_otus<-read.table(paste("~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/OTU_",file,".txt",sep=""),sep=";")
			category.names<-c("Query","Size","Ref","Kingdoms","Divisions","Classes","Orders","Families","Genera","Species")
			plots<-list()
			otu_names[file,1]<-unite_otus[1,1]
			otu_names[file,2]<-unite_otus[1,2]
			for (category in c(4:10))
			{
			#p1<-ggplot(unite_otus,aes(x=factor(unite_otus[,category]),y=log(unite_otus[,19]),colour=unite_otus[,category])) + geom_boxplot(outlier.size = 0) + geom_jitter() + theme(axis.text.x=element_text(angle=90)) + guides(colour=FALSE)+ labs(x="" , y="log E value", title=paste ("Adscription of OTU",file,"to",category.names[category],sep=" "))
			#print(p1)
			foo_mean_E<-mat.or.vec(length(levels(unite_otus[,category])),1)
			names(foo_mean_E)<-levels(unite_otus[,category])
			foo_min_E<-foo_mean_E
			for (i in names(foo_mean_E))
				{
					if (sum(unite_otus[,category]==i,na.rm=T)==1)
					{
						foo_mean_E[i]<-log(unite_otus[unite_otus[,category]==i,19])[is.na(log(unite_otus[unite_otus[,category]==i,19]))==F]
						foo_min_E[i]<-foo_mean_E[i]
					}
					if (sum(unite_otus[,category]==i,na.rm=T)!=1)
					{
						foo_mean_E[i]<-mean(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)
						foo_min_E[i]<-min(log(unite_otus[unite_otus[,category]==i,19]),na.rm=T)				
					}
				}
				if (length(foo_mean_E)!=1)
					{
						if (length(foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))])!=0)
								{
									foo_mean_E<-foo_mean_E[-which(names(foo_mean_E)%in%c("uncultured","NA","unidentified","unknown"))]
									foo_min_E<-foo_min_E[-which(names(foo_min_E)%in%c("uncultured","NA","unidentified","unknown"))]
								}
					}
				foo_bar<-c(names(foo_mean_E)[foo_mean_E==min(foo_mean_E)],min(foo_mean_E))	
				if (length(foo_bar)==2)
					{
						values_mean[file,((category-3)*2)-1]<-foo_bar[1]
						values_mean[file,(category-3)*2]<-foo_bar[2]
					}
				if (length(foo_bar)>2)
					{
						values_mean[file,(category-3)*2]<-foo_bar[length(foo_bar)]
						values_mean[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
					}
				foo_bar<-c(names(foo_min_E)[foo_min_E==min(foo_min_E)],min(foo_min_E))	
				if (length(foo_bar)==2)
					{
						values_min[file,((category-3)*2)-1]<-foo_bar[1]
						values_min[file,(category-3)*2]<-foo_bar[2]
						
					}
				if (length(foo_bar)>2)
					{
						values_min[file,(category-3)*2]<-foo_bar[length(foo_bar)]
						values_min[file,((category-3)*2)-1]<-paste(foo_bar[-length(foo_bar)],collapse="/")
					}
				}
		}
	}
}
save(values_mean,values_min,otu_names,file="~/Desktop/ANTONIA/2_Not_collapsed/6B_Blast_OTUs_UNITE/Workspace_full.Rdata")
#dev.off()