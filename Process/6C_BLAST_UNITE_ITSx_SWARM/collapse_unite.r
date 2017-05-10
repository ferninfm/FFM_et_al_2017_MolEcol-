#!/usr/bin/env Rscript
library(ggplot2)
pdf("~/Dropbox/ANTONIA/2_Not_collapsed/6C_BLAST_UNITE_ITSx_SWARM/boxplots_corrected1.pdf", onefile = TRUE)
category.names<-c("query", "size", "ref", "Kingdom", "Division", "Class", "Order", "Family", "Genus", "Pseudospecies")
load("~/Dropbox/ANTONIA/2_Not_collapsed/6C_BLAST_UNITE_ITSx_SWARM/seq_otus.Rdata")
foo<-read.table("~/Dropbox/ANTONIA/2_Not_collapsed/6C_BLAST_UNITE_ITSx_SWARM/ITS1_OTUs_Blast_Processed.txt",sep=";",header=T)
unite_min1<-foo[!duplicated(foo[,1]),]
unite_mean<-unite_min<-unite_min1[,-c(18:20)]
unite_mean[,4:17]<-0
unite_min[,4:17]<-0
cuboti<-cuboti[as.character(unite_min1[,1]),]
for (file in 1:dim(unite_min1)[1])
	{
	unite_otus<-foo[foo[,1]==unite_min1[file,1],]
	if (dim(unite_otus)[1]==1)
		{
			unite_otus[19]<-log(unite_otus[19])
			unite_min[file,4:17]<-unite_mean[file,4:17]<-as.matrix(unite_otus[c(4,19,5,19,6,19,7,19,8,19,9,19,10,19)])
		}
	if (dim(unite_otus)[1]>1)
		{
		for (category in c(4:10))
			{
				p1<-ggplot(unite_otus,aes(x=factor(unite_otus[,category]),y=log(as.numeric(unite_otus[,19])),colour=unite_otus[,category])) + geom_boxplot(outlier.size = 0) + geom_jitter() + theme(axis.text.x=element_text(angle=90)) + guides(colour=FALSE)+ labs(x="" , y="log E value", title=paste ("Ascription of OTU",cuboti[file,2],"to",category.names[category],sep=" "))
				print(p1)
				foo_mean_E<-as.character(unite_otus[!duplicated(unite_otus[,category]),category])
				foo_mean_E<-foo_mean_E[!is.na(foo_mean_E)]
				names(foo_mean_E)<-foo_mean_E
				cats.foo<-foo_min_E<-foo_mean_E
				for (i in names(foo_mean_E))
					{
						if (sum(unite_otus[,category]==i,na.rm=T)==1)
						{
							foo_mean_E[i]<-log(as.numeric(unite_otus[unite_otus[,category]==i,19]))[is.na(log(unite_otus[unite_otus[,category]==i,19]))==F]
							foo_min_E[i]<-foo_mean_E[i]
						}
						if (sum(unite_otus[,category]==i,na.rm=T)!=1)
						{
							foo_mean_E[i]<-mean(log(as.numeric(unite_otus[unite_otus[,category]==i,19])),na.rm=T)
							foo_min_E[i]<-min(log(as.numeric(unite_otus[unite_otus[,category]==i,19])),na.rm=T)				
						}
					}
					foo_mean_E<-as.numeric(foo_mean_E)
					foo_min_E<-as.numeric(foo_min_E)
					names(foo_min_E)<-names(foo_mean_E)<-cats.foo
					rm(cats.foo)
				if (length(foo_mean_E)==1)
						{
						unite_mean[file,((category-3)*2)+2]<-names(foo_mean_E)
					    unite_mean[file,((category-3)*2)+3]<-foo_mean_E
						unite_min[file,((category-3)*2)+2]<-names(foo_min_E)
					    unite_min[file,((category-3)*2)+3]<-foo_min_E
					 	}						
				if (length(foo_mean_E)!=1)
					{
					foo_mean_E<-foo_mean_E[!(names(foo_mean_E)%in%c("0","uncultured","unidentified","unknown"))]
					foo_min_E<-foo_min_E[!(names(foo_min_E)%in%c("0","uncultured","unidentified","unknown"))]
					unite_mean[file,((category-3)*2)+2]<-names(sort(foo_mean_E,decreasing=F)[1])
				    unite_mean[file,((category-3)*2)+3]<-sort(foo_mean_E,decreasing=F)[1]
					unite_min[file,((category-3)*2)+2]<-names(sort(foo_min_E,decreasing=F)[1])
				    unite_min[file,((category-3)*2)+3]<-sort(foo_min_E,decreasing=F)[1]
				 	}
			}
		}
	}
rm(category,file,foo,foo_mean_E,foo_min_E,unite_otus,i)
names(unite_mean)<-names(unite_min)<-c("query", "size", "ref", "kingdom", "e_k", "division", "e_d", "classe", "e_c", "order", "e_o", "family","e_f", "genus", "e_g","species","e_s")
#
mv <- function (a, b) { 
		anm <- deparse(substitute(a)) 
     	bnm <- deparse(substitute(b)) 
     	if (!exists(anm,where=1,inherits=FALSE)) 
         	stop(paste(anm, "does not exist.\n")) 
     	if (exists(bnm,where=1,inherits=FALSE)) { 
         	ans <- readline(paste("Overwrite ", bnm, "? (y/n) ", sep =   "")) 
         	if (ans != "y") 
             	return(invisible()) 
     	} 
     	assign(bnm, a, pos = 1) 
     	rm(list = anm, pos = 1) 
     	invisible() 
        } 
mv(unite_mean,unite_mean_new)
mv(unite_min,unite_min_new)
mv(unite_min1,unite_min1_new)
#save.image(file="~/Desktop/ANTONIA/2_Not_collapsed/6C_BLAST_UNITE_ITSx_SWARM/Workspace_full.Rdata")
dev.off()