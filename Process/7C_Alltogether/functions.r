#
#
# PRINT STACKED BARCHARTS
#
#
tabulate.and.print<-function(data=seq.subset,samples=1,groups=19,counts=3, category="Class",titulo=paste("Ascription of seqences to",category, sep=" "),n.cols=3, normalize=TRUE,reorder=FALSE,niveles,flip=FALSE,legend.pos="right",print=TRUE,extra=scale_fill_discrete())
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
if (normalize==TRUE)
	{
	table.clases[,]<-normalize.rows(table.clases)
	}
dc<-data.frame(table.clases)
dc<-dc[dc$Freq!=0,]#$
if (reorder==TRUE)
{
dc$Var2<-factor(dc$Var2, levels=niveles) 
}
if (print==TRUE)
	{
		require(ggplot2)
		if (flip==TRUE)
		{
		p1<-ggplot(dc, aes(x = factor(Var1), y=Freq, fill=factor(Var2), order=as.numeric(factor(Var2)))) + 		geom_bar(stat="identity",position=("stack")) + coord_flip() + theme(legend.position=legend.pos, legend.text=element_text(size=6)) + guides(fill=guide_legend(ncol=n.cols)) + labs(x="" , y="", title=titulo) + extra
		}else{
		if (flip==FALSE)
		{
		p1<-ggplot(dc, aes(x = factor(Var1), y=Freq, fill=factor(Var2),order=as.numeric(factor(Var2)))) + 	geom_bar(stat="identity",position=("stack"))+ theme(axis.text.x=element_text(angle=90), 	legend.position=legend.pos, legend.text=element_text(size=6)) + guides(fill=guide_legend(ncol=n.cols))+ 	labs(x="" , y="", title=titulo) + extra
			}
		}
	}
	if (print==FALSE)
		{
			p1<-NULL
		}	
return(list(table.clases,p1))
}
#
#
# PRINT PIE CHARTS
#
#
tabulate.and.pie<-function(data=seq.subset,samples=1,groups=19,counts=3, category="Class",titulo=paste("Ascription of seqences to",category, sep=" "),n.cols=3, normalize=TRUE,flip=FALSE, extra=scale_fill_manual(values=c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")),reorder=FALSE,niveles)
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
if (normalize==TRUE)
	{
	table.clases[,]<-normalize.rows(table.clases)
	}
dc<-data.frame(table.clases)
dc<-dc[dc$Freq!=0,]#$
if (reorder==TRUE)
{
dc$Var2<-factor(dc$Var2, levels=niveles) 
}
require(ggplot2)
p1<-ggplot(dc,aes(x = factor(Var1), y=Freq,fill=as.factor(Var2), order=as.numeric(factor(Var2))))+ geom_bar(stat="identity", width=1000) + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="right", legend.title=element_blank(), legend.text=element_text(size=6)) + facet_wrap(~ Var1) + extra + coord_polar(theta = "y")
return(list(table.clases,p1))
}
#
#
# MULTIPLOT FUNCTION
#
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#
#EXTRACT LEGEND
#
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}
#
# COLORES FER
#
#
colores.fer<-function(x=dc1$Var2,y=dc1$Freq,others="Other Orders",ui="unidentified")
{
	table.foo<-table(as.character(x))
	for (i in names(table.foo))
	{
		table.foo[i]<-sum(y[x==i])
	}
	table.foo<-names(sort(table.foo,decreasing=T))
	table.foo<-c(table.foo[table.foo!=ui&table.foo!=others],others, ui)
	n.colors<-length(table.foo)
	#
	col1<-c("#12309E","#6078D0","#3451B9","#0C2274","#051345")
	col2<-c("#E9A200","#FFD063","#FFC030","#AB7700","#664700")
	col3<-c("#E70006","#FE6367","#FD2F35","#A90004","#650003")
	col4<-c("#09BC00","#5FE358","#31D528","#078A00","#045300")
	col5<-c("#61099C","#9F58CE","#7E2BB7","#470673","#2A0245")
	grey<-c("#999999","#F0F0F0","#A8A8A8","#707070","#505050")
	#
	first<-c(col1[1],col2[1],col1[2],col2[2],col2[4],col1[4])
	second<-c(col1[1],col2[1],col3[1],col1[2],col2[2],col3[2],col1[4],col2[4],col3[4],col1[3],col2[3],col3[3])
	third<-c(col1[1],col2[1],col3[1],col4[1],col1[2],col2[2],col3[2],col4[2],col1[4],col2[4],col3[4],col4[4],col1[3],col2[3],col3[3],col4[3])
	fourth<-c(col1[1],col2[1],col3[1],col4[1],col1[2],col2[2],col3[2],col4[2],col1[4],col2[4],col3[4],col4[4],col1[3],col2[3],col3[3],col4[3],grey[3],grey[4],grey[5])
	fifth<-c(col1[1],col2[1],col3[1],col4[1],col5[1],col1[2],col2[2],col3[2],col4[2],col5[2],col1[4],col2[4],col3[4],col4[4],col5[4],col1[3],col2[3],col3[3],col4[3],col5[3],col1[5],col2[5],col3[5],col4[5],col5[5],grey[3],grey[4],grey[5])
	
	if(n.colors==2)
	{
		out<-c(col1[1],col2[1])
		}else{
			if(n.colors%in%c(3:7))
				{
				out<-c(first[1:(n.colors-1)],grey[1])			
				}else{
					if(n.colors%in%c(8:14))
					{
					out<-c(second[1:(n.colors-1)],grey[2],grey[1])
					}else{
						if(n.colors%in%c(15:18))
						{
						out<-c(third[1:(n.colors-1)],grey[2],grey[1])
						}else{
							if(n.colors%in%c(19:21))
								{
								out<-c(fourth[1:(n.colors-1)],grey[2],grey[1])
								}else{
									if(n.colors%in%c(22:30))
										{
										out<-c(fifth[1:(n.colors-1)],grey[2],grey[1])
										}else{
											out<-palette(rainbow(n.colors))
										}
									}
								}
							}
						}
					}
names(out)<-table.foo
return(out)
}
#
##This function creates true individual based rarefaction curves				###
##Author: Joshua Jacobs   email@joshuajacobs.org  www.joshuajacobs.org/R/rarefaction.html	###
##Usage rarefaction(x,y), x is the community matrix with species as columns			###

rarefaction<-function(x,subsample=5, symbol=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
	{
		require(vegan)
		options(warn=-1)
		x <-as.matrix(x)
		y1<-rowSums(x)
        rare.data<-x                                   
		select<-unique(sort(c((apply(x, 1, sum)), (seq(0,(max(y1)), by=subsample)), recursive=TRUE))) # NUMBER OF SUBSAMPLES
		# PREPARE OUTPUT
        storesummary.e<-matrix(data=NA, ncol=length(rare.data[,1]),nrow=length(select))
        rownames(storesummary.e)<-c(select)
        colnames(storesummary.e)<-rownames(x)
        storesummary.se<-matrix(data=NA, ncol=length(rare.data[,1]),nrow=length(select))
        rownames(storesummary.se)<-c(select)
        colnames(storesummary.se)<-rownames(x)
		#
		for(i in 1:length(select))  #the for loop
               	{
               	select.c<-select[i]  #assigns the 'i'th element of select to select.c
				foo<-rarefy(x,select.c, se=T)  #use whatever vegan fn you want
				storesummary.e[i,]<-foo[1,]
				storesummary.se[i,]<-foo[2,]                
               	}
		storesummary.e<-as.data.frame(storesummary.e)               
        richness.error<<-storesummary.se      
		for (i in 1:(length(storesummary.e)))
			{
				storesummary.e[,i]<-ifelse(select>sum(x[i,]), NA, storesummary.e[,i])
			}
		
		alpha<-list("richness"= storesummary.e, "SE"=richness.error, "subsample"=select)
		return(alpha)        
}
        