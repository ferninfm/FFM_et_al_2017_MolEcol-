#Read seqs per sample
library(ape)
out<-NULL
for (file in c("A032","A138","A172","A194","A227","A229","A243","A280","A360","A361","A368","A405","A418","A420","A434","A440","A476","A482","A608","A622","A623","A636","A670","A792","A809","A832"))
{
foo<-read.dna(paste("/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/1_Input_unique_filtered_fasta/",file,".unique.filtered.fa",sep=""),format="fasta")
foo<-cbind(names(foo),file)
out<-rbind(out,foo)
}
#
# Get a single reference per cluster and grep it to the cluster table
#
otus<-readLines("~/Desktop/ANTONIA/2_Not_collapsed/5B_Cluster_UCLUST/uclust_results.txt")
out2<-out
for (i in 1:dim(out2)[1])
{
out2[i,1]<-grep(unlist(strsplit(out[i]," "))[1],otus)
}
# Process the cluster table
parsed.otus<-mat.or.vec(length(otus),7)
colnames(parsed.otus)<-c("reference","size","sequences","classification","identity","OTU_A","OTU_B")
for (i in 1:length(otus))
	{
	foo<-unlist(strsplit(otus[[i]],"\t"))
	if (length(foo)==5)
		{
		parsed.otus[i,c(4,5,6)]<-c(foo[c(2,3,5)])
		}
	if (length(foo)==6)
		{
		parsed.otus[i,c(4,5,6,7)]<-c(foo[c(2,3,6,5)])
		}
	foo<-unlist(strsplit(foo[1],";"))
	parsed.otus[i,c(1,2,3)]<-c(foo[1],substring(foo[2],6),foo[3])
	}
	
#Sample matrix
Samples<-cbind(c("A032","A138","A172","A194","A227","A229","A243","A280","A360","A361","A368","A405","A418","A420","A434","A440","A476","A482","A608","A622","A623","A636","A670","A792","A809","A832"),
c("Umb_cyl","Can_vit","Rhi_geo","Rhi_geo","Leca_schw","Psori_con","Leca_poly","Teph_atr","Leca_intr","Teph_atr","Leca_bic","Rhi_geo","Leca_poly","Aspi_myr","Leca_poly","Teph_atr","Psori_con","Leca_poly","Aspi_myr","Vari_lac","Aca_fus","Leci_lap","Leca_poly","Leci_lap","Teph_atr","Leca_bic"),
c("N","N","N","Endo_macro","N","N","N","Sky_teph","N","N","N","Mue_pyg","Lich_lec","N","Lich_lec","Mue_atr","N","Cer_epi","Sag_fis","Stig_eucl","N","Mue_pyg","Mue_pyg","N","Tae_atr","Arth_var"))

#Read Megan taxonomy
foo<-readLines("~/Desktop/ANTONIA/2_Not_collapsed/3A_Megan/2_Megan_taxonomy/All.txt")
megan.taxonomy<-mat.or.vec(length(foo),2)
for (i in 1:length(foo))
	{
		megan.taxonomy[i,]<-unlist(strsplit(foo[i],"; ;"))
	}
# Unrecorded stuff

# Parse Megan Taxonomy
colnames(final.dataset)<-c("sample","seq_length","read_repeats","megan.taxonomy","classification","identity","OTU_A","OTU_B","Kingdom","kev","Division","dev","Class","cev","Order","oev","Family","fev","Genus","gev","Species","spev")
     