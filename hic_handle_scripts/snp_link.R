library(tidyverse)
library(strawr)

#43 is cardiac
#hic.data.frame <- strawr::straw("KR", "RH_BR_43.hic", "11", "11", "BP", 10000)

#35 is VE
#hic.data.frame <- strawr::straw("KR", "RH_BR_35.hic", "11", "11", "BP", 10000)

#53 is Epi
#hic.data.frame <- strawr::straw("KR", "RH_BR_53.hic", "11", "11", "BP", 10000)


#mutations<-read.table('results_GroupD5-D15_CM-v-rest_logfc1_redo2.AF_intersect.bed')
#mutations<-read.table('results_GroupD5-D15_VE-v-rest_logfc1_redo2.AF_intersect.bed')
#mutations<-read.table('results_GroupD5-D20_Epi-v-rest_logfc1_redo2.AF_intersect.bed')

mutations<-read.table('/data/SNP_Hi-C_check.txt',header=TRUE)

search_distance=2e6

gene_tbl<-read.table('/data/refGene_hg38.tsv')

extract_tss<-function(row){
        if (row[4]=='+')
        {
                tss<-row[5]
        }else{
                tss<-row[6]
        }
        tss<-as.numeric(gsub(" ","",tss, fixed=TRUE))
        return(tss)
}

all_tss<-apply(gene_tbl,1,extract_tss)

df_list<-vector("list",dim(mutations)[1])
best_hits<-rep(NA,dim(mutations)[1])
reads_supporting_best_hit<-rep(NA,dim(mutations)[1])
distance_to_best_hits<-rep(NA,dim(mutations)[1])

#make room for at least 400 hits per mutation, should be enough
complete_hits<-data.frame(matrix(NA,nrow = 400*dim(mutations)[1], ncol = 9))
names(complete_hits)<-c(names(mutations),"RefSeq","symbol","Count","Distance","diff_ex_score")

#index for complete table
k=1
shared_keys<-names(mutations)

for ( i in 1:dim(mutations)[1])
{
        chr_chr=mutations$Chr[i]
        chr=substring(chr_chr,4,nchar(chr_chr))
        query_start=mutations$Start[i]-search_distance
        query_end=mutations$End[i]+search_distance

	mut_type=substr(mutations$name[i],1,2)
	mut_loc=mutations$Start[i]

	query_str<-paste(chr,as.character(query_start),as.character(query_end), sep=':')
	if (mut_type=="CM") {
	   hic.data.frame<-straw("NONE","/data/RH_BR_43.hic",query_str,query_str,"BP",10000)
	} else if (mut_type=="VE") {
	  hic.data.frame<-straw("NONE","/data/RH_BR_35.hic",query_str,query_str,"BP",10000)
	} else if (mut_type=="EP"){
	  hic.data.frame<-straw("NONE","/data/RH_BR_53.hic",query_str,query_str,"BP",10000)
	} else {
	  print("UNKNOWN MUTATION TYPE")
	}

	bin_borders<-unique(c(hic.data.frame[,1],hic.data.frame[,2]))
	bin_borders[which(bin_borders<mut_loc)][-1]
	center_bin<-tail(bin_borders[which(bin_borders<mut_loc)],n=1)

	hits<-gene_tbl[which(gene_tbl[,3]==chr_chr & all_tss>query_start & all_tss<query_end),]
	hit_tss<-all_tss[which(gene_tbl[,3]==chr_chr & all_tss>query_start & all_tss<query_end)]

	read_support<-rep(0,dim(hits)[1])
	distance_to_hit<-rep(NA,dim(hits)[1])
	hit_names<-hits$V13

	for (j in 1:dim(hits)[1])
	{
		#gene name
		complete_hits[k,shared_keys]<-mutations[i,shared_keys]
		complete_hits[k,'symbol']<-hits$V13[j]
		complete_hits[k,'RefSeq']<-hits$V2[j]

		hit_bin<-tail(bin_borders[which(bin_borders<hit_tss[j])],n=1)
		try(read_support[j]<-hic.data.frame[which(hic.data.frame[,1]==center_bin & hic.data.frame[,2]==hit_bin),3],silent=TRUE)
		complete_hits[k,'Count']<-read_support[j]
		distance_to_hit[j]<-abs(hit_tss[j]-mut_loc)
		complete_hits[k,'Distance']<-abs(hit_tss[j]-mut_loc)
		
		k=k+1
	}
	
	best_hit<-hits$V13[which(read_support==max(read_support))]
	distanct_to_best_hit<-distance_to_hit[which(read_support==max(read_support))]
	reads_supporting_best_hit[[i]]<-max(read_support)
	hit_df=data.frame(gene_names=hits$V13,read_support=read_support,refseq=hits$V2)
	best_hits[[i]]<-best_hit[1]
	distance_to_best_hits[[i]]<-distanct_to_best_hit[1]
}

df_list<-mutations

#save.image('cardiac_spec_muts.Rds')
#save.image('ve_spec_muts.Rds')
#save.image('epi_spec_muts.Rds')
save.image('SNP_Hi-C_check_obj.Rdata')

df_list$top_hit=best_hits
df_list$Count=reads_supporting_best_hit
df_list$Distance=distance_to_best_hits

#best_hits_df<-data.frame(mutations=mutations$V16,best_hits=best_hits, reads_supporting=reads_supporting_best_hit)

#write.table(best_hits_df,file='cardiac_best_hits.txt',quote=FALSE)
#write.table(best_hits_df,file='VE_best_hits.txt',quote=FALSE)
#write.table(best_hits_df,file='epi_best_hits.txt',quote=FALSE)

write.table(df_list, file='/data/SNP_Hi-C_check_with_best_hits.txt', quote=FALSE, sep='\t',row.names=FALSE)

cm_rna<-read.csv('/data/one_vs_other3Lineages.DeSeq2.rna.anno/LineageCM_vs_Other3Lineages_dge_anno.csv')
ep_rna<-read.csv('/data/one_vs_other3Lineages.DeSeq2.rna.anno/LineageEpi_vs_Other3Lineages_dge_anno.csv')
ve_rna<-read.csv('/data/one_vs_other3Lineages.DeSeq2.rna.anno/LineageVE_vs_Other3Lineages_dge_anno.csv')

tab_list<-list(cm_rna,ep_rna,ve_rna)
names(tab_list)<-c('CM', 'EP', 'VE')

diffex_score<-rep(NA,dim(mutations)[1])

for ( i in 1:dim(mutations)[1])
{
	gene<-df_list[i,]$top_hit
	mut_type=substr(df_list$name[i],1,2)
	diffex_score[i]<-NA
	try(diffex_score[i]<-tab_list[[mut_type]]$log2FoldChange[which(tab_list[[mut_type]]$GeneSymbol==gene)][1])

}

df_list$diff_ex_score<-diffex_score

first_na<-min(which(is.na(complete_hits$name)))
complete_hits<-complete_hits[1:(first_na-1),]

for ( i in 1:dim(complete_hits)[1])
{
        gene<-complete_hits[i,'symbol']
        mut_type=substr(complete_hits[i,'name'],1,2)
        try(complete_hits[i,'diff_ex_score']<-tab_list[[mut_type]]$log2FoldChange[which(tab_list[[mut_type]]$GeneSymbol==gene)][1])
}


non_zero_hits<-complete_hits[which(complete_hits$Count>0),]

non_redundant_hits<-non_zero_hits[,c(1:5,7:10)]
non_redundant_hits<-unique(non_redundant_hits)

write.table(non_redundant_hits,file='/data/SNP_Hi-C_check_all_nonredundant_wDiffex.txt',quote=FALSE, sep='\t', row.names=FALSE)