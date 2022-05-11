
library(strawr)

#43 is cardiac
#hic.data.frame <- strawr::straw("KR", "RH_BR_43.hic", "11", "11", "BP", 10000)

#35 is VE
#hic.data.frame <- strawr::straw("KR", "RH_BR_35.hic", "11", "11", "BP", 10000)

#53 is Epi
#hic.data.frame <- strawr::straw("KR", "RH_BR_53.hic", "11", "11", "BP", 10000)


#mutations<-read.table('results_GroupD5-D15_CM-v-rest_logfc1_redo2.AF_intersect.bed')
#mutations<-read.table('results_GroupD5-D15_VE-v-rest_logfc1_redo2.AF_intersect.bed')
mutations<-read.table('results_GroupD5-D20_Epi-v-rest_logfc1_redo2.AF_intersect.bed')


search_distance=2e6

gene_tbl<-read.table('refGene_hg38.tsv')

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

for ( i in 1:dim(mutations)[1])
{
        chr_chr=mutations[i,1]
        chr=substring(chr_chr,4,nchar(chr_chr))
        query_start=mutations[i,2]-search_distance
        query_end=mutations[i,3]+search_distance

	mut_loc=mutations[i,13]

	query_str<-paste(chr,as.character(query_start),as.character(query_end), sep=':')
	
	#hic.data.frame<-straw("NONE","RH_BR_43.hic",query_str,query_str,"BP",10000)
	#hic.data.frame<-straw("NONE","RH_BR_35.hic",query_str,query_str,"BP",10000)
	hic.data.frame<-straw("NONE","RH_BR_53.hic",query_str,query_str,"BP",10000)

	bin_borders<-unique(c(hic.data.frame[,1],hic.data.frame[,2]))
	bin_borders[which(bin_borders<mut_loc)][-1]
	center_bin<-tail(bin_borders[which(bin_borders<mut_loc)],n=1)

	hits<-gene_tbl[which(gene_tbl[,3]==chr_chr & all_tss>query_start & all_tss<query_end),]
	hit_tss<-all_tss[which(gene_tbl[,3]==chr_chr & all_tss>query_start & all_tss<query_end)]

	read_support<-rep(0,dim(hits)[1])

	for (j in 1:dim(hits)[1])
	{
		#gene name
		hits$V13[j]
		hit_bin<-tail(bin_borders[which(bin_borders<hit_tss[j])],n=1)
		try(read_support[j]<-hic.data.frame[which(hic.data.frame[,1]==center_bin & hic.data.frame[,2]==hit_bin),3],silent=TRUE)
	}
	best_hit<-hits$V13[which(read_support==max(read_support))]
	reads_supporting_best_hit[[i]]<-max(read_support)
	hit_df=data.frame(gene_names=hits$V13,read_support=read_support,refseq=hits$V2)
	df_list[[i]]<-hit_df
	best_hits[[i]]<-best_hit[1]
}

names(df_list)<-mutations$V16

#save.image('cardiac_spec_muts.Rds')
#save.image('ve_spec_muts.Rds')
save.image('epi_spec_muts.Rds')

best_hits_df<-data.frame(mutations=mutations$V16,best_hits=best_hits, reads_supporting=reads_supporting_best_hit)

#write.table(best_hits_df,file='cardiac_best_hits.txt',quote=FALSE)
#write.table(best_hits_df,file='VE_best_hits.txt',quote=FALSE)
write.table(best_hits_df,file='epi_best_hits.txt',quote=FALSE)
