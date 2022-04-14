#following  https://bioconductor.org/packages/devel/bioc/vignettes/HiCBricks/inst/doc/IntroductionToHiCBricks.html#23_Creating_Brick_objects_from_mcool_binary_files

library('HiCBricks')

mcool_out_dir<-'~/data'

mcool_path=file.path(mcool_out_dir, "RH_BR_43.mcool")

out_dir <- file.path(tempdir(), "mcool_to_Brick_test")
dir.create(out_dir)

Create_many_Bricks_from_mcool(output_directory = out_dir,
file_prefix = "mcool_to_Brick_test", 
mcool = mcool_path, 
resolution = 10000,
experiment_name = "Testing mcool creation",
remove_existing = TRUE)

out_dir <- file.path(tempdir(), "mcool_to_Brick_test")
My_BrickContainer <- load_BrickContainer(project_dir = out_dir)

#load brick

Brick_load_data_from_mcool(Brick = My_BrickContainer,
mcool = mcool_path,
resolution = 10000,
cooler_read_limit = 10000000,
matrix_chunk = 2000,
remove_prior = TRUE,
norm_factor = NULL)

mutations<-read.table('results_GroupD5-D15_CM-v-rest_logfc1_redo2.AF_intersect.bed')

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

for ( i in 1:dim(mutations)[1])
{
	chr_chr=mutations[i,1]
	chr=substring(chr_chr,4,nchar(chr_chr))
	query_start=mutations[i,2]-search_distance
	query_end=mutations[i,3]+search_distance

	query_range=Brick_fetch_range_index(Brick = My_BrickContainer,
	chr = chr, start = query_start, end = query_end, resolution = 10000)

	hits<-gene_tbl[which(gene_tbl[,3]==chr_chr & all_tss>query_start & all_tss<query_end),]

	#x coord is the hit
	query_str<-paste(chr,as.character(query_start),as.character(query_end), sep=':')

	bins<-Brick_return_region_position(Brick = My_BrickContainer,
    	region = query_str,
	resolution = 10000)

	Sub_matrix <- Brick_get_matrix_within_coords(Brick = My_BrickContainer,
    	x_coords=query_str,
    	force = TRUE,
    	resolution = 10000,
    	y_coords = query_str)

}


#We will use promoters of Refseq genes extended to a width of 5kb
#from UCSC browser: http://genome.ucsc.edu/cgi-bin/hgTables

gene_tbl<-read.table('refGene_hg38.tsv')

refseq.db <- makeTxDbFromUCSC(genome="hg38", table="refGene")
refseq.genes = genes(refseq.db)
refseq.transcripts = transcriptsBy(refseq.db, by="gene")

refseq.transcripts = refseq.transcripts[ names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) ]
refseq_promoters <- promoters(refseq.transcripts, 2500,2500)
