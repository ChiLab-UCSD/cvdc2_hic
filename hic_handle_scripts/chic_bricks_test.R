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

mcool_out_dir <- file.path(tempdir(), "mcool_out_dir")
mcool_path=file.path(mcool_out_dir, "RH_BR_43.mcool")

Brick_load_data_from_mcool(Brick = My_BrickContainer,
mcool = mcool_path,
resolution = 10000,
cooler_read_limit = 10000000,
matrix_chunk = 2000,
remove_prior = TRUE,
norm_factor = "Knight-Ruitz")

mutations<-read.table('results_GroupD5-D15_CM-v-rest_logfc1_redo2.AF_intersect.bed')

search_distance=5e6

for ( i in 1:dim(mutations)[1])
{
	chr=mutations[i,1]
	start=mutations[i,2]-search_distance
	end=mutations[i,3]+search_distance

	Brick_fetch_range_index(Brick = My_BrickContainer,
	chr = chr, start = start, end = end, resolution = 10000)

}
