library(optparse)

option_list <- list(
  make_option(opt_str = c("-p","--in_pre"),
              type = "character",
              default = NULL,
              help = "Input ATAC preprocessed directory",
              metavar = "character"),
  make_option(opt_str = c("-t","--in_tenx"),
              type = "character",
              default = NULL,
              help = "Input 10x cellranger-atac outs/ directory",
              metavar = "character"),
  make_option(opt_str = c("-k", "--in_key"),
              type = "character",
              default = NULL,
              help = "Input SampleSheet.csv file",
              metavar = "character"),
  make_option(opt_str = c("-s","--in_sample"),
              type = "character",
              default = NULL,
              help = "sample",
              metavar = "character"),
  make_option(opt_str = c("-j","--in_sum"),
              type = "character",
              default = NULL,
              help = "summary files",
              metavar = "character"),
  make_option(opt_str = c("-g","--genome"),
              type = "character",
              default = "hg38",
              help = "Genome (hg38 or hg19)",
              metavar = "character"),
  make_option(opt_str = c("-q", "--qc_only"),
              type = "character",
              default = "FALSE",
              help = "Generate only QC metrics, not .arrow or .h5 files.",
              metavar = "character"),
  make_option(opt_str = c("-r","--references"),
              type = "character",
              default = "all",
              help = "Reference sets to use, comma separated (all [default] or some of altius,gene_bodies,great,grr,peaks,tss)",
              metavar = "character"),
  make_option(opt_str = c("-x","--window_sizes"),
              type = "character",
              default = "all",
              help = "Window sizes to use, comma separated (all [default] or some of window_5k,window_20k,window_100k)",
              metavar = "character"),
  make_option(opt_str = c("-d","--out_dir"),
              type = "character",
              default = NULL,
              help = "Output directory",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_html"),
              type = "character",
              default = NULL,
              help = "Output HTML run summary file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

if(is.null(args$out_html)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

rmd_loc <- "02_crossplatform_atac_qc.Rmd"

file.copy(system.file("rmarkdown/02_crossplatform_atac_qc.Rmd", package = "scATAC_Seq_Pipeline"),
          rmd_loc,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_loc,
  params = list(in_pre = args$in_pre,
                in_tenx = args$in_tenx,
                in_key = args$in_key,
                in_sum = args$in_sum,
                in_sample = args$in_sample,
                genome = args$genome,
                qc_only = args$qc_only,
                refs = args$references,
                window_sizes = args$window_sizes,
                out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_loc)