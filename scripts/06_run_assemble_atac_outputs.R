library(optparse)

option_list <- list(
  make_option(opt_str = c("-p","--in_post"),
              type = "character",
              default = NULL,
              help = "Input postprocessed ATAC directory",
              metavar = "character"),
  make_option(opt_str = c("-f","--in_frag"),
              type = "character",
              default = NULL,
              help = "Input ATAC filtered fragments.tsv.gz",
              metavar = "character"),
  make_option(opt_str = c("-m","--in_meta"),
              type = "character",
              default = NULL,
              help = "Input ATAC filtered metadata.csv.gz",
              metavar = "character"),
  make_option(opt_str = c("-g","--genome"),
              type = "character",
              default = "hg38",
              help = "Genome (hg38 or hg19)",
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
  make_option(opt_str = c('-n', "--n_cores"),
              type = "character",
              default = "auto",
              help = "Number of cores to use for parallel processing",
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

if(is.null(args$out_dir)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

if(!dir.exists(args$out_dir)) {
  dir.create(args$out_dir)
}

rmd_loc <- file.path(args$out_dir, "assemble_atac_outputs.Rmd")

file.copy(system.file("rmarkdown/assemble_atac_outputs.Rmd", package = "scATACSeqPipeline"),
          rmd_loc,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_loc,
  params = list(in_post = args$in_post,
                in_frag = args$in_frag,
                in_meta = args$in_meta,
                genome = args$genome,
                refs = args$references,
                window_sizes = args$window_sizes,
                n_cores = args$n_cores,
                out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_loc)
