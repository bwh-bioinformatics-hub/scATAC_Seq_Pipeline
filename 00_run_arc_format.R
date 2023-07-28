library(optparse)

option_list <- list(
  make_option(opt_str = c("-t","--in_tenx"),
              type = "character",
              default = NULL,
              help = "Input 10x cellranger-atac outs/ directory",
              metavar = "character"),
  make_option(opt_str = c("-d","--out_dir"),
              type = "character",
              default = NULL,
              help = "Output Directory",
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

rmd_loc <- "00_arc_format.rmd"

file.copy(system.file("rmarkdown/00_arc_format.rmd", package = "ATACSeqPipeline"),
          rmd_loc,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_loc,
  params = list(in_tenx = args$in_tenx,out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_loc)
