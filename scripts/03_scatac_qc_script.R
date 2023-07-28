library(optparse)

option_list <- list(
  make_option(opt_str = c("-s","--in_sample"),
              type = "character",
              default = NULL,
              help = "Batch identifier",
              metavar = "character"),
  make_option(opt_str = c("-m","--in_method"),
              type = "character",
              default = NULL,
              help = "Input batch pipeline modality string",
              metavar = "character"),
  make_option(opt_str = c("-i","--in_dir"),
              type = "character",
              default = NULL,
              help = "Input directory containing h5 and json files",
              metavar = "character"),
  make_option(opt_str = c("-k","--in_key"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
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

if(is.null(args$in_sample)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

if(!dir.exists(args$out_dir)) {
  dir.create(args$out_dir)
}

rmd_path <- file.path(args$out_dir,
                      paste0(args$in_sample,
                             "_scatac_seq_qc_report_generator.Rmd"))

file.copy(system.file("rmarkdown/scatac_seq_qc_report_generator.Rmd", package = "scATACSeqPipeline"),
          rmd_path,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_path,
  params = list(in_sample = args$in_sample,
                in_method = args$in_method,
                in_dir  = args$in_dir,
                in_key  = args$in_key,
                out_dir = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_path)
