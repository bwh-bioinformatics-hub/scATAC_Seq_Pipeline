library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)

## UCSC Genome Browser chrom.sizes files:
genomes <- c("hg19","hg38","mm10","mm9")

for(genome in genomes) {
  out_file <- paste0("inst/reference/", genome, ".chrom.sizes")
  download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".chrom.sizes"),
                out_file)
}

## UCSC Genome Browser chain files:
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
              "inst/reference/hg38ToHg19.over.chain.gz")
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
              "inst/reference/hg19ToHg38.over.chain.gz")

## Altius DHS Index from ENCODE:
temp_file <- tempfile(fileext = ".tsv")
download.file("https://www.encodeproject.org/files/ENCFF503GCK/@@download/ENCFF503GCK.tsv",
              temp_file)

alt_idx_raw <- fread(temp_file)

alt_idx_bed <- alt_idx_raw[,c("seqname","start","end","identifier")]
fwrite(alt_idx_bed,
       "inst/reference/hg38_altius.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

alt_idx_gr <- GRanges(seqnames = alt_idx_raw$seqname,
                      IRanges(start = alt_idx_raw$start,
                              end = alt_idx_raw$end),
                      identifier = alt_idx_raw$identifier)

saveRDS(alt_idx_gr,
        "inst/reference/hg38_altius_gr.rds")

# Conversion to hg19
temp_chain <- tempfile(fileext = ".over.chain")
R.utils::gunzip("inst/reference/hg38ToHg19.over.chain.gz",
                temp_chain,
                remove = FALSE)

ch_to_38 <- import.chain(temp_chain)

hg19_alt_idx_lo <- liftOver(alt_idx_gr, ch_to_38)
n_hg19_lo <- elementNROWS(hg19_alt_idx_lo)

keep_hg19_lo <- n_hg19_lo == 1

hg19_alt_idx <- unlist(hg19_alt_idx_lo[keep_hg19_lo])
hg19_alt_idx <- sort(hg19_alt_idx, ignore.strand = TRUE)

saveRDS(hg19_alt_idx,
        "inst/reference/hg19_altius_gr.rds")

hg19_alt_idx_bed <- as.data.frame(hg19_alt_idx)
hg19_alt_idx_bed <- hg19_alt_idx_bed[,c("seqnames","start","end","identifier")]

fwrite(hg19_alt_idx_bed,
       "inst/reference/hg19_altius.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

## Buenrostro/GSE123577 peak set
temp_chain <- tempfile(fileext = ".over.chain")
R.utils::gunzip("inst/reference/hg19ToHg38.over.chain.gz",
                temp_chain,
                remove = FALSE)

ch_to_19 <- import.chain(temp_chain)

temp_file <- tempfile(fileext = ".bed.gz")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123577&format=file&file=GSE123577%5Fpbmc%5Fpeaks%2Ebed%2Egz",
              temp_file)

hg19_peaks <- fread(temp_file)
names(hg19_peaks) <- c("chr","start","end")
hg19_gr <- convert_fragments_gr(list(hg19_peaks))[[1]]

hg38_lo <- liftOver(hg19_gr, ch_to_19)
n_lo <- elementNROWS(hg38_lo)

keep_lo <- n_lo == 1

hg38_gr <- unlist(hg38_lo[keep_lo])
hg38_gr <- sort(hg38_gr, ignore.strand = TRUE)
hg19_gr <- hg19_gr[keep_lo]
hg19_gr <- sort(hg19_gr, ignore.strand = TRUE)

saveRDS(hg19_gr, "inst/reference/hg19_peaks_gr.rds")
saveRDS(hg38_gr, "inst/reference/hg38_peaks_gr.rds")

hg19_bed <- as.data.frame(hg19_gr)[,1:3]
fwrite(hg19_bed,
       "inst/reference/hg19_peaks.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

hg38_bed <- as.data.frame(hg38_gr)[,1:3]
fwrite(hg38_bed,
       "inst/reference/hg38_peaks.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

## ENSEMBL Hg38 v93
## Couldn't download directly, so downloaded through browser from:
## ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/

hg38_chrom_sizes <- read_chrom_sizes("hg38")

gtf <- fread("C:/Users/lucasg/Downloads/Homo_sapiens.GRCh38.93.gtf.gz")
names(gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
gtf <- gtf[feature == "gene"]

gtf$attribute <- gsub(" ?[a-z|_]+ \"([A-Za-z0-9|_|-]+)\"",
                      "\\1",
                      gtf$attribute)

gtf <- separate(gtf,
                attribute,
                sep = ";",
                into = c("gene_id","gene_version","gene_name","gene_source","gene_biotype"))

gtf$gene_name <- sub('gene_name "([^"]+)"', "\\1", gtf$gene_name)
gtf <- gtf[seqname != "MT"]
gtf$seqname <- paste0("chr", gtf$seqname)
gtf <- gtf[seqname %in% hg38_chrom_sizes$chr]

tenx_feat <- fread(system.file("reference/GRCh38_10x_gene_metadata.csv.gz",
                               package = "H5weaver"))
tenx_ensembl <- tenx_feat$id

keep_gtf <- gtf[gene_id %in% tenx_ensembl]

fwrite(keep_gtf,
       "inst/reference/hg38_ensemble93_tenx_genes.tsv.gz")

# Gene Bodies
gene_gr <- GRanges(seqnames = keep_gtf$seqname,
                   ranges = IRanges(start = keep_gtf$start,
                           end = keep_gtf$end),
                   strand = keep_gtf$strand,
                   gene_id = keep_gtf$gene_id,
                   gene_name = keep_gtf$gene_name)
gene_gr <- sort(gene_gr, ignore.strand = TRUE)

saveRDS(gene_gr,
        "inst/reference/hg38_gene_bodies_gr.rds")

gene_bed <- as.data.frame(gene_gr)[,c(1:3, 6)]
fwrite(gene_bed,
       "inst/reference/hg38_gene_bodies.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# TSS Regions
tss_2kb_gr <- resize(gene_gr,
                     width = 2e3,
                     fix = "start")
tss_2kb_gr <- resize(tss_2kb_gr,
                     width = 4e3,
                     fix = "end")

tss_2kb_gr <- GenomicRanges::sort(tss_2kb_gr, ignore.strand = TRUE)
# Rename S4Vectors DFrame to DataFrame for backwards compatibility
class(tss_2kb_gr@elementMetadata) <- "DataFrame"

saveRDS(tss_2kb_gr,
        "inst/reference/hg38_tss_gr.rds")

tss_2kb_bed <- as.data.frame(tss_2kb_gr)[,c(1:3, 6)]
fwrite(tss_2kb_bed,
       "inst/reference/hg38_tss.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# Gene Regulatory Regions
# +/- 20kb from TSS, but exclude the TSS +/- 1kb
hg38_chrom_sizes <- read_chrom_sizes("hg38")

plus_genes <- gene_gr[strand(gene_gr) == "+"]
minus_genes <- gene_gr[strand(gene_gr) == "-"]

up_20kb_plus <- plus_genes
start(up_20kb_plus) <- sapply(start(plus_genes), function(x) max(x - 2e4, 1))
end(up_20kb_plus) <- sapply(start(plus_genes), function(x) max(x - 1e3, 2))

dn_20kb_plus <- plus_genes
plus_chr_sizes <- hg38_chrom_sizes$size[match(as.character(seqnames(plus_genes)), hg38_chrom_sizes$chr)]
end(dn_20kb_plus) <- mapply(min,
                            start(plus_genes) + 2e4,
                            plus_chr_sizes)
start(dn_20kb_plus) <- mapply(min,
                              start(plus_genes) + 1e3,
                              plus_chr_sizes - 1)

up_20kb_minus <- minus_genes
minus_chr_sizes <- hg38_chrom_sizes$size[match(as.character(seqnames(minus_genes)), hg38_chrom_sizes$chr)]
end(up_20kb_minus) <- mapply(min, end(minus_genes) + 2e4, minus_chr_sizes)
start(up_20kb_minus) <- mapply(min, end(minus_genes) + 1e3, minus_chr_sizes)

dn_20kb_minus <- minus_genes
start(dn_20kb_minus) <- sapply(end(minus_genes), function(x) max(x - 2e4, 1))
end(dn_20kb_minus) <- sapply(end(minus_genes), function(x) max(x - 1e3, 2))

up_regions <- c(up_20kb_plus, up_20kb_minus)
up_regions$gene_id <- paste0(up_regions$gene_id, "-up")
dn_regions <- c(dn_20kb_plus, dn_20kb_minus)
dn_regions$gene_id <- paste0(dn_regions$gene_id, "-dn")

grr_gr <- c(up_regions, dn_regions)
grr_gr <- sort(grr_gr, ignore.strand = TRUE)

saveRDS(grr_gr,
        "inst/reference/hg38_grr_gr.rds")

grr_bed <- as.data.frame(grr_gr)[,c(1:3, 6)]
fwrite(grr_bed,
       "inst/reference/hg38_grr.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# GREAT-like regions
# Core region -5kb to +1kb; Extended up to 1Mb unless intersecting another core
core_gr <- promoters(gene_gr, upstream = 5e3, downstream = 1e3)
core_gr <- sort(core_gr, ignore.strand = TRUE)

core_chrs <- unique(seqnames(core_gr))
ext_chr_grs <- split(core_gr, seqnames(core_gr))

ext_chr_grs <- lapply(core_chrs,
                      function(x) {
                        core_chr_gr <- ext_chr_grs[[x]]
                        chr_length <- hg38_chrom_sizes$size[match(x, hg38_chrom_sizes$chr)]

                        lag_core_dist <- start(core_chr_gr) - c(0, end(core_chr_gr)[1:(length(core_chr_gr) - 1)])
                        lag_core_dist[lag_core_dist > 1e6] <- 1e6
                        lag_core_dist[lag_core_dist < 0] <- 0
                        start(core_chr_gr) <- start(core_chr_gr) - lag_core_dist
                        start(core_chr_gr)[start(core_chr_gr) < 1] <- 1

                        lead_core_dist <- c(start(core_chr_gr)[2:length(core_chr_gr)], chr_length) - end(core_chr_gr)
                        lead_core_dist[lead_core_dist > 1e6] <- 1e6
                        lead_core_dist[lead_core_dist < 0] <- 0
                        end(core_chr_gr) <- end(core_chr_gr) + lead_core_dist
                        end(core_chr_gr)[end(core_chr_gr) > chr_length] <- chr_length

                        core_chr_gr
                      })

great_gr <- do.call("c", ext_chr_grs)
great_gr <- sort(great_gr, ignore.strand = TRUE)

saveRDS(great_gr,
        "inst/reference/hg38_great_gr.rds")

great_bed <- as.data.frame(great_gr)[,c(1:3, 6)]
fwrite(great_bed,
       "inst/reference/hg38_great.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

## Same process for hg19/GRCh37 build 87
hg19_chrom_sizes <- read_chrom_sizes("hg19")

gtf <- fread("C:/Users/lucasg/Downloads/Homo_sapiens.GRCh37.87.gtf.gz", skip = 5)
names(gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
gtf <- gtf[feature == "gene"]

gtf$attribute <- gsub(" ?[a-z|_]+ \"([A-Za-z0-9|_|-]+)\"",
                      "\\1",
                      gtf$attribute)

gtf <- separate(gtf,
                attribute,
                sep = ";",
                into = c("gene_id","gene_version","gene_name","gene_source","gene_biotype"))

gtf$gene_name <- sub('gene_name "([^"]+)"', "\\1", gtf$gene_name)
gtf <- gtf[seqname != "MT"]
gtf$seqname <- paste0("chr", gtf$seqname)
gtf <- gtf[seqname %in% hg19_chrom_sizes$chr]

temp_file <- tempfile(fileext = ".tar.gz")
download.file("http://cf.10xgenomics.com/samples/cell-exp/1.0.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
              temp_file, mode = "wb")

untar(temp_file, list = TRUE)
temp_dir <- tempdir()

feat_con <- untar(temp_file,files="filtered_gene_bc_matrices/hg19/genes.tsv", exdir = temp_dir)

tenx_feat <- read.table(file.path(temp_dir, "filtered_gene_bc_matrices/hg19/genes.tsv"))
names(tenx_feat) <- c("id","name")

keep_gtf <- gtf[gene_id %in% tenx_feat$id]

fwrite(keep_gtf,
       "inst/reference/hg19_ensemble87_tenx_genes.tsv.gz")

# Gene Bodies
gene_gr <- GRanges(seqnames = keep_gtf$seqname,
                   ranges = IRanges(start = keep_gtf$start,
                                    end = keep_gtf$end),
                   strand = keep_gtf$strand,
                   gene_id = keep_gtf$gene_id,
                   gene_name = keep_gtf$gene_name)
gene_gr <- sort(gene_gr, ignore.strand = TRUE)
class(gene_gr@elementMetadata) <- "DataFrame"

saveRDS(gene_gr,
        "inst/reference/hg19_gene_bodies_gr.rds")

gene_bed <- as.data.frame(gene_gr)[,c(1:3, 6)]
fwrite(gene_bed,
       "inst/reference/hg19_gene_bodies.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# TSS Regions
tss_2kb_gr <- resize(gene_gr,
                     width = 2e3,
                     fix = "start")
tss_2kb_gr <- resize(tss_2kb_gr,
                     width = 4e3,
                     fix = "end")

tss_2kb_gr <- GenomicRanges::sort(tss_2kb_gr, ignore.strand = TRUE)
# Rename S4Vectors DFrame to DataFrame for backwards compatibility
class(tss_2kb_gr@elementMetadata) <- "DataFrame"

saveRDS(tss_2kb_gr,
        "inst/reference/hg19_tss_gr.rds")

tss_2kb_bed <- as.data.frame(tss_2kb_gr)[,c(1:3, 6)]
fwrite(tss_2kb_bed,
       "inst/reference/hg19_tss.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# Gene Regulatory Regions
# +/- 20kb from TSS, but exclude the TSS +/- 1kb
hg19_chrom_sizes <- read_chrom_sizes("hg19")

plus_genes <- gene_gr[strand(gene_gr) == "+"]
minus_genes <- gene_gr[strand(gene_gr) == "-"]

up_20kb_plus <- plus_genes
start(up_20kb_plus) <- sapply(start(plus_genes), function(x) max(x - 2e4, 1))
end(up_20kb_plus) <- sapply(start(plus_genes), function(x) max(x - 1e3, 2))

dn_20kb_plus <- plus_genes
plus_chr_sizes <- hg19_chrom_sizes$size[match(as.character(seqnames(plus_genes)), hg19_chrom_sizes$chr)]
end(dn_20kb_plus) <- mapply(min,
                            start(plus_genes) + 2e4,
                            plus_chr_sizes)
start(dn_20kb_plus) <- mapply(min,
                              start(plus_genes) + 1e3,
                              plus_chr_sizes - 1)

up_20kb_minus <- minus_genes
minus_chr_sizes <- hg19_chrom_sizes$size[match(as.character(seqnames(minus_genes)), hg19_chrom_sizes$chr)]
end(up_20kb_minus) <- mapply(min, end(minus_genes) + 2e4, minus_chr_sizes)
start(up_20kb_minus) <- mapply(min, end(minus_genes) + 1e3, minus_chr_sizes)

dn_20kb_minus <- minus_genes
start(dn_20kb_minus) <- sapply(end(minus_genes), function(x) max(x - 2e4, 1))
end(dn_20kb_minus) <- sapply(end(minus_genes), function(x) max(x - 1e3, 2))

up_regions <- c(up_20kb_plus, up_20kb_minus)
up_regions$gene_id <- paste0(up_regions$gene_id, "-up")
dn_regions <- c(dn_20kb_plus, dn_20kb_minus)
dn_regions$gene_id <- paste0(dn_regions$gene_id, "-dn")

grr_gr <- c(up_regions, dn_regions)
grr_gr <- sort(grr_gr, ignore.strand = TRUE)
class(grr_gr@elementMetadata) <- "DataFrame"

saveRDS(grr_gr,
        "inst/reference/hg19_grr_gr.rds")

grr_bed <- as.data.frame(grr_gr)[,c(1:3, 6)]
fwrite(grr_bed,
       "inst/reference/hg19_grr.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

# GREAT-like regions
# Core region -5kb to +1kb; Extended up to 1Mb unless intersecting another core
core_gr <- promoters(gene_gr, upstream = 5e3, downstream = 1e3)
core_gr <- sort(core_gr, ignore.strand = TRUE)

core_chrs <- unique(seqnames(core_gr))
ext_chr_grs <- split(core_gr, seqnames(core_gr))

ext_chr_grs <- lapply(core_chrs,
                      function(x) {
                        core_chr_gr <- ext_chr_grs[[x]]
                        chr_length <- hg19_chrom_sizes$size[match(x, hg19_chrom_sizes$chr)]

                        lag_core_dist <- start(core_chr_gr) - c(0, end(core_chr_gr)[1:(length(core_chr_gr) - 1)])
                        lag_core_dist[lag_core_dist > 1e6] <- 1e6
                        lag_core_dist[lag_core_dist < 0] <- 0
                        start(core_chr_gr) <- start(core_chr_gr) - lag_core_dist
                        start(core_chr_gr)[start(core_chr_gr) < 1] <- 1

                        lead_core_dist <- c(start(core_chr_gr)[2:length(core_chr_gr)], chr_length) - end(core_chr_gr)
                        lead_core_dist[lead_core_dist > 1e6] <- 1e6
                        lead_core_dist[lead_core_dist < 0] <- 0
                        end(core_chr_gr) <- end(core_chr_gr) + lead_core_dist
                        end(core_chr_gr)[end(core_chr_gr) > chr_length] <- chr_length

                        core_chr_gr
                      })

great_gr <- do.call("c", ext_chr_grs)
great_gr <- sort(great_gr, ignore.strand = TRUE)

saveRDS(great_gr,
        "inst/reference/hg19_great_gr.rds")

great_bed <- as.data.frame(great_gr)[,c(1:3, 6)]
fwrite(great_bed,
       "inst/reference/hg19_great.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)


### Seurat RNA-seq reference
library(Seurat)
library(H5weaver)

tenx_genes <- fread("hg38_ensemble93_tenx_genes.tsv.gz")
tenx_genes$symbols <- make.unique(tenx_genes$gene_name)

if(!file.exists("pbmc_10k_v3.rds")) {
  download.file("https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1",
                "pbmc_10k_v3.rds")
}

pbmc_so <- readRDS("pbmc_10k_v3.rds")
pbmc_mat <- pbmc_so@assays$RNA@counts

pbmc_genes <- rownames(pbmc_mat)
keep_genes <- pbmc_genes %in% tenx_genes$symbols
pbmc_mat <- pbmc_mat[keep_genes,]

pbmc_genes <- rownames(pbmc_mat)
pbmc_ensembl <- tenx_genes$gene_id[match(pbmc_genes, tenx_genes$symbols)]
#rownames(pbmc_mat) <- pbmc_ensembl

pbmc_meta <- pbmc_so@meta.data
pbmc_meta$original_barcodes <- rownames(pbmc_meta)

pbmc_meta <- pbmc_meta[pbmc_meta$celltype != "Platelets",]
pbmc_mat <- pbmc_mat[,pbmc_meta$original_barcodes]

pbmc_meta <- lapply(pbmc_meta,
                    function(x) {
                      if(class(x) == "factor") {
                        as.character(x)
                      } else {
                        x
                      }
                    })

name_conversion <- c(B.Naive = "pre-B cell",
                     B.Activated = "B cell progenitor",
                     T.CD4.Naive = "CD4 Naive",
                     T.CD8.Naive = "CD8 Naive",
                     T.CD4.Memory = "CD4 Memory",
                     T.CD8.Effector = "CD8 effector",
                     T.DoubleNegative = "Double negative T cell",
                     NK = "NK cell",
                     DC.Plasmacytoid = "pDC",
                     DC.Myeloid = "Dendritic cell",
                     Mono.CD14 = "CD14+ Monocytes",
                     Mono.CD16 = "CD16+ Monocytes")

pbmc_meta$celltype <- names(name_conversion[match(pbmc_meta$celltype, name_conversion)])

pbmc_list <- list(matrix_dgCMatrix = pbmc_mat,
                  matrix = list(observations = pbmc_meta,
                                features = list(name = pbmc_genes))
                  )

pbmc_list <- h5_list_convert_from_dgCMatrix(pbmc_list, "matrix")

write_h5_list(pbmc_list, "seurat_rna_pbmc_ref.h5")

file.remove("pbmc_10k_v3.rds")
