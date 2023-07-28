library(purrr)
library(ATAComb)
options(stringsAsFactors = FALSE)

dataset_names <- c("tenx_nextgem",
                   "leukopak",
                   "high_neutrophil")
dataset_well_id <- c("X000-AP0C3W1",
                     "X005-AP0C1W4",
                     "X027-AP0C1W1")
dataset_unique_in_qc_pass <- c(0.704, 0.8386, 0.3262)
dataset_filtered_frip <- c(0.69, 0.752, 0.492)
dataset_signal_frac <- dataset_unique_in_qc_pass * dataset_filtered_frip

summary_csv_files <- file.path("inst/reference/saturation",
                               paste0(dataset_well_id, "_summary.csv"))
saturation_files <- file.path("inst/reference/saturation",
                              paste0(dataset_well_id, "_saturation_curve.tsv.gz"))

summaries <- map(summary_csv_files, read.csv)
saturations <- map(saturation_files, read.table, header = F, col.names = c("count", "freq"))

saturations <- map(saturations, as.matrix)

metrics <- map(1:length(summaries),
               function(x) {
                 summary <- summaries[[x]]
                 saturation <- saturations[[x]]
                 signal_frac <- dataset_signal_frac[x]
                 data.frame(total_reads = summary$num_fragments,
                            total_umis = summary$total_usable_fragments,
                            total_counts = sum(saturation[,"count"] * saturation[,"freq"]))
               })

projections <- map(1:length(summaries),
                   function(x) {
                     saturation <- saturations[[x]]
                     metric <- metrics[[x]]
                     suppressWarnings(diversity_projection(saturation, metric, max_val = 2e9))
                   })

results <- map(1:length(projections),
               function(x) {
                 projection <- projections[[x]]
                 signal_frac <- dataset_signal_frac[x]
                 projection$signal_umis <- floor(projection$expected_umis * signal_frac)
                 projection$dataset <- dataset_names[x]
                 projection
               })

results <- do.call(rbind, results)

write.csv(results,
          "inst/reference/saturation/reference_projections.csv")

 library(ggplot2)

ggplot() +
  geom_line(data = results,
            aes(x = n_raw_reads,
                y = expected_umis,
                group = dataset,
                color = dataset),
            size = 2)

ggplot() +
  geom_line(data = results,
            aes(x = n_raw_reads,
                y = signal_umis,
                group = dataset,
                color = dataset),
            size = 2)
