library(tidyverse)

jaccard_calc <- function(data, pipelines) {
    item_mat <- matrix(data = NA, nrow = length(pipelines), ncol = length(pipelines))
    row.names(item_mat) <- pipelines
    colnames(item_mat) <- pipelines
    for (pipeline1 in pipelines) {
        for (pipeline2 in pipelines) {
            item_mat[pipeline1, pipeline2] <- bayesbio::jaccardSets(data[[pipeline1]], data[[pipeline2]])
        }
    }
    return(item_mat)
}


jaccard_analysis <- function(merged_data) {
    pipelines <- merged_data %>%
        distinct(pipeline) %>%
        pull()
    cor_results <- merged_data %>%
        group_by(statparam, resource, bio_context, status, main_dataset) %>%
        group_split() %>%
        purrr::map(., function(x) {
            to_analyse <- x %>%
                select(pipeline, scores, items) %>%
                pivot_wider(
                    names_from = pipeline,
                    values_from = scores,
                    values_fn = {
                        mean
                    }
                ) %>%
                pivot_longer(
                    cols = -items,
                    names_to = "pipeline",
                    values_to = "scores"
                )
            resource <- x %>%
                distinct(resource) %>%
                pull()
            if (resource == "progeny") {
                n_extr <- 3
            } else if (resource == "dorothea") {
                n_extr <- 15
            } else if (resource == "msigdb_hallmarks") {
                n_extr <- 5
            }
            item_collector <- list()
            for (pipeline in pipelines) {
                max_items <- to_analyse %>%
                    filter(pipeline == !!pipeline) %>%
                    slice_max(scores, n = n_extr) %>%
                    pull(items)
                min_items <- to_analyse %>%
                    filter(pipeline == !!pipeline) %>%
                    slice_min(scores, n = n_extr) %>%
                    pull(items)
                extreme_items <- c(max_items, min_items)
                item_collector[[pipeline]] <- extreme_items
            }

            cor_results <- jaccard_calc(item_collector, pipelines) %>%
                as.data.frame() %>%
                rownames_to_column(var = "feature_1") %>%
                pivot_longer(-feature_1) %>%
                mutate(
                    statparam = unique(x$statparam),
                    bio_context = unique(x$bio_context),
                    resource = unique(x$resource),
                    status = unique(x$status),
                    main_dataset = unique(x$main_dataset)
                )

            return(cor_results)
        }) %>%
        bind_rows() %>%
        subset(feature_1 != name) %>%
        rowwise() %>%
        mutate(
            id = paste0(
                sort(c(feature_1, name))[1],
                " - ",
                sort(c(feature_1, name))[2]
            )
        ) %>%
        distinct(
            id,
            statparam,
            bio_context,
            resource,
            status,
            main_dataset,
            .keep_all = TRUE
        )
    return(cor_results)
}

args <- commandArgs(trailingOnly = FALSE)
datafile <- args[grep("--file", args) + 1][2]
dataset_id <- args[grep("--dataset", args) + 1]

merged_data <- read_tsv(datafile)
jaccard_results <- jaccard_analysis(merged_data)
write_tsv(cor_results, file = paste0(dataset_id, "__jaccard.tsv"))

