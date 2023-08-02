
#' @title Jaccard index calculation
#' @description This function calculates the Jaccard index between two lists of items for each pairwise comparison of pipelines
#' @param data A list of dataframes containing n items for each pipeline
#' @param pipelines A vector of pipeline names
#' @return A matrix of Jaccard indices
#' @export
#' @examples
#' jaccard_calc(data, pipelines)
jaccard_calc <- function(data, pipelines) {
    jaccard_mat <- matrix(data = NA, nrow = length(pipelines), ncol = length(pipelines))
    agreement_mat <- matrix(data = NA, nrow = length(pipelines), ncol = length(pipelines))
    row.names(jaccard_mat) <- pipelines
    colnames(jaccard_mat) <- pipelines
    row.names(agreement_mat) <- pipelines
    colnames(agreement_mat) <- pipelines
    for (pipeline1 in pipelines) {
        for (pipeline2 in pipelines) {
            min_n_items <- min(length(data[[pipeline1]]), length(data[[pipeline2]]))
            jaccard_mat[pipeline1, pipeline2] <- bayesbio::jaccardSets(data[[pipeline1]][1:min_n_items], data[[pipeline2]][1:min_n_items])
            agreement_mat[pipeline1, pipeline2] <- length(intersect(data[[pipeline1]], data[[pipeline2]]))*2/length(c(data[[pipeline1]][1:min_n_items],data[[pipeline2]][1:min_n_items]))
        }
    }
    return(list(jaccard_mat, agreement_mat))
}

#' @title Jaccard analysis
#' @description This function calculates the Jaccard index for the top and bottom n scores for each pipeline, depending on the prior knowledge source.
#' The resulting matrix is transformed into a long-format dataframe.
#' @param merged_data A dataframe containing the scores and items for each pipeline
#' @return A long-format dataframe containing the Jaccard indices for every combination of pipelines
#' @export
#' @examples
#' jaccard_analysis(merged_data)
jaccard_analysis <- function(merged_data, p_cutoff = 1) {
    merged_data <- merged_data %>%
        mutate(status_pipeline = paste0(status, "-", pipeline))
    jaccard_results <- merged_data %>%
        group_by(statparam, resource, bio_context, main_dataset, subset) %>%
        group_split() %>%
        purrr::map(., function(x) {
            pipelines <- x %>%
                distinct(status_pipeline) %>%
                arrange() %>%
                pull()
            to_analyse <- x %>%
                select(status_pipeline, act, items, padj)
                # pivot_wider(
                #     names_from = status_pipeline,
                #     values_from = act,
                #     values_fn = {
                #         mean
                #     }
                # ) %>%
                # pivot_longer(
                #     cols = -items,
                #     names_to = "status_pipeline",
                #     values_to = "scores"
                # )
            resource <- x %>%
                distinct(resource) %>%
                pull()
            if (grepl("progeny", resource)) {
                n_extr <- 3
            } else if (grepl("dorothea", resource)) {
                n_extr <- 15
            } else if (grepl("msigdb_hallmarks", resource)) {
                n_extr <- 5
            } else {
                n_extr <- round(length(x$items)*0.05)
            }
            item_collector <- list()
            n_items <- tibble()
            for (pipeline in pipelines) {
                max_items <- to_analyse %>%
                    filter(status_pipeline == !!pipeline) %>%
                    slice_max(act, n = n_extr) %>%
                    filter(padj < !!p_cutoff) %>%
                    pull(items)
                min_items <- to_analyse %>%
                    filter(status_pipeline == !!pipeline) %>%
                    slice_min(act, n = n_extr) %>%
                    filter(padj < !!p_cutoff) %>%
                    pull(items)
                extreme_items <- c(max_items, min_items)
                # item_collector[['max']][[pipeline]] <- max_items
                # item_collector[['min']][[pipeline]] <- min_items
                item_collector[[pipeline]] <- extreme_items
                sub_n_items <- tibble(pipeline = pipeline, n_extr = length(extreme_items))
                n_items <- bind_rows(n_items, sub_n_items)
            }

            n_items_df <- c(outer(pipelines,pipelines,FUN=paste,sep="__")) %>%
            as_tibble() %>%
            separate(
                col=value,
                into=c("feature_1", "name"),
                sep="__"
            ) %>% 
            rowwise() %>% 
            mutate(nitems = min(n_items[n_items[, "pipeline"]==feature_1, ][["n_extr"]], n_items[n_items[, "pipeline"]==name, ][["n_extr"]]))
            
            jaccard_results <- tibble()
            agreement_results <- tibble()
            # extremes <- c('max', 'min')
            # for(extreme in extremes){
                jaccard_results_sub <- jaccard_calc(item_collector, pipelines) %>% #item_collector[[extreme]]
                    .[[1]] %>%
                    as.data.frame() %>%
                    rownames_to_column(var = "feature_1") %>%
                    pivot_longer(-feature_1) %>%
                    mutate(
                        statparam = unique(x$statparam),
                        bio_context = unique(x$bio_context),
                        resource = unique(x$resource),
                        main_dataset = unique(x$main_dataset),
                        subset = unique(x$subset),
                        # extreme = !!extreme,
                        type = 'jaccard'
                    ) %>% left_join(n_items_df, by = c('feature_1', 'name'))

                agreement_results_sub <- jaccard_calc(item_collector, pipelines) %>%
                    .[[2]] %>%
                    as.data.frame() %>%
                    rownames_to_column(var = "feature_1") %>%
                    pivot_longer(-feature_1) %>%
                    mutate(
                        statparam = unique(x$statparam),
                        bio_context = unique(x$bio_context),
                        resource = unique(x$resource),
                        main_dataset = unique(x$main_dataset),
                        subset = unique(x$subset),
                        # extreme = !!extreme,
                        type = 'agreement'
                    ) %>% left_join(n_items_df, by = c('feature_1', 'name'))

                jaccard_results <- bind_rows(jaccard_results, jaccard_results_sub)
                agreement_results <- bind_rows(agreement_results, agreement_results_sub)
            # }

            
            results_df <- bind_rows(jaccard_results, agreement_results)

            return(results_df)
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
            main_dataset,
            subset,
            type,
            # extreme,
            .keep_all = TRUE
        )
    return(jaccard_results)
}

