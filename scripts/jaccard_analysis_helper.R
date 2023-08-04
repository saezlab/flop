#' @title Jaccard index calculation
#' @description This function calculates the Jaccard index between two lists of items for each pairwise comparison of pipelines
#' @param data A list of dataframes containing n items for each pipeline
#' @param pipelines A vector of pipeline names
#' @return A matrix of Jaccard indices
#' @export
#' @examples
#' jaccard_calc(data, pipelines)
jaccard_calc <- function(data, pipelines) {
    # creates the results matrices for storage
    jaccard_mat <- matrix(data = NA, nrow = length(pipelines), ncol = length(pipelines))
    agreement_mat <- matrix(data = NA, nrow = length(pipelines), ncol = length(pipelines))
    row.names(jaccard_mat) <- pipelines
    colnames(jaccard_mat) <- pipelines
    row.names(agreement_mat) <- pipelines
    colnames(agreement_mat) <- pipelines
    for (pipeline1 in pipelines) {
        for (pipeline2 in pipelines) {
            n_items <- data %>%
                filter(status_pipeline %in% c(pipeline1, pipeline2)) %>% 
                group_by(status_pipeline, extreme) %>%
                summarise(n_items = n())
            # chooses the minimum number of items between the two pipelines for comparison
            n_items_df <- n_items %>% group_by(extreme) %>% summarise(n_items = min(n_items))

            data_df <- data %>% filter(status_pipeline %in% c(pipeline1, pipeline2)) %>% 
                mutate(extreme = factor(extreme, levels = c('top', 'bottom'))) %>%
                group_by(status_pipeline, extreme) %>% arrange(padj, .by_group = TRUE) %>%
                group_split() %>%
                purrr::map(., function(x){
                    extreme <- x %>% pull(extreme) %>% unique()
                    n_item <- n_items_df %>% filter(extreme == !!extreme) %>% pull(n_items)
                    pipeline <- x %>% pull(status_pipeline) %>% unique()
                    subsetted_items <- x %>% slice_head(n=n_item) %>% mutate(extreme = as.character(extreme))
                    if(extreme == 'bottom'){
                        subsetted_items <- subsetted_items %>% arrange(desc(act))
                    }
                    return(subsetted_items)
                }) %>% bind_rows() 

            # if there are no items for one of the pipelines, the result value is set to NA
            if(nrow(data_df) == 0){
                jaccard_mat[pipeline1, pipeline2] <- NA
                agreement_mat[pipeline1, pipeline2] <- NA
                next
            }

            data_list <- list()
            data_list[[pipeline1]] <- data_df %>% filter(status_pipeline == pipeline1) %>% pull(items)
            data_list[[pipeline2]] <- data_df %>% filter(status_pipeline == pipeline2) %>% pull(items)

            jaccard_mat[pipeline1, pipeline2] <- bayesbio::jaccardSets(data_list[[pipeline1]], data_list[[pipeline2]])
            agreement_mat[pipeline1, pipeline2] <- length(intersect(data_list[[pipeline1]], data_list[[pipeline2]]))*2/length(c(data_list[[pipeline1]],data_list[[pipeline2]]))
        }
    }
    return(list(jaccard_mat, agreement_mat))
}

#' @title Jaccard analysis
#' @description This function calculates the Jaccard index for the top and bottom n scores for each pipeline, depending on the prior knowledge source.
#' The resulting matrix is transformed into a long-format dataframe.
#' @param merged_data A dataframe containing the scores and items for each pipeline
#' @param p_cutoff The p-value cutoff for the scores
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
            resource <- x %>%
                distinct(resource) %>%
                pull()
            # Cutoffs for the different PK sources
            if (grepl("progeny", resource)) {
                n_extr <- 3
            } else if (grepl("dorothea", resource)) {
                n_extr <- 15
            } else if (grepl("msigdb_hallmarks", resource)) {
                n_extr <- 5
            } else {
                n_extr <- round(length(unique(x$items))*0.05)
            }
            item_collector <- list()
            n_items <- tibble()
            extreme_items <- tibble()
            for (pipeline in pipelines) {
                max_items <- to_analyse %>%
                    filter(status_pipeline == !!pipeline) %>%
                    slice_max(act, n = n_extr) %>%
                    filter(padj < !!p_cutoff) %>% # If applicable, filter by p-value
                    mutate(extreme = 'top')
                min_items <- to_analyse %>%
                    filter(status_pipeline == !!pipeline) %>%
                    slice_min(act, n = n_extr) %>%
                    filter(padj < !!p_cutoff)  %>%
                    mutate(extreme = 'bottom')
                extreme_items_sub <- rbind(max_items, min_items)
                extreme_items <- rbind(extreme_items, extreme_items_sub)
                sub_n_items <- tibble(pipeline = pipeline, n_extr = nrow(extreme_items_sub))
                n_items <- bind_rows(n_items, sub_n_items)
            }

            n_items_df <- c(outer(pipelines, pipelines, FUN = paste, sep = "__")) %>%
            as_tibble() %>%
            separate(
                col=value,
                into=c("feature_1", "name"),
                sep="__"
            ) %>% 
            rowwise() %>% 
            mutate(nitems = min(n_items[n_items[, "pipeline"] == feature_1, ][["n_extr"]], n_items[n_items[, "pipeline"] == name, ][["n_extr"]]))
            
            jaccard_results <- tibble()
            agreement_results <- tibble()
            module_results <- jaccard_calc(extreme_items, pipelines)
            # The results are transformed into a long-format dataframe and merged
            jaccard_results_sub <- module_results %>%
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
                    type = 'jaccard'
                ) %>% left_join(n_items_df, by = c('feature_1', 'name'))

            agreement_results_sub <- module_results %>%
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
                    type = 'agreement'
                ) %>% left_join(n_items_df, by = c('feature_1', 'name'))

            jaccard_results <- bind_rows(jaccard_results, jaccard_results_sub)
            agreement_results <- bind_rows(agreement_results, agreement_results_sub)

            
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
            .keep_all = TRUE
        )
    return(jaccard_results)
}

