#' @title Correlation analysis
#' @description This function performs a pairwise rank correlation analysis, for a specific statistical parameter
#' @param data A dataframe containing the data to be analysed
#' @param corrparam A string containing the name of the statistical parameter to be used for the correlation analysis
#' @return A dataframe containing the correlation results
#' @importFrom dplyr %>% select add_column as_tibble rownames_to_column column_to_rownames
#' @importFrom purrr map
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stats cor
#' @export
#' @examples
#' corr_analysis(merged_data, "stats")
corr_analysis <- function(merged_data) {
  cor_results <- merged_data %>%
    mutate(status_pipeline = paste0(status, "-", pipeline)) %>%
    group_by(statparam, resource, bio_context, main_dataset) %>%
    group_split() %>% 
    purrr::map(., function(x) {
      to_cor <- x %>%
        select(status_pipeline, scores, items) %>%
        pivot_wider(
          names_from = status_pipeline,
          values_from = scores,
          values_fn = {mean}
        ) %>%
        column_to_rownames(var = "items") %>%
        select(order(colnames(.))) %>%
        as.matrix()

      cor_results_mat <- cor(to_cor, method = "spearman", use = 'complete.obs') 
      test = cor(to_cor[complete.cases(to_cor),], method = "spearman", use = 'complete.obs') 
      all(test == cor_results_mat)
      occ_mat <- ifelse(is.na(to_cor), 0, 1)
      mat_crossprod <- crossprod(occ_mat)
      
      mat_crossprod[upper.tri(mat_crossprod)] <- NA

      mat_crossprod_df <- mat_crossprod %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "feature_1") %>% 
        pivot_longer(-feature_1) %>% 
        filter(!is.na(value)) %>%
        mutate(
          type = "occurrence",
          statparam = unique(x$statparam),
          bio_context = unique(x$bio_context),
          resource = unique(x$resource),
          main_dataset = unique(x$main_dataset)
        )

      cor_results <- cor_results_mat %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature_1") %>%
        pivot_longer(-feature_1) %>%
        mutate(
          type = "correlation",
          statparam = unique(x$statparam),
          bio_context = unique(x$bio_context),
          resource = unique(x$resource),
          main_dataset = unique(x$main_dataset)
        )
      
      full_df <- bind_rows(mat_crossprod_df, cor_results)

      return(full_df)

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
      type,
      statparam,
      bio_context,
      resource,
      main_dataset,
      .keep_all = TRUE)
  return(cor_results)
}
