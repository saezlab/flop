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
corr_analysis <- function(merged_data, corrparam) {
  cor_results <- merged_data %>%
    group_by(statparam, resource, bio_context, status, main_dataset) %>%
    group_split() %>% 
    purrr::map(., function(x) {
      to_cor <- x %>%
        select(!!corrparam, scores, items) %>%
        pivot_wider(
          names_from = !!corrparam,
          values_from = scores,
          values_fn = {mean}
        ) %>%
        column_to_rownames(var = "items") %>%
        select(order(colnames(.))) %>%
        as.matrix()

      cor_results <- cor(to_cor, method = "spearman") %>%
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
      .keep_all = TRUE)
  return(cor_results)
}
