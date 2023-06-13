#Source: https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
#Functions partly written by Atrebas

clustering_k <- function(merged_data, k, resource, status_i) {
  cluster_results <- list()
  for (pipeline in pipelines) {
    filt_data <- merged_data %>%
      filter(
        pipeline == !!pipeline,
        statparam == !!statparam,
        resource == !!resource,
        status == !!status_i
        ) %>%
      mutate(id = sub(pattern = "__.*", "", obs_id))
    numeric_data <- filt_data %>%
      pivot_wider(
        names_from = items,
        values_from = scores,
        id_cols = id,
        values_fn = {mean}
      ) %>%
      column_to_rownames("id")
    cluster <- numeric_data %>%
      dist() %>%
      hclust(., "ave")
    hcdata <- hc_calculator(cluster, k)
    cluster_results[[status_i]][[resource]][[pipeline]] <- hcdata
  }
  return(cluster_results)
}

#' @title Hierarchical clustering calculator
#' @description Creates k number of clusters and adds information to the hierarchical clustering object
#' @param hc A hierarchical clustering object
#' @param k number of clusters
#' @importFrom ggdendro dendro_data
#' @importFrom stats cutree
#' @return A list containing the hierarchical clustering object with additional information
#' @export
#' @examples
#' hc_calculator(cluster_info, k)
hc_calculator <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}
