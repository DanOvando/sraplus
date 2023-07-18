#' provide summary of sraplus fit
#'
#' @param fit 
#' @param output 
#'
#' @return a summary of key sraplus outputs
#' @export
#'
summarize_sralpus <- function(fit, output = "table"){
  
    result_summary <- fit$results %>% 
      dplyr::group_by(variable) %>% 
      dplyr::filter(year == max(year)) 
    
    result_summary_plot <- result_summary %>% 
      ggplot2::ggplot(aes(variable, mean, ymin = lower, ymax = upper)) + 
      ggplot2::geom_pointrange() + 
      ggplot2::facet_wrap(~variable, scales = "free") + 
      ggplot2::labs(x = '', y = "Estimate", caption = "Point is mean estimated value, bars are 95% interval. Values are in the most recent year") + 
      theme_sraplus() + 
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
      )  
    
    if (output == "table"){
      
      out <- result_summary
      
    } else {
      out <- result_summary_plot
      
    }
    
  
  return(out)
}