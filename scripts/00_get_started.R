ordered_cp_dotplot <- function(ck, marker_term_select_combine){ 
  ck_select=ck %>% dplyr::filter(Description %in% marker_term_select_combine)
  ck_select@compareClusterResult$Description = factor(ck_select@compareClusterResult$Description, levels=rev(marker_term_select_combine))
  
  default_labeller <- function(n) {
    function(str){
      str <- gsub("_", " ", str)
      ep_str_wrap(str, n)
    }
  }
  
  object = ck_select@compareClusterResult
  x= "Cluster"
  colorBy="p.adjust"
  showCategory=15
  
  by="count"
  size="count"
  
  #by="geneRatio"
  #size="geneRatio"
  
  
  split=NULL
  includeAll=TRUE
  font.size=12
  title=""
  label_format = 30
  group = FALSE
  shape = FALSE
  color <- NULL
  
  df <- fortify(object, showCategory=showCategory, by=by,
                includeAll=includeAll, split=split)
  
  label_func <- default_labeller(label_format)
  
  by2 <- switch(size, rowPercentage = "Percentage", 
                count         = "Count", 
                geneRatio     = "GeneRatio")    
  
  df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
  
  p <- ggplot(df, aes_string(x = 'cluster', y = "Description", size = by2))      
  
  p = p +  geom_point(aes_string(color = colorBy)) +
    scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) + 
    ylab(NULL) + ggtitle(title) + 
    DOSE::theme_dose(font.size) + 
    scale_size_continuous(range=c(3, 8)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) 
  
  return(p)
}
