#' Create cell count table from Seurat object
#'
#' @author Subrata Paul
#'
#' @param seurat.obj Seurat object
#' @param cell.type.col Column name on the seurat meta.data containing cell type annotation
#' @param cluster.col Column name on the seurat meta.data containing the cluster, default 'seurat_clusters'
#'
#'

cell_count_table<-function(seurat.obj, cell.type.col, cluster.col = 'seurat_clusters'){
  out = table(seurat.obj[[cell.type.col]][[cell.type.col]], seurat.obj[[cluster.col]][[cluster.col]])
  out = as.data.frame.matrix(out)
  names(out) = paste0('Cluster.', names(out))
  out$Total = rowSums(out)
  total_row = data.frame(t(colSums(out)))
  rownames(total_row)<-'Total'
  out = rbind(out, total_row)


  by_all = format(100*out/out['Total','Total'], digits = 2, zero.print = T, drop0trailing = T)
  by_cols = format(100*sweep(out, MARGIN = 1, STATS = unlist(out$Total), FUN = '/'),
                   digits = 2, zero.print = T, drop0trailing = T)


  less_than_one<-function(x){
    ifelse(as.numeric(x)<1 & as.numeric(x)>0, '<1', x)
  }

  formatted<-out
  for(i in 1:nrow(out)){
    for(j in 1:ncol(out)){
      added_percent = paste0(out[i,j], '(', less_than_one(by_cols[i,j]), ', ', less_than_one(by_all[i,j]),')')
      formatted[i,j]<-ifelse(as.numeric(out[i,j])>0, added_percent, '0')
    }
  }
  return(formatted)
}
