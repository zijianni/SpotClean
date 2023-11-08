#' @title Filter and return highly expressed or highly variable genes
#'
#' @description This function is for filtering genes in the expression count
#' matrix based on their average expression and variability. Usually genes
#' with low expressions and variability are less interesting and do not
#' contribute too much to downstream analyses, but rather bring technical noise.
#' It is always recommended to pre-filter gene expression matrix before
#' any analyses.
#'
#' We apply the
#' function \code{FindVariableFeatures} with \code{selection.method = "mvp"}
#' from package \code{Seurat} on log transformed expression matrix to
#' detect high variable genes. This method does not require a pre-specified
#' number of high variable genes.
#'
#' @param count_mat (matrix of num) Input count matrix to be filtered.
#' Can be either standard matrix format or sparse matrix format.
#'
#' @param top_high (int) Only look for highly expressed and variable genes
#' within this number of top expressed genes. Default: 5000.
#'
#' @param mean_cutoff (num) Genes with average expressions among all
#' spots exceeding this cutoff are kept as highly expressed genes.
#'
#' @param return_matrix (logical) Whether return filtered matrix instead
#' of gene names. Default: \code{FALSE}
#'
#' @param verbose (logical) Whether print progress information.
#' Default: \code{TRUE}
#'
#' @return A vector of gene names or a filtered expression count
#' matrix with the same class as \code{count_mat}.
#'
#' @examples
#'
#' data(mbrain_raw)
#' dim(mbrain_raw)
#'
#' mbrain_raw_f <- keepHighGene(mbrain_raw, mean_cutoff=100,
#'                              return_matrix=TRUE)
#' dim(mbrain_raw_f)
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Matrix rowMeans
#'
#'
#' @export


keepHighGene <- function(count_mat, top_high=5000,
                         mean_cutoff=1, return_matrix=FALSE,
                         verbose=TRUE){

    mean_exp <- rowMeans(count_mat)
    # keep at most this number of highly expressed genes.
    top_genes <- rank(-mean_exp)<=top_high
    count_mat <- count_mat[top_genes,,drop=FALSE]
    mean_exp <- mean_exp[top_genes]

    high_exp_genes <- mean_exp>=mean_cutoff

    S_vf <- NormalizeData(CreateSeuratObject(count_mat), verbose = FALSE)
    S_vf <- FindVariableFeatures(S_vf,
                                 selection.method = "mvp", verbose = FALSE)
    
    # accommodate changes in Seurat v5 objects
    if(as.integer(gsub("\\<(\\d+)\\.\\d+\\.\\d+", "\\1", S_vf@version))>=5){
        high_variable_genes <- S_vf@assays$RNA@meta.data$vf_mvp_data_variable
    }else{
        high_variable_genes <- S_vf@assays$RNA@meta.features$mvp.variable
    }
    
    gene_tokeep <- high_variable_genes | high_exp_genes

    if(verbose){
        message("Kept ",sum(gene_tokeep),
                " highly expressed or highly variable genes.")
    }

    if(return_matrix){
        return(count_mat[gene_tokeep,, drop=FALSE])
    }else{
        return(names(which(gene_tokeep)))
    }

}
