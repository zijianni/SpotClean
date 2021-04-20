#' @title Filter and return highly expressed or highly variable genes
#'
#' @description This function is for filtering genes in the expression count
#' matrix based on their average expression and variability. Usually genes
#' with low expressions and variability are less interesting and do not
#' contribute too much to downstream analyses, but rather bring technical noise.
#' It is always recommended to pre-filter gene expression matrix before
#' any analyses.
#'
#' Variability is defined using coefficient of variation (CV). We apply the
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
#' @return A filtered expression count matrix with the same
#' class as \code{count_mat}.
#'
#' @examples
#' data(MbrainSmall)
#' dim(mbrain_raw)
#' mbrain_raw_f <- KeepHighGene(mbrain_raw, mean_cutoff=100)
#' dim(mbrain_raw_f)
#'
#' @import dplyr
#' @importFrom Seurat CreateSeuratObject FindVariableFeatures
#' @import Matrix
#'
#'
#' @export


KeepHighGene <- function(count_mat, top_high=5000, mean_cutoff=1){
    # Keep highly expressed and highly variable genes.
    # When calling high variable genes, input data will be log transformed.
    # Args:
    #   count_mat (matrix-like): gene expression data matrix
    #   top_high (int): Only look for highly expressed genes within these top expressed genes.
    #   mean_cutoff (num): mean expression cutoff. Genes reach this cutoff will be kept.
    #   max_cutoff (num): max expression (among all spots) cutoff. Genes reach this cutoff will be kept.
    #   top_variable (int): Only look for highly variable genes within these top variable genes.
    # Returns:
    #   (named vector of int): indexes of resulting genes in input data.
    #        Names of each entry are the corresponding gene names.

    mean_exp <- rowMeans(count_mat)
    top_genes <- rank(-mean_exp)<=top_high # keep at most 2000 highly expressed genes.
    count_mat <- count_mat[top_genes,]
    mean_exp <- mean_exp[top_genes]

    high_exp_genes <- mean_exp>=mean_cutoff

    S_vf <- CreateSeuratObject(log1p(count_mat)) %>%
        FindVariableFeatures(selection.method = "mvp", verbose = FALSE)

    high_variable_genes <- S_vf@assays$RNA@meta.features$mvp.variable

    gene_tokeep <- high_variable_genes | high_exp_genes

    message("Kept ",sum(gene_tokeep),
            " highly expressed or highly variable genes.")

    return(count_mat[gene_tokeep,])
}
