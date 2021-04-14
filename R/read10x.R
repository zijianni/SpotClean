#' @name Read10xRaw
#' @rdname Read10xRaw
#' 
#' @title Read 10x output data
#'
#' @description \code{Read10xRaw} is a one-line handy function for reading 
#' 10x Cell Ranger output data, producing a count matrix for input to 
#' \code{CB2FindCell}. \code{Read10xRawH5} is for reading 10x Cell Ranger output
#' HDF5 file (ended with .h5).  Works under both old (<3) and new (>=3) 
#' Cell Ranger version.
#'
#' @param dir The directory of 10x output data. For Cell Ranger version <3,
#' directory should include three files: barcodes.tsv, genes.tsv, matrix.mtx.
#' For Cell Ranger version >=3, directory should include three
#' files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#' 
#' @param h5file The path of 10x output HDF5 file (ended with .h5).
#' 
#' @param row.name Specify either using gene symbols 
#' (\code{row.name = "symbol"}) or gene Ensembl IDs (\code{row.name = "id"})
#' as row names of the count matrix. Default is \code{row.name = "symbol"}.
#' 
#' @param meta Logical. If \code{TRUE}, returns a list containing both the 
#' count matrix and metadata of genes (features). Metadata includes feature 
#' names, IDs and other additional information depending on Cell Ranger 
#' output. If \code{FALSE} (default), only returns the count matrix.
#'
#' @return If \code{meta = TRUE}, returns a list of two elements: a 
#' "dgCMatrix" sparse matrix containing expression counts and a data 
#' frame containing metadata of genes (features). For the count matrix, 
#' each row is a gene (feature) and each column is a barcode.  If 
#' \code{meta = FALSE}, only returns the count matrix.
#'
#' @examples
#' 
#' # simulate 10x output files
#' data(mbrainSub)
#' data_dir <- file.path(tempdir(),"CB2example")
#' dir.create(data_dir)
#' gene_name <- rownames(mbrainSub)
#' 
#' # For simplicity, use gene names to generate gene IDs to fit the format.
#' gene_id <- paste0("ENSG_fake_",gene_name)
#' barcode_id <- colnames(mbrainSub)
#' Matrix::writeMM(mbrainSub,file = file.path(data_dir,"matrix.mtx"))
#' write.table(barcode_id,file = file.path(data_dir,"barcodes.tsv"),
#'     sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#' write.table(cbind(gene_id,gene_name),file = file.path(data_dir,"genes.tsv"),
#'     sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#' 
#' # read files
#' list.files(data_dir)
#' mbrainSub_new <- Read10xRaw(data_dir)
#' str(mbrainSub_new)
#' identical(mbrainSub, mbrainSub_new)
#'
#' @importFrom utils read.delim
#' @importFrom Matrix readMM
#' @importFrom Matrix sparseMatrix
#' 
#'
#' @export

Read10xRaw <- function(dir = NULL,
                    row.name = "symbol",
                    meta = FALSE) {
    if(is.null(dir)){
        dir <- getwd()
    }
    if (!row.name %in% c("symbol", "id"))
        stop("row.name should be either \"symbol\" or \"id\".")
    if (!dir.exists(dir))
        stop("Directory does not exist.")
    dir <- gsub("/$", "", dir)
    fname <- list.files(dir)
    
    #check if Cell Ranger version >=3
    V3 <- "features.tsv.gz" %in% fname
    if (V3) {
        if(!all(c("barcodes.tsv.gz","features.tsv.gz",
                "matrix.mtx.gz")%in%fname)){
            stop("No 10x output file detected.")
        }
        Barcode <- file.path(dir, "barcodes.tsv.gz")
        Gene <- file.path(dir, "features.tsv.gz")
        CountMat <- file.path(dir, "matrix.mtx.gz")
    } else{
        if(!all(c("barcodes.tsv","genes.tsv",
                "matrix.mtx")%in%fname)){
            stop("No 10x output file detected.")
        }
        Barcode <- file.path(dir, "barcodes.tsv")
        Gene <- file.path(dir, "genes.tsv")
        CountMat <- file.path(dir, "matrix.mtx")
    }
    barcode <- readLines(Barcode)
    gene.meta <-
        read.delim(Gene, header = FALSE, colClasses = "character")
    if (V3) {
        colnames(gene.meta) <- c("id", "symbol", "type")
    } else{
        colnames(gene.meta) <- c("id", "symbol")
    }
    if (row.name == "symbol") {
        #gene names as row names of count matrix
        gene <- gene.meta$symbol
    } else{
        #gene ids as row names of count matrix
        gene <- gene.meta$id
    }
    countmat <- as(readMM(CountMat), "dgCMatrix")
    colnames(countmat) <- barcode
    rownames(countmat) <- make.unique(gene)
    if (meta) {
        #return a list including count matrix and gene metadata
        return(list(CountMatrix = countmat, Metadata = gene.meta))
    } else{
        #only return count matrix
        return(countmat)
    }
}



#' @rdname Read10xRaw
#' @importFrom rhdf5 h5ls
#' @importFrom rhdf5 h5read
#'
#' @export

Read10xRawH5 <- function(h5file,
                    row.name = "symbol",
                    meta = FALSE) {
    if (!row.name %in% c("symbol", "id"))
        stop("row.name should be either \"symbol\" or \"id\".")
    if (!file.exists(h5file))
        stop("File does not exist.")
    fname <- h5ls(h5file)
    #check if Cell Ranger version >=3
    V3 <- "/matrix" %in% fname$group
    if (V3) {
        data.temp <-  h5read(h5file, "/matrix")
        barcode <- data.temp$barcodes
        gene.meta <- data.temp$features
        if (row.name == "symbol") {
            #gene names as row names of count matrix
            gene <- gene.meta$name
        } else{
            #gene ids as row names of count matrix
            gene <- gene.meta$id
        }
        countmat <-
            sparseMatrix(
                i = data.temp$indices,
                p = data.temp$indptr,
                x = as.numeric(data.temp$data),
                index1 = FALSE,
                dims = data.temp$shape
            )
        colnames(countmat) <- as.vector(barcode)
        rownames(countmat) <- make.unique(gene)
        
        if (meta) {
            #return a list including count matrix and gene metadata
            return(list(CountMatrix = countmat, Metadata = gene.meta))
        } else{
            #only return count matrix
            return(countmat)
        }
    } else{
        subname <- fname$name[fname$group == "/"]
        #Cell Ranger V2 data might include multiple count matrices
        CountMatList <- vector("list", length(subname))
        names(CountMatList) <- subname
        for (i in seq_along(subname)) {
            data.temp <- h5read(h5file, paste0("/", subname[i]))
            barcode <- data.temp$barcodes
            gene.meta <-
                data.frame(symbol = data.temp$gene_names, id = data.temp$genes)
            if (row.name == "symbol") {
                #gene names as row names of count matrix
                gene <- gene.meta$symbol
            } else{
                #gene ids as row names of count matrix                          
                gene <- gene.meta$id
            }
            countmat <-
                sparseMatrix(
                    i = data.temp$indices,
                    p = data.temp$indptr,
                    x = data.temp$data,
                    dims = data.temp$shape,
                    index1 = FALSE,
                    giveCsparse = TRUE
                )
            colnames(countmat) <- as.vector(barcode)
            rownames(countmat) <- make.unique(gene)
            if (meta) {
                #return a list including count matrix and gene metadata
                CountMatList[[i]] <-
                    list(CountMatrix = countmat, Metadata = gene.meta)
            } else{
                #only return count matrix
                CountMatList[[i]] <- countmat
            }
        }
        if (length(unique(fname$group)) > 2) {
            warning("Detected multiple datasets. 
                    Returning a list with multiple matrices.")
            return(CountMatList)
        } else{
            return(CountMatList[[1]])
        }
    }
}
