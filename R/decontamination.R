#' @title Decontaminate spot swapping effect in spatial transcriptomics data
#'
#' @description This is the main function implementing the SpotClean method
#' for decontaminating spot swapping effect in spatial transcriptomics data.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{createSlide()}, or a \code{SpatialExperiment} object created from
#' \code{SpatialExperiment::read10xVisium()}.
#'
#' @param ... Arguments passed to other methods
#'
#'
#' @return For slide object created from \code{createSlide()}, returns a
#' slide object where the decontaminated expression matrix is in the
#' "decont" assay slot and the contamination statistics are in
#' metadata slots. Contamination statistics include ambient RNA contamination
#' (ARC) score, bleeding rate, distal rate, contamination radius,
#' contamination kernel weight matrix, log-likelihood value in each iteration,
#' estimated proportion of contamination in each tissue spot in observed data.
#' Since decontaminated and raw data have different number of columns, they can
#' not be stored in a single object.
#'
#' For \code{SpatialExperiment} object created from
#' \code{SpatialExperiment::read10xVisium()}, returns a
#' \code{SpatialExperiment} object where the decontaminated expression matrix
#' is in the "decont" assay slot and the contamination statistics are in
#' metadata slots. Raw expression matrix is also stored in the "counts" assay
#' slot. Genes are filtered based on \code{gene_cutoff}.
#'
#' @details Briefly, the contamination level for the slide is estimated based on
#' the total counts of all spots. UMI counts travelling around the slide are
#' assumed to follow Poisson distributions and modeled by a mixture of
#' Gaussian (proximal) and uniform (distal) kernels. The underlying
#' uncontaminated gene expressions are estimated by EM algorithm to
#' maximize the data likelihood. Detailed derivation can be found in our
#' manuscript.
#'
#'
#' @examples
#'
#' data(mbrain_raw)
#' spatial_dir <- system.file(file.path("extdata",
#'                                      "V1_Adult_Mouse_Brain_spatial"),
#'                            package = "SpotClean")
#' mbrain_slide_info <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
#'                                        "tissue_positions_list.csv"),
#'              tissue_img_file = file.path(spatial_dir,
#'                                        "tissue_lowres_image.png"),
#'              scale_factor_file = file.path(spatial_dir,
#'                                        "scalefactors_json.json"))
#' mbrain_obj <- createSlide(mbrain_raw,
#'                           mbrain_slide_info)
#'
#' mbrain_decont_obj <- spotclean(mbrain_obj, tol=10, candidate_radius=20)
#' mbrain_decont_obj


#' @rdname spotclean
#'
#' @export

spotclean <- function(slide_obj, ...) {
    if(!class(slide_obj)%in%c("SummarizedExperiment","SpatialExperiment")){
        stop("Invalid slide object.")
    }
    UseMethod(generic = "spotclean", object = slide_obj)
}

#' @param gene_keep (vector of chr) Gene names to keep for decontamination.
#' We recommend not decontaminating lowly expressed and lowly variable genes
#' in order to save computation time. Even if user include them, their
#' decontaminated expressions will not change too much from raw expressions.
#' When setting to \code{NULL}, \code{keepHighGene()} will be
#' automatically called to filter out lowly expressed and lowly variable genes
#' before decontamination.
#' Default: \code{NULL}.
#'
#' @param maxit (int) Maximum iteration for EM parameter updates. Default: 30.
#'
#' @param tol (num) Tolerance to define convergence in EM parameter updates.
#' When the element-wise maximum difference between current and updated
#' parameter matrix is less than \code{tol}, parameters are considered
#' converged. Default: 1
#'
#' @param candidate_radius (vector of num) Candidate contamination radius.
#' A series of radius to try when estimating contamination parameters.
#' Default: {c(5, 10, 15, 20, 25, 30)}
#'
#' @param kernel (chr): name of kernel to use to model local contamination.
#' Supports "gaussian", "linear", "laplace", "cauchy". Default: "gaussian".
#'
#' @param verbose (logical) Whether print progress information.
#' Default: \code{TRUE}

#' @import Matrix
#' @importFrom SummarizedExperiment assays colData SummarizedExperiment assays<-
#' @importMethodsFrom S4Vectors metadata metadata<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom dplyr filter select rename
#' @importFrom stats dist quantile optim cor lm coef
#' @importFrom methods as
#' @importFrom SpatialExperiment scaleFactors spatialCoords
#' @importFrom rlang .data
#'
#' @method spotclean SummarizedExperiment
#' @rdname spotclean
#'
#' @export

spotclean.SummarizedExperiment <- function(slide_obj, gene_keep=NULL,
                                           maxit=30, tol=1,
                                           candidate_radius=5*seq_len(6),
                                           kernel="gaussian",
                                           verbose=TRUE, ...){

    # validate arguments
    .check_arguments(gene_keep, gene_cutoff=Inf, maxit, tol,
                     candidate_radius, kernel, verbose)

    raw_data <- assays(slide_obj)$raw  # raw data matrix
    if(is.null(raw_data)){
        stop("Cannot find raw data in input slide object.")
    }
    slide <- metadata(slide_obj)$slide  # slide info

    # run SpotClean
    res <- .SpotClean(raw_data=raw_data, slide=slide,
                      gene_keep=gene_keep,
                      maxit=maxit, tol=tol,
                      candidate_radius=candidate_radius,
                      kernel=kernel,
                      verbose=verbose)

    # create output
    metadata(slide_obj)$slide <- slide[slide$tissue==1,]
    decont_obj <- createSlide(res$decont,
                              c(metadata(slide_obj),res$meta),
                              gene_cutoff = 0, verbose=FALSE)
    names(decont_obj@assays) <- "decont"

    return(decont_obj)
}

#' @param gene_cutoff (num) Filter out genes with average expressions
#' among tissue spots below or equal to this cutoff.
#' Only applies to \code{SpatialExperiment} object.
#' Default: 0.1.
#'
#' @method spotclean SpatialExperiment
#' @rdname spotclean
#'
#' @export

spotclean.SpatialExperiment <- function(slide_obj, gene_keep=NULL,
                                        gene_cutoff=0.1,
                                        maxit=30, tol=1,
                                        candidate_radius=5*seq_len(6),
                                        kernel="gaussian",
                                        verbose=TRUE, ...){
    # validate arguments
    .check_arguments(gene_keep, gene_cutoff, maxit, tol,
                     candidate_radius, kernel, verbose)

    # raw data matrix
    raw_data <- as(object = assays(slide_obj)$counts, Class = 'CsparseMatrix')
    if(is.null(raw_data)){
        stop("Cannot find raw data in input slide object.")
    }

    # collect spot info from the SpatialExperiment object
    slide <- data.frame(colData(slide_obj))
    slide <- rename(slide, tissue="in_tissue", row="array_row", col="array_col")
    slide$barcode <- rownames(slide)
    slide$tissue <- factor(as.integer(slide$tissue))
    image_pos <- spatialCoords(slide_obj)*scaleFactors(slide_obj)
    slide$imagecol <- image_pos[,"pxl_col_in_fullres"]
    slide$imagerow <- image_pos[,"pxl_row_in_fullres"]

    # filter genes
    gene_cutoff <- max(gene_cutoff,0)
    good_gene <- rowMeans(raw_data[,slide$tissue==1])>gene_cutoff
    if(verbose){
        message("Filtered out ",sum(!good_gene)
                ," genes with average expressions below or equal to ",
                gene_cutoff, ".")
    }
    raw_data <- raw_data[good_gene,]

    # run SpotClean
    res <- .SpotClean(raw_data=raw_data, slide=slide,
                      gene_keep=gene_keep,
                      maxit=maxit, tol=tol,
                      candidate_radius=candidate_radius,
                      kernel=kernel,
                      verbose=verbose)

    # create output
    decont_obj <- slide_obj[rownames(res$decont),
                            slide$barcode[slide$tissue==1]]
    assays(decont_obj)$decont <- res$decont
    metadata(decont_obj) <- c(metadata(decont_obj), res$meta)
    return(decont_obj)
}

.check_arguments <- function(gene_keep,
                             gene_cutoff,
                             maxit, tol,
                             candidate_radius,
                             kernel,
                             verbose){
    bad_args <- c()
    if(!(is.null(gene_keep) | is.character(gene_keep))){
        bad_args <- c(bad_args, "gene_keep")
    }
    if(!is.numeric(gene_cutoff)){
        bad_args <- c(bad_args, "gene_cutoff")
    }
    if(!is.numeric(maxit)) {
        bad_args <- c(bad_args, "maxit")
    }
    if(!is.numeric(candidate_radius)){
        bad_args <- c(bad_args, "candidate_radius")
    }
    if(!kernel%in%c("gaussian", "linear", "laplace", "cauchy")){
        bad_args <- c(bad_args, "kernel")
    }
    if(!is.logical(verbose)){
        bad_args <- c(bad_args, "verbose")
    }
    if(length(bad_args)>0){
        stop("invalid argument(s) (",paste(bad_args,collapse = ", "),")")
    }

    if(gene_cutoff<0) {
        stop("gene_cutoff should be non-negative")
    }
    if(maxit<=1 | maxit%%1!=0) {
        stop("maxit should be an integer greater than 1")
    }
    if(!all(candidate_radius>0)){
        stop("candidate_radius should be positive")
    }
}

.SpotClean <- function(raw_data, slide,
                       gene_keep=NULL,
                       maxit=30, tol=1,
                       candidate_radius=5*seq_len(6),
                       kernel="gaussian",
                       verbose=TRUE){

    if(verbose){
        message(Sys.time(), " Start.")
    }

    #############################
    # Step 0: setup
    #############################

    # some universal variables
    n_spots <- ncol(raw_data)  # number of spots
    ts_idx <- which(slide$tissue==1)  # tissue index
    raw_ts_data <- raw_data[,ts_idx, drop=FALSE]  # raw tissue data matrix

    # Euclidean distance matrix for spots
    slide_distance <- .calculate_euclidean_weight(
        select(slide, "imagerow", "imagecol")
    )
    rownames(slide_distance) <- colnames(slide_distance) <- slide$barcode

    # some constant matrices
    I1_y <- matrix(1,n_spots,length(ts_idx))
    I1_yy <- matrix(1,length(ts_idx),length(ts_idx))
    I_yy <- diag(length(ts_idx))

    # calculate ARC score
    ARC_score <- arcScore(raw_data,
                          background_bcs=slide$barcode[slide$tissue==0])

    # estimate spot distance in pixels
    # Note: imagecol may not correspond to col, but correspond to row
    # when the image is rotated
    if(abs(cor(slide$imagecol, slide$col))>0.99){
        lm_tmp <- lm(slide$imagecol ~ slide$col)
    }else{
        lm_tmp <- lm(slide$imagerow ~ slide$col)
    }
    spot_distance <- abs(coef(lm_tmp)[2])

    # Keep highly expressed and highly variable genes
    if(is.null(gene_keep)){
        gene_keep <- keepHighGene(raw_ts_data, verbose=verbose)

    }else{
        gene_keep <- intersect(gene_keep, rownames(raw_data))
        if(length(gene_keep)==0){
            stop("Specified genes not found in raw data.")
        }
    }

    #############################
    # Step 1: Estimate contamination parameters
    #############################

    total_counts <- colSums(raw_data[gene_keep,])
    cont_out <- c()

    if(verbose){
        message(Sys.time(), " Estimating contamination parameters...")
        pb <- txtProgressBar(max = length(candidate_radius),
                             style = 3, file = stderr())
    }

    for(i in seq_along(candidate_radius)){

        slide_weight <- .local_kernel(
            slide_distance,
            .points_to_sdv(candidate_radius[i], spot_distance),
            kernel
        )
        W_y <- t(slide_weight[ts_idx,]/rowSums(slide_weight[ts_idx,]))
        W_yy <- unname(W_y[ts_idx,])
        Wyy_tWyy <- W_yy+t(W_yy)
        WtW <- unname(crossprod(W_y))

        cont_out <- c(cont_out, .estimate_cont(obs_exp=total_counts,
                                               ts_idx=ts_idx,
                                               n_spots=n_spots,
                                               W_y=W_y, W_yy=W_yy, WtW=WtW,
                                               Wyy_tWyy=Wyy_tWyy, I_yy=I_yy,
                                               I1_y=I1_y, I1_yy=I1_yy))

        if(verbose){
            setTxtProgressBar(pb, i)
        }
    }
    gc()

    # Find the best solution
    best_radius <- which.min(unlist(lapply(cont_out,function(x) x$value)))
    bleed_rate <- cont_out[[best_radius]]$par[1]
    distal_rate <- cont_out[[best_radius]]$par[2]
    cont_radius <- candidate_radius[best_radius]

    # contamination weight matrix
    slide_weight <- .local_kernel(slide_distance,
                                  .points_to_sdv(cont_radius, spot_distance),
                                  kernel)
    slide_weight <- slide_weight[ts_idx,]/rowSums(slide_weight[ts_idx,])
    weight_mat <- distal_rate/n_spots+(1-distal_rate)*slide_weight

    #############################
    # Step 2: Decontaminate expression matrix
    #############################

    if(verbose){
        message("\n",Sys.time()," Decontaminating genes ...")
    }

    scale_factor <- rowSums(raw_data[gene_keep,])/
        rowSums(raw_ts_data[gene_keep,])
    decont_data_init <- raw_ts_data[gene_keep,]*scale_factor

    decont_out <- .decont_EM(raw_data=raw_data[gene_keep,],
                             decont_data_init = decont_data_init,
                             bleed_rate=bleed_rate,
                             distal_rate=distal_rate,
                             slide_weight = slide_weight,
                             ts_idx = ts_idx,
                             maxit=maxit,
                             tol=tol,
                             verbose=verbose)

    if(verbose){
        message("\n",Sys.time()," Scaling genes...")
    }

    decont_data <- decont_out$decont_data

    # scale up the rest genes
    gene_notkeep <- setdiff(rownames(raw_data),gene_keep)
    if(length(gene_notkeep)>0){

        scale_ratio <- rowSums(raw_data[gene_notkeep,, drop=FALSE])/
            rowSums(raw_ts_data[gene_notkeep,, drop=FALSE])
        scale_ratio[is.na(scale_ratio)] <- 1 # 0/0=NaN
        decont_data <- rbind(decont_data,
                             raw_ts_data[gene_notkeep,, drop=FALSE]*scale_ratio)
    }

    decont_data <- as(decont_data[rownames(raw_data),],"CsparseMatrix")

    # Calculation contamination rate in each tissue spot
    decont_total_counts <- colSums(decont_data)
    cont_rate <- .calculate_cont_rate(decont_total_counts,
                                      bleed_rate, distal_rate,
                                      weight_mat, n_spots)

    # write results
    meta <- list(
        bleeding_rate=bleed_rate,
        distal_rate=distal_rate,
        contamination_radius=cont_radius,
        weight_matrix=weight_mat,
        loglh=decont_out$loglh,
        decontaminated_genes=gene_keep,
        contamination_rate=cont_rate,
        ARC_score=ARC_score
    )

    if(verbose){
        message(Sys.time(), " All finished.")
    }

    return(list(decont=decont_data, meta=meta))

}

.calculate_euclidean_weight <- function(x){
    # Calculate Gaussian distances matrix of a set of 2-d points
    # Args:
    #   x (matrix of num): Each row is the coordinates of one point.
    # Returns:
    #   (matrix of num) Euclidean distance matrix

    if(ncol(x)!=2){
        stop("Incorrect dimension of input coordinates.")
    }

    # Calculate Euclidean distances
    d_mat <- dist(x, method = "euclidean", diag = TRUE)
    d_mat <- as.matrix(d_mat)
    return(d_mat)
}


.gaussian_kernel <- function(x, sigma){
    # Gaussian kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Gaussian distance

    exp(x^2/(-2*sigma^2))
}

.linear_kernel <- function(x, sigma){
    # linear kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): nonzero radius
    # Returns:
    #   (num) linear weights
    abs(sigma-pmin(x,sigma))
}

.laplace_kernel <- function(x, sigma){
    # Laplac kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Laplace weights
    exp(-abs(x)/sigma)
}

.cauchy_kernel <- function(x, sigma){
    # Cauchy kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Cauchy weights
    1/(1+x^2/sigma^2)
}

.kernel_list <- list(gaussian=.gaussian_kernel,
                     linear=.linear_kernel,
                     laplace=.laplace_kernel,
                     cauchy=.cauchy_kernel)

.local_kernel <- function(x, sigma, kernel){
    if(kernel=="linear"){
        # Linear kernel has different contamination radius
        # calculation from Gaussian,etc.
        sigma <- 2*sigma
    }
    .kernel_list[[kernel]](x,sigma)
}

.points_to_sdv <- function(cont_radius, d){
    # Transform points radius to standard deviation (the bandwidth) in
    # Gaussian kernel. This bandwidth assures that around 95% of bled-out
    # expressions are within the circle of specified spots radius
    # Args:
    #   cont_radius (int): number of points as radius where 95% contamination
    #       goes into
    #   d (num): pixel distance between two adjacent spots in a same row
    # Returns:
    #   (num): standard deviation in the Gaussian weight function

    r <- cont_radius*d # radius in pixels
    return(r/2) # sigma=r/2 -> radius=2*standard deviation -> mean +- 2SD: 95%
}

.decont_EM <- function(raw_data, decont_data_init,
                       bleed_rate, distal_rate,
                       slide_weight, ts_idx,
                       maxit=30, tol=1,
                       verbose=TRUE){
    # EM parameter updates
    # Args:
    #   raw_data (matrix of num): raw expression matrix
    #   decont_data_init (matrix of num): initial value of decontaminated
    #       expression matrix
    #   bleed_rate (num): bleeding rate
    #   distal_rate (num): distal rate
    #   slide_weight (matrix of num): Gaussian weight matrix
    #   ts_idx (vector of int): indices of tissue spots
    #   maxit (int): maximum number of iteration
    #   tol (num): tolerance for convergence
    #   verbose (logical): whether print progress messages
    # Returns:
    #   a list containing decontaminated expression matrix and a series
    #   of log-likelihood values during iteration


    # Initialization
    decont_data <- decont_data_init
    raw_ts_data <- raw_data[,ts_idx]

    Loglh <- c()
    n_it <- 1
    N_gene <- nrow(raw_data)
    N_spot <- ncol(raw_data)
    N_ts_spot <- length(ts_idx)
    bleed_weight_mat <- bleed_rate*(
        distal_rate/N_spot+(1-distal_rate)*slide_weight
    )

    # Iteration
    repeat{
        if(verbose){
            message("\nIteration: ",n_it)
        }

        # E-step

        Eta <- decont_data%*%bleed_weight_mat
        Eta[,ts_idx] <- Eta[,ts_idx]+decont_data*(1-bleed_rate)

        # stayed
        S_ <- (1-bleed_rate)*raw_ts_data*decont_data/Eta[,ts_idx]
        # contamination
        C_ <- (raw_data/Eta)%*%t(slide_weight)*
            decont_data*bleed_rate*(1-distal_rate)+
            decont_data*bleed_rate*distal_rate/N_spot *
            rowSums(raw_data/Eta)


        # M-step

        decont_data_new <- S_+C_
        loglh <- sum(raw_data*log(Eta)-Eta, na.rm = TRUE)

        # Difference between two adjacent iterations
        Lambda_maxdiff <- max(abs(decont_data_new-decont_data))

        decont_data <- decont_data_new
        Loglh <- c(Loglh, loglh)

        if(verbose){
            message("Log-likelihood: ",round(loglh,3),
                    "\nMax difference of decontaminated expressions: ",
                    round(Lambda_maxdiff,3))
        }

        if(n_it>1){
            if(Lambda_maxdiff < tol){
                if(verbose) message("Parameter converged.")
                break
            }else if(n_it>=maxit){
                if(verbose) message("Reached maximum iteration.")
                break
            }
        }

        n_it <- n_it+1

    }

    return(list(decont_data=decont_data,
                loglh=Loglh))
}


.estimate_cont <- function(obs_exp, ts_idx, n_spots,
                           W_y, W_yy, WtW, Wyy_tWyy, I_yy, I1_y, I1_yy){
    # Estimate contamination parameters using graident descent
    # Args:
    #   obs_exp (vector of num): observed expressions for all spots
    #   ts_idx (vector of int): indices of tissue spots
    #   n_spots (int): total number of spots
    #   W_y (matrix of num): contamination weight matrix
    #   W_yy (matrix of num): contamination weight matrix from tissue to tissue
    #   WtW (matrix of num): W_transpose dot-product W
    #   Wyy_tWyy (matrix of num): Wyy_transpose dot-product Wyy
    #   I_yy (matrix of num): diagonal matrix of tissue spots
    #   I1_y (matrix of num): matrix of ones
    #   I1_yy (matrix of num): matrix of ones for tissue spots
    # Returns:
    #   a list of optimization outputs from optim()

    bg_idx <- setdiff(seq_along(obs_exp),ts_idx)
    nonzero_pos <- which(obs_exp[ts_idx]>0)

    # estimate initial values
    bleed_rate_init <- sum(obs_exp[bg_idx])/length(bg_idx)*
        n_spots/sum(obs_exp)

    bg_quant <- quantile(obs_exp[bg_idx], c(0.25, 0.5))
    bg_trim <- obs_exp[bg_idx]>=bg_quant[1] & obs_exp[bg_idx]<=bg_quant[2]
    uniform_cont <- mean(obs_exp[bg_idx][bg_trim])
    distal_rate_init <- uniform_cont/sum(obs_exp[bg_idx])*length(bg_idx)
    # trim the initial value
    distal_rate_init <- min(distal_rate_init, 0.5)

    mu_init <- obs_exp[ts_idx][nonzero_pos]
    mu_init <- mu_init/sum(mu_init)*sum(obs_exp)

    # assign parameter bounds
    lower_bounds <- rep(0,length(mu_init)+2)
    lower_bounds[c(1,2)] <- c(sum(obs_exp[bg_idx])/sum(obs_exp),0.1)
    upper_bounds <- rep(Inf,length(lower_bounds))
    upper_bounds[c(1,2)] <- 1

    # calculate other matrices
    I1tZ <- crossprod(I1_y,obs_exp)
    WtZ <- crossprod(W_y,obs_exp)

    # other candidate initial values
    cand_init <- list(pmax(c(bleed_rate_init,distal_rate_init),0.1),
                      c(0.3,0.3),c(0.5,0.3))

    # minimize RSS using L-BFGS-B
    opt_list <- list()
    for(init in seq_along(cand_init)){
        x_init <- unname(c(cand_init[[init]],mu_init))
        opt_list[[init]] <- optim(x_init, .fn_optim, .gr_optim,
                                  method = "L-BFGS-B",
                                  obs_exp=obs_exp, ts_idx=ts_idx,
                                  nonzero_pos = nonzero_pos, n_spots=n_spots,
                                  W_yy=W_yy, WtW=WtW,
                                  Wyy_tWyy=Wyy_tWyy,I_yy=I_yy,
                                  I1_yy=I1_yy, WtZ=WtZ, I1tZ=I1tZ,
                                  lower=lower_bounds,
                                  upper=upper_bounds,
                                  control = list(maxit=100))
    }

    best_par <- which.min(unlist(lapply(opt_list, function(x) x$value)))
    opt <- opt_list[[best_par]]

    x_final <- numeric(length(ts_idx))
    x_final[nonzero_pos] <- opt$par[-c(1,2)]
    opt$par <- c(opt$par[c(1,2)],x_final)

    return(list(opt))
}

.calculate_cont_rate <- function(decont_total_counts,
                                 bleed_rate, distal_rate,
                                 weight_mat, n_spots){
    # Estimate proportion of contamination in each tissue spot in observed data

    # expression originated from the spot = non-bled expression +
    # contamination going to itself

    stayed_counts <- decont_total_counts*(1-bleed_rate)+
        decont_total_counts*bleed_rate*
        diag(weight_mat[names(decont_total_counts),
                        names(decont_total_counts)])

    # fitted total expression = stayed + received
    fitted_total_counts <- decont_total_counts*(1-bleed_rate)+
        ((decont_total_counts*bleed_rate)%*%
             weight_mat)[,names(decont_total_counts)]

    received_counts <- fitted_total_counts-stayed_counts

    cont_rate <- received_counts/fitted_total_counts

    return(cont_rate)
}

.fn_optim <- function(x, obs_exp, ts_idx, nonzero_pos,
                      W_yy, WtW, I_yy, Wyy_tWyy, I1_yy, n_spots, WtZ, I1tZ){
    # objective function: residual sum of squares

    x_coef <- x[-c(1,2)]
    x_c_rate <- x[1] # bleed_rate
    x_g_rate <- x[2] # distal_rate

    x_coef%*%.AtA(x_c_rate, x_g_rate,
                  WtW[nonzero_pos,nonzero_pos],
                  I_yy[nonzero_pos,nonzero_pos],
                  Wyy_tWyy[nonzero_pos,nonzero_pos],
                  I1_yy[nonzero_pos,nonzero_pos],
                  n_spots)%*%x_coef-
        2*crossprod(.AtZ(x_c_rate, x_g_rate,
                         WtZ[nonzero_pos], I1tZ[nonzero_pos],
                         obs_exp[ts_idx[nonzero_pos]],n_spots),
                    x_coef)


}


.gr_optim <- function(x, obs_exp, ts_idx, nonzero_pos,
                      W_yy, WtW, I_yy, Wyy_tWyy, I1_yy, n_spots, WtZ, I1tZ){
    # Gradient of .fn_optim

    # recover full dimension of x
    N <- length(ts_idx)
    x_coef <- numeric(N)
    x_coef[nonzero_pos] <- x[-c(1,2)]
    x_c_rate <- x[1] # bleed_rate
    x_g_rate <- x[2] # distal_rate


    x_coef1 <- x[-c(1,2)]
    x_c_rate <- x[1] # bleed_rate
    x_g_rate <- x[2] # distal_rate


    # gradient of spots expressions
    g_coef <- 2*crossprod(.AtA(x_c_rate, x_g_rate,
                               WtW[nonzero_pos,nonzero_pos],
                               I_yy[nonzero_pos,nonzero_pos],
                               Wyy_tWyy[nonzero_pos,nonzero_pos],
                               I1_yy[nonzero_pos,nonzero_pos],
                               n_spots),x_coef1)-
        2*.AtZ(x_c_rate, x_g_rate, WtZ[nonzero_pos], I1tZ[nonzero_pos],
               obs_exp[ts_idx[nonzero_pos]], n_spots)


    # gradient of bleeding rate
    g_c_rate <- x_coef1%*%.dAtA_dr(x_c_rate, x_g_rate,
                                   WtW[nonzero_pos,nonzero_pos],
                                   I_yy[nonzero_pos,nonzero_pos],
                                   Wyy_tWyy[nonzero_pos,nonzero_pos],
                                   I1_yy[nonzero_pos,nonzero_pos],
                                   n_spots) %*%x_coef1-
        2*crossprod(.dAtZ_dr(x_c_rate, x_g_rate, WtZ[nonzero_pos],
                             I1tZ[nonzero_pos], n_spots,
                             obs_exp[ts_idx[nonzero_pos]]),x_coef1)

    # gradient of distal rate
    g_g_rate <- x_coef1%*%.dAtA_dc(x_c_rate, x_g_rate,
                                   WtW[nonzero_pos,nonzero_pos],
                                   Wyy_tWyy[nonzero_pos,nonzero_pos],
                                   n_spots,
                                   I1_yy[nonzero_pos,nonzero_pos]) %*%x_coef1-
        2*crossprod(.dAtZ_dc(x_c_rate, WtZ[nonzero_pos],
                             I1tZ[nonzero_pos],
                             n_spots),x_coef1)

    c(g_c_rate,g_g_rate,g_coef)
}


# below are other internal functions for gradient calculation

.AtA <- function(bleed_rate, distal_rate, WtW,
                 I_yy, Wyy_tWyy,I1_yy, n_spots){
    (1-bleed_rate)^2*I_yy+
        bleed_rate^2*(1-distal_rate)^2*WtW+
        bleed_rate*(1-bleed_rate)*(1-distal_rate)*Wyy_tWyy+
        bleed_rate*distal_rate*(2-bleed_rate*distal_rate)/n_spots*I1_yy
}

.AtZ <- function(bleed_rate, distal_rate, WtZ, I1tZ, obs_ts_exp, n_spots){
    (1-bleed_rate)*obs_ts_exp+bleed_rate*(1-distal_rate)*WtZ+
        bleed_rate*distal_rate/n_spots*I1tZ
}

.dAtA_dr <- function(bleed_rate, distal_rate,
                     WtW, I_yy, Wyy_tWyy,I1_yy, n_spots){
    2*(bleed_rate-1)*I_yy+2*bleed_rate*(1-distal_rate)^2*WtW+
        (1-2*bleed_rate)*(1-distal_rate)*Wyy_tWyy+
        2*(distal_rate-bleed_rate*distal_rate^2)/n_spots*I1_yy
}

.dAtZ_dr <- function(bleed_rate, distal_rate,
                     WtZ, I1tZ, n_spots, obs_ts_exp){
    (1-distal_rate)*WtZ+distal_rate/n_spots*I1tZ-obs_ts_exp

}

.dAtA_dc <- function(bleed_rate, distal_rate,
                     WtW, Wyy_tWyy, n_spots, I1_yy){
    2*bleed_rate^2*(distal_rate-1)*WtW+
        bleed_rate*(bleed_rate-1)*Wyy_tWyy+
        2*(bleed_rate-bleed_rate^2*distal_rate)/n_spots*I1_yy
}

.dAtZ_dc <- function(bleed_rate, WtZ, I1tZ, n_spots){
    bleed_rate/n_spots*I1tZ-bleed_rate*WtZ

}
