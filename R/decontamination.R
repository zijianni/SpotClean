#' @title Decontaminate spot swapping effect in spatial transcriptomics data
#'
#' @description This is the main function implementing the TBD method
#' for decontaminating spot swapping effect in spatial transcriptomics data.
#'
#' @param slide_obj A slide object created or inherited from
#' \code{CreateSlide()}.
#'
#' @param gene_keep (vectof of chr) Gene names to keep for decontamination.
#' We recommend not decontaminating lowly expressed and lowly variable genes
#' in order to save computation time. Even if user include them, their
#' decontaminated expressions will not change too much from raw expressions.
#' When setting to \code{NULL}, \code{KeepHighGene()} will be
#' automatically called to filter out lowly expressed and lowly variable genes
#' before decontamination.
#' Default: \code{NULL}.
#'
#' @param maxit (int) Maximum iteration for EM parameter updates. Default: 30.
#'
#' @param tol (num) Tolerance to define convergence in EM parameter updates.
#' When the element-wise maximum difference between current and updated
#' parameter matrix is less than \code{tol}, parameters are considered
#' converged. Default: 1.
#'
#' @param candidate_radius (vector of num) Candidate contamination radius.
#' A series of radius to try when estimating contamination parameters.
#' Default: {5, 10, 15, 20, 25, 30}
#'
#' @param verbose (logical) Whether print progress information.
#' Default: \code{TRUE}
#'
#' @return A slide object where the decontaminated expression matrix is in the
#' "decont" assay slot and the contamination statistics are in
#' metadata slots. Contamination statistics include ambient RNA contamination
#' (ARC) score, bleeeding rate, global contamination rate, contamination radius,
#' contamination kernel weight matrix, log-likelihood value in each iteration.
#' Since decontaminated and raw data have different columns, they can
#' not be stored in a single SummarizedExperiment object.
#' Do not overwrite the raw slide object.
#'
#' @details Briefly, the contamination level for the slide is estimated based on
#' the total counts of all spots. UMI counts travelling around the slide are
#' assumed to follow Poisson distributions and modeled by a mixture of
#' Gaussian and uniform kernels. The underlying uncontaminated gene expressions
#' are estimated by EM algorithm to maximize the data likelihood. Detailed
#' derivation can be found in our manuscript.
#'
#'
#' @examples
#'
#' data(MbrainSmall)
#' mbrain_obj <- CreateSlide(mbrain_raw,
#'                           mbrain_slide_info)
#' mbrain_decont_obj <- TBD(mbrain_obj, candidate_radius=20)
#' mbrain_decont_obj




#' @import Matrix
#' @import stats
#' @import SummarizedExperiment
#' @importFrom utils txtProgressBar
#' @importFrom dplyr filter
#'
#' @export

TBD <- function(slide_obj, gene_keep=NULL,
                maxit=30, tol=1,
                candidate_radius=5*1:6,
                verbose=TRUE){

    if(verbose){
        message(Sys.time(), " Start.")
    }

    #############################
    # Step 0: setup
    #############################

    # some universal variables
    raw_data <- assays(slide_obj)$raw  # raw data matrix
    slide <- metadata(slide_obj)$slide  # slide info
    n_spots <- ncol(raw_data)  # number of spots
    ts_idx <- which(slide$tissue==1)  # tissue index
    raw_ts_data <- raw_data[,ts_idx, drop=FALSE]  # raw tissue data matrix

    # Euclidean distance matrix for spots
    slide_distance <- .calculate_euclidean_weight(
        select(slide, imagerow, imagecol)
    )
    rownames(slide_distance) <- colnames(slide_distance) <- slide$barcode

    # some constant matrices
    I1_y <- matrix(1,n_spots,length(ts_idx))
    I1_yy <- matrix(1,length(ts_idx),length(ts_idx))
    I_yy <- diag(length(ts_idx))

    # calculate ARC score
    ARC_score <- ARCScore(raw_data, which(slide$tissue==0))

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
    mean_exp <- rowSums(raw_data)/length(ts_idx)
    if(is.null(gene_keep)){
        gene_keep <- rownames(KeepHighGene(raw_ts_data, verbose=verbose))

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

        slide_weight <- .gaussian_kernel(
            slide_distance, .points_to_sdv(candidate_radius[i], spot_distance)
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
    best_radius <- which.min(lapply(cont_out,function(x) x$value) %>% unlist)
    cont_rate <- cont_out[[best_radius]]$par[1]
    global_rate <- cont_out[[best_radius]]$par[2]
    cont_radius <- candidate_radius[best_radius]

    # contamination weight matrix
    slide_weight <- .gaussian_kernel(slide_distance, .points_to_sdv(cont_radius, spot_distance))
    slide_weight <- slide_weight[ts_idx,]/rowSums(slide_weight[ts_idx,])
    bleed_weight_mat <- cont_rate*(global_rate/n_spots+(1-global_rate)*slide_weight)

    #############################
    # Step 2: Decontaminate expression matrix
    #############################

    if(verbose){
        message("\n",Sys.time()," Decontaminating genes ...")
    }

    scale_factor <- rowSums(raw_data[gene_keep,])/rowSums(raw_ts_data[gene_keep,])
    decont_data_init <- raw_ts_data[gene_keep,]*scale_factor

    decont_out <- .decont_EM(raw_data=raw_data[gene_keep,],
                             decont_data_init = decont_data_init,
                             cont_rate=cont_rate,
                             global_rate=global_rate,
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

        scale_ratio <- rowSums(raw_data[gene_notkeep,, drop=F])/
            rowSums(raw_ts_data[gene_notkeep,, drop=F])
        scale_ratio[is.na(scale_ratio)] <- 1 # 0/0=NaN
        decont_data <- rbind(decont_data,
                             raw_ts_data[gene_notkeep,, drop=F]*scale_ratio)
    }

    decont_data <- as(decont_data[rownames(raw_data),],"dgCMatrix")

    # write results
    meta <- c(metadata(slide_obj), list(
        ARC_score=ARC_score,
        bleeding_rate=cont_rate,
        global_contamination_rate=global_rate,
        contamination_radius=cont_radius,
        weight_matrix=bleed_weight_mat,
        loglh=decont_out$loglh
    ))
    meta$slide <- filter(meta$slide, tissue==1)
    decont_obj <- CreateSlide(decont_data, meta,
                              gene_cutoff = 0, verbose=FALSE)
    assayNames(decont_obj) <- "decont"

    if(verbose){
        message(Sys.time(), " All finished.")
    }

    return(decont_obj)
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
                       cont_rate, global_rate,
                       slide_weight, ts_idx,
                       maxit=30, tol=1,
                       verbose=TRUE){
    # EM parameter updates
    # Args:
    #   raw_data (matrix of num): raw expression matrix
    #   decont_data_init (matrix of num): initial value of decontaminated
    #       expression matrix
    #   cont_rate (num): bleeding rate
    #   global_rate (num): global contamination rate
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
    bleed_weight_mat <- cont_rate*(global_rate/N_spot+(1-global_rate)*slide_weight)

    # Iteration
    repeat{
        if(verbose){
            message("\nIteration: ",n_it)
        }

        # E-step

        Eta <- decont_data%*%bleed_weight_mat
        Eta[,ts_idx] <- Eta[,ts_idx]+decont_data*(1-cont_rate)

        # stayed
        S <- (1-cont_rate)*raw_ts_data*decont_data/Eta[,ts_idx]
        # local contamination
        M <- (raw_data/Eta)%*%t(slide_weight)*
            decont_data*cont_rate*(1-global_rate)
        # global contamination
        N <- decont_data*cont_rate*global_rate/N_spot * rowSums(raw_data/Eta)

        # M-step

        decont_data_new <- S+M+N
        loglh <- sum(raw_data*log(Eta)-Eta, na.rm = T)

        # Difference between two adjacent iterations
        Lambda_maxdiff <- max(abs(decont_data_new-decont_data))

        decont_data <- decont_data_new
        Loglh <- c(Loglh, loglh)

        if(verbose){
            message("Log-likelihood: ",round(loglh,3),
                    "\nMax difference of decontaminated expressions: ",round(Lambda_maxdiff,3))
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
    cont_rate_init <- sum(obs_exp[bg_idx])/length(bg_idx)*
        n_spots/sum(obs_exp)

    bg_quant <- quantile(obs_exp[bg_idx], c(0.25, 0.5))
    bg_trim <- obs_exp[bg_idx]>=bg_quant[1] & obs_exp[bg_idx]<=bg_quant[2]
    global_cont <- mean(obs_exp[bg_idx][bg_trim])
    global_rate_init <- global_cont/sum(obs_exp[bg_idx])*length(bg_idx)
    # trim the initial value
    global_rate_init <- min(global_rate_init, 0.5)

    mu_init <- obs_exp[ts_idx][nonzero_pos]
    mu_init <- mu_init/sum(mu_init)*sum(obs_exp)
    x_init <- unname(c(cont_rate_init,global_rate_init,mu_init))
    x_init[1:2] <- pmax(x_init[1:2], 0.1)

    # calculate other matrices
    I1tZ <- crossprod(I1_y,obs_exp)
    WtZ <- crossprod(W_y,obs_exp)

    # assign parameter bounds
    lower_bounds <- rep(0,length(x_init))
    lower_bounds[1:2] <- c(sum(obs_exp[bg_idx])/sum(obs_exp),0.1)
    upper_bounds <- rep(Inf,length(x_init))
    upper_bounds[1:2] <- 1

    # minimize RSS using L-BFGS-B
    opt <- optim(x_init, .fn_optim, .gr_optim, method = "L-BFGS-B",
                 obs_exp=obs_exp, ts_idx=ts_idx,
                 nonzero_pos = nonzero_pos, n_spots=n_spots,
                 W_yy=W_yy, WtW=WtW, Wyy_tWyy=Wyy_tWyy,I_yy=I_yy,
                 I1_yy=I1_yy, WtZ=WtZ, I1tZ=I1tZ,
                 lower=lower_bounds,
                 upper=upper_bounds,
                 control = list(maxit=100))

    x_final <- numeric(length(ts_idx))
    x_final[nonzero_pos] <- opt$par[-(1:2)]
    opt$par <- c(opt$par[1:2],x_final)

    return(list(opt))
}


.fn_optim <- function(x, obs_exp, ts_idx, nonzero_pos,
                       W_yy, WtW, I_yy, Wyy_tWyy, I1_yy, n_spots, WtZ, I1tZ){
    # objective function: residual sum of squares

    x_coef <- x[-(1:2)]
    x_c_rate <- x[1] # cont_rate
    x_g_rate <- x[2] # global_rate

    x_coef%*%.AtA(x_c_rate, x_g_rate,
                  WtW[nonzero_pos,nonzero_pos],
                  I_yy[nonzero_pos,nonzero_pos],
                  Wyy_tWyy[nonzero_pos,nonzero_pos],
                  I1_yy[nonzero_pos,nonzero_pos],
                  n_spots)%*%x_coef-
        2*crossprod(.AtZ(x_c_rate, x_g_rate, WtZ[nonzero_pos], I1tZ[nonzero_pos],
                         obs_exp[ts_idx[nonzero_pos]],n_spots),x_coef)


}


.gr_optim <- function(x, obs_exp, ts_idx, nonzero_pos,
                       W_yy, WtW, I_yy, Wyy_tWyy, I1_yy, n_spots, WtZ, I1tZ){
    # Gradient of .fr_optim

    # recover full dimension of x
    N <- length(ts_idx)
    x_coef <- numeric(N)
    x_coef[nonzero_pos] <- x[-(1:2)]
    x_c_rate <- x[1] # cont_rate
    x_g_rate <- x[2] # global_rate


    x_coef1 <- x[-(1:2)]
    x_c_rate <- x[1] # cont_rate
    x_g_rate <- x[2] # global_rate


    # gradient of spots expressions
    g_coef <- 2*crossprod(.AtA(x_c_rate, x_g_rate, WtW[nonzero_pos,nonzero_pos],
                               I_yy[nonzero_pos,nonzero_pos],
                               Wyy_tWyy[nonzero_pos,nonzero_pos],
                               I1_yy[nonzero_pos,nonzero_pos], n_spots),x_coef1)-
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

    # gradient of global contamination rate
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

.AtA <- function(cont_rate, global_rate, WtW, I_yy, Wyy_tWyy,I1_yy, n_spots){
    (1-cont_rate)^2*I_yy+
        cont_rate^2*(1-global_rate)^2*WtW+
        cont_rate*(1-cont_rate)*(1-global_rate)*Wyy_tWyy+
        cont_rate*global_rate*(2-cont_rate*global_rate)/n_spots*I1_yy
}

.AtZ <- function(cont_rate, global_rate, WtZ, I1tZ, obs_ts_exp, n_spots){
    (1-cont_rate)*obs_ts_exp+cont_rate*(1-global_rate)*WtZ+cont_rate*global_rate/n_spots*I1tZ
}

.dAtA_dr <- function(cont_rate, global_rate, WtW, I_yy, Wyy_tWyy,I1_yy, n_spots){
    2*(cont_rate-1)*I_yy+2*cont_rate*(1-global_rate)^2*WtW+
        (1-2*cont_rate)*(1-global_rate)*Wyy_tWyy+
        2*(global_rate-cont_rate*global_rate^2)/n_spots*I1_yy
}

.dAtZ_dr <- function(cont_rate, global_rate, WtZ, I1tZ, n_spots, obs_ts_exp){
    (1-global_rate)*WtZ+global_rate/n_spots*I1tZ-obs_ts_exp

}

.dAtA_dc <- function(cont_rate, global_rate, WtW, Wyy_tWyy, n_spots, I1_yy){
    2*cont_rate^2*(global_rate-1)*WtW+
        cont_rate*(cont_rate-1)*Wyy_tWyy+
        2*(cont_rate-cont_rate^2*global_rate)/n_spots*I1_yy
}

.dAtZ_dc <- function(cont_rate, WtZ, I1tZ, n_spots){
    cont_rate/n_spots*I1tZ-cont_rate*WtZ

}
