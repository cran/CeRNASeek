ceRNA.cernia <-
function(miRtar,targetce=NULL, geneexp,miRexp, mres, numMIR = 3, cor_cutoff = 0, s_cutoff = 0.5) {
# functions (parMM, graphWeights, recommendation, dtHybrid) of cernia method are from
# the website: https://github.com/dsardina/cernia 
#Copyright 2016 Rosalba Giugno Licensed under
# the Apache License, Version 2.0 (the 'License'); you may not use this file except in
# compliance with the License. You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed under the
# License is distributed on an 'AS IS' BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing permissions and
# limitations under the License.

    dtHybrid <- function(miRtar) {
    mir <- unique(miRtar[, 1])
    tar <- unique(miRtar[, 2])
    A <- matrix(nrow = length(tar), ncol = length(mir), data = 0)
    colnames(A) <- mir
    rownames(A) <- tar
    for (i in seq_len(nrow(miRtar))) {
        A[which(tar %in% as.character(miRtar[i, 2])),
            which(mir %in% as.character(miRtar[i, 1]))] <- 1
    }
    cl <- makeCluster(getOption("cl.cores", 2)) 
    #cl <- makeCluster(detectCores() - 2)
    M <- recommendation(A, cl = cl)
    W <- graphWeights(nrow(M), ncol(M), M, cl = cl)

    stopCluster(cl)

    return(W)
}

graphWeights <- function(n, m, A, lambda = 0.5, alpha = 0.5, S = NA, S1 = NA, cl = NA) {

    if (nrow(A) != n || ncol(A) != m) {
        stop("The matrix A should be an n by m matrix.")
    }

    has.similarity <- (!all(is.na(S)) && is.matrix(S) && !all(is.na(S1)) && is.matrix(S1))

    if (has.similarity) {
        if (nrow(S1) != m || ncol(S1) != m) {
            stop("The matrix S1 should be an m by m matrix.")
        }
        if (nrow(S) != n || ncol(S) != n) {
            stop("The matrix S should be an n by n matrix.")
        }
    }

    Ky <- diag(1/colSums(A))
    Ky[is.infinite(Ky) | is.na(Ky)] <- 0  

    kx <- rowSums(A)
    Nx <- 1/(matrix(kx, nrow = n, ncol = n, byrow = TRUE)^(lambda) * matrix(kx, nrow = n, ncol = n,
        byrow = FALSE)^(1 - lambda))

    Nx[is.infinite(Nx) | is.na(Nx)] <- 0  
    kx[is.infinite(kx) | is.na(kx)] <- 0 

    W <- t(parMM(cl, A, Ky))
    W <- parMM(cl, A, W)
    W <- Nx * W
    rownames(W) <- rownames(A)
    colnames(W) <- rownames(A)

    if (has.similarity) {
        X5 <- parMM(cl, A, S1)
        X6 <- parMM(cl, X5, t(A))
        X7 <- parMM(cl, A, matrix(1, nrow = m, ncol = m))
        X8 <- parMM(cl, X7, t(A))
        S2 <- X6/X8
        W <- W * (1 + (alpha * S) + ((1 - alpha) * S2))
    }

    W[is.nan(W)] <- 0 
    return(W)
}

recommendation <- function(A, lambda = 0.5, alpha = 0.5, S = NA, S1 = NA, cl = NA) {

    n <- nrow(A)
    m <- ncol(A)
    W <- graphWeights(n = n, m = m, A = A, lambda = lambda, alpha = alpha, S = S, S1 = S1, cl = cl)

    R <- parMM(cl, W, A)
    return(R)
}

parMM <- function(cl, A, B) {

    if (!all(is.na(cl)) && is.object(cl)) {

        nA <- nrow(A)
        ncl <- length(cl)
        i <- seq_len(nA)
        if (ncl == 1) {
            splitIndices <- i
        } else {
            fuzz <- min((nA - 1)/1000, 0.4 * nA/ncl)
            breaks <- seq(1 - fuzz, nA + fuzz, length = ncl + 1)
            splitIndices <- structure(split(i, cut(i, breaks)), names = NULL)
        }
        splitRows <- lapply(splitIndices, function(i) A[i, , drop = FALSE])

        fun <- get(as.character("rbind"), envir = .GlobalEnv, mode = "function")
        R <- do.call("fun", lapply(clusterApply(cl = cl, x = splitRows, get("%*%"), B), enquote))
    } else {
        R <- A %*% B
    }
    return(R)
}

if(is.null(targetce)){
    r1 <- nrow(miRtar)
    n1 <- ncol(miRtar)

    geneexpNames <- row.names(geneexp)
	miRexpNames <- row.names(miRexp)
    miRtar <- as.matrix(miRtar)
    miRtar <- miRtar[intersect(which(miRtar[, 1] %in% miRexpNames),which(miRtar[, 2] %in% geneexpNames)), ]
      
    miRNA <- miRtar[, 1]
    target <- miRtar[, 2]

    miR_uni <- unique(miRNA)
    tar_uni <- unique(target)

     t1<- length(tar_uni)

    tmp <- matrix(NA, t1 * (t1 - 1)/2, 2)
    tmp1 <- matrix(NA, t1 * (t1 - 1)/2, 10)

    for (i in seq_len(t1 - 1)) {
        for (j in seq(i + 1, t1)) {

            miRNA1 <- miRtar[which(miRtar[, 2] %in% tar_uni[i]), 1]
            miRNA2 <- miRtar[which(miRtar[, 2] %in% tar_uni[j]), 1]
            inter <- intersect(miRNA1, miRNA2)
            comres <- mres[mres[, 2] %in% c(tar_uni[i], tar_uni[j]) & mres[, 1] %in% inter, ]

            M1 <- length(miRNA1)
            M2 <- length(miRNA2)
            M3 <- length(intersect(miRNA1, miRNA2))
            
            

            index1 <- which(geneexpNames %in% tar_uni[i])
            index2 <- which(geneexpNames %in% tar_uni[j])

            M6 <- cor.test(as.numeric(geneexp[index1, ]),as.numeric( geneexp[index2, ]))$estimate

            if (M3 >= numMIR  & M6 > cor_cutoff &
                nrow(comres) > 0) {

                tmp[(i - 1) * t1 + j - sum(seq_len(i)), 1] <- tar_uni[i]
                tmp[(i - 1) * t1 + j - sum(seq_len(i)), 2] <- tar_uni[j]

                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 1] <- M3
                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 2] <- log(M3/min(M1, M2))
                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 3] <- sum(sapply(inter, function(x) {
                  MREs <- comres[comres[, 1] == x, ]
                  if (nrow(MREs) <= 0) return(1)
                  gap_score <- MREs[, c("gap_l", "gap_r")]
                  L <- abs(max(gap_score[, 1]) - min(gap_score[, 2]))
                  return(log(nrow(MREs)/L))
                }))
                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 4] = sum(sapply(inter, function(x) {
                  pos <- comres[comres[, 1] == x, c("gap_l", "gap_r")]
                  if (nrow(pos) <= 0) return(1)
                  return(log(abs(max(pos[, 1]) - min(pos[, 2]))^2/sum((pos[, 2] -
                    pos[, 1])^2)))
                }))
                nr <- nrow(comres)
                if (nr == length(unique(comres[, 1]))) {
                  tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 5] <- log(1/nr)
                } else {
                  tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 5] <- log((nr - length(unique(comres[,
                    1])) + 1)/nr)
                }

                comres <- mres[mres[, 2] %in% c(tar_uni[i], tar_uni[j]) & mres[, 1] %in%
                  inter, ]
                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 6] <- sum(sapply(inter, function(x) {
                  MREs <- comres[comres[, 1] == x, ]
                  if (nrow(MREs) <= 0) return(1)
                  gap_score <- MREs[, c("gap_l", "gap_r")]
                  L <- abs(max(gap_score[, 1]) - min(gap_score[, 2]))
                  return(log(sum(abs(MREs[, 3]))/L))
                }))
                recomm_score <- dtHybrid(miRtar)
                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 7] <- recomm_score[tar_uni[i],
                  tar_uni[j]]

                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 8] <- log(M6)

                tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 9] <- tmp1[(i - 1) * t1 + j - sum(seq_len(i)),
                  2] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 3] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)),
                  4] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 5] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)),
                  6] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)), 7] + tmp1[(i - 1) * t1 + j - sum(seq_len(i)),
                  8]
            }
        }
    }

    tmp <- tmp[apply(tmp, 1, function(x) !all(is.na(x))), ]
    tmp1 <- tmp1[apply(tmp1, 1, function(x) !all(is.na(x))), ]
    tmp1[, 10] <- (tmp1[, 9] - min(tmp1[, 9]))/(max(tmp1[, 9]) - min(tmp1[, 9]))
	tmp1[, 10]<-tmp1[, 9]
    tmp <-  tmp[(tmp1[,10] > s_cutoff), ]

    tmp1 <- tmp1[( tmp1[, 10] > s_cutoff) , ]

     if (is.vector(tmp1)) {
        result <- c(tmp, tmp1)
        names(result) <- c("targetce", "anotherce", "miRNAs_num",
            "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
            "Normalized score")
    } else {
        result <- cbind(tmp, tmp1)
        colnames(result) <- c("targetce", "anotherce", "miRNAs_num",
            "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
            "Normalized score")
	
    }
	
	return(result)
  }
  
  
else{
    r1 <- nrow(miRtar)
    n1 <- ncol(miRtar)

    geneexpNames <- row.names(geneexp)
	miRexpNames <- row.names(miRexp)
    miRtar <- as.matrix(miRtar)
    miRtar <- miRtar[intersect(which(miRtar[, 1] %in% miRexpNames),which(miRtar[, 2] %in% geneexpNames)), ]
      
    miRNA <- miRtar[, 1]
    target <- miRtar[, 2]

    miR_uni <- unique(miRNA)
    tar_uni <- unique(target)

     t1<- length(tar_uni)

    tmp <- matrix(NA, t1 , 2)
    tmp1 <- matrix(NA, t1, 10)
	
    for (j in seq_len(t1)) {

            miRNA1 <- miRtar[which(miRtar[, 2] %in% targetce), 1]
            miRNA2 <- miRtar[which(miRtar[, 2] %in% tar_uni[j]), 1]
            inter <- intersect(miRNA1, miRNA2)
            comres <- mres[mres[, 2] %in% c(targetce, tar_uni[j]) & mres[, 1] %in% inter, ]

            M1 <- length(miRNA1)
            M2 <- length(miRNA2)
            M3 <- length(intersect(miRNA1, miRNA2))
            
            

            index1 <- which(geneexpNames %in% targetce)
            index2 <- which(geneexpNames %in% tar_uni[j])

            M6 <- cor.test(as.numeric(geneexp[index1, ]),as.numeric( geneexp[index2, ]))$estimate

            if (M3 >= numMIR  & M6 > cor_cutoff &
                nrow(comres) > 0) {

                tmp[j, 1] <- targetce
                tmp[j, 2] <- tar_uni[j]

                tmp1[j, 1] <- M3
                tmp1[j, 2] <- log(M3/min(M1, M2))
                tmp1[j, 3] <- sum(sapply(inter, function(x) {
                  MREs <- comres[comres[, 1] == x, ]
                  if (nrow(MREs) <= 0) return(1)
                  gap_score <- MREs[, c("gap_l", "gap_r")]
                  L <- abs(max(gap_score[, 1]) - min(gap_score[, 2]))
                  return(log(nrow(MREs)/L))
                }))
                tmp1[j, 4] = sum(sapply(inter, function(x) {
                  pos <- comres[comres[, 1] == x, c("gap_l", "gap_r")]
                  if (nrow(pos) <= 0) return(1)
                  return(log(abs(max(pos[, 1]) - min(pos[, 2]))^2/sum((pos[, 2] -
                    pos[, 1])^2)))
                }))
                nr <- nrow(comres)
                if (nr == length(unique(comres[, 1]))) {
                  tmp1[j, 5] <- log(1/nr)
                } else {
                  tmp1[j, 5] <- log((nr - length(unique(comres[,
                    1])) + 1)/nr)
                }

                comres <- mres[mres[, 2] %in% c(targetce, tar_uni[j]) & mres[, 1] %in%
                  inter, ]
                tmp1[j, 6] <- sum(sapply(inter, function(x) {
                  MREs <- comres[comres[, 1] == x, ]
                  if (nrow(MREs) <= 0) return(1)
                  gap_score <- MREs[, c("gap_l", "gap_r")]
                  L <- abs(max(gap_score[, 1]) - min(gap_score[, 2]))
                  return(log(sum(abs(MREs[, 3]))/L))
                }))
                recomm_score <- dtHybrid(miRtar)
                tmp1[j, 7] <- recomm_score[targetce,tar_uni[j]]

                tmp1[j, 8] <- log(M6)

                tmp1[j, 9] <- tmp1[j,2] + tmp1[j, 3] + tmp1[j,4] + tmp1[j, 5] + tmp1[j,6] + tmp1[j, 7] + tmp1[j,8]
            }
        }
    

    tmp <- tmp[apply(tmp, 1, function(x) !all(is.na(x))), ]
    tmp1 <- tmp1[apply(tmp1, 1, function(x) !all(is.na(x))), ]
    tmp1[, 10] <- (tmp1[, 9] - min(tmp1[, 9]))/(max(tmp1[, 9]) - min(tmp1[, 9]))
    tmp <-  tmp[(tmp1[,10] > s_cutoff), ]

    tmp1 <- tmp1[( tmp1[, 10] > s_cutoff) , ]

     if (is.vector(tmp1)) {
        result <- c(tmp, tmp1)
        names(result) <- c("targetce", "anotherce", "miRNAs_num",
            "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
            "Normalized score")
    } else {
        result <- cbind(tmp, tmp1)
        colnames(result) <- c("targetce", "anotherce", "miRNAs_num",
            "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
            "Normalized score")
	
    }
	
	    return(result)
  }

}
