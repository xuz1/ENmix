## Adopted from R package DunedinPACE. 
##Daniel W Belsky, Avshalom Caspi, David L Corcoran,et al. DunedinPACE, a DNA methylation biomarker of the pace of aging, eLIfe, 2022

DunedinPACE <-function (betas, proportionOfProbesRequired = 0.8) 
{
#    require(preprocessCore)
    mPOA_Models=list()
    load(system.file("mPOA_Models.RData",package="ENmix"))
        if (!is.numeric(as.matrix(betas))) {
            stop("betas matrix/data.frame is not numeric!")
        }
        probeOverlap <- length(which(rownames(betas) %in% mPOA_Models$model_probes))/length(mPOA_Models$model_probes)
        probeOverlap_background <- length(which(rownames(betas) %in% 
            mPOA_Models$gold_standard_probes))/length(mPOA_Models$gold_standard_probes)
        if (probeOverlap < proportionOfProbesRequired | probeOverlap_background < 
            proportionOfProbesRequired) {
            result <- rep(NA, ncol(betas))
            names(result) <- colnames(betas)
            model_results=result
            model_results
        }
        else {
            betas.mat <- as.matrix(betas[which(rownames(betas) %in% 
                mPOA_Models$gold_standard_probes), 
                ])
            probesNotInMatrix <- mPOA_Models$gold_standard_probes[which(mPOA_Models$gold_standard_probes %in% 
                rownames(betas.mat) == F)]
            if (length(probesNotInMatrix) > 0) {
                for (probe in probesNotInMatrix) {
                  tmp.mat <- matrix(0, nrow = 1, ncol = ncol(betas.mat))
                  rownames(tmp.mat) <- probe
                  colnames(tmp.mat) <- colnames(betas.mat)
                  tmp.mat[probe, ] <- rep(mPOA_Models$gold_standard_means[probe], 
                    ncol(tmp.mat))
                  betas.mat <- rbind(betas.mat, tmp.mat)
                }
            }
            samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 
                2, function(x) {
                  1 - (length(which(is.na(x)))/length(x)) < proportionOfProbesRequired
                }))]
            if (length(samplesToRemove) > 0) {
                betas.mat <- betas.mat[, -which(colnames(betas.mat) %in% 
                  samplesToRemove)]
            }
            if (ncol(betas.mat) > 0){
                pctValuesPresent <- apply(betas.mat, 1, function(x) {
                  1 - (length(which(is.na(x)))/length(x))
                })
                probesToAdjust <- which(pctValuesPresent < 1 & 
                  pctValuesPresent >= proportionOfProbesRequired)
                if (length(probesToAdjust) > 0) {
                  if (length(probesToAdjust) > 1) {
                    betas.mat[probesToAdjust, ] <- t(apply(betas.mat[probesToAdjust, 
                      ], 1, function(x) {
                      x[is.na(x)] = mean(x, na.rm = TRUE)
                      x
                    }))
                  }
                  else {
                    betas.mat[probesToAdjust, which(is.na(betas.mat[probesToAdjust, 
                      ]))] <- mean(betas.mat[probesToAdjust, 
                      ], na.rm = T)
                  }
                }
                if (length(which(pctValuesPresent < proportionOfProbesRequired)) > 
                  0) {
                  probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < 
                    proportionOfProbesRequired)]
                  for (probe in probesToReplaceWithMean) {
                    betas.mat[probe, ] <- rep(mPOA_Models$model_means[probe], 
                      ncol(betas.mat))
                  }
                }
                betas.norm <- normalize.quantiles.use.target(betas.mat, 
                  target = mPOA_Models$gold_standard_means)
                rownames(betas.norm) <- rownames(betas.mat)
                colnames(betas.norm) <- colnames(betas.mat)
                score = mPOA_Models$model_intercept + 
                  rowSums(t(betas.norm[mPOA_Models$model_probes, 
                    ]) %*% diag(mPOA_Models$model_weights))
                names(score) <- colnames(betas.norm)
                if (length(samplesToRemove) > 0) {
                  score.tmp <- rep(NA, length(samplesToRemove))
                  names(score.tmp) <- samplesToRemove
                  score <- c(score, score.tmp)
                }
                model_results <- score[colnames(betas)]
                model_results
            }
            else {
                model_results <- rep(NA, ncol(betas.mat))
                names(model_results) <- colnames(betas.mat)
                model_results
            }
        }
}
