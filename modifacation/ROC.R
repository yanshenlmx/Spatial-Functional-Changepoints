miss <- c()
for (phi in c(2,5)) {
  for (rho in c(1,1.5)) {
    for (null in c(12,25,38)) {
      for (seed in c(1:100)) {
        dir <- paste('/Users/yanshen/Library/CloudStorage/Box-Box/detection_simulation/d', phi, rho, null, seed, sep = '_')
        setwd(dir = dir)
        if (!file.exists("index_alt.RData")) {
          name <- paste('d', phi, rho, null, seed, sep = '_')
          miss <- c(miss, name)
        }
      }
    }
  }
}
## F1score
res <- NULL
for (phi in c(2,5)) {
  for (rho in c(1,1.5)) {
    for (null in c(12,25,38)) {
      scenarios <- c()
      for (seed in c(1:100)) {
        dir <- paste('/Users/yanshen/Library/CloudStorage/Box-Box/detection_simulation/d', phi, rho, null, seed, sep = '_')
        setwd(dir = dir)
        if (file.exists("index_alt.RData")) {
          name <- paste('d', phi, rho, null, seed, sep = '_')
          iters = 20000; thin_step = 10
          r_thin = seq(iters*3/4+1, iters, thin_step) 
          ns = 50
          
          load("index_alt.RData")
          load("pvalue.RData")
          p_temp = read.csv("p.csv", header = TRUE)[,-1]
          p_mean = apply(as.matrix(p_temp[r_thin,]), 2, mean)
          
          # Get the null locations
          index_null = NULL
          for(i in 1:ns){
            if(!i %in% index_alt){index_null = c(index_null, i)}
          }
          
          #### ROC for BH ####
          BH_ROC = data.frame(matrix(nrow = 0, ncol = 5))
          colnames(BH_ROC) = c("p_threshold", "sensitivity", "specificity", "precision", "recall")
          for(th_BH in seq(0, 1, by = 0.01)){
            alt_BH = which(p_mean > th_BH)
            null_BH = NULL
            for(i in 1:ns){
              if(!i %in% alt_BH){null_BH = c(null_BH, i)}
            }
            
            TN_BH = FP_BH = FN_BH = TP_BH = 0
            for(i in index_null){
              if(i %in% null_BH){TN_BH = TN_BH + 1
              }else{FP_BH = FP_BH + 1}
            }
            for(i in index_alt){
              if(i %in% null_BH){FN_BH = FN_BH + 1
              }else{TP_BH = TP_BH + 1}
            }
            BH_ROC = rbind(BH_ROC, data.frame(p_threshold = th_BH, 
                                              sensitivity = TP_BH/(TP_BH+FN_BH), 
                                              specificity = TN_BH/(FP_BH+TN_BH),
                                              precision = TP_BH/(TP_BH+FP_BH), 
                                              recall = TP_BH/(TP_BH+FN_BH),
                                              accuracy = (TP_BH + TN_BH)/50,
                                              F1score = 2*TP_BH/(2*TP_BH+FP_BH+FN_BH)))
          }
          
          #### ROC for FF ####
          FF_ROC = data.frame(matrix(nrow = 0, ncol = 5))
          colnames(FF_ROC) = c("p_threshold", "sensitivity", "specificity", "precision", "recall")
          for(th_FDR in seq(0, 1, by = 0.01)){
            alt_FF_FDR = which(p.adjust(pvalue, "BH") <= th_FDR)
            
            null_FF_FDR = NULL
            for(i in 1:ns){
              if(!i %in% alt_FF_FDR){null_FF_FDR = c(null_FF_FDR, i)}
            }
            
            TN_FF_FDR = FP_FF_FDR = FN_FF_FDR = TP_FF_FDR = 0
            for(i in index_null){
              if(i %in% null_FF_FDR){TN_FF_FDR = TN_FF_FDR + 1
              }else{FP_FF_FDR = FP_FF_FDR + 1}
            }
            for(i in index_alt){
              if(i %in% null_FF_FDR){FN_FF_FDR = FN_FF_FDR + 1
              }else{TP_FF_FDR = TP_FF_FDR + 1}
            }
            FF_ROC = rbind(FF_ROC, data.frame(p_threshold = th_FDR, 
                                              sensitivity = TP_FF_FDR/(TP_FF_FDR+FN_FF_FDR), 
                                              specificity = TN_FF_FDR/(FP_FF_FDR+TN_FF_FDR),
                                              precision = TP_FF_FDR/(TP_FF_FDR+FP_FF_FDR), 
                                              recall = TP_FF_FDR/(TP_FF_FDR+FN_FF_FDR),
                                              accuracy = (TP_FF_FDR + TN_FF_FDR)/50,
                                              F1score = 2*TP_FF_FDR/(2*TP_FF_FDR+FP_FF_FDR+FN_FF_FDR)))
          }
          
          BH_ROC_range <- BH_ROC[BH_ROC$p_threshold >= 0.4 & BH_ROC$p_threshold <= 0.6,]
          FF_ROC_range <- FF_ROC[FF_ROC$p_threshold >= 0 & FF_ROC$p_threshold <= 0.1,]
          
          F1score <- data.frame(FDR_BH = BH_ROC$F1score[which(BH_ROC$p_threshold == 0.4)], 
                                FDR_FF = FF_ROC$F1score[which(FF_ROC$p_threshold == 0.05)],
                                nature_BH = BH_ROC$F1score[which(BH_ROC$p_threshold == 0.5)],
                                nature_FF = FF_ROC$F1score[which(FF_ROC$p_threshold == 0.05)],
                                max_BH = max(BH_ROC_range$F1score),
                                max_FF = max(FF_ROC_range$F1score),
                                p_max_BH = BH_ROC_range$p_threshold[which.max(BH_ROC_range$F1score)][1],
                                p_max_FF = FF_ROC_range$p_threshold[which.max(FF_ROC_range$F1score)][1],
                                median_BH = median(BH_ROC_range$F1score),
                                median_FF = median(FF_ROC_range$F1score),
                                p_median_BH = BH_ROC_range$p_threshold[which(BH_ROC_range$F1score == median(BH_ROC_range$F1score))][1],
                                p_median_FF = FF_ROC_range$p_threshold[which(FF_ROC_range$F1score == median(FF_ROC_range$F1score))][1],
                                mean_BH = mean(BH_ROC$F1score[BH_ROC$p_threshold >= 0.4 & BH_ROC$p_threshold <= 0.6]),
                                mean_FF = mean(FF_ROC$F1score[FF_ROC$p_threshold >= 0 & FF_ROC$p_threshold <= 0.1]))
          scenarios <- rbind(scenarios, F1score)
        }
      }
      mean_var <- rbind(apply(scenarios, 2, mean), 
                        apply(scenarios, 2, sd)/(sqrt(nrow(scenarios))), 
                        apply(scenarios, 2, min),
                        apply(scenarios, 2, max))
      rownames(mean_var) <- c(paste(phi, rho, null, 'mean',sep = '_'), 
                              paste(phi, rho, null, 'se',sep = '_'), 
                              paste(phi, rho, null, 'pmin',sep = '_'),
                              paste(phi, rho, null, 'pmax',sep = '_'))
      res <- rbind(res, mean_var)
    }
  }
}

## ROC

res <- NULL
for (phi in c(2,5)) {
  for (rho in c(1,1.5)) {
    for (null in c(12,25,38)) {
      scenarios_BH_speci <- c()
      scenarios_BH_sensi <- c()
      scenarios_FF_speci <- c()
      scenarios_FF_sensi <- c()
      for (seed in c(1:100)) {
        dir <- paste('/Users/yanshen/Library/CloudStorage/Box-Box/detection_simulation/d', phi, rho, null, seed, sep = '_')
        setwd(dir = dir)
        if (file.exists("index_alt.RData")) {
          name <- paste('d', phi, rho, null, seed, sep = '_')
          iters = 20000; thin_step = 10
          r_thin = seq(iters*3/4+1, iters, thin_step) 
          ns = 50
          
          load("index_alt.RData")
          load("pvalue.RData")
          p_temp = read.csv("p.csv", header = TRUE)[,-1]
          p_mean = apply(as.matrix(p_temp[r_thin,]), 2, mean)
          
          # Get the null locations
          index_null = NULL
          for(i in 1:ns){
            if(!i %in% index_alt){index_null = c(index_null, i)}
          }
          
          #### ROC for BH ####
          BH_ROC = data.frame(matrix(nrow = 0, ncol = 5))
          colnames(BH_ROC) = c("p_threshold", "sensitivity", "specificity")
          for(th_BH in seq(0, 1, by = 0.01)){
            alt_BH = which(p_mean > th_BH)
            null_BH = NULL
            for(i in 1:ns){
              if(!i %in% alt_BH){null_BH = c(null_BH, i)}
            }
            
            TN_BH = FP_BH = FN_BH = TP_BH = 0
            for(i in index_null){
              if(i %in% null_BH){TN_BH = TN_BH + 1
              }else{FP_BH = FP_BH + 1}
            }
            for(i in index_alt){
              if(i %in% null_BH){FN_BH = FN_BH + 1
              }else{TP_BH = TP_BH + 1}
            }
            BH_ROC = rbind(BH_ROC, data.frame(p_threshold = th_BH, 
                                              sensitivity = TP_BH/(TP_BH+FN_BH), 
                                              specificity = TN_BH/(FP_BH+TN_BH)))
          }
          BH_sensi <- BH_ROC$sensitivity
          BH_speci <- BH_ROC$specificity
          #### ROC for FF ####
          FF_ROC = data.frame(matrix(nrow = 0, ncol = 5))
          colnames(FF_ROC) = c("p_threshold", "sensitivity", "specificity")
          for(th_FDR in seq(0, 1, by = 0.01)){
            alt_FF_FDR = which(p.adjust(pvalue, "BH") <= th_FDR)
            
            null_FF_FDR = NULL
            for(i in 1:ns){
              if(!i %in% alt_FF_FDR){null_FF_FDR = c(null_FF_FDR, i)}
            }
            
            TN_FF_FDR = FP_FF_FDR = FN_FF_FDR = TP_FF_FDR = 0
            for(i in index_null){
              if(i %in% null_FF_FDR){TN_FF_FDR = TN_FF_FDR + 1
              }else{FP_FF_FDR = FP_FF_FDR + 1}
            }
            for(i in index_alt){
              if(i %in% null_FF_FDR){FN_FF_FDR = FN_FF_FDR + 1
              }else{TP_FF_FDR = TP_FF_FDR + 1}
            }
            FF_ROC = rbind(FF_ROC, data.frame(p_threshold = th_FDR, 
                                              sensitivity = TP_FF_FDR/(TP_FF_FDR+FN_FF_FDR), 
                                              specificity = TN_FF_FDR/(FP_FF_FDR+TN_FF_FDR)))
          }
          FF_sensi <- FF_ROC$sensitivity
          FF_speci <- FF_ROC$specificity
          scenarios_BH_sensi <- cbind(scenarios_BH_sensi, BH_sensi)
          scenarios_BH_speci <- cbind(scenarios_BH_speci, BH_speci)
          scenarios_FF_sensi <- cbind(scenarios_FF_sensi, FF_sensi)
          scenarios_FF_speci <- cbind(scenarios_FF_speci, FF_speci)
        }
      }
      mean_BH <- cbind(apply(scenarios_BH_sensi, 1, mean), 
                        apply(scenarios_BH_speci, 1, mean))
      mean_FF <- cbind(apply(scenarios_FF_sensi, 1, mean), 
                       apply(scenarios_FF_speci, 1, mean))
      setwd("~/Desktop/Postdoc/changepoint detection/modifacation")
      png(paste("roc_", phi, "_", rho, "_", null, "_", ".png", sep = ""))
      plot(1 - mean_BH[,2], mean_BH[,1], type = "l", col = "orange", 
           xlab = "1-specificity", ylab = "sensitivity", lwd = 2, 
           main = paste("roc_", phi, "_", rho, "_", null, sep = ""))
      lines(1 - mean_FF[,2], mean_FF[,1], type = "l", col = "blue", lwd = 2)
      legend('bottomright', legend = c('BH', 'FF'), col = c('orange', 'blue'), lty = c(1))
      dev.off()
    }
  }
}


res_round <- round(res, 3)

res_se <- sapply(seq(1,23,2), function(xx){
  for (j in c(1:10)) {
    res_round[xx,j] <- paste(res_round[xx,j],'(', res_round[xx+1,j], ')', sep = "")
  }
  return(res_round[xx,])
})

res_se <- t(res_se)

library(xtable)
print(xtable(res_se, type = "latex", digits = c(0,0,0,0,0,0,0,0,0,0,0)), file = "F1score.tex")


iters = 20000; thin_step = 10
r_thin = seq(iters*3/4+1, iters, thin_step) 
ns = 50

load("index_alt.RData")
load("pvalue.RData")
p_temp = read.csv("p.csv", header = TRUE)[,-1]
p_mean = apply(as.matrix(p_temp[r_thin,]), 2, mean)

# Get the null locations
index_null = NULL
for(i in 1:ns){
  if(!i %in% index_alt){index_null = c(index_null, i)}
}

#### ROC for BH ####
BH_ROC = data.frame(matrix(nrow = 0, ncol = 5))
colnames(BH_ROC) = c("p_threshold", "sensitivity", "specificity", "precision", "recall")
for(th_BH in seq(0, 1, by = 0.01)){
  alt_BH = which(p_mean > th_BH)
  null_BH = NULL
  for(i in 1:ns){
    if(!i %in% alt_BH){null_BH = c(null_BH, i)}
  }
  
  TN_BH = FP_BH = FN_BH = TP_BH = 0
  for(i in index_null){
    if(i %in% null_BH){TN_BH = TN_BH + 1
    }else{FP_BH = FP_BH + 1}
  }
  for(i in index_alt){
    if(i %in% null_BH){FN_BH = FN_BH + 1
    }else{TP_BH = TP_BH + 1}
  }
  BH_ROC = rbind(BH_ROC, data.frame(p_threshold = th_BH, 
                                    sensitivity = TP_BH/(TP_BH+FN_BH), 
                                    specificity = TN_BH/(FP_BH+TN_BH),
                                    precision = TP_BH/(TP_BH+FP_BH), 
                                    recall = TP_BH/(TP_BH+FN_BH),
                                    accuracy = (TP_BH + TN_BH)/50,
                                    F1score = 2*TP_BH/(2*TP_BH+FP_BH+FN_BH)))
}

#### ROC for FF ####
FF_ROC = data.frame(matrix(nrow = 0, ncol = 5))
colnames(FF_ROC) = c("p_threshold", "sensitivity", "specificity", "precision", "recall")
for(th_FDR in seq(0, 1, by = 0.01)){
  alt_FF_FDR = which(p.adjust(pvalue, "BH") <= th_FDR)
  
  null_FF_FDR = NULL
  for(i in 1:ns){
    if(!i %in% alt_FF_FDR){null_FF_FDR = c(null_FF_FDR, i)}
  }
  
  TN_FF_FDR = FP_FF_FDR = FN_FF_FDR = TP_FF_FDR = 0
  for(i in index_null){
    if(i %in% null_FF_FDR){TN_FF_FDR = TN_FF_FDR + 1
    }else{FP_FF_FDR = FP_FF_FDR + 1}
  }
  for(i in index_alt){
    if(i %in% null_FF_FDR){FN_FF_FDR = FN_FF_FDR + 1
    }else{TP_FF_FDR = TP_FF_FDR + 1}
  }
  FF_ROC = rbind(FF_ROC, data.frame(p_threshold = th_FDR, 
                                    sensitivity = TP_FF_FDR/(TP_FF_FDR+FN_FF_FDR), 
                                    specificity = TN_FF_FDR/(FP_FF_FDR+TN_FF_FDR),
                                    precision = TP_FF_FDR/(TP_FF_FDR+FP_FF_FDR), 
                                    recall = TP_FF_FDR/(TP_FF_FDR+FN_FF_FDR),
                                    accuracy = (TP_FF_FDR + TN_FF_FDR)/50,
                                    F1score = 2*TP_FF_FDR/(2*TP_FF_FDR+FP_FF_FDR+FN_FF_FDR)))
}

F1score <- rbind(F1score, data.frame(FDR_BH = BH_ROC$F1score[which(BH_ROC$p_threshold == 0.4)], 
                                     FDR_FF = FF_ROC$F1score[which(FF_ROC$p_threshold == 0.05)],
                                     nature_BH = BH_ROC$F1score[which(BH_ROC$p_threshold == 0.5)],
                                     nature_FF = FF_ROC$F1score[which(FF_ROC$p_threshold == 0.05)],
                                     max_BH = max(BH_ROC$F1score[BH_ROC$p_threshold >= 0.4 & BH_ROC$p_threshold <= 0.6]),
                                     max_FF = max(FF_ROC$F1score[FF_ROC$p_threshold >= 0 & FF_ROC$p_threshold <= 0.1]),
                                     median_BH = median(BH_ROC$F1score[BH_ROC$p_threshold >= 0.4 & BH_ROC$p_threshold <= 0.6]),
                                     median_FF = median(FF_ROC$F1score[FF_ROC$p_threshold >= 0 & FF_ROC$p_threshold <= 0.1]),
                                     mean_BH = mean(BH_ROC$F1score[BH_ROC$p_threshold >= 0.4 & BH_ROC$p_threshold <= 0.6]),
                                     mean_FF = mean(FF_ROC$F1score[FF_ROC$p_threshold >= 0 & FF_ROC$p_threshold <= 0.1])))

F1score
library(xtable)
print(xtable(F1score, type = "latex", digits = c(0,4,4,4,4,4,4,4,4,4,4)), file = "F1score.tex")







