#####   Helper functions
library(ggplot2)

save_LOVE_result <- function(clusters, delta, lbd, only.pure = F, data_name = "") {
  K <- length(clusters)
  if (only.pure)
    savePath = paste("../output/", data_name, format(Sys.time(),"/love_%m.%d_K_"),K,"_delta_",delta,"_lbd_",lbd,"_pure",sep='')
  else
    savePath = paste("../output/", data_name, format(Sys.time(),"/love_%m.%d_K_"),K,"_delta_",delta,"_lbd_",lbd,sep='')
  
  dir.create(file.path(getwd(),savePath))

  foreach (k = 1:K) %do% {
    groupk = sort(unlist(clusters[[k]], use.names = F))
    name = file.path(getwd(), savePath, paste("group",k,"-size-",length(groupk),
                                            ".txt", sep=''))
    write.table(groupk, file=name, sep="\t", col.name=FALSE, row.names=FALSE)
  }
  return(savePath)
}



Violin_plot <- function(methods, CV_result, MSE_flag, ylab_name = NULL, ylims = NULL, 
                        label_names = NULL) {
  col_ind <- 1:length(methods)
  if (MSE_flag) 
    data_wide_CV <- data.frame(1 - CV_result[ ,col_ind])
  else
    data_wide_CV <- data.frame(CV_result[ ,col_ind])
  
  names(data_wide_CV) <- methods
  
  data_long_CV <- reshape(data_wide_CV, direction = "long", timevar = "Method",
                          varying = methods, v.names = "MSE", times = methods)
  
  if (MSE_flag) {
    if (is.null(ylims))
      ylims <- c(0.75, max(data_long_CV$MSE))
    ylab_name <- ifelse(is.null(ylab_name), "1 - mean square error", ylab_name)
  } else {
    if (is.null(ylims))
      ylims <- c(0.25, 1)
    ylab_name <- ifelse(is.null(ylab_name), "Accuracy", ylab_name)
  }
  
  if (is.null(label_names))
    label_names <- methods
  
  data_long_CV$Method <- factor(data_long_CV$Method, levels = label_names)
  
  p1 <- ggplot(data_long_CV, aes(Method, MSE, color = Method, fill = Method)) +
    geom_violin(alpha = 0.5) + ylab(ylab_name) + ylim(ylims) + theme_classic()
 
  p1 
}


Fitted_plot <- function(Y, fitted_mat, methods, order) {
  y_ranges <- range(as.vector(fitted_mat), Y)
  par(mfrow=c(2,3))
  order_inc <- order(Y)
  corr <- c()
  for (k in 1:length(methods)) {
    i <- order[k]
    corr[i] <- cor(Y, fitted_mat[,i], method = "spearman")
    plot(Y[order_inc], fitted_mat[order_inc,i], ylim = range(range(Y), range(fitted_mat[,i])), col = i, pch = 1, xlim = range(Y),
         ylab = "predicted Y", xlab = "observed Y", main = paste(methods[i]," (corr = ", round(corr[i],2), ")", sep = ""))
  }
  names(corr) <- methods
}


  

Plot_2D <- function(Y_label, data_mat, label_names = c(1, 0), 
                    loc_legend = "topleft", axis_names = c("CR1", "CR2")) {

  res_svd <- svd(data_mat, 2, 2)
  p_comp <- res_svd$u %*% diag(res_svd$d[1:2])
  ind <- Y_label == 1
  
  par(mar = c(2.8, 2.8, 1, 1), mgp = c(1.8, 0.5, 0))
  
  if (length(label_names) == 3) {
    plot(p_comp[Y_label == 0,1], p_comp[Y_label == 0,2], col = 2, xlab = axis_names[1], 
         ylab = axis_names[2], main = "", type = "p", xlim = range(p_comp[,1]), 
         ylim = range(p_comp[,2]), pch = 16, cex = 1.2)
    points(p_comp[Y_label == 1, 1], p_comp[Y_label == 1, 2], col = 3, pch = 16, cex = 1.2)
    points(p_comp[Y_label == 2, 1], p_comp[Y_label == 2, 2], col = 4, pch = 16, cex = 1.2)
    legend(loc_legend, legend = label_names, col = 2:4, pch = 16, bty = "o", cex = 1.2,
           lty = 0)
  } else {
    plot(p_comp[ind,1], p_comp[ind,2], col = 2, xlab = axis_names[1], ylab = axis_names[2], main = "",
         type = "p", xlim = range(p_comp[,1]), ylim = range(p_comp[,2]), pch = 16, cex = 1.2)
    points(p_comp[!ind, 1], p_comp[!ind, 2], col = 3, pch = 16, cex = 1.2)
    legend(loc_legend, legend = label_names, col = 2:3, pch = 16, bty = "o", cex = 1.2,
           lty = 0)
  }
}



