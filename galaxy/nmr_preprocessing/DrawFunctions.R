
require(ggplot2)
require(gridExtra)
require(reshape2)


Draw <- function(Signal_data, type.draw = c("signal", "pca"), output = c("default", 
                                                                         "window", "png", "pdf"), dirpath = ".", filename = "%003d", height = 480, 
                 width = 640, pdf.onefile = TRUE, ...) {
  
  # Data initialisation and checks ----------------------------------------------
  type.draw <- match.arg(type.draw)
  output <- match.arg(output)
  fullpath <- paste(file.path(dirpath, filename), output, sep = ".")
  createFile <- TRUE
  createWindow <- FALSE
  
  # Drawing --------------------------------------------------------------------
  # output
  switch(output, default = {
    createFile <- FALSE
  }, window = {
    createWindow <- TRUE
    createFile <- FALSE
  }, png = {
    grDevices::png(fullpath, width, height)
  }, pdf = {
    grDevices::pdf(fullpath, width = width/72, height = height/72, 
                   onefile = pdf.onefile)
  }, {
    stop("Unknown output type.")
  })
  
  # Drawing type (signal/spectrum or PCA)
  funs <- list(signal = DrawSignal, pca = DrawPCA)
  if (type.draw %in% names(funs)) {
    fun <- funs[[type.draw]]
  } else {
    stop(paste("Unknown type:", type.draw))
  }
  
  # Plot finalisation ----------------------------------------------
  if (is.vector(Signal_data)) {
    Signal_data <- vec2mat(Signal_data)
  }
  fun(Signal_data, createWindow = createWindow, ...)
  if (createFile) {
    grDevices::dev.off()
  }
}



DrawSignal <- function(Signal_data, subtype = c("stacked", "together", 
                                                "separate", "diffmean", "diffmedian", "diffwith"), 
                       ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
                       xlab = "rowname", RowNames = NULL, row = 1, num.stacked = 4, 
                       main.title = NULL, createWindow) {
  
  # nticks
  
  # Data initialisation and checks ----------------------------------------------
  
  subtype <- match.arg(subtype)
  
  
  n <- nrow(Signal_data)
  m <- ncol(Signal_data)
  
  
  scale <- colnames(Signal_data)
  
  num.plot <- sum(ReImModArg)
  
  Var <- rowname <- value <- NULL  # only for R CMD check
  
  # Drawing array
  if (num.plot <= 0) {
    stop("Nothing selected in ReImModArg.")
  } else if (num.plot <= 2) {
    if (vertical)  {
      nrow <- num.plot
      ncol <- 1
    } else  {
      nrow <- 1
      ncol <- num.plot
    }
  } else {
    nrow <- 2
    ncol <- 2
  }
  
  # RowNames 
  if (is.null(RowNames))  {
    RowNames <- rownames(Signal_data)
    if (is.null(RowNames))  {
      RowNames <- 1:n
    }
  } else {
    if (!is.vector(RowNames)) {
      stop("RowNames is not a vector")
    }
    if (length(RowNames) != n)  {
      stop(paste("RowNames has length", length(RowNames), "and there are", n, "FIDs."))
    }
  }
  
  if (n == 1) {
    RowNames <- deparse(substitute(Signal_data))
  }
  
  elements <- list()
  if (ReImModArg[1]) {
    elements[["Re"]] <- Re(Signal_data)
  }
  if (ReImModArg[2]) {
    elements[["Im"]] <- Im(Signal_data)
  }
  if (ReImModArg[3]) {
    elements[["Mod"]] <- Mod(Signal_data)
  }
  if (ReImModArg[4]) {
    elements[["Arg"]] <- Arg(Signal_data)
  }
  
  
  
  
  # Drawing --------------------------------------------------------------------
  
  y = x = NULL # only for R CMD check
  
  
  # SEPARATE or STACKED ===============
  if (subtype == "separate" | subtype == "stacked")  {
    
    i <- 1
    while (i <= n)  {
      if (createWindow)  {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      if (subtype == "separate")  {
        # The other uses gridExtra to do that
        graphics::par(mfrow = c(nrow, ncol))
      }
      plots <- list()
      if (subtype == "separate")  {
        last <- i
      } else  {
        last <- min(i + num.stacked - 1, n)
      }
      for (name in names(elements))  {
        if (subtype == "separate")   {
          if (n == 1) {
            df <- data.frame(x = as.numeric(scale), y = elements[[name]])
          } else {df <- data.frame(x = as.numeric(scale), y = elements[[name]][i, ])
          }
          
          plots[[name]] <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y)) + 
            ggplot2::geom_line() + 
            ggplot2::theme(legend.position = "none") + 
            ggplot2::labs(x = xlab, y = name) +
            ggplot2::ggtitle(RowNames[i]) +
            ggplot2::theme_bw()
          
          if ((df[1, "x"] - df[(dim(df)[1]), "x"]) > 0) {
            plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
          }
          
        } else   {
          
          if (n == 1) {
            melted <- data.frame(rowname = rep(name, m), 
                                 Var = as.numeric(scale), value = elements[[name]][i,])
          } else {melted <- reshape2::melt(elements[[name]][i:last, ], 
                                           varnames = c("rowname", "Var"))
          }
          
          
          plots[[name]] <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value)) + 
            ggplot2::geom_line() + 
            ggplot2::facet_grid(rowname ~ ., scales = "free_y") + 
            ggplot2::theme(legend.position = "none") + 
            ggplot2::labs(x = xlab, y = name) +
            ggplot2::ggtitle(main.title) +
            ggplot2::theme_bw()
          
          if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 0) {
            plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
          }
        }
      }
      
      if (subtype == "stacked")  {
        do.call(gridExtra::grid.arrange, c(plots, list(nrow = nrow, ncol = ncol)))
      }
      i <- last + 1
    }
  } else if (subtype %in% c("together", "diffmean", "diffmedian", "diffwith")) {
    
    # TOGHETER or DIFFMEAN or DIFFMEDIAN or DIFFWITH ===============
    
    rainbow_colors <- grDevices::rainbow(n)
    
    if (createWindow) {
      grDevices::dev.new(noRStudioGD = TRUE)
    }
    graphics::par(mfrow = c(nrow, ncol))
    
    plots <- list()
    
    # Loop for Re, Im, Mod and Arg
    for (name in names(elements)) {
      # Get this part of the signal
      element <- elements[[name]]
      
      # Express the signal according to a reference if asked by `subtype'
      if (subtype == "diffmean")  {
        element <- sweep(element, MARGIN = 2, colMeans(element),  `-`)
      } else if (subtype == "diffmedian") {
        element <- sweep(element, MARGIN = 2, matrixStats::colMedians(element), `-`)
      } else if (subtype == "diffwith")  {
        element <- sweep(element, MARGIN = 2, element[row, ], `-`)
        if (row == 1 & n > 1)  {
          # Since we use plot on the first row and lines on the following, the y
          # scale is calculated at the first row so if the first row is all 0, it
          # causes problems
          tmp <- element[1, ]
          element[1, ] <- element[2, ]
          element[2, ] <- tmp
        }
      }
      
      
      melted <- reshape2::melt(elements[[name]], varnames = c("rowname", "Var"))
      
      
      plots[[name]] <- ggplot2::ggplot(melted, ggplot2::aes(x = Var, 
                                                            y = value, group = rowname, colour = rowname)) + ggplot2::geom_line() + 
        ggplot2::labs(x = xlab, y = name) + ggplot2::scale_colour_discrete(name = NULL) + 
        ggplot2::ggtitle(main.title)
      
      if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 
          0)  {
        plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
      }
      
      do.call(gridExtra::grid.arrange, c(plots, list(nrow = nrow, 
                                                            ncol = ncol)))
    }
  }
}

DrawPCA <- function(Signal_data, drawNames = TRUE, main = "PCA score plot", 
                    Class = NULL, axes = c(1, 2), type = c("scores", "loadings"), 
                    loadingstype = c("l", "p"), num.stacked = 4, xlab = "rowname", 
                    createWindow){
  
  
  # Data initialisation and checks ----------------------------------------------
  loadingstype <- match.arg(loadingstype)
  type <- match.arg(type)
  
  checkArg(main, "str", can.be.null = TRUE)
  
  # Class
  if (!is.vector(Class, mode = "any") & !is.null(Class)) {
    stop("Class is not a numeric vector")
  }
  if (is.vector(Class, mode = "numeric") & length(Class) != nrow(Signal_data)) {
    stop("the length of Class is not equal to the nrow of Signal_data")
  }
  
  # axes
  if (!is.vector(axes, mode = "numeric")) {
    stop("axes is not a numeric vector")
  }
  
  if (nrow(Signal_data) < 2) {
    stop("At least 2 spectra are needed for PCA.")
  }
  
  if (0 %in% Class) {
    Class <- Class + 1
  }
  
  # axes for scores plot
  Xax <- axes[1]
  Yax <- axes[2]
  
  
  # PCA ----------------------------------------------
  
  pca <- stats::prcomp(Re(Signal_data))
  
  # Eigenvalues
  eig <- (pca$sdev)^2
  
  # Variances in percentage
  variance <- eig * 100/sum(eig)
  
  # scores
  scores <- as.data.frame(pca$x)
  
  # loadings
  loadings <- matrix(data = NA, nrow = nrow(pca$rotation), ncol = ncol(pca$rotation))
  for (i in 1:length(eig)) {
    loadings[, i] <- pca$rotation[, i] * pca$sdev[i]
  }
  rownames(loadings) <- colnames(Signal_data)
  colnames(loadings) <- paste0("Loading", c(1:length(eig)))
  loadings <- as.data.frame(loadings)
  
  
  
  # Drawing ----------------------------------------------
  
  plots <- list()
  
  Var <- rowname <- value <- NULL  # only for R CMD check
  
  # SCORES ===============
  if (type == "scores") {
    
    if (!length(axes) == 2) {
      stop("the length of axes is not equal to 2 for scores plot")
    }
    
    if (createWindow)  {
      grDevices::dev.new(noRStudioGD = TRUE)
    }
    Xlim <- c(min(pca$x[, Xax]) * 1.4, max(pca$x[, Xax]) * 1.4)
    Ylim <- c(min(pca$x[, Yax]) * 1.4, max(pca$x[, Yax]) * 1.4)
    
    plots <- ggplot2::ggplot(scores, ggplot2::aes(get(colnames(scores)[Xax]), 
                                                  get(colnames(scores)[Yax]))) + ggplot2::xlim(Xlim) + ggplot2::ylim(Ylim)
    
    if (is.null(Class))  {
      plots <- plots + ggplot2::geom_jitter()
    } else  {
      plots <- plots + ggplot2::geom_jitter(ggplot2::aes(colour = Class, 
                                                         shape = Class))
    }
    
    plots <- plots + ggplot2::ggtitle(main) + ggplot2::geom_vline(xintercept = 0, 
                                                                  size = 0.1) + ggplot2::geom_hline(yintercept = 0, size = 0.1) + 
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60", 
                                                                                    size = 0.2), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "gray98")) + 
      ggplot2::labs(x = paste0("PC", Xax, " (", round(variance[Xax], 
                                                      2), "%)"), y = paste0("PC", Yax, " (", round(variance[Yax], 
                                                                                                   2), "%)"))
    
    if (drawNames)  {
      if (is.null(Class))  {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, 
                                                                    Xax], y = scores[, Yax], label = rownames(Signal_data)), 
                                            hjust = 0, nudge_x = (Xlim[2]/25), show.legend = FALSE, 
                                            size = 2)
      } else  {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, 
                                                                    Xax], y = scores[, Yax], label = rownames(Signal_data), 
                                                         colour = Class, shape = Class), hjust = 0, nudge_x = (Xlim[2]/25), 
                                            show.legend = F, size = 2)
      }
    }
    
    print(ggplot2::last_plot())
    
    
  } else {
    # LOADINGS ===============
    
    loadings <- loadings[, axes]
    
    if (is.vector(loadings)) {
      n <- 1
    } else {
      n <- ncol(loadings)
    }
    
    
    i <- 1
    while (i <= n)  {
      if (createWindow)  {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      
      last <- min(i + num.stacked - 1, n)
      
      if (n == 1)   {
        melted <- reshape2::melt(t(loadings), varnames = c("rowname",  "Var"))
      } else {
        melted <- reshape2::melt(t(loadings[, i:last]), varnames = c("rowname", 
                                                                     "Var"))
      }
      
      
      plots <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value))
      
      if (loadingstype == "l")  {
        plots <- plots + ggplot2::geom_line()
      } else if (loadingstype == "p")  {
        plots <- plots + ggplot2::geom_point(size = 0.5)
      } else  {
        warning("loadingstype is misspecified")
      }
      
      plots <- plots + ggplot2::ggtitle(main) + ggplot2::facet_grid(rowname ~  ., scales = "free_y") + 
        ggplot2::theme(legend.position = "none") + 
        ggplot2::labs(x = xlab, y = "Loadings") + 
        ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60") + 
        ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("(", round(variance[i:last], 
                                                                               2), "%)"), vjust = 1, hjust = 1)
      
      if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 0)  {
        plots <- plots + ggplot2::scale_x_reverse()
      }
      
      
      # Plot finalisation ----------------------------------------------
      
      i <- last + 1
      print(ggplot2::last_plot())
    }
    
  }
  
  
}
