#' Plot method (S4) for [Enrichment-class] object
#'
#' [plot] is a generic function for plotting of R objects.
#' The function invokes particular [methods] which depend on the [class] of the first argument.
#'
#' @param x An object of class [Enrichment-class] which the Z statistics or p-values have to be drawn.
#' @param what Default `"Genome"`. Plot Z statistics or p-values for the defined `what`,
#'   *i.e.*, `"All"`, `"Genome"` or numeric vector.
#' @param type Plot the selected analysis for `"eSNP"` and/or `"xSNP"`.
#' @param ggplot Use ggplot (default `FALSE`) instead of classic plot method.
#' @param pvalue If `TRUE`, p-value convergense is plotted. Otherwise, Z statistic is plotted.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see [par]).
#'
#' @examples
#' if (interactive()) {
#'   data(toyEnrichment)
#'   reSample(toyEnrichment, 10)
#'   plot(toyEnrichment)
#' }
#'
#' @name plot
#' @exportMethod plot
methods::setGeneric(name = "plot", def = function (x, y, ...) standardGeneric("plot"))
#' @name plot
#' @rdname plot
#' @aliases plot,Enrichment-method
methods::setMethod(f = "plot", signature = "Enrichment", definition = function(x, what = "Genome", type = c("eSNP", "xSNP"), ggplot = FALSE, pvalue = TRUE, ...) {
  if (is.null(unlist(x@Call$reSample, use.names = FALSE))) {
    stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
  }

  if (is.null(what) | any(!what %in% c("All", "Genome", seq_len(22)))) {
    stop('[Enrichment:plot] "what" must be: "All", "Genome" or numeric value (atomic or vector).', call. = FALSE)
  }
  if (is.null(type) | any(!type %in% c("eSNP", "xSNP"))) {
    stop('[Enrichment:plot] "type" must be: "eSNP" and/or "xSNP".', call. = FALSE)
  }
  if (any(type %in% "xSNP") & length(x["xSNP"]["List"]) == 0) {
    type <- "eSNP"
  }

  .computeER4plot <- function(EnrichSNPObject) {
    ER <- EnrichSNPObject@EnrichmentRatio
    if (nrow(EnrichSNPObject@Resampling) == 0) {
      stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
    }
    resampling <- EnrichSNPObject@Resampling[, 5]
    ERsample <- NULL
    size <- length(resampling)
    if (size >= 1000) {
      interv <- unique(c(seq(from = min(1000, floor(0.1 * size)), to = size, by = floor(size / 1000)), size))
    } else {
      interv <- unique(c(seq(from = max(floor(0.1 * size), 3), to = size, by = 1), size))
    }
    ERsample <- sapply(interv, function(k) {
      resamplingInterv <- resampling[1:k]
      resamplingClean <- resamplingInterv[!(is.na(resamplingInterv) | is.infinite(resamplingInterv))]
      mu <- mean(resamplingClean)
      sigma <- sqrt(stats::var(resamplingClean))

      if (sigma == 0 | is.na(sigma)) {
        if (mu == 0) 0 else (ER - mu)
      } else {
        (ER - mu) / sigma
      }
    })
    names(ERsample) <- interv
    as.matrix(ERsample)
  }
  .computeEmpP4plot <- function(EnrichSNPObject) {
    ER <- EnrichSNPObject@EnrichmentRatio
    if (nrow(EnrichSNPObject@Resampling) == 0) {
      stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
    }
    resampling <- EnrichSNPObject@Resampling[, 5]
    ERsample <- NULL
    size <- length(resampling)
    if (size >= 1000) {
      interv <- unique(c(seq(from = min(1000, floor(0.1 * size)), to = size, by = floor(size / 1000)), size))
    } else {
      interv <- unique(c(seq(from = max(floor(0.1 * size), 3), to = size, by = 1), size))
    }
    ERsample <- sapply(interv, function(k) {
      resamplingInterv <- resampling[1:k]
      resamplingClean <- resamplingInterv[!(is.na(resamplingInterv) | is.infinite(resamplingInterv))]
      sum(EnrichSNPObject@EnrichmentRatio < resamplingClean) / length(resamplingClean)
    })
    names(ERsample) <- interv
    as.matrix(ERsample)
  }

  if (x@Call$reSample$empiricPvalue) {
    matrixER <- list(eSNP = NULL, xSNP = NULL)
    for (iType in type) {
      if (length(what) == 1) {
        switch(EXPR = as.character(what),
          "Genome" = {
            matrixER[[iType]] <- .computeEmpP4plot(x[iType])
            colnames(matrixER[[iType]]) <- "Genome"
          },
          "All" = {
            matrixER[[iType]] <- .computeEmpP4plot(x[iType])
            for (j in seq_len(22)) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- c("Genome", paste0("Chrom", seq_len(22)))
          }, {
            for (j in what) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- paste0("Chrom", what)
          }
        )
      } else {
        for (j in what) {
          matrixER[[iType]] <- cbind(matrixER[[iType]], .computeEmpP4plot(x["Chromosomes", j][iType]))
        }
        colnames(matrixER[[iType]]) <- paste0("Chrom", what)
      }
    }
  } else {
    matrixER <- list(eSNP = NULL, xSNP = NULL)
    for (iType in type) {
      if (length(what) == 1) {
        switch(EXPR = as.character(what),
          "Genome" = {
            matrixER[[iType]] <- .computeER4plot(x[iType])
            colnames(matrixER[[iType]]) <- "Genome"
          },
          "All" = {
            matrixER[[iType]] <- .computeER4plot(x[iType])
            for (j in seq_len(22)) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- c("Genome", paste0("Chrom", seq_len(22)))
          }, {
            for (j in what) {
              matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
            }
            colnames(matrixER[[iType]]) <- paste0("Chrom", what)
          }
        )
      } else {
        for (j in what) {
          matrixER[[iType]] <- cbind(matrixER[[iType]], .computeER4plot(x["Chromosomes", j][iType]))
        }
        colnames(matrixER[[iType]]) <- paste0("Chrom", what)
      }
    }
  }
  if (ggplot) {
    if (!require("ggplot2") | !require("grid")) {
      stop('[Enrichment:plot] "ggPlot2" and "grid" packages must be installed with "ggplot=TRUE".', call. = FALSE)
    }
    multiplot <- function(..., plotlist = NULL, cols = 1, rows = 1, layout = NULL) {
      plots <- c(list(...), plotlist)
      numPlots <- length(plots)
      if (is.null(layout)) {
        layout <- matrix(seq_len(cols * ceiling(numPlots / cols)), ncol = cols, nrow = rows, byrow = TRUE)
      }
      if (numPlots == 1) {
        print(plots[[1]])
      } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
        }
      }
    }
    .ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
      if ((diff(h) %% 360) < 1) {
        h[2] <- h[2] - 360 / n
      }
      grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    listPlots <- list()
    for (iType in type) {
      if (pvalue) {
        if (x@Call$reSample$empiricPvalue) {
          ylab <- "P-Value (Empirical)"
        } else {
          matrixER[[iType]] <- apply(matrixER[[iType]], 2, stats::pnorm, lower.tail = FALSE)
          ylab <- "P-Value (From Z-statistic)"
        }
      } else {
        ylab <- "Z statistic"
      }
      if (ncol(matrixER[[iType]]) > 1) {
        matrixER[[iType]] <- apply(matrixER[[iType]], 2, scale)
      }
      tmpDF <- as.data.frame(t(matrixER[[iType]]))
      cnames <- colnames(tmpDF)
      colnames(tmpDF) <- paste0("R", colnames(tmpDF))
      tmpDF$IID <- factor(colnames(matrixER[[iType]]), levels = c("Genome", paste0("Chrom", seq_len(22))), labels = c("Genome", paste0("Chrom", seq_len(22))))
      tmp <- stats::reshape(tmpDF, idvar = "IID", direction = "long", varying = list(grep("R", colnames(tmpDF))), times = cnames)
      colnames(tmp) <- c("IID", "Resampling", "Z")
      tmp[, "Resampling"] <- as.numeric(tmp[, "Resampling"])

      p <- ggplot2::ggplot(data = tmp, ggplot2::aes_string(x = "Resampling", y = "Z", colour = "IID")) + ggplot2::geom_line()
      noGridColour <- "transparent"
      base_size <- 12
      base_family <- ""
      p <- p + ggplot2::theme(
        line = ggplot2::element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"),
        rect = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
        text = ggplot2::element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
        axis.text = ggplot2::element_text(size = ggplot2::rel(0.8), colour = "black"),
        strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(vjust = 1),
        axis.text.y = ggplot2::element_text(hjust = 1),
        axis.ticks = ggplot2::element_line(colour = "black"),
        axis.title.x = ggplot2::element_text(),
        axis.title.y = ggplot2::element_text(angle = 90),
        axis.ticks.length = ggplot2::unit(0.15, "cm"),
        axis.ticks.margin = ggplot2::unit(0.1, "cm"),
        legend.background = ggplot2::element_rect(fill = "white", colour = "black"),
        legend.margin = ggplot2::unit(0.2, "cm"),
        legend.key = ggplot2::element_rect(fill = "white", colour = "black"),
        legend.key.size = ggplot2::unit(1.2, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
        legend.text.align = NULL,
        legend.title = ggplot2::element_text(size = ggplot2::rel(0.8), face = "bold", hjust = 0),
        legend.title.align = NULL,
        legend.position = "right",
        legend.direction = NULL,
        legend.justification = "center",
        legend.box = NULL,
        panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = noGridColour[1]),
        panel.grid.minor = ggplot2::element_line(colour = noGridColour[length(noGridColour)], size = 0.25),
        panel.margin = ggplot2::unit(0.25, "lines"),
        strip.background = ggplot2::element_rect(fill = "black", colour = "black"),
        strip.text.x = ggplot2::element_text(colour = "white"),
        strip.text.y = ggplot2::element_text(angle = -90, colour = "white"),
        plot.background = ggplot2::element_rect(colour = "white"),
        plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)),
        plot.margin = ggplot2::unit(c(1, 1, 0.5, 0.5), "lines"),
        complete = TRUE
      )
      p <- p + ggplot2::xlab(paste(iType, "Resampling"))
      if (ncol(matrixER[[iType]]) > 1) {
        p <- p + ggplot2::ylab(paste(ylab, "(scale and center)"))
      } else {
        p <- p + ggplot2::ylab(ylab)
      }
      p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
      if (length(what) >= 5 | what == "All") {
        p <- p + ggplot2::theme(legend.background = ggplot2::element_rect(fill = "gray90", linetype = "dotted"))
      } else {
        quarter <- floor(3 / 4 * nrow(tmp)):nrow(tmp)
        rangeZtmp <- range(tmp[, "Z"])
        if (rangeZtmp[1] == rangeZtmp[2]) {
          if (all(rangeZtmp == 0)) {
            rangeZtmp[1] <- -1
            rangeZtmp[2] <- 1
          } else {
            rangeZtmp[1] <- rangeZtmp[1] * 0.90
            rangeZtmp[2] <- rangeZtmp[2] * 1.10
          }
        }
        rangeZ <- seq(rangeZtmp[1], rangeZtmp[2], by = diff(rangeZtmp) * 1 / 3)
        names(rangeZ) <- c("0%", "33%", "66%", "100%")
        rangeQuarter <- range(tmp[quarter, "Z"])
        inf <- apply(sapply(rangeQuarter, function(lim) lim < rangeZ), 1, all)
        sup <- apply(sapply(rangeQuarter, function(lim) lim > rangeZ), 1, all)
        if (sum(inf) <= sum(sup)) {
          p <- p + ggplot2::theme(legend.justification = c(1, 0), legend.position = c(1, 0))
        } else {
          p <- p + ggplot2::theme(legend.justification = c(1, 1), legend.position = c(1, 1))
        }
      }

      if ("Genome" %in% unique(tmp$IID)) {
        p <- p + ggplot2::scale_colour_manual(values = c("black", .ggplotColours(ifelse(length(unique(tmp$IID)) - 1 > 0, length(unique(tmp$IID)) - 1, 1))))
      }
      listPlots[[iType]] <- p
    }
    multiplot(plotlist = listPlots, cols = length(listPlots))
    invisible(listPlots)
  } else {
    graphics::par(mfrow = c(1, length(type)))
    for (iType in type) {
      if (pvalue) {
        if (x@Call$reSample$empiricPvalue) {
          matrixER[[iType]] <- .computeEmpP4plot(x[iType])
          ylab <- "P-Value (Empirical)"
        } else {
          matrixER[[iType]] <- apply(matrixER[[iType]], 2, stats::pnorm, lower.tail = FALSE)
          ylab <- "P-Value (From Z-statistic)"
        }
      } else {
        ylab <- "Z statistic"
      }
      nbCol <- ncol(matrixER[[iType]])
      ylim <- range(stats::na.exclude(matrixER[[iType]]))
      xNames <- rownames(matrixER[[iType]])
      colors <- grDevices::rainbow(nbCol)
      if (nbCol > 1) {
        matrixER[[iType]] <- apply(matrixER[[iType]], 2, scale)
        ylab <- paste(ylab, "(scale and center)")
        graphics::plot(x = xNames, y = matrixER[[iType]][, 1], ylab = ylab, xlab = iType, type = "l", ylim = ylim)
        res <- sapply(seq_len(ncol(matrixER[[iType]][, -1])), function(iER) {
          graphics::lines(x = xNames, y = matrixER[[iType]][, iER + 1], iType = "l", ylim = ylim, col = colors[iER + 1])
        })
      } else {
        graphics::plot(x = xNames, y = matrixER[[iType]][, 1], ylab = ylab, xlab = iType, type = "l", ylim = ylim)
      }
    }
    invisible()
  }
})
