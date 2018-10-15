printf <- function (...) { print(sprintf(...)) }

plotMyloess <-
function (lo, add = T, col = "#990066", pred.int = TRUE, conf.int = TRUE, 
    ...) 
{
    if (!add) {
        plot(lo$x, lo$y, ...)
    }
    x <- c(lo$x, rev(lo$x))
    if (conf.int) 
        polygon(x, c(lo$ci.lower, rev(lo$ci.upper)), border = F, 
            col = paste0(col, "66"))
    if (pred.int) 
        polygon(x, c(lo$pr.lower, rev(lo$pr.upper)), border = F, 
            col = paste0(col, "33"))
    lines(lo$x, lo$y, lwd = 2, col = col)
}


myloess.sd <-
function (x, y = NULL, level = 0.95, ...) 
{
    level <- 1 - (1 - level)/2
    xy <- xy.coords(x, y)
    x <- xy$x
    x0 <- unique(sort(x))
    y <- xy$y
    sel <- !is.na(x) & !is.na(y)
    if (sum(!sel) > 0) {
        warning(sprintf("Removing %d pairs containing NA", sum(!sel)))
        x <- x[sel]
        y <- y[sel]
        x0 <- x0[sel]
    }
    mod <- loess(y ~ x, ...)
    yfit <- predict(mod, data.frame(x = x0))
    r <- residuals(mod)
    modr <- loess(I(r^2) ~ x, ...)
    sd <- sqrt(pmax(0, predict(modr, data.frame(x = x0))))
    lo.conf <- predict(mod, newdata = x0, se = T)
    lo.conf$upper <- lo.conf$fit + qt(level, lo.conf$df) * lo.conf$se.fit
    lo.conf$lower <- lo.conf$fit - qt(level, lo.conf$df) * lo.conf$se.fit
    list(model = mod, x = x0, y = yfit, sd = sd, pr.upper = yfit + 
        qnorm(level) * sd, pr.lower = yfit - qnorm(level) * sd, 
        ci.upper = lo.conf$upper, ci.lower = lo.conf$lower)
}


fig_label <-
function (text, region = "figure", pos = "topleft", cex = NULL, 
    ...) 
{
    region <- match.arg(region, c("figure", "plot", "device"))
    pos <- match.arg(pos, c("topleft", "top", "topright", "left", 
        "center", "right", "bottomleft", "bottom", "bottomright"))
    if (region %in% c("figure", "device")) {
        ds <- dev.size("in")
        x <- grconvertX(c(0, ds[1]), from = "in", to = "user")
        y <- grconvertY(c(0, ds[2]), from = "in", to = "user")
        if (region == "figure") {
            fig <- par("fig")
            dx <- (x[2] - x[1])
            dy <- (y[2] - y[1])
            x <- x[1] + dx * fig[1:2]
            y <- y[1] + dy * fig[3:4]
        }
    }
    if (region == "plot") {
        u <- par("usr")
        x <- u[1:2]
        y <- u[3:4]
    }
    sw <- strwidth(text, cex = cex) * 60/100
    sh <- strheight(text, cex = cex) * 60/100
    x1 <- switch(pos, topleft = x[1] + sw, left = x[1] + sw, 
        bottomleft = x[1] + sw, top = (x[1] + x[2])/2, center = (x[1] + 
            x[2])/2, bottom = (x[1] + x[2])/2, topright = x[2] - 
            sw, right = x[2] - sw, bottomright = x[2] - sw)
    y1 <- switch(pos, topleft = y[2] - sh, top = y[2] - sh, topright = y[2] - 
        sh, left = (y[1] + y[2])/2, center = (y[1] + y[2])/2, 
        right = (y[1] + y[2])/2, bottomleft = y[1] + sh, bottom = y[1] + 
            sh, bottomright = y[1] + sh)
    old.par <- par(xpd = NA)
    on.exit(par(old.par))
    text(x1, y1, text, cex = cex, ...)
    return(invisible(c(x, y)))
}


rocformat <-
function (x) 
{
    ret <- sprintf("AUC=%.2f (95%% CI=%.2f-%.2f)", x$auc, x$ci[1], 
        x$ci[3])
    return(ret)
}


rocplot <- 
function (ret, title = "ROC", file = NULL, confmat = TRUE, ci = TRUE, 
    positive.name = NULL, add = FALSE, ...) 
{
    require(pROC)
    if (!"data.frame" %in% class(ret)) 
        stop("ret must be a data.frame object")
    if (!all(c("reality", "prediction") %in% colnames(ret))) 
        stop("ret must have the following column names: reality, prediction")
    if (!"decision" %in% colnames(ret)) {
        if (confmat) 
            warning("No decision column; setting confmat=FALSE")
        confmat <- FALSE
    }
    if (is.null(positive.name)) 
        positive.name <- ret$reality[1]
    rr <- roc(response = (ret$reality == positive.name), predictor = ret$prediction, 
        plot = T, ci = T, main = title, add = add, ...)
    if (ci && !add) {
        cat(sprintf("# AUC= %.2f, CI= %.2f - %.2f\\n", rr$auc, 
            rr$ci[1], rr$ci[3]))
        legend("topleft", sprintf("AUC= %.2f\\nCI= %.2f - %.2f", 
            rr$auc, rr$ci[1], rr$ci[3]), bty = "n")
    }
    if (confmat && !add) {
        pp <- par("family")
        par(family = "mono")
        legend("bottomright", confuMat(ret$reality, ret$decision, 
            as.text = T), bty = "n", inset = 0.05)
        par(family = pp)
    }
    if (!is.null(file)) 
        copyan(file = file)
    invisible(rr)
}


mypalette <-
function (n = NULL, transparent = "99", alpha = NULL) 
{
    pal <- "E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7 999999 E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7"
    pal <- unlist(strsplit(pal, " "))
    if (!is.null(alpha)) {
        if (alpha > 1) 
            alpha <- 1
        if (alpha < 0) 
            alpha <- 0
        transparent <- sprintf("%02X", as.integer(alpha * 255))
    }
    pal <- paste("#", pal, transparent, sep = "")
    if (!is.null(n)) {
        if (n > length(pal)) {
            pal <- rep(pal, ceiling(n/length(pal)))
        }
        else {
            pal <- pal[1:n]
        }
    }
    return(pal)
}


plotROCseries <- 
function (prs, legend = TRUE, legend.names = NULL, add.stats = TRUE, 
    col = NULL, lwd = 3, lty = NULL, ...) 
{
    if (!is(prs, "list")) 
        stop("prs must be a list")
    if (is.null(col)) 
        col <- mypalette()
    N <- length(prs)
    add <- FALSE
    if (is.null(lwd)) 
        lwd <- 3
    if (length(lwd) == 1) 
        lwd <- rep(lwd, N)
    if (is.null(lty)) 
        lty <- 1
    if (length(lty) == 1) 
        lty <- rep(lty, N)
    if (length(col) == 1) 
        col <- rep(col, N)
    if (is.null(legend.names)) 
        legend.names <- names(prs)
    if (is.null(legend.names)) 
        legend.names <- 1:N
    rocs <- list()
    for (i in 1:N) {
        rocs[[i]] <- rocplot(prs[[i]], col = col[i], add = add, 
            lwd = lwd, lty = lty[i], confmat = FALSE, ci = FALSE, 
            ...)
        if (!add) 
            add <- TRUE
    }
    ciinf <- sapply(rocs, rocformat)
    if (add.stats) {
        legend.names <- paste(legend.names, ciinf, sep = ", ")
    }
    if (legend) {
        legend("bottomright", legend.names, lwd = lwd, col = col, 
            bty = "n", lty = lty)
    }
    return(invisible(rocs))
}


topTableAll <- 
function (fit, contrasts = NULL, genelist = NULL, idcoll = NULL, 
    sign.fc.col = FALSE, sig.threshold = 1e-05, msd = FALSE, 
    ci = FALSE, rawp = FALSE) 
{
    require(limma)
    if (is.null(contrasts)) {
        contrasts <- colnames(fit$coefficients)
    }
    common <- colnames(fit$genes)
    if (!missing(genelist)) {
        common <- union(common, colnames(genelist))
    }
    else {
        genelist <- fit$genes
    }
    ret <- NULL
    for (c in contrasts) {
        tt <- topTable(fit, coef = c, genelist = genelist, sort.by = "none", 
            number = Inf, confint = msd)
        if (msd | ci) {
            if (any(grepl("CI.025", colnames(tt)))) {
                tt$CI.L <- tt$CI.025
                tt$CI.R <- tt$CI.975
            }
            tt$msd <- apply(sign(tt$logFC) * cbind(tt$CI.L, tt$CI.R), 
                1, min)
        }
        if (is.null(ret)) {
            print(colnames(tt))
            ret <- tt[, common]
        }
        ret[, paste("logFC", c, sep = ".")] <- tt$logFC
        if (sign.fc.col) {
            ret[, paste("logFC.s", c, sep = ".")] <- tt$logFC
            ret[tt$adj.P.Val > sig.threshold, paste("logFC.s", 
                c, sep = ".")] <- 0
        }
        if (msd) {
            ret[, paste("msd", c, sep = ".")] <- tt$msd
        }
        if (ci) {
            ret[, paste("cil", c, sep = ".")] <- tt$CI.L
            ret[, paste("cir", c, sep = ".")] <- tt$CI.R
        }
        if (rawp) {
            ret[, paste("P.Value", c, sep = ".")] <- tt$P.Value
        }
        ret[, paste("qval", c, sep = ".")] <- tt$adj.P.Val
    }
    if (!is.null(idcoll)) {
        rownames(ret) <- ret[, idcoll]
    }
    return(ret)
}


pval2star <- 
function (pv, ns = "", thresholds = c(0.05, 0.001, 1e-04)) 
{
    ret <- as.character(pv)
    thresholds <- c(thresholds, 0)
    strings <- c(ns, "*", "**", "***")
    for (i in length(thresholds):1) ret[pv > thresholds[i]] <- strings[i]
    ret
}


## older version of stratify
stratify.old <- function(df, ...) {
  strata <- list(...)

  f <- strata[[1]]
  f <- factor(f)
  catf("Stratifying by %s\n", names(strata)[1])
  ret <- sapply(levels(f), function(x) df[ f == x, ], simplify=FALSE)

  if(length(strata) > 1) {
    ret <- sapply(ret, function(r) {
      params <- c(list(df=r), strata[-1])
      do.call(stratify.old, params)
    }, simplify=FALSE)

    ret <- unlist(ret, recursive=FALSE)
  }
  ret
}

#' Stratify a vector, crossing each pair of strata
#' @param x a vector 
#' @param strata a list of strata. Each stratum is a list explicitly defining the stratum either by enumeration or as subscripts / logical values
#' @param by.which if true, then the strata are defined by subscripts / logical values
#' @example
#' st <- list(one=list(AJ=1:15, KT=16:20), two=list(odd=seq(1,20,2), foo=5:15, even=seq(2,20,2)))
#' stratify(LETTERS[1:30], st)
stratify <- function(x, strata, by.which=TRUE) {

  catf("Stratifying by %s\n", names(strata)[1])

  if(by.which) {
    strata <- sapply(strata, function(s) {
      sapply(s, function(ss) x[ss], simplify=FALSE)
    }, simplify=FALSE)
  }

  f <- strata[[1]]

  ret <- sapply(f, function(y) intersect(x, y), simplify=FALSE)

  if(length(strata) > 1) {
    ret <- sapply(ret, function(r) {
      stratify(r, strata[-1], by.which=FALSE)
    }, simplify=FALSE)
    ret <- unlist(ret, recursive=FALSE)
  }
  ret
}


## shorthand to return the AUC calculated by pROC
auc <- function(x) {
  require(pROC)
  pval <- auc.pval(x)
  ret <- c(roc(response=x$reality, predictor=x$prediction, ci=TRUE)$ci[1:3], pval)
}

auc.pval <- function(x, only.pval=TRUE) {
  ret <- wilcox.test(x$prediction ~ x$reality)
  if(only.pval) ret <- ret$p.value
  ret
}

auc.txt <- function(x, type=1) {

  sprintf("%.2f (95%% CI: %.2f–%.2f)", x[2], x[1], x[3])
}

auctxt <- function(x, f="%.2f; 95%% CI: %.2f–%.2f") {
    vals <- auc(x)
      sprintf(f, vals[2], vals[1], vals[3])
}

auctxtb1 <- function(x, f="%.2f (95%% CI: %.2f–%.2f)") {
    vals <- auc(x)
      sprintf(f, vals[2], vals[1], vals[3])
}

auctxtb2 <- function(x, f="%.2f [95%% CI: %.2f–%.2f]") {
    vals <- auc(x)
      sprintf(f, vals[2], vals[1], vals[3])
}

auctxtb0 <- function(x, f="%.2f") {
    vals <- auc(x)
      sprintf(f, vals[2])
}




## remove rows which have less (or equal) than N columns with value higher than minimal
## if 0 < N < 1, treat N as a proportion
filter.by.min <- function(x, d=NULL, n=0.2) {
  if(n < 1 && n > 0) {
    n <- round(n * ncol(x))
  }
  counts <- apply(x, 1, function(xx) sum(xx > min(xx)))

  if(!is.null(d)) {
    return(d[ counts >= n, ])
  } 
  x[ counts >= n, ]
}

## remove rows which contain only one value
remove.identical.rows <- function(x, d=NULL) {
  counts <- apply(x, 1, function(xx) length(unique(xx)))

  if(!is.null(d)) {
    return(d[ counts > 1,, drop=F ])
  } 

  x[ counts > 1,, drop=F ]
}

## rank-normalize an EList
EList.ranknorm <- function(rg) {
  ret <- rg
  ret$E <- t(apply(rg$E, 1, rank))
  ret
}

## create an EList, make sure col/row names are OK
make.EList <- function(E, genes, targets) {
  require(limma)
  rg <- new("EList",
    list(E=E,
      genes=genes[rownames(E),],
      targets=targets[colnames(E),]))
  colnames(rg) <- rownames(rg$targets) <- colnames(E)
  rownames(rg) <- rownames(rg$genes) <- rownames(E)
  rg
}

## plot lines grouped by a factor
tcplot <- function(x, y, f, col, lwd=2, ...) {

  uids <- unique(f)

  for(id in uids) {
    sel <- f == id
    xx <- x[sel]
    yy <- y[sel]
    o  <- order(xx)
    scol <- col[sel][1]
    lines(xx[o], yy[o], lwd=lwd, col=scol)
  }

}


fcfc.plot <- function(fit, coef1, coef2, xlab=NULL, ylab=NULL,
  pch=19, col="#33333333", show.conf=FALSE, add=FALSE, xlim=NULL, ylim=NULL, ...) {
  if(!all(c(coef1, coef2) %in% colnames(fit$coefficients))) stop("Incorrect coefficients")
  tt <- topTableAll(fit, msd=T, ci=TRUE)

  if(is.null(xlab)) xlab <- coef1
  if(is.null(ylab)) ylab <- coef2


  l1 <- tt[,paste0("logFC.", coef1)]
  l2 <- tt[,paste0("logFC.", coef2)]

  if(is.null(xlim)) xlim <- range(l1)
  if(is.null(ylim)) ylim <- range(l2)


  if(!add) {
    plot( NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    all_ablines()
  }

  if(show.conf) {
    segments(
      l1, tt[,paste0("cil.", coef2)],
      l1, tt[,paste0("cir.", coef2)],
      col="#33333333")

    segments(
      tt[,paste0("cil.", coef1)], l2,
      tt[,paste0("cir.", coef1)], l2,
      col="#33333333")
  }
  
  points(l1, l2, pch=pch, col=col) 
}

## fit a standardized model to data
model.fit <- function(rg, selection, formula, contrasts) {
  ret <- list()

  ret$rg <- rg[ , selection]
  ret$d  <- model.matrix(as.formula(formula), data=ret$rg$targets)
  if(!is.null(contrasts)) {
    ret$c  <- makeContrasts(contrasts=contrasts, levels=ret$d)
    colnames(ret$c) <- names(contrasts)
    ret$fit <- eBayes(contrasts.fit(lmFit(ret$rg, ret$d), ret$c))
  } else {
    ret$fit <- eBayes(lmFit(ret$rg, ret$d))
  }
  ret$tt  <- topTableAll(ret$fit, ci=T, msd=T)
  ret
}

## shortcut for topTable
topTT <- function(fit, coef=1, number=Inf, p.value=0.05, ...) {
  topTable(fit, coef=coef, p.value=p.value, number=number, ...)
}

## shorthand to quickly display a biochemical
loess.plot.bioch <- function(i, rg, data=NULL, ...) {
  if(!is.null(data)) y <- data[i,]
  else               y <- rg$E[i,]
  main <- all.biochemicals[i, "BIOCHEMICAL"]
  substr(main, 1, 1) <- toupper(substr(main, 1, 1))

  loess.plot(rg, y, main=main, ...)
}

## plot function for p-val effects plot
makePFfunc <- function(E, P) {
  colP <- paste0(colorRampPalette(c("grey", "red"))(5)[2:5], "99", sep="")
  colN <- paste0(colorRampPalette(c("grey", "blue"))(5)[2:5], "99", sep="")

  steps <- -log10(c(1, 0.5, 0.01, 0.001))


  E2 <- abs(E)
  E2 <- (E2 - min(E2)) / (max(E2) - min(E2))
  #E2 <- log(E2 + 2*min(E2))

  pf <- function(row, col, x, y, col.w, row.h, e, p) { 
    #printf("r=%s c=%s x=%.2f y=%.2f w=%.2f e=%.2f P=%.2f", row, col, x, y, col.w, e, p) ;
    if(is.na(e)) {
      color <- "#33333333"
    } else {
      if(E[row, col] > 0) color <- colP[ findInterval(p, steps) ]
      else color <- colN[ findInterval(p, steps) ]
      print(color)
    }
    points(x, y, col=color, pch=19, cex=1 + E2[row, col] * 2) 

    }

pf

}

loess.plot <- function(rg, y, 
  ylab="Relative abundance", xlab="Time to TB", bty="n", 
  show.points=TRUE,
  plot.ctrls=FALSE,
  plot.cases=TRUE,
  ylim=NULL, legend=TRUE,
  ...) {
  x <- -rg$targets$TIME_TO_TB
  cases <- rg$targets$CASE_CONTROL == "case"
  ctrls <- rg$targets$CASE_CONTROL == "control"

  #x <- -rg$targets$TIME_TO_TB
  #y <- mod[i, ]
  #y <- rg$E[i,]
  xyml.cases <- myloess.sd(x[cases], y[cases], degree=1)
  xyml.ctrls <- myloess.sd(x[ctrls], y[ctrls], degree=1)

  if(is.null(ylim)) {
    ylim <- c( min(c(xyml.cases$pr.lower, xyml.ctrls$pr.lower)),
             max(c(xyml.cases$pr.upper, xyml.ctrls$pr.upper)))
  }
  plot(NULL, ylim=ylim, xlim=range(x), xlab=xlab, ylab=ylab, bty=bty, ...)

  if(show.points) {
    points(x[ctrls], y[ctrls], pch=19, col="#00999933")
    points(x[cases], y[cases], pch=19, col="#99009933")
  }

  #plotMyloess(xyml.ctrls, col="#009999", pred.int=FALSE)
  plotMyloess(xyml.cases, col="#990099", pred.int=FALSE)

  abline(h=median(y[ctrls]), col="#009999", lwd=2)
  sdctrl <- IQR(y[ctrls])
  abline(h=quantile(y[ctrls], c(0.25, 0.75)), col="#009999", lwd=1, lty=2)

  if(legend) {
    legend("topleft", c( "Progressors", "Controls"), pch=19, col=c("#990099", "#009999"), bty="n")
  }

}



## row-wise wilcoxon test with correction for multiple testing
run.wilcox <- function(d, group1, group2, genelist=NULL, adj.method="fdr") {

  pvals <- apply(d, 1, function(x) {
    wilcox.test(x[group1], x[group2])$p.value
  })

  padj <- p.adjust(pvals, method=adj.method)

  ret <- data.frame(P.Value=pvals, adj.P.Val=padj)
  if(!is.null(genelist)) {
    ret <- data.frame(genelist, ret)
  }
  ret
}

## Create a ROC curve with additional data, based on pROC roc() function
rocfigure <- function( ret, title= "ROC", file= NULL, confmat= TRUE, ci= TRUE ) {

  rr <- roc( response= ( ret$reality == "case" ), predictor= ret[,1], plot= T, ci= T, main= title )

  if( ci )  {
    cat( sprintf( "# AUC= %.2f, CI= %.2f - %.2f\n", rr$auc, rr$ci[1], rr$ci[3] ) )
    legend( "topleft", sprintf( "AUC= %.2f\nCI= %.2f - %.2f", rr$auc, rr$ci[1], rr$ci[3] ), bty= "n" )
  }

  if( confmat ) {
    pp <- par( "family" )
    par( family= "mono" )
    legend( "bottomright", confuMat( ret[,4], ret[,3], as.text= T ), bty= "n", inset= 0.05 )
    par( family= pp )
  }

  if( ! is.null( file ) ) copyan( file= file )

  invisible( rr )
}

## like xpredict, but constructs a secondary, reduced model and reports
## only AUC for each of the 100 replicates
xpredict.nvar <- function(x, y, x.response, y.response, rep=100,
  nvar=c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50, 75, 100, 125))  {

  set.x <- rownames(x)
  set.y <- rownames(y)

  common <- intersect(set.x, set.y)

  x <- x[common,,drop=F]
  y <- y[common,,drop=F]
  print(dim(x))
  print(dim(y))

  res <- sapply(1:rep, function(i) {
    catf("\r%d", i)
    rf <- randomForest(t(x), factor(x.response), importance=T)
    var.l <- rownames(rf$importance)[order(-rf$importance[,4])]
    sapply(nvar, function(n) {
      sel <- var.l[1:n]
      x2 <- x[sel,]
      rf2 <- randomForest(t(x2), factor(x.response))
      pr <- predict(rf2, newdata=t(y), type="prob")
      roc(response=y.response, predictor=pr[,1])$auc[1]
    })
  })

  cat("\n")
  res
}


## given two matrices, train the model on the first, apply to the second
xpredict <- function(x, y, x.response, y.response, nvar=NULL, 
  cleanx=FALSE, cleany=FALSE, cutoff.method="MaxProdSpSe", 
  return.list=FALSE,
  rf=NULL,
  ...) {
  require(randomForest)
  require(OptimalCutpoints)

  set.x <- rownames(x)
  set.y <- rownames(y)

  common <- intersect(set.x, set.y)

  x <- x[common,,drop=F]
  y <- y[common,,drop=F]
  n.vars <- length(common)

  if(cleanx) 
    x <- x[ , !colnames(x) %in% colnames(y),drop=F]

  if(cleany)
    y <- y[ , !colnames(y) %in% colnames(x),drop=F]

  if(!is.factor(x.response))
    x.response <- factor(x.response)

  if(is.null(rf)) 
    rf <- randomForest(t(x), x.response, importance=T)

  if(!is.null(nvar)) {
    sel <- order(-rf$importance[,4])[1:nvar]
    x <- x[sel,]
    rf <- randomForest(t(x), x.response, importance=T)
    n.vars <- nvar

  }



  ret <- data.frame(predict(rf, newdata=t(y), type="prob"))

  cutpoints <-
    optimal.cutpoints("pr", "re", 
    data=data.frame(pr=rf$votes[,1], re=x.response),
    tag.healthy=rf$classes[2],
    methods=cutoff.method)
    
  cutoff <- 
    cutpoints[[cutoff.method]][[1]]$optimal.cutoff$cutoff[1]

  attr(ret, "cutoff") <- cutoff
  attr(ret, "nvar") <- n.vars
  ret$decision <- colnames( ret )[ ( ret[,1] < ret[,2] ) + 1 ]
  ret$decision2 <- colnames(ret)[ (ret[,1] < cutoff) + 1 ]
  ret$reality  <- y.response
  ret$prediction <- ret[,1]
  if(return.list) {
    return(list(ret=ret, rf=rf, cutoff=cutoff, cutpoints=cutpoints))
  } else {
    return(ret)
  }
}

## merge matrices using common biochemicals
mergesets <- function(x, ...) {
  inp <- list(...)

  for(y in inp) {
    if(ncol(y) > 0) {
      common <- intersect(rownames(x), rownames(y))
      x <- x[common,]
      y <- y[common,]
      x <- cbind(x, y)
    }
  }

  x
}

## plot a series of ROC prediction objects as boxplot
seriesboxplot <- function(x, xlim=c(0.25, 1), rcol="#ccddcc", bty="n", mm=0.4, ...) {
  require(pROC)
  oldpar <- par(mar= c(5,1,4,9)+0.1)
  on.exit(par(oldpar))

  s <- sapply(x, function(xx) as.vector(roc(response=xx$reality, predictor=xx$prediction, ci=T)$ci))
  print(s)
  n <- ncol(s)
  y <- -1:-n

  plot(x=NULL, y=NULL, xlim=xlim, ylim=c(-(n+0.5), -0.5), bty=bty, yaxt="n", ylab="", xlab="AUC", ...) 
  segments(s[1,], y, s[3,], y, col="grey", lwd=4)

  #rect((1:n)-mm, s[1,], (1:n)+mm, s[3,], col=rcol)
  abline(v=0.5, col="#993333", lwd=3)
  abline(v=(6:10)/10, col="grey", lwd=1)
  #segments((1:n)-mm, s[2,], (1:n)+mm, s[2,], lwd=2)
  axis(4, at=y, labels=names(x), tick=F, las=2, cex.axis=1.5, mgp=c(0,0,0))
  points(s[2,], y, pch=18, cex=1.4)
}

## based on prediction objects, returns a table with auc, ci and p-value
## method: method to adjust for multiple testing
# Prv: prevalence assumed as in Suliman et al. (2%)
# Spc: specificity according to TIPP (75%)
rocseriestable <- function(x, method="fdr", levels=c("control", "case"), 
	Spc=0.75, Prv=0.02) {
  require(pROC)

  rocs <- sapply(x, function(xx) as.vector(roc(response=xx$reality, predictor=xx$prediction, levels=levels, ci=T)),
		simplify=FALSE)
	s <- sapply(rocs, function(x) x$ci[1:3])
  #s <- sapply(x, function(xx) as.vector(roc(response=xx$reality, predictor=xx$prediction, ci=T)$ci))
  s <- t(s)
  pvals <- sapply(x, auc.pval)

  ret <- data.frame(name=names(x), AUC=s[,2], CI.L=s[,1], CI.R=s[,3], p.value=pvals)
  ret$q.value <- p.adjust(ret$p.value, method=method)
	
  # calculating ppv, npv and confidence intervals
	foo <- sapply(rocs, function(x) {
		#bar <- ci.coords(x, x="best", ret=c("ppv", "npv"), best.method="c", best.policy="random")
		#bar0 <- coords(x, x="best", ret=c("ppv", "npv"), best.method="c")
		thr <- coords(x, x="best", best.method="y", ret="threshold")[1]
		sens <- coords(x, x=Spc, input="specificity", ret="sensitivity")[1]
		tab <- table(x$predictor < thr, x$response)
		bar <- summary(epi.tests(tab))[ c("ppv", "npv"), c(2,1,3) ]
		
		npv2 <- Spc * (1 - Prv) / ( (1-sens) * Prv + Spc * (1- Prv) )
		ppv2 <- sens * Prv / (sens * Prv + (1 - Spc) * (1- Prv))
		return(c(PPV=bar[1,2], PPV.CI.L=bar[1,1], PPV.CI.R=bar[1,3],
		         NPV=bar[2,2], NPV.CI.L=bar[2,1], NPV.CI.R=bar[2,3], SENS=sens, 
						 PPV.2=ppv2, NPV.2=npv2))
	})
	
	ret <- cbind(ret, t(foo))


  ret$star <- pval2star(ret$q.value, thresholds=c(0.05, 0.01, 0.001))
  ret
}


rocseriestable.fmt <- function(tab) {
	tab$star <- NULL
	tab$AUC <- sprintf("%.2f (%.2f-%.2f)", tab$AUC, tab$CI.L, tab$CI.R)
	tab$PPV <- sprintf("%.2f (%.2f-%.2f)", tab$PPV, tab$PPV.CI.L, tab$PPV.CI.R)
	tab$NPV <- sprintf("%.2f (%.2f-%.2f)", tab$NPV, tab$NPV.CI.L, tab$NPV.CI.R)

	tab <- tab[ , !colnames(tab) %in% c("CI.L", "CI.R", "PPV.CI.L", "PPV.CI.R", "NPV.CI.L", "NPV.CI.R") ]
	tab





}



## given a factor, sample the overrepresented groups such that each factor
## level is represented the same number of times
balanceset <- function( x ) {

  if( ! is.factor( x ) ) x <- factor( x )

  counts <- summary( x )
  min.c  <- min( counts )
  min.c.n <- names( counts )[ which.min( counts ) ]

  sel <- which( x == min.c.n )
  num.min <- length( sel )
  # catf( "\n%s %d %d\n", min.c.n, min.c, num.min )

  for( l in levels( x )[ levels( x ) != min.c.n ] ) {
    sel <- c( sel, sample( which( x == l ), num.min ) )
  }

  return( sel )
}


## kfold cross-validation with additional features
mykfold <- function(x, response, grouping_var=NULL, folds=10) {


  require(randomForest)
  require(caret)
  if(is.null(grouping_var)) grouping_var <- 1:length(response)

  gvuq <- unique(grouping_var)
  kf <- createFolds(gvuq, k=folds)
  kf <- sapply(kf, function(f) which(grouping_var %in% gvuq[f]), simplify=F)

  ret <- sapply(kf, function(f) {
    cat(".")
    sel <- setdiff(1:length(response), f)
    rf.m <- randomForest(t(x[,sel]), factor(response[sel]))
    predict( rf.m, newdata= t(x[,f,drop=F]), type= 'prob' ) 
  }, simplify=F)
  cat("\n")

  ret <- Reduce(rbind, ret)


  ret <- data.frame( ret )
  ret$decision <- colnames( ret )[ ( ret[,1] < ret[,2] ) + 1 ]
  ret$reality  <- response[Reduce(c, kf)]
  ret$prediction <- ret[,1]
  return( ret )
}

## kfold cross-validation with additional features
mykfold2 <- function(x, response, grouping_var=NULL, folds=10, var=NULL) {


  require(randomForest)
  require(caret)
  if(is.null(grouping_var)) grouping_var <- 1:length(response)

  gvuq <- unique(grouping_var)
  kf <- createFolds(gvuq, k=folds)
  kf <- sapply(kf, function(f) which(grouping_var %in% gvuq[f]), simplify=F)

  ret <- sapply(kf, function(f) {
    cat(".")
    sel <- setdiff(1:length(response), f)
    rf.m <- randomForest(t(x[,sel]), factor(response[sel]), importance=TRUE)

		if(!is.null(var)) {
			selvar <- order(-rf.m$importance[,4])[1:var]
    	rf.m <- randomForest(t(x[selvar,sel]), factor(response[sel]), importance=TRUE)
		}

    predict( rf.m, newdata= t(x[,f,drop=F]), type= 'prob' ) 
  }, simplify=F)
  cat("\n")

  ret <- Reduce(rbind, ret)


  ret <- data.frame( ret )
  ret$decision <- colnames( ret )[ ( ret[,1] < ret[,2] ) + 1 ]
  ret$reality  <- response[Reduce(c, kf)]
  ret$prediction <- ret[,1]
  return( ret )
}




## leave one out cross-validation with additional functions
myloo2 <- function(x, response, grouping_var=NULL, 
                      secondary_model= FALSE, nvar= 20, balancing= TRUE, 
                      mask= NULL, test_mask= NULL, bestsample_sel= FALSE) {

  require( randomForest )
  .importance <- secondary_model

  ret <- NULL
  nn  <- ncol(x)

  ## loop over all samples
  for( i in 1:nn ) {
    cat( sprintf( "\r%d %.0f%%", i, 100*i/nn ) )

    # from the training set, omit all the samples with the same grouping
    # var (e.g. all samples from the same person)
    if( !is.null( grouping_var ) ) train.sel <- grouping_var != grouping_var[i]
    else                           train.sel <- 1:nn != i
    
    if( ! is.null( mask ) ) train.sel <- train.sel & mask

    train.x <- x[,train.sel, drop=F]
    train.r <- response[train.sel]

    if(balancing) {
      train.sel <- balanceset(train.r)
      train.x <- train.x[,train.sel, drop=F]
      train.r <- train.r[train.sel]
    }

    #if( balancing ) train <- balanceset( tmp$targets[ , response ] )
    #else            train <- 1:ncol( tmp )

    rf.m <- randomForest( t(train.x), factor(train.r), importance= .importance, mtry= nrow(train.x) )

    if(secondary_model) {
      if( nvar > nrow( rf.m$importance ) ) nvar <- nrow( rf.m$importance )
      sel_var <- order( rf.m$importance[,4], decreasing= TRUE )[1:nvar]

      if(bestsample_sel) {
        s_sel <- rf.m$predicted == train.r 
        train.r <- train.r[s_sel]
        train.x <- train.x[s_sel]
      }

      rf.m <- randomForest(t(train.x[sel_var,]), factor(train.r), mtry= nvar )
    }

    ret <- rbind( ret, predict( rf.m, newdata= t(x[,i,drop=F]), type= 'prob' ) )
  }

  cat( "\n" )

  ret <- data.frame( ret )
  ret$decision <- colnames( ret )[ ( ret[,1] < ret[,2] ) + 1 ]
  ret$reality  <- response
  ret$prediction <- ret[,1]
  return( ret )
}


## leave one out cross-validation with additional functions
myloo <- function( rg, match_sample_type= FALSE, grouping_var= NULL, 
                       secondary_model= FALSE, nvar= 20, balancing= TRUE, 
                       mask= NULL, test_mask= NULL, bestsample_sel= FALSE,
                       response= "GROUP" ) {
  require( randomForest )


  .importance <- secondary_model

  ret <- NULL
  nn  <- ncol( rg )

  for( i in 1:ncol( rg ) ) {
    cat( sprintf( "\r%d %.0f%%", i, 100*i/ncol(rg) ) )
    s <- rg[ , i ]

    # from the training set, omit all the samples with the same grouping
    # var (e.g. all samples from the same person)
    if( !is.null( grouping_var ) ) tmp.sel <- grouping_var != grouping_var[i]
    else                           tmp.sel <- rep( TRUE, nn )
    
    if( ! is.null( mask ) ) tmp.sel <- tmp.sel & mask

    if( match_sample_type )
      tmp.sel <- tmp.sel & rg$targets$SAMPLE_TYPE == s$targets$SAMPLE_TYPE

    tmp <- rg[ , tmp.sel ]

    if( balancing ) train <- balanceset( tmp$targets[ , response ] )
    else            train <- 1:ncol( tmp )

    tmp <- tmp[ , train ]

    rf.m <- randomForest( t( tmp$E ), factor( tmp$targets[ , response ] ), importance= .importance, mtry= nrow( tmp$E ) )

    if( secondary_model ) {

      if( nvar > nrow( rf.m$importance ) ) nvar <- nrow( rf.m$importance )
      sel_var <- order( rf.m$importance[,4], decreasing= TRUE )[1:nvar]

      if( bestsample_sel ) {
        s_sel <- rf.m$predicted == tmp$targets[ , response ] 
        tmp <- tmp[ , s_sel ]
        tmp <- tmp[ , balanceset( tmp$targets[ , response ] ) ]

      }


      rf.m <- randomForest( t( tmp$E[ sel_var, ] ), factor( tmp$targets[ , response ] ), mtry= nvar )
    }

    ret <- rbind( ret, predict( rf.m, newdata= t( s$E ), type= 'prob' ) )

  }

  cat( "\n" )

  ret <- data.frame( ret )
  ret$decision <- colnames( ret )[ ( ret[,1] < ret[,2] ) + 1 ]
  ret$reality  <- rg$targets[ , response ]
  ret$prediction <- ret[,1]
  return( ret )
}


summarize_rets <- function( l, meta ) {

  ll <- length( l )

  ret <- data.frame( matrix( 0, nrow= length( unique( meta$DONOR_ID ) ), ncol= 2 ) )
  rownames( ret ) <- unique( meta$DONOR_ID )
  colnames( ret ) <- c( "case", "control" )
  print( ret )

  for( i in 1:ll ) {
    r <- l[[i]]

    map <- as.character( meta[ rownames( r ), ]$DONOR_ID )
    case <- map[ r$decision == "case" ]
    control <- map[ r$decision == "control" ]

    ret[ case, ]$case <- ret[ case, ]$case + 1
    ret[ control, ]$control <- ret[ control, ]$control + 1

  }

  return( invisible( ret ) )
}


## use caret package for creating folds; make sure each fold does not
## contain samples from the same subjects as the training set
createFoldsBySubject <- function(subject, foldfun=createFolds, ...) {

  subjectID <- unique(subject)
  folds <- foldfun(subjectID, ...)

  ret <- lapply(folds, function(x) which(subject %in% subjectID[x]))
  ret
}


plotRocResample <- function(x, ...) {

  add <- FALSE
  for(i in unique(x$Resample)) {
    tmp <- x[ x$Resample == i,]
    roc(response=tmp$obs, predictor=tmp$control, plot=T, add=add, ...)
    if(!add) add <- TRUE
  }


}

## show diagnostic information about an RF model
rf.diagn <- function(rf) { 
  require(beeswarm)
  print(rf) ; rf.pr <- rf2pr(rf) ; 
  oldpar <- par(mfrow=c(1,2))
  rocplot(rf.pr)
  on.exit(par(oldpar))
  boxplot(rf.pr$case ~ rf.pr$reality)
  beeswarm(rf.pr$case ~ rf.pr$reality, add=T) 
}

## make a data frame with new predictions
makepred <- function(name, rf, rg, thr=NULL) {
  x1 <- predict(rf, newdata=t(rg$E), type="prob")[,1]

  if(is.null(thr)) {
    pr <- rf$votes[,1]
    re <- as.character(rf$y)
    tmp <- roc(response=re, predictor=pr)
    thr <- coords(tmp, "b", "t")[1]
    warning(sprintf("threshold set to %.2f", thr))
  }

  x2 <- c( "control", "case" )[ (x1 > thr) + 1 ]
  ret <- data.frame(x1, x2)
  colnames(ret) <- c(name, paste0(name, "CLASS"))
  return(ret)
}

## merge two ELists, intersecting the gene IDs
fusesets <- function(x1, x2, geneset=NULL) {
  require(limma)
  common <- intersect(rownames(x1), rownames(x2))
  if(!is.null(geneset))
    common <- intersect(common, geneset)
  if(length(common) == 0) stop("No common genes")
  commont <- intersect(colnames(x1$targets), colnames(x2$targets))
  if(length(commont) == 0) stop("No common target columns")
  x1 <- x1[match(common, rownames(x1)),]
  x2 <- x2[match(common, rownames(x2)),]
  new("EList",
    list(E=cbind(x1$E, x2$E),
         genes=x1$genes,
         targets=rbind(x1$targets[,commont], x2$targets[,commont])
         ))
}

## calculate the p-value in fishers exact test
fisher.exact <- function(x) {
  n <- length(x)
  1-pchisq(-2 * sum(log(x)), 2 * n)
}


calcRevTP <- function(tp, donors) {

  revtp <- rep(NA, length(tp))
  for(d in unique(donors)) {
    sel <- donors == d
    x <- factor(tp[sel], levels=c("BL", "M06", "M18"))
    x <- as.numeric(droplevels(x))
    revtp[sel] <- gsub("minus0", "last", paste0("minus", max(x) - x))
  }
  revtp
}


fisherpred <- function(x, donors) {
  donors <- as.character(donors)
  sapply(unique(donors), function(d) {
    fisher.exact(x[donors==d])
  })
}

mkfisherpred <- function(name, x, donors, thr) {

  x1 <- fisherpred(x, donors)
  x2 <- c( "control", "case" )[ (x1 > thr) + 1 ]
  ret <- data.frame(x1, x2)
  colnames(ret) <- c(name, paste0(name, "CLASS"))
  ret
}

mkfishercons <- function(name, df, thr) {
  x1 <- apply(df, 1, function(x) {
    fisher.exact(x)})
  x2 <- c( "control", "case" )[ (x1 > thr) + 1 ]
  ret <- data.frame(x1, x2)
  colnames(ret) <- c(name, paste0(name, "CLASS"))
  return(ret)
}


## make a series of predictions based on different models
testTrainSets <- function(train, test, thr=NULL) {

  ret <- lapply(names(train), function(xn) {
    print(xn)
    x <- train[[xn]]

    common <- intersect(rownames(test$genes), rownames(x$genes))
    rf <- randomForest(t(x$E[common,]), factor(x$targets$group))
    rf.diagn(rf)

    makepred(xn, rf, test, thr)
  })

  Reduce(cbind, ret)
}
 



doubleplot <- function( m, myylim= NULL ) {

  heal <- rgm$targets$group == "A_TST_NEG" 
  lat  <- rgm$targets$group == "B_TST_POS" 
  tb   <- rgm$targets$group == "C_TB_ACTIVE"
  case <- rg2$targets$CASE_CONTROL == "case"

  hh <- strheight( "Qp" ) * 1.8
  layout( mat= matrix( c( 1, 2 ), nrow= 1 ), widths= c( 3/5, 2/5 ) )
  if( is.null( myylim ) ) 
    myylim <- range( c( rgm$E[ m, ], rg2$E[ m, ] ) ) ; myylim[2] <- myylim[2] + 2 * hh

  boxplot( rgm$E[ m, ] ~  rgm$targets$group, names= c( "healthy", "infected", "TB" ), frame= F, ylim= myylim )
  beeswarm( rgm$E[ m, ] ~  rgm$targets$group, pwcol= c( "#33333333", "#33333333", "#66000066" )[ factor( rgm$targets$group ) ], add= T, pch= 19 )
  p1 <- wilcox.test( rgm$E[ m, heal ], rgm$E[ m, tb ] )$p.value
  if( p1 < 0.05 ) signbox( myylim[2], 1, 3, p1 )
  p1 <- wilcox.test( rgm$E[ m, lat ], rgm$E[ m, tb ] )$p.value
  if( p1 < 0.05 ) signbox( myylim[2] - hh, 2, 3, p1 )

  boxplot( rg2$E[ m, ] ~  factor( rg2$targets$CASE_CONTROL, levels= c( "control", "case" ) ), frame= F, ylim= myylim )
  beeswarm( rg2$E[ m, ] ~  factor( rg2$targets$CASE_CONTROL, levels= c( "control", "case" ) ), pwcol= c( "#66000066", "#33333333" )[ factor( rg2$targets$CASE_CONTROL ) ], add= T, pch= 19 )
  p1 <- wilcox.test( rg2$E[ m, case ], rg2$E[ m, ! case ] )$p.value
  if( p1 < 0.05 ) signbox( myylim[2], 1, 2, p1 )

  dev.copy2pdf( file= gsub( "[^a-zA-Z0-9_.]*", "", gsub( " ", "_", paste( "example_", rgm$genes[ m, "BIOCHEMICAL.NAME" ], ".pdf", sep= "" ) ) ) )
}


signbox <- function( h1, x1, x2, p1, width= 0.75, height= NULL ) {
  if( is.null( height) ) height <- strheight( "Qp" ) * 1.5
  h1 <- rep( h1, 2 )
  mid <- ( x1 + x2 ) / 2
  segments( c( x1, mid + width / 2 ), h1, c( mid - width / 2, x2 ), h1 )
  segments( c( x1, x2 ), h1 - height / 2, c( x1, x2 ), h1 )
  rect( mid - width / 2, h1[2] - height / 2, mid + width / 2, h1[2] + height / 2 )
  text( mid, h1[2], pval2star( p1, thresholds= c( 0.05, 0.01, 0.001 ) ), cex= 1.5 )
}  



myvioplot <- function (x, group, omit.empty=TRUE, range = 1.5, h = NULL, ylim = NULL, names = NULL, 
    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
    at, add = FALSE, wex = 1, drawRect = TRUE, bty="n", 
    xlab="", ylab="", ...) {

  require(sm)

  if(length(x) != length(group)) stop("Length of group must be equal to length of x")
  sel <- !is.na(x) & !is.na(group)

  x     <- x[sel]
  group <- group[sel]

  if(!is.factor(group)) group <- factor(group)

  if(omit.empty) {
    group <- factor(group)
  }


  n <- length(levels(group))

  if(length(col) == 1) col <- rep(col, n)

  # datas <- list(x, ...)
  # n <- length(datas)

  if (missing(at)) at <- 1:n

  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)

  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med  <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)

  args <- list(display = "none")

  if (!(is.null(h))) args <- c(args, h = h)

  for (i in 1:n) {
      # data <- datas[[i]]
      data <- x[ group == levels(group)[i] ]
      if(length(data) == 0) next ;
      data.min <- min(data)
      data.max <- max(data)
      q1[i] <- quantile(data, 0.25)
      q3[i] <- quantile(data, 0.75)
      med[i] <- median(data)
      iqd <- q3[i] - q1[i]
      upper[i] <- min(q3[i] + range * iqd, data.max)
      lower[i] <- max(q1[i] - range * iqd, data.min)
      est.xlim <- c(min(lower[i], data.min), max(upper[i], data.max))
      smout <- do.call("sm.density", c(list(data, xlim = est.xlim), args))
      hscale <- 0.4/max(smout$estimate) * wex
      base[[i]] <- smout$eval.points
      height[[i]] <- smout$estimate * hscale
      t <- range(base[[i]])
      baserange[1] <- min(baserange[1], t[1])
      baserange[2] <- max(baserange[2], t[2])
  }

  if (!add) {
      xlim <- if (n == 1) 
          at + c(-0.5, 0.5)
      else range(at) + min(diff(at))/2 * c(-1, 1)
      if (is.null(ylim)) {
          ylim <- baserange
      }
  }

  if (is.null(names)) label <- levels(group)
  else label <- names

  boxwidth <- 0.05 * wex

  if (!horizontal) {
      if (!add) {
          plot(NULL, type="n", xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, bty=bty, ...)
          #plot.window(xlim = xlim, ylim = ylim)
          axis(2)
          axis(1, at = at, label = label)
      }
      
      for (i in 1:n) {
          polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
          if (drawRect) {
              lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
              rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
                q3[i], col = rectCol)
              points(at[i], med[i], pch = pchMed, col = colMed)
          }
      }
  } else {
      if (!add) {
          plot(NULL, type="n", xlim=ylim, ylim=xlim, xaxt="n", yaxt="n", xlab=ylab, ylab=xlab, bty=bty, ...)
          #plot.window(xlim = ylim, ylim = xlim)
          axis(1)
          axis(2, at = at, label = label)
      }
  
      for (i in 1:n) {
          polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
          if (drawRect) {
              lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
                lty = lty)
              rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
                boxwidth/2, col = rectCol)
              points(med[i], at[i], pch = pchMed, col = colMed)
          }
      }
  }

  invisible(list(upper = upper, lower = lower, median = med, q1 = q1, q3 = q3))
}




hockeyL <- function(bp, x, y) {

  N   <- length(x)
  pm <- pmax(x, bp)
  fit <- lm(y ~ pm )
  coef <- coefficients(fit)
  Mu <- coef[1] + coef[2] * pm

  SSE <- sum(fit$residuals^2)
  sigma2 <- SSE/N
  L <- sum(dnorm(y, mean = Mu, sd = sqrt(sigma2), log=TRUE))
  return(L)
}




## hockey fit
## Fit a "hockey stick" linear regression, with hockey slope=0 before the breakpoint
## based on "piecewise.linear" from the SiZer package
hockey <- function (x, y, CI = FALSE, R = 1000, p = 0.05) {

  bp <- optimize(hockeyL, range(x), x=x, y=y, maximum=T)$maximum
  bp <- hockey0(x, y)
  model <- lm(y ~ pmax(x,bp))

  res <- list()
  res$breakpoint <- bp

  res$model <- model
  res$data <- data.frame(x=x, y=y)
  res$pred.x <- seq(min(x), max(x), length = 200)
  res$pred.y <- predict(model, data.frame(x = res$pred.x))

  res$CI <- NULL

  if (CI) {

    my.cp <- function(df, sel) {
      x <- df[sel, 1]
      y <- df[sel, 2]

      bp <- optimize(hockeyL, range(x), x=x, y=y, maximum=T)$maximum

      model <- lm(y ~ pmax(x,bp) )
      res <- c(bp, model$coefficients[2] )
      return(res)
    }

    boot.result <- boot(data.frame(x=x, y=y), my.cp, R = R)
    res$CI <- t(apply(boot.result$t, 2, quantile, probs = c(p/2, 1 - p/2)))
    rownames(res$CI) <- c("Breakpoint", "Slope")
    #res$CI <- t(res$CI)
  }

  class(res) <- "hockey"
  return(res)
}




plot.hockey <- function(x, ...) {

  plot(NULL, xlim=range(x$data$x), ylim=range(x$data$y), ...)
  if(!is.null(x$CI)) {
    ci.l <- x$CI[ "Breakpoint", 1]
    ci.r <- x$CI[ "Breakpoint", 2]
    rect(ci.l, min(x$data$y), ci.r, max(x$data$y), col="#cccccc", border=NA)
  }
  points(x$data$x, x$data$y)
  n <- length(x$pred.x)
  abline(v=x$breakpoint, lty=2)
  lines(c(x$pred.x[1], x$breakpoint, x$pred.x[n]), 
        c(x$pred.y[1], x$pred.y[1],  x$pred.y[n]))


}


## capitalize first letter in a string
capit <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}



## Manuscript helper: register figures

fig_register <- function(x, prefix=NULL) {

  env <- parent.env(environment())

  i <- 1
  x <- paste0(prefix, x)
  for(i in 1:length(x)) {
    env[[x[i]]] <- i
  }

}



## make lme4 model

# run.lmer <- function(rgC, model1, model2) {
# 	message("run.lmer here")
# 
# 	message("running lmer.0")
# 	lmer.0 <- apply(rgC$E, 1, function(y) lmer(as.formula(model1), data=rgC$targets, REML=F))
# 	message("running lmer.M")
# 	lmer.M <- apply(rgC$E, 1, function(y) lmer(as.formula(model2), data=rgC$targets, REML=F))
# 	message("done")
# 	lmer.res <- t(sapply(names(lmer.0), function(n) { xx <- anova(lmer.0[[n]], lmer.M[[n]]) ; Reduce(rbind, xx) }))
# 	colnames(lmer.res) <- c("df.1", "AIC.1", "BIC.1", "LL.1", "DEV.1", "x", "x1", "x2", "df.2", "AIC.2", "BIC.2", "LL.2", "DEV.2", "Chisq", "Chi_df", "p")
# 	lmer.res <- lmer.res[ , -c(6:8) ]
# 	lmer.res <- data.frame(rgC$genes[ , c("ID", "BIOCHEMICAL", "HMDB")], lmer.res[,c("AIC.1", "AIC.2", "p")])
# 	lmer.res$p.adj <- p.adjust(lmer.res$p, method="fdr")
# 
# 	lmer.res
# }

run.lmer.2 <- function(rgC, mod1, mod2) {

	require(piecewiseSEM)
	lmer.res <- apply(rgC$E, 1, function(y) {
		m1 <- lmer(as.formula(mod1), data=rgC$targets, REML=F)
		m2 <- lmer(as.formula(mod2), data=rgC$targets, REML=F)
		foo <- sem.model.fits(list(m1, m2))$Conditional
		f2 <- (foo[2] - foo[1])/(1 - foo[2]) # Cohen's f^2 measure
		res <- Reduce(rbind, anova(m1, m2))
		c(res, f2)

	})
	lmer.res <- t(lmer.res)

	colnames(lmer.res) <- c("df.1", "AIC.1", "BIC.1", "LL.1", "DEV.1", "x", "x1", "x2", "df.2", "AIC.2", "BIC.2", "LL.2", "DEV.2", "Chisq", "Chi_df", "P.Value", "F2")
	lmer.res <- lmer.res[ , -c(6:8) ]
	lmer.res <- data.frame(rgC$genes[ , c("ID", "BIOCHEMICAL", "HMDB")], lmer.res[,c("BIC.1", "BIC.2", "F2", "P.Value")])
	lmer.res$adj.P.Val <- p.adjust(lmer.res$P.Value, method="fdr")

	lmer.res
}


## runs a binomial regression model
## in the formula, the metabolite levels for a given metabolite are coded as "x"
## e.g. run.glm(rg, f="group ~ x")
## rg is an EList object
run.glm <- function(rg, mod1, mod2) {

	res <- apply(rg$E, 1, function(x) {
		m  <- glm(as.formula(mod1), data=rg$targets, family=binomial(logit))
		m2 <- glm(as.formula(mod2), data=rg$targets, family=binomial(logit))

		coef <- summary(m)$coefficients["x", ]
		p <- 1 - pchisq(m2$deviance - m$deviance, m2$df.residual - m$df.residual)

		or   <- exp(coef["Estimate"])
		names(or) <- NULL
		ci.l <- or / exp(qnorm(0.975) * coef["Std. Error"])
		ci.r <- or * exp(qnorm(0.975) * coef["Std. Error"])

		res <- c(or, ci.l, ci.r, p)
		names(res) <- c("OR", "CI.L", "CI.R", "P.Value")
		res
	})

	res <- data.frame(t(res), stringsAsFactors=FALSE)
	res$adj.P.Val <- p.adjust(res$P.Value, method="fdr")
	res
}

## format a p value into pandoc (with exponent)
pvf <- function(p, digits=2) {
	ret <- as.character(signif(p, digits=digits))
	ret <- gsub("e([-]?[0-9]+)", "⋅10^\\1^", ret)
	ret
}


myheatmap <- function(x, coords=NULL, 
  labCol=NULL, labRow=NULL,
  col=heat.colors(12), 
	cex.row=NULL, cex.col=NULL, breaks=NULL,
  showValues=FALSE,
  colLabVert=FALSE) {

  if(is.null(coords)) coords <- par("usr")

  nc <- ncol(x)
  nr <- nrow(x)

  x0 <- coords[1]
  y0 <- coords[3]
  tot.w <- coords[2] - coords[1]
  tot.h <- coords[4] - coords[3]

  pos <- cbind(as.vector(x), as.vector(col(x)), nr - as.vector(row(x)) + 1)

  if(is.null(breaks)) 
    breaks <- seq(min(pos[,1]), max(pos[,1]), length.out=length(col) + 1)

  col.xy <- col[findInterval(pos[,1], breaks, rightmost.closed=TRUE)]

  xp0 <- x0 + (pos[,2] - 1)/nc * tot.w
  yp0 <- y0 + (pos[,3] - 1)/nr * tot.h
  xp1 <- x0 + pos[,2]/nc * tot.w
  yp1 <- y0 + pos[,3]/nr * tot.h

  rect(xp0,
       yp0,
       xp1,
       yp1,
       border=NA,
       col=col.xy)

  if(showValues) 
    text((xp0+xp1)/2, (yp0+yp1)/2, sprintf("%.2f", (pos[,1])))

  sep.w <- strwidth("x")/2
  sep.h <- strheight("Xj")/2
  print(sep.h)

  if(!is.null(labCol)) {
    xt <- x0 + seq(0.5, nc - 0.5)/nc * tot.w
    yt <- rep(y0 - sep.h, nc)
    if(colLabVert) {
      text(xt, yt, labCol, srt=90, pos=2)
    } else {
      text(xt, yt - sep.h, labCol, cex=cex.col)
    }
  }

  if(!is.null(labRow)) {
    xt <- x0 + tot.w + sep.w
    yt <- y0 + seq(nr - 0.5, 0.5, by=-1)/nr * tot.h 
    text(xt, yt, labRow, pos=4, cex=cex.row)
  }


}
