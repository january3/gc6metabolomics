## ----------------------------------------------------------------------
## Note: For the sake of keeping the original predictions untouched, this
## script has been left as is. Some of the variables are outdated or renamed
## in the latter versions of the script for result analysis.
## Contact january@mpiib-berlin.mpg.de if you wish to run this script
## yourself.
## ----------------------------------------------------------------------


## ----------------------------------------------------------------------
## SET4 predictions
## ----------------------------------------------------------------------

set.seed(12345)

## The Metabo models
pred4l <- list()
pred4l$metabo <- list()
x          <- metabo_export$E
x.response <- metabo_export$targets$GROUP

y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$metabo$plasma <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$metabo$serum <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$metabo$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)

## 15 best variables
pred4l$metabo15 <- list()
y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$metabo15$plasma15 <- xpredict(x, y, x.response, y.response=NULL, nvar=15)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$metabo15$serum15 <- xpredict(x, y, x.response, y.response=NULL, nvar=15)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$metabo15$plasmarpmi15 <- xpredict(x, y, x.response, y.response=NULL, nvar=15)


## serum samples, < 6M, SUN baseline excluded
## "LS"  - "late serum"
pred4l$ls <- list()
sel <- with(set123, SAMPLE_TYPE == "serum" &
                    (TIME_TO_TB < 6 | group == "control") &
                    (SITE != "SUN" | TIMEPOINT != "BL")) 
x <- serum[ , colnames(serum) %in% set123$SAMPLE_ID[sel] ]
x.response <- set123[ colnames(x), "group" ]

y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$ls$plasma <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$ls$serum <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$ls$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)

## serum samples, BL
## "BS"  - "baseline serum"
pred4l$bs <- list()
sel <- with(set123, SAMPLE_TYPE == "serum" & 
                    (TIMEPOINT == "BL" | group == "control"))
x <- serum[ , colnames(serum) %in% set123$SAMPLE_ID[sel] ]
x.response <- set123[ colnames(x), "group" ]

y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$bs$plasma <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$bs$serum <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$bs$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)


## plasma samples, < 6M, SUN baseline excluded
## "LP"  - "late plasma"
pred4l$lp <- list()
sel <- with(set123, SAMPLE_TYPE == "plasma" &
                    (TIME_TO_TB < 6 | group == "control" ) &
                    (SITE != "SUN" | TIMEPOINT != "BL")) 
x <- plasma[ , colnames(plasma) %in% set123$SAMPLE_ID[sel] ]
x.response <- set123[ colnames(x), "group" ]

y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$lp$plasma <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$lp$serum <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$lp$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)

## plasma samples, BASELINE
## "BP"  - "baseline plasma"
pred4l$bp <- list()
sel <- with(set123, SAMPLE_TYPE == "plasma" & 
                    (TIMEPOINT == "BL" | group == "control"))

x <- plasma[ , colnames(plasma) %in% set123$SAMPLE_ID[sel] ]
x.response <- set123[ colnames(x), "group" ]

y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$bp$plasma <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$bp$serum <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$bp$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL) #, nvar=10)


## TOTAL LATE
## combined predictions. all models, all samples < 6M and not sun baseline
sel <- with(set123,
                    (TIME_TO_TB < 6 | group == "control") &
                    (SITE != "SUN" | TIMEPOINT != "BL"))
common <- intersect(rownames(rg$genes), rownames(metabo_export$genes))
x <- cbind(rg$E[common,sel], metabo_export$E[common, ])
x.response <- c(set123$group[sel], metabo_export$targets$GROUP)

rf <- randomForest(t(x), factor(x.response), importance=TRUE)
rf.total <- rf
save(rf.total, file="rf.total.rda")

pred4l$tot <- list()
y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$tot$plasma <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$tot$serum <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$tot$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)


## TOTAL BASELINE
## combined predictions. all models, all samples < 6M and not sun baseline
sel <- with(set123, TIMEPOINT == "BL")
common <- intersect(rownames(rg$genes), rownames(metabo_export$genes))
x <- cbind(rg$E[common,sel], metabo_export$E[common, ])
x.response <- c(set123$group[sel], metabo_export$targets$GROUP)

rf <- randomForest(t(x), factor(x.response), importance=TRUE)
rf.total.bl <- rf
save(rf.total.bl, file="rf.total.bl.rda")

pred4l$tot.bl <- list()
y          <- plasma[, colnames(plasma) %in% set4$ID]
pred4l$tot.bl$plasma <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
pred4l$tot.bl$serum <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
pred4l$tot.bl$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)


## GG ratio
sel <- with(rg$targets, SAMPLE_TYPE == "serum" &
                    (TIME_TO_TB < 6 | group == "control") &
                    (SITE != "SUN" | TIMEPOINT != "BL"))
x <- t(data.frame(ratioGG=rg$E[ "M.53", sel]/rg$E["M.57", sel], cortisol=rg$E["M.1712", sel]))
x.response <- rg$targets$group[sel]
rf <- randomForest(t(x), factor(x.response), importance=TRUE)

pred4l$ggc <- list()
y          <- plasma[, colnames(plasma) %in% set4$ID]
y <- t(data.frame(ratioGG=y[ "M.53", ]/y["M.57", ], cortisol=y["M.1712", ]))
pred4l$ggc$plasma <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- serum[, colnames(serum) %in% set4$ID]
y <- t(data.frame(ratioGG=y[ "M.53", ]/y["M.57", ], cortisol=y["M.1712", ]))
pred4l$ggc$serum <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)
y          <- plasmarpmi[, colnames(plasmarpmi) %in% set4$ID]
y <- t(data.frame(ratioGG=y[ "M.53", ]/y["M.57", ], cortisol=y["M.1712", ]))
pred4l$ggc$plasmarpmi <- xpredict(x, y, x.response, y.response=NULL, rf=rf) #, nvar=10)


extrd <- function(d, field) {
  Reduce(c, 
    sapply(d, function(x) setNames(x[[field]], rownames(x))))
}

pred4vec <- sapply(pred4l, function(x) extrd(x, "case"))
pred4 <- data.frame(ID=rownames(pred4vec), pred4vec, stringsAsFactors=FALSE)
colnames(pred4) <- c("ID", "Metabo", "Metabo15", "Late.serum", "BL.serum", "Late.plasma", "BL.plasma", "Total", "Total.BL", "GGC")
write.csv(pred4, file="set4_predictions_weiner.csv")
write.csv(pred4, file="data/set4_predictions_weiner.csv")

pred4vec.class <- sapply(pred4l, function(x) extrd(x, "decision2"))
pred4.class <- data.frame(ID=rownames(pred4vec.class), pred4vec.class, stringsAsFactors=FALSE)
colnames(pred4.class) <- c("ID", "Metabo", "Metabo15", "Late.serum", "BL.serum", "Late.plasma", "BL.plasma", "Total", "Total.BL", "GGC")
pred4.class$consensus <- c( "control", "case" )[ 
  (apply(pred4.class[,-1], 1, function(x) sum(x == "case")) > (ncol(pred4.class)-1)/2) + 1
  ]


write.csv(pred4.class, file="set4_predictions_weiner_classes.csv")
write.csv(pred4.class, file="data/set4_predictions_weiner_classes.csv")

cutoffs <- sapply(pred4l, function(x) sapply(x, function(y) attr(y, "cutoff")))
colnames(cutoffs) <- c("Metabo", "Metabo15", "Late.serum", "BL.serum", "Late.plasma", "BL.plasma", "Total", "Total.BL", "GGC")
write.csv(cutoffs, file="set4_predictions_weiner_cutoffs.csv")
write.csv(cutoffs, file="data/set4_predictions_weiner_cutoffs.csv")


nvars <- sapply(pred4l, function(x) min(sapply(x, function(xx) attr(xx, "nvar"))))
nvars <- data.frame(ID=c("Metabo", "Metabo15", "Late.serum", "BL.serum", "Late.plasma", "BL.plasma", "Total", "Total.BL", "GGC"),
                    N=nvars)
write.csv(nvars, file="set4_predictions_weiner_nvars.csv")
write.csv(nvars, file="data/set4_predictions_weiner_nvars.csv")

