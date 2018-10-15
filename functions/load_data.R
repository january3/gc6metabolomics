## Loads all data necessary to generate the manuscript
cat( "load_data.R: Loading data\n" )

## set13: the training set
## set4: the test set
## set: the full set (training + test)
set    <- read.csv("data/full_unblinded_metadata_with_smoking_tst.csv", row.names=1)
set4   <- set[ set$SET == "SET4", ]
set123 <- set[ set$SET %in% c("SET1", "SET2", "SET3"), ]

sites  <- unique(set$SITE)
donors <- set[ !duplicated(set$DONOR_ID), ]
donors$TIMEPOINT <- donors$TIME_TO_TB <- donors$revTP <- NULL

## sample data
plasma <- read.csv("data/measurements_plasma_full.csv", row.names=1, stringsAsFactors=FALSE)
serum  <- read.csv("data/measurements_serum_full.csv", row.names=1, stringsAsFactors=FALSE)
rpmi   <- read.csv("data/measurements_plasmarpmi_full.csv", row.names=1, stringsAsFactors=FALSE)
all.biochemicals <- read.csv("data/biochemicals_full_list_4.csv", row.names=1, stringsAsFactors=FALSE)

## [Omitted]
## Filter by minimal number of values higher than the imputed minimum
#plasma <- filter.by.min(plasma, n=50)
#serum  <- filter.by.min(serum, n=50)
#rpmi   <- filter.by.min(rpmi, n=25)

## Create limma-like models and rank-based normalization
## Create limma's EList objects from data frames
r.serum  <- make.EList(data.matrix(serum), all.biochemicals, set)
r.plasma <- make.EList(data.matrix(plasma), all.biochemicals, set)
r.rpmi   <- make.EList(data.matrix(rpmi), all.biochemicals, set)


## rank-normalize the Elists
r.serumR  <- EList.ranknorm(r.serum)
r.plasmaR <- EList.ranknorm(r.plasma)
r.rpmiR   <- EList.ranknorm(r.rpmi)

## Comparison with Weiner et al. 2012 data
metabo <- read.csv("data/weiner_2012.csv", row.names=1, stringsAsFactors=F)
mm <- metabo[,14:ncol(metabo)]
tt <- setNames(data.frame(Reduce(rbind, strsplit(colnames(mm), split="_")), row.names=NULL), c("ID", "CLASS", "sex"))
tt <- data.frame(tt, row.names=tt[,"ID"], group=c("control", "case")[ (tt$CLASS == "TB") + 1])
colnames(mm) <- tt$ID
r.metabo <- make.EList(mm, metabo[,1:13], tt)
rm(tt, mm, metabo)

## Additional metabolic profiling data: TB patients vs other respiratory
## diseases (ORD)
tbord <- read.csv("data/TB_ORD_Gambia_Sutherland.csv")
rownames(tbord) <- tbord$SAMPLE_NAME
tbord$group <- tbord$CASE_CONTROL

tbord.g <- read.csv("data/TB_ORD_Gambia_Sutherland_biochemicals.csv", sep="\t")[,1:14]
tbord.d <- read.csv("data/TB_ORD_Gambia_Sutherland_data.csv", sep=",", comment.char="")[,-c(1:16)]
rownames(tbord.d) <- rownames(tbord.g) <- tbord.g$ID <- paste0("M.", tbord.g$COMP_ID)


## predictions.
## predictions were generated using code shown in the file "predictions.R", only once.
## the files have been then saved, locked and disseminated in the consortium
## this code is not re-run when the manuscript is compiled!

pred          <- read.csv("data/set4_predictions_weiner.csv", row.names=1)
pred.decision <- read.csv("data/set4_predictions_weiner_classes.csv", row.names=1)
pred.nvars    <- read.csv("data/set4_predictions_weiner_nvars.csv", row.names=1)
pred$ID <- NULL

## remove Metabo and Metabo15; these models are not necessary for blinded
## predictions since they
## are wholy independent of the training set
## instead, they are applied to the full set1234 data set
pred$Metabo   <- NULL
pred$Metabo15 <- NULL

colnames(pred)[colnames(pred)=="Total"] <- "TotalF"
colnames(pred.decision)[colnames(pred.decision)=="Total"] <- "TotalF"

pred.nvars$ID[ pred.nvars$ID == "Total" ] <- "TotalF"
pred.nvars <- rbind(pred.nvars, list(ID="Compact", N=25))

## predictions from Fergal Duffy

fd_univ <- read.table("data/FD_gbm_universal_loocv_predictions.csv", sep=",", quote='"', header=T)
colnames(fd_univ) <- c("X", "prediction", "reality", "sampleID", "Site", "timeToTB", "donorID", "type", "TIMEPOINT", "MetabDonorID", "DonorID", "RnaID")
fd_univ$decision <- c("case", "control")[ (fd_univ$prediction < 0.5) + 1 ]
fd_blind <- read.table("data/FD_gbm_universal_blind_predictions.csv", sep=",", quote='"', header=T)
colnames(fd_blind) <- c("X", "prediction", "sampleID", "MetabDonorID.x", "Site", "reality", "timepoint", "sampleType", "MetabDonorID.y", "DonorID", "RnaID")
fd_blind$decision <- c("case", "control")[ (fd_blind$prediction < 0.5) + 1 ]
fd_blind$timeToTB <- set$TIME_TO_TB[ match(fd_blind$sampleID, set$SAMPLE_ID) ]

pred$Compact          <- fd_blind$prediction[ match(rownames(pred), fd_blind$sampleID) ]
pred.decision$Compact <- fd_blind$decision[ match(rownames(pred.decision), fd_blind$sampleID)]

pred.models <- colnames(pred)
cat("load_data.R: done\n")




