## create subsets of data stratified by SITE and proximity to diagnosis
## set 1-3: the training set

cat( "make_strata.R: generating stratified subsets\n")

## stratification by site and proximity to disease

## for each site, the distal, proximate and full data sets, merged from
## plasma and rpmi
ss3 <- list()
cutoff <- 5 ## cutoff defining the proximate / distal data

## creating sample type stratification, sets 1-3 only
## "ss" stands for "sample stratification"
stype.ss3 <- sapply(c("plasma", "serum", "rpmi"), function(stype) {

	x <- get(stype)

	x <- x[, colnames(x) %in% set123$SAMPLE_ID ]
	xmeta <- set123[ colnames(x), ]
	list(all=x,
	     distal=x[ , xmeta$TIME_TO_TB >= cutoff | xmeta$group == "control" ],
			 proximate=x[ , xmeta$TIME_TO_TB < cutoff | xmeta$group == "control" ])
}, simplify=FALSE)
stype.ss3 <- unlist(stype.ss3, recursive=FALSE)

stype.ss3[["rpmi.proximate"]] <- NULL # only one case

## creating site stratification, sets 1-3 only
## "ss" stands for "sample stratification"
site.sel <- c(sapply(sites, function(s) set123$SITE == s, simplify=FALSE), list(TOT=TRUE))

ss3 <- sapply(site.sel, function(sel) {
  se <- serum[ , colnames(serum) %in% set123$SAMPLE_ID[sel] ]
  pl <- plasma[ , colnames(plasma) %in% set123$SAMPLE_ID[sel] ]
  rp <- rpmi[ , colnames(rpmi) %in% set123$SAMPLE_ID[sel] ]

  foo.all <- mergesets(se, pl, rp)
  tmp <- set123[ colnames(foo.all), ]
  list(all=foo.all, 
       proximate=foo.all[ , tmp$TIME_TO_TB < cutoff  | tmp$group == "control" ],
       distal=foo.all[ , tmp$TIME_TO_TB >= cutoff | tmp$group == "control"]
			 )
}, simplify=FALSE)

ss3 <- unlist(ss3, recursive=FALSE)

ss3$MAK.proximate <- NULL ## number of cases not sufficient





## ----------------- Creating site stratification --------------------------
## creating site stratification, sets 1-4 (full)
## "ss" stands for "sample stratification"

## for each site, the distal, proximate and full data sets, merged from
## plasma and rpmi
cutoff <- 5
sites <- unique(set$SITE)
site.sel <- c(sapply(sites, function(s) set$SITE == s, simplify=FALSE), list(TOT=TRUE))
ss4 <- list()

ss4 <- sapply(site.sel, function(sel) {
  se <- serum[ , colnames(serum) %in% set$SAMPLE_ID[sel] ]
  pl <- plasma[ , colnames(plasma) %in% set$SAMPLE_ID[sel] ]
  rp <- rpmi[ , colnames(rpmi) %in% set$SAMPLE_ID[sel] ]

  foo.all <- mergesets(se, pl, rp)
  tmp <- set[ colnames(foo.all), ]
  list(all=foo.all, 
       distal=foo.all[ , tmp$TIME_TO_TB >= cutoff | tmp$group == "control"],
       proximate=foo.all[ , tmp$TIME_TO_TB < cutoff  | tmp$group == "control" ])
}, simplify=FALSE)

ss4 <- unlist(ss4, recursive=FALSE)

cat("make_strata.R: done\n")
