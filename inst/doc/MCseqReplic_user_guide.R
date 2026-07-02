## ----include=FALSE------------------------------------------------------------
library(MCseqReplic)

## ----data---------------------------------------------------------------------
exdata <- read.table(text="
                a a b b
                a a b b
                b b a a
                a c c b
                b b a c
                b b a c
                ")
weights=rep(1, nrow(exdata))
s.exdata <- TraMineR::seqdef(exdata, weights = weights, id=paste("id",1:nrow(exdata), sep=""))


## ----repdata------------------------------------------------------------------
## 3 altered sequence datasets
set.seed(3)
(altseq.list <- MCseqReplicate(s.exdata, J=1, R=3, include.obs=TRUE))


## ----repdiss------------------------------------------------------------------
(dist.list <- MCdisslist(altseq.list, method="LCS"))


## -----------------------------------------------------------------------------
MCpj(Emean=1.05, pzero=.5)

## ----MCseqdistSE--------------------------------------------------------------

(MCdistSE <- MCseqdistSE(dist.list))

## ----MCratios-----------------------------------------------------------------
MCratios(MCdistSE)

## ----MCdisscorr---------------------------------------------------------------
MCdisscorr(dist.list)


## ----repcompdiss-dissassoc----------------------------------------------------
sex <- c("f","f","f","m","m","m")
assoc.list <- lapply(dist.list, 
                     TraMineR::dissassoc, 
                     group=sex)
assoc.list[[1]]


## ----repcompdiss-dissCompare, message=FALSE-----------------------------------
library(TraMineRextras) ## for function dissCompare
comp.list <- suppressMessages(lapply(dist.list, 
                                     TraMineRextras::dissCompare, 
                                     group=sex, 
                                     squared=FALSE, 
                                     s=0))
comp.list[[1]]

## ----compgrp-stat-table-------------------------------------------------------
comptab <- MCcompgrp(dist.list, group=sex, 
                     dissassoc.args=list(R=1000),
                     dissCompare.args=list(squared=FALSE))
round(comptab,3)


## ----compgrp-summary----------------------------------------------------------
summary(comptab[-nrow(comptab),])


## ----clusqual-----------------------------------------------------------------
clqual <- MCclustqual(dist.list, clustmeth="ward.D", ncluster=4, verbose=FALSE)
round(clqual$qual.tab[[2]],3)


## ----qual.max-----------------------------------------------------------------
clqual


## ----qual.mfreq---------------------------------------------------------------
clqual$max.freq


## ----plotCQI, out.width="70%", fig.width=6, fig.height=4----------------------
ggplotMCcqi(clqual, cqi="PBC") 

## ----clusters3----------------------------------------------------------------
clust.list <- lapply(dist.list, WeightedCluster::wcKMedoids, k=3, cluster.only=TRUE)
clust.list


## ----clustcomp, echo=TRUE, warning=FALSE, message=FALSE-----------------------
(res <- MCclustcomp(clust.list, AMI=TRUE))


## ----MCmdscorr, message=FALSE-------------------------------------------------
MCmdscorr(dist.list, verbose=FALSE, core=1)


## ----MCmds-scores, out.width="70%", fig.width=8, fig.height=7-----------------
MCmdsboth <- MCmdscorr(dist.list, what="both") 
MCmds <- MCmdsboth[[2]]
nset <- length(MCmds)
title <- paste0("MC",1:nset)
if (attr(dist.list,"obs")) title[nset] <- "Obs"

layout(matrix(c(1:4,rep(5,4)),nrow=2,byrow=TRUE), heights = c(.75,.25))
for (i in 1:length(altseq.list)) {
  seqIplot(altseq.list[[i]],sortv=MCmds[[i]], with.legend=FALSE, 
           ylab=NA, yaxis=FALSE, main=title[i])
}
seqlegend(altseq.list[[1]],ncol=3, cex=1.5 )


## ----2mdsfactors--------------------------------------------------------------

mds2.list <- lapply(dist.list, cmdscale, k=2)
mds2.list[[1]]


