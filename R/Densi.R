
library(phangorn)
library(PhyloOrchard)
library(datelife2) # getAges

#
# mrp.supertree slightly modified from Liam Revell's phytools package  
# densiTree based on overlayPlot
#
mrp.supertree <- function (phy, weights = NULL) 
{
  if (!class(phy) == "multiPhylo") 
    stop("phy must be object of class 'multiPhylo.'")
  X <- list()
  characters <- 0
  for (i in 1:length(phy)) {
    temp <- prop.part(phy[[i]])
    X[[i]] <- matrix(0, nrow = length(phy[[i]]$tip), ncol = length(temp) - 1)
    for (j in 1:ncol(X[[i]])) X[[i]][c(temp[[j + 1]]), j] <- 1
    rownames(X[[i]]) <- attr(temp, "labels")
    if (i == 1) species <- phy[[i]]$tip.label
    else species <- union(species, phy[[i]]$tip.label)
    characters <- characters + ncol(X[[i]])
  }
  XX <- matrix(data = "?", nrow = length(species), ncol = characters, dimnames = list(species))
  j <- 1
  for (i in 1:length(X)) {
    XX[rownames(X[[i]]), c(j:((j - 1) + ncol(X[[i]])))] <- X[[i]][1:nrow(X[[i]]), 1:ncol(X[[i]])]
    j <- j + ncol(X[[i]])
  }
  contrast <- matrix(data = c(1, 0, 0, 1, 1, 1), 3, 2, dimnames = list(c("0", "1", "?"), c("0", "1")), byrow = TRUE)
  XX <- phyDat(XX, type = "USER", contrast = contrast)  
#  print(min(fitch(phy, XX)))
  supertree <- pratchet(XX, trace = 0)
#  print(min(fitch(supertree, XX)))
  return(supertree)
}

# we want a rooted supertree
rootedSuperTree = function(x){
  fun = function(x){
    x=reorder(x, "postorder")
    nTips = length(x$tip)
    x$edge[x$edge>nTips] = x$edge[x$edge>nTips] + 2L
    l=nrow(x$edge)
    oldroot = x$edge[l,1L]
    x$edge=rbind(x$edge,matrix(c(rep(nTips+2,2),oldroot,nTips+1),2L,2L))
    x$edge.length=c(x$edge.length, 100, 100)
    x$tip.label=c(x$tip.label, "ZZZ")
    x$Nnode=x$Nnode+1L
    x
  }
  if(!is.null(attr(x, "TipLabel")))x = .uncompressTipLabel(x)
  x = unclass(x)
  x = lapply(x, fun)    
  class(x)="multiPhylo"
  res = mrp.supertree(x)
  res = root(res, "ZZZ")
  res$edge.length = rep(.1, nrow(res$edge))
  res = drop.tip(res, "ZZZ")
  reorder(res, "postorder")
  res
}


densiTree <- function(x, type="cladogram", alpha=1/length(x), optim=FALSE, scaleX=FALSE, col=1, width=1, cex=.8, ...) {
  if(class(x)!="multiPhylo")stop("x must be of class multiPhylo")
  compressed <- ifelse(is.null(attr(x, "TipLabel")), FALSE, TRUE)
  if(optim | !compressed)consensus <- rootedSuperTree(x)
  else consensus=x[[1]]
  consensus = reorder(consensus, "postorder")
  e2 = reorder(consensus)$edge[,2]
  nTip = as.integer(length(consensus$tip))
  tiporder = e2[e2<=nTip]   
  maxBT = max(getAges(x))
  if(scaleX) maxBT=1.0
  label = rev(pretty(c(maxBT,0)))
  maxBT = max(label)
  xy = plotPhyloCoor(consensus, ...)
  yy = xy[,2]
  plot.new() 
  tl = which.max(nchar(consensus$tip.label))
  sw <- strwidth(consensus$tip.label[tl],cex=cex) * 1.1
  plot.window(xlim=c(0, 1.0+sw), ylim=c(0, nTip+1))
  axis(side=1,at=seq(0,1.0, length.out=length(label)), labels=label)
  text(x=rep(1.0,Ntip(consensus)),y=yy[1:nTip],labels=consensus$tip.label,pos=4,cex=cex)  
  tip.order = yy[1:nTip]
  for (treeindex in 1:length(x)) {
    tmp <- reorder(x[[treeindex]], "postorder")
    xy <- plotPhyloCoor(tmp, tip.order=tiporder, ...)
    xx = xy[,1]
    yy = xy[,2]
    if(scaleX) xx <- xx/max(xx)
    else xx <- xx/maxBT 
    xx <- xx + (1.0 - max(xx))
    e1=tmp$edge[,1]
    e2=tmp$edge[,2]
    if(type=="cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color=adjustcolor(col, alpha=alpha), edge.width=width, edge.lty=1)
    if(type=="phylogram"){
      Ntip <- min(e1)-1L 
      Nnode <- tmp$Nnode 
      phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, TRUE, edge.color=adjustcolor(col, alpha=alpha), edge.width=width, 1) 
    }
  }  
}



data(Laurasiatherian)
set.seed(1)
bs <- bootstrap.phyDat(Laurasiatherian, FUN = function(x)upgma(dist.hamming(x)), bs=100, multicore=FALSE)
class(bs) <- 'multiPhylo'
bs = .compressTipLabel(bs)
#cnet = consensusNet(bs, .3)
#plot(cnet, show.edge.label=TRUE)

# cladogram nice to show conflict
densiTree(bs, optim=TRUE, type="cladogram", col="blue")

# cladogram are nice to show different ages
data(BinindaEmondsEtAl2007)
BinindaEmondsEtAl2007 <- .compressTipLabel(BinindaEmondsEtAl2007) 
densiTree(BinindaEmondsEtAl2007, type="phylogram", col="red")



