#' Date a tree with initial branch lengths with treePL.
#' @inheritParams use_calibrations_bladj
#' @return A phylo object
#' @details
#' This function uses treePL as described in Smith, S. A., & O’Meara, B. C. (2012).
#' \doi{10.1093/bioinformatics/bts492}, with
#' the function `treePL.phylo`. It attempts to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the
#' range of the minimum age and maximum age and try again. And repeat.
#' If expand = 0, it uses the summarized calibrations.
#' In some cases, it returns edge lengths in relative time (with maximum tree depth = 1)
#' instead of absolute time, as given by calibrations. In this case, the function returns NA.
#' This is an issue from PATHd8.
#' @references
#' Smith, S. A., & O’Meara, B. C. (2012).
#' "treePL: divergence time estimation using penalized likelihood for large phylogenies".
#' Bioinformatics, 28(20), 2689-2690, \doi{10.1093/bioinformatics/bts492}.
#' @export
use_calibrations_treePL <- function(phy, calibrations) {
  message("... Using secondary calibrations with treePL")
  phy <- input_process(phy)
  if (!inherits(phy, "phylo")) {
    message("phy is not a phylo object")
    return(NA)
  }
  if (is.null(phy$edge.length)) {
    message("phy does not have branch lengths, consider using a dating method that does not require data, such as BLADJ or MrBayes.")
    return(NA)
  }
  if (any(phy$edge.length < 0)) {
    message("phy has negative branch lengths. treePL fixes this to run, asigning them to 0 right now.")
    phy <- tree_fix_brlen(phy)
  }
  # fix any negative branch lengths?
  # phy <- tree_fix_brlen(phy, fixing_criterion = "negative", fixing_method = 0, ultrametric = FALSE)
  # following line checks whether all calibrations are present in phy, and returns that in calibs$present_calibrations
  calibs <- match_all_calibrations(phy, calibrations)
  if (nrow(calibs$present_calibrations) < 1) {
    message("Nodes in calibrations (determined by taxon pairs) do not match any nodes in phy.")
    message("Dating analysis is not possible with this set of calibrations.")
    return(NA)
  }
  chronogram <- NA
  # matches calibrations pairs to nodes in phy and summarizes min and max ages:
  used_calibrations <- calibs$matched_calibrations[, c("NodeNames", "MaxAge", "MinAge", "taxonA", "taxonB")]
  names(used_calibrations)[grep("NodeNames", names(used_calibrations))] <- "MRCA"
  # make sure that the max age of the most inclusive calibration (the deepest calibrated node) is at least the same as the oldest internal calibration:
  # first, find the most inclusive node that has a calibration:
  # this is a vector with the number of tips included by each calibrated node:
  pp <- sapply(ape::prop.part(phy), length)[as.numeric(names(calibs$phy$calibration_distribution)) - ape::Ntip(calibs$phy)]
  # the most inclusive is the one with the most extant tips,
  # but this might not be true
  # we might have two most inclusive nodes, set up a better test to find these
  maxage_deepestnode <- sapply(calibs$phy$calibration_distribution, max)[which.max(pp)]
  maxage_allnodes <- max(sapply(calibs$phy$calibration_distribution, max))
  if (maxage_deepestnode < maxage_allnodes) {
    message("max age of deepest node calibrations is smaller than max age of internal calibrations")
    message(paste("setting max age to", maxage_allnodes, "(instead of", maxage_deepestnode, ")"))
    used_calibrations[which.max(pp), "MaxAge"] <- maxage_allnodes
  }
  chronogram <- tryCatch(treePL.phylo(calibs$phy, used_calibrations,
    base = ".tmp_treePL", rm = FALSE, opts = list(
      smooth = 100, nthreads = 1, optad = 0, opt = 1,
      cvstart = 1000, cviter = 3, cvend = 0.1, thorough = TRUE
    )
  ),
  error = function(e) NULL
  )
  if (inherits(chronogram, "phylo")) {
    problem <- NULL
    if (is.null(chronogram$edge.length) | all(is.na(chronogram$edge.length))) {
      chronogram$edge.length <- NULL
      problem <- "TreePL returned a tree with no branch lengths."
    } else {
      if (all(chronogram$edge.length == 0)) {
        chronogram$edge.length <- NULL
        problem <- "TreePL returned a tree with branch lengths equal to 0."
      }
      if (is.null(problem)) { # then fix negative branch lengths again
        # chronogram$edge.length[which(chronogram$edge.length<0)] <- 0
        chronogram <- tree_fix_brlen(
          tree = chronogram, fixing_criterion =
            "negative", fixing_method = 0, ultrametric = TRUE
        )
      }
      if (round(max(ape::node.depth.edgelength(chronogram)), digits = 3) == 1) {
        problem <- "Edge lengths seem to be relative to maximum age = 1 (and not absolute to time given by calibrations)."
      }
    }
    if (!is.null(problem)) {
      message(paste(problem, "\nThis is an issue from TreePL; returning tree with a $problem."))
    }
    chronogram$problem <- problem
    chronogram$dating_method <- "treePL"
    chronogram$calibration_distribution <- calibs$phy$calibration_distribution
    chronogram$used_calibrations <- used_calibrations
    chronogram$present_calibrations <- calibs$present_calibrations
  } else {
    message("Dating analysis with TreePL failed with this tree and set of calibrations.")
    return(list(phy = calibs$phy, used_calibrations = used_calibrations))
  }
  return(chronogram)
}
# i=1
# phy <- tax_phyloall_bold[[3]][[1]]
# phy <- make_bold_otol_tree(tax_phyloallall[[3]][[1]], chronogram = FALSE)
# phy <- tax_phyloall_bold2[[3]][[1]] # is NA
# calibrations <- tax_othercalall[[3]][[1]]
# phy$edge.length
# calibrations$MaxAge
# calibs <- match_all_calibrations(tax_phyloall_bold[[1]][[i]], tax_othercalall[[1]][[i]])
# calibs$phy$edge.length
# used_calibrations <- calibs$matched_calibrations
# used_calibrations$MaxAge
# chronogram <- geiger::PATHd8.phylo(calibs$phy, used_calibrations)
# chronogram$edge.length
# plot(chronogram, main = i)
# ape::axisPhylo()
# chr <- use_calibrations_treePL(phy = catsanddogs_phyloall[[1]], calibrations = catsanddogs_calibrations)

#'
#'
#'
#'
#'
#'
treePL.phylo <- function(phy, calibrations = NULL, base = "", rm = TRUE, ...) {
  phy$node.label <- NULL
  if (!is.null(calibrations)) {
    infile <- geiger::write.treePL(phy = phy, calibrations = calibrations, base = base, ...)
  } else {
    infile <- paste(base, "infile", sep = ".")
    ape::write.tree(phy, infile)
  }
  smooth.file <- paste(base, "dated.tre", sep = ".")
  outfile <- paste(base, "treePL.orig.out", sep = ".")
  if (file.exists(outfile)) unlink(outfile)
  if (!system("which treePL", ignore.stdout = TRUE) == 0) stop("Install 'treePL' before proceeding.")
  # infile <- "~/Desktop/Evolution19Dating/TreePL_examples/ASNTreePLConfig2.txt"
  # outfile <- "~/Desktop/Evolution19Dating/TreePL_examples/ASNTreePLConfig2.treePL.orig.out"
  # infile <- ".tmp_treePL.infile"
  # outfile <- ".tmp_treePL.orig.out"
  system(paste("treePL ", infile, " >", outfile, sep = " "))
  # system(paste("grep \"tree\" ", outfile, ">", parsed.outfile, sep=" "))
  smoothed <- ape::read.tree(smooth.file)
  # smoothed=read.tree(".tmp_treePL.dated.tre.r8s")
  if (rm & base == "") {
    unlink(smooth.file)
    unlink(outfile)
    unlink(infile)
  }
  return(smoothed)
}

write.treePL <- function(phy, calibrations, nsites = 10000, min = 0.0001, base = "", opts = list(smooth = 100, nthreads = 8, optad = 0, opt = 1, cvstart = 1000, cviter = 3, cvend = 0.1, thorough = TRUE)) {
  # 	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
  # 	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB
  # 	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis
  # 	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia
  # 	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens

  if (file.exists(inp <- paste(base, "infile", sep = "."))) unlink(inp)
  if (file.exists(int <- paste(base, "intree", sep = "."))) unlink(inp)

  poss <- list(
    cv = "numeric",
    collapse = "boolean",
    checkconstraints = "boolean",
    cvstart = "numeric",
    cvstop = "numeric",
    cvmultstep = "numeric",
    verbose = "boolean",
    lftemp = "numeric",
    pltemp = "numeric",
    plcool = "numeric",
    lfstoptemp = "numeric",
    plstoptemp = "numeric",
    lfrtstep = "numeric",
    plrtstep = "numeric",
    thorough = "boolean",
    lfiter = "integer",
    pliter = "integer",
    cviter = "integer",
    ldfsimaniter = "integer",
    plsimaniter = "integer",
    cvsimaniter = "integer",
    calcgrad = "numeric",
    paramverbose = "boolean",
    prime = "boolean",
    opt = "boolean",
    optad = "boolean",
    optcvad = "boolean",
    moredetail = "boolean",
    moredetailad = "boolean",
    moredetailcvad = "boolean",
    randomcv = "boolean",
    ftol = "numeric",
    xtol = "numeric",
    mapspace = "boolean",
    nthreads = "integer"
  )
  if (length(opts) == 0) {
    print(poss)
    stop("No 'opts' specified")
  }

  # correct small branch lengths
  z <- phy$edge.length[which(phy$edge.length > 0)]
  if (any(z < min)) {
    scl <- min / min(z)
    phy$edge.length <- phy$edge.length * scl
  }
  ape::write.tree(phy, file = int)

  ## 	check appropriateness of constraints ##
  # 	check for 'calibrations' and 'phy' congruence
  # 	if(!is.null(phy)){
  # 		check=function(t, phy) all(t%in%phy$tip.label)
  # 		a=check(calibrations$taxonA, phy)
  # 		b=check(calibrations$taxonB, phy)

  # 		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
  # 	}

  ## 	build r8s file
  # 	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
  constraints <- constraintnames <- character(nrow(calibrations))
  for (i in 1:nrow(calibrations)) {
    cal <- calibrations[i, ]
    taxon <- cal$MRCA
    desc <- c(cal$taxonA, cal$taxonB)

    txt1 <- ifelse(!is.na(cal$MinAge), paste("min =", taxon, cal$MinAge, sep = " "), "")
    txt2 <- ifelse(!is.na(cal$MaxAge), paste("max =", taxon, cal$MaxAge, sep = " "), "")
    txt <- paste(txt1, txt2, sep = "\n")

    constraints[i] <- txt
    constraintnames[i] <- paste("mrca =", taxon, desc[1], desc[2], sep = " ")
  }
  infile <- list(
    tree = paste("treefile = ", int, sep = ""),
    ns = paste("numsites = ", nsites, sep = ""),
    names = paste(unlist(constraintnames), collapse = "\n"),
    mrca = paste(unlist(constraints), collapse = "\n"),
    out = paste("outfile = ", paste(base, "dated", "tre", sep = "."), sep = ""),
    opt = paste(names(opts), opts, sep = "=", collapse = "\n")
  )

  inp <- paste(base, "infile", sep = ".")
  writeLines(paste(infile, collapse = "\n\n"), con = inp)
  attr(inp, "dating_method") <- "treePL"
  return(inp)
}
