relevant_age_quantiles <- function(ages, probs = c(0.5, 0, 0.025, 0.975, 1)) {
  # just utility wrapper function with different defaults
  return(stats::quantile(ages, probs))
}

numeric_vector_to_html_row <- function(x, digits = 2) {
  return(paste(paste0("<td>", round(x, digits)), "</td>", sep = "", collapse = ""))
}

patristic_matrix_sample <- function(patristic_matrix_array, uncertainty) {
  # if (dim(patristic_matrix_array)[3] == 1) {
  # 	patristic_matrix<-patristic_matrix_array[,,1]
  # 	#need order of node depths, from just the upper triangular and diagonal part of the matrix
  # 	element.order<-order(patristic_matrix[upper.tri(patristic_matrix,diag = FALSE)],decreasing = TRUE)
  # 	new.patristic_matrix<-patristic_matrix*0
  # 	cur.val<-patristic_matrix[upper.tri(patristic_matrix,diag = FALSE)][element.order[1]]
  #   new.patristic_matrix[upper.tri(new.patristic_matrix,diag = FALSE)][element.order[1]] <- cur.val + runif(1, -cur.val*uncertainty/100, cur.val*uncertainty/100)
  # 	element.order<-element.order[-1]
  #  	for (i in sequence(length(element.order))) {
  #  		cur.val<-patristic_matrix[upper.tri(patristic_matrix,diag = FALSE)][element.order[i]]
  #  		new.patristic_matrix[upper.tri(new.patristic_matrix,diag = FALSE)][element.order[i]] <- cur.val + runif(1, -cur.val*uncertainty/100, min(cur.val*uncertainty/100, min( ))
  #  	}
  #  }
  #  else {
  return(patristic_matrix <- patristic_matrix_array[, , sample.int(1, size = dim(patristic_matrix_array)[3], replace = TRUE)])
  # }
}

patristic_matrix_subset <- function(patristic_matrix, taxa, phy4 = NULL) {
  # gets a subset of the patristic_matrix. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic_matrix.new <- patristic_matrix[rownames(patristic_matrix) %in% taxa, colnames(patristic_matrix) %in% taxa]
  problem.new <- "none"
  final.size <- sum(rownames(patristic_matrix.new) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem.new <- "some of the queried taxa are not on this chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic_matrix.new <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic_matrix = patristic_matrix.new, problem = problem.new))
}

#' Return the relevant authors for a set of studies
#' @param results.index A vector from datelife_result_study_index() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each author, with names equal to author names
#' @export
datelife_authors_tabulate <- function(results.index,
                                      cache = "opentree_chronograms") {
  if ("opentree_chronograms" %in% cache) {
    utils::data("opentree_chronograms", package = "datelife")
    cache <- get("opentree_chronograms")
  }
  authors <- cache$authors[results.index]
  return(table(unlist(authors)))
}

#' Return the relevant curators for a set of studies
#' @param results.index A vector from datelife_result_study_index() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each curator, with names equal to curator names
#' @export
relevant_curators_tabulate <- function(results.index,
                                       cache = "opentree_chronograms") {
  if ("opentree_chronograms" %in% cache) {
    utils::data("opentree_chronograms", package = "datelife")
    cache <- get("opentree_chronograms")
  }
  curators <- cache$curators[results.index]
  return(table(unlist(curators)))
}

#
# ReadDistance <- function(file) {
# 	data <- readLines(file, n=-1)[-1] #read in phylip distance, perhaps from SDM
# 	data <- strsplit(data, " +")
# 	data[sapply(data, length)>0] #trim trailing
# 	final_matrix <- matrix(nrow = length(data), ncol = length(data))
# 	all.names <- rep(NA, length(data))
# 	for (data.index in sequence(length(data))) {
# 		local.line <- data[[data.index]]
# 		all.names[data.index] <- local.line[1]
# 		local.line <- as.numeric(local.line[-1])
#
# 	}
# }

#
# SDM <- function(datelife_result, weights = NULL)) {
#
# 	patristic.array <- patristic_matrix_list_to_array(datelife_result)
# 	k <- length(datelife_result)
# 	if(is.null(weights)) {
# 		weights <- rep(1/k, k)
# 	} else {
# 		weights <- weights/sum(weights) #just to make sure total weight is 1
# 	}
#
# 	#Use their appendix and stick it in solve
# 	#a*x = b
# 	#The a matrix has a vector of alpha values, then a_i for matrix 1, a_i for matrix 2..., then
#
# }
#
# #rewrite of code in ape by Andrei Popescu niteloserpopescu@gmail.com to allow for debugging
# apeSDM <- function(...) {
# 	st <- list(...)
# 	k <- length(st)/2
# 	ONEtoK <- seq_len(k)
# 	for (i in ONEtoK) st[[i]] <- as.matrix(st[[i]])
# 	ROWNAMES <- lapply(st[ONEtoK], rownames)
# 	NROWS <- sapply(ROWNAMES, length)
# 	tot <- sum(NROWS)
# 	labels <- unique(unlist(ROWNAMES))
# 	sp <- unlist(st[k + ONEtoK])
# 	astart <- numeric(tot)
# 	astart[1] <- k
# 	for (i in 2:k) astart[i] <- astart[i - 1] + NROWS[i - 1]
# 	n <- length(labels)
# 	miustart <- k + tot
# 	niustart <- miustart + n
# 	lambstart <- niustart + k - 1
# 	X <- matrix(0, n, n, dimnames = list(labels, labels))
# 	V <- w <- X
# 	tmp <- 2 * k + tot + n
# 	col <- numeric(tmp)
# 	for (i in 1:(n - 1)) {
# 		for (j in (i + 1):n) {
# 			for (p in ONEtoK) {
# 				if (is.element(labels[i], ROWNAMES[[p]]) && is.element(labels[j],
# 					ROWNAMES[[p]])) {
# 						w[i, j] <- w[j, i] <- w[i, j] + sp[p]
# 					}
# 				}
# 			}
# 		}
# 		ONEtoN <- seq_len(n)
# 		Q <- matrix(0, tmp, tmp)
# 		for (p in ONEtoK) {
# 			d_p <- st[[p]]
# 			for (l in ONEtoK) {
# 				d <- st[[l]]
# 				sum <- 0
# 				dijp <- -1
# 				if (l == p) {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							if (all(!is.na(pos))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								sum <- sum + dij * dij - sp[l] * dij *
# 								dij/w[i, j]
# 								tmp2 <- dij - sp[l] * dij/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] + tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] + tmp2
# 							}
# 						}
# 					}
# 				}
# 				else {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							posp <- match(labels[c(i, j)], ROWNAMES[[p]])
# 							if (all(!is.na(pos)) && all(!is.na(posp))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								dijp <- d_p[posp[1L], posp[2L]]
# 								sum <- sum - sp[l] * dij * dijp/w[i, j]
# 								tmp2 <- sp[l] * dijp/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] - tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] - tmp2
# 							}
# 						}
# 					}
# 				}
# 				Q[p, l] <- sum
# 			}
# 			Q[p, lambstart + 1] <- 1
# 		}
# 		r <- k
# 		for (p in ONEtoK) {
# 			dp <- st[[p]]
# 			for (i in ONEtoN) {
# 				if (is.element(labels[i], ROWNAMES[[p]])) {
# 					r <- r + 1
# 					for (l in ONEtoK) {
# 						d <- st[[l]]
# 						if (l == p) {
# 							ipos <- match(labels[i], ROWNAMES[[p]])
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								jpos <- match(labels[j], ROWNAMES[[p]])
# 								if (!is.na(jpos)) {
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] + dij - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- 1 - sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] + tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] + tmp2
# 								}
# 							}
# 						}
# 						else {
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								if (!is.element(labels[j], rownames(dp)))
# 								next
# 								pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 								if (all(!is.na(pos))) {
# 									ipos <- pos[1L]
# 									jpos <- pos[2L]
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] - tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] - tmp2
# 								}
# 							}
# 						}
# 					}
# 					if (p < k)
# 					Q[r, ] <- Q[r, ] * sp[p]
# 					Q[r, miustart + i] <- 1
# 					if (p < k)
# 					Q[r, niustart + p] <- 1
# 				}
# 			}
# 		}
# 		r <- r + 1
# 		col[r] <- k
# 		Q[r, ONEtoK] <- 1
# 		for (i in ONEtoN) {
# 			r <- r + 1
# 			for (p in ONEtoK) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		for (p in 1:(k - 1)) {
# 			r <- r + 1
# 			for (i in ONEtoN) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		a <- solve(Q, col, 1e-19)
# 		for (i in ONEtoN) {
# 			for (j in ONEtoN) {
# 				if (i == j) {
# 					X[i, j] <- V[i, j] <- 0
# 					next
# 				}
# 				sum <- 0
# 				sumv <- 0
# 				for (p in ONEtoK) {
# 					d <- st[[p]]
# 					pos <- match(labels[c(i, j)], ROWNAMES[[p]])
# 					if (all(!is.na(pos))) {
# 						ipos <- pos[1L]
# 						jpos <- pos[2L]
# 						dij <- d[ipos, jpos]
# 						sum <- sum + sp[p] * (a[p] * dij + a[astart[p] +
# 							ipos] + a[astart[p] + jpos])
# 							sumv <- sumv + sp[p] * (a[p] * dij)^2
# 						}
# 					}
# 					X[i, j] <- sum/w[i, j]
# 					V[i, j] <- sumv/(w[i, j])^2
# 				}
# 	}
# 	list(X, V)
# }
#
#
#
#
#
