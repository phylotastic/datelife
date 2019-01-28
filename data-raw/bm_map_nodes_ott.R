# benchmarking map_nodes_ott
load(data-raw/opentree_chronograms_oct2018.rda)
# preliminary assessment:
tot_time_list <- vector(mode = "list", length = length(opentree_chronograms$trees))
for (i in seq(opentree_chronograms$trees)){
    print(i)
    start_time <- Sys.time()
    map_nodes_ott(tree = opentree_chronograms$trees[[i]])
    tot_time_list[[i]] <- Sys.time() - start_time
}
save(tot_time_list, file = "data-raw/benchmark/opentree_chronograms/bm_map_nodes_ott-tot_time_list.rda")
# plot preliminary assessment:
load(file = "data-raw/bm_map_nodes_ott-tot_time_list.rda")
ii <- sapply(seq(opentree_chronograms$trees), function(x)
        any(grepl("not.mapped", opentree_chronograms$trees[[x]]$tip.label)))
length(opentree_chronograms$trees)
oo <- opentree_chronograms$trees[!ii]
length(oo)
tt <- tot_time_list[!ii]
table(unname(sapply(oo, ape::Ntip)))
dimnames(table(unname(sapply(oo, ape::Ntip))))
tip_number <- as.numeric(names(table(unname(sapply(oo, ape::Ntip)))))
frequency <- table(unname(sapply(oo, ape::Ntip)))
plot(x = tip_number,
        y = frequency, main = "Frequency of tip number in cached Opentree chronograms")
plot(x = log(tip_number),
        y = frequency, main = "Frequency of tip number in cached Opentree chronograms")
hist(table(unname(sapply(oo, ape::Ntip))))  # 18 chronograms have 72 tips
unique(names(oo[sapply(oo, ape::Ntip) == 72]))  # belonging to only 2 studies
table(unname(sapply(oo, ape::Ntip)))
table(tt)
plot(x = unname(sapply(oo, ape::Ntip)), y = tt,
        main = "map_nodes_ott running time by tip number",
        xlab = "Tip number", ylab = "Run time")
plot(x = log(unname(sapply(oo, ape::Ntip))), y = tt,
        main = "map_nodes_ott running time by tip number",
        xlab = "log(Tip number)", ylab = "Run time")
tt[which(tt > 20)]
oo[which(tt > 20)]
oo[which(unname(sapply(oo, ape::Ntip)) > 4000)]
oo[unname(sapply(oo, ape::Ntip)) > 1000]
# Preliminary assessment summary:
# Three chronograms with ~4k tips took the most time (around 50s each to finish running map_nodes_ott)
# this is interesting bc chronograms with waymore tips (~80k) did not take as much.
# I decided to microbenchmark with one chrrnogram per tip number cohort.
# For most cohorts there's only one chronogram.
# The following cohorts have more than one chronogram. In these cases only the first one has chosen for microbenchmarking.
# Tip number cohorts with more than one chronogram:
onesies <- oo[!duplicated(unname(sapply(oo, ape::Ntip)))]  # leave one chronogram for each tip number cohort
onesies <- onesies[order(sapply(onesies, ape::Ntip))]
length(onesies) # 138
library(microbenchmark)
mbm <- c()
for(i in seq(onesies)){
    print(i)
    mm <- microbenchmark(map_nodes_ott(tree = onesies[[i]]), times = 100L)
    levels(mm$expr)[1] <- paste(ape::Ntip(onesies[[i]]), "tips")
    mbm <- rbind(mbm, mm)
}
save(mbm, file = "data-raw/bm_map_nodes_ott-mbm.rda")
x <- runif(100)
benchmark <- microbenchmark(
  sqrt(x),
  x ^ 0.5
)
names(benchmark)
library(microbenchmark)
library(ggplot2)

#' Plot a microbenchmark object with ggplot, with time in y axis
#' @param file
#' @examples
# plot_runtime(file = "data-raw/benchmark/opentree_chronograms/bm_map_nodes_ott-mbm.pdf", benchmark = mbm,
xlabel = "Number of tips in chronogram", ylabel = "Time (milliseconds)",
xticks_angle = 90, xticks_breaks = , xticks_labels = ,
yticks_breaks = c(1e+03, 1e+04, 3e+04, 6e+04), yticks_labels = )
levels(mbm$expr)
class(mbm$time)
attributes(mbm)
is.numeric(NULL)
microbenchmark:::autoplot.microbenchmark(mbm)
plot_runtime <- function(file = NULL, width = 6, height = 3, benchmark = NULL,
                            ymin = 0, ymax = 2000,
                            xlabel = "Independent var", ylabel = "Dependent var",
                            xticks_angle = NULL, xticks_labels = NULL, xticks_breaks = NULL,
                            yticks_angle = NULL, yticks_labels = NULL, yticks_breaks = NULL){
    if (!is.null(file)) {
        pdf(file= file, width = width, height = height)
    }
    plt <- ggplot2::ggplot(benchmark, ggplot2::aes_string(x = names(benchmark)[1], y = names(benchmark)[2]))
    # to change the span of time when plotting several things, not needed for now:
    if (!"expr" %in% names(benchmark)){
        plt <- plt + ggplot2::coord_cartesian(ylim = c(ymin, ymax))
    }
    plt <- plt + ggplot2::stat_ydensity()
    plt <- plt + ggplot2::scale_x_discrete(name = xlabel, if(!is.null(xticks_breaks)){breaks = xticks_breaks},
                                            if(!is.null(xticks_labels)){labels = xticks_labels})
    if(is.numeric(xticks_angle)){
        plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle= xticks_angle, size=5), legend.position = "none")
    }
    plt <- plt + ggplot2::scale_y_log10(name = ylabel)
    if (!is.null(file)) {
        print(plt)
        dev.off()
    } else {
        return(plt)
    }
}
