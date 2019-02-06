# test_that("summary works for a set of different taxa", {
#     skip("test taxa batch")
#     skip_on_cran()
#     skip_on_travis()
    # utils::data(names_subset2)
    # taxa_list <- list(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis"))
    # taxa_list <- c(taxa_list, lapply(1:50, function(x) sample(names_subset2, size = 20, replace = FALSE)))
    taxa_list <- list(c("Yixianornis grabaui", "Amphibia", "Amphibia", "Amphibia",
      "Amphibia", "Sauropsida", "Bucerotiformes", "Struthioniformes",
      "Forpus passerinus", "Anhimidae", "Deinonychus", "Vorona", "Vegavis",
      "Meleagris gallopavo", "Coenocorypha", "Megapodius", "Podargidae",
      "Melopsittacus undulatus", "Notornis", "Falco sparverius", "Testudines",
      "Animalia", "Acryllium", "Circaetus gallicus", "Leptocardii",
      "Oryx", "Oryx", "Coliiformes", "Caprimulgiformes", "Falconiformes",
      "Parotia", "Carduelis tristis", "Coelurosauria", "Anatalavis",
      "Diomedea exulans", "Burhinus", "Casuariiformes", "Archilochus colubris",
      "Fregata minor", "Bilateria", "Linnaeus", "Steatornithidae",
      "Confuciusornis sanctus", "Cyclostomata", "Cyclostomata", "Sphenisciformes",
      "Chionis", "Dipnoi", "Sophia", "Sophia", "Sophia", "Myxini",
      "Caracara", "Patricia", "Patricia", "Hydrophasianus", "Podicipediformes",
      "Geococcyx californianus", "Anatidae", "Archaeopteryx", "Archaeopteryx",
      "Piciformes", "Accipitriformes", "Spheniscus magellanicus", "Numididae",
      "Sialia sialis", "Jeholornis", "Anseranatidae", "Passeriformes",
      "Jixiangornis", "Aquila chrysaetos", "Hyperoartia", "Microraptor",
      "Apatornis", "Gaviiformes", "Vanellus", "Anseriformes", "Canaria",
      "Diapsida", "Phoenicopteriformes", "Anseranas", "Apterygiformes",
      "Psittacidae", "Ciconiiformes", "Malurus coronatus", "Megapodiidae",
      "Geospiza scandens", "Chauna", "Chauna", "Taeniopygia guttata",
      "Scansoriopterygidae", "Bubulcus ibis", "Talegalla", "Rahonavis",
      "Parus major", "Dinornithiformes", "Dromaeosauridae", "Asio",
      "Myiopsitta monachus", "Synthliboramphus antiquus", "Diomedea immutabilis",
      "Fratercula arctica", "Archelosauria", "Mammalia", "Merganetta",
      "Pteroclidiformes", "Puffinus griseus", "Ascidiacea", "Chondrichthyes",
      "Austinornis lentus", "Coraciiformes", "Pelecaniformes", "Grus grus",
      "Anhima", "Galliformes", "Eulipoa", "Eulipoa", "Coliidae", "Psittaciformes",
      "Agnatha", "Leipoa", "Caprimulgidae", "Aurornis xui", "Xiaotingia",
      "Crax", "Actinopterygii", "Aythya valisineria", "Strigiformes",
      "Cephalochordata", "Trogoniformes", "Passer domesticus", "Paradisaea raggiana",
      "Loriculus", "Turdoides bicolor", "Troodontidae", "Nyctea scandiaca",
      "Reptilia", "Protopteryx", "Accipitridae", "Hongshanornithidae",
      "Charadriiformes", "Alectura", "Alectura", "Crocodilia", "Uria aalge",
      "Alexandra", "Alexandra", "Alexandra", "Archaeopteryx lithographica",
      "Tinamiformes", "Nyctibiidae", "Squamata", "Caprimulgus tristigma",
      "Probosciger aterrimus", "Xiaotingia zhengi", "Apus apus", "Sarcopterygii",
      "Cracidae", "Opisthocomiformes", "Tachyeres", "Chordata", "Anchiornis huxleyi",
      "Aepypodius", "Anapsida", "Rheiformes", "Palaeognathae", "Apodiformes",
      "Theropoda", "Ichthyornis", "Dinosauria", "Actophilornis", "Neognathae",
      "Gruiformes", "Apsaravis", "Pitohui", "Alectura lathami", "Jacana",
      "Deuterostomia", "Macrocephalon", "Gansus", "Pachyptila belcheri",
      "Vertebrata", "Vertebrata", "Vertebrata", "Sapeornis", "Cathartes aura",
      "Gnathostomata", "Gnathostomata", "Maina", "Patagopteryx", "Cuculiformes",
      "Urochordata", "Songlingornithidae", "Plectropterus", "Phalaropus lobatus",
      "Sylvia", "Sylvia", "Irediparra", "Aves", "Phasianidae", "Lepidosauria",
      "Osteichthyes", "Musophagiformes", "Eurypyga helias", "Archosauria",
      "Oxyura vittata", "Columbiformes", "Procellariiformes", "Caprimulgus ruficollis",
      "Galloanserae", "Appendicularia", "Appendicularia", "Appendicularia",
      "Appendicularia", "Appendicularia", "Latina", "Oreortyx", "Odontophoridae",
      "Rhynchortyx", "Phasianinae", "Baso", "Agelastes", "Callipepla",
      "Callipepla", "Callipepla", "Tetraoninae", "Cyrtonyx", "Colinus",
      "Colinus", "Dactylortyx", "Basa", "Perdicinae", "Columbea", "Odontophorus",
      "Odontophorus", "Dendrortyx", "Europaea", "Guttera", "Philortyx",
      "Numida", "Meleagridinae"),
      c("felis", "canidae"))
    n <- 1
    med_nj <- med_upgma <- #vector(,length(taxa_list))
    sdm_nj <- sdm_upgma <-
    med_nj_cm <- med_upgma_cm <-
    sdm_nj_cm <- sdm_upgma_cm <-
    maxbr_all <- true_all <- good_names <- vector(,length(taxa_list))
    phylo_all <- phylo_med_nj <- phylo_med_upgma <- phylo_sdm_nj <-
    phylo_sdm_upgma <- vector(mode="list", length(taxa_list))

    for (taxa in taxa_list){
        query <- make_datelife_query(input = taxa, get_spp_from_taxon = TRUE)
        res <- get_datelife_result(input = query)
        print(n)
        if(length(res) == 0) {
          med_nj[n] <- med_upgma[n] <- med_nj_cm[n] <-
          med_upgma_cm[n] <- sdm_nj[n] <- sdm_upgma[n] <-
          sdm_nj_cm[n] <- sdm_upgma_cm[n] <- maxbr_all[n] <- "0 trees"
          n <- n + 1
          next()
        }
        if(length(res) == 1){
          med_nj[n] <- med_upgma[n] <- med_nj_cm[n] <- med_upgma_cm[n] <-
          sdm_nj[n] <- sdm_upgma[n] <- sdm_nj_cm[n] <- sdm_upgma_cm[n] <-
          maxbr_all[n] <- "1 tree"
          n <- n + 1
          next()
        }
        good_names[n] <- TRUE
        mrcas <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "mrca")
        names(mrcas) <- NULL
        phylo_all[[n]] <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_all")
        maxbr <- sapply(phylo_all[[n]], function(x) max(ape::branching.times(x)))
        names(maxbr) <- NULL
        maxbr_all[n] <- max(maxbr)
        #the following tests that both the source patristic matrix and the derived phylo object have the same maximum age.
        true_all[n] <- all(round(sort(mrcas), digits = 5) == round(sort(maxbr), digits = 5))
        ###############
        median.result <- NULL
        overlap <- 2
        while(is.null(median.result)){
          best_grove <- datelife::filter_for_grove(datelife_result = res,
                        criterion = "taxa", n = overlap)
                        patristic.array <- patristic_matrix_list_to_array(best_grove)
                        median.matrix <- summary_patristic_matrix_array(patristic.array)
          median.result <- tryCatch(patristic_matrix_to_phylo(median.matrix,
                            clustering_method = "nj", fix_negative_brlen = FALSE)
                            , error = function(e) NULL)
          overlap <- overlap + 1
        }
        phylo_med_nj[[n]] <- median.result
        med_nj[n] <- tryCatch(round(max(ape::branching.times(median.result)), digits = 3),
            error = function(e) NaN)
        med_nj_cm[n] <- tryCatch(median.result$clustering_method,
            error = function(e) NaN)
        median.result_upgma <- tryCatch(patristic_matrix_to_phylo(median.matrix,
                          clustering_method = "upgma", fix_negative_brlen = FALSE)
                          , error = function(e) NULL)
        phylo_med_upgma[[n]] <- median.result_upgma
        med_upgma[n] <- tryCatch(round(max(ape::branching.times(median.result_upgma)), digits = 3),
            error = function(e) NaN)
        med_upgma_cm[n] <- tryCatch(median.result_upgma$clustering_method,
            error = function(e) NaN)
        ################
        weighting = "flat"
        verbose = FALSE
        phy <- NA
        # used.studies <- names(datelife_result)
        unpadded.matrices <- lapply(best_grove, patristic_matrix_unpad)
        good.matrix.indices <- c()
        for(i in sequence(length(unpadded.matrices))) {
          test.result <- NA
          try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent = TRUE)
          if (verbose){
            message(cat(i, "out of", length(unpadded.matrices), "chronograms tried: "), appendLF = FALSE)
          }
          if(is.finite(test.result)) {
            good.matrix.indices <- append(good.matrix.indices,i)
            if (verbose){
              message(cat(" Ok."))
            }
          } else {
            if (verbose){
              message(cat(" Failed."))
            }
          }
        }
        if(length(good.matrix.indices) > 1) {
          unpadded.matrices <- unpadded.matrices[good.matrix.indices]
          # used.studies <- used.studies[good.matrix.indices]
          weights = rep(1, length(unpadded.matrices))
          if (weighting=="taxa") {
            weights = unname(sapply(unpadded.matrices, dim)[1,])
          }
          if (weighting=="inverse") {
            weights = 1/unname(sapply(unpadded.matrices, dim)[1,])
          }
          if (verbose){
            message(cat("\n", "Synthesizing", length(unpadded.matrices), "chronograms with SDM"))
          }
          SDM.result <- do.call(ape::SDM, c(unpadded.matrices, weights))[[1]]
        } else {
            if(length(good.matrix.indices) == length(best_grove)) {
              warning("There are not enough input chronograms to run SDM. This is not your fault.")
            } else {
              warning("All input chronograms throw an error when running SDM. This is not your fault.")
            }
            SDM.result <- NA
        }
        sdm.result <- tryCatch(patristic_matrix_to_phylo(SDM.result, clustering_method = "nj",
                      fix_negative_brlen = FALSE),
                      error = function(e) NA)
        phylo_sdm_nj[[n]] <- sdm.result
        sdm_nj[n] <- tryCatch(round(max(ape::branching.times(sdm.result)), digits = 3),
          error = function(e) NaN)
        sdm_nj_cm[n] <- tryCatch(sdm.result$clustering_method,
            error = function(e) "none")

        sdm.result_upgma <- tryCatch(datelife_result_sdm(best_grove, clustering_method = "upgma",
                      fix_negative_brlen = FALSE),
                      error = function(e) NA)
        phylo_sdm_upgma[[n]] <- sdm.result_upgma
        sdm_upgma[n] <- tryCatch(round(max(ape::branching.times(sdm.result_upgma)), digits = 3),
            error = function(e) NaN)
        sdm_upgma_cm[n] <- tryCatch(sdm.result_upgma$clustering_method,
          error = function(e) "none")
        n <- n + 1
    }
#end of test that
# })
# utils::data(names_subset2)
# taxa_list <- list(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis"))
# taxa_list <- c(taxa_list, lapply(1:50, function(x) sample(names_subset2, size = 20, replace = FALSE)))
# source("~/Desktop/datelife/data-raw/test_taxa_batch.R")
x <- cbind(med_nj, med_upgma, sdm_nj, sdm_upgma, maxbr_all, true_all, sapply(phylo_all, length))
final_test <- taxa_list[good_names]
# final_data <- x[good_names,]
final_data <- x[c(1,2,10:13,22,34,37,41,52),]
final_trees <- phylo_all[good_names]
save(final_test, final_data, final_trees, file = "test_taxa.RData")
lapply(final_trees[c(4,5,8,9,10)], names)
final_test[c(4,5,8,9,10)]

# names(med_nj) <- med_nj_cm
# names(med_upgma) <- med_upgma_cm
# names(sdm_nj) <- sdm_nj_cm
# names(sdm_upgma) <- sdm_upgma_cm
# med_nj
# med_upgma
# sdm_nj
# sdm_upgma
