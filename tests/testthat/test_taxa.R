# in this file we:
 # debug by running the main functions over and over but with different set of taxa
 # register inexpected results from tnrs services to be hopefully corrected in the future

test_that("birds from wikipedia work", {
  taxa <- c("Yixianornis grabaui", "Amphibia", "Amphibia", "Amphibia",
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
    "Numida", "Meleagridinae")
  expect_true(inherits(datelife_search(taxa, summary_format="phylo_median"), "phylo"))
  expect_true(inherits(datelife_search(taxa, summary_format="phylo_sdm"), "phylo"))
})

test_that("Mus higher-taxon search is giving species back", {
  skip("skipping Mus")
  expect_silent(make_datelife_query("Echinus", get_spp_from_taxon = TRUE))
  expect_silent(make_datelife_query("Mus", get_spp_from_taxon = TRUE))
  expect_true(length(rphylotastic::taxon_get_species("Mus")) > 0)
})

test_that("patristic_matrix_to_phylo is accurate", {
  skip_on_travis()
  skip("skipping taxa test")
  utils::data(names_subset2)
  taxa_list <- list(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis"))
  taxa_list <- c(taxa_list, lapply(1:50, function(x) sample(names_subset2, size = 20, replace = FALSE)))
  n <- 0
  phylocm_med <- c()
  phylocm_sdm <- c()
  phylobr_sdm <- c()
  phylobr_med <- c()
  maxbr_all <- c()
  for (taxa in taxa_list){
    query <- make_datelife_query(input = taxa, get_spp_from_taxon = TRUE)
    res <- get_datelife_result(input = query)
    if(length(res) == 0) {
      phylocm_med <- c(phylocm_med, NA)
      phylocm_sdm <- c(phylocm_sdm, NA)
      phylobr_med <- c(phylobr_med, NA)
      phylobr_sdm <- c(phylobr_sdm, NA)
      maxbr_all <- c(maxbr_all, NA)
      n <- n + 1
      print(n)
      next()
    }
    mrcas <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "mrca")
    names(mrcas) <- NULL
    phyloall <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_all")
    maxbr <- sapply(phyloall, function(x) max(ape::branching.times(x)))
    names(maxbr) <- NULL
    #the following tests that both the source patristic matrix and the derived phylo object have the same maximum age.
    expect_true(all(round(sort(mrcas), digits = 5) == round(sort(maxbr), digits = 5)))
    phylosdm <- tryCatch(summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_sdm"),
      error = function(e) NA)
    # expect_true(max(maxbr)*1.2 > max(ape::branching.times(phylosdm)) & max(ape::branching.times(phylosdm)) > max(maxbr)*0.8)
    phylomed <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_median")
    # expect_true(max(maxbr)*1.2 > max(ape::branching.times(phylomed)) & max(ape::branching.times(phylomed)) > max(maxbr)*0.8)
    n <- n + 1
    phylocm_med <- c(phylocm_med, phylomed$clustering_method)
    phylocm_sdm <- c(phylocm_sdm, tryCatch(phylosdm$clustering_method,
      error = function(e) NaN))
    phylobr_med <- c(phylobr_med, max(ape::branching.times(phylomed)))
    phylobr_sdm <- c(phylobr_sdm, tryCatch(max(ape::branching.times(phylosdm)),
      error = function(e) NaN))
    maxbr_all <- c(maxbr_all, max(maxbr))
    print(n)
  }
})
