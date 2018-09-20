test_that("felis/canidae divergence is accurate", {
    query <- make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
    cats_and_dogs_results <- get_datelife_result(input = query)
    matrix_max_ages <- sapply(cats_and_dogs_results, max)
    taxa <- c("felis", "canidae")
    cats_and_dogs <- datelife_search(input = taxa, get_spp_from_taxon = TRUE,
      summary_format = "phylo_all")
    phylo_max_ages <- sapply(cats_and_dogs, function(x) max(ape::branching.times(x)))
    expect_true(all(names(matrix_max_ages) == names(phylo_max_ages)))
    # names(matrix_max_ages) <- names(phylo_max_ages)<- NULL
    ns <- 20
    format(round(sort(matrix_max_ages/2)), nsmall = ns) == format(round(sort(phylo_max_ages)), nsmall = ns)
    # ages from our cache range from 54.9 to 70.9, this includes upper limit confidence interval chronograms
    # timetree study-derived ages range from 39.7 to 67.1. This excludes confidence intervals
    median_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
      summary_format = "phylo_median")
    sdm_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
        summary_format = "phylo_sdm")
})

test_that("Mus higher-taxon search is giving species back"){
  expect_silent(make_datelife_query("Echinus", get_spp_from_taxon = TRUE))
  expect_silent(make_datelife_query("Mus", get_spp_from_taxon = TRUE))
  expect_true(length(rphylotastic::taxon_get_species("Mus")) > 0)
})

test_that("birds and cat sdm is super young"){
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis")
  query <- make_datelife_query(input = taxa, get_spp_from_taxon = TRUE)
  res <- get_datelife_result(input = query)
  all_phylo <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_all")
  sdm_phylo <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "phylo_sdm")
  mrcas <- summarize_datelife_result(datelife_query= query, datelife_result = res, summary_format = "mrca")
  names(mrcas) <- NULL
  max(ape::branching.times(sdm_phylo))
})

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
  expect_equal(class(datelife_search(taxa, summary_format="phylo_median")), "phylo")
  expect_equal(class(datelife_search(taxa, summary_format="phylo_sdm")), "phylo")
})

test_that("median and sdm work ok with very variable source chronograms"){
    utils::data(names_subset2)
    spp_query <- make_datelife_query(names_subset2)
    spp_dl_result <- get_datelife_result(spp_query)
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_median")), "phylo")
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_sdm")), "phylo")
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_median")), "character")
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_sdm")), "character")
}
