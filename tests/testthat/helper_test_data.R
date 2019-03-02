# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- make_datelife_query(threebirds)
threebirds_result <- get_datelife_result(threebirds_query)

utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
utils::data(subset2_search)
subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result)$best_grove
# trees <- lapply(subset2_bestgrove, patristic_matrix_to_phylo)

subset2_bestgrove_phyloall <- summarize_datelife_result(subset2_bestgrove, summary_format = "phylo_all")
# subset2_bestgrove_consensus <- ape::consensus()
subset2_sdm_matrix <- get_sdm_matrix(subset2_bestgrove)

birds_wiki <- c("Yixianornis grabaui", "Amphibia", "Amphibia", "Amphibia",
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
