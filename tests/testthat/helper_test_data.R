# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- make_datelife_query(threebirds)
threebirds_result <- get_datelife_result(threebirds_query)
threebirds_median <- summarize_datelife_result(threebirds_result, threebirds_query,
    summary_format = "phylo_median")

utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
utils::data(subset2_search)
subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result)$best_grove
# subset2_bestgrove_phyloall <- summarize_datelife_result(subset2_bestgrove, summary_format = "phylo_all")
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

birds_yellowstone <- c("Caracara", "Chen caerulescens", "Gavia immer", "Chen rossii",
"Branta canadensis", "Cygnus columbianus", "Aix sponsa", "Podiceps grisegena",
"Anas strepera", "Podiceps nigricollis", "Anas americana", "Anas platyrhynchos",
"Cygnus buccinator", "Anas discors", "Podilymbus podiceps", "Podiceps auritus",
"Aechmophorus occidentalis", "Aechmophorus clarkii", "Anas cyanoptera", "Anas aclypeata",
"Anas acuta", "Anas crecca", "Aythya valisineria", "Aythya americana", "Aythya collaris",
"Aythya marila", "Aythya affinis", "Histrionicus histrionicus", "Bucephala albeola",
"Bucephala clangula", "Bucephala islandica", "Lophodytes cucullatus", "Mergus merganser",
"Mergus serrator", "Oxyura jamaicensis", "Phalacrocorax auritus", "Buteo jamaicensis",
"Buteo regalis", "Buteo lagopus", "Aquila chrysaetos", "Rallus limicola", "Porzana carolina",
"Fulica americana", "Grus canadensis", "Pluvialis squatarola", "Charadrius vociferus",
"Ardea herodias", "Plegadis chihi", "Cathartes aura", "Pelecanus erythrorhynchos",
"Nycticorax nycticorax", "Pandion haliaetus", "Haliaeetus leucophalus", "Circus cyaneus",
"Accipiter striatus", "Accipiter cooperii", "Bonasa umbellus", "Dendragapus obscurus",
"Charadrius semipalmatus", "Accipiter gentilis", "Buteo platypterus", "Buteo swainsoni",
"Himantopus mexicanus", "Recurvirostra americana", "Actitis macularia", "Tringa solitaria",
"Tringa semipalmata", "Tringa flavipes", "Tringa melanoleuca", "Bartramia longicauda",
"Numenius americanus", "Limosa fedoa", "Arenaria interpres", "Calidris alba",
"Calidris pusilla", "Calidris mauri", "Calidris minutilla", "Calidris fusciollis",
"Calidris bairdii", "Calidris melanotos", "Limnodromus griseus", "Limnodromus scolopaceus",
"Gallinago delicata", "Phalaropus tricolor", "Phalaropus lobatus", "Empidonax hammondii",
"Empidonax oberholseri", "Empidonax occidentalis", "Tyrannus verticalis", "Tyrannus tyrannus",
"Lanius ludovicianus", "Lanius excubitor", "Vireo gilvus", "Larus pipixcan",
"Perisoreus canadensis", "Larus delawarensis", "Larus californicus", "Larus argentatus",
"Cyanocitta stelleri", "Sterna caspia", "Cyanocitta cristata", "Sterna hirundo",
"Nucifraga columbiana", "Sterna forsteri", "Pica hudsonia", "Gymnorhinus cyanocephalus",
"Columba livia", "Zenaida macroura", "Eremophila alpestris", "Tachycineta bicolor",
"Tachycineta thalassina", "Streptopelia decaocto", "Megascops kennicottii",
"Bubo virginianus", "Corvus brachyrhynchos", "Corvus corax", "Strix nebulosa",
"Riparia riparia", "Asio otus", "Glaucidium gnoma", "Asio flammeus", "Aegolius funereus",
"Aegolius acadicus", "Chordeiles minor", "Aeronautes saxatalis", "Selasphorus platycercus",
"Selasphorus rufus", "Stellula calliope", "Stelgidopteryx serripennis", "Petrochelidon pyrrhonota",
"Hirundo rustica", "Spizella", "Spizella breweri", "Oreothlypis ruficapilla",
"Geothlypis tolmiei", "Geothlypis trichas", "Setophaga ruticilla", "Setophaga fusca",
"Setophaga petechia", "Setophaga striata", "Setophaga coronata", "Setophaga townsendi",
"Cardellina pusilla", "Pipilo chlorurus", "Pipilo maculatus", "Spizella arborea",
"Pooecetes gramineus", "Chondestes grammacus", "Ammodramus savannarum", "Passerella iliaca",
"Melospiza melodia", "Poecile atricapillus", "Poecile", "Melospiza lincolnii",
"Zonotrichia leucophrys", "Certhia americana", "Salpinctes obsoletus", "Troglodytes aedon",
"Cistothorus palustris", "Sitta canadensis", "Sitta carolinensis", "Junco hyemalis",
"Picoides dorsalis", "Piranga ludoviciana", "Pheucticus melanocephalus", "Passerina amoena",
"Agelaius phoeniceus", "Sturnella neglecta", "Xanthocephalus xanthocephalus", "Regulus satrapa",
"Euphagus cyanocephalus", "Molothrus ater", "Icterus bullockii", "Sialia mexicana",
"Leucosticte atrata", "Pinicola enucleator", "Carpodacus cassinii", "Loxia curvirostra",
"Acanthis flammea", "Regulus calendula", "Sialia currucoides", "Myadestes townsendi",
"Catharus fuscescens", "Catharus ustulatus", "Leucosticte tephroctis", "Loxia leucoptera",
"Catharus guttatus", "Turdus migratorius", "Spinus pinus", "Dumetella carolinensis",
"Spinus tristis", "Orescoptes montanus", "Sturnus vulgarus", "Passer domesticus",
"Falco sparverius", "Falco columbarius", "Falco peregrinus", "Falco mexicanus",
"Anthus rubescens", "Empidonax traillii", "Contopus sordidulus", "Contopus cooperi",
"Picoides villosus", "Dryocopus pileatus", "Oreothlypis peregrina", "Colaptes auratus",
"Passerculus sandwichensis", "Picoides arcticus", "Picoides pubescens", "Mniotilta varia",
"Cinclus mexicanus", "Sphyrapicus nuchalis", "Sphyrapicus thyroideus", "Seiurus noveboracensis",
"Melanerpes erythrocephalus", "Ceryle alcyon", "Melanerpes lewis", "Plectrophenax nivalis",
"Oreothlypis celata", "Larus philadelphia", "Coccothraustes vespertinus",
"Bombycilla", "Bombycilla cedrorum")

birds_yellowstone_otoltree <- get_otol_synthetic_tree(birds_yellowstone)
names(birds_yellowstone_otoltree)
data.frame(birds_yellowstone_otoltree$tip.label, birds_yellowstone_otoltree$ott_ids)
birds_yellowstone_query <- make_datelife_query(input = birds_yellowstone_otoltree)
birds_yellowstone_result <- get_datelife_result(birds_yellowstone_query)
birds_yellowstone_phyloall <- summarize_datelife_result(birds_yellowstone_result, birds_yellowstone_query, "phylo_all")
# unname(ape::is.ultrametric(phylo_all))
# unname(sapply(phylo_all, function(x) max(ape::branching.times(x))))
# unname(sapply(datelife_result, function(x) max(x)))
# plot_phylo_all(birds_yellowstone_phyloall)
birds_yellowstone_phylomedian <- summarize_datelife_result(birds_yellowstone_result, birds_yellowstone_query, "phylo_median")
plot(birds_yellowstone_phylomedian, cex = 0.25)
ape::axisPhylo(cex = 0.5)
