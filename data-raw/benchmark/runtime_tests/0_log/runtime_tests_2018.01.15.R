# Lunes 15 de enero 2018

################################
devtools::install_github("phylotastic/datelife", force=TRUE)
source('~/Desktop/datelife/R/datelife.R')
library(datelife)
library(ape)
library(geiger)
library(testthat)
################################

load('~/Google Drive/datelife/runtime_tests/1_name_samples/random_sample_10_aves_spp.RData')
make_bold_otol_tree(input=random_sample_10_aves_spp[[3]])
# Error in make_bold_otol_tree(input = get(xname)[[j]]) :
  # Not enough sequence coverage in BOLD to perform analysis of this set of taxa
j <- 12
i <- 10
input <- random_sample_10_aves_spp[[12]]
input <- input_check(input)  # ok
input <- input$cleaned.names

# fix another error in make_bold_otol_tree:
# Error in check_ott_ids(ott_ids) : NAs are not allowed
# In addition: Warning messages:
# 1: In make_bold_otol_tree(input = get(xname)[[j]]) :
  # Not enough sequences available in BOLD. No tree was constructed.
# 2: In check_tnrs(res) :
  # Mesopicos elliotii, Primobucco olsoni, Motacilla feldegg are not matched
# 3: In rotl::tnrs_match_names(names = input) : NAs introduced by coercion
j <- 14
i <- 10
input <- random_sample_10_aves_spp[[14]]
rr <- rotl::tnrs_match_names(names = input)
rr <- rr[!is.na(rr$unique_name),]
taxon.index <- which(grepl(i, sequences$species_name))  # what happens here if there are no sequences from the taxon????
# taxon.index is empty in this case
# I thought it was taking species with different names, so was trying to find a way to homogenize names, using ott_ids from rotl:
?bold::bold_seqspec
# use ott_id from rr to avoid name mispecifications:
names(sequences)
sequences$species_taxID
rr$ott_id
for(i in rr$ott_id){
	y <- i %in% sequences$species_taxID
	print(y)
}

# yet another error in make_bold_otol_tree:
# Error: $ operator is invalid for atomic vectors
j <- 83
i <- 10
input <- random_sample_10_aves_spp[[83]]
# Here's the error, input names do not match names from sequences obtained with bold::bold_seqspec(taxon = input, marker = marker)
	sequences <- bold::bold_seqspec(taxon = input, marker = marker)
# so just have to move the conditional step length(sequences) == 1 before

# another error in make_bold_otol_tree:
# Error in data[, tree$tip.label] : subscript out of bounds
j <- 34
i <- 100
load('~/Google Drive/datelife/runtime_tests/1_name_samples/random_sample_100_aves_spp.RData')
input <- random_sample_100_aves_spp[[34]]
# I think input names in phy do not match any of the seqeunces found in bold
rownames(final.sequences) %in% phy$tip.label
# nope, they all match but one:
ww <- which(rownames(final.sequences) %in% phy$tip.label)
rownames(final.sequences)[ww]
length(phy$tip.label)  # 94
length(rr$unique_name)  # 99
length(mm)  # 94
phy$tip.label[34]  # "Asthenes"
?tol_induced_subtree
 [1] "Serinus_melanochrous"        "Emberiza_striolata"
 [3] "Emberiza_sulphurata"         "Ammodramus_savannarum"
 [5] "Icterus_chrysater"           "Sporophila_leucoptera"
 [7] "Passerina_leclancherii"      "Ploceus_victoriae"
 [9] "Ploceus_castanops"           "Vidua_funerea"
[11] "Pyrenestes_sanguineus"       "Turdus_smithi"
[13] "Turdus_fumigatus"            "Cyornis_hoevelli"
[15] "Polioptila_melanura"         "Randia_pseudozosterops"
[17] "Prinia_erythroptera"         "Pnoepyga_albiventer"
[19] "Zosterops_vellalavella"      "Cyanistes_semilarvatus"
[21] "Rhipidura_brachyrhyncha"     "Rhipidura_lepida"
[23] "Hypothymis_azurea"           "Pyrrhocorax_graculus"
[25] "Cyanocorax_chrysops"         "Dicrurus_ludwigii"
[27] "Bias_musicus"                "Laniarius_nigerrimus"
[29] "Sphecotheres_vieilloti"      "Pachycephala_melanura"
[31] "Archboldia_papuensis"        "Cranioleuca_pallida"
[33] "Thripophaga_berlepschi"      "Asthenes"
[35] "Margarornis_stellatus"       "Cinclodes_palliatus"
[37] "Grallaria_alleni"            "Hylopezus_dives"
[39] "Rhytipterna_holerythra"      "Contopus_pallidus"
[41] "Ochthoeca_jelskii"           "Fluvicola_albiventer"
[43] "Myiophobus_lintoni"          "Psephotellus_pulcherrimus"
[45] "Micropsitta_meeki"           "Poicephalus_flavifrons"
[47] "Pionus_cyanescens"           "Amazona_autumnalis"
[49] "Falco_vespertinus"           "Megalaima_rafflesii"
[51] "Selenidera_spectabilis"      "Picoides_villosus"
[53] "Colaptes_rivolii"            "Chrysocolaptes_validus"
[55] "Prodotiscus_zambesiae"       "Todiramphus_lazuli"
[57] "Bycanistes_brevis"           "Rhinopomastus_cyanomelas"
[59] "Mascarenotus_murivorus"      "Glaucidium_brasilianum"
[61] "Circaetus_fasciolatus"       "Podiceps_auritus"
[63] "Calidris_subminuta"          "Actitis_macularia"
[65] "Numenius_minutus"            "Charadrius_javanicus"
[67] "Burhinus_oedicnemus"         "Bubulcus_ibis"
[69] "Mesembrinibis_cayennensis"   "Phalacrocorax_varius"
[71] "Spheniscus_magellanicus"     "Aptenodytes_australis"
[73] "Hydrobates_furcatus"         "Aenigmatolimnas_marginalis"
[75] "Gallirallus_philippensis"    "Geophaps_smithii"
[77] "Geophaps_scripta"            "Ptilinopus_chalcurus"
[79] "Ptilinopus_huttoni"          "Leptotila_rufaxilla"
[81] "Tauraco_ruspolii"            "Coua_verreauxi"
[83] "Coua_coquereli"              "Centropus_celebensis"
[85] "Anthracothorax_prevostii"    "Colibri_coruscans"
[87] "Eriocnemis_alinae"           "Sephanoides_fernandensis"
[89] "Chrysolophus_amherstiae"     "Malacorhynchus_membranaceus"
[91] "Cygnus_buccinator"           "Virago_castanea"
[93] "Pyropia_pulchra"             "Dictyopteris_gracilis"
# found the problem. We're using species names instead of ott_ids to match input names and rotl_induced_subtree.
# fixed it in datelife now
