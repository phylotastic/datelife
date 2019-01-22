# Miércoles 21 de Febrero 2018

# 1. Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# first with doML=FALSE
# finished up to 1500 spp
# # ran the test from 300 again, and got the following: 
[1] 300
[1] 400
[1] 500
[1] 700
[1] 1000
[1] 1500
Error in b_GET(paste0(bbase(), "API_Public/combined"), args, ...) :
  <style type="text/css">
  #kohana_error { background: #ddd; font-size: 1em; font-family:sans-serif; text-align: left; color: #111; }
  #kohana_error h1,
  #kohana_error h2 { margin: 0; padding: 1em; font-size: 1em; font-weight: normal; background: #911; color: #fff; }
  #kohana_error h1 a,
  #kohana_error h2 a { color: #fff; }
  #kohana_error h2 { background: #222; }
  #kohana_error h3 { margin: 0; padding: 0.4em 0 0; font-size: 1em; font-weight: normal; }
  #kohana_error p { margin: 0; padding: 0.2em 0; }
  #kohana_error a { color: #1b323b; }
  #kohana_error pre { overflow: auto; white-space: pre-wrap; }
  #kohana_error table { width: 100%; display: block; margin: 0 0 0.4em; padding: 0; border-collapse: collapse; background: #fff; }
#kohana_error table td { border: solid 1px #ddd; text-align: left; vertical-align: top; padding: 0.4em; }
#kohana_error div.content { padding: 0.4em 1em 1em; overflow: hidden; }
#kohana_error pre.source { margin: 0 0 1em; padding: 0.4em; background: #fff; border: d
In addition: There were 50 or more warnings (use warnings() to see the first 50)
# warning came out again, check it:
# 50: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
#               EOF within quoted string
> warnings()
Warning messages:
  1: In make_bold_otol_tree(input = get(xname)[[1]]) :
  Negative branch lengths in BOLD chronogram.

2: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
             EOF within quoted string
3: In make_bold_otol_tree(input = get(xname)[[j]]) :
             Negative branch lengths in BOLD chronogram.
           
4: In make_bold_otol_tree(input = get(xname)[[j]]) :
             Negative branch lengths in BOLD chronogram.
           
5: In make_bold_otol_tree(input = get(xname)[[j]]) :
             Negative branch lengths in BOLD chronogram.
           
6: In make_bold_otol_tree(input = get(xname)[[j]]) :
             Negative branch lengths in BOLD chronogram.
           
7: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
           
8: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
           
9: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
10: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
11: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
12: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
13: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
14: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
15: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                      
16: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
17: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                  
18: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                  
19: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                  
20: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
21: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                              
22: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                              
23: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
24: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                          
25: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
26: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                      
27: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                      
28: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
29: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                  
30: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                  
31: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                  
32: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                  
33: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
34: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
35: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                                          
36: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                                          
37: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
38: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
                                                                                                                      
39: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
                                                                                                                                  40: In make_bold_otol_tree(input = get(xname)[[j]]) :
 Negative branch lengths in BOLD chronogram.
41: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
42: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
43: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
44: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
45: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
46: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
47: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string
48: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
49: In make_bold_otol_tree(input = get(xname)[[j]]) :
Negative branch lengths in BOLD chronogram.
50: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
EOF within quoted string

i <- 2000
xname <- paste0("random_sample_",i, "_aves_spp")
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
load(file=paste0(xname,".RData"))
make_bold_otol_tree(input=random_sample_2000_aves_spp[[1]])
# 
# Phylogenetic tree with 475 tips and 474 internal nodes.
# 
# Tip labels:
#   Lepidothrix coronata, Chrysolampis mosquitus, Oenanthe moesta, Northiella haematogaster, Numida meleagris, Spinus magellanicus, ...
# Node labels:
#   Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...
# 
# Rooted; includes branch lengths.
# Warning messages:
#   1: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#                EOF within quoted string
#              2: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#                           EOF within quoted string
#                         3: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#                                      EOF within quoted string
#                                    4: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[1]]) :
#                                      Negative branch lengths in BOLD chronogram.
make_bold_otol_tree(input=random_sample_2000_aves_spp[[2]])
# Phylogenetic tree with 547 tips and 546 internal nodes.

# Tip labels:
        # Sayornis nigricans, Buteo brachyurus, Aphelocephala leucopsis, Penelope marail, Porzana fusca, Columba rupestris, ...
# Node labels:
        # , Aves ott81461, Neognathae ott241846, , , , ...

# Rooted; includes branch lengths.
# Warning messages:
# 1: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string
# 2: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string
# 3: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string
# 4: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[2]]) :
        # Negative branch lengths in BOLD chronogram.
make_bold_otol_tree(input=random_sample_2000_aves_spp[[3]])
# Phylogenetic tree with 518 tips and 517 internal nodes.

# Tip labels:
        # Heliobletus contaminatus, Hylophylax punctulatus, Tragopan blythii, Picumnus spilogaster, Lamprotornis pulcher, Megascops colombianus, ...
# Node labels:
        # Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...

# Rooted; includes branch lengths.
# Warning messages:
# 1: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string
# 2: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string
# 3: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[3]]) :
        # Negative branch lengths in BOLD chronogram.
for (i in 4:10){
	cat("number = ", i, "\n")
	x <- make_bold_otol_tree(input=random_sample_2000_aves_spp[[i]])
	print(x)
}
#retrieved exactly the same warnings
# number =  4

# Phylogenetic tree with 527 tips and 526 internal nodes.

# Tip labels:
        # Mycteria americana, Hylocharis eliciae, Botaurus lentiginosus, Veniliornis sanguineus, Cacatua galerita, Sylvia curruca, ...
# Node labels:
        # Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...

# Rooted; includes branch lengths.
# number =  5

# Phylogenetic tree with 535 tips and 534 internal nodes.

# Tip labels:
        # Spizella pallida, Aphrodroma brevirostris, Phleocryptes melanops, Myioborus castaneocapilla, Sitta przewalskii, Dacnis cayana, ...
# Node labels:
        # Aves ott81461, Neognathae ott241846, , , , , ...

# Rooted; includes branch lengths.
# number =  6

# Phylogenetic tree with 499 tips and 498 internal nodes.

# Tip labels:
        # Turdus merula, Poecile carolinensis, Calonectris leucomelas, Cyanistes caeruleus, Branta bernicla, Yuhina diademata, ...
# Node labels:
        # , Aves ott81461, Neognathae ott241846, , , , ...

# Rooted; includes branch lengths.
# number =  7

# Phylogenetic tree with 570 tips and 569 internal nodes.

# Tip labels:
        # Campethera caroli, Atthis heloisa, Cyanocompsa brissonii, Spinus cucullatus, Bulweria bulwerii, Chamaeza ruficauda, ...
# Node labels:
        # , Aves ott81461, Neognathae ott241846, , , , ...

# Rooted; includes branch lengths.
# number =  8

# Phylogenetic tree with 484 tips and 483 internal nodes.

# Tip labels:
        # Picoides nuttallii, Otus mindorensis, Amazilia edward, Gallirallus striatus, Phalacrocorax capillatus, Cinnyris bouvieri, ...
# Node labels:
        # , Aves ott81461, Neognathae ott241846, , , , ...

# Rooted; includes branch lengths.
# number =  9

# Phylogenetic tree with 531 tips and 530 internal nodes.

# Tip labels:
        # Laniarius luehderi, Rhipidura phasiana, Dumetella carolinensis, Geotrygon violacea, Clytolaema rubricauda, Pionus tumultuosus, ...
# Node labels:
        # Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...

# Rooted; includes branch lengths.
# number =  10

# Phylogenetic tree with 526 tips and 525 internal nodes.

# Tip labels:
        # Malimbus nitens, Tetraogallus altaicus, Tringa solitaria, Coracias garrulus, Emberiza caesia, Fulica armillata, ...
# Node labels:
        # , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

# Rooted; includes branch lengths.

for (i in 11:100){
	cat("number = ", i, "\n")
	x <- make_bold_otol_tree(input=random_sample_2000_aves_spp[[i]])
	print(x)
}
number =  11

Phylogenetic tree with 553 tips and 552 internal nodes.

Tip labels:
  Campylorhamphus procurvoides, Ardea cinerea, Bleda canicapillus, Prunella atrogularis, Lophozosterops dohertyi, Crax rubra, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  12

Phylogenetic tree with 496 tips and 495 internal nodes.

Tip labels:
  Oceanodroma leucorhoa, Brachyramphus perdix, Hypocryptadius cinnamomeus, Amalocichla incerta, Bleda canicapillus, Oenanthe isabellina, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  13

Phylogenetic tree with 547 tips and 546 internal nodes.

Tip labels:
  Erpornis zantholeuca, Conopophaga peruviana, Ictinia plumbea, Notopholia corrusca, Aplonis brunneicapillus, Knipolegus poecilurus, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  14

Phylogenetic tree with 488 tips and 487 internal nodes.

Tip labels:
  Megaceryle torquata, Elaenia mesoleuca, Cyanoloxia glaucocaerulea, Pelecanus rufescens, Tachycineta bicolor, Trogon melanurus, ...
Node labels:
  Aves ott81461, Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.
number =  15

Phylogenetic tree with 482 tips and 481 internal nodes.

Tip labels:
  Eurostopodus argus, Urosticte ruficrissa, Tangara florida, Bonasa umbellus, Phyllomyias griseiceps, Cinnyris venustus, ...
Node labels:
  , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  16

Phylogenetic tree with 473 tips and 472 internal nodes.

Tip labels:
  Callipepla gambelii, Asthenes pyrrholeuca, Pyrilia pulchra, Catharus fuscater, Microcerculus marginatus, Acridotheres cristatellus, ...
Node labels:
  Eukaryota ott304358, Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  17

Phylogenetic tree with 580 tips and 579 internal nodes.

Tip labels:
  Acrocephalus agricola, Coeligena helianthea, Dromas ardeola, Falco amurensis, Synallaxis brachyura, Botaurus pinnatus, ...
Node labels:
  Aves ott81461, Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.
number =  18

Phylogenetic tree with 536 tips and 535 internal nodes.

Tip labels:
  Melanitta fusca, Calidris pusilla, Anhinga rufa, Riparia diluta, Accipiter soloensis, Cyanoramphus auriceps, ...
Node labels:
  Eukaryota ott304358, Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  19

Phylogenetic tree with 554 tips and 553 internal nodes.

Tip labels:
  Dromas ardeola, Compsothraupis loricata, Zosterops everetti, Chlorospingus tacarcunae, Sclateria naevia, Psittacara holochlorus, ...
Node labels:
  , Amniota ott229560, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  20

Phylogenetic tree with 496 tips and 495 internal nodes.

Tip labels:
  Knipolegus hudsoni, Sylviorthorhynchus desmursii, Dicaeum australe, Icterus parisorum, Amazilia fimbriata, Spinus magellanicus, ...
Node labels:
  , Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.
number =  21

Phylogenetic tree with 586 tips and 585 internal nodes.

Tip labels:
  Larus vegae, Drymornis bridgesii, Haematopus longirostris, Polysticta stelleri, Calidris alba, Hypocolius ampelinus, ...
Node labels:
  , Amniota ott229560, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  22

Phylogenetic tree with 553 tips and 552 internal nodes.

Tip labels:
  Drymophila ochropyga, Psittacara wagleri, Sclerurus mexicanus, Coeligena torquata, Mesophoyx intermedia, Amazona barbadensis, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  23

Phylogenetic tree with 485 tips and 484 internal nodes.

Tip labels:
  Glareola nordmanni, Phoenicurus erythrogastrus, Bleda canicapillus, Elvira chionura, Todiramphus gambieri, Leucopsar rothschildi, ...
Node labels:
  Eukaryota ott304358, Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  24

Phylogenetic tree with 543 tips and 542 internal nodes.

Tip labels:
  Erythrura trichroa, Pipra fasciicauda, Dendrocincla merula, Hylophilus semicinereus, Pyrrhura melanura, Doryfera johannae, ...
Node labels:
  Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  25

Phylogenetic tree with 569 tips and 568 internal nodes.

Tip labels:
  Leptoptilos javanicus, Eupodotis afra, Cepphus grylle, Phoenicurus ochruros, Vidua macroura, Calidris fuscicollis, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  26

Phylogenetic tree with 401 tips and 400 internal nodes.

Tip labels:
  Turdoides malcolmi, Topaza pyra, Streptoprocne zonaris, Aplonis brunneicapillus, Brotogeris sanctithomae, Poecile carolinensis, ...
Node labels:
  Eukaryota ott304358, Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  27

Phylogenetic tree with 519 tips and 518 internal nodes.

Tip labels:
  Patagioenas leucocephala, Eleoscytalopus indigoticus, Saxicola insignis, Macronectes giganteus, Patagioenas maculosa, Phalaropus tricolor, ...
Node labels:
  Eukaryota ott304358, , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , ...

Rooted; includes branch lengths.
number =  28

Phylogenetic tree with 569 tips and 568 internal nodes.

Tip labels:
  Limnodromus semipalmatus, Tachyeres pteneres, Forpus passerinus, Agropsar philippensis, Larus glaucescens, Eumyias indigo, ...
Node labels:
  , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  29

Phylogenetic tree with 581 tips and 580 internal nodes.

Tip labels:
  Pogoniulus subsulphureus, Falco subbuteo, Melaenornis pammelaina, Oenanthe xanthoprymna, Drymornis bridgesii, Oenanthe picata, ...
Node labels:
  Eukaryota ott304358, , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , ...

Rooted; includes branch lengths.
number =  30

Phylogenetic tree with 508 tips and 507 internal nodes.

Tip labels:
  Myiothlypis leucophrys, Toxostoma bendirei, Icterus mesomelas, Sporophila pileata, Platalea ajaja, Sporophila caerulescens, ...
Node labels:
  Opisthokonta ott332573, , Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  31

Phylogenetic tree with 505 tips and 504 internal nodes.

Tip labels:
  Phoenicopterus chilensis, Vireo bellii, Phrygilus gayi, Apteryx australis, Tringa incana, Anser indicus, ...
Node labels:
  , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  32

Phylogenetic tree with 567 tips and 566 internal nodes.

Tip labels:
  Empidonax wrightii, Sturnia pagodarum, Phaethornis yaruqui, Formicarius analis, Phylloscopus orientalis, Pholia sharpii, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  33

Phylogenetic tree with 497 tips and 496 internal nodes.

Tip labels:
  Cyornis umbratilis, Buteogallus meridionalis, Charadrius forbesi, Porzana spiloptera, Falco vespertinus, Iduna caligata, ...
Node labels:
  Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  34

Phylogenetic tree with 566 tips and 565 internal nodes.
Warning message:
  In make_bold_otol_tree(input = asd.names) :
  Negative branch lengths in BOLD chronogram.

> get_datelife_result(asd, process_input=TRUE)

Tip labels:
  Arborophila gingica, Veniliornis sanguineus, Phaethornis malaris, Locustella ochotensis, Euphonia fulvicrissa, Psittacara leucophthalmus, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  94

Phylogenetic tree with 546 tips and 545 internal nodes.

Tip labels:
  Tinamus major, Pselliophorus luteoviridis, Onychognathus salvadorii, Amazona arausiaca, Myrmotherula menetriesii, Aythya australis, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  95

Phylogenetic tree with 562 tips and 561 internal nodes.

Tip labels:
  Aethia pusilla, Turdus hortulorum, Zentrygon costaricensis, Cuculus optatus, Pachyptila turtur, Euphonia pectoralis, ...
Node labels:
  Aves ott81461, Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.
number =  96

Phylogenetic tree with 506 tips and 505 internal nodes.

Tip labels:
  Aethopyga siparaja, Eumyias panayensis, Phaetusa simplex, Streptoprocne phelpsi, Micropsitta finschii, Herpsilochmus sticturus, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  97

Phylogenetic tree with 430 tips and 429 internal nodes.

Tip labels:
  Pluvialis fulva, Eurillas curvirostris, Copsychus luzoniensis, Ficedula semitorquata, Rollandia rolland, Apteryx australis, ...
Node labels:
  Eukaryota ott304358, , Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
number =  98

Phylogenetic tree with 571 tips and 570 internal nodes.

Tip labels:
  Sarcoramphus papa, Actenoides lindsayi, Myrmotherula schisticolor, Leucophaeus atricilla, Mino kreffti, Glareola pratincola, ...
Node labels:
  , Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  99

Phylogenetic tree with 467 tips and 466 internal nodes.

Tip labels:
  Phaethornis bourcieri, Euphonia chlorotica, Geositta punensis, Pedionomus torquatus, Pomatorhinus ruficollis, Zosterops citrinella, ...
Node labels:
  Eukaryota ott304358, Aves ott81461, Neognathae ott241846, , , , ...

Rooted; includes branch lengths.
number =  100

Phylogenetic tree with 391 tips and 390 internal nodes.

Tip labels:
  Turdoides gularis, Amazona festiva, Sporophila leucoptera, Falculea palliata, Aplonis pelzelni, Picoides scalaris, ...
Node labels:
  , Euteleostomi ott114654, Aves ott81461, Neognathae ott241846, , , ...

Rooted; includes branch lengths.
There were 50 or more warnings (use warnings() to see the first 50)
> warnings()
Warning messages:
  1: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
               EOF within quoted string
             2: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                          EOF within quoted string
                        3: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                          Negative branch lengths in BOLD chronogram.
                        
                        4: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                     EOF within quoted string
                                   5: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                EOF within quoted string
                                              6: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                           EOF within quoted string
                                                         7: In collapse_singles(tr) :
                                                           Dropping singleton nodes with labels: Spheniscidae ott494367
                                                         8: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                           Negative branch lengths in BOLD chronogram.
                                                         
                                                         9: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                      EOF within quoted string
                                                                    10: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                  EOF within quoted string
                                                                                11: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                  Negative branch lengths in BOLD chronogram.
                                                                                
                                                                                12: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                              EOF within quoted string
                                                                                            13: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                          EOF within quoted string
                                                                                                        14: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                      EOF within quoted string
                                                                                                                    15: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                  EOF within quoted string
                                                                                                                                16: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                  Negative branch lengths in BOLD chronogram.
                                                                                                                                
                                                                                                                                17: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                              EOF within quoted string
                                                                                                                                            18: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                          EOF within quoted string
                                                                                                                                                        19: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                    20: In collapse_singles(tr) :
                                                                                                                                                                      Dropping singleton nodes with labels: Spheniscidae ott494367
                                                                                                                                                                    21: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                    
                                                                                                                                                                    22: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                  EOF within quoted string
                                                                                                                                                                                23: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                              EOF within quoted string
                                                                                                                                                                                            24: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                          EOF within quoted string
                                                                                                                                                                                                        25: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                          Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                        
                                                                                                                                                                                                        26: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                                                                    27: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    28: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                  EOF within quoted string
                                                                                                                                                                                                                                29: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                              EOF within quoted string
                                                                                                                                                                                                                                            30: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                          EOF within quoted string
                                                                                                                                                                                                                                                        31: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                                                                                                                    32: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                    33: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                  EOF within quoted string
                                                                                                                                                                                                                                                                                34: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                              EOF within quoted string
                                                                                                                                                                                                                                                                                            35: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                          EOF within quoted string
                                                                                                                                                                                                                                                                                                        36: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                                                                                                                                                                    37: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                  EOF within quoted string
                                                                                                                                                                                                                                                                                                                                38: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                                                                                  Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                39: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                              EOF within quoted string
                                                                                                                                                                                                                                                                                                                                            40: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                          EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                        41: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                                    42: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                    43: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                    44: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                                  EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                                                45: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                                                                                                                                  Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                46: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                                              EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                                                            47: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                                                          EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                                                                        48: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                                                                      EOF within quoted string
                                                                                                                                                                                                                                                                                                                                                                                                                    49: In make_bold_otol_tree(input = random_sample_2000_aves_spp[[i]]) :
                                                                                                                                                                                                                                                                                                                                                                                                                      Negative branch lengths in BOLD chronogram.
                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                    50: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
                                                                                                                                                                                                                                                                                                                                                                                                                                  EOF within quoted string

# there appears to be no problem: rerun from 1500 names on a new screen 

