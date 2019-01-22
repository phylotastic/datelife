library(drake)
library(microbenchmark)
library(testthat)

plan <- drake_plan(
	setwd("~/Desktop/datelife/data-raw/"),
	ninput = c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000),
	res = c(),
	for(i in ninput){
		xname = paste0("runtime_tests/2_tests/1_same_spp_names/aves",i,".1.gfr.runtime_2017.12.29")
		x = paste0(xname, ".RData")
		load(x)
		res = rbind(res, get(xname))
	},
	et1 = expect_true(length(res) == 2),
	et2 = expect_true(inherits(res, "microbenchmark")),
	load('runtime_tests/2_tests/0_all_names/aves.all.gfr.runtime_2017.12.29.RData'),
	expect_true(inherits(aves.all.gfr.runtime_2017.12.29, "microbenchmark")),
	resf = rbind(res, aves.all.gfr.runtime_2017.12.29),
	plt = microbenchmark:::autoplot.microbenchmark(resf),
	report1 = rmarkdown::render(
    		knitr_in("runtime_tests/report1.Rmd"),
    		output_file = file_out("report1.html"),
    		quiet = TRUE
  	)	
)
config <- drake_config(plan)
vis_drake_graph(config)
make(plan)

