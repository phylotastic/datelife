library(drake)
library(microbenchmark)
library(testthat)
library(here)
setwd("~/Desktop/datelife/data-raw/benchmark/runtime_tests")
set_here()
plan <- drake_plan(
	ninput = c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000),
	res = c(),
	ff = for(i in ninput){
		xname = paste0("2_tests/1_same_spp_names/aves",i,".1.gfr.runtime_2017.12.29")
		x = paste0(xname, ".RData")
		load(x)
		res = rbind(res, get(xname))
	},
	# expect_true(length(res) == 2),
	# expect_true(inherits(res, "microbenchmark")),
	load("2_tests/0_all_names/aves.all.gfr.runtime_2017.12.29.RData"),
	# expect_true(inherits(aves.all.gfr.runtime_2017.12.29, "microbenchmark")),
	resf = rbind(res, aves.all.gfr.runtime_2017.12.29),
	plt = microbenchmark:::autoplot.microbenchmark(resf)
)
config <- drake_config(plan)
vis_drake_graph(config)
make(plan)
config <- drake_config(plan)
vis_drake_graph(config)
outdated(config)
diagnose(drake_target_1)
