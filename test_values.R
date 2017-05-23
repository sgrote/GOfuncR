
# check that results match the results created previously or store new results
# needs to be executed in /FuncBlocks_package/   (to use on any computer)

Sys.setlocale("LC_COLLATE","C")
library(FuncBlocks)
setwd("./tests")
source("test_parallel.R")  # this does check on the fly when sourced
source("test_go_enrich_values.R")
source("test_get_functions_values.R")

# one optional Command Line argument "set_values" if not, then default *test values*
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
	set_values = FALSE  # set values (or test values)
	if (!("saved_results.RData") %in% dir()){
		stop("saved_results.RData not found")
	}
} else if (args[1] == "set_values"){
	set_values = TRUE
} else {
	stop ("Please use 'set_values' or nothing as command line argument (nothing = test-values)")
}

### compute values
# get go_enrich test-set values (lists of 3)
go_enrich_results = go_enrich_values()
# get get-functions test-set values (dataframes or vectors)
get_results = get_values()

#### save or check values
test_results = c(go_enrich_results, get_results)
if (set_values){
	saved_results = test_results
	session_info = sessionInfo()
	save(list=c("saved_results", "session_info"), file="saved_results.RData")
	message(paste("\nWrote saved_results.RData to", getwd()))
} else {
	message("checking results...")
	load("saved_results.RData")
	failed = 0
	for (i in 1:length(test_results)){
		if (!(isTRUE(all.equal(saved_results[[i]], test_results[[i]])))){
			failed = failed+1
			message(paste("\ntest_results[[",i,"]] is failing (", names(saved_results)[i],").",sep=""))
			if(is.data.frame(test_results[[i]])){
				message("saved_results")
				print(head(saved_results[[i]]))
				message("test_results")
				print(head(test_results[[i]]))
			} else { # assuming it's a go_enrich list then
				# NEW: loop again: print only those that dont match from res[[1]],res[[2]],res[[3]]
				for (j in 1:length(saved_results[[i]])){
						if(!(isTRUE(all.equal(saved_results[[i]][[j]], test_results[[i]][[j]])))){
						message(paste("saved_results[[",j,"]]",sep=""))
						print(head(saved_results[[i]][[j]]))
						message(paste("test_results[[",j,"]]",sep=""))
						print(head(test_results[[i]][[j]]))
					}
				}
			}	
		}
	}
	message(paste("\nNumer of failed tests:",failed))		
}

