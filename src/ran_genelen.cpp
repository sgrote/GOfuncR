
// draw unique integers [1,total length of genes] and select genes dependent on gene length

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include "structures.h"

#include <Rcpp.h>

std::set<int> rannum_genelen(int n_candidate, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos, long total_length){
	
	//Rcpp::Rcout << std::endl << "ran_genelen option:" << std::endl << std::endl;
	
	// get random genes and their index, without duplicates
	std::set<int> random_numbers;
	// draw as many unique numbers as candidate genes
	while (random_numbers.size() < n_candidate) {	
		// get random number [1,total length of genes]
		long ran = R::runif(0,1) * (total_length) + 1; 
		// select gene
		int k = 0;
		while (ran > genes_pos[k].cumu_len){
			k++;
		}	
		random_numbers.insert(genename_to_index.find(genes_pos[k].name)->second);
	}	
	//Rcpp::Rcout << "Candidate genes: " << n_candidate << ", Random genes: " << random_numbers.size() << std::endl;
	return(random_numbers);
}
