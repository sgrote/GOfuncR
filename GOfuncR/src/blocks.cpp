
// given input and background regions choose random regions on any chromosome
// find genes in these regions and add their index to random set

#include <set>
//#include <stdlib.h>     // srand, rand 
//#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#include "structures.h"

#include <Rcpp.h>

// take input regions, background_regions and number of randomsets as input

std::set<int> rannum_blocks(std::vector<bed_str> candidate_bed, std::vector<bed_str> background_bed, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos){
	
	//Rcpp::Rcout << std::endl << "Blocks option:" << std::endl;

	// create randomset  
	//srand (time(NULL)); // initialize random seed
	std::vector<bed_str> background; // copy of background_bed to modify
	std::set<int> random_numbers; // indices of selected genes	
	int sum_genes = 0;

	// loop over candidate regions
	bool incomplete = true;
	int trials = 0; //for placing all candidate regions without overlap
	while (incomplete){
		trials++ ; 
		background = background_bed; // leave background_bed unmodified for future trials
		for (int j=0; j < candidate_bed.size(); j++){
			
			// get cumulative length for each region (without unstartable ends)			
			long total_len = 0;
			for (int k=0; k < background.size(); k++){  
				int this_len = background[k].len - candidate_bed[j].len;
				if (this_len < 0) this_len = 0;
				total_len += this_len;
				background[k].cumu_len = total_len;
			}	
			//Rcpp::Rcout << std::endl << "Residual background regions for candidate_region " << j+1 << std::endl;
			//for (int i=0; i < background.size(); i++){
				//Rcpp::Rcout << background[i].chrom << " " << background[i].start << " " << background[i].end  << " " << background[i].len << " " << background[i].cumu_len << std::endl;
			//}
			// candidate region does not fit anymore, try again placing randomly
			if(total_len == 0){
				incomplete = true;
				Rcpp::Rcout << "The candidate does not fit - try again..." << std::endl;
				Rcpp::Rcout << "This was trial " << trials << std::endl;
				if(trials == 10){
					Rcpp::Rcout << "Error: " << trials << " times in a row the candidate regions could not be placed randomly without forcing them to overlap. Consider using larger background regions." << std::endl;
					Rcpp::stop("Background regions too small."); 
				}
				break;		
			} else {
				incomplete = false;		
			}
			// choose random number [1, total length] 
			// int runif(0,1)*10 = [0,9]
			long ran = R::runif(0,1) * (total_len) + 1; 
			//long ran = rand() % total_len + 1; 
			long last_cumu = 0;
			int k = 0;
			// choose chromosome and region
			while (ran > background[k].cumu_len){  
				last_cumu = background[k].cumu_len;
				k++;
			}		
			long ran_start = background[k].start + (ran - last_cumu);			
			long ran_end = ran_start + candidate_bed[j].len;
			std::string ran_chrom = background[k].chrom; 
			
			//Rcpp::Rcout << std::endl << "Choosen random region " << j+1 << ":" << std::endl;
			//Rcpp::Rcout << ran_chrom << " " << ran_start << " " << ran_end << " " << ran_end-ran_start << std::endl;	

			// go through genes positions and select those that overlap randomly chosen region			
			for (int g=0; g<genes_pos.size(); g++){
				if(genes_pos[g].chrom == ran_chrom &&
				((genes_pos[g].start >= ran_start && genes_pos[g].start < ran_end) ||
				(genes_pos[g].end >= ran_start && genes_pos[g].end < ran_end) ||
				(genes_pos[g].start <= ran_start && genes_pos[g].end >= ran_end))){
					// add to set of randomly chosen test genes 
					// get index for every random gene name and add to random numbers
					random_numbers.insert(genename_to_index.find(genes_pos[g].name)->second);  
					sum_genes ++;
					//Rcpp::Rcout << genes_pos[g].name << " " << genes_pos[g].chrom << " " << genes_pos[g].start << " " << genes_pos[g].end << " " << genename_to_index.find(genes_pos[g].name)->second << std::endl;			
				}
			}
			// update background regions to avoid overlapping random regions
			// split chosen background region into two parts
			// second part 
			bed_str tail;
			tail.chrom = background[k].chrom;
			tail.start = ran_end;
			tail.end = background[k].end;
			tail.len = tail.end - tail.start;
			background.push_back(tail);
			// first part
			background[k].end = ran_start;
			background[k].len = ran_start - background[k].start;					
		} // end candidate regions	
	} // end while incomplete
		
	//Rcpp::Rcout << std::endl << "sum of random genes: " << sum_genes << std::endl;
	//Rcpp::Rcout << "sum of unique random genes: " << random_numbers.size() << std::endl;

	return(random_numbers);		
}
