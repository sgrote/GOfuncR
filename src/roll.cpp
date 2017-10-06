
// given input and background regions choose random regions on the same chromosome like input regions
// find genes in these regions and add their index to random set

// removing used parts from background regions
// case a) candidate starts in one and ends in the second bg-region 
// case b) like (a) with blocks inbetween which are used completely
// case c) canidate region fits in just one block 
// case d) there is just one bg-block, candidate region starts at end and ends at beginning (handled like case a) 
	
#include <set>
//#include <stdlib.h>     // srand, rand 
//#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#include "structures.h"

#include <Rcpp.h>

// take input regions, background_regions, gene_name-to-index-map and genes positions as arguments 

std::set<int> rannum_roll(std::vector<bed_str> candidate_bed, std::vector<bed_str> background_bed, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos){
	
	//Rcpp::Rcout << std::endl << "Circ_chrom option:" << std::endl;
	
	//srand (time(NULL)); // initialize random seed
	std::set<int> random_numbers; // indices of selected genes	
	int sum_genes = 0;

	// loop over candidate regions
	for (int j=0; j < candidate_bed.size(); j++){		
		std::string candidate_chrom = candidate_bed[j].chrom;
		// candidate.cumulen: length of background regions on same chromosome	
		// go through mappable regions, sum up length from chrom and store last value from chromosome
		// last value from chrom gets assigned	
		long chrom_len = 0;
		for (int i=0; i < background_bed.size(); i++){
			if (background_bed[i].chrom == candidate_chrom){
				chrom_len += background_bed[i].len;
				background_bed[i].cumu_len = chrom_len;
				candidate_bed[j].cumu_len = chrom_len;
			} 
		}
		//Rcpp::Rcout << std::endl << "Candidate region " << j+1 << " modified: " << std::endl;
		//Rcpp::Rcout << candidate_bed[j].chrom << " " << candidate_bed[j].start << " " << candidate_bed[j].end  << " " << candidate_bed[j].len << " " << candidate_bed[j].cumu_len << std::endl;
		//Rcpp::Rcout << std::endl << "Background regions modified:" << std::endl;
		//for (int i=0; i < background_bed.size(); i++){
			//Rcpp::Rcout << background_bed[i].chrom << " " << background_bed[i].start << " " << background_bed[i].end  << " " << background_bed[i].len << " " << background_bed[i].cumu_len << std::endl;
		//}	
		// find chromosome in mappable regions
		int k = 0;
		while (background_bed[k].chrom != candidate_chrom){  
			 k++;
		}
		int chrom_start = k; // is needed below to go back to beginning of chrom when rolling
		//Rcpp::Rcout << std::endl << "Start of background region for candidate region " << j+1 << ":" << std::endl;
		//Rcpp::Rcout << background_bed[k].chrom << " " << background_bed[k].start << " " << background_bed[k].end  << " " << background_bed[k].len << " " << background_bed[k].cumu_len << std::endl;
	
		// choose random number [1, total length of chrom] 
		// int runif(0,1)*10 = [0,9]
		long ran = R::runif(0,1) * candidate_bed[j].cumu_len + 1; 
		//long ran = rand() % candidate_bed[j].cumu_len + 1; 
		long last_cumu = 0;
		// choose region (k is first index of background on the current chromosome)
		while (ran > background_bed[k].cumu_len){  
			last_cumu = background_bed[k].cumu_len;
			k++;
		}		
		long ran_start = background_bed[k].start + (ran - last_cumu);			
		long ran_end = ran_start + candidate_bed[j].len;
		std::string ran_chrom = background_bed[k].chrom; 				
		
		// let it roll!
		//Rcpp::Rcout << "Chosen random regions " << j+1 << " of length " << candidate_bed[j].len << ":" << std::endl;
		long overhang = 0;
		bool need_more = true;	
		int used_bg_blocks = 0; //number of background regions used, also 2 if only 1 block is started again
		while (true){
			used_bg_blocks ++;
			// if current bg-block is not long enough:
			if (ran_end > background_bed[k].end){
				overhang = ran_end - background_bed[k].end;
				ran_end = background_bed[k].end;
				// shorten this block if it is the first one (cases a,b,d)
				if (used_bg_blocks == 1){
					background_bed[k].end = ran_start;
					background_bed[k].len = background_bed[k].end - background_bed[k].start;
				} else { // delete this block (used completely, case b)					
					background_bed.erase(background_bed.begin() + k);
					k--; // account for deleted object					
				}									
			} else {
				need_more = false;
				overhang = 0;
				// shorten if this block is not the first (cases a,b,d) (else it needs to be split in two)
				if (used_bg_blocks != 1){
					background_bed[k].start = ran_end;
					background_bed[k].len = background_bed[k].end - background_bed[k].start;
				}
			}	
			//Rcpp::Rcout << "need more: " << need_more << ", k: " << k << ", chrom: " << ran_chrom << ", start: " << ran_start << ", end: " << ran_end << ", overhang: " << overhang << ", used bg_blocks: " << used_bg_blocks << std::endl;
							
			// go through genes positions and select those that overlap randomly chosen region			
			for (int g=0; g<genes_pos.size(); g++){
				if(genes_pos[g].chrom == ran_chrom &&
				((genes_pos[g].start >= ran_start && genes_pos[g].start < ran_end) ||
				(genes_pos[g].end >= ran_start && genes_pos[g].end < ran_end) ||
				(genes_pos[g].start <= ran_start && genes_pos[g].end >= ran_end))){
					// add to set of randomly chosen test genes	
					random_numbers.insert(genename_to_index.find(genes_pos[g].name)->second);  
					sum_genes ++;
					//Rcpp::Rcout << genes_pos[g].name << " " << genes_pos[g].chrom << " " << genes_pos[g].start << " " << genes_pos[g].end << " " << genename_to_index.find(genes_pos[g].name)->second << std::endl;
				}				
			}
			if (!need_more) break;				
			// next block 
			k ++;
			// roll chromosome
			if ((k == background_bed.size()) || (background_bed[k].chrom != ran_chrom)){
				k = chrom_start;
			}
			ran_start = background_bed[k].start;
			ran_end = ran_start + overhang;				
		}
		// update background_bed, case (c): only one block was used, split in two and insert tail after head
		if (used_bg_blocks == 1){
			// tail	
			bed_str tail;
			tail.chrom = background_bed[k].chrom;
			tail.start = ran_end;
			tail.end = background_bed[k].end;
			tail.len = tail.end - tail.start;
			background_bed.insert(background_bed.begin()+k+1, tail); // insert after index of current block	
			// head
			background_bed[k].end = ran_start;
			background_bed[k].len = background_bed[k].end - background_bed[k].start;
		}
					
	} // end candidate regions	
	
	//Rcpp::Rcout << std::endl << "sum of random genes: " << sum_genes << std::endl;
	//Rcpp::Rcout << "sum of unique random genes: " << random_numbers.size() << std::endl;

	return(random_numbers);		
}
