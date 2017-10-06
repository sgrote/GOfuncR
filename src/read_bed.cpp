
// Read bedfile and store in vector of regions
	
#include "structures.h"
#include <fstream>
#include <sstream>
#include <vector>


std::vector< bed_str > read_bed(std::string bed_file){
	std::vector< bed_str > region_bed;
	std::ifstream region (bed_file.c_str());
	std::string line;
	while(std::getline( region, line )){
		bed_str bed;
		std::istringstream is( line.c_str()) ;
		is >> bed.chrom >> bed.start >> bed.end;	
		bed.len = bed.end - bed.start; // lengths of regions			
		region_bed.push_back(bed);  
	}
	region.close();
	return(region_bed);
}
