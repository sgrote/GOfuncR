#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>

struct gen_pos_str {
	std::string name; std::string chrom; long start; long end; long cumu_len;
};	
struct bed_str {
	std::string chrom; long start; long end; long len; long cumu_len;
};	
	
#endif
