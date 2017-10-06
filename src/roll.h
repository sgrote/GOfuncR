#ifndef ROLL_H
#define ROLL_H

#include <set>
#include <vector>
#include <map>
#include "structures.h"

std::set<int> rannum_roll(std::vector<bed_str> candidate_bed, std::vector<bed_str> background_bed, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos);

#endif
