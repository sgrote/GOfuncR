#ifndef GENELEN_H
#define GENELEN_H

#include <set>
#include <vector>
#include <map>
#include "structures.h"

std::set<int> rannum_genelen(int n_candidate, const std::map<std::string,int> &genename_to_index, std::vector<gen_pos_str> genes_pos, long total_length);

#endif
