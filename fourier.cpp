#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>
#include "fourier.h"
using namespace std;
namespace SubstringMatching{
	std::vector<size_t> FindSubstrings(const std::string& str, const std::string& pattern){
		std::vector<size_t> res;
		for (size_t i = 0; i < str.size(); i++){
			if (!pattern.compare(str.substr(i, pattern.size()))){
				res.push_back(i);
			} 
		}
		return res;
	}
}
