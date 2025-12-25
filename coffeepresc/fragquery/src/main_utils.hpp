#include "common.hpp"
#include "Point3d.hpp"
#include "Vector3d.hpp"
#include "InterEnergyGrid.hpp"
#include "Molecule.hpp"
#include "Fragment.hpp"
#include "OBMol.hpp"

#include <cmath>

#ifndef MAIN_UTILS_H_
#define MAIN_UTILS_H_

namespace main_utils {
	void makeFolder(std::string folderName);
	std::string getDate();
	fragdock::Point3d<int> to_score_num(int k, 
                                      const fragdock::Point3d<int>& score_num, 
                                      const fragdock::Point3d<int>& search_num, 
                                      const fragdock::Point3d<int>& ratio);
	fragdock::Vector3d operator/(fragdock::Vector3d v, fragdock::Point3d<fltype> p);
	fragdock::Point3d<int> round(const fragdock::Vector3d& v);
	// not restretto
	std::vector<fragdock::Fragment> convert_fragments(std::vector<OpenBabel::OBMol>& obmols);
	std::string option_desc(const std::string main, const std::vector<std::pair<std::string, std::string> >& options);
	
  struct progress_bar {
    const int bar_size = 100;
    const int total;
    const std::string desc;
    const bool verbose = false; // default is not to show progress bar
    progress_bar(int total, bool verbose) : total(total), desc(""), verbose(verbose) { display(0); }
    progress_bar(int total, std::string desc, bool verbose) : total(total), desc(desc+":"), verbose(verbose) { display(0); }
    void display(int progress) const;
    void clear() const;
  };
} // namespace main_utils

namespace fragdock {
	InterEnergyGrid makeDistanceGrid(const Point3d<fltype>& center,
                                   const Point3d<fltype>& pitch,
                                   const Point3d<int>& num,
                                   const Molecule& receptor_mol);
} // namespace fragdock

#endif