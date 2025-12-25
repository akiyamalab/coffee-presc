#include "main_utils.hpp"

#include <sys/stat.h>

namespace main_utils {
	void makeFolder(std::string folderName){
    struct stat st;
    if(stat(folderName.c_str(), &st)==-1){
      std::cout << "there isn't folder. Try create: " << folderName << std::endl;
      if(mkdir(folderName.c_str(), 0775) == -1){
	std::cerr << "Fail to make folder: " << folderName << std::endl;
	abort();
      }
    }else if(!S_ISDIR(st.st_mode)){
      std::cerr << "Fail to make folder. There is same name file (not folder): " << folderName << std::endl;
      abort();
    }
  }

	std::string getDate() {
    time_t timer = time(NULL);
    tm* date = localtime(&timer);
    std::string dateStr = (boost::format("%02d_%02d_%02d_%02d_%02d")
			   % date->tm_mon % date->tm_mday
			   % date->tm_hour % date->tm_min % date-> tm_sec).str();
    return dateStr;
  }

	fragdock::Point3d<int> to_score_num(int k, const fragdock::Point3d<int>& score_num, const fragdock::Point3d<int>& search_num, const fragdock::Point3d<int>& ratio) {
    return score_num / 2 + (-search_num / 2 + k) * ratio;
  }

	fragdock::Vector3d operator/(fragdock::Vector3d v, fragdock::Point3d<fltype> p){
    return fragdock::Vector3d(v.x/p.x, v.y/p.y, v.z/p.z);
  }
	
  fragdock::Point3d<int> round(const fragdock::Vector3d& v) {
    return fragdock::Point3d<int>(std::round(v.x), std::round(v.y), std::round(v.z));
  }

	std::vector<fragdock::Fragment> convert_fragments(std::vector<OpenBabel::OBMol>& obmols) {
    std::vector<fragdock::Fragment> fragments_frag(obmols.size()); /* a vector of fragment objects */
    for (uint f_ind = 0; f_ind < obmols.size(); ++f_ind) {
      OpenBabel::OBMol& ob_fragment = obmols[f_ind];
      ob_fragment.AddHydrogens();
      fragments_frag[f_ind] = format::toFragment(ob_fragment, f_ind);
      ob_fragment.DeleteHydrogens();
      fragments_frag[f_ind].translate(-fragments_frag[f_ind].getCenter());
      fragments_frag[f_ind].setSmiles(OpenBabel::canonicalSmiles(ob_fragment));
    }
    return fragments_frag;
  }

  std::string option_desc(const std::string main, const std::vector<std::pair<std::string, std::string> >& options) {
    // align options
    size_t max_size = 0;
    for (auto& option : options) {
      max_size = std::max(max_size, option.first.size());
    }
    std::string desc = main + "\n";
    for (auto& option : options) {
      desc += "  " + option.first;
      for (size_t i = 0; i < max_size - option.first.size(); ++i) {
        desc += " ";
      }
      desc += ": " + option.second + "\n";
    }
    // last "\n" is not needed
    desc.pop_back();
    return desc;
  }

  void progress_bar::display(int progress) const {
    if (!verbose) return;
    float fraction = static_cast<float>(progress) / total;
    int progress_size = static_cast<int>(fraction * bar_size);
    std::string frac_per = std::to_string(static_cast<int>(fraction * 100));

    std::cout << desc;
    for (int i = 0; i < 3 - frac_per.size(); ++i) {
      std::cout << " ";
    }
    std::cout << frac_per << "%[";
    for (int i = 0; i < bar_size; ++i) {
      if (i < progress_size)
        std::cout << "#";
      else
        std::cout << ".";
    }
    std::cout << "] " << progress << "/" << total << "\r";
    std::cout.flush();
  }

  void progress_bar::clear() const {
    if (!verbose) return;
    std::cout << "\r";
    int del_size = desc.size() +  bar_size + std::to_string(total).size()*2 + 10;
    for (int i = 0; i < del_size; ++i) {
      std::cout << " ";
    }
    std::cout << "\r";
    std::cout.flush();
  }
} // namespace main_utils

namespace fragdock {
  InterEnergyGrid makeDistanceGrid(const Point3d<fltype>& center,
                              const Point3d<fltype>& pitch,
                              const Point3d<int>& num,
                              const Molecule& receptor_mol) {

    InterEnergyGrid grid(center, pitch, num);
    for(int x=0; x<grid.getNum().x; x++) {
      for(int y=0; y<grid.getNum().y; y++) {
        for(int z=0; z<grid.getNum().z; z++) {

          fltype mindist = 1e4;
          Vector3d pos = grid.convert(x, y, z);

          for (const Atom& a : receptor_mol.getAtoms()) {
            if (a.getXSType() != XS_TYPE_H) {
              mindist = std::min(mindist, (pos - a).abs());
            }
          }
          grid.setInterEnergy(x, y, z, mindist);
        }
      }
    }
    return grid;
  }

} // namespace fragdock