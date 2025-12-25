#include "OpenDx.hpp"

namespace opendx {
  template<typename T>
  std::string point3d_to_string(const fragdock::Point3d<T>& p) {
    using namespace std;
    return to_string(p.x) + " " + to_string(p.y) + " " + to_string(p.z);
  }

  void OpenDx::writeDx(std::ofstream& ofs) const {
    using namespace std;

    fragdock::Point3d<int> num = fgrid.getGrid().getNum();
    fragdock::Point3d<fltype> center = fgrid.getGrid().getCenter();
    fragdock::Point3d<fltype> pitch = fgrid.getGrid().getPitch();
    fragdock::Point3d<fltype> origin = center - pitch * num / 2;

    ofs << "# Data from Compound Selection" << endl;
    ofs << "# Score: -energy" << endl;
    ofs << "# !! Energy is reciprocal !! " << endl;
    
    ofs << "object 1 class gridpositions counts " << point3d_to_string<int>(num) << endl;
    ofs << "origin " << point3d_to_string<fltype>(origin) << endl;
    ofs << "delta " << pitch.x << " 0 0" << endl;
    ofs << "delta 0 " << pitch.y << " 0" << endl;
    ofs << "delta 0 0 " << pitch.z << endl;
    ofs << "object 2 class gridconnections counts " << point3d_to_string<int>(num) << endl;
    ofs << "object 3 class array type \"double\" rank 0 items " << num.x*num.y*num.z << " data follows" << endl;
    
    for (int x = 0; x < num.x; x++) {
      for (int y = 0; y < num.y; y++) {
        for (int z = 0; z < num.z; z++) {
          ofs << -fgrid.getGrid().getInterEnergy(x, y, z) << " ";
        }
        ofs << endl;
      }
    }

    ofs << "attribute \"dep\" string \"positions\"" << endl;
    ofs << "object \"positions\" class field" << endl;
    ofs << "component \"positions\" value 1" << endl;
    ofs << "component \"connections\" value 2" << endl;
    ofs << "component \"data\" value 3" << endl;
    ofs << "end" << endl;
  }
  void OpenDx::write(const std::string& filename) const {
    using namespace std;

    ofstream ofs(filename.c_str());
    if(!ofs){
      cerr << "OpenDx::write() : file could not open. " << filename << endl;
      return ;
    }
    writeDx(ofs);
    ofs.close();
  }
}