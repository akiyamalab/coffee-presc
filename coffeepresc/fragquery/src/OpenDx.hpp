#include "common.hpp"
#include "Fragment.hpp"
#include "FragmentInterEnergyGrid.hpp"
#include "Point3d.hpp"

#ifndef OPEN_DX_H_
#define OPEN_DX_H_
namespace opendx {
  class OpenDx {
  private:
    // fragdock::Fragment fragment;
    fragdock::FragmentInterEnergyGrid fgrid;
  public:
    OpenDx(const fragdock::FragmentInterEnergyGrid& fgrid) : fgrid(fgrid) {}
    void writeDx(std::ofstream& ofs) const;
    void write(const std::string& filename) const;
  };

  template<typename T>
  std::string point3d_to_string(const fragdock::Point3d<T>& p);
}
#endif