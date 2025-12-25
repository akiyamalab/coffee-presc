#include "common.hpp"
#include "InterEnergyGrid.hpp"
// #include "EnergyCalculator.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "Fragment.hpp"
#include "Point3d.hpp"
#include "infile_reader.hpp"

#ifndef FRAGMENT_ENERGY_GRID_H_
#define FRAGMENT_ENERGY_GRID_H_
namespace fragdock {
  struct FragPose {
    Vector3d coord;
    fltype score;
    FragPose() {}
    FragPose(const Vector3d& coord, fltype score) : coord(coord), score(score) {}
  };

  class FragmentInterEnergyGrid /* : public InterEnergyGrid */ {
    InterEnergyGrid grid;
    // void parse(const std::string& filename, int rot_size);
  public:
    int frag_idx;
    FragmentInterEnergyGrid() { frag_idx = -1; }
    // FragmentInterEnergyGrid(int frag_id, const std::string& filename, int rot_size)
    // : frag_id(frag_id) { parse(filename, rot_size); }
    FragmentInterEnergyGrid(const Fragment& orig_frag,
                       const std::vector<Vector3d>& rot_angles,
                       const std::vector<AtomInterEnergyGrid>& atom_grids,
                       const InterEnergyGrid& distance_grid);
    ~FragmentInterEnergyGrid() {}
    const InterEnergyGrid& getGrid() const { return grid; }
    // void writeFile(const std::string& filename) const;    
    fltype getSearchScore(int x, int y, int z) const { return getGrid().getInterEnergy(x, y, z); }
    std::vector<FragPose> getBestPoses(const Point3d<fltype>& search_width, const format::QueryParams& query_params) const;
  };
}
#endif

