#include "common.hpp"
#include "InterEnergyGrid.hpp"
#include <string>
#include <vector>

#ifndef INFILE_READER_H_
#define INFILE_READER_H_

namespace format {
  struct SearchGrid {
    fragdock::Point3d<fltype> center, outer_width, inner_width, search_pitch, score_pitch;
  };

  struct QueryParams {
    int64_t promising_pose = 40;
    fltype cluster_size = 1.0;
    fltype distance_min = 0.0;
    fltype distance_max = 100.0;
    fltype distance_width = 0.1;
  };

  struct DockingConfiguration {
    SearchGrid grid;
    QueryParams query_params;
    std::vector<std::string> ligand_files;
    std::string receptor_file, output_file;
    std::string log_file, grid_folder, dx_folder;
    std::string rotangs_file;
    std::string fragments_file;
    bool single_fragment;
    bool no_oberrorlog;
    int verbosity;

    void checkConfigValidity() const;
  };

  DockingConfiguration ParseInFile(const char *filename);
} // namespace format

#endif