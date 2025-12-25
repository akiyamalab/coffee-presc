#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "infile_reader.hpp"


namespace format {
  DockingConfiguration ParseInFile(const char *filename){
    std::ifstream ifs(filename);
    if (ifs.fail()){
      std::cerr << "opening grid file failed:" << filename << std::endl;
      abort();
    }
    DockingConfiguration conf;
    std::string buffer;
    while(!ifs.eof()){
      std::getline(ifs, buffer);
      try {
        if (boost::algorithm::starts_with(buffer, "OUTERBOX ")) {
          std::string str = buffer.substr(9);
          std::vector<std::string> vals;
          boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
          boost::algorithm::trim(vals[0]);
          boost::algorithm::trim(vals[1]);
          boost::algorithm::trim(vals[2]);
          conf.grid.outer_width = fragdock::Point3d<fltype>(
            boost::lexical_cast<fltype>(vals[0]),
            boost::lexical_cast<fltype>(vals[1]),
            boost::lexical_cast<fltype>(vals[2]));
        }
        else if (boost::algorithm::starts_with(buffer, "INNERBOX ")) {
          std::string str = buffer.substr(9);
          std::vector<std::string> vals;
          boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
          boost::algorithm::trim(vals[0]);
          boost::algorithm::trim(vals[1]);
          boost::algorithm::trim(vals[2]);
          conf.grid.inner_width = fragdock::Point3d<fltype>(
            boost::lexical_cast<fltype>(vals[0]),
            boost::lexical_cast<fltype>(vals[1]),
            boost::lexical_cast<fltype>(vals[2]));
        }
        else if (boost::algorithm::starts_with(buffer, "BOX_CENTER ")) {
          std::string str = buffer.substr(11);
          std::vector<std::string> vals;
          boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
          boost::algorithm::trim(vals[0]);
          boost::algorithm::trim(vals[1]);
          boost::algorithm::trim(vals[2]);
          conf.grid.center = fragdock::Point3d<fltype>(
            boost::lexical_cast<fltype>(vals[0]),
            boost::lexical_cast<fltype>(vals[1]),
            boost::lexical_cast<fltype>(vals[2]));
        }
        else if (boost::algorithm::starts_with(buffer, "SCORING_PITCH ")) {
          std::string str = buffer.substr(14);
          std::vector<std::string> vals;
          boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
          boost::algorithm::trim(vals[0]);
          boost::algorithm::trim(vals[1]);
          boost::algorithm::trim(vals[2]);
          conf.grid.score_pitch = fragdock::Point3d<fltype>(
            boost::lexical_cast<fltype>(vals[0]),
            boost::lexical_cast<fltype>(vals[1]),
            boost::lexical_cast<fltype>(vals[2]));
        }
        else if (boost::algorithm::starts_with(buffer, "RECEPTOR ")) {
          conf.receptor_file = buffer.substr(9);
        }
        else if (boost::algorithm::starts_with(buffer, "OUTPUT ")) {
          conf.output_file = buffer.substr(7);
        }
        else if (boost::algorithm::starts_with(buffer, "LOG ")) {
          conf.log_file = buffer.substr(4);
        }
        else if (boost::algorithm::starts_with(buffer, "GRID_FOLDER ")) {
          conf.grid_folder = buffer.substr(12);
        }
        else if (boost::algorithm::starts_with(buffer, "ROTANGS ")) {
          conf.rotangs_file = buffer.substr(8);
        }
        else if (boost::algorithm::starts_with(buffer, "FRAGMENTS ")) {
          // FRAGMENTS [FRAGMENTS]
          conf.fragments_file = buffer.substr(10);
        }
        else if (boost::algorithm::starts_with(buffer, "DX_FOLDER ")) {
          conf.dx_folder = buffer.substr(10);
        }
        else if (boost::algorithm::starts_with(buffer, "PROMISING_POSE ")) {
          std::string str = buffer.substr(15);
          boost::algorithm::trim(str);
          conf.query_params.promising_pose = boost::lexical_cast<int64_t>(str);
        }
        else if (boost::algorithm::starts_with(buffer, "CLUSTER_SIZE ")) {
          std::string str = buffer.substr(13);
          boost::algorithm::trim(str);
          conf.query_params.cluster_size = boost::lexical_cast<fltype>(str);
        }
        else if (boost::algorithm::starts_with(buffer, "DISTANCE_WIDTH ")) {
          std::string str = buffer.substr(15);
          boost::algorithm::trim(str);
          conf.query_params.distance_width = boost::lexical_cast<fltype>(str);
        }
        else if (boost::algorithm::starts_with(buffer, "DISTANCE_MIN ")) {
          std::string str = buffer.substr(13);
          boost::algorithm::trim(str);
          conf.query_params.distance_min = boost::lexical_cast<fltype>(str);
        }
        else if (boost::algorithm::starts_with(buffer, "DISTANCE_MAX ")) {
          std::string str = buffer.substr(13);
          boost::algorithm::trim(str);
          conf.query_params.distance_max = boost::lexical_cast<fltype>(str);
        }
        else if (boost::algorithm::starts_with(buffer, "SINGLE_FRAGMENT ")) {
          std::string str = buffer.substr(16);
          boost::algorithm::trim(str);
          if (str == "true" || str == "1") conf.single_fragment = true;
          else conf.single_fragment = false;
        }
        else if (boost::algorithm::starts_with(buffer, "NO_OBERRORLOG ")) {
          std::string str = buffer.substr(14);
          boost::algorithm::trim(str);
          if (str == "true" || str == "1") conf.no_oberrorlog = true;
          else conf.no_oberrorlog = false;
        }
        else if (boost::algorithm::starts_with(buffer, "VERBOSITY ")) {
          std::string str = buffer.substr(10);
          boost::algorithm::trim(str);
          conf.verbosity = boost::lexical_cast<int>(str);
        }
      } catch (const boost::bad_lexical_cast &e) {
        std::cerr << "DockingConfiguration ParseInFile(): lexical_cast error for line: " << buffer << std::endl;
        abort();
      }
    }

    return conf;
  }
} // namespace format