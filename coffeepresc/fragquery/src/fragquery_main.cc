#include "common.hpp"
#include "main_utils.hpp"
#include "utils.hpp"
#include "OBMol.hpp"
#include "Point3d.hpp"
#include "infile_reader.hpp"
#include "log_writer_stream.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "FragmentInterEnergyGrid.hpp"
#include "EnergyCalculator.hpp"
#include "QueryGenerator.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <chrono>
#include <stdexcept>

namespace {
  format::DockingConfiguration parseArgs(int argc, char **argv){
    using namespace boost::program_options;
    options_description options("Options");
    options_description hidden("Hidden options");
    positional_options_description pos_desc;
    hidden.add_options()
      ("conf-file", value<std::string>(), "configuration file");
    pos_desc.add("conf-file", 1);
    options.add_options()
      ("help,h", "show help")
      ("output,o", value<std::string>(), "output file (.txt file)")
      ("receptor,r", value<std::string>(), "receptor file (.pdb file)")
      ("grid,g", value<std::string>(), "grid folder")
      ("log", value<std::string>(), "log file")
      ("fragments,f", value<std::string>(), "fragments file (.sdf file)")
      ("promising_pose", value<int64_t>(), "the number of promising pose each fragment (default: 40)")
      ("cluster_size", value<fltype>(), "cluster size (Å) of promising poses (default: 1.0)")
      ("distance_width", value<fltype>(), "discretization width (Å) of relative distance (default: 0.1)")
      ("distance_min", value<fltype>(), "minimum distance for query (default: 0.0)")
      ("distance_max", value<fltype>(), "maximum distance for query (default: 100.0)")
      ("single_fragment", value<bool>()->implicit_value(true)->default_value(false), "also create single fragment queries (distance is 0.0)")
      ("no_oberrorlog", value<bool>()->implicit_value(true)->default_value(false), "stop OBEroorLog")
      ("verbosity,v", value<int>()->implicit_value(1)->default_value(1), "verbosity level (0=errors only, 1=info+, 2=debug+)");
    options_description desc;
    desc.add(options).add(hidden);
    variables_map vmap;
    store(command_line_parser(argc, argv).
	  options(desc).positional(pos_desc).run(), vmap);
    notify(vmap);

    if (!vmap.count("conf-file") || vmap.count("help")){
      if (!vmap.count("conf-file") && !vmap.count("help")){
	std::cout << "too few arguments" << std::endl;
      }
      std::cout << "Usage: ligandock conf-file [options]\n"
		<< options << std::endl;
      std::exit(0);
    }
    format::DockingConfiguration conf = format::ParseInFile(vmap["conf-file"].as<std::string>().c_str());
    if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    if (vmap.count("fragments")) conf.fragments_file = vmap["fragments"].as<std::string>();
    if (vmap.count("output")) conf.output_file = vmap["output"].as<std::string>();
    if (vmap.count("grid")) conf.grid_folder = vmap["grid"].as<std::string>();
    if (vmap.count("log")) conf.log_file = vmap["log"].as<std::string>();
    if (vmap.count("promising_pose")) conf.query_params.promising_pose = vmap["promising_pose"].as<int64_t>();
    if (vmap.count("cluster_size")) conf.query_params.cluster_size = vmap["cluster_size"].as<fltype>();
    if (vmap.count("distance_width")) conf.query_params.distance_width = vmap["distance_width"].as<fltype>();
    if (vmap.count("distance_min")) conf.query_params.distance_min = vmap["distance_min"].as<fltype>();
    if (vmap.count("distance_max")) conf.query_params.distance_max = vmap["distance_max"].as<fltype>();
    conf.single_fragment = vmap["single_fragment"].as<bool>();
    conf.no_oberrorlog = vmap["no_oberrorlog"].as<bool>();
    if (vmap.count("verbosity")) conf.verbosity = vmap["verbosity"].as<int>();
    return conf;
  }

  void logConfig(const format::DockingConfiguration config){
    logs::lout << "Receptor file name  : "+config.receptor_file   << std::endl;
    logs::lout << "Fragments file name : "+config.fragments_file  << std::endl;
    logs::lout << "Output file name    : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name    : "+config.grid_folder     << std::endl;
  }

  struct query_record {
    std::string smi1, smi2;
    fltype dist_min, dist_max;
    int rank1, rank2;
    fltype score1, score2, pair_score;
    query_record(const std::string& smi1, const std::string& smi2, query::output_query oq) 
      : smi1(smi1), smi2(smi2), dist_min(oq.dist_min), dist_max(oq.dist_max), rank1(oq.rank1), rank2(oq.rank2), score1(oq.score1), score2(oq.score2), pair_score(oq.pair_score) {}
    query_record(const std::string& smi1, const std::string& smi2, fltype dist_min, fltype dist_max, int rank1, int rank2, fltype score1, fltype score2, fltype pair_score)
      : smi1(smi1), smi2(smi2), dist_min(dist_min), dist_max(dist_max), rank1(rank1), rank2(rank2), score1(score1), score2(score2), pair_score(pair_score) {}
    
    // 1st: pair_score, 2nd: smiles, 3rd: dist
    bool operator<(const query_record& o) {
      if (std::abs(pair_score - o.pair_score) < EPS) {
        if (smi1+smi2 == o.smi1+o.smi2) return dist_min < o.dist_min;
        else return smi1+smi2 < o.smi1+o.smi2;
      }
      return pair_score < o.pair_score;
    }
    static std::string getHeader() {
      return "f_1,f_2,start,end,rank_1,rank_2,score_1,score_2,pair_score";
    }
    friend std::ostream& operator<<(std::ostream& os, const query_record& qp);
  };

  std::ostream& operator<<(std::ostream& os, const query_record& qp) {
    os << qp.smi1 << "," << qp.smi2 << "," << qp.dist_min << "," << qp.dist_max << "," << qp.rank1 << "," << qp.rank2 << "," << qp.score1 << "," << qp.score2 << "," << qp.pair_score;
    return os;
  }

  void DebugQueryParams(const format::QueryParams& query_params) {
    logs::lout << logs::debug << "config.query_params.promising_pose: " << query_params.promising_pose  << std::endl;
    logs::lout << logs::debug << "config.query_params.cluster_size  : " << query_params.cluster_size    << std::endl;
    logs::lout << logs::debug << "config.query_params.distance_min  : " << query_params.distance_min    << std::endl;
    logs::lout << logs::debug << "config.query_params.distance_max  : " << query_params.distance_max    << std::endl;
    logs::lout << logs::debug << "config.query_params.distance_width: " << query_params.distance_width  << std::endl;
  }
} // namespace

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;
  using namespace main_utils;

  std::chrono::milliseconds whole_time(0);
  auto t0 = std::chrono::system_clock::now();

  format::DockingConfiguration config = parseArgs(argc, argv);

  if(config.log_file == ""){
    config.log_file = config.output_file + "__" + getDate() + ".log";
  }
  logs::log_init(config.log_file, config.verbosity);
  logConfig(config);

  if (config.no_oberrorlog) {
    OpenBabel::obErrorLog.StopLogging();
    logs::lout << "obErrorLog.StopLogging" << endl;
  }

  // parse receptor file
  const OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  const Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse fragments file
  vector<OpenBabel::OBMol> fragments = format::ParseFileToOBMol(config.fragments_file);

  int frags_sz = fragments.size(); /* the number of fragments */
  logs::lout << "number of fragments: " << frags_sz << endl;

  // ================================================================
  // prepare atomgrids and rotations
  // ================================================================
  logs::lout << logs::info << "[start] read energy grids" << endl;
  vector<AtomInterEnergyGrid> atom_grids = AtomInterEnergyGrid::readAtomGrids(config.grid_folder);
  logs::lout << logs::info << "[ end ] read energy grids" << endl;
  // logs::lout << "atom grid size: " << atom_grids.size() << endl;


  logs::lout << logs::info << "[TIME STAMP] START MOLECULE OBJECT CONVERSION" << endl;
   vector<Fragment> fragments_frag = convert_fragments(fragments);
  logs::lout << logs::info << "[TIME STAMP] END   MOLECULE OBJECT CONVERSION" << endl;


  InterEnergyGrid distance_grid = makeDistanceGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol);

  // ---------------------------------------

  DebugQueryParams(config.query_params);

  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;

  std::chrono::milliseconds query_time(0);
  std::chrono::milliseconds pose_time(0);
  int query_cnt = 0;

  auto t1 = std::chrono::system_clock::now();

  // vector of (position, score) of each fragment
  vector<vector<FragPose> > best_poses(frags_sz, vector<FragPose>());

  progress_bar grid_prog(frags_sz, "Searching fragment best poses", config.verbosity == 2);

  for (int f_ind = 0; f_ind < frags_sz; ++f_ind) {
    FragmentInterEnergyGrid frag_grid(fragments_frag[f_ind], makeRotations60(), atom_grids, distance_grid);
    best_poses[f_ind] = frag_grid.getBestPoses(config.grid.inner_width, config.query_params);
    grid_prog.display(f_ind+1);
  }
  grid_prog.clear();

  logs::lout << logs::info << "[TIME STAMP] END   CALCULATING BY FRAGGRID" << endl;


  auto t2 = std::chrono::system_clock::now();

  pose_time += std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1);
  logs::lout << logs::info << "[TIME STAMP] pose_time  : " << pose_time.count() << endl;

  std::chrono::milliseconds scoring_time(0);


  logs::lout << logs::info << "[TIME STAMP] START CALCULATING QUERY SCORES" << endl;

  vector<query_record> records;

  if (config.single_fragment) {
    logs::lout << logs::debug << "Single fragment queries: True" << endl;
    for (int i = 0; i < frags_sz; ++i) {
      // pair_score is the best score of the fragment
      fltype score = best_poses[i][0].score;
      records.push_back(query_record(fragments_frag[i].getSmiles(),
                                    fragments_frag[i].getSmiles(),
                                    0.0, 0.0, // distance is 0.0
                                    -1, -1,
                                    0, 0,
                                    score));
      query_cnt++;
    }
  }

  progress_bar output_prog(frags_sz*frags_sz, "Calculating query scores", config.verbosity == 2);
  for (int i = 0; i < frags_sz; ++i) {
    for (int j = i; j < frags_sz; ++j) {
      const vector<FragPose>& poses1 = best_poses[i];
      const vector<FragPose>& poses2 = best_poses[j];
      query::QueryGenerator query_generator(config.query_params);
      for (int k = 0; k < poses1.size(); ++k) {
        for (int l = (i==j)?k+1:0; l < poses2.size(); ++l) {
          // distance between two fragments
          fltype distance = (poses1[k].coord - poses2[l].coord).abs();
          fltype pair_score = poses1[k].score + poses2[l].score;
          query_generator.append(distance, k, l, poses1[k].score, poses2[l].score, pair_score);
        }
      }
      vector<query::output_query> oqueries = query_generator.getPairInfoVec();
      for (int d = 0; d < oqueries.size(); ++d) {
        records.push_back(query_record(fragments_frag[i].getSmiles(),
                                      fragments_frag[j].getSmiles(),
                                      oqueries[d]));
        query_cnt++;
      }
      output_prog.display(i*frags_sz+j+1);
    }
  }
  output_prog.clear();

  std::sort(records.begin(), records.end());
  ofstream outputcsv(config.output_file);
  int query_size = (NUM_QUERIES == -1) ? records.size() : min(NUM_QUERIES, (int)records.size());
  outputcsv << query_record::getHeader() << endl;
  for (int i = 0; i < query_size; ++i) {
    outputcsv << records[i] << endl;
  }
  outputcsv.close();

  auto t3 = std::chrono::system_clock::now();

  scoring_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t2);
  query_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t1);

  logs::lout << logs::info << "[TIME STAMP] END   CALCULATING QUERY SCORES" << endl;
  logs::lout << logs::info << "[TIME STAMP] scoring_time  : " << scoring_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] query_time : " << query_time.count() << endl;
  logs::lout << logs::debug << "query_cnt  : " << query_cnt << endl;

  auto t4 = std::chrono::system_clock::now();
  whole_time += std::chrono::duration_cast< std::chrono::milliseconds >(t4 - t0);

  logs::lout << logs::info << "################ Program end ################" << endl;
  logs::lout << logs::info << "[TIME STAMP] whole_time : " << whole_time.count() << endl;

  logs::close();
  return 0;
}