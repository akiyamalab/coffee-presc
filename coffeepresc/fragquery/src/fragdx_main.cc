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
#include "OpenDx.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <chrono>
#include <stdexcept>

#include <unordered_map>
#include <thread>

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
      ("output,o", value<std::string>(), "output file (.sh file)")
      ("receptor,r", value<std::string>(), "receptor file (.pdb file)")
      ("grid,g", value<std::string>(), "grid folder")
      ("log", value<std::string>(), "log file")
      ("fragments,f", value<std::string>(), "fragments file (.sdf file)")
      ("dx", value<std::string>(), "dx folder")
      ("no_oberrorlog", value<bool>()->implicit_value(true)->default_value(false), "stop OBEroorLog")
      ("verbose,v", value<bool>()->implicit_value(true)->default_value(false), "enable verbose output");
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
    if (vmap.count("dx")) conf.dx_folder = vmap["dx"].as<std::string>();
    conf.no_oberrorlog = vmap["no_oberrorlog"].as<bool>();
    conf.verbose = vmap["verbose"].as<bool>();
    return conf;
  }

  void logConfig(const format::DockingConfiguration config){
    logs::lout << "Receptor file name  : "+config.receptor_file   << std::endl;
    logs::lout << "Fragments file name  : "+config.fragments_file << std::endl;
    logs::lout << "Output file name    : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name    : "+config.grid_folder     << std::endl;
    logs::lout << "DX folder name      : "+config.dx_folder       << std::endl;
  }

  std::string escape_filename(const std::string& filename) {
    std::string result;
    for (char c : filename) {
      switch (c) {
        case '\\':
          result += "Yen";
          break;
        case '/':
          result += "Sla";
          break;
        case ':':
          result += "Col";
          break;
        case '*':
          result += "Aste";
          break;
        case '?':
          result += "Ques";
          break;
        case '\"':
          result += "Dquo";
          break;
        case '>':
          result += "Gth";
          break;
        case '<':
          result += "Lth";
          break;
        case '|':
          result += "Vbar";
          break;
        default:
          result += c;
      }
    }
    return result;
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
    config.log_file = config.output_file + "fraggdx__" + getDate() + ".log";
  }
  logs::log_init(config.log_file, 1);
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

  const Point3d<int>& score_num = atom_grids[0].getNum(); // # of grid points of score grids (atom grids, fragment grids)


  // vector<Fragment> frag_library; /* list of unique fragments which poses are normalized*/
  // vector<int> frag_importance; /* # fragment heavy atoms * fragment occurrence */
  

  // logs::lout << logs::info << "start pre-calculate energy" << endl;
  // EnergyCalculator ene_calculator(1.0); 
  // logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << logs::info << "[TIME STAMP] START MOLECULE OBJECT CONVERSION" << endl;
   vector<Fragment> fragments_frag = convert_fragments(fragments);
  logs::lout << logs::info << "[TIME STAMP] END   MOLECULE OBJECT CONVERSION" << endl;


  InterEnergyGrid distance_grid = makeDistanceGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol);

  // // Amount of calculation cost reduction by reusing fragment grid
  // int reduces = 0;

  // ---------------------------------------


  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;

  std::string fraggrid_folder = config.dx_folder;
  makeFolder(fraggrid_folder);

  progress_bar frag_progress(frags_sz, "Writing fragment grids", config.verbose);

  for (int f_ind = 0; f_ind < frags_sz; ++f_ind) {
    FragmentInterEnergyGrid frag_grid(fragments_frag[f_ind], makeRotations60(), atom_grids, distance_grid);
    opendx::OpenDx open_dx(frag_grid);
    open_dx.write(fraggrid_folder + "/" + escape_filename(fragments_frag[f_ind].getSmiles()) + ".dx");
    frag_progress.display(f_ind+1);
  }
  frag_progress.clear();

  logs::lout << logs::info << "[TIME STAMP] END CALCULATING BY FRAGGRID" << endl;
  auto t3 = std::chrono::system_clock::now();
  whole_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t0);

  logs::lout << logs::info << "################ Program end ################" << endl;
  logs::lout << logs::info << "[TIME STAMP] whole_time : " << whole_time.count() << endl;

  logs::close();
  return 0;
}