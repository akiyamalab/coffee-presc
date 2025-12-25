#include "log_writer_stream.hpp"

#include <fstream>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/format.hpp>

namespace {
  const std::string getDateString() {
    time_t timer = time(NULL);
    tm* date = localtime(&timer);
    std::string dateStr = (boost::format("[%02d/%02d %02d:%02d:%02d]") 
			   % (date->tm_mon+1) % date->tm_mday
			   % date->tm_hour % date->tm_min % date-> tm_sec).str();
    return dateStr;
  }
  // file-local verbosity level (0=errors only, 1=info+, 2=debug+)
  int verbosity_level = 1;
}

namespace logs{
  io::filtering_ostream lout;
  LogType info("INFO");
  LogType debug("DEBUG");
  LogType warn("WARN");
  LogType error("ERROR");
  void set_verbosity(int v) {
    verbosity_level = std::max(0, std::min(2, v));
  }
  void log_init(const std::string& filename, int verbosity) {
    set_verbosity(verbosity);

    static std::ofstream ofs(filename.c_str(), std::ios::app); // append to existing log file
    static io::tee_filter<std::ostream>  fileFilt(ofs); // tee all passed data to logfile
    logs::lout.push(fileFilt); // 1st, tee off any data to the file (raw boost XML)
    // If verbosity > 0 then also write to stdout
    if (verbosity > 0) {
      logs::lout.push(std::cout); // 2nd, tee off any data to cout (raw boost XML)
    }
  }
  void close() {
    io::close(logs::lout);
  }
  std::ostream& operator<<(std::ostream& os, const LogType& lt) {
    // suppress output according to verbosity level
    if ((lt.type == "DEBUG" && verbosity_level < 2) ||
        ((lt.type == "INFO" || lt.type == "WARN") && verbosity_level < 1)) {
      static std::ofstream null_ofs("/dev/null");
      return null_ofs;
    }
    os << getDateString() << " [time=(" << time(NULL) << ")] [" << lt.type << "] ";
    return os;
  }
}