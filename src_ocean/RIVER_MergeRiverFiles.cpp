#include "River.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
  try {
    if (argc <= 3) {
      std::cerr << "RIVER_MergeRiverFiles [TargetFile] [file1] [file2] .... "
                   "[fileN]\n";
      std::cerr << "with file.nml the file describing the chosen options\n";
      return -1;
    }
    std::string RiverFile = argv[1];
    std::vector<std::string> ListRiverFile;
    for (int i = 2; i < argc; i++) {
      std::string eFile = argv[i];
      ListRiverFile.push_back(eFile);
    }
    MergeRiverFile(RiverFile, ListRiverFile);
    std::cerr << "Normal termination of RIVER_MergeRiverFiles\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in RIVER_MergeRiverFiles\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "runtime = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
}
