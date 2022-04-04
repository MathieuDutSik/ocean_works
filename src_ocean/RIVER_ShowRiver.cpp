#include "River.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
  try {
    FullNamelist eFull = NAMELIST_RetrieveData();
    if (argc != 2) {
      std::cerr << "RIVER_ShowRiver [file.nml]\n";
      std::cerr << "with file.nml the file describing the river description "
                   "and the list of times\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    PrintRiverInformation(eFull);
    std::cerr << "Normal termination of RIVER_ShowRiver\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in RIVER_ShowRiver\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "runtime = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
}
