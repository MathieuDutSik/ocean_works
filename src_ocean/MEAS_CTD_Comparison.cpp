#include "Satellite.h"

int main(int argc, char *argv[])
{
  srand_random_set();
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_GetStandardCTD_COMPARISON();
    if (argc != 2) {
      std::cerr << "MEAS_CTD_Comparison [file.nml]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_ctd_Comparison_Request(eFull);
    std::cerr << "Normal termination of MEAS_CTD_Comparison\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in MEAS_CTD_Comparison\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
