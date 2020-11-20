#include "ROMSboundary.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    std::cerr << "Before NAMELIST_GetStandardPLOT_BOUNDARY\n";
    FullNamelist eFull = NAMELIST_GetStandardPLOT_BOUNDARY();
    std::cerr << " After NAMELIST_GetStandardPLOT_BOUNDARY\n";
    if (argc != 2) {
      std::cerr << "PLOT_bound is used as\n";
      std::cerr << "PLOT_bound [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    BOUND_Plotting_Function(eFull);
    std::cerr << "Normal termination of PLOT_bound\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_bound\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
