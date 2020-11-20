#include "PointOutput.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_GetStandard_PlotBuoy();
    if (argc != 2) {
      std::cerr << "PLOT_buoy is used as\n";
      std::cerr << "PLOT_buoy [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    BUOY_Plot(eFull);
    std::cerr << "Normal termination of PLOT_buoy\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_buoy\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
