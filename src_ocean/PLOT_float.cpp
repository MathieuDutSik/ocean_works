#include "Floats.h"

int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_GetStandard_PlotRomsFloats();
    if (argc != 2) {
      std::cerr << "PLOT_float [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
        NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    PLOT_ROMS_float(eFull);
    std::cerr << "Normal termination of PLOT_float\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_float\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
