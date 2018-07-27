#include "River.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    FullNamelist eFull = NAMELIST_PLOT_River();
    if (argc != 2) {
      std::cerr << "PLOT_river is used as\n";
      std::cerr << "PLOT_river [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    PlotRiverInformation(eFull);
    std::cerr << "Normal termination of PLOT_river\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_river\n";
    exit(e.eVal);
  }
}
