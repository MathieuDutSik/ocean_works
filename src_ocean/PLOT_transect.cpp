#include "Transect.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_PlotTransect();
    if (argc != 2) {
      std::cerr << "PLOT_transect is used as\n";
      std::cerr << "PLOT_transect [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    TRANSECT_Plot(eFull);
    std::cerr << "Normal termination of PLOT_transect\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_transect\n";
    exit(e.eVal);
  }
}
