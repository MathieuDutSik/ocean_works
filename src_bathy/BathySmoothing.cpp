#include "SmoothingBathy.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_Bathymetry_Smoothing();
    if (argc != 2) {
      std::cerr << "PLOT_transect is used as\n";
      std::cerr << "PLOT_transect [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    DoFullSmoothing(eFull);
    std::cerr << "Normal termination of BathySmoothing\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
