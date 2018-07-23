#include "ROMSboundary.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
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
    exit(e.eVal);
  }
}
