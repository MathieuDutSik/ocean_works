#include "Plotting_fct.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_GetStandard_PlotRoutine_pair();
    if (argc != 2) {
      std::cerr << "PLOT_diff_results is used as\n";
      std::cerr << "PLOT_diff_results [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    PAIR_Plotting_Function(eFull);
    std::cerr << "Normal termination of PLOT_diff_results\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
