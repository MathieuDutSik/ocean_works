#include "Plotting_fct.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    std::cerr << "Before NAMELIST_GetStandard_PlotRoutine_single\n";
    FullNamelist eFull = NAMELIST_GetStandard_PlotRoutine_single();
    std::cerr << " After NAMELIST_GetStandard_PlotRoutine_single\n";
    if (argc != 2) {
      std::cerr << "PLOT_results is used as\n";
      std::cerr << "PLOT_results [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SINGLE_Plotting_Function(eFull);
    std::cerr << "Normal termination of PLOT_results\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
