#include "PointOutput.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull = NAMELIST_GetStandard_MultipleVarPlot();
    if (argc != 2) {
      std::cerr << "PLOT_PointOutput is used as\n";
      std::cerr << "PLOT_PointOutput [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    PointOutputPlot(eFull);
    std::cerr << "Normal termination of PLOT_PointOutput\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in PLOT_PointOutput\n";
    exit(e.eVal);
  }
}
