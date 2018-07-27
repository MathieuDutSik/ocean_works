#include "Floats.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull = NAMELIST_GetStandard_ComputeFloatTrajectories();
    if (argc != 2) {
      std::cerr << "COMP_drifters is used as\n";
      std::cerr << "COMP_drifters [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    ProcessFloatComputation(eFull);
    std::cerr << "Normal termination of COMP_drifters\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in COMP_drifters\n";
    exit(e.eVal);
  }
}
