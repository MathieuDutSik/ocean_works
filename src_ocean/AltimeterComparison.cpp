#include "Satellite.h"

#ifdef DEBUG_TRAP_NAN
#include <fenv.h>
#endif


int main(int argc, char *argv[])
{
#ifdef DEBUG_TRAP_NAN
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  try {
    FullNamelist eFull=NAMELIST_GetStandardALTIMETRY_COMPARISON();
    if (argc != 2) {
      std::cerr << "AltimeterComparison is used as\n";
      std::cerr << "AltimeterComparison [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Altimetry_Comparison_Request(eFull);
    std::cerr << "Normal termination of AltimeterComparison\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
