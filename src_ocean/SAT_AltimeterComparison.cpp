#include "Satellite.h"

#ifdef DEBUG_TRAP_NAN
#include <fenv.h>
#endif


int main(int argc, char *argv[])
{
#ifdef DEBUG_TRAP_NAN
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_GetStandardALTIMETRY_COMPARISON();
    if (argc != 2) {
      std::cerr << "SAT_AltimeterComparison [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Altimetry_Comparison_Request(eFull);
    std::cerr << "Normal termination of SAT_AltimeterComparison\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in SAT_AltimeterComparison\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
