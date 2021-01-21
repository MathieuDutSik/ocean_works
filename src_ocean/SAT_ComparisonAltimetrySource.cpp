#include "Satellite.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_Comparison_Altimetry_Source();
    if (argc != 2) {
      std::cerr << "SAT_ComparisonAltimetrySource is used as\n";
      std::cerr << "SAT_ComparisonAltimetrySource [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Comparison_Altimetry_Sources(eFull);
    std::cerr << "Normal termination of SAT_ComparisonAltimetrySource\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in SAT_ComparisonAltimetrySource\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
