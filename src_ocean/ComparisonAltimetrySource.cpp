#include "Satellite.h"
int main(int argc, char *argv[])
{
  try {
    FullNamelist eFull=NAMELIST_Comparison_Altimetry_Source();
    if (argc != 2) {
      std::cerr << "ComparisonAltimetrySource is used as\n";
      std::cerr << "ComparisonAltimetrySource [alti.nml]\n";
      std::cerr << "with alti.nml the file describing the chosen options\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Process_Comparison_Altimetry_Sources(eFull);
    std::cerr << "Normal termination of ComparisonAltimetrySource\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ComparisonAltimetrySource\n";
    exit(e.eVal);
  }
}
