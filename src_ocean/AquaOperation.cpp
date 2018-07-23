#include "AquaSatellite.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    FullNamelist eFull=NAMELIST_GetStandardAQUA();
    if (argc != 2) {
      std::cerr << "AquaOperation is used as\n";
      std::cerr << "AquaOperation [file.nml]\n";
      std::cerr << "with file.nml the file describing the input and plotting used\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    AquaDownloading(eFull);
    std::cerr << "Normal termination of AquaOperation\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
