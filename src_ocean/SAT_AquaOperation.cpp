#include "AquaSatellite.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
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
    std::cerr << "Error in AquaOperation\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
