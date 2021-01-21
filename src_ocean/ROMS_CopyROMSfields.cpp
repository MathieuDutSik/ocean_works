#include "ROMSfunctionality.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull = NAMELIST_CREATE_DEFAULT_SETTING();
    if (argc != 2) {
      std::cerr << "ROMS_CopyROMSfields [file.nml]\n";
      std::cerr << "with file.nml the file describing the default variable chosen\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    CreateDefaultInputFiles(eFull);
    std::cerr << "Normal termination of ROMS_CopyROMSfields\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ROMS_CopyROMSfields\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
