#include "Model_interpolation.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_GetStandard_CREATE_sflux();
    if (argc != 2) {
      std::cerr << "DATA_CreateSflux [file.nml]\n";
      std::cerr << "with file.nml the file describing the choices made\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    CREATE_sflux_files(eFull);
    std::cerr << "Normal termination of CREATE_sflux\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in DATA_CreateSflux\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
