#include "Model_interpolation.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    FullNamelist eFull=NAMELIST_GetStandard_CREATE_sflux();
    if (argc != 2) {
      std::cerr << "CREATE_sflux is used as\n";
      std::cerr << "CREATE_sflux [file.nml]\n";
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
    exit(e.eVal);
  }
}
