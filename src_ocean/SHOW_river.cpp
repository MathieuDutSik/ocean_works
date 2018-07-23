#include "River.h"
int main(int argc, char *argv[])
{
  std::cerr << std::fixed;
  std::cerr << std::setprecision(9);
  try {
    FullNamelist eFull = NAMELIST_RetrieveData();
    if (argc != 2) {
      std::cerr << "SHOW_river is used as\n";
      std::cerr << "SHOW_river [file.nml]\n";
      std::cerr << "with file.nml the file describing the river description and the list of times\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    PrintRiverInformation(eFull);
    std::cerr << "Normal termination of SHOW_river\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
