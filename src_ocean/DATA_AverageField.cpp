#include "Model_interpolation.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    FullNamelist eFull=NAMELIST_InfileAveraging();
    if (argc != 2) {
      std::cerr << "DATA_AverageField [file.nml]\n";
      std::cerr << "with file.nml the file describing the interpolation process\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    Average_field_Function(eFull);
    std::cerr << "Normal termination of DATA_AverageField\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in DATA_AverageField\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
