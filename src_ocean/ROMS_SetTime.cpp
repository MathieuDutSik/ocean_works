#include "Plotting_fct.h"
int main(int argc, char *argv[])
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  try {
    if (argc != 3) {
      std::cerr << "ROMS_SetTime is used as\n";
      std::cerr << "ROMS_SetTime [NetcdfInitialFile] [Date]\n";
      std::cerr << "with NetcdfInitialFile the file considered and [Date] the date chosen\n";
      return -1;
    }
    std::string eFileName=argv[1];
    std::string eDate=argv[2];
    if (!IsExistingFile(eFileName)) {
      std::cerr << "The file eFileName=" << eFileName << " is missing\n";
      throw TerminalException{1};
    }
    netCDF::NcFile dataFile(eFileName, netCDF::NcFile::write, netCDF::NcFile::nc4);
    netCDF::NcVar eVar = dataFile.getVar("ocean_time");
    eVar.putAtt("units", "seconds since 1968-05-23 00:00:00");
    double eDate2=DATE_ConvertString2mjd(eDate);
    double eDate1=DATE_ConvertSix2mjd({1968,5,23,0,0,0});
    double eTime = (eDate2 - eDate1) * 86400;
    eVar.putVar(&eTime);
    std::cerr << "Normal termination of ROMS_SetTime\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Error in ROMS_SetTime\n";
    exit(e.eVal);
  }
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "runtime = " << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
}
