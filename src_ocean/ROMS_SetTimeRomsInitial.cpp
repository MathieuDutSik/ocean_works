// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_netcdf.h"
#include "mjdv2.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "ROMS_SetTimeRomsInitial [file.nc] timepos dateStrFile\n";
      std::cerr << "Date format is as in 20160120.000000\n";
      return -1;
    }
    //
    std::string eFileName = argv[1];
    int timepos;
    sscanf(argv[2], "%d", &timepos);
    std::string dateStrFile = argv[3];
    //
    double eTimeDay = DATE_ConvertString2mjd(dateStrFile);
    netCDF::NcFile dataFile(eFileName, netCDF::NcFile::write);
    if (dataFile.isNull()) {
      std::cerr << "Error with dataFile reading\n";
      throw TerminalException{1};
    }
    netCDF::NcVar data = dataFile.getVar("ocean_time");
    if (data.isNull()) {
      std::cerr << "Error with data reading\n";
      throw TerminalException{1};
    }
    double eTimeDayRef = DATE_ConvertSix2mjd({1968, 5, 23, 0, 0, 0});
    double eTimeWrite = 86400 * (eTimeDay - eTimeDayRef);
    //
    std::vector<size_t> start{size_t(timepos)};
    std::vector<size_t> count{1};
    data.putVar(start, count, &eTimeWrite);
    data.putAtt("long_name", "time since initialization");
    data.putAtt("units", "seconds since 1968-05-23 00:00:00");
    data.putAtt("calendar", "julian");
    data.putAtt("field", "time, scalar, series");
    //
    netCDF::NcVar data_day = dataFile.getVar("ocean_time_day");
    if (data_day.isNull()) {
      std::cerr << "No need to write ocean_time_day\n";
    } else {
      double eTimeWrite_day = eTimeDay - eTimeDayRef;
      std::cerr << "Need to write ocean_time_day\n";
      data_day.putVar(start, count, &eTimeWrite_day);
      data_day.putAtt("long_name", "time since initialization");
      data_day.putAtt("units", "days since 1968-05-23 00:00:00");
      data_day.putAtt("calendar", "julian");
      data_day.putAtt("field", "time, scalar, series");
    }
    //
    netCDF::NcVar data_str = dataFile.getVar("ocean_time_str");
    if (data_str.isNull()) {
      std::cerr << "No need to write ocean_time_str\n";
    } else {
      std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
      std::cerr << "Need to write ocean_time_str\n";
      std::vector<size_t> start2{size_t(timepos), 0};
      std::vector<size_t> count2{1, 19};
      data_str.putVar(start2, count2, strPres.c_str());
    }
    std::cerr << "Normal termination of ROMS_SetTimeRomsInitial\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ROMS_SetTimeRomsInitial\n";
    exit(e.eVal);
  }
  runtime(time1);
}
