#ifndef SRC_OCEAN_DATA_ACCESS_H_
#define SRC_OCEAN_DATA_ACCESS_H_

#include "Basic_netcdf.h"
#include "Basic_string.h"
#include "MAT_Matrix.h"
#include "Namelist.h"
#include "SphericalGeom.h"
#include "mjdv2.h"
#include <map>
#include <string>
#include <vector>

struct RadarMeasLine {
  double lon;
  double lat;
  double U;
  double V;
  double QCcur;
  double Wheight;
  double Wdir;
  double QCwave;
};

struct OneArrayRadarMeas {
  double date;
  std::vector<RadarMeasLine> ListMeas;
};

std::vector<OneArrayRadarMeas> ReadRadarFile(std::string const &eFileName) {
  std::cerr << "eFileName = " << eFileName << "\n";
  if (!IsExistingFile(eFileName)) {
    std::cerr << "Error in ReadRadarFile\n";
    std::cerr << "Missing file eFileName = " << eFileName << "\n";
    throw TerminalException{1};
  }
  std::vector<OneArrayRadarMeas> ListRet;
  std::vector<std::vector<int>> ListTimePos;
  std::ifstream INfs(eFileName);
  int pos = -1;
  auto FuncInsert = [&](double const &date, std::vector<int> const &eTimePos,
                        RadarMeasLine const &eMeas) -> void {
    bool IsCorrect;
    if (pos == -1) {
      IsCorrect = false;
    } else {
      if (ListTimePos[pos] == eTimePos)
        IsCorrect = true;
      else
        IsCorrect = false;
    }
    if (IsCorrect) {
      ListRet[pos].ListMeas.push_back(eMeas);
    } else {
      int nbRet = ListRet.size();
      for (int iRet = 0; iRet < nbRet; iRet++)
        if (ListTimePos[iRet] == eTimePos) {
          ListRet[iRet].ListMeas.push_back(eMeas);
          return;
        }
      ListRet.push_back({date, {eMeas}});
      ListTimePos.push_back(eTimePos);
      pos = nbRet;
    }
  };
  auto MySplit = [&](std::string const &eStr) -> std::vector<std::string> {
    int len = eStr.size();
    //    std::cerr << "len=" << len << "\n";
    if (len == 0)
      return {};
    //    std::cerr << "MySplit, step 1\n";
    std::vector<std::string> ListStr = STRING_Split_Strict(eStr, ";");
    //    std::cerr << "MySplit, step 2\n";
    std::string LastChar = eStr.substr(len - 1, 1);
    //    std::cerr << "MySplit, step 3\n";
    if (LastChar != std::string(";"))
      return ListStr;
    //    std::cerr << "MySplit, step 4\n";
    ListStr.push_back("");
    //    std::cerr << "MySplit, step 5\n";
    return ListStr;
  };
  bool IsFirst = true;
  auto ExtractDouble = [](std::string const &eS) -> double {
    if (eS.size() == 0)
      return -double(999);
    double eVal;
    std::istringstream(eS) >> eVal;
    return eVal;
  };
  int nbLine = 0;
  while (!INfs.eof()) {
    std::string str;
    std::getline(INfs, str);
    if (!IsFirst) {
      //      std::cerr << "Passing here\n";
      //      std::cerr << "str=" << str << "\n";
      std::vector<std::string> ListStr = MySplit(str);
      //      std::cerr << "Step 1 : |ListStr|=" << ListStr.size() << "\n";
      if (ListStr.size() == 13) {
        nbLine++;
        int day, month, year, hour, min, sec = 0;
        std::istringstream(ListStr[0]) >> day;
        std::istringstream(ListStr[1]) >> month;
        std::istringstream(ListStr[2]) >> year;
        std::istringstream(ListStr[3]) >> hour;
        std::istringstream(ListStr[4]) >> min;
        std::vector<int> eDate{year, month, day, hour, min, sec};
        double date = DATE_ConvertSix2mjd(eDate);
        double lat, lon, U, V, QCcur, Wheight, Wdir, QCwave;
        lat = ExtractDouble(ListStr[5]);
        lon = ExtractDouble(ListStr[6]);
        U = ExtractDouble(ListStr[7]);
        V = ExtractDouble(ListStr[8]);
        QCcur = ExtractDouble(ListStr[9]);
        Wheight = ExtractDouble(ListStr[10]);
        Wdir = ExtractDouble(ListStr[11]);
        QCwave = ExtractDouble(ListStr[12]);
        /*
        std::istringstream(ListStr[5])  >> lat;
        std::istringstream(ListStr[6])  >> lon;
        std::istringstream(ListStr[7])  >> U;
        std::istringstream(ListStr[8])  >> V;
        std::istringstream(ListStr[9])  >> QCcur;
        std::istringstream(ListStr[10]) >> Wheight;
        std::istringstream(ListStr[11]) >> Wdir;
        std::istringstream(ListStr[12]) >> QCwave;*/
        RadarMeasLine eRad{lon, lat, U, V, QCcur, Wheight, Wdir, QCwave};
        FuncInsert(date, eDate, eRad);
      }
    }
    IsFirst = false;
  }
  std::cerr << "Return |ListRet|=" << ListRet.size() << " nbLine=" << nbLine
            << "\n";
  return ListRet;
}

struct DATAaqua {
  bool IsCorrect;
  double minTime;
  double maxTime;
  MyMatrix<double> LON;
  MyMatrix<double> LAT;
  MyMatrix<int> MSK;
  MyMatrix<double> VAR;
};

DATAaqua ReadAquaSatellite(std::string const &DataFile,
                           std::string const &VarName) {
  std::cerr << "Beginning of ReadAquaSatellite DataFile=" << DataFile << "\n";
  if (!IsExistingFile(DataFile)) {
    std::cerr << "ReadAquaSatellite : The file DataFile = " << DataFile
              << " is missing\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(DataFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "ReadAquaSatellite : null entry in the file\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  //
  // Reading variables
  //
  std::cerr << "Reading variables\n";
  netCDF::NcDim dimPixel = dataFile.getDim("pixels_per_line");
  if (dimPixel.isNull()) {
    std::cerr << "ReadAquaSatellite : null entry for dimPixel\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  int pixels_per_line = dimPixel.getSize();
  //
  netCDF::NcDim dimLine = dataFile.getDim("number_of_lines");
  if (dimLine.isNull()) {
    std::cerr << "ReadAquaSatellite : null entry for dimLine\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  int number_of_lines = dimLine.getSize();
  std::cerr << "pixels_per_line=" << pixels_per_line
            << "   number_of_lines=" << number_of_lines << "\n";
  //
  // Reading groups
  //
  //  std::cerr << "Reading groups\n";
  netCDF::NcGroup GrpGeophys =
      dataFile.getGroup("geophysical_data", netCDF::NcGroup::ChildrenGrps);
  if (GrpGeophys.isNull()) {
    std::cerr << "ReadAquaSatellite : error when accessing geophysical_data\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "We have GrpGeophys\n";
  //
  netCDF::NcGroup GrpLine =
      dataFile.getGroup("scan_line_attributes", netCDF::NcGroup::ChildrenGrps);
  if (GrpLine.isNull()) {
    std::cerr
        << "ReadAquaSatellite : error when accessing scan_line_attributes\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "We have GrpLine\n";
  //
  netCDF::NcGroup GrpNav =
      dataFile.getGroup("navigation_data", netCDF::NcGroup::ChildrenGrps);
  if (GrpLine.isNull()) {
    std::cerr << "ReadAquaSatellite : error when accessing navigation_data\n";
    std::cerr << "DataFile = " << DataFile << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "We have GrpNav\n";
  //
  // Printing the group names
  //
  bool DoPrintGroupName = false;
  if (DoPrintGroupName) {
    std::multimap<std::string, netCDF::NcGroup> ListGroup =
        dataFile.getGroups(netCDF::NcGroup::ChildrenGrps);
    auto iter = ListGroup.begin();
    while (iter != ListGroup.end()) {
      std::cerr << "Group=" << iter->first << "\n";
      iter++;
    }
  }
  //
  // Reading the longitude and latitudes.
  //
  std::cerr << "Reading longitude/latitude\n";
  netCDF::NcVar data_slon = GrpLine.getVar("slon");
  std::cerr << "We have data_slon\n";
  MyVector<double> slon = NC_Read1Dvariable_data(data_slon);
  //  std::cerr << "Reading longitude/latitude, step 1\n";
  //
  netCDF::NcVar data_clon = GrpLine.getVar("clon");
  MyVector<double> clon = NC_Read1Dvariable_data(data_clon);
  //  std::cerr << "Reading longitude/latitude, step 2\n";
  //
  netCDF::NcVar data_elon = GrpLine.getVar("elon");
  MyVector<double> elon = NC_Read1Dvariable_data(data_elon);
  //  std::cerr << "Reading longitude/latitude, step 3\n";
  //
  netCDF::NcVar data_slat = GrpLine.getVar("slat");
  MyVector<double> slat = NC_Read1Dvariable_data(data_slat);
  //  std::cerr << "Reading longitude/latitude, step 4\n";
  //
  netCDF::NcVar data_clat = GrpLine.getVar("clat");
  MyVector<double> clat = NC_Read1Dvariable_data(data_clat);
  //  std::cerr << "Reading longitude/latitude, step 5\n";
  //
  netCDF::NcVar data_elat = GrpLine.getVar("elat");
  MyVector<double> elat = NC_Read1Dvariable_data(data_elat);
  //  std::cerr << "Reading longitude/latitude, step 6\n";
  //
  // Reading the times
  //
  std::cerr << "Reading times\n";
  netCDF::NcVar data_year = GrpLine.getVar("year");
  MyVector<int> ListYear = NC_Read1Dvariable_int_data(data_year);
  //
  netCDF::NcVar data_day = GrpLine.getVar("day");
  MyVector<int> ListDay = NC_Read1Dvariable_int_data(data_day);
  //
  netCDF::NcVar data_msec = GrpLine.getVar("msec");
  MyVector<int> ListMSec = NC_Read1Dvariable_int_data(data_msec);
  //
  int nbTime = ListYear.size();
  std::vector<double> ListTime(nbTime);
  for (int iTime = 0; iTime < nbTime; iTime++) {
    int iYear = ListYear[iTime];
    int iDay = ListDay[iTime];
    int MSec = ListMSec[iTime];
    double eTimeFirst = DATE_ConvertSix2mjd({iYear, 1, 1, 0, 0, 0});
    double eTimeDay1 = eTimeFirst + double(iDay - 1);
    double eTimeFinal = eTimeDay1 + double(MSec) / double(1000 * 86400);
    /*    if (iTime == 0 || iTime == nbTime-1) {
      std::cerr << "iTime=" << iTime << "\n";
      std::cerr << "  iYear=" << iYear << " iDay=" << iDay << "  MSec=" << MSec
      << "\n"; std::cerr << "  eTimeFirst=" << eTimeFirst << "\n"; std::cerr <<
      "  eTimeDay1=" << eTimeDay1 << "\n"; std::cerr << "  eTimeFinal=" <<
      eTimeFinal << "\n";
      }*/
    ListTime[iTime] = eTimeFinal;
  }
  double minTime = VectorMin(ListTime);
  double maxTime = VectorMax(ListTime);
  std::cerr << "minTime=" << minTime << " maxTime=" << maxTime << "\n";
  //
  // Reading number of missing values
  //
  /*
  std::vector<std::string> ListVarName{"chlor_a", "chl_ocx", "pic", "poc",
  "ipar", "par"}; for (auto & eVar : ListVarName) { netCDF::NcVar
  data_fill=GrpGeophys.getVar(eVar); int
  nbFill=NC_ReadVariable_NbFillValue_data(data_fill); std::cerr << "eVar=" <<
  eVar << " nbFill=" << nbFill << "\n";
    }*/
  //
  // Reading data
  //
  std::cerr << "Reading data\n";
  netCDF::NcVar data_read = GrpGeophys.getVar(VarName);
  MyMatrix<double> VarField;
  VarField = NC_Read2Dvariable_data(data_read);
  int nbRow = VarField.rows();
  int nbCol = VarField.cols();
  if (nbRow != number_of_lines || nbCol != pixels_per_line) {
    std::cerr << "We have nbRow=" << nbRow << " and nbCol=" << nbCol << "\n";
    std::cerr << "But pixels_per_line=" << pixels_per_line
              << " number_of_lines=" << number_of_lines << "\n";
    throw TerminalException{1};
  }
  //
  // Reading Mask from the FillValue
  //
  MyMatrix<int> MSK = NC_Read2Dvariable_Mask_data(data_read);
  std::cerr << "|MSK|=" << MSK.rows() << " / " << MSK.cols() << "\n";
  //
  // Construction of the LON/LAT arrays
  //
  std::cerr << "Setting up LON/LAT arrays\n";
  bool ComputeLLanal = false;
  if (ComputeLLanal) {
    MyMatrix<double> LONanal(nbRow, nbCol);
    MyMatrix<double> LATanal(nbRow, nbCol);
    for (int iLine = 0; iLine < number_of_lines; iLine++) {
      double sLon = slon(iLine);
      double eLon = elon(iLine);
      double sLat = slat(iLine);
      double eLat = elat(iLine);
      std::vector<double> sXYZ = GetXYZcoordinateLL(sLon, sLat);
      std::vector<double> eXYZ = GetXYZcoordinateLL(eLon, eLat);
      for (int iPixel = 0; iPixel < pixels_per_line; iPixel++) {
        double sCoeff = double(iPixel);
        double eCoeff = double(pixels_per_line - 1 - iPixel);
        std::vector<double> TheXYZ = GetXYZaverage(sXYZ, eXYZ, sCoeff, eCoeff);
        std::vector<double> TheLL = GetLLcoordinateXYZ(TheXYZ);
        LONanal(iLine, iPixel) = TheLL[0];
        LATanal(iLine, iPixel) = TheLL[1];
      }
    }
  }
  std::vector<int> ListStatus(number_of_lines);
  std::vector<int> ListCorrect;
  for (int iLine = 0; iLine < number_of_lines; iLine++) {
    int eStatus = 1;
    if (elon(iLine) < -200 || elon(iLine) > 200)
      eStatus = 0;
    if (slon(iLine) < -200 || slon(iLine) > 200)
      eStatus = 0;
    if (elat(iLine) < -100 || elat(iLine) > 100)
      eStatus = 0;
    if (slat(iLine) < -100 || elat(iLine) > 100)
      eStatus = 0;
    ListStatus[iLine] = eStatus;
    if (eStatus)
      ListCorrect.push_back(iLine);
  }
  int iFirstLine = VectorMin(ListCorrect);
  int iLastLine = VectorMax(ListCorrect);
  int nbLineRed = iLastLine + 1 - iFirstLine;
  //
  netCDF::NcVar data_lon = GrpNav.getVar("longitude");
  MyMatrix<double> LONfile = NC_Read2Dvariable_data(data_lon);
  //
  netCDF::NcVar data_lat = GrpNav.getVar("latitude");
  MyMatrix<double> LATfile = NC_Read2Dvariable_data(data_lat);
  //
  // Now we determine if there is need for doing some plots
  //
  DATAaqua eDATA;
  eDATA.minTime = minTime;
  eDATA.maxTime = maxTime;
  std::cerr << "|LONfile|=" << LONfile.rows() << " / " << LONfile.cols()
            << "\n";
  MyMatrix<double> LONred(nbLineRed, pixels_per_line);
  MyMatrix<double> LATred(nbLineRed, pixels_per_line);
  MyMatrix<double> VARred(nbLineRed, pixels_per_line);
  MyMatrix<int> MSKred(nbLineRed, pixels_per_line);
  int idx = 0;
  for (int iLine = iFirstLine; iLine <= iLastLine; iLine++) {
    for (int i = 0; i < pixels_per_line; i++) {
      LONred(idx, i) = LONfile(iLine, i);
      LATred(idx, i) = LATfile(iLine, i);
      VARred(idx, i) = VarField(iLine, i);
      MSKred(idx, i) = MSK(iLine, i);
    }
    idx++;
  }
  eDATA.LON = LONred;
  eDATA.LAT = LATred;
  eDATA.MSK = MSKred;
  eDATA.VAR = VARred;
  eDATA.IsCorrect = true;
  return eDATA;
}

std::vector<int> ConvertAquaNameToTime(std::string const &CoreString) {
  //  std::cerr << "CoreString=" << CoreString << "\n";
  std::string strYear = CoreString.substr(1, 4);
  std::string strDay = CoreString.substr(5, 3);
  int eYear, nbDay;
  std::istringstream(strYear) >> eYear;
  std::istringstream(strDay) >> nbDay;
  double FirstDayYear = DATE_ConvertSix2mjd({eYear, 1, 1, 0, 0, 0});
  double RelTimeDay = FirstDayYear + nbDay - 1;
  std::vector<int> LTime = DATE_ConvertMjd2six(RelTimeDay);
  return LTime;
}

std::vector<std::string> GetAquaListLinkDownload_Specified(
    FullNamelist const &eFull, double const &MjdFirst, double const &MjdLast) {
  std::string rndString = random_string(20);
  std::string TheDir = "/tmp/AQUA_" + rndString;
  CreateDirectory(TheDir);
  //
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  std::string AquaLink = BlPROC.ListStringValues.at("AquaLink");
  std::vector<std::string> ListStringLink;
  auto GetSpecificDay = [&](double const &eDay) -> void {
    std::vector<int> dateVect = DATE_ConvertMjd2six(eDay);
    std::string strPres = DATE_ConvertMjd2mystringPres(eDay);
    std::cerr << "GetSpecificDay date=" << strPres << "\n";
    int eYear = dateVect[0];
    double FirstDayYear = DATE_ConvertSix2mjd({eYear, 1, 1, 0, 0, 0});
    int nbDay = int(eDay - FirstDayYear);
    std::string strDay = StringNumber(nbDay, 3);
    std::string TmpStdout =
        TheDir + "tmp_" + StringNumber(eYear, 4) + "_" + strDay + "_stdout";
    std::string TmpStderr =
        TheDir + "tmp_" + StringNumber(eYear, 4) + "_" + strDay + "_stderr";
    std::string eComm1 = "(cd " + TheDir + " && wget " + AquaLink +
                         StringNumber(eYear, 4) + "/" + strDay + " > " +
                         TmpStdout + " 2> " + TmpStderr + ")";
    int iret = system(eComm1.c_str());
    std::cerr << "iret=" << iret << "\n";
    if (iret == 0) {
      std::cerr << "File correctly downloaded |ListStringLink|="
                << ListStringLink.size() << "\n";
      std::cerr << "eComm1=" << eComm1 << "\n";
      std::string TheFile = TheDir + "/" + strDay;
      if (!IsExistingFile(TheFile)) {
        std::cerr << "We had iret=" << iret << "\n";
        std::cerr << "Which seems to indicate successful download\n";
        std::cerr << "but the file TheFile=" << TheFile << " is missing\n";
        std::cerr << "Please debug\n";
        std::cerr << "eComm1=" << eComm1 << "\n";
        std::cerr << "Good luck\n";
        throw TerminalException{1};
      }
      std::vector<std::string> ListLines = ReadFullFile(TheFile);
      for (auto &eLine : ListLines) {
        std::vector<std::string> LStr = STRING_Split(eLine, "getfile");
        if (LStr.size() == 2) {
          std::vector<std::string> LStrB = STRING_Split(eLine, "'");
          std::string eLink = LStrB[1];
          ListStringLink.push_back(eLink);
        }
      }
      RemoveFile(TheFile);
    } else {
      std::cerr << "The File could not be downloaded\n";
    }
  };
  for (double eDay = MjdFirst; eDay <= MjdLast; eDay++)
    GetSpecificDay(eDay);
  RemoveEmptyDirectory(TheDir);
  std::cerr
      << "Now leaving GetAquaListLinkDownload_Specified, |ListStringLink|="
      << ListStringLink.size() << "\n";
  return ListStringLink;
}

std::vector<std::string> GetAquaListLinkDownload(FullNamelist const &eFull) {
  std::time_t t = std::time(nullptr);
  struct tm *now = std::gmtime(&t);
  int eYear = now->tm_year + 1900;
  int eMon = now->tm_mon + 1;
  int eDay = now->tm_mday;
  double TheDay = DATE_ConvertSix2mjd({eYear, eMon, eDay, 0, 0, 0});
  //
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  int EarliestDay_shift = BlPROC.ListIntValues.at("EarliestDay_shift");
  double MjdFirst = TheDay + double(EarliestDay_shift);
  double MjdLast = TheDay + double(2);
  return GetAquaListLinkDownload_Specified(eFull, MjdFirst, MjdLast);
}

#endif  // SRC_OCEAN_DATA_ACCESS_H_
