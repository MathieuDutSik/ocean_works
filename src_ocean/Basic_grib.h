// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_BASIC_GRIB_H_
#define SRC_OCEAN_BASIC_GRIB_H_

#include "ArrHistory.h"
#include "Basic_Ocean_types.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "MAT_Matrix.h"
#include "SphericalGeom.h"
#include "Timings.h"
#include "grib_api.h"
#include "mjdv2.h"
#include <string>
#include <utility>
#include <vector>

struct CosmoGridInfo {
  double latitudeOfSouthernPoleInDegrees;
  double longitudeOfSouthernPoleInDegrees;
  double angleOfRotationInDegrees;
  //
  double longitudeOfFirstGridPointInDegrees;
  double latitudeOfFirstGridPointInDegrees;
  double longitudeOfLastGridPointInDegrees;
  double latitudeOfLastGridPointInDegrees;
  //
  double iDirectionIncrementInDegrees;
  double jDirectionIncrementInDegrees;
};

void GRIB_CheckAllowedKeywords(std::string const &eModelName) {
  std::vector<std::string> LStr = STRING_Split(eModelName, ":");
  std::vector<std::string> ListAllowed{
      "optimaltime", "shifttime", "timestartfromfilename", "retrieveallstates"};
  for (size_t i = 1; i < LStr.size(); i++) {
    std::vector<std::string> LStrB = STRING_Split(LStr[i], "_");
    std::string eStr = LStrB[0];
    if (PositionVect(ListAllowed, eStr) == -1) {
      std::cerr << "Found keyword = " << eStr << "\n";
      std::cerr << "ListAllowed =";
      for (auto &fStr : ListAllowed) {
        std::cerr << " " << fStr;
      }
      throw TerminalException{1};
    }
  }
}

double phirot2phi(double const &phirot, double const &rlarot,
                  double const &polphi, double const &pollam,
                  double const &polgam) {
  double pi = 3.1415926535;
  double d180 = 180;
  double eMult = pi / d180;
  double eMultInv = d180 / pi;
  double zsinpol = sin(eMult * polphi);
  double zcospol = cos(eMult * polphi);
  double zphis = eMult * phirot;
  double zrlas;
  if (rlarot > d180) {
    zrlas = rlarot - 360;
  } else {
    zrlas = rlarot;
  }
  zrlas = eMult * zrlas;
  double zarg;
  if (fabs(polgam) > 0) {
    double zgam = eMult * polgam;
    zarg = zsinpol * sin(zphis) +
           zcospol * cos(zphis) *
               (cos(zrlas) * cos(zgam) - sin(zgam) * sin(zrlas));
  } else {
    zarg = zcospol * cos(zphis) * cos(zrlas) + zsinpol * sin(zphis);
  }
  double phirot2phi = eMultInv * asin(zarg);
  return phirot2phi;
}

double rlarot2rla(double const &phirot, double const &rlarot,
                  double const &polphi, double const &pollam,
                  double const &polgam) {
  double pi = 3.1415926535;
  double d180 = 180;
  double eMult = pi / d180;
  double eMultInv = d180 / pi;
  double zsinpol = sin(eMult * polphi);
  double zcospol = cos(eMult * polphi);
  double zphis = eMult * phirot;
  double zrlas;
  if (rlarot > d180) {
    zrlas = rlarot - 360;
  } else {
    zrlas = rlarot;
  }
  zrlas = eMult * zrlas;
  double zlampol = eMult * pollam;
  double zarg1, zarg2;
  if (fabs(polgam) > 0) {
    double zgam = eMult * polgam;
    zarg1 =
        sin(zlampol) * (-zsinpol * cos(zphis) *
                            (cos(zrlas) * cos(zgam) - sin(zrlas) * sin(zgam)) +
                        zcospol * sin(zphis)) -
        cos(zlampol) * cos(zphis) *
            (sin(zrlas) * cos(zgam) + cos(zrlas) * sin(zgam));
    zarg2 =
        cos(zlampol) * (-zsinpol * cos(zphis) *
                            (cos(zrlas) * cos(zgam) - sin(zrlas) * sin(zgam)) +
                        zcospol * sin(zphis)) +
        sin(zlampol) * cos(zphis) *
            (sin(zrlas) * cos(zgam) + cos(zrlas) * sin(zgam));
  } else {
    zarg1 = sin(zlampol) *
                (-zsinpol * cos(zrlas) * cos(zphis) + zcospol * sin(zphis)) -
            cos(zlampol) * sin(zrlas) * cos(zphis);
    zarg2 = cos(zlampol) *
                (-zsinpol * cos(zrlas) * cos(zphis) + zcospol * sin(zphis)) +
            sin(zlampol) * sin(zrlas) * cos(zphis);
  }
  //  if (zarg2 == 0) zarg2=1.0e-20;
  double rlarot2rla = eMultInv * atan2(zarg1, zarg2);
  return rlarot2rla;
}

void Apply_COSMO_Transformation(MyMatrix<double> &LON, MyMatrix<double> &LAT,
                                CosmoGridInfo const &eCosmoGrid) {
  double pollat_sp = eCosmoGrid.latitudeOfSouthernPoleInDegrees;
  double pollon_sp = eCosmoGrid.longitudeOfSouthernPoleInDegrees;
  double polgam = eCosmoGrid.angleOfRotationInDegrees;
  double zstartlon_tot = eCosmoGrid.longitudeOfFirstGridPointInDegrees;
  double zstartlat_tot = eCosmoGrid.latitudeOfFirstGridPointInDegrees;
  double zendlon_tot = eCosmoGrid.longitudeOfLastGridPointInDegrees;
  double zendlat_tot = eCosmoGrid.latitudeOfLastGridPointInDegrees;
  double dlon = eCosmoGrid.iDirectionIncrementInDegrees;
  double dlat = eCosmoGrid.jDirectionIncrementInDegrees;
  double tolLL = static_cast<double>(1) / static_cast<double>(100000);
  if (fabs(zendlon_tot - zstartlon_tot) < tolLL ||
      fabs(zendlat_tot - zstartlat_tot) < tolLL) {
    std::cerr << "Error of consistency in zstartlat / zendlat\n";
    throw TerminalException{1};
  }
  //
  int eta_rho = LON.rows();
  int xi_rho = LON.cols();
  double pollat = -pollat_sp;
  double pollon = pollon_sp - static_cast<double>(180);
  // For zstartlon_tot / zstartlat_tot we are unsure. Maybe there is a shift
  double startlon_tot = zstartlon_tot;
  double startlat_tot = zstartlat_tot;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      double eLonR = startlon_tot + i * dlon;
      double eLatR = startlat_tot + j * dlat;
      double eLat = phirot2phi(eLatR, eLonR, pollat, pollon, polgam);
      double eLon = rlarot2rla(eLatR, eLonR, pollat, pollon, polgam);
      LON(i, j) = eLon;
      LAT(i, j) = eLat;
    }
}

GridArray GRIB_ReadGridArray(std::string const &FileName,
                             std::string const &eModelName) {
  //  grib_context *c;
  grib_handle *h = NULL;
  int err;
  FILE *in = NULL;
  in = fopen(FileName.c_str(), "r");
  //  unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;

  while ((h = grib_handle_new_from_file(0, in, &err)) != NULL) {
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err, 0);
    //
    long Ni, Nj, numberOfDataPoints;
    GRIB_CHECK(grib_get_long(h, "Ni", &Ni), 0);
    //    std::cerr << "Ni=" << Ni << "\n";
    GRIB_CHECK(grib_get_long(h, "Nj", &Nj), 0);
    //    std::cerr << "Nj=" << Nj << "\n";
    GRIB_CHECK(grib_get_long(h, "numberOfDataPoints", &numberOfDataPoints), 0);
    CosmoGridInfo eCosmoGrid;
    if (eModelName == "GRIB_COSMO") {
      double latitudeOfSouthernPoleInDegrees, longitudeOfSouthernPoleInDegrees,
          angleOfRotationInDegrees;
      double latitudeOfFirstGridPointInDegrees,
          longitudeOfFirstGridPointInDegrees, latitudeOfLastGridPointInDegrees,
          longitudeOfLastGridPointInDegrees;
      double iDirectionIncrementInDegrees, jDirectionIncrementInDegrees;
      //
      // The southern pole coordinates
      //
      GRIB_CHECK(grib_get_double(h, "latitudeOfSouthernPoleInDegrees",
                                 &latitudeOfSouthernPoleInDegrees),
                 0);
      //      std::cerr << "latitudeOfSouthernPoleInDegrees=" <<
      //      latitudeOfSouthernPoleInDegrees << "\n";
      eCosmoGrid.latitudeOfSouthernPoleInDegrees =
          latitudeOfSouthernPoleInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "longitudeOfSouthernPoleInDegrees",
                                 &longitudeOfSouthernPoleInDegrees),
                 0);
      //      std::cerr << "longitudeOfSouthernPoleInDegrees=" <<
      //      longitudeOfSouthernPoleInDegrees << "\n";
      eCosmoGrid.longitudeOfSouthernPoleInDegrees =
          longitudeOfSouthernPoleInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "angleOfRotationInDegrees",
                                 &angleOfRotationInDegrees),
                 0);
      //      std::cerr << "angleOfRotationInDegrees=" <<
      //      angleOfRotationInDegrees << "\n";
      eCosmoGrid.angleOfRotationInDegrees = angleOfRotationInDegrees;
      //
      // The first and last longitudes
      //
      GRIB_CHECK(grib_get_double(h, "latitudeOfFirstGridPointInDegrees",
                                 &latitudeOfFirstGridPointInDegrees),
                 0);
      //      std::cerr << "latitudeOfFirstGridPointInDegrees=" <<
      //      latitudeOfFirstGridPointInDegrees << "\n";
      eCosmoGrid.latitudeOfFirstGridPointInDegrees =
          latitudeOfFirstGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "longitudeOfFirstGridPointInDegrees",
                                 &longitudeOfFirstGridPointInDegrees),
                 0);
      //      std::cerr << "longitudeOfFirstGridPointInDegrees=" <<
      //      longitudeOfFirstGridPointInDegrees << "\n";
      eCosmoGrid.longitudeOfFirstGridPointInDegrees =
          longitudeOfFirstGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "latitudeOfLastGridPointInDegrees",
                                 &latitudeOfLastGridPointInDegrees),
                 0);
      //      std::cerr << "latitudeOfLastGridPointInDegrees=" <<
      //      latitudeOfLastGridPointInDegrees << "\n";
      eCosmoGrid.latitudeOfLastGridPointInDegrees =
          latitudeOfLastGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "longitudeOfLastGridPointInDegrees",
                                 &longitudeOfLastGridPointInDegrees),
                 0);
      //      std::cerr << "longitudeOfLastGridPointInDegrees=" <<
      //      longitudeOfLastGridPointInDegrees << "\n";
      eCosmoGrid.longitudeOfLastGridPointInDegrees =
          longitudeOfLastGridPointInDegrees;
      //
      // the increments
      //
      GRIB_CHECK(grib_get_double(h, "iDirectionIncrementInDegrees",
                                 &iDirectionIncrementInDegrees),
                 0);
      //      std::cerr << "iDirectionIncrementInDegrees=" <<
      //      iDirectionIncrementInDegrees << "\n";
      eCosmoGrid.iDirectionIncrementInDegrees = iDirectionIncrementInDegrees;
      //
      GRIB_CHECK(grib_get_double(h, "jDirectionIncrementInDegrees",
                                 &jDirectionIncrementInDegrees),
                 0);
      //      std::cerr << "jDirectionIncrementInDegrees=" <<
      //      jDirectionIncrementInDegrees << "\n";
      eCosmoGrid.jDirectionIncrementInDegrees = jDirectionIncrementInDegrees;
      //
    }
    //    std::cerr << "NumberOfDataPoints=" << numberOfDataPoints << "\n";
    size_t size = numberOfDataPoints;
    std::vector<double> lats(size), lons(size), values(size);
    err = grib_get_data(h, lats.data(), lons.data(), values.data());
    grib_handle_delete(h);
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err, 0);
    int eta_rho = Ni;
    int xi_rho = Nj;
    MyMatrix<double> LON(eta_rho, xi_rho);
    MyMatrix<double> LAT(eta_rho, xi_rho);
    MyMatrix<uint8_t> MSK(eta_rho, xi_rho);
    int idx = 0;
    for (int j = 0; j < xi_rho; j++)
      for (int i = 0; i < eta_rho; i++) {
        LON(i, j) = lons[idx];
        LAT(i, j) = lats[idx];
        MSK(i, j) = 1;
        idx++;
      }
    if (eModelName == "GRIB_COSMO")
      Apply_COSMO_Transformation(LON, LAT, eCosmoGrid);
    std::cerr << "GRIB: [0,0]            lon=" << LON(0, 0)
              << " lat=" << LAT(0, 0) << "\n";
    std::cerr << "      [eta_rho,0]      lon=" << LON(eta_rho - 1, 0)
              << " lat=" << LAT(eta_rho - 1, 0) << "\n";
    std::cerr << "      [eta_rho,xi_rho] lon=" << LON(eta_rho - 1, xi_rho - 1)
              << " lat=" << LAT(eta_rho - 1, xi_rho - 1) << "\n";
    std::cerr << "      [0,xi_rho]       lon=" << LON(0, xi_rho - 1)
              << " lat=" << LAT(0, xi_rho - 1) << "\n";
    MyMatrix<double> ANG = CreateAngleMatrix(LON, LAT);
    GridArray GrdArr;
    GrdArr.ModelName = eModelName;
    GrdArr.IsFE = 0;
    GrdArr.IsSpherical = true;
    GrdArr.GrdArrRho.LON = LON;
    GrdArr.GrdArrRho.LAT = LAT;
    GrdArr.GrdArrRho.MSK = MSK;
    GrdArr.GrdArrRho.ANG = ANG;
    GrdArr.ARVD.IsAssigned = false;
    GrdArr.ARVD.Zcoordinate = false;
    return GrdArr;
  }
  std::cerr << "Failed to find the variable. Error in GRIB_ReadGridArray\n";
  throw TerminalException{1};
}

void PrintGRIBmessageInfo(std::ostream &os, GRIB_MessageInfo const &eMesg) {
  os << "FileName=" << eMesg.FileName << "\n";
  os << "  shortName=" << eMesg.shortName
     << "  |shortName|=" << eMesg.shortName.size() << "\n";
  os << "  cfVarName=" << eMesg.cfVarName << "\n";
  os << "  cfVarNameECMWF=" << eMesg.cfVarNameECMWF << "\n";
  os << "  name=" << eMesg.name << "\n";
  os << "  units=" << eMesg.units << "\n";
  os << "  time=" << eMesg.time
     << " string=" << DATE_ConvertMjd2mystringPres(eMesg.time) << "\n";
  os << "  timeStart=" << eMesg.timeStart
     << " string=" << DATE_ConvertMjd2mystringPres(eMesg.timeStart) << "\n";
  os << "  stepRange=" << eMesg.stepRange << "\n";
}

void PrintVectorGRIBmessageInfo(std::ostream &os,
                                std::vector<GRIB_MessageInfo> const &ListMesg) {
  int nbMesg = ListMesg.size();
  for (int iMesg = 0; iMesg < nbMesg; iMesg++) {
    os << "MESSAGE info iMesg=" << iMesg << " / " << nbMesg << "\n";
    PrintGRIBmessageInfo(os, ListMesg[iMesg]);
  }
}

/* Some runs have some shift in time by say 1H.
   Using this option, we can correct for it.
*/
double GetShiftTime(std::string const &eModelName) {
  std::vector<std::string> LStr = STRING_Split(eModelName, ":");
  for (size_t i = 1; i < LStr.size(); i++) {
    std::vector<std::string> LStrB = STRING_Split(LStr[i], "_");
    if (LStrB[0] == "shifttime") {
      double eValHour;
      std::istringstream(LStrB[1]) >> eValHour;
      double eValDay = eValHour / static_cast<double>(24);
      return eValDay;
    }
  }
  return 0;
}

double ExtractTimeStartFromName(double const &PreTimeStart, double const &eTime,
                                std::string const &eModelName,
                                std::string const &FileName) {
  std::vector<std::string> LStr = STRING_Split(eModelName, ":");
  for (size_t i = 1; i < LStr.size(); i++) {
    if (LStr[i] == "timestartfromfilename") {
      std::vector<std::string> LStrB = STRING_Split(FileName, "+");
      if (LStrB.size() != 2) {
        std::cerr << "We cannot extract timestart from filename. Heuristic for "
                     "ALADIN does not apply\n";
        throw TerminalException{1};
      }
      std::string str = LStrB[1];
      int len1 = str.size();
      int WeMatch = false;
      std::string strB;
      if (len1 == 8) {
        strB = str.substr(0, 4);
        WeMatch = true;
      }
      if (len1 == 6) {
        strB = str.substr(0, 2);
        WeMatch = true;
      }
      if (!WeMatch) {
        std::cerr << "str=" << str
                  << " but it should have length 6 or 8, i.e. be of the form "
                     "0048.grb or 48.grb\n";
        throw TerminalException{1};
      }
      double eValHour;
      std::istringstream(strB) >> eValHour;
      double eValDay = eValHour / static_cast<double>(24);
      return eTime - eValDay;
    }
  }
  return PreTimeStart;
}

std::vector<GRIB_MessageInfo>
GRIB_GetAllMessagesFromFile(std::string const &FileName,
                            std::string const &eModelName) {
  grib_handle *h = NULL;
  int err;
  FILE *in = NULL;
  in = fopen(FileName.c_str(), "r");
  unsigned long key_iterator_filter_flags = GRIB_KEYS_ITERATOR_ALL_KEYS;
  std::vector<GRIB_MessageInfo> ListInfo;
  int idx = 0;
  while ((h = grib_handle_new_from_file(0, in, &err)) != NULL) {
    //    std::cerr << "idx=" << idx << "\n";
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err, 0);
    //
    //    std::string eStr="ls";
    char name_space[3] = "ls";
    grib_keys_iterator *kiter = NULL;
    kiter = grib_keys_iterator_new(h, key_iterator_filter_flags, name_space);
    if (!kiter) {
      printf("ERROR: Unable to create keys iterator\n");
      throw TerminalException{1};
    }
    std::string ShortNameValue;
    std::string NameValue;
    std::string cfVarName;
    std::string cfVarNameECMWF;
    std::string units;
    std::string DataDateValue = "unset";
    std::string DataTimeValue = "unset";
    std::string StepRangeValue;
    //    std::cerr << "              Before the while loop over keys\n";
    while (grib_keys_iterator_next(kiter)) {
      const size_t MAX_VAL_LEN = 1024;
      size_t eff_len;
      char value[MAX_VAL_LEN];
      const char *name = grib_keys_iterator_get_name(kiter);
      bzero(value, MAX_VAL_LEN);
      int err = grib_get_string(h, name, value, &eff_len);
      GRIB_CHECK(err, 0);
      std::string nameStr = name;
      std::string valueStr = value;
      //      std::cerr << "nameStr=" << nameStr << " valueStr=" << valueStr <<
      //      "\n";
      if (nameStr == "shortName")
        ShortNameValue = valueStr;
      if (nameStr == "name")
        NameValue = valueStr;
      if (nameStr == "cfVarName")
        cfVarName = valueStr;
      if (nameStr == "cfVarNameECMWF")
        cfVarNameECMWF = valueStr;
      if (nameStr == "units")
        units = valueStr;
      if (nameStr == "dataDate")
        DataDateValue = valueStr;
      if (nameStr == "dataTime")
        DataTimeValue = valueStr;
      if (nameStr == "stepRange")
        StepRangeValue = valueStr;
    }
    grib_keys_iterator_delete(kiter);
    int stepRange;
    if (StepRangeValue == "0") {
      stepRange = 0;
    } else {
      std::vector<std::string> LStr = STRING_Split(StepRangeValue, "-");
      int siz = LStr.size();
      if (siz == 1) {
        stepRange = stoi(StepRangeValue);
      } else {
        if (siz != 2) {
          std::cerr << "Inconsistency in our assumptions\n";
          throw TerminalException{1};
        }
        stepRange = stoi(LStr[1]);
      }
    }
    //    std::cerr << "StepRangeValue=" << StepRangeValue << " stepRange=" <<
    //    stepRange << "\n";
    //
    //    std::cerr << "DataDateValue=" << DataDateValue << "\n";
    int year, month, day;
    if (DataDateValue == "unset") {
      long dataDate;
      GRIB_CHECK(grib_get_long(h, "dataDate", &dataDate), 0);
      //      std::cerr << "dataDate=" << dataDate << "\n";
      day = dataDate % 100;
      int res1 = (dataDate - day) / 100;
      month = res1 % 100;
      int res2 = (res1 - month) / 100;
      year = res2;
    } else {
      std::string yearStr = DataDateValue.substr(0, 4);
      //    std::cerr << "yearStr=" << yearStr << "\n";
      std::string monthStr = DataDateValue.substr(4, 2);
      //    std::cerr << "monthStr=" << monthStr << "\n";
      std::string dayStr = DataDateValue.substr(6, 2);
      //    std::cerr << "dayStr=" << dayStr << "\n";
      year = stoi(yearStr);
      month = stoi(monthStr);
      day = stoi(dayStr);
    }
    double shiftTime = GetShiftTime(eModelName);
    //
    if (DataTimeValue == "unset") {
      long dataTime;
      GRIB_CHECK(grib_get_long(h, "dataTime", &dataTime), 0);
      //      std::cerr << "dataTime=" << dataTime << "\n";
      if (dataTime == 0) {
        DataTimeValue = "0000";
      } else {
        DataTimeValue = std::to_string(dataTime);
      }
    }
    grib_handle_delete(h);
    //    std::cerr << "DataTimeValue=" << DataTimeValue << "\n";
    auto GetStringHourMin = [&]() -> std::pair<std::string, std::string> {
      int siz = DataTimeValue.size();
      if (siz == 4)
        return {DataTimeValue.substr(0, 2), DataTimeValue.substr(2, 2)};
      if (siz == 3)
        return {DataTimeValue.substr(0, 1), DataTimeValue.substr(1, 2)};
      if (siz == 1)
        return {std::string("0"), DataTimeValue};
      std::cerr << "Failed to find siz matching siz=" << siz << "\n";
      throw TerminalException{1};
    };
    std::pair<std::string, std::string> ePair = GetStringHourMin();
    std::string HourStr = ePair.first;
    std::string MinStr = ePair.second;
    //    std::cerr << "HourStr=" << HourStr << "\n";
    //    std::cerr << "MinStr=" << MinStr << "\n";
    int hour = stoi(HourStr);
    int min = stoi(MinStr);
    //
    double PreTimeStart = DATE_ConvertSix2mjd({year, month, day, hour, min, 0});
    double eTime = PreTimeStart +
                   static_cast<double>(stepRange) / static_cast<double>(24) +
                   shiftTime;
    double eTimeStart =
        ExtractTimeStartFromName(PreTimeStart, eTime, eModelName, FileName);
    std::string strDate = DATE_ConvertMjd2mystringPres(eTime);
    //    std::cerr << "FileName=" << FileName << " year=" << year << " month="
    //    << month << " day=" << day << " hour=" << hour << " min=" << min << "
    //    strDate=" << strDate << "\n";
    //
    //    std::cerr << "Inserting message FileName=" << FileName << "\n";
    //    std::cerr << "  idx=" << idx << "\n";
    //    std::cerr << "  ShortNameValue=" << ShortNameValue << "\n";
    //    std::cerr << "  eTime=" << eTime << "\n";
    GRIB_MessageInfo eInfo{
        ShortNameValue, cfVarName,  cfVarNameECMWF, NameValue, units, idx,
        eTime,          eTimeStart, stepRange,      FileName};
    ListInfo.push_back(eInfo);
    idx++;
  }
  fclose(in);
  return ListInfo;
}

MyMatrix<double> GRIB_ReadFromMessageInfo(GRIB_MessageInfo const &eMesg) {
  grib_handle *h = NULL;
  int err;
  FILE *in = NULL;
  in = fopen(eMesg.FileName.c_str(), "r");
  int idx = 0;
  while ((h = grib_handle_new_from_file(0, in, &err)) != NULL) {
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err, 0);
    if (idx == eMesg.idx) {
      long Ni, Nj, numberOfDataPoints;
      GRIB_CHECK(grib_get_long(h, "Ni", &Ni), 0);
      GRIB_CHECK(grib_get_long(h, "Nj", &Nj), 0);
      GRIB_CHECK(grib_get_long(h, "numberOfDataPoints", &numberOfDataPoints),
                 0);
      size_t size = numberOfDataPoints;
      std::vector<double> lats(size), lons(size), values(size);
      err = grib_get_data(h, lats.data(), lons.data(), values.data());
      grib_handle_delete(h);
      if (err != GRIB_SUCCESS)
        GRIB_CHECK(err, 0);
      int eta_rho = Ni;
      int xi_rho = Nj;
      MyMatrix<double> VAL(eta_rho, xi_rho);
      int idxPos = 0;
      for (int j = 0; j < xi_rho; j++)
        for (int i = 0; i < eta_rho; i++) {
          VAL(i, j) = values[idxPos];
          idxPos++;
        }
      fclose(in);
      return VAL;
    }
    grib_handle_delete(h);
    idx++;
  }
  fclose(in);
  //
  std::cerr << "Failed to find the matching GRIB_MessageInfo\n";
  throw TerminalException{1};
}

bool GRIB_HasVariable(std::vector<GRIB_MessageInfo> const &ListMessage,
                      std::string const &nameSearch,
                      std::string const &VarName) {
  for (auto &eMesg : ListMessage) {
    if (nameSearch == "shortName" && eMesg.shortName == VarName)
      return true;
    if (nameSearch == "cfVarName" && eMesg.cfVarName == VarName)
      return true;
    if (nameSearch == "cfVarNameECMWF" && eMesg.cfVarNameECMWF == VarName)
      return true;
  }
  return false;
}

MyMatrix<double>
GRIB_Read2Dvariable(std::vector<GRIB_MessageInfo> const &ListMessage,
                    std::string const &nameSearch, std::string const &VarName) {
  /*
  std::cerr << "|ListMessage|=" << ListMessage.size() << "\n";
  for (auto & eMesg : ListMessage) {
    std::cerr << "  eMesg : shortName = " << eMesg.shortName << " eMesg.FileName
  = " << eMesg.FileName << "\n";
    }*/
  for (auto &eMesg : ListMessage) {
    if (nameSearch == "shortName" && eMesg.shortName == VarName)
      return GRIB_ReadFromMessageInfo(eMesg);
    if (nameSearch == "cfVarName" && eMesg.cfVarName == VarName)
      return GRIB_ReadFromMessageInfo(eMesg);
    if (nameSearch == "cfVarNameECMWF" && eMesg.cfVarNameECMWF == VarName)
      return GRIB_ReadFromMessageInfo(eMesg);
  }
  std::cerr << "|ListMessage|=" << ListMessage.size() << "\n";
  for (auto &eMesg : ListMessage) {
    std::cerr << "Message:\n";
    std::cerr << "  eMesg.shortName = " << eMesg.shortName << "\n";
    std::cerr << "  eMesg.cfVarName = " << eMesg.cfVarName << "\n";
    std::cerr << "  eMesg.cfVarNameECMWF = " << eMesg.cfVarNameECMWF << "\n";
    std::cerr << "  eMesg.FileName = " << eMesg.FileName << "\n";
    std::cerr << "  eMesg.idx = " << eMesg.idx << "\n";
  }
  std::cerr << "Error in GRIB_Read2Dvariable\n";
  std::cerr << "Failed to find the variable =" << VarName << "\n";
  std::cerr << "Exiting\n";
  throw TerminalException{1};
}

MyMatrix<double>
GRID_Get2DVariableTimeDifferentiate(TotalArrGetData const &TotalArr,
                                    std::string const &VarName,
                                    double const &eTimeDay) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  double tolDay = static_cast<double>(1) / static_cast<double>(1000000);
  int nbTime = TotalArr.eArr.ListTime.size();
  if (nbTime > 1) {
    double deltTimeEst =
        (TotalArr.eArr.ListTime[nbTime - 1] - TotalArr.eArr.ListTime[0]) /
        static_cast<double>(nbTime - 1);
    tolDay = deltTimeEst / static_cast<double>(100);
  }
  //
  int nbTimeStart = TotalArr.eArr.ListStartTime.size();
  //  int nbTime=TotalArr.eArr.ListStartTime.size();
  //
  auto search = TotalArr.eArr.MatchingByVariable.find(VarName);
  if (search == TotalArr.eArr.MatchingByVariable.end()) {
    std::cerr << "Error in GRIB_Get2DvariableTimeDifferentiate\n";
    std::cerr << "The variable VarName = " << VarName << "\n";
    std::cerr << "is absent of the list of allowed variables\n";
    throw TerminalException{1};
  }
  //
  struct ShootSolution {
    int iMesgLow;
    int iMesgUpp;
    double eTimeLow;
    double eTimeUpp;
    double eStartTime;
    double DeltaTimeDay;
  };
  std::vector<ShootSolution> ListShootSolution;
  int TotalNbMessage = TotalArr.eArr.ListIStartTime.size();
#ifdef TIMINGS
  SingletonTime time2;
#endif
  std::vector<int> ListITimeStart;
  for (int iTimeStart = 0; iTimeStart < nbTimeStart; iTimeStart++) {
    double eStartTime = TotalArr.eArr.ListStartTime[iTimeStart];
    double eEndTime = TotalArr.eArr.ListEndTime[iTimeStart];
    if (eStartTime <= eTimeDay && eTimeDay <= eEndTime)
      ListITimeStart.push_back(iTimeStart);
  }
  //  std::cerr << "|ListITimeStart|=" << ListITimeStart.size() << "\n";
  for (int &iTimeStart : ListITimeStart) {
    double eStartTime = TotalArr.eArr.ListStartTime[iTimeStart];
    std::vector<int> ListIMesg;
    std::vector<double> ListTime;
    for (int iMesg : TotalArr.eArr.ListListIMesg[iTimeStart]) {
      if (TotalArr.eArr.ListAllMessage[iMesg].shortName == VarName) {
        int iTime = TotalArr.eArr.ListITime[iMesg];
        ListIMesg.push_back(iMesg);
        ListTime.push_back(TotalArr.eArr.ListTime[iTime]);
      }
    }
    int nbEnt = ListIMesg.size();
    for (int iEntUpp = 1; iEntUpp < nbEnt; iEntUpp++) {
      int iEntLow = iEntUpp - 1;
      int iMesgLow = ListIMesg[iEntLow];
      int iMesgUpp = ListIMesg[iEntUpp];
      int iTimeLow = TotalArr.eArr.ListITime[iMesgLow];
      int iTimeUpp = TotalArr.eArr.ListITime[iMesgUpp];
      double eTimeLow = TotalArr.eArr.ListTime[iTimeLow];
      double eTimeUpp = TotalArr.eArr.ListTime[iTimeUpp];
      if (eTimeLow - tolDay < eTimeDay && eTimeDay < eTimeUpp + tolDay) {
        double DeltaTimeDay = eTimeUpp - eTimeLow;
        ListShootSolution.push_back(
            {iMesgLow, iMesgUpp, eTimeLow, eTimeUpp, eStartTime, DeltaTimeDay});
      }
    }
  }
#ifdef TIMINGS
  SingletonTime time3;
#endif
  if (ListShootSolution.size() == 0) {
    std::cerr << "Printing debug information\n";
    std::vector<int> ListNBEnt;
    for (int iTimeStart = 0; iTimeStart < nbTimeStart; iTimeStart++) {
      double eStartTime = TotalArr.eArr.ListStartTime[iTimeStart];
      double eEndTime = TotalArr.eArr.ListEndTime[iTimeStart];
      std::string strStartTime = DATE_ConvertMjd2mystringPres(eStartTime);
      std::string strEndTime = DATE_ConvertMjd2mystringPres(eEndTime);
      std::vector<int> ListIMesg;
      std::vector<double> ListTime;
      std::string retString;
      for (int iMesg = 0; iMesg < TotalNbMessage; iMesg++) {
        if (TotalArr.eArr.ListIStartTime[iMesg] == iTimeStart &&
            TotalArr.eArr.ListAllMessage[iMesg].shortName == VarName) {
          ListIMesg.push_back(iMesg);
          int iTime = TotalArr.eArr.ListITime[iMesg];
          double eTime = TotalArr.eArr.ListTime[iTime];
          std::string strTime = DATE_ConvertMjd2mystringPres(eTime);
          retString += " " + std::to_string(iMesg) + ":" + strTime;
        }
      }
      std::cerr << iTimeStart << "/" << nbTimeStart
                << " start_time=" << strStartTime << " end_time=" << strEndTime
                << " list=" << retString << "\n";
      ListNBEnt.push_back(static_cast<int>(ListIMesg.size()));
    }
    CollectedResult<int> eColl = Collected(ListNBEnt);
    for (size_t u = 0; u < eColl.LVal.size(); u++) {
      int eVal = eColl.LVal[u];
      std::cerr << "u=" << u << " eVal=" << eVal << " eMult=" << eColl.LMult[u]
                << "\n";
    }
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "VarName=" << VarName << "\n";
    std::cerr << "eTimeDay=" << eTimeDay << "  strPres=" << strPres << " \n";
    std::cerr << "TotalNbMessage=" << TotalNbMessage
              << "  nbTimeStart=" << nbTimeStart << "\n";
    std::cerr
        << "|ListShootSolution| = 0 so we did not find any possible scenario\n";
    std::cerr << "for the differentiation\n";
    throw TerminalException{1};
  }
  ShootSolution eSol = ListShootSolution[0];
  for (auto &fSol : ListShootSolution) {
    double eDelta_Start = eTimeDay - eSol.eStartTime;
    double fDelta_Start = eTimeDay - fSol.eStartTime;
    if (fDelta_Start < eDelta_Start - tolDay ||
        (fabs(eDelta_Start - fDelta_Start) < tolDay &&
         fSol.DeltaTimeDay < eSol.DeltaTimeDay - tolDay))
      eSol = fSol;
  }
#ifdef TIMINGS
  SingletonTime time4;
#endif
  //
  MyMatrix<double> Flow =
      GRIB_ReadFromMessageInfo(TotalArr.eArr.ListAllMessage[eSol.iMesgLow]);
#ifdef TIMINGS
  SingletonTime time5;
#endif
  MyMatrix<double> Fupp =
      GRIB_ReadFromMessageInfo(TotalArr.eArr.ListAllMessage[eSol.iMesgUpp]);
#ifdef TIMINGS
  SingletonTime time6;
#endif
  double DeltaTimeSec = eSol.DeltaTimeDay * static_cast<double>(86400);
  std::cerr << "DeltaTimeSec=" << DeltaTimeSec << "\n";
  MyMatrix<double> Fret = (Fupp - Flow) / DeltaTimeSec;
#ifdef TIMINGS
  SingletonTime time7;
  std::cerr << "|Paperwork|=" << ms(time1, time2) << "\n";
  std::cerr << "|TimeStart Loop|=" << ms(time2, time3) << "\n";
  std::cerr << "|Selecting eSol|=" << ms(time3, time4) << "\n";
  std::cerr << "|Flow|=" << ms(time4, time5) << "\n";
  std::cerr << "|Fupp|=" << ms(time5, time6) << "\n";
  std::cerr << "|Fret|=" << ms(time6, time7) << "\n";
#endif
  return Fret;
}

MyMatrix<double> GRIB_Get2DvariableSpecTime(TotalArrGetData const &TotalArr,
                                            std::string const &nameSearch,
                                            std::string const &VarName,
                                            double const &eTimeDay) {
  //  std::cerr << "GRIB_Get2DvariableSpecTime, step 1 date=" <<
  //  DATE_ConvertMjd2mystringPres(eTimeDay) << "\n"; std::cerr <<
  //  "GRIB_Get2DvariableSpecTime, step 1 VarName=" << VarName << " nameSearch="
  //  << nameSearch << "\n";
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
#endif
  auto search = TotalArr.eArr.MatchingByVariable.find(VarName);
  if (search == TotalArr.eArr.MatchingByVariable.end()) {
    std::cerr << "Error in GRIB_Get2DvariableSpecTime\n";
    std::cerr << "The variable VarName = " << VarName << "\n";
    std::cerr << "is absent of the list of allowed variables\n";
    throw TerminalException{1};
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
#endif
  std::vector<int> const &RelListIndex =
      TotalArr.eArr.MatchingByVariable.at(VarName);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 =
      std::chrono::system_clock::now();
#endif
  int nbTime = RelListIndex.size();
  auto f = [&](int const &pos) -> double {
    return TotalArr.eArr.ListTime[RelListIndex[pos]];
  };
  InterpInfo eInterpInfo = GetTimeInterpolationInfo_F(nbTime, f, eTimeDay);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 =
      std::chrono::system_clock::now();
  std::cerr << "|MatchingByVariable1|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                     time1)
                   .count()
            << "\n";
  std::cerr << "|MatchingByVariable2|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time3 -
                                                                     time2)
                   .count()
            << "\n";
  std::cerr << "|eInterpInfo|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time4 -
                                                                     time3)
                   .count()
            << "\n";
  std::cerr << "UseSingleEntry=" << eInterpInfo.UseSingleEntry << "\n";
#endif
  if (eInterpInfo.UseSingleEntry) {
    int iTime = eInterpInfo.iTimeLow;
    int iTimeReal = RelListIndex[iTime];
    //    std::cerr << "GRIB_Get2DvariableSpecTime, step 2\n";
    return GRIB_Read2Dvariable(TotalArr.eArr.ListListMessages[iTimeReal],
                               nameSearch, VarName);
  }
  double alphaLow = eInterpInfo.alphaLow;
  int iTimeLow = eInterpInfo.iTimeLow;
  double alphaUpp = eInterpInfo.alphaUpp;
  int iTimeUpp = eInterpInfo.iTimeUpp;
  //
  int iTimeRealLow = RelListIndex[iTimeLow];
  int iTimeRealUpp = RelListIndex[iTimeUpp];
  //  std::cerr << "GRIB_Get2DvariableSpecTime, step 3\n";
  MyMatrix<double> eVarLow = GRIB_Read2Dvariable(
      TotalArr.eArr.ListListMessages[iTimeRealLow], nameSearch, VarName);
  //  std::cerr << "GRIB_Get2DvariableSpecTime, step 4\n";
  MyMatrix<double> eVarUpp = GRIB_Read2Dvariable(
      TotalArr.eArr.ListListMessages[iTimeRealUpp], nameSearch, VarName);
  //  std::cerr << "GRIB_Get2DvariableSpecTime, step 5\n";
  return alphaLow * eVarLow + alphaUpp * eVarUpp;
}

void PrintGribMessage(std::string const &FileName,
                      MyMatrix<double> const &TheField,
                      GRIB_MessageInfo const &eMesg, GridArray const &GrdArr) {
  FILE *of = NULL;
  long step = 0;
  grib_handle *h = NULL;
  grib_multi_handle *mh = NULL;
  const int start_section = 4; /* Grib2 Product Definition Section */

  /* create a new empty multi field handle */
  mh = grib_multi_handle_new(0);
  if (!mh) {
    std::cerr << "unable to create multi field handle\n";
    throw TerminalException{1};
  }

  grib_set_long(h, "step", step);
  grib_multi_handle_append(h, start_section, mh);
  if (IsExistingFile(FileName)) {
    of = fopen(FileName.c_str(), "a");
  } else {
    of = fopen(FileName.c_str(), "w");
  }
  if (!of) {
    std::cerr << "ERROR: unable to open file " << FileName << "\n";
    throw TerminalException{1};
  }
  grib_multi_handle_write(mh, of);
  fclose(of);
  grib_handle_delete(h);
}

// clang-format off
#endif  // SRC_OCEAN_BASIC_GRIB_H_
// clang-format on
