// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_SATELLITE_H_
#define SRC_OCEAN_SATELLITE_H_

#include "Basic_netcdf.h"
#include "Basic_file.h"
#include "Data_Access.h"
#include "Interpolation.h"
#include "Model_data_loading.h"
#include "Model_grids.h"
#include "NamelistExampleOcean.h"
#include "Plotting_fct.h"
#include "SphericalGeom.h"
#include "Statistics.h"
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

struct SingleEntryMeasurement {
  double Time, Lon, Lat;
  double WindSpeed, WindSpeed_cor, WindSpeed_used;
  double Sigma0, Sigma0_cal, Sigma0std, Sigma0second, Sigma0secondstd;
  double Swh, SwhStd, Swh_cor, Swh_used;
  std::vector<double> Swh_model, WindSpeed_model;
  double WaveDir;
  int Satellite;
  double DistToCoast; // not yet set
  double Altitude;
  double Ucurr, Vcurr;
  double QCcur, QCwave;
  double Footprint;
  double gradientLL;
};

SingleEntryMeasurement GetSingleEntryMeasurement() {
  double v = std::nan("1");
  return {v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          {},
          {},
          v,
          std::numeric_limits<int>::max(),
          v,
          v,
          v,
          v,
          v,
          v,
          v,
          v};
}

struct SatelliteSerInfo {
  double BeginTime;
  double EndTime;
  int FirstDayTime;
  int LastDayTime;
  std::string strBegin;
  std::string strEnd;
};

std::vector<std::string> GetAllNamesOfSatelliteAltimeter() {
  return {"ERS1", "ERS2",   "ENVISAT", "TOPEX", "POSEIDON",   "JASON1",
          "GFO",  "JASON2", "CRYOSAT", "SARAL", "RADAR_SPLIT", "JASON3", "SENTINEL_3A", "SENTINEL_3B"};
}

std::string GetNameOfSatelliteAltimeter(int iSat) {
  return GetAllNamesOfSatelliteAltimeter()[iSat - 1];
}

SatelliteSerInfo GetTimeSerInfo_From_BeginEnd(double const &BeginTime,
                                              double const &EndTime) {
  //
  int FirstDayTime = static_cast<int>(floor(BeginTime));
  int LastDayTime = static_cast<int>(ceil(EndTime) - 1);
  //
  std::string strBegin = DATE_ConvertMjd2mystringPres(BeginTime);
  std::string strEnd = DATE_ConvertMjd2mystringPres(EndTime);
  //
  return {BeginTime, EndTime, FirstDayTime, LastDayTime, strBegin, strEnd};
}

SatelliteSerInfo
RetrieveTimeInformation(std::vector<TotalArrGetData> const &ListTotalArr,
                        FullNamelist const &eFull) {
  int nbGrid = ListTotalArr.size();
  std::vector<double> ListFirstTime(nbGrid);
  std::vector<double> ListLastTime(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    ListFirstTime[iGrid] = MinimumTimeHistoryArray(ListTotalArr[iGrid].eArr);
    ListLastTime[iGrid] = MaximumTimeHistoryArray(ListTotalArr[iGrid].eArr);
  }
  //
  SingleBlock eBlSELECT = eFull.get_block("SELECT");
  std::string strBEGTC = eBlSELECT.get_string("BEGTC");
  std::string strENDTC = eBlSELECT.get_string("ENDTC");
  double BeginTime = 0, EndTime = 0;
  if (strBEGTC == "earliest") {
    BeginTime = VectorMax(ListFirstTime);
  } else {
    BeginTime = CT2MJD(strBEGTC);
  }
  if (strENDTC == "latest") {
    EndTime = VectorMin(ListLastTime);
  } else {
    EndTime = CT2MJD(strENDTC);
  }
  return GetTimeSerInfo_From_BeginEnd(BeginTime, EndTime);
}

SatelliteSerInfo RetrieveTimeInformation_Begin_End(FullNamelist const &eFull) {
  SingleBlock eBlSELECT = eFull.get_block("SELECT");
  std::string strBEGTC = eBlSELECT.get_string("BEGTC");
  std::string strENDTC = eBlSELECT.get_string("ENDTC");
  double BeginTime = CT2MJD(strBEGTC);
  double EndTime = CT2MJD(strENDTC);
  return GetTimeSerInfo_From_BeginEnd(BeginTime, EndTime);
}

//
//
//
std::vector<SingleEntryMeasurement>
AssignGradient(std::vector<SingleEntryMeasurement> const &ListEnt) {
  std::map<int, std::vector<SingleEntryMeasurement>> MapSatName;
  for (auto &eEnt : ListEnt) {
    auto iter = MapSatName.find(eEnt.Satellite);
    if (iter == MapSatName.end()) {
      MapSatName[eEnt.Satellite] = {eEnt};
    } else {
      MapSatName[eEnt.Satellite].push_back(eEnt);
    }
  }
  std::vector<SingleEntryMeasurement> ListRet;
  for (auto &ePair : MapSatName) {
    sort(ePair.second.begin(), ePair.second.end(),
         [](SingleEntryMeasurement const &eEnt1,
            SingleEntryMeasurement const &eEnt2) -> bool {
           if (eEnt1.Time < eEnt2.Time)
             return true;
           if (eEnt2.Time < eEnt1.Time)
             return false;
           std::cerr << "Error in Comparison. Dates seem to be identical for "
                        "the same satellite. Impossible\n";
           throw TerminalException{1};
         });
    int len = ePair.second.size();
    std::vector<double> ListDeltaTime(len - 1);
    std::vector<double> ListDistKM(len - 1);
    for (int i = 0; i < len - 1; i++) {
      double DeltaTime = ePair.second[i + 1].Time - ePair.second[i].Time;
      ListDeltaTime[i] = DeltaTime;
      double eLon1 = ePair.second[i].Lon;
      double eLat1 = ePair.second[i].Lat;
      double eLon2 = ePair.second[i + 1].Lon;
      double eLat2 = ePair.second[i + 1].Lat;
      double eDistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
      ListDistKM[i] = eDistKM;
    }
    if (len > 1) {
      double FactorMult = 4.0;
      double minDeltaTime = VectorMin(ListDeltaTime);
      double minDistKM = VectorMin(ListDistKM);
      std::vector<bool> ListStatus(len - 1, true);
      for (int i = 0; i < len - 1; i++) {
        bool test1 = ListDeltaTime[i] > (minDeltaTime * FactorMult);
        bool test2 = ListDistKM[i] > (minDistKM * FactorMult);
        if (test1 != test2) {
          std::cerr << "We have difference in conclusion at i=" << i
                    << " test1=" << test1 << " test2=" << test2 << "\n";
        }
        if (test1 || test2)
          ListStatus[i] = false;
      }
      std::vector<double> ListGrad(len - 1);
      for (int i = 0; i < len - 1; i++) {
        double eLon1 = ePair.second[i].Lon;
        double eLat1 = ePair.second[i].Lat;
        double eLon2 = ePair.second[i + 1].Lon;
        double eLat2 = ePair.second[i + 1].Lat;
        double eGrad = (eLon1 - eLon2) / (eLat1 - eLat2);
        ListGrad[i] = eGrad;
      }
      for (int i = 0; i < len; i++) {
        double sum = 0;
        int sumI = 0;
        if (i < len - 1) {
          if (ListStatus[i]) {
            sum += ListGrad[i];
            sumI++;
          }
        }
        if (i > 0) {
          if (ListStatus[i - 1]) {
            sum += ListGrad[i - 1];
            sumI++;
          }
        }
        if (sumI > 0) {
          double eGrad = sum / static_cast<double>(sumI);
          ePair.second[i].gradientLL = eGrad;
        }
      }
    }
    ListRet.insert(ListRet.end(), ePair.second.begin(), ePair.second.end());
  }
  return ListRet;
}

std::vector<SingleEntryMeasurement>
READ_ALTI_FILE_IFREMER(std::string const &eFileAlti) {
  std::cerr << "Read Altimetry file : " << eFileAlti << "\n";
  std::vector<double> ListTime = NC_ReadTimeFromFile(eFileAlti, "time");
  // other variables
  MyVector<double> ListLon = NC_Read1Dvariable(eFileAlti, "lon");
  MyVector<double> ListLat = NC_Read1Dvariable(eFileAlti, "lat");
  MyVector<double> ListWindSpeed = NC_Read1Dvariable(eFileAlti, "wind_speed");
  MyVector<double> ListWindSpeed_cor =
      NC_Read1Dvariable(eFileAlti, "wind_speed_cor");
  MyVector<double> ListSigma0 = NC_Read1Dvariable(eFileAlti, "sigma0");
  MyVector<double> ListSigma0cal = NC_Read1Dvariable(eFileAlti, "sigma0_cal");
  MyVector<double> ListSigma0std = NC_Read1Dvariable(eFileAlti, "sigma0std");
  MyVector<double> ListSigma0second =
      NC_Read1Dvariable(eFileAlti, "sigma0second");
  MyVector<double> ListSigma0secondstd =
      NC_Read1Dvariable(eFileAlti, "sigma0secondstd");
  MyVector<double> ListSwh = NC_Read1Dvariable(eFileAlti, "swh");
  MyVector<double> ListSwhstd = NC_Read1Dvariable(eFileAlti, "swhstd");
  MyVector<double> ListSwhcor = NC_Read1Dvariable(eFileAlti, "swhcor");
  MyVector<int> ListSatellite = NC_Read1Dvariable_int(eFileAlti, "satellite");
  //
  int nbMeas = ListTime.size();
  std::vector<SingleEntryMeasurement> ListEnt;
  for (int iMeas = 0; iMeas < nbMeas; iMeas++) {
    SingleEntryMeasurement eEnt = GetSingleEntryMeasurement();
    eEnt.Satellite = ListSatellite[iMeas];
    eEnt.Time = ListTime[iMeas];
    eEnt.Lon = ListLon(iMeas);
    eEnt.Lat = ListLat(iMeas);
    eEnt.WindSpeed = ListWindSpeed(iMeas);
    eEnt.WindSpeed_cor = ListWindSpeed_cor(iMeas);
    eEnt.Sigma0 = ListSigma0(iMeas);
    eEnt.Sigma0_cal = ListSigma0cal(iMeas);
    eEnt.Sigma0std = ListSigma0std(iMeas);
    eEnt.Sigma0second = ListSigma0second(iMeas);
    eEnt.Sigma0secondstd = ListSigma0secondstd(iMeas);
    eEnt.Swh = ListSwh(iMeas);
    eEnt.SwhStd = ListSwhstd(iMeas);
    eEnt.Swh_cor = ListSwhcor(iMeas);
    // unset values
    eEnt.Swh_used = -1;
    eEnt.WindSpeed_used = -1;
    eEnt.Footprint = -1;
    eEnt.QCwave = 255;
    eEnt.gradientLL = std::nan("1");
    ListEnt.push_back(eEnt);
  }
  return ListEnt;
}

std::vector<SingleEntryMeasurement>
READ_ALTI_FILE_EUMETCAST_SINGLE(std::string const &eFileAlti) {
  std::cerr << "Eumetcast, Read Altimetry file : " << eFileAlti << "\n";
  std::vector<std::string> F_LStr = STRING_Split(eFileAlti, "/");
  std::string eFileNaked = F_LStr[F_LStr.size() - 1];
  std::vector<std::string> F_LStr2 = STRING_Split(eFileNaked, "_");
  std::string SatName = F_LStr2[0];
  std::string SwhString, WindSpeedString, Sigma0String;
  std::string eNameCanonical;
  if (SatName == "JA2") {
    eNameCanonical = "JASON2";
    SwhString = "swh_ku_mle3";
    WindSpeedString = "wind_speed_alt_mle3";
    Sigma0String = "sig0_ku_mle3";
  }
  if (SatName == "SRL") {
    eNameCanonical = "SARAL";
    SwhString = "swh";
    WindSpeedString = "wind_speed_alt";
    Sigma0String = "sig0";
  }
  std::vector<double> ListTime = NC_ReadTimeFromFile(eFileAlti, "time");
  // other variables
  MyVector<double> ListLon = NC_Read1Dvariable(eFileAlti, "lon");
  MyVector<double> ListLat = NC_Read1Dvariable(eFileAlti, "lat");
  MyVector<double> ListWindSpeed =
      NC_Read1Dvariable(eFileAlti, WindSpeedString);
  MyVector<double> ListSwh = NC_Read1Dvariable(eFileAlti, SwhString);
  MyVector<double> ListSigma0 = NC_Read1Dvariable(eFileAlti, Sigma0String);
  MyVector<double> ListAltitude =
      NC_Read1Dvariable(eFileAlti, "mean_sea_surface");
  //
  // Determining the satellite identifier
  //
  int eSatellite = -1;
  std::vector<std::string> ListName = GetAllNamesOfSatelliteAltimeter();
  for (int iSat = 0; iSat < static_cast<int>(ListName.size()); iSat++)
    if (ListName[iSat] == eNameCanonical)
      eSatellite = iSat + 1;
  if (eSatellite == -1) {
    std::cerr << "We have SatName = " << SatName << "\n";
    std::cerr << "But the only allowed ones are JA2 and SRL\n";
    throw TerminalException{1};
  }
  //
  // Building return array
  //
  int nbMeas = ListTime.size();
  std::vector<SingleEntryMeasurement> ListEnt(nbMeas);
  for (int iMeas = 0; iMeas < nbMeas; iMeas++) {
    SingleEntryMeasurement eEnt = GetSingleEntryMeasurement();
    eEnt.Time = ListTime[iMeas];
    eEnt.Lon = ListLon(iMeas);
    eEnt.Lat = ListLat(iMeas);
    eEnt.WindSpeed = ListWindSpeed(iMeas);
    eEnt.WindSpeed_cor = ListWindSpeed(iMeas);
    eEnt.Sigma0 = ListSigma0(iMeas);
    eEnt.Swh = ListSwh(iMeas);
    eEnt.Swh_cor = ListSwh(iMeas);
    eEnt.Altitude = ListAltitude(iMeas);
    eEnt.Satellite = eSatellite;
    eEnt.gradientLL = std::nan("1");
    ListEnt[iMeas] = eEnt;
  }
  return ListEnt;
}

std::vector<SingleEntryMeasurement>
READ_ALTI_FILE_EUMETCAST(int const &year, int const &month, int const &day,
                         std::string const &eDirData) {
  std::string FullDir =
      eDirData + StringNumber(year, 4) + "/" + StringNumber(month, 2) + "/";
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(FullDir);
  std::vector<SingleEntryMeasurement> ListEntReturn;
  std::string eStrDate =
      StringNumber(year, 4) + StringNumber(month, 2) + StringNumber(day, 2);
  for (auto &eFile : ListFile) {
    //    std::vector<std::string> LStr=STRING_Split(eFile, "/");
    //    std::string eFileNaked=LStr[LStr.size()-1];
    std::vector<std::string> LStr2 = STRING_Split(eFile, "_");
    if (LStr2[1] == "OPR" && LStr2[4] == eStrDate) {
      std::string FullFile = FullDir + eFile;
      std::vector<SingleEntryMeasurement> eV =
          READ_ALTI_FILE_EUMETCAST_SINGLE(FullFile);
      ListEntReturn.insert(ListEntReturn.end(), eV.begin(), eV.end());
    }
  }
  return ListEntReturn;
}

struct SingleSearchEntry {
  size_t iEntry;
  int eEta;
  int eXi;
  double TimeCoeff;
  double SpatialCoeff;
};

std::vector<TotalArrGetData> RealAllArrayHistory(FullNamelist const &eFull) {
  SingleBlock eBlPROC = eFull.get_block("PROC");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  int nbGrid = ListGridFile.size();
  std::vector<TotalArrGetData> ListTotalArr(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    std::string eModelName = ListModelName[iGrid];
    std::string GridFile = ListGridFile[iGrid];
    std::string HisPrefix = ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    ListTotalArr[iGrid] = RetrieveTotalArr(eTriple);
  }
  return ListTotalArr;
}

// The format of SingleRecInterp is just too nice not to use it
// for all methods.
//
std::vector<SingleRecInterp> Function_Interpolation_Nearest_Etc(
    GridArray const &GrdArr, MyMatrix<double> const &ListXY,
    MyVector<double> const &ListRadius, std::string const &eMethod) {
  if (eMethod == "linear_interpolation") {
    return General_FindInterpolationWeight(GrdArr, ListXY, false);
  }
  if (eMethod == "nearest") {
#ifdef USE_OPENCV_LIBRARY
    return NearestInterpolation_FindWeight(GrdArr, ListXY);
#else
    std::cerr
        << "Need to compile with opencv in order to have the nearest option\n";
    throw TerminalException{1};
#endif
  }

  std::cerr << "Function_Interpolation_Nearest_Etc accepts as method\n";
  std::cerr << "linear_interpolation, nearest, \n";
  throw TerminalException{1};
}

std::vector<double>
ConvertListStringValueToVector_SAT(std::vector<std::string> const &ListStr,
                                   double const &eDefault) {
  std::vector<std::string> ListSatName = GetAllNamesOfSatelliteAltimeter();
  int nbSat = ListSatName.size();
  // no assignment put default
  if (ListStr.size() == 0) {
    return std::vector<double>(nbSat, eDefault);
  }
  // Case of a blanket value for everything
  if (ListStr.size() == 1) {
    if (IsFullyNumeric(ListStr[0])) {
      double eVal;
      std::istringstream(ListStr[0]) >> eVal;
      //      std::cerr << "eVal=" << eVal << "\n";
      return std::vector<double>(nbSat, eVal);
    }
  }
  //  std::cerr << "nbSat=" << nbSat << " |ListStr|=" << ListStr.size() << "\n";
  double eNAN = std::nan("1");
  std::vector<double> eVect(nbSat, eNAN);
  for (auto &ePair : ListStr) {
    std::vector<std::string> LSpl = STRING_Split(ePair, ":");
    if (LSpl.size() != 2) {
      std::cerr << "The splitting of ePair by : failed\n";
      std::cerr
          << "Error is in ListMinHS_meas (or ListMaxHS_meas or similar)\n";
      std::cerr << "Legal format are for example\n";
      std::cerr << "ListMinHS_meas = \"CRYOSAT:0.121\", \"JASON2:0\",\n";
      std::cerr << "or\n";
      std::cerr << "ListMinHS_meas = \"0\",";
      throw TerminalException{1};
    }
    bool WeAss = false;
    for (int iSat = 0; iSat < nbSat; iSat++) {
      if (ListSatName[iSat] == LSpl[0]) {
        double eVal;
        std::istringstream(LSpl[1]) >> eVal;
        eVect[iSat] = eVal;
        WeAss = true;
      }
    }
    if (!WeAss) {
      std::cerr << "String is not recognized\n";
      std::cerr << "LSpl[0]=" << LSpl[0] << "\n";
      throw TerminalException{1};
    }
  }
  return eVect;
}

void InterpolateAltimeterData(std::vector<SingleEntryMeasurement> &ListEntry,
                              std::vector<TotalArrGetData> const &ListTotalArr,
                              FullNamelist const &eFull) {
  std::cerr << "InterpolateAltimeterData, step 1\n";
  SingleBlock eBlPROC = eFull.get_block("PROC");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  //
  SingleBlock eBlSELECT = eFull.get_block("SELECT");
  std::vector<std::string> ListFootprintKM_str =
    eBlSELECT.get_list_string("ListRadiusFootprintKM");
  std::vector<double> ListFootprintKM = ConvertListStringValueToVector_SAT(
      ListFootprintKM_str, 40); // 40km footprint
  std::string TheMethodInterp =
    eBlSELECT.get_string("InterpolationMethod");
  //
  int nbGrid = ListGridFile.size();
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    std::string eModelName = ListModelName[iGrid];
    std::string GridFile = ListGridFile[iGrid];
    std::string HisPrefix = ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
    int nbRow = GrdArr.GrdArrRho.MSK.rows();
    int nbCol = GrdArr.GrdArrRho.MSK.cols();
    int sumCoeff = GrdArr.GrdArrRho.MSK.sum();
    int deltaMSK = nbRow * nbCol - sumCoeff;
    std::cerr << "LON(min/max)=" << GrdArr.GrdArrRho.LON.minCoeff() << " / "
              << GrdArr.GrdArrRho.LON.maxCoeff() << "\n";
    std::cerr << "LAT(min/max)=" << GrdArr.GrdArrRho.LAT.minCoeff() << " / "
              << GrdArr.GrdArrRho.LAT.maxCoeff() << "\n";
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol
              << " sum(MSK)=" << sumCoeff << " deltaMSK=" << deltaMSK << "\n";
    //    std::cerr << "Before call to ReadArrayHistory B, iGrid=" << iGrid <<
    //    "\n"; ArrayHistory eArr=ReadArrayHistory(eTriple);
    std::cerr << "Before call to ARR_PrintHistoryArray\n";
    ARR_PrintHistoryArray(std::cerr, ListTotalArr[iGrid].eArr);
    std::cerr << "After call to ARR_PrintHistoryArray\n";
    std::cerr << "InterpolateAltimeterData, step 2\n";
    //
    // Since it is about altimeter data, we have to select first the needed
    // indexes and the needed
    //
    size_t nbEntry = ListEntry.size();
    MyMatrix<double> ListXY(2, nbEntry);
    MyVector<double> ListRadius(nbEntry);
    std::vector<double> ListTime(nbEntry);
    std::vector<int> ListStatusTime(nbEntry, 0);
    std::vector<double> ListLON(nbEntry), ListLAT(nbEntry);
    for (size_t iEntry = 0; iEntry < nbEntry; iEntry++) {
      double eLon = ListEntry[iEntry].Lon;
      double eLat = ListEntry[iEntry].Lat;
      double eFootprint = ListEntry[iEntry].Footprint;
      ListXY(0, iEntry) = eLon;
      ListXY(1, iEntry) = eLat;
      ListRadius(iEntry) = eFootprint;
      ListLON[iEntry] = eLon;
      ListLAT[iEntry] = eLat;
    }
    std::cerr << "ListXY, LON(min/max)=" << VectorMin(ListLON) << " / "
              << VectorMax(ListLON) << " nbEntry=" << nbEntry << "\n";
    std::cerr << "ListXY, LAT(min/max)=" << VectorMin(ListLAT) << " / "
              << VectorMax(ListLAT) << "\n";
    std::cerr << "InterpolateAltimeterData, step 3\n";
    int nbTime = 0;
    for (size_t iEntry = 0; iEntry < nbEntry; iEntry++) {
      double eTime = ListEntry[iEntry].Time;
      InterpInfo eInterpInfo =
          GetTimeInterpolationInfoGeneralized(ListTotalArr[iGrid].eArr, eTime);
      int iTimeLow = eInterpInfo.iTimeLow;
      if (iTimeLow >= nbTime) {
        nbTime = iTimeLow + 1;
      }
      int iTimeUpp = eInterpInfo.iTimeUpp;
      if (iTimeUpp >= nbTime) {
        nbTime = iTimeUpp + 1;
      }
    }
    std::cerr << "InterpolateAltimeterData, step 4\n";
    std::vector<SingleRecInterp> LRec = Function_Interpolation_Nearest_Etc(
        GrdArr, ListXY, ListRadius, TheMethodInterp);
    std::cerr << "InterpolateAltimeterData, step 5 |LRec|=" << LRec.size()
              << "\n";
    std::vector<std::vector<SingleSearchEntry>> ListListCases(nbTime);
    std::cerr << "InterpolateAltimeterData, step 5.1 nbEntry=" << nbEntry
              << "\n";
    for (size_t iEntry = 0; iEntry < nbEntry; iEntry++) {
      double eTime = ListEntry[iEntry].Time;
      ListTime[iEntry] = eTime;
      InterpInfo eInterpInfo =
          GetTimeInterpolationInfoGeneralized(ListTotalArr[iGrid].eArr, eTime);
      SingleRecInterp eRec = LRec[iEntry];
      int iTimeLow = eInterpInfo.iTimeLow;
      int iTimeUpp = eInterpInfo.iTimeUpp;
      double alphaLow = eInterpInfo.alphaLow;
      double alphaUpp = eInterpInfo.alphaUpp;
      //
      double eValAssign = -1000;
      if (eRec.status == 1) {
        eValAssign = 0;
        for (auto &ePart : eRec.LPart) {
          int eEta = ePart.eEta;
          int eXi = ePart.eXi;
          double eSpatialCoeff = ePart.eCoeff;
          //
          SingleSearchEntry eEnt = {iEntry, eEta, eXi, alphaLow, eSpatialCoeff};
          SingleSearchEntry fEnt = {iEntry, eEta, eXi, alphaUpp, eSpatialCoeff};
          ListListCases[iTimeLow].push_back(eEnt);
          ListListCases[iTimeUpp].push_back(fEnt);
        }
      }
      ListEntry[iEntry].Swh_model.push_back(eValAssign);
      ListEntry[iEntry].WindSpeed_model.push_back(eValAssign);
    }
    std::cerr << "InterpolateAltimeterData, step 6\n";
    std::vector<int> ListNbMatch(nbEntry, 0);
    std::vector<double> ListSumWeight_wnd(nbEntry, 0);
    std::vector<double> ListSumWeight_swh(nbEntry, 0);
    SingleBlock const& BlockPROCESS = eFull.get_block("PROCESS");
    bool DoWnd = BlockPROCESS.get_bool("DO_WNDMAG");
    bool DoSwh = BlockPROCESS.get_bool("DO_HS");
    for (int iTime = 0; iTime < nbTime; iTime++) {
      int nbCase = ListListCases[iTime].size();
      double eTimeDay = ARR_GetTime(ListTotalArr[iGrid].eArr, iTime);
      if (nbCase > 0) {
        if (DoWnd) {
          RecVar eRecVar = ModelSpecificVarSpecificTime(ListTotalArr[iGrid],
                                                        "WINDMAG", eTimeDay);
          for (auto &eSingEntry : ListListCases[iTime]) {
            int iEntry = eSingEntry.iEntry;
            int eEta = eSingEntry.eEta;
            int eXi = eSingEntry.eXi;
            double TimeCoeff = eSingEntry.TimeCoeff;
            double SpatialCoeff = eSingEntry.SpatialCoeff;
            ListEntry[iEntry].WindSpeed_model[iGrid] +=
                TimeCoeff * SpatialCoeff * eRecVar.F(eEta, eXi);
            ListSumWeight_wnd[iEntry] += TimeCoeff * SpatialCoeff;
          }
        }
        if (DoSwh) {
          RecVar eRecVar = ModelSpecificVarSpecificTime(ListTotalArr[iGrid],
                                                        "Hwave", eTimeDay);
          for (auto &eSingEntry : ListListCases[iTime]) {
            int iEntry = eSingEntry.iEntry;
            int eEta = eSingEntry.eEta;
            int eXi = eSingEntry.eXi;
            double TimeCoeff = eSingEntry.TimeCoeff;
            double SpatialCoeff = eSingEntry.SpatialCoeff;
            ListEntry[iEntry].Swh_model[iGrid] +=
                TimeCoeff * SpatialCoeff * eRecVar.F(eEta, eXi);
            ListNbMatch[iEntry]++;
            ListSumWeight_swh[iEntry] += TimeCoeff * SpatialCoeff;
          }
        }
      }
    }
    double SumErrorWeight = 0;
    if (DoSwh) {
      for (size_t iEntry = 0; iEntry < nbEntry; iEntry++) {
        //    std::cerr << "iEntry=" << iEntry << " Swh_model=" <<
        //    ListEntry[iEntry].Swh_model << "\n";
        double eSumWeight = ListSumWeight_swh[iEntry];
        //    std::cerr << "  nbMatch=" << ListNbMatch[iEntry] << " SumWeight="
        //    << eSumWeight << "\n";
        SumErrorWeight += fabs(eSumWeight - static_cast<double>(1));
      }
    }
    if (DoWnd) {
      for (size_t iEntry = 0; iEntry < nbEntry; iEntry++) {
        //    std::cerr << "iEntry=" << iEntry << " Swh_model=" <<
        //    ListEntry[iEntry].Swh_model << "\n";
        double eSumWeight = ListSumWeight_wnd[iEntry];
        //    std::cerr << "  nbMatch=" << ListNbMatch[iEntry] << " SumWeight="
        //    << eSumWeight << "\n";
        SumErrorWeight += fabs(eSumWeight - static_cast<double>(1));
      }
    }
    std::cerr << "SumErrorWeight = " << SumErrorWeight << "\n";
    std::cerr << "InterpolateAltimeterData, step 7\n";
  }
}

void CheckSatelliteList(FullNamelist const &eFull) {
  SingleBlock const& BlockSELECT = eFull.get_block("SELECT");
  std::vector<std::string> AllSatNames = GetAllNamesOfSatelliteAltimeter();
  std::vector<std::string> AllowedSatNames = BlockSELECT.get_list_string("AllowedSatellites");
  for (auto &eSatName : AllowedSatNames) {
    std::cerr << "eSatName=" << eSatName << "\n";
    bool WeMatch = false;
    for (auto &fSatName : AllSatNames)
      if (fSatName == eSatName)
        WeMatch = true;
    if (!WeMatch) {
      std::cerr << "Error in input.\n";
      std::cerr << "AllowedSatellites was assigned to:\n";
      WriteStdVector(std::cerr, AllowedSatNames);
      std::cerr << "On the other the only satellites that you can put are "
                   "following:\n";
      WriteStdVector(std::cerr, AllSatNames);
      throw TerminalException{1};
    }
  }
}

std::vector<int> GetListSatelliteId_vect(FullNamelist const &eFull) {
  SingleBlock const& BlockSELECT = eFull.get_block("SELECT");
  std::set<int> SatelliteId;
  CheckSatelliteList(eFull);
  std::vector<std::string> AllSatNames = GetAllNamesOfSatelliteAltimeter();
  int nbSatellite = AllSatNames.size();
  std::vector<int> ListStatusSatellite(nbSatellite, 0);
  std::vector<std::string> AllowedSatNames = BlockSELECT.get_list_string("AllowedSatellites");
  for (auto &eSatName : AllowedSatNames)
    for (int iSat = 0; iSat < nbSatellite; iSat++)
      if (AllSatNames[iSat] == eSatName)
        ListStatusSatellite[iSat] = 1;
  return ListStatusSatellite;
}

std::vector<SingleEntryMeasurement>
READ_RADAR_CSV_FILE(int const &year, int const &month, int const &day,
                    std::string const &eDirData) {
  //  std::cerr << "READ_RADAR_CSV_FILE, step 1\n";
  std::string strYear = StringNumber(year, 4);
  std::string strMonth = StringNumber(month, 2);
  std::string strDay = StringNumber(day, 2);
  std::string TheFileComplete = eDirData + strYear + "/" + strMonth + "/wm_" +
                                strYear + strMonth + strDay + ".csv";
  //  std::cerr << "READ_RADAR_CSV_FILE, step 2\n";
  std::vector<OneArrayRadarMeas> ListRadarMeas = ReadRadarFile(TheFileComplete);
  //  std::cerr << "READ_RADAR_CSV_FILE, step 3\n";
  std::vector<SingleEntryMeasurement> ListEnt;
  int idRadarSplit = -1;
  std::vector<std::string> ListStr = GetAllNamesOfSatelliteAltimeter();
  for (int iPos = 0; iPos < static_cast<int>(ListStr.size()); iPos++) {
    if (ListStr[iPos] == "RADAR_SPLIT") {
      idRadarSplit = iPos + 1;
    }
  }
  std::vector<double> ListQCWave;
  for (auto &eRec : ListRadarMeas) {
    double date = eRec.date;
    for (auto &eEnt : eRec.ListMeas) {
      SingleEntryMeasurement NewEnt = GetSingleEntryMeasurement();
      NewEnt.Lon = eEnt.lon;
      NewEnt.Lat = eEnt.lat;
      NewEnt.Time = date;
      NewEnt.QCcur = eEnt.QCcur;
      NewEnt.QCwave = eEnt.QCwave;
      NewEnt.WaveDir = eEnt.Wdir;
      NewEnt.Swh = eEnt.Wheight;
      NewEnt.Swh_cor = eEnt.Wheight;
      NewEnt.Ucurr = eEnt.U;
      NewEnt.Vcurr = eEnt.V;
      NewEnt.Satellite = idRadarSplit;
      NewEnt.gradientLL = std::nan("1");
      ListEnt.push_back(NewEnt);
      ListQCWave.push_back(eEnt.QCwave);
    }
  }
  std::cerr << "READ_RADAR_CSV_FILE, step 4 |ListEnt|=" << ListEnt.size()
            << "\n";
  if (ListQCWave.size() > 0) {
    std::cerr << "ListQCWave(min/max) = " << VectorMin(ListQCWave) << " / "
              << VectorMax(ListQCWave) << "\n";
  }
  return ListEnt;
}

std::vector<SingleEntryMeasurement>
RETRIEVE_RELEVANT_ALTI_DATA(SatelliteSerInfo const &eRecSer,
                            FullNamelist const &eFull) {
  SingleBlock const& BlockPROC = eFull.get_block("PROC");
  std::vector<std::string> ListTypeData = BlockPROC.get_list_string("ListTypeData");
  std::vector<std::string> ListDirData = BlockPROC.get_list_string("ListDirData");
  if (ListTypeData.size() != ListDirData.size()) {
    std::cerr << "|ListTypeData| = " << ListTypeData.size() << "\n";
    std::cerr << "|ListDirData|  = " << ListDirData.size() << "\n";
    std::cerr << "ListTypeData has not the same size as ListDirData\n";
    throw TerminalException{1};
  }
  if (ListTypeData.size() == 0) {
    std::string ePrefixAlti = STRING_GETENV("ALTIMETER_DIRECTORY");
    ListTypeData.push_back("IFREMER");
    ListDirData.push_back(ePrefixAlti);
  }
  int nbType = ListTypeData.size();
  //
  // Now the selection information
  //
  SingleBlock eBlSELECT = eFull.get_block("SELECT");
  double MinWind = eBlSELECT.get_double("MinWIND");
  double MaxWind = eBlSELECT.get_double("MaxWIND");
  std::vector<std::string> ListMinHs_measStr =
    eBlSELECT.get_list_string("ListMinHS_meas");
  std::vector<double> ListMinHs_meas =
      ConvertListStringValueToVector_SAT(ListMinHs_measStr, 0);
  //
  std::vector<std::string> ListMaxHs_measStr =
    eBlSELECT.get_list_string("ListMaxHS_meas");
  std::vector<double> ListMaxHs_meas =
      ConvertListStringValueToVector_SAT(ListMaxHs_measStr, 998);
  //
  std::vector<std::string> ListFootprintKM_str =
    eBlSELECT.get_list_string("ListRadiusFootprintKM");
  std::vector<double> ListFootprintKM = ConvertListStringValueToVector_SAT(
      ListFootprintKM_str, 40); // 40km footprint
  //
  int GEOSELECTION = eBlSELECT.get_int("GEOSELECTION");
  double MinLON = eBlSELECT.get_double("MinLON");
  double MaxLON = eBlSELECT.get_double("MaxLON");
  double MinLAT = eBlSELECT.get_double("MinLAT");
  double MaxLAT = eBlSELECT.get_double("MaxLAT");
  double minQCwave = eBlSELECT.get_double("minQCwave");
  double maxQCwave = eBlSELECT.get_double("maxQCwave");
  std::vector<double> LonPoly = eBlSELECT.get_list_double("LONPOLY");
  std::vector<double> LatPoly = eBlSELECT.get_list_double("LATPOLY");
  std::vector<int> ListStatusSatellite = GetListSatelliteId_vect(eFull);
  SingleBlock eBlPROCESS = eFull.get_block("PROCESS");
  bool USE_CORRECTED = eBlPROCESS.get_bool("USE_CORRECTED");
  bool SelectingHours = eBlSELECT.get_bool("SelectingHours");
  std::vector<int> ListAllowedHours =
    eBlSELECT.get_list_int("ListAllowedHours");
  double LargeValue = 1.0e31;
  //
  // The insertion stuff.
  //
  std::vector<SingleEntryMeasurement> eRetList;
  auto FuncInsertEntry = [&](SingleEntryMeasurement const &eEnt) -> bool {
    double eLon = eEnt.Lon;
    double eLat = eEnt.Lat;
    //
    // Geographical selection
    //
    bool eStatusGeo = true;
    if (GEOSELECTION == 1) {
      if (!(eLon >= MinLON && eLon <= MaxLON && eLat >= MinLAT &&
            eLat <= MaxLAT))
        eStatusGeo = false;
    }
    if (GEOSELECTION == 2) {
      bool IsInside = IsPointInside(eLon, eLat, LonPoly, LatPoly);
      if (!IsInside)
        eStatusGeo = false;
    }
    if (!eStatusGeo)
      return false;
    if (eEnt.QCwave < minQCwave || eEnt.QCwave > maxQCwave)
      return false;
    //
    // Time selection
    //
    double eTime = eEnt.Time;
    bool eStatusTime = false;
    if (eTime >= eRecSer.BeginTime && eTime <= eRecSer.EndTime)
      eStatusTime = true;
    if (!eStatusTime)
      return false;
    //
    // selection by allowed satellites
    //
    int eSatellite = eEnt.Satellite;
    if (ListStatusSatellite[eSatellite - 1] == 0)
      return false;
    //
    // Now copying and specifying values for WIND
    //
    SingleEntryMeasurement fEnt = eEnt;
    double WindSpeed_used;
    if (USE_CORRECTED)
      WindSpeed_used = eEnt.WindSpeed_cor;
    else
      WindSpeed_used = eEnt.WindSpeed;
    if (WindSpeed_used > MaxWind || WindSpeed_used < MinWind)
      WindSpeed_used = LargeValue;
    fEnt.WindSpeed_used = WindSpeed_used;
    //    std::cerr << "FuncInsertEntry WindSpeed_used = " << WindSpeed_used <<
    //    "\n";
    //
    // Now copying and specifying values for WIND
    //
    double Swh_used;
    if (USE_CORRECTED)
      Swh_used = eEnt.Swh_cor;
    else
      Swh_used = eEnt.Swh;
    double eMinHs_meas = ListMinHs_meas[eSatellite - 1];
    double eMaxHs_meas = ListMaxHs_meas[eSatellite - 1];
    if (std::isnan(eMinHs_meas) || std::isnan(eMaxHs_meas)) {
      std::cerr << "eMinHs_meas=" << eMinHs_meas << "\n";
      std::cerr << "eMaxHs_meas=" << eMaxHs_meas << "\n";
      std::cerr << "None of them should be a NaN\n";
      throw TerminalException{1};
    }
    if (Swh_used > eMaxHs_meas || Swh_used < eMinHs_meas)
      Swh_used = LargeValue;
    fEnt.Swh_used = Swh_used;
    //
    // Now setting up the footprint values
    //
    double eFootprint = ListFootprintKM[eSatellite - 1];
    fEnt.Footprint = eFootprint;
    //
    // Now putting the limitation on the hours.
    //
    if (SelectingHours) {
      std::vector<int> eDate = DATE_ConvertMjd2six(eEnt.Time);
      int eHour = eDate[3];
      if (PositionVect(ListAllowedHours, eHour) == -1)
        return false;
    }
    //
    // We passed everything. Now inserting
    //
    eRetList.push_back(fEnt);
    return true;
  };
  int len = eRecSer.LastDayTime - eRecSer.FirstDayTime;
  for (int eDay = eRecSer.FirstDayTime; eDay <= eRecSer.LastDayTime; eDay++) {
    double eDayDoubl = static_cast<double>(eDay);
    int pos = eDay - eRecSer.FirstDayTime;
    std::vector<int> eDate = DATE_ConvertMjd2six(eDayDoubl);
    std::string strPres = DATE_ConvertMjd2mystringPres(eDayDoubl);
    int year = eDate[0];
    int month = eDate[1];
    int day = eDate[2];
    std::cerr << "pos=" << pos << " / " << len << "    eDay=" << eDay
              << " date=" << strPres << "\n";
    for (int iType = 0; iType < nbType; iType++) {
      std::string eTypeData = ListTypeData[iType];
      std::string eDirData = ListDirData[iType];
      auto ReadAltiFile = [&]() -> std::vector<SingleEntryMeasurement> {
        if (eTypeData == "IFREMER") {
          std::string eFileAlti =
              eDirData + StringNumber(year, 4) + "/" + StringNumber(month, 2) +
              "/wm_" + StringNumber(year, 4) + StringNumber(month, 2) +
              StringNumber(day, 2) + ".nc";
          if (!IsExistingFile(eFileAlti)) {
            std::cerr << "Error in IFREMER case\n";
            std::cerr << "Missing file " + eFileAlti + "\n";
            std::cerr << "eDirData=" << eDirData << "\n";
            std::cerr << "Please use perl script DownloadAltimeterIfremer\n";
            std::cerr << "for downloading the data\n";
            throw TerminalException{1};
          }
          return READ_ALTI_FILE_IFREMER(eFileAlti);
        }
        if (eTypeData == "EUMETCAST") {
          return READ_ALTI_FILE_EUMETCAST(year, month, day, eDirData);
        }
        if (eTypeData == "RADAR_SPLIT") {
          return READ_RADAR_CSV_FILE(year, month, day, eDirData);
        }
        std::cerr << "Did not find a matching entry for eTypeData\n";
        std::cerr << "Allowed values are IFREMER, EUMETCAST, RADAR_SPLIT\n";
        throw TerminalException{1};
      };
      std::vector<SingleEntryMeasurement> ListEnt = ReadAltiFile();
      std::cerr << "After the data READS\n";
      int nbEnt = ListEnt.size();
      std::cerr << "GEOSELECTION = " << GEOSELECTION << "\n";
      int nbIns = 0;
      for (auto &eEnt : ListEnt) {
        // std::cerr << "Before FuncInsert\n";
        bool res = FuncInsertEntry(eEnt);
        // std::cerr << "After FuncInsert\n";
        if (res)
          nbIns++;
      }
      std::cerr << "  eTypeData=" << eTypeData << " nbEnt=" << nbEnt
                << " nbIns=" << nbIns << "\n";
    }
  }
  std::cerr << "Before returning eRetList\n";
  std::vector<double> ListMJD;
  for (auto &eEnt : eRetList) {
    double eMJD = eEnt.Time;
    ListMJD.push_back(eMJD);
  }
  double minMJD = VectorMin(ListMJD);
  double maxMJD = VectorMax(ListMJD);
  std::string dateMin = DATE_ConvertMjd2mystringPres(minMJD);
  std::string dateMax = DATE_ConvertMjd2mystringPres(maxMJD);
  std::cerr << "DATE min = " << dateMin << " max=" << dateMax
            << " |ListMJD|=" << ListMJD.size()
            << " |eRetList|=" << eRetList.size() << "\n";
  return eRetList;
}

struct PairListWindWave {
  int eSat;
  std::vector<double> ListTimeWind;
  std::vector<double> ListTimeWave;
  std::vector<PairMM> ListPairWind;
  std::vector<PairMM> ListPairWave;
  std::vector<PairLL> ListPairLLwind;
  std::vector<PairLL> ListPairLLwave;
};

struct SatelliteListTrack {
  int eSat;
  std::vector<std::vector<SingleEntryMeasurement>> ListListEntAltimeter;
  double avgDistKM;
};

std::vector<int> GetEliminationWrongModelInterpolation(
    std::vector<SingleEntryMeasurement> const &eVectEnt,
    FullNamelist const &eFull) {
  SingleBlock const& BlPROCESS = eFull.get_block("PROCESS");
  int siz = eVectEnt.size();
  std::vector<int> ListStatus(siz, 1);
  for (int i = 0; i < siz; i++) {
    SingleEntryMeasurement eEnt = eVectEnt[i];
    int eStatus = 1;
    if (BlPROCESS.get_bool("DO_WNDMAG")) {
      for (auto &eVal : eEnt.WindSpeed_model)
        if (eVal < -700)
          eStatus = 0;
    }
    if (BlPROCESS.get_bool("DO_HS")) {
      for (auto &eVal : eEnt.Swh_model)
        if (eVal < -700)
          eStatus = 0;
    }
    ListStatus[i] = eStatus;
  }
  return ListStatus;
}

std::vector<int>
GetListStatusTrackLength(std::vector<SingleEntryMeasurement> const &eVectEnt,
                         FullNamelist const &eFull) {
  int siz = eVectEnt.size();
  std::vector<int> ListStatus(siz, 1);
  auto eBlSELECT = eFull.get_block("SELECT");
  double MaxDistTrackPointKM =
    eBlSELECT.get_double("MaxDistTrackPointKM");
  double MaxDeltaTimeTrackPointDay =
    eBlSELECT.get_double("MaxDeltaTimeTrackPointDay");
  int MinimalTrackSize = eBlSELECT.get_int("MinimalTrackSize");
  std::set<int> SatelliteId;
  for (auto &eEnt : eVectEnt)
    SatelliteId.insert(eEnt.Satellite);
  for (auto &eSat : SatelliteId) {
    std::vector<SingleEntryMeasurement> ListEnt;
    std::vector<int> ListIL;
    int nbL = eVectEnt.size();
    for (int iL = 0; iL < nbL; iL++) {
      SingleEntryMeasurement eEnt = eVectEnt[iL];
      if (eEnt.Satellite == eSat) {
        ListEnt.push_back(eEnt);
        ListIL.push_back(iL);
      }
    }
    int nbEnt = ListEnt.size();
    for (int iEnt = 0; iEnt < nbEnt - 1; iEnt++) {
      double eDiff = ListEnt[iEnt + 1].Time - ListEnt[iEnt].Time;
      if (eDiff < 0) {
        std::cerr << "We have a decrease in time when it should increase\n";
        throw TerminalException{1};
      }
    }
    std::vector<int> ListSep;
    for (int iEnt = 0; iEnt < nbEnt - 1; iEnt++) {
      double eLon1 = ListEnt[iEnt].Lon;
      double eLat1 = ListEnt[iEnt].Lat;
      double eLon2 = ListEnt[iEnt + 1].Lon;
      double eLat2 = ListEnt[iEnt + 1].Lat;
      double eDistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
      double DeltaTimeDay = fabs(ListEnt[iEnt + 1].Time - ListEnt[iEnt].Time);
      if (eDistKM < MaxDistTrackPointKM &&
          DeltaTimeDay < MaxDeltaTimeTrackPointDay)
        ListSep.push_back(iEnt);
    }
    int nbScene = ListSep.size() + 1;
    for (int iScene = 0; iScene < nbScene; iScene++) {
      int iBegin, iEnd;
      if (nbScene == 0) {
        iBegin = 0;
        iEnd = nbEnt - 1;
      } else {
        if (iScene == 0) {
          iBegin = 0;
          iEnd = ListSep[0];
        } else {
          if (iScene == nbScene - 1) {
            iBegin = ListSep[nbScene - 2] + 1;
            iEnd = nbEnt - 1;
          } else {
            iBegin = ListSep[iScene - 1] + 1;
            iEnd = ListSep[iScene];
          }
        }
      }
      int len = 1 + iEnd - iBegin;
      if (len < MinimalTrackSize) {
        for (int i = iBegin; i <= iEnd; i++) {
          int iL = ListIL[i];
          ListStatus[iL] = 0;
        }
      }
    }
  }
  return ListStatus;
}

struct SmoothArr {
  int LenTotal;
  std::vector<int> ListShift;
  std::vector<double> ListWeight;
};

SmoothArr GetSmoothingArray(double const &avgDistKM_model,
                            double const &avgDistKM_track) {
  double TheSizeReal = avgDistKM_model / avgDistKM_track;
  double w = (TheSizeReal - 1) / 2;
  //  std::cerr << "TheSizeReal=" << TheSizeReal << " w=" << w << "\n";
  int wlow = static_cast<int>(floor(w));
  int LenTotal = 1 + 2 * wlow + 2;
  //  std::cerr << "wlow=" << wlow << " LenTotal=" << LenTotal << "\n";
  std::vector<int> ListShift;
  std::vector<double> ListWeight;
  double eOne = 1.0;
  ListShift.push_back(0);
  ListWeight.push_back(eOne);
  for (int widx = 1; widx <= wlow; widx++) {
    ListShift.push_back(widx);
    ListWeight.push_back(eOne);
    ListShift.push_back(-widx);
    ListWeight.push_back(eOne);
  }
  double wres = (TheSizeReal - eOne - static_cast<double>(2 * wlow)) / 2;
  ListShift.push_back(wlow + 1);
  ListWeight.push_back(wres);
  ListShift.push_back(-1 - wlow);
  ListWeight.push_back(wres);
  //  for (int i=0; i<LenTotal; i++)
  //    std::cerr << "i=" << i << " shift=" << ListShift[i] << " weight=" <<
  //    ListWeight[i] << "\n";
  return {LenTotal, ListShift, ListWeight};
}

std::vector<SingleEntryMeasurement>
SpatialAveragingTrack(std::vector<SingleEntryMeasurement> const &eVectEnt,
                      SmoothArr const &eSmoothArr) {
  int len = eVectEnt.size();
  std::vector<SingleEntryMeasurement> eVectRet(len);
  int LenTotal = eSmoothArr.LenTotal;
  for (int i = 0; i < len; i++) {
    SingleEntryMeasurement eRefEnt = eVectEnt[i];
    double Time = eRefEnt.Time;
    double Lon = eRefEnt.Lon;
    double Lat = eRefEnt.Lat;
    std::vector<double> WindSpeed_model = eRefEnt.WindSpeed_model;
    std::vector<double> Swh_model = eRefEnt.Swh_model;
    //
    double SumSwh = 0;
    double SumSwh_cor = 0;
    double SumSwh_used = 0;
    double SumWindSpeed = 0;
    double SumWindSpeed_cor = 0;
    double SumWindSpeed_used = 0;
    double SumWeight = 0;
    for (int iShift = 0; iShift < LenTotal; iShift++) {
      int iNew = i + eSmoothArr.ListShift[iShift];
      double eWeight = eSmoothArr.ListWeight[iShift];
      if (iNew >= 0 && iNew < len) {
        auto eRef = eVectEnt[iNew];
        SumWeight += eWeight;
        SumSwh += eRef.Swh;
        SumSwh_cor += eRef.Swh_cor;
        SumSwh_used += eRef.Swh_used;
        SumWindSpeed += eRef.WindSpeed;
        SumWindSpeed_cor += eRef.WindSpeed_cor;
        SumWindSpeed_used += eRef.WindSpeed_used;
      }
    }
    SingleEntryMeasurement eNewEnt;
    eNewEnt.Time = Time;
    eNewEnt.Lon = Lon;
    eNewEnt.Lat = Lat;
    eNewEnt.WindSpeed_model = WindSpeed_model;
    eNewEnt.Swh_model = Swh_model;
    eNewEnt.Satellite = eRefEnt.Satellite;
    //
    eNewEnt.Swh = SumSwh / SumWeight;
    eNewEnt.Swh_cor = SumSwh_cor / SumWeight;
    eNewEnt.Swh_used = SumSwh_used / SumWeight;
    eNewEnt.WindSpeed = SumWindSpeed / SumWeight;
    eNewEnt.WindSpeed_cor = SumWindSpeed_cor / SumWeight;
    eNewEnt.WindSpeed_used = SumWindSpeed_used / SumWeight;
    eVectRet[i] = eNewEnt;
  }
  return eVectRet;
}

std::vector<SatelliteListTrack>
GetListTrackAltimeter(std::vector<SingleEntryMeasurement> const &eVectEnt,
                      double const &avgDistKM_model,
                      FullNamelist const &eFull) {
  std::cerr << "GetListTrackAltimeter |eVectEnt|=" << eVectEnt.size() << "\n";
  std::vector<SatelliteListTrack> RetList;
  SingleBlock const& eBlSELECT = eFull.get_block("SELECT");
  double MaxDistTrackPointKM =
    eBlSELECT.get_double("MaxDistTrackPointKM");
  double MaxDeltaTimeTrackPointDay =
    eBlSELECT.get_double("MaxDeltaTimeTrackPointDay");
  int MethodTrackSmoothing = eBlSELECT.get_int("MethodTrackSmoothing");
  std::vector<std::string> ListFootprintKM_str =
    eBlSELECT.get_list_string("ListRadiusFootprintKM");
  std::vector<double> ListFootprintKM = ConvertListStringValueToVector_SAT(
      ListFootprintKM_str, 40); // 40km footprint
  //
  std::set<int> SatelliteId;
  for (auto &eEnt : eVectEnt)
    SatelliteId.insert(eEnt.Satellite);
  for (auto &eSat : SatelliteId) {
    std::vector<SingleEntryMeasurement> ListEnt;
    for (auto &eEnt : eVectEnt)
      if (eEnt.Satellite == eSat)
        ListEnt.push_back(eEnt);
    int nbEnt = ListEnt.size();
    for (int iEnt = 0; iEnt < nbEnt - 1; iEnt++) {
      double eDiff = ListEnt[iEnt + 1].Time - ListEnt[iEnt].Time;
      if (eDiff < 0) {
        std::cerr << "We have a decrease in time when it should increase\n";
        throw TerminalException{1};
      }
    }
    std::vector<int> ListSep;
    double SumDistKM = 0;
    int nbPair = 0;
    for (int iEnt = 0; iEnt < nbEnt - 1; iEnt++) {
      double eLon1 = ListEnt[iEnt].Lon;
      double eLat1 = ListEnt[iEnt].Lat;
      double eLon2 = ListEnt[iEnt + 1].Lon;
      double eLat2 = ListEnt[iEnt + 1].Lat;
      double eDistKM = GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
      double DeltaTimeDay = fabs(ListEnt[iEnt + 1].Time - ListEnt[iEnt].Time);
      if (eDistKM < MaxDistTrackPointKM &&
          DeltaTimeDay < MaxDeltaTimeTrackPointDay) {
        SumDistKM += eDistKM;
        nbPair++;
      } else {
        ListSep.push_back(iEnt);
      }
    }
    double avgDistKM_track = SumDistKM / static_cast<double>(nbPair);
    int nbScene = ListSep.size() + 1;
    std::cerr << "nbScene=" << nbScene << "\n";
    std::vector<std::vector<SingleEntryMeasurement>> ListListEntAltimeter;
    for (int iScene = 0; iScene < nbScene; iScene++) {
      int iBegin, iEnd;
      if (nbScene == 1) {
        iBegin = 0;
        iEnd = nbEnt - 1;
      } else {
        if (iScene == 0) {
          iBegin = 0;
          iEnd = ListSep[0];
        } else {
          if (iScene == nbScene - 1) {
            iBegin = ListSep[nbScene - 2] + 1;
            iEnd = nbEnt - 1;
          } else {
            iBegin = ListSep[iScene - 1] + 1;
            iEnd = ListSep[iScene];
          }
        }
      }
      int len = 1 + iEnd - iBegin;
      std::vector<SingleEntryMeasurement> ListEntAltimeter(len);
      for (int i = iBegin; i <= iEnd; i++)
        ListEntAltimeter[i - iBegin] = ListEnt[i];
      double avgDistKM_satellite = ListFootprintKM[eSat - 1];
      //      std::cerr << "avgDistKM_model=" << avgDistKM_model << "
      //      avgDistKM_track=" << avgDistKM_track << " avgDistKM_satellite=" <<
      //      avgDistKM_satellite << " MethodTrackSmoothing=" <<
      //      MethodTrackSmoothing << "\n";
      auto GetTrackValues = [&]() -> std::vector<SingleEntryMeasurement> {
        if (MethodTrackSmoothing == 0)
          return ListEntAltimeter;
        if (MethodTrackSmoothing == 1) {
          SmoothArr eSmoothArr =
              GetSmoothingArray(avgDistKM_model, avgDistKM_track);
          return SpatialAveragingTrack(ListEntAltimeter, eSmoothArr);
        }
        if (MethodTrackSmoothing == 2) {
          SmoothArr eSmoothArr =
              GetSmoothingArray(avgDistKM_satellite, avgDistKM_track);
          return SpatialAveragingTrack(ListEntAltimeter, eSmoothArr);
        }
        std::cerr << "MethodTrackSmoothing=" << MethodTrackSmoothing << "\n";
        throw TerminalException{1};
      };
      ListListEntAltimeter.push_back(GetTrackValues());
    }
    RetList.push_back({eSat, ListListEntAltimeter, avgDistKM_track});
    std::cerr << "eSat=" << eSat
              << " |ListListEntAltimeter|=" << ListListEntAltimeter.size()
              << "\n";
  }
  std::cerr << "GetListTrackAltimeter |RetList|=" << RetList.size() << "\n";
  return RetList;
}

std::set<int>
GetListSatelliteId_set(std::vector<SingleEntryMeasurement> const &eVectEnt,
                       FullNamelist const &eFull) {
  std::set<int> PreSatelliteId;
  for (auto &eEnt : eVectEnt)
    PreSatelliteId.insert(eEnt.Satellite);
  CheckSatelliteList(eFull);
  SingleBlock const& BlockSELECT = eFull.get_block("SELECT");
  std::vector<std::string> AllSatNames = GetAllNamesOfSatelliteAltimeter();
  std::vector<std::string> AllowedSatNames = BlockSELECT.get_list_string("AllowedSatellites");
  std::set<int> SatelliteId;
  for (auto &eEnt : PreSatelliteId) {
    std::string eSatName = AllSatNames[eEnt - 1];
    bool WeMatch = false;
    for (auto &fSatName : AllowedSatNames)
      if (eSatName == fSatName)
        WeMatch = true;
    if (WeMatch)
      SatelliteId.insert(eEnt);
  }
  return SatelliteId;
}

void RAW_STATISTICS_ALTIMETER(std::vector<PairListWindWave> const &eSS,
                              SatelliteSerInfo const &eRecSer,
                              PermanentInfoDrawing const &ePerm,
                              int const &iGrid) {
  std::cerr << "Running RAW_STATISTICS_ALTIMETER\n";
  std::string FileName =
      ePerm.eDir + "Statistics_iGrid" + StringNumber(iGrid, 1) + ".txt";
  std::ofstream os(FileName);
  os << "Comparison of model results with altimeter\n";
  os << "Beginning time = " << eRecSer.strBegin << "\n";
  os << "Ending    time = " << eRecSer.strEnd << "\n";
  SingleBlock eBlPROC = ePerm.eFull.get_block("PROC");
  SingleBlock eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  os << "HisPrefix = " << ListHisPrefix[iGrid] << "\n";
  os << "ModelName = " << ListModelName[iGrid] << "\n";
  os << "GridFile = " << ListGridFile[iGrid] << "\n";
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  os << "Do wind statistic = " << DO_WNDMAG << "\n";
  os << "Do wave statistic = " << DO_HS << "\n";
  if (DO_WNDMAG) {
    std::vector<PairMM> TotalListPairWind;
    for (auto &eRec : eSS)
      TotalListPairWind.insert(TotalListPairWind.end(),
                               eRec.ListPairWind.begin(),
                               eRec.ListPairWind.end());
    std::string eNameWind = "All satellites wind statistics";
    T_stat eStatWind = ComputeStatistics_Pair(TotalListPairWind);
    Print_Down_Statistics(os, eNameWind, eStatWind);
    for (auto &eRec : eSS) {
      std::string SatName = GetNameOfSatelliteAltimeter(eRec.eSat);
      std::string eNameWindB = SatName + " wind statistics";
      T_stat eStatWindB = ComputeStatistics_Pair(eRec.ListPairWind);
      Print_Down_Statistics(os, eNameWindB, eStatWindB);
    }
  }
  if (DO_HS) {
    std::vector<PairMM> TotalListPairWave;
    for (auto &eRec : eSS)
      TotalListPairWave.insert(TotalListPairWave.end(),
                               eRec.ListPairWave.begin(),
                               eRec.ListPairWave.end());
    std::string eNameWave = "All satellites wave statistics";
    //    std::cerr << "Measurements and model of HWAVE\n";
    //    for (auto & ePair : TotalListPairWave) {
    //      std::cerr << "meas=" << ePair.Meas << " model=" << ePair.Model <<
    //      "\n";
    //    }
    T_stat eStatWave = ComputeStatistics_Pair(TotalListPairWave);
    Print_Down_Statistics(os, eNameWave, eStatWave);
    for (auto &eRec : eSS) {
      std::string SatName = GetNameOfSatelliteAltimeter(eRec.eSat);
      std::string eNameWaveB = SatName + " wave statistics";
      T_stat eStatWaveB = ComputeStatistics_Pair(eRec.ListPairWave);
      Print_Down_Statistics(os, eNameWaveB, eStatWaveB);
    }
  }
}

void RAW_STATISTICS_BREAKDOWN_ALTIMETER(
    std::vector<PairListWindWave> const &eSS, SatelliteSerInfo const &eRecSer,
    PermanentInfoDrawing const &ePerm, int const &iGrid,
    std::string const &OptPart) {
  std::cerr << "Running RAW_STATISTICS_BREAKDOWN_ALTIMETER OptPart=" << OptPart
            << "\n";
  std::string FileName = ePerm.eDir + "Statistics_" + OptPart + "_iGrid" +
                         StringNumber(iGrid, 1) + ".txt";
  std::ofstream os(FileName);
  os << "Comparison of model results with altimeter\n";
  os << "Beginning time = " << eRecSer.strBegin << "\n";
  os << "Ending    time = " << eRecSer.strEnd << "\n";
  SingleBlock eBlPROC = ePerm.eFull.get_block("PROC");
  SingleBlock eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  SingleBlock eBlSELECT = ePerm.eFull.get_block("SELECT");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  os << "HisPrefix = " << ListHisPrefix[iGrid] << "\n";
  os << "ModelName = " << ListModelName[iGrid] << "\n";
  os << "GridFile = " << ListGridFile[iGrid] << "\n";
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  os << "Do wind statistic = " << DO_WNDMAG << "\n";
  os << "Do wave statistic = " << DO_HS << "\n";
  std::vector<double> ListLineLonStart =
    eBlSELECT.get_list_double("ListLineLonStart");
  std::vector<double> ListLineLatStart =
    eBlSELECT.get_list_double("ListLineLatStart");
  std::vector<double> ListLineLonEnd =
    eBlSELECT.get_list_double("ListLineLonEnd");
  std::vector<double> ListLineLatEnd =
    eBlSELECT.get_list_double("ListLineLatEnd");
  std::set<size_t> ListLen{ListLineLonStart.size(), ListLineLatStart.size(),
                           ListLineLonEnd.size(), ListLineLatEnd.size()};
  if (ListLen.size() != 1) {
    std::cerr << "All the ListLineLonStart, LatStart, LonEnd, LatEnd should "
                 "have the same length\n";
    throw TerminalException{1};
  }
  int nbLine = ListLineLonStart.size();
  auto FuncPrint = [&](int const &opt) -> void {
    std::vector<double> TotalListTime;
    std::vector<PairMM> TotalListPair;
    std::vector<PairLL> TotalListPairLL;
    for (auto &eRec : eSS) {
      if (opt == 1) {
        TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWind.begin(),
                             eRec.ListTimeWind.end());
        TotalListPair.insert(TotalListPair.end(), eRec.ListPairWind.begin(),
                             eRec.ListPairWind.end());
        TotalListPairLL.insert(TotalListPairLL.end(),
                               eRec.ListPairLLwind.begin(),
                               eRec.ListPairLLwind.end());
      }
      if (opt == 2) {
        TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWave.begin(),
                             eRec.ListTimeWave.end());
        TotalListPair.insert(TotalListPair.end(), eRec.ListPairWave.begin(),
                             eRec.ListPairWave.end());
        TotalListPairLL.insert(TotalListPairLL.end(),
                               eRec.ListPairLLwave.begin(),
                               eRec.ListPairLLwave.end());
      }
    }
    std::set<std::string> SetDisc;
    std::vector<std::string> ListDisc;
    int nbEnt = TotalListTime.size();
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      double eTime = TotalListTime[iEnt];
      std::vector<int> eDate = DATE_ConvertMjd2six(eTime);
      int eYear = eDate[0];
      int eMonth = eDate[1];
      //      int eDay=eDate[2];
      int eHour = eDate[3];
      std::string eDisc;
      if (OptPart == "monthly") {
        eDisc = StringNumber(eYear, 4) + " / " + StringNumber(eMonth, 2);
      }
      if (OptPart == "hourly") {
        eDisc = StringNumber(eHour, 2);
      }
      if (OptPart == "geography") {
        for (int iLine = 0; iLine < nbLine; iLine++) {
          double eLon = TotalListPairLL[iEnt].eLon;
          double eLat = TotalListPairLL[iEnt].eLat;
          double eLonStart = ListLineLonStart[iLine];
          double eLatStart = ListLineLatStart[iLine];
          double eLonEnd = ListLineLonEnd[iLine];
          double eLatEnd = ListLineLatEnd[iLine];
          double scal = (eLon - eLonStart) * (eLatEnd - eLatStart) -
                        (eLat - eLatStart) * (eLonEnd - eLonStart);
          if (scal > 0)
            eDisc += " 1";
          else
            eDisc += " 0";
        }
      }
      SetDisc.insert(eDisc);
      ListDisc.push_back(eDisc);
    }
    //
    os << "Monthly statistics: Mean Error / RMSE / Correlation / ScatterIndex "
          ": nbMeas\n";
    for (auto &eDisc : SetDisc) {
      std::vector<PairMM> SelectListPair;
      int len = TotalListTime.size();
      for (int i = 0; i < len; i++)
        if (ListDisc[i] == eDisc)
          SelectListPair.push_back(TotalListPair[i]);
      T_stat eStat = ComputeStatistics_Pair(SelectListPair);
      os << eDisc << " : " << eStat.MeanError << " / " << eStat.RMSE << " / "
         << eStat.Correlation << " / " << eStat.ScatterIndex << "      :  "
         << eStat.nbMeas << "\n";
    }
  };
  if (DO_WNDMAG) {
    os << "Wind speed statistics\n";
    FuncPrint(1);
  }
  if (DO_HS) {
    os << "Wave height statistics\n";
    FuncPrint(2);
  }
}

void RAW_BIN_RANGING_STATISTICS_BREAKDOWN(
    std::vector<PairListWindWave> const &eSS, SatelliteSerInfo const &eRecSer,
    PermanentInfoDrawing const &ePerm, int const &iGrid) {
  std::cerr << "Running RAW_BIN_RANGING_STATISTICS_BREAKDOWN\n";
  std::string FileName = ePerm.eDir + "Statistics_bin_ranging_iGrid" +
                         StringNumber(iGrid, 1) + ".txt";
  std::ofstream os(FileName);
  os << "Comparison of model results with altimeter\n";
  os << "Beginning time = " << eRecSer.strBegin << "\n";
  os << "Ending    time = " << eRecSer.strEnd << "\n";
  SingleBlock const& eBlPROC = ePerm.eFull.get_block("PROC");
  SingleBlock const& eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  SingleBlock const& eBlSELECT = ePerm.eFull.get_block("SELECT");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  os << "HisPrefix = " << ListHisPrefix[iGrid] << "\n";
  os << "ModelName = " << ListModelName[iGrid] << "\n";
  os << "GridFile = " << ListGridFile[iGrid] << "\n";
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  os << "Do wind statistic = " << DO_WNDMAG << "\n";
  os << "Do wave statistic = " << DO_HS << "\n";
  std::vector<int> ListSat = {-1};
  for (auto &eRec : eSS)
    ListSat.push_back(eRec.eSat);
  auto FuncPrint = [&](int const &eSatIn, int const &opt,
                       std::vector<double> const ListBinValues) -> void {
    std::vector<double> TotalListTime;
    std::vector<PairMM> TotalListPair;
    std::vector<PairLL> TotalListPairLL;
    for (auto &eRec : eSS) {
      if (eSatIn == -1 || eSatIn == eRec.eSat) {
        if (opt == 1) {
          TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWind.begin(),
                               eRec.ListTimeWind.end());
          TotalListPair.insert(TotalListPair.end(), eRec.ListPairWind.begin(),
                               eRec.ListPairWind.end());
          TotalListPairLL.insert(TotalListPairLL.end(),
                                 eRec.ListPairLLwind.begin(),
                                 eRec.ListPairLLwind.end());
        }
        if (opt == 2) {
          TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWave.begin(),
                               eRec.ListTimeWave.end());
          TotalListPair.insert(TotalListPair.end(), eRec.ListPairWave.begin(),
                               eRec.ListPairWave.end());
          TotalListPairLL.insert(TotalListPairLL.end(),
                                 eRec.ListPairLLwave.begin(),
                                 eRec.ListPairLLwave.end());
        }
      }
    }
    int nbBin = ListBinValues.size();
    if (nbBin == 0) {
      os << "ListBinValues is empty. So, we cannot compute the bin sorting\n";
      return;
    }
    std::vector<int> ListMatchMeas(nbBin + 1, 0);
    std::vector<int> ListMatchModel(nbBin + 1, 0);
    auto fInsert = [&](std::vector<int> &eVectMatch,
                       double const &eVal) -> void {
      if (eVal < ListBinValues[0]) {
        eVectMatch[0]++;
        return;
      }
      for (int iBin = 1; iBin < nbBin; iBin++) {
        if (ListBinValues[iBin - 1] <= eVal && eVal < ListBinValues[iBin]) {
          eVectMatch[iBin]++;
          return;
        }
      }
      eVectMatch[nbBin]++;
    };

    for (auto &ePair : TotalListPair) {
      fInsert(ListMatchMeas, ePair.Meas);
      fInsert(ListMatchModel, ePair.Model);
    }
    //
    os << "Value by meas / model |TotalListPair|=" << TotalListPair.size()
       << "\n";
    os << "Int = (0, " << ListBinValues[0] << ") nb=" << ListMatchMeas[0]
       << " / " << ListMatchModel[0] << "\n";
    for (int iBin = 1; iBin < nbBin; iBin++) {
      os << "Int = (" << ListBinValues[iBin - 1] << ", " << ListBinValues[iBin]
         << ") nb=" << ListMatchMeas[iBin] << " / " << ListMatchModel[iBin]
         << "\n";
    }
    os << "Int = (" << ListBinValues[nbBin - 1]
       << ", inf) nb=" << ListMatchMeas[nbBin] << " / " << ListMatchModel[nbBin]
       << "\n";
  };
  auto GetRelSatName = [&](int const &eSat) -> std::string {
    if (eSat == -1)
      return "all satellites";
    return GetNameOfSatelliteAltimeter(eSat);
  };

  for (auto &eSat : ListSat) {
    os << "----------------- Analysis for : " << GetRelSatName(eSat)
       << " ------------------\n";
    if (DO_WNDMAG) {
      os << "Wind speed statistics\n";
      std::vector<double> ListBinWindValues =
        eBlSELECT.get_list_double("ListBinWindValues");
      FuncPrint(eSat, 1, ListBinWindValues);
    }
    if (DO_HS) {
      os << "Wave height statistics\n";
      std::vector<double> ListBinWaveValues =
        eBlSELECT.get_list_double("ListBinWaveValues");
      FuncPrint(eSat, 2, ListBinWaveValues);
    }
  }
}

void BREAKDOWN_GEOG_POINT(std::vector<PairListWindWave> const &eSS,
                          SatelliteSerInfo const &eRecSer,
                          NCLcaller<GeneralType> &eCall,
                          PermanentInfoDrawing const &ePerm, int const &iGrid) {
  std::cerr << "Running BREAKDOWN_GEOG_POINT\n";
  if (eSS.size() == 0) {
    std::cerr << "eSS array is empty. Probably no valid measurement. Exiting "
                 "the subroutine without output\n";
    return;
  }
  SingleBlock const& eBlPLOT = ePerm.eFull.get_block("PLOT");
  SingleBlock const& eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  SingleBlock const& eBlSELECT = ePerm.eFull.get_block("SELECT");
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  double tolLL = eBlSELECT.get_double("tolLL");
  bool PlotIndividualPoint =
    eBlPROCESS.get_bool("PlotIndividualPoints");
  bool PlotGeographicStat = eBlPROCESS.get_bool("PlotGeographicStat");
  int nbBlock = eBlPLOT.get_int("nbBlock");
  auto FuncPrint = [&](int const &opt, std::string const &varSymb,
                       std::string const &unitStr) -> void {
    std::vector<double> TotalListTime;
    std::vector<PairMM> TotalListPairMM;
    std::vector<PairLL> TotalListPairLL;
    for (auto &eRec : eSS) {
      if (opt == 1) {
        TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWind.begin(),
                             eRec.ListTimeWind.end());
        TotalListPairMM.insert(TotalListPairMM.end(), eRec.ListPairWind.begin(),
                               eRec.ListPairWind.end());
        TotalListPairLL.insert(TotalListPairLL.end(),
                               eRec.ListPairLLwind.begin(),
                               eRec.ListPairLLwind.end());
      }
      if (opt == 2) {
        TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWave.begin(),
                             eRec.ListTimeWave.end());
        TotalListPairMM.insert(TotalListPairMM.end(), eRec.ListPairWave.begin(),
                               eRec.ListPairWave.end());
        TotalListPairLL.insert(TotalListPairLL.end(),
                               eRec.ListPairLLwave.begin(),
                               eRec.ListPairLLwave.end());
      }
    }
    std::vector<PairLL> ListPointLL;
    std::vector<std::vector<double>> ListListTime;
    std::vector<std::vector<PairMM>> ListListMeas;
    auto FuncInsert = [&](PairLL const &ePairLL, PairMM const &ePairMM,
                          double const &eTime) -> void {
      int len = ListPointLL.size();
      for (int u = 0; u < len; u++) {
        double err_LL = fabs(ListPointLL[u].eLon - ePairLL.eLon) +
                        fabs(ListPointLL[u].eLat - ePairLL.eLat);
        if (err_LL < tolLL) {
          ListListTime[u].push_back(eTime);
          ListListMeas[u].push_back(ePairMM);
          return;
        }
      }
      ListPointLL.push_back(ePairLL);
      ListListTime.push_back({eTime});
      ListListMeas.push_back({ePairMM});
    };
    std::cerr << "FuncPrint, step 1\n";

    int nbEnt = TotalListPairLL.size();
    for (int iEnt = 0; iEnt < nbEnt; iEnt++)
      FuncInsert(TotalListPairLL[iEnt], TotalListPairMM[iEnt],
                 TotalListTime[iEnt]);
    int nbPointLL = ListPointLL.size();
    std::cerr << "FuncPrint, step 2 nbPointLL=" << nbPointLL << "\n";
    //
    // Now finding the squares
    //
    auto FuncInsertDouble = [&](std::set<double> &TheList,
                                double const &eVal) -> void {
      for (auto &fVal : TheList) {
        double eDist = fabs(fVal - eVal);
        if (eDist < tolLL)
          return;
      }
      TheList.insert(eVal);
    };
    auto GetPositionDouble = [&](std::vector<double> &TheList,
                                 double const &eVal) -> int {
      int eLen = TheList.size();
      for (int u = 0; u < eLen; u++) {
        double eDist = fabs(TheList[u] - eVal);
        if (eDist < tolLL)
          return u;
      }
      std::cerr << "We should not reach that stage\n";
      return -1;
    };
    std::set<double> SetLon;
    std::set<double> SetLat;
    for (int iPoint = 0; iPoint < nbPointLL; iPoint++) {
      double eLon = ListPointLL[iPoint].eLon;
      double eLat = ListPointLL[iPoint].eLat;
      FuncInsertDouble(SetLon, eLon);
      FuncInsertDouble(SetLat, eLat);
    }
    std::vector<double> ListLon;
    std::vector<double> ListLat;
    for (auto &eVal : SetLon)
      ListLon.push_back(eVal);
    for (auto &eVal : SetLat)
      ListLat.push_back(eVal);

    std::vector<int> ListIdxLon(nbPointLL);
    std::vector<int> ListIdxLat(nbPointLL);
    for (int iPoint = 0; iPoint < nbPointLL; iPoint++) {
      double eLon = ListPointLL[iPoint].eLon;
      double eLat = ListPointLL[iPoint].eLat;
      int idxLon = GetPositionDouble(ListLon, eLon);
      int idxLat = GetPositionDouble(ListLat, eLat);
      ListIdxLon[iPoint] = idxLon;
      ListIdxLat[iPoint] = idxLat;
    }
    std::cerr << "FuncPrint, step 3\n";
    int nbLon = ListLon.size();
    int nbLat = ListLat.size();
    MyMatrix<double> LON(nbLon, nbLat);
    MyMatrix<double> LAT(nbLon, nbLat);
    for (int iLon = 0; iLon < nbLon; iLon++)
      for (int iLat = 0; iLat < nbLat; iLat++) {
        LON(iLon, iLat) = ListLon[iLon];
        LAT(iLon, iLat) = ListLat[iLat];
      }
    std::cerr << "FuncPrint, step 4 nbLon=" << nbLon << " nbLat=" << nbLat
              << "\n";
    std::cerr << "LON (min/max)=" << LON.minCoeff() << " / " << LON.maxCoeff()
              << "\n";
    std::cerr << "LAT (min/max)=" << LAT.minCoeff() << " / " << LAT.maxCoeff()
              << "\n";
    MyMatrix<uint8_t> MSK = ZeroMatrix<uint8_t>(nbLon, nbLat);
    for (int iPoint = 0; iPoint < nbPointLL; iPoint++) {
      int idxLon = ListIdxLon[iPoint];
      int idxLat = ListIdxLat[iPoint];
      MSK(idxLon, idxLat) = 1;
    }
    std::cerr << "MSK=\n";
    for (int iLon = 0; iLon < nbLon; iLon++) {
      for (int iLat = 0; iLat < nbLat; iLat++)
        std::cerr << (int)MSK(iLon, iLat);
      std::cerr << "\n";
    }
    MyMatrix<double> ANG = CreateAngleMatrix(LON, LAT);
    if (PlotGeographicStat) {
      GridArray GrdArr;
      GrdArr.ModelName = "Statistic";
      GrdArr.IsFE = 0;
      GrdArr.IsSpherical = true;
      GrdArr.GrdArrRho.LON = LON;
      GrdArr.GrdArrRho.LAT = LAT;
      GrdArr.GrdArrRho.MSK = MSK;
      GrdArr.GrdArrRho.ANG = ANG;
      GrdArr.ARVD.IsAssigned = false;
      GrdArr.ARVD.Zcoordinate = false;
      std::cerr << "FuncPrint, step 6\n";
      //
      int nbTypeStat = 4;
      std::vector<std::string> ListTypeStat{"Correlation", "MeanError", "RMSE",
                                            "nbMeas"};
      std::vector<std::vector<double>> ListValStat(nbTypeStat);
      for (int iPointLL = 0; iPointLL < nbPointLL; iPointLL++) {
        T_stat eStat = ComputeStatistics_Pair(ListListMeas[iPointLL]);
        ListValStat[0].push_back(eStat.Correlation);
        ListValStat[1].push_back(eStat.MeanError);
        ListValStat[2].push_back(eStat.RMSE);
        ListValStat[3].push_back(static_cast<double>(eStat.nbMeas));
      }
      std::cerr << "FuncPrint, step 7\n";
      auto GetRelUnit = [&](std::string const &typeStat) -> std::string {
        if (typeStat == "Correlation" || typeStat == "nbMeas")
          return "1";
        if (typeStat == "RMSE" || typeStat == "MeanError")
          return unitStr;
        throw TerminalException{1};
      };
      for (int iTypeStat = 0; iTypeStat < nbTypeStat; iTypeStat++) {
        std::string eTypeStat = ListTypeStat[iTypeStat];
        std::cerr << "iTypeStat=" << iTypeStat << " eTypeStat=" << eTypeStat
                  << "\n";
        //
        MyMatrix<double> VAL = ZeroMatrix<double>(nbLon, nbLat);
        for (int iPoint = 0; iPoint < nbPointLL; iPoint++) {
          int idxLon = ListIdxLon[iPoint];
          int idxLat = ListIdxLat[iPoint];
          VAL(idxLon, idxLat) = ListValStat[iTypeStat][iPoint];
        }
        std::string FileName = ePerm.eDir + eTypeStat + "_iGrid" +
                               StringNumber(iGrid, 1) + "_" + varSymb;
        //
        PairMinMax ePair = ComputeMinMax(GrdArr, VAL);
        RecVar eRecVar;
        eRecVar.RecS.strAll = eTypeStat;
        eRecVar.RecS.VarName1 = eTypeStat;
        eRecVar.RecS.VarName2 = eTypeStat;
        eRecVar.RecS.minval = ePair.TheMin;
        eRecVar.RecS.maxval = ePair.TheMax;
        eRecVar.RecS.Unit = GetRelUnit(eTypeStat);
        eRecVar.F = VAL;
        std::cerr << "Unit=" << eRecVar.RecS.Unit << "\n";
        //
        double deltaLL = 0.1;
        QuadArray eQuad = GetQuadArray(LON, LAT, deltaLL);
        std::cerr << "eQuad=" << eQuad << "\n";
        // BasicArrayDraw(eQuad);
        DrawArr eDrawArr = ePerm.eDrawArr;
        eDrawArr.eQuadFrame = eQuad;
        //      eDrawArr.DoColorBar=true;
        //      eDrawArr.nbLevelSpa=20;
        eDrawArr.TitleStr = eTypeStat + " for " + varSymb;
        eDrawArr.DoTitle = true;
        eDrawArr.VarNameUF = eTypeStat + "_for_" + varSymb;
        //
        PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
      }
    }
    if (PlotIndividualPoint) {
      for (int iPoint = 0; iPoint < nbPointLL; iPoint++) {
        double eLon = ListPointLL[iPoint].eLon;
        double eLat = ListPointLL[iPoint].eLat;
        int idxLon = ListIdxLon[iPoint];
        int idxLat = ListIdxLat[iPoint];
        int nbTime = ListListMeas[iPoint].size();
        std::vector<int> ListPos = DivideListPosition(nbTime, nbBlock);
        //
        for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
          std::cerr << "iPoint=" << iPoint << " idxLon=" << idxLon
                    << " idxLat=" << idxLat << " eLon=" << eLon
                    << " eLat=" << eLat << " iBlock=" << iBlock << "\n";
          DrawLinesArr eDrawArr;
          eDrawArr.DoTitle = true;
          eDrawArr.TitleStr = "Time series at lon=" + DoubleToString(eLon) +
                              " lat=" + DoubleToString(eLat);
          eDrawArr.IsTimeSeries = true;
          eDrawArr.PairComparison = true;
          eDrawArr.DoExplicitLabel = false;
          eDrawArr.VarName = "Point_" + IntToString(idxLon) + "_" +
                             IntToString(idxLat) + "_" + IntToString(iBlock);
          eDrawArr.ListName_plot = {"model", "meas"};
          eDrawArr.YAxisString = varSymb + "(" + unitStr + ")";
          //
          int pos1 = ListPos[iBlock];
          int pos2 = ListPos[iBlock + 1];
          int len = pos2 - pos1;
          MyVector<double> ListTime(len);
          MyVector<double> ListMeas(len);
          MyVector<double> ListModel(len);
          for (int pos = pos1; pos < pos2; pos++) {
            int posA = pos - pos1;
            ListTime(posA) = ListListTime[iPoint][pos];
            ListMeas(posA) = ListListMeas[iPoint][pos].Meas;
            ListModel(posA) = ListListMeas[iPoint][pos].Model;
          }
          eDrawArr.ListX = ListTime;
          eDrawArr.ListListVect = {ListModel, ListMeas};
          double TheMax = std::max(ListModel.maxCoeff(), ListMeas.maxCoeff());
          double TheMin = std::max(ListModel.minCoeff(), ListMeas.minCoeff());
          eDrawArr.TheMax = TheMax;
          eDrawArr.TheMin = TheMin;
          std::string FileName =
              ePerm.eDir + "TimeSeries_" + IntToString(idxLon) + "_" +
              IntToString(idxLat) + "_" + IntToString(iBlock) + "_" + varSymb;
          LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
        }
      }
    }
  };
  if (DO_WNDMAG) {
    std::cerr << "Wind speed statistics\n";
    FuncPrint(1, "wind", "m/s");
  }
  if (DO_HS) {
    std::cerr << "Wave height statistics\n";
    FuncPrint(2, "hs", "m");
  }
}

void PRINT_RAW_DATA(std::vector<std::vector<PairListWindWave>> const &ListSS,
                    std::vector<SatelliteListTrack> const &LTrack,
                    SatelliteSerInfo const &eRecSer,
                    NCLcaller<GeneralType> &eCall,
                    PermanentInfoDrawing const &ePerm) {
  std::cerr << "Running PRINT_RAW_DATA (enter)\n";
  SingleBlock eBlPROC = ePerm.eFull.get_block("PROC");
  SingleBlock eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  SingleBlock eBlSELECT = ePerm.eFull.get_block("SELECT");
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  auto FuncPrint = [&](int const &opt, std::string const &OutFileTxt) -> void {
    std::vector<MyVector<double>> ListVectorModel;
    MyVector<double> eVectorMeas;
    bool IsFirst = true;
    for (auto &eSS : ListSS) {
      std::vector<double> TotalListTime;
      std::vector<PairMM> TotalListPairMM;
      std::vector<PairLL> TotalListPairLL;
      for (auto &eRec : eSS) {
        if (opt == 1) {
          TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWind.begin(),
                               eRec.ListTimeWind.end());
          TotalListPairMM.insert(TotalListPairMM.end(),
                                 eRec.ListPairWind.begin(),
                                 eRec.ListPairWind.end());
          TotalListPairLL.insert(TotalListPairLL.end(),
                                 eRec.ListPairLLwind.begin(),
                                 eRec.ListPairLLwind.end());
        }
        if (opt == 2) {
          TotalListTime.insert(TotalListTime.end(), eRec.ListTimeWave.begin(),
                               eRec.ListTimeWave.end());
          TotalListPairMM.insert(TotalListPairMM.end(),
                                 eRec.ListPairWave.begin(),
                                 eRec.ListPairWave.end());
          TotalListPairLL.insert(TotalListPairLL.end(),
                                 eRec.ListPairLLwave.begin(),
                                 eRec.ListPairLLwave.end());
        }
      }
      int len = TotalListPairMM.size();
      MyVector<double> eVectorModel(len);
      for (int i = 0; i < len; i++)
        eVectorModel(i) = TotalListPairMM[i].Model;
      ListVectorModel.push_back(eVectorModel);
      if (IsFirst) {
        MyVector<double> eVect(len);
        for (int i = 0; i < len; i++)
          eVect(i) = TotalListPairMM[i].Meas;
        eVectorMeas = eVect;
        IsFirst = false;
      }
    }
    int nbModel = ListSS.size();
    std::string eFile = ePerm.eDir + OutFileTxt;
    std::ofstream os(eFile);
    os << "Writing measurement\n";
    int len = eVectorMeas.size();
    for (int i = 0; i < len; i++) {
      os << "meas=" << eVectorMeas(i) << " model=";
      for (int iModel = 0; iModel < nbModel; iModel++) {
        if (iModel > 0)
          os << " ";
        os << ListVectorModel[iModel](i);
      }
      os << "\n";
    }
  };
  if (DO_WNDMAG) {
    std::cerr << "Wind speed raw data printing\n";
    FuncPrint(1, "RawDataWind.txt");
  }
  if (DO_HS) {
    std::cerr << "Wave height raw data printing\n";
    FuncPrint(2, "RawDataWave.txt");
  }
}

void WRITE_MEASUREMENT_GEOG_POSITION(std::vector<PairListWindWave> const &eSS,
                                     SatelliteSerInfo const &eRecSer,
                                     PermanentInfoDrawing const &ePerm,
                                     int const &iGrid) {
  std::cerr << "Beginning of WRITE_MEASUREMENT_GEOG_POSITION\n";
  SingleBlock eBlPROCESS = ePerm.eFull.get_block("PROCESS");
  SingleBlock eBlSELECT = ePerm.eFull.get_block("SELECT");
  bool DO_WNDMAG = eBlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = eBlPROCESS.get_bool("DO_HS");
  double tolLL = eBlSELECT.get_double("tolLL");
  std::cerr << "tolLL=" << tolLL << "\n";
  auto FuncPrint = [&](int const &opt) -> void {
    std::set<PairLL> TotalListPairLL;
    for (auto &eRec : eSS) {
      if (opt == 1)
        TotalListPairLL.insert(eRec.ListPairLLwind.begin(),
                               eRec.ListPairLLwind.end());
      if (opt == 2)
        TotalListPairLL.insert(eRec.ListPairLLwave.begin(),
                               eRec.ListPairLLwave.end());
    }
    std::cerr << "|TotalListPairLL|=" << TotalListPairLL.size() << "\n";
    std::vector<PairLL> TotalListPairLL_b;
    auto FuncInsert = [&](PairLL const &eLL) -> void {
      for (auto &fLL : TotalListPairLL_b) {
        double err_LL = fabs(fLL.eLon - eLL.eLon) + fabs(fLL.eLat - eLL.eLat);
        if (err_LL < tolLL)
          return;
      }
      TotalListPairLL_b.push_back(eLL);
    };
    for (auto &eLL : TotalListPairLL)
      FuncInsert(eLL);
    int nbPoint = TotalListPairLL_b.size();
    std::cerr << "|TotalListPairLL_b|=" << TotalListPairLL_b.size() << "\n";
    //

    std::cerr << "Writing coordinates into files opt=" << opt << "\n";
    std::string eFileNC = ePerm.eDir + "Geog_position_iGrid" +
                          StringNumber(iGrid, 1) + "_opt" +
                          StringNumber(opt, 1) + ".nc";
    std::cerr << "Coordinates : eFileNC=" << eFileNC << "\n";
    if (!FILE_IsFileMakeable(eFileNC)) {
      std::cerr << "Request to create file eFileNC=" << eFileNC << "\n";
      std::cerr << "but the directory does not exist\n";
      throw TerminalException{1};
    }
    netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                            netCDF::NcFile::nc4);
    netCDF::NcDim nbDim = dataFile.addDim("nbPoint", nbPoint);
    std::vector<std::string> LDim{"nbPoint"};
    //
    netCDF::NcVar eVAR_lon = dataFile.addVar("lon", "double", LDim);
    netCDF::NcVar eVAR_lat = dataFile.addVar("lat", "double", LDim);
    std::vector<double> A(nbPoint);
    for (int u = 0; u < nbPoint; u++)
      A[u] = TotalListPairLL_b[u].eLon;
    eVAR_lon.putVar(A.data());
    for (int u = 0; u < nbPoint; u++)
      A[u] = TotalListPairLL_b[u].eLat;
    eVAR_lat.putVar(A.data());
    std::cerr << "After writing of the geography file\n";
  };
  if (DO_WNDMAG) {
    std::cerr << "Wind speed positions\n";
    FuncPrint(1);
  }
  if (DO_HS) {
    std::cerr << "Wave height positions\n";
    FuncPrint(2);
  }
}

void RAW_SCATTER_ALTIMETER(std::ostream &os,
                           std::vector<PairListWindWave> const &eSS,
                           NCLcaller<GeneralType> &eCall,
                           PermanentInfoDrawing const &ePerm,
                           int const &iGrid) {
  std::cerr << "Running RAW_SCATTER_ALTIMETER\n";
  int idWind = 1;
  int idWave = 2;
  SingleBlock eBlPLOT = ePerm.eFull.get_block("PLOT");
  bool UseDynamicRangeInScatter =
    eBlPLOT.get_bool("UseDynamicRangeInScatter");
  auto fPlot = [&](std::vector<PairMM> const &LPair, int const &idWindWave,
                   int const &idSat) -> void {
    DrawScatterArr eDrw;
    int siz = LPair.size();
    MyVector<double> eVectA(siz);
    MyVector<double> eVectB(siz);
    for (int i = 0; i < siz; i++) {
      PairMM ePair = LPair[i];
      eVectA[i] = ePair.Meas;
      eVectB[i] = ePair.Model;
    }
    std::vector<double> data_rangeA(2);
    std::vector<double> data_rangeB(2);
    //
    std::string SatName, SatNameFile;
    if (idSat == -1) {
      SatName = "all satellites";
      SatNameFile = "allsat";
    } else {
      SatName = GetNameOfSatelliteAltimeter(idSat);
      SatNameFile = SatName;
    }
    std::string VarType, eUnit;
    double TheMax;
    if (idWindWave == idWind) {
      VarType = "Wind";
      eUnit = "m/s";
      TheMax = 20;
    } else {
      VarType = "Wave";
      eUnit = "m";
      TheMax = 4;
    }
    if (UseDynamicRangeInScatter && siz > 0) {
      double maxA = eVectA.maxCoeff();
      double maxB = eVectB.maxCoeff();
      TheMax = std::max(maxA, maxB);
    }
    data_rangeA[0] = 0;
    data_rangeA[1] = TheMax;
    data_rangeB[0] = 0;
    data_rangeB[1] = TheMax;
    eDrw.VarNameAB_file = "Scatter_iGrid" + StringNumber(iGrid, 1) + "_Sat" +
                          SatNameFile + "_" + VarType;
    //
    eDrw.DoTitle = false;
    eDrw.AddStatMeasModel = true;
    eDrw.NameA_plot = "Data (" + eUnit + ")";
    eDrw.NameB_plot = "Model (" + eUnit + ")";
    eDrw.data_rangeA = data_rangeA;
    eDrw.data_rangeB = data_rangeB;
    eDrw.eVectA = eVectA;
    eDrw.eVectB = eVectB;
    eDrw.aSize = 100;
    eDrw.bSize = 100;
    PLOT_SCATTER(eDrw, eCall, ePerm);
  };
  SingleBlock const& BlPROCESS = ePerm.eFull.get_block("PROCESS");
  if (BlPROCESS.get_bool("DO_WNDMAG")) {
    std::vector<PairMM> TotalListPairWind;
    for (auto &eRec : eSS)
      TotalListPairWind.insert(TotalListPairWind.end(),
                               eRec.ListPairWind.begin(),
                               eRec.ListPairWind.end());
    fPlot(TotalListPairWind, idWind, -1);
    for (auto &eRec : eSS)
      fPlot(eRec.ListPairWind, idWind, eRec.eSat);
  }
  if (BlPROCESS.get_bool("DO_HS")) {
    std::vector<PairMM> TotalListPairWave;
    for (auto &eRec : eSS)
      TotalListPairWave.insert(TotalListPairWave.end(),
                               eRec.ListPairWave.begin(),
                               eRec.ListPairWave.end());
    fPlot(TotalListPairWave, idWave, -1);
    for (auto &eRec : eSS)
      fPlot(eRec.ListPairWave, idWave, eRec.eSat);
  }
}

void RAW_PLOT_ALTIMETER_TRACKS(std::ostream &os,
                               std::vector<SatelliteListTrack> const &LTrack,
                               bool const &PlotAllTrack,
                               bool const &PlotIndividualTrack,
                               NCLcaller<GeneralType> &eCall,
                               PermanentInfoDrawing const &ePerm) {
  std::cerr << "Running RAW_PLOT_ALTIMETER_TRACKS\n";
  double MinLon, MinLat, MaxLon, MaxLat;
  SingleBlock eBlSEL = ePerm.eFull.get_block("SELECT");
  if (eBlSEL.get_int("GEOSELECTION") == 1) {
    MinLon = eBlSEL.get_double("MinLON");
    MaxLon = eBlSEL.get_double("MaxLON");
    MinLat = eBlSEL.get_double("MinLAT");
    MaxLat = eBlSEL.get_double("MaxLAT");
  } else {
    std::vector<double> ListLon = eBlSEL.get_list_double("LONPOLY");
    std::vector<double> ListLat = eBlSEL.get_list_double("LATPOLY");
    int siz = ListLon.size();
    MinLon = ListLon[0];
    MaxLon = ListLon[0];
    MinLat = ListLat[0];
    MaxLat = ListLat[0];
    for (int i = 1; i < siz; i++) {
      double eLon = ListLon[i];
      double eLat = ListLat[i];
      if (MinLon > eLon)
        MinLon = eLon;
      if (MaxLon < eLon)
        MaxLon = eLon;
      if (MinLat > eLat)
        MinLat = eLat;
      if (MaxLat < eLat)
        MaxLat = eLat;
    }
  }
  QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};
  int nbSplitLon = 100;
  int nbSplitLat = 100;
  GridArray GrdArr = RECTANGULAR_GRID_ARRAY(eQuad, nbSplitLon, nbSplitLat);
  RecVar eRecVar = GetTrivialArrayPlot(GrdArr);
  //
  auto fPlot = [&](std::vector<SeqLineSegment> const &TheList,
                   int const &idSat) -> void {
    std::string SatName, SatNameFile;
    if (idSat == -1) {
      SatName = "all satellites";
      SatNameFile = "allsat";
    } else {
      SatName = GetNameOfSatelliteAltimeter(idSat);
      SatNameFile = SatName;
    }
    std::string TitleStr = "Tracks for " + SatName;
    std::string FileName = ePerm.eDir + "Tracks_for_" + SatNameFile;
    //
    DrawArr eDrw = ePerm.eDrawArr;
    eDrw.eQuadFrame = eQuad;
    eDrw.ColorMap = "WhBlGrYeRe";
    eDrw.DoColorBar = false;
    eDrw.TitleStr = TitleStr;
    eDrw.ListLineSegment = TheList;
    eDrw.VarNameUF = IntToString(idSat) + "_" + SatNameFile;
    //
    PLOT_PCOLOR(FileName, GrdArr, eDrw, eRecVar, eCall, ePerm);
  };
  struct SatelliteListTrackWrite {
    int eSat;
    std::vector<SeqLineSegment> ListLineSegment;
  };
  std::vector<SatelliteListTrackWrite> ListListLineSegment;
  struct SatelliteListTrackWriteIndividual {
    int eSat;
    double minTime;
    double maxTime;
    SeqLineSegment eLineSegment;
  };
  std::vector<SatelliteListTrackWriteIndividual> ListLineSegmentIndividual;
  for (auto &eRecTrack : LTrack) {
    std::vector<SeqLineSegment> eListLineSegment;
    for (auto &eListEnt : eRecTrack.ListListEntAltimeter) {
      std::vector<PairLL> ListPairLL;
      std::vector<double> ListTime;
      for (auto &eEnt : eListEnt) {
        double eLon = eEnt.Lon;
        double eLat = eEnt.Lat;
        double eTime = eEnt.Time;
        PairLL ePairLL{eLon, eLat};
        ListPairLL.push_back(ePairLL);
        ListTime.push_back(eTime);
      }
      SeqLineSegment eSeq{ListPairLL, false};
      eListLineSegment.push_back(eSeq);
      double minTime = VectorMin(ListTime);
      double maxTime = VectorMax(ListTime);
      ListLineSegmentIndividual.push_back(
          {eRecTrack.eSat, minTime, maxTime, eSeq});
    }
    ListListLineSegment.push_back({eRecTrack.eSat, eListLineSegment});
  }
  std::vector<SeqLineSegment> ListTotal;
  for (auto &eRec : ListListLineSegment)
    ListTotal.insert(ListTotal.end(), eRec.ListLineSegment.begin(),
                     eRec.ListLineSegment.end());
  if (PlotAllTrack) {
    fPlot(ListTotal, -1);
    for (auto &eRec : ListListLineSegment)
      fPlot(eRec.ListLineSegment, eRec.eSat);
  }
  if (PlotIndividualTrack) {
    std::vector<std::string> ListPossAlt = GetAllNamesOfSatelliteAltimeter();
    int nbPossAlt = ListPossAlt.size();
    std::cerr << "nbPossAlt=" << nbPossAlt << "\n";
    std::vector<int> PosAltimeter(nbPossAlt + 1, 0);
    std::cerr << "|ListLineSegmentIndividual|="
              << ListLineSegmentIndividual.size() << "\n";
    for (auto &eRec : ListLineSegmentIndividual) {
      std::string SatName = GetNameOfSatelliteAltimeter(eRec.eSat);
      double avgPairLength = AveragePairLength(eRec.eLineSegment);
      //      std::cerr << "eRec.eSat=" << eRec.eSat << "\n";
      PosAltimeter[eRec.eSat]++;
      int ePos = PosAltimeter[eRec.eSat];
      std::string strPos = StringNumber(ePos, 4);
      std::string strFirst = DATE_ConvertMjd2mystringPres(eRec.minTime);
      std::string strLast = DATE_ConvertMjd2mystringPres(eRec.maxTime);
      std::string TitleStr = "Tracks " + strPos + " for " + SatName + " from " +
                             strFirst + " to " + strLast +
                             " distPair=" + DoubleTo4dot1f(avgPairLength);
      std::string FileName = ePerm.eDir + "Track_" + SatName + "_at" + strPos;
      //
      DrawArr eDrw = ePerm.eDrawArr;
      eDrw.eQuadFrame = eQuad;
      eDrw.ColorMap = "WhBlGrYeRe";
      eDrw.DoColorBar = false;
      eDrw.TitleStr = TitleStr;
      eDrw.ListLineSegment = {eRec.eLineSegment};
      eDrw.VarNameUF =
          "Track_" + strPos + "_" + IntToString(eRec.eSat) + "_" + SatName;
      //
      PLOT_PCOLOR(FileName, GrdArr, eDrw, eRecVar, eCall, ePerm);
    }
    std::cerr << "Finished the individual plots\n";
  }
  std::cerr << "Leaving RAW_PLOT_ALTIMETER_TRACKS\n";
}

void RAW_PLOT_VALUE_TRACKS(std::ostream &os,
                           std::vector<SatelliteListTrack> const &LTrack,
                           std::vector<TotalArrGetData> const &ListTotalArr,
                           NCLcaller<GeneralType> &eCall,
                           PermanentInfoDrawing const &ePerm) {
  std::cerr << "Running RAW_PLOT_VALUE_TRACKS\n";
  int idWind = 1;
  int idWave = 2;
  SingleBlock const& BlPROCESS = ePerm.eFull.get_block("PROCESS");
  int MinEntryTrackPlot = BlPROCESS.get_int("MinEntryTrackPlot");
  bool PlotAddiWind = BlPROCESS.get_bool("PlotAddiWind");
  bool PlotAddiWave = BlPROCESS.get_bool("PlotAddiWave");
  PlotBound ePlotBound = ReadPlotBound(ePerm.eFull);
  int nbGrid = ListTotalArr.size();
  std::string ExtensionReal =
      GetRealInfoNclExtension(ePerm.Extension).eExtensionReal;
  auto fPlot = [&](MyVector<double> const &ListLat,
                   std::vector<PairLL> const &ListPairLL,
                   MyVector<double> const &ListMeas,
                   std::vector<MyVector<double>> const &ListListModel,
                   int const &idWindWave, int const &eSat, int const &iTrack,
                   double const &eTimeDay) -> void {
    std::string eVarName;
    double TheMin, TheMax;
    if (idWindWave == idWind) {
      eVarName = "Wind";
      TheMin = 0;
      TheMax = 20;
    } else {
      eVarName = "Wave";
      TheMin = 0;
      TheMax = 4;
    }
    //    std::cerr << "min(ListLat)=" << VectorMin(ListLat) << "\n";
    //    std::cerr << "max(ListLat)=" << VectorMax(ListLat) << "\n";
    std::string SatName = GetNameOfSatelliteAltimeter(eSat);
    std::string SatNameFile = SatName;
    //    std::cerr << "eTimeDay=" << eTimeDay << "\n";
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    //    std::cerr << "strPres=" << strPres << "\n";
    std::string strFile = DATE_ConvertMjd2mystringFile(eTimeDay);
    DrawLinesArr eDrawArr;
    eDrawArr.DoTitle = true;
    eDrawArr.TitleStr = "Track Nr" + IntToString(iTrack) + " of " + eVarName +
                        " for " + SatName + " at " + strPres;
    std::string fVarName =
        SatName + "_" + eVarName + "_Track" + StringNumber(iTrack, 4);
    eDrawArr.IsTimeSeries = false;
    eDrawArr.PairComparison = true;
    eDrawArr.DoExplicitLabel = false;
    eDrawArr.VarName = fVarName;
    eDrawArr.TheMax = TheMax;
    eDrawArr.TheMin = TheMin;
    eDrawArr.ListX = ListLat;
    std::vector<MyVector<double>> LLMeasModel{ListMeas};
    for (auto &eList : ListListModel)
      LLMeasModel.push_back(eList);
    eDrawArr.ListListVect = LLMeasModel;
    std::vector<std::string> LLStr{"meas."};
    for (int iGrid = 0; iGrid < static_cast<int>(ListListModel.size()); iGrid++)
      LLStr.push_back("model " + StringNumber(iGrid, 1));
    eDrawArr.ListName_plot = LLStr;
    std::string FileNameLinePlot = ePerm.eDir + "LINEPLOT_" + SatNameFile +
                                   "_" + eVarName + "_Track" +
                                   StringNumber(iTrack, 4) + "_at_" + strFile;
    std::string FileNameLinePlotExt = FileNameLinePlot + "." + ExtensionReal;
    LINES_PLOT(FileNameLinePlot, eDrawArr, eCall, ePerm);
    //
    // Now plotting if needed the wind
    //
    SeqLineSegment eSeq{ListPairLL, false};
    DrawArr eDrawArrQP = ePerm.eDrawArr;
    eDrawArrQP.ListLineSegment = {eSeq};

    //
    if (PlotAddiWind || PlotAddiWave) {
      std::vector<BashOper> ListBashOper;
      VarQuery eQuery{eTimeDay, -1, std::string("instant"),
                      static_cast<double>(0), std::string("direct")};
      for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
        eDrawArrQP.eQuadFrame = GetQuadArray(ListTotalArr[iGrid].GrdArr);
        auto fPlotVAR = [&](std::string const &eNature) -> void {
          // eNature is WIND10 or Hwave
          eDrawArrQP.VarNameUF = eNature;
          RecVar eRecVar = ModelSpecificVarSpecificTimeGeneral(
              ListTotalArr[iGrid], eNature, eQuery, ePlotBound);
          std::string FileName = ePerm.eDir + "LinePlotCompl_iGrid" +
                                 IntToString(iGrid) + "_" + SatNameFile + "_" +
                                 eNature + "_" + eRecVar.RecS.strAll;
          std::string FileNameMergeExt =
              ePerm.eDir + "LinePlotPcolor_iGrid" + IntToString(iGrid) + "_" +
              SatNameFile + "_" + eNature + "_" + eRecVar.RecS.strAll +
              "_merge." + ExtensionReal;
          std::cerr << "fPlotVAR FileName=" << FileName << "\n";
          std::string FileNameExt = FileName + "." + ExtensionReal;
          std::string TheCommand = "convert +append " + FileNameLinePlotExt +
                                   " " + FileNameExt + " " + FileNameMergeExt;
          int retValCorrect = 0; // convert operation returns 0 if all is ok.
          BashOper eBashOper{
              TheCommand, {FileNameLinePlotExt, FileNameExt}, retValCorrect};
          ListBashOper.push_back(eBashOper);
          //
          PLOT_PCOLOR_OR_QUIVER(FileName, ListTotalArr[iGrid].GrdArr,
                                eDrawArrQP, eRecVar, eCall, ePerm);
        };
        //
        if (PlotAddiWind)
          fPlotVAR("WIND10");
        if (PlotAddiWave)
          fPlotVAR("Hwave");
      }
      std::cerr << "|ListBashOper|=" << ListBashOper.size() << "\n";
      for (auto &eBashOper : ListBashOper) {
        std::cerr << "Submitting one eBashOper\n";
        GeneralType eGen(eBashOper);
        eCall.SubmitJob(eGen);
      }
    }
  };
  for (auto &eRec : LTrack) {
    int eSat = eRec.eSat;
    int iTrack = 0;
    for (auto &SingBlock : eRec.ListListEntAltimeter) {
      int TheSize = SingBlock.size();
      if (TheSize > MinEntryTrackPlot) {
        MyVector<double> ListLat(TheSize);
        std::vector<PairLL> ListPairLL(TheSize);
        MyVector<double> ListMeasWind(TheSize);
        MyVector<double> ListMeasWave(TheSize);
        MyVector<double> eV(TheSize);
        std::vector<MyVector<double>> ListListModelWind(nbGrid, eV);
        std::vector<MyVector<double>> ListListModelWave(nbGrid, eV);
        int idx = 0;
        double SumTimeDay = 0;
        for (auto &eEnt : SingBlock) {
          ListLat[idx] = eEnt.Lat;
          ListPairLL[idx] = {eEnt.Lon, eEnt.Lat};
          ListMeasWind(idx) = eEnt.WindSpeed_used;
          ListMeasWave(idx) = eEnt.Swh_used;
          for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
            ListListModelWind[iGrid](idx) = eEnt.WindSpeed_model[iGrid];
            ListListModelWave[iGrid](idx) = eEnt.Swh_model[iGrid];
          }
          SumTimeDay += eEnt.Time;
          // std::cerr << "idx=" << idx << " time=" << eEnt.Time << "\n";
          // std::cerr << "  SumTimeDay=" << SumTimeDay << "\n";
          idx++;
        }
        double eTimeDay = SumTimeDay / static_cast<double>(idx);
        // std::cerr << "eTimeDay=" << eTimeDay << "\n";
        SingleBlock const& BlPROCESS = ePerm.eFull.get_block("PROCESS");
        if (BlPROCESS.get_bool("DO_WNDMAG"))
          fPlot(ListLat, ListPairLL, ListMeasWind, ListListModelWind, idWind,
                eSat, iTrack, eTimeDay);
        if (BlPROCESS.get_bool("DO_HS"))
          fPlot(ListLat, ListPairLL, ListMeasWave, ListListModelWave, idWave,
                eSat, iTrack, eTimeDay);
        iTrack++;
      }
    }
  }
}

std::vector<std::vector<PairListWindWave>>
MergeTracksForRawStatistics(std::vector<SatelliteListTrack> const &LTrack,
                            FullNamelist const &eFull) {
  std::cerr << "MergeTracksForRawStatistics, |LTrack|=" << LTrack.size()
            << "\n";
  SingleBlock const& eBlSELECT = eFull.get_block("SELECT");
  SingleBlock const& eBlPROC = eFull.get_block("PROC");
  double MinWind = eBlSELECT.get_double("MinWIND");
  double MaxWind = eBlSELECT.get_double("MaxWIND");
  int nbGrid = eBlPROC.get_list_string("ListMODELNAME").size();
  std::vector<std::string> ListMinHs_measStr =
    eBlSELECT.get_list_string("ListMinHS_meas");
  std::vector<std::string> ListMaxHs_measStr =
    eBlSELECT.get_list_string("ListMaxHS_meas");
  std::vector<std::string> ListFootprintKM_str =
    eBlSELECT.get_list_string("ListRadiusFootprintKM");
  std::cerr << "ListMinHs_MeasStr =";
  for (auto &eStr : ListMinHs_measStr)
    std::cerr << " " << eStr;
  std::cerr << "\n";
  std::cerr << "ListMaxHs_MeasStr =";
  for (auto &eStr : ListMaxHs_measStr)
    std::cerr << " " << eStr;
  std::cerr << "\n";

  std::vector<double> ListMinHs_meas =
      ConvertListStringValueToVector_SAT(ListMinHs_measStr, 0);
  std::vector<double> ListMaxHs_meas =
      ConvertListStringValueToVector_SAT(ListMaxHs_measStr, 998);

  std::cerr << "ListMinHs_meas =";
  for (auto &eStr : ListMinHs_meas)
    std::cerr << " " << eStr;
  std::cerr << "\n";
  std::cerr << "ListMaxHs_meas =";
  for (auto &eStr : ListMaxHs_meas)
    std::cerr << " " << eStr;
  std::cerr << "\n";

  std::vector<double> ListFootprintKM =
      ConvertListStringValueToVector_SAT(ListFootprintKM_str, 40);

  // Default value is for a fooprint of 40KM

  double MinHs_model = eBlSELECT.get_double("MinHS_model");
  double MaxHs_model = eBlSELECT.get_double("MaxHS_model");

  std::vector<std::vector<PairListWindWave>> RetList(nbGrid);
  for (auto &eRecTrack : LTrack) {
    std::vector<std::vector<PairMM>> ListPairWind(nbGrid);
    std::vector<std::vector<PairMM>> ListPairWave(nbGrid);
    std::vector<std::vector<double>> ListTimeWind(nbGrid);
    std::vector<std::vector<double>> ListTimeWave(nbGrid);
    std::vector<std::vector<PairLL>> ListPairLLwind(nbGrid);
    std::vector<std::vector<PairLL>> ListPairLLwave(nbGrid);
    int nbMatch = 0;
    for (auto &eListEnt : eRecTrack.ListListEntAltimeter)
      for (auto &eEnt : eListEnt) {
        nbMatch++;
        int eSatellite = eEnt.Satellite;
        // std::cerr << "eSatellite=" << eSatellite << "\n";
        double eMinHs_meas = ListMinHs_meas[eSatellite - 1];
        double eMaxHs_meas = ListMaxHs_meas[eSatellite - 1];
        if (std::isnan(eMinHs_meas) || std::isnan(eMaxHs_meas)) {
          std::cerr << "eMinHs_meas=" << eMinHs_meas << "\n";
          std::cerr << "eMaxHs_meas=" << eMaxHs_meas << "\n";
          std::cerr << "None of them should be a NaN\n";
          throw TerminalException{1};
        }
        if (eEnt.Swh_used < eMaxHs_meas && eEnt.Swh_used > eMinHs_meas) {
          bool IsOK = true;
          for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
            double eSwh = eEnt.Swh_model[iGrid];
            if (eSwh > MaxHs_model || eSwh < MinHs_model)
              IsOK = false;
          }
          if (IsOK) {
            for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
              PairMM ePairHs{eEnt.Swh_used, eEnt.Swh_model[iGrid]};
              ListTimeWave[iGrid].push_back(eEnt.Time);
              ListPairWave[iGrid].push_back(ePairHs);
              ListPairLLwave[iGrid].push_back({eEnt.Lon, eEnt.Lat});
            }
          }
        }
        if (eEnt.WindSpeed_used < MaxWind && eEnt.WindSpeed_used > MinWind) {
          for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
            PairMM ePairWind{eEnt.WindSpeed_used, eEnt.WindSpeed_model[iGrid]};
            ListTimeWind[iGrid].push_back(eEnt.Time);
            ListPairWind[iGrid].push_back(ePairWind);
            ListPairLLwind[iGrid].push_back({eEnt.Lon, eEnt.Lat});
          }
        }
      }
    int nbCorrWind = ListPairWind[0].size();
    int nbCorrWave = ListPairWave[0].size();
    std::string SatName = GetNameOfSatelliteAltimeter(eRecTrack.eSat);
    std::cerr << "Sat=" << SatName << " nbMatch=" << nbMatch
              << "  nbCorr(wind/wave)=" << nbCorrWind << "/" << nbCorrWave
              << "\n";
    for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
      PairListWindWave eSS{eRecTrack.eSat,       ListTimeWind[iGrid],
                           ListTimeWave[iGrid],  ListPairWind[iGrid],
                           ListPairWave[iGrid],  ListPairLLwind[iGrid],
                           ListPairLLwave[iGrid]};
      RetList[iGrid].push_back(eSS);
    }
  }
  std::cerr << "MergeTracksForRawStatistics, |LTrack|=" << LTrack.size()
            << "\n";
  return RetList;
}

std::vector<PairLL> ReadLonLatDiscFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in ReadLonLatDiscFile\n";
    std::cerr << "The file eFile=" << eFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  int nbLine = FILE_GetNumberLine(eFile);
  std::ifstream is;
  is.open(eFile);
  std::vector<PairLL> ListPt(nbLine);
  for (int iLine = 0; iLine < nbLine; iLine++) {
    double eLon, eLat;
    is >> eLon;
    is >> eLat;
    PairLL ePt{eLon, eLat};
    ListPt[iLine] = ePt;
  }
  return ListPt;
}

std::vector<int>
FilterByMinDistCoast(std::vector<SingleEntryMeasurement> const &eList,
                     FullNamelist const &eFull) {
  int nbEntry = eList.size();
  std::vector<int> ListStatus(nbEntry, 1);
  SingleBlock BlSELECT = eFull.get_block("SELECT");
  std::string eFileCoast = BlSELECT.get_string("LonLatDiscFile");
  double MinDistCoastKM = BlSELECT.get_double("MinDistCoastKM");
  std::cerr << "eFileCoast=" << eFileCoast << "\n";
  std::vector<PairLL> ListPtCoast = ReadLonLatDiscFile(eFileCoast);
  std::cerr << "|ListPtCoast|=" << ListPtCoast.size() << "\n";
  std::vector<PairLL> ListPt(nbEntry);
  for (int iEntry = 0; iEntry < nbEntry; iEntry++) {
    double eLon = eList[iEntry].Lon;
    double eLat = eList[iEntry].Lat;
    PairLL ePt{eLon, eLat};
    ListPt[iEntry] = ePt;
  }
  std::vector<double> ListMinDist =
      GetListMinimalDistances(ListPtCoast, ListPt);
  std::vector<SingleEntryMeasurement> RetList;
  for (int iEntry = 0; iEntry < nbEntry; iEntry++)
    if (ListMinDist[iEntry] < MinDistCoastKM)
      ListStatus[iEntry] = 0;
  return ListStatus;
}

std::vector<int>
FilterByOvervalues(std::vector<SingleEntryMeasurement> const &eList,
                   FullNamelist const &eFull) {
  SingleBlock const& BlPROCESS = eFull.get_block("PROCESS");
  bool DO_WNDMAG = BlPROCESS.get_bool("DO_WNDMAG");
  bool DO_HS = BlPROCESS.get_bool("DO_HS");
  std::cerr << "FilterByOvervalues DO_WNDMAG=" << DO_WNDMAG
            << " DO_HS=" << DO_HS << "\n";
  int nbEntry = eList.size();
  std::vector<int> ListStatus(nbEntry, 0);
  double LargeValue = 30000;
  int nbCorrWind = 0;
  int nbCorrWave = 0;
  int nbCorrWindWave = 0;
  for (int iEntry = 0; iEntry < nbEntry; iEntry++) {
    SingleEntryMeasurement eEnt = eList[iEntry];
    if (eEnt.WindSpeed_used < LargeValue)
      nbCorrWind++;
    if (eEnt.Swh_used < LargeValue)
      nbCorrWave++;
    if (eEnt.Swh_used < LargeValue && eEnt.WindSpeed_used < LargeValue)
      nbCorrWindWave++;
    bool IsCorrect = true;
    if (DO_WNDMAG) {
      if (eEnt.WindSpeed_used > LargeValue)
        IsCorrect = false;
    }
    if (DO_HS) {
      if (eEnt.Swh_used > LargeValue)
        IsCorrect = false;
    }
    if (IsCorrect)
      ListStatus[iEntry] = 1;
  }
  std::cerr << "nbCorrWind=" << nbCorrWind << " nbCorrWave=" << nbCorrWave
            << " nbCorrWindWave=" << nbCorrWindWave << "\n";
  return ListStatus;
}

std::vector<int>
MergeListListStatus(std::vector<std::vector<int>> const &ListListStatus,
                    std::vector<SingleEntryMeasurement> const &eVectEnt) {
  int siz = eVectEnt.size();
  std::vector<int> ListStatus(siz, 1);
  for (auto &eVect : ListListStatus)
    for (int i = 0; i < siz; i++)
      ListStatus[i] *= eVect[i];
  return ListStatus;
}

std::vector<SingleEntryMeasurement>
SelectByStatus(std::vector<std::vector<int>> const &ListListStatus,
               std::vector<SingleEntryMeasurement> const &eVectEnt) {
  int siz = eVectEnt.size();
  std::vector<SingleEntryMeasurement> RetVect;
  int nbStatus = ListListStatus.size();
  for (int i = 0; i < siz; i++) {
    int eStatus = 1;
    for (int iStatus = 0; iStatus < nbStatus; iStatus++)
      eStatus *= ListListStatus[iStatus][i];
    if (eStatus == 1)
      RetVect.push_back(eVectEnt[i]);
  }
  return RetVect;
}

std::vector<SingleEntryMeasurement>
FilteringData(std::vector<SingleEntryMeasurement> const &PreListSingleEntry,
              FullNamelist const &eFull) {
  std::vector<std::vector<int>> ListListStatus;
  auto InsertListStatus = [&](std::vector<int> const &eListStatus,
                              std::string const &meth) -> void {
    int siz = eListStatus.size();
    int eSum = 0;
    for (int i = 0; i < siz; i++)
      eSum += eListStatus[i];
    std::cerr << "meth = " << meth << "  siz=" << siz << " sum=" << eSum
              << "\n";
    ListListStatus.push_back(eListStatus);
  };
  SingleBlock const& eBlSELECT = eFull.get_block("SELECT");
  InsertListStatus(FilterByOvervalues(PreListSingleEntry, eFull),
                   "max over value");
  if (eBlSELECT.get_bool("DoMinDistCoast"))
    InsertListStatus(FilterByMinDistCoast(PreListSingleEntry, eFull),
                     "min dist coast");
  if (eBlSELECT.get_bool("EliminationShortTrack"))
    InsertListStatus(GetListStatusTrackLength(PreListSingleEntry, eFull),
                     "elim short track");
  if (eBlSELECT.get_bool("EliminationWrongModelInterpolationValues"))
    InsertListStatus(
        GetEliminationWrongModelInterpolation(PreListSingleEntry, eFull),
        "wrong model interpolation");
  return SelectByStatus(ListListStatus, PreListSingleEntry);
}

FullNamelist NAMELIST_GetStandardALTIMETRY_COMPARISON() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListListStringValues1["ListMODELNAME"] = {};
  ListListStringValues1["ListGridFile"] = {};
  ListListStringValues1["ListHisPrefix"] = {};
  ListStringValues1["BoundFile"] = "unset";
  ListStringValues1["PicPrefix"] = "unset PicPrefix";
  ListStringValues1["Extension"] = "png";
  ListStringValues1["__NaturePlot"] = "ALTIMETRY";
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Lines_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  ListListStringValues1["ListTypeData"] = {};
  ListListStringValues1["ListDirData"] = {};
  ListIntValues1["NPROC"] = 1;
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues1);
  BlockPROC.setListBoolValues(ListBoolValues1);
  BlockPROC.setListDoubleValues(ListDoubleValues1);
  BlockPROC.setListListDoubleValues(ListListDoubleValues1);
  BlockPROC.setListStringValues(ListStringValues1);
  BlockPROC.setListListStringValues(ListListStringValues1);
  ListBlock["PROC"] = BlockPROC;
  // SELECT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<int>> ListListIntValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["GEOSELECTION"] = 1;
  ListDoubleValues2["MinLON"] = -7;
  ListDoubleValues2["MaxLON"] = 37;
  ListDoubleValues2["MinLAT"] = 30;
  ListDoubleValues2["MaxLAT"] = 46;
  ListDoubleValues2["MaxDistTrackPointKM"] = 20;
  ListDoubleValues2["MaxDeltaTimeTrackPointDay"] = 0.001;
  ListDoubleValues2["MinWIND"] = 0;
  ListDoubleValues2["MaxWIND"] = 998;
  ListDoubleValues2["MinHS_model"] = 0;
  ListDoubleValues2["MaxHS_model"] = 998;
  ListDoubleValues2["minQCwave"] = -10000;
  ListDoubleValues2["maxQCwave"] = 10000;
  ListDoubleValues2["tolLL"] = 0;
  ListStringValues2["BEGTC"] = "20110915.000000";
  ListStringValues2["ENDTC"] = "20110925.000000";
  ListListDoubleValues2["LONPOLY"] = {10, 10, 10};
  ListListDoubleValues2["LATPOLY"] = {10, 10, 10};
  ListListDoubleValues2["ListBinWindValues"] = {};
  ListListDoubleValues2["ListBinWaveValues"] = {};
  ListListStringValues2["ListMinHS_meas"] = {};
  ListListStringValues2["ListMaxHS_meas"] = {};
  ListListStringValues2["ListRadiusFootprintKM"] = {};
  ListListStringValues2["AllowedSatellites"] =
      GetAllNamesOfSatelliteAltimeter();
  ListIntValues2["MethodTrackSmoothing"] =
      0; // 0 for no smoothing, 1 for smoothing with length given by model
         // resolution and 2 for smoothing with length given by sat footprint
  ListBoolValues2["EliminationShortTrack"] = false;
  ListBoolValues2["EliminationWrongModelInterpolationValues"] = true;
  ListBoolValues2["SelectingHours"] = false;
  ListListIntValues2["ListAllowedHours"] = {0,  1,  2,  3,  4,  5,  6,  7,
                                            8,  9,  10, 11, 12, 13, 14, 15,
                                            16, 17, 18, 19, 20, 21, 22, 23};
  ListIntValues2["MinimalTrackSize"] = 30;
  ListBoolValues2["DoMinDistCoast"] = false;
  ListDoubleValues2["MinDistCoastKM"] = 60;
  ListStringValues2["LonLatDiscFile"] = "LonLatDisc.txt";
  ListStringValues2["InterpolationMethod"] = "linear_interpolation";
  ListListDoubleValues2["ListLineLonStart"] = {};
  ListListDoubleValues2["ListLineLatStart"] = {};
  ListListDoubleValues2["ListLineLonEnd"] = {};
  ListListDoubleValues2["ListLineLatEnd"] = {};
  SingleBlock BlockSELECT;
  BlockSELECT.setListIntValues(ListIntValues2);
  BlockSELECT.setListBoolValues(ListBoolValues2);
  BlockSELECT.setListDoubleValues(ListDoubleValues2);
  BlockSELECT.setListListDoubleValues(ListListDoubleValues2);
  BlockSELECT.setListListIntValues(ListListIntValues2);
  BlockSELECT.setListStringValues(ListStringValues2);
  BlockSELECT.setListListStringValues(ListListStringValues2);
  ListBlock["SELECT"] = BlockSELECT;
  // PROCESS
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListBoolValues3["DO_WNDMAG"] = false;
  ListBoolValues3["DO_HS"] = false;
  ListBoolValues3["DO_STAT"] = false;
  ListBoolValues3["DO_MONTHLY_STAT"] = false;
  ListBoolValues3["DO_HOURLY_STAT"] = false;
  ListBoolValues3["DO_BIN_RANGING_STAT"] = false;
  ListBoolValues3["DO_GEOGRAPHIC_STAT"] = false;
  ListBoolValues3["WRITE_MEASUREMENT_GEOG_POSITION"] = false;
  ListBoolValues3["PLOT_STAT_GEOGRAPHIC_BREAKDOWN"] = false;
  ListBoolValues3["PlotIndividualPoints"] = false;
  ListBoolValues3["PlotGeographicStat"] = false;
  ListBoolValues3["DO_TXTRAW"] = false;
  ListBoolValues3["DO_NCOUT"] = false;
  ListBoolValues3["PLOT_ALL_TRACKS"] = false;
  ListBoolValues3["PLOT_INDIVIDUAL_TRACKS"] = false;
  ListBoolValues3["PLOT_TRACKS"] = false;
  ListBoolValues3["USE_CORRECTED"] = true;
  ListBoolValues3["DO_SCATTERPLOT"] = true;
  ListBoolValues3["SPATIALAVER"] = false;
  ListBoolValues3["PlotAddiWind"] = false;
  ListBoolValues3["PlotAddiWave"] = false;
  ListIntValues3["MinEntryTrackPlot"] = 30;
  ListStringValues3["FILE_SAVE_TXT"] = "alldata.txt";
  ListBoolValues3["DO_SAVE_NC"] = false;
  ListStringValues3["FILE_SAVE_NC"] = "NeededFileName.nc";
  SingleBlock BlockPROCESS;
  BlockPROCESS.setListIntValues(ListIntValues3);
  BlockPROCESS.setListBoolValues(ListBoolValues3);
  BlockPROCESS.setListDoubleValues(ListDoubleValues3);
  BlockPROCESS.setListStringValues(ListStringValues3);
  BlockPROCESS.setListListStringValues(ListListStringValues3);
  ListBlock["PROCESS"] = BlockPROCESS;
  // PLOT
  std::map<std::string, int> ListIntValues4;
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, double> ListDoubleValues4;
  std::map<std::string, std::string> ListStringValues4;
  std::map<std::string, std::vector<std::string>> ListListStringValues4;
  std::map<std::string, std::vector<double>> ListListDoubleValues4;
  ListIntValues4["nbBlock"] = 1;
  ListIntValues4["nbLevelSpa"] = 20;
  ListIntValues4["nbLabelStride"] = 10;
  ListBoolValues4["DoColorBar"] = true;
  ListBoolValues4["DoTitle"] = true;
  ListBoolValues4["FillLand"] = true;
  ListBoolValues4["PrintMMA"] = false;
  ListBoolValues4["DrawRiver"] = false;
  ListBoolValues4["DrawContourBathy"] = false;
  ListBoolValues4["DrawAnnotation"] = false;
  ListBoolValues4["DoTitleString"] = true;
  ListBoolValues4["cnSmoothingOn"] = true;
  ListBoolValues4["UseNativeGrid"] = true;
  ListBoolValues4["UseDynamicRangeInScatter"] = false;
  ListBoolValues4["VariableRange"] = false;
  ListListStringValues4["BoundSingle_var"] = {};
  ListListDoubleValues4["BoundSingle_min"] = {};
  ListListDoubleValues4["BoundSingle_max"] = {};
  ListDoubleValues4["vcRefLengthF"] = 0.02;
  ListListStringValues2["ListSubstitution"] = {};
  ListDoubleValues4["AnnotationLon"] = -400;
  ListDoubleValues4["AnnotationLat"] = -400;
  ListStringValues4["AnnotationText"] = "unset_text";
  ListStringValues4["FileDirectNCLins"] = "irrelevant";
  ListStringValues4["GridResolution"] = "HighRes";
  ListStringValues4["ColorMap"] = "BlAqGrYeOrReVi200";
  ListStringValues4["cnFillMode"] = "RasterFill";
  ListBoolValues4["cnFillOn"] = true;
  ListBoolValues4["cnLinesOn"] = false;
  ListBoolValues4["cnLineLabelsOn"] = false;
  ListStringValues4["LandPortr"] = "Landscape";
  ListStringValues4["optStatStr"] = "double";
  SingleBlock BlockPLOT;
  BlockPLOT.setListIntValues(ListIntValues4);
  BlockPLOT.setListBoolValues(ListBoolValues4);
  BlockPLOT.setListDoubleValues(ListDoubleValues4);
  BlockPLOT.setListStringValues(ListStringValues4);
  BlockPLOT.setListListStringValues(ListListStringValues4);
  BlockPLOT.setListListDoubleValues(ListListDoubleValues4);
  ListBlock["PLOT"] = BlockPLOT;
  // Merging all data
  return FullNamelist(ListBlock);
}

void Process_Altimetry_Comparison_Request(FullNamelist const &eFull) {
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);

  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  SingleBlock const& eBlPROC = eFull.get_block("PROC");
  std::vector<std::string> ListModelName =
    eBlPROC.get_list_string("ListMODELNAME");
  std::vector<std::string> ListGridFile =
    eBlPROC.get_list_string("ListGridFile");
  std::vector<std::string> ListHisPrefix =
    eBlPROC.get_list_string("ListHisPrefix");
  int nbGrid = ListModelName.size();
  std::vector<double> ListAvgDist_model(nbGrid);
  std::vector<TotalArrGetData> ListTotalArr = RealAllArrayHistory(eFull);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    double avgDistKM_model = GetGridSpacing(ListTotalArr[iGrid].GrdArr);
    std::cerr << "iGrid=" << iGrid << " avgDistKM_model=" << avgDistKM_model
              << "\n";
    ListAvgDist_model[iGrid] = avgDistKM_model;
  }
  double avgDistKM_model = VectorMax(ListAvgDist_model);
  std::cerr << "Final avgDistKM_model=" << avgDistKM_model << "\n";
  //
  std::vector<std::string> AllSatNames = GetAllNamesOfSatelliteAltimeter();
  std::cerr << "AllSatNames =";
  for (auto &eStr : AllSatNames)
    std::cerr << " " << eStr;
  std::cerr << "\n";
  //
  std::cerr << "Before RETRIEVE_RELEVANT_ALTI_DATA\n";
  SatelliteSerInfo eRecSer = RetrieveTimeInformation(ListTotalArr, eFull);
  std::cerr << "After RetrieveTimeInformation\n";
  std::vector<SingleEntryMeasurement> PreListSingleEntry =
      RETRIEVE_RELEVANT_ALTI_DATA(eRecSer, eFull);
  std::cerr << "After RETRIEVE_RELEVANT_Alti_DATA\n";
  InterpolateAltimeterData(PreListSingleEntry, ListTotalArr, eFull);
  std::cerr << "|PreListSingleEntry| = " << PreListSingleEntry.size() << "\n";
  std::set<int> SatelliteId = GetListSatelliteId_set(PreListSingleEntry, eFull);
  std::cerr << "Ater GetListSatelliteId_set\n";
  SingleBlock eBlSELECT = eFull.get_block("SELECT");
  std::vector<SingleEntryMeasurement> ListSingleEntry =
      FilteringData(PreListSingleEntry, eFull);
  std::cerr << "|ListSingleEntry|=" << ListSingleEntry.size() << "\n";
  std::vector<SatelliteListTrack> ListTrackInfo =
      GetListTrackAltimeter(ListSingleEntry, avgDistKM_model, eFull);
  std::vector<std::vector<PairListWindWave>> ListSS =
      MergeTracksForRawStatistics(ListTrackInfo, eFull);
  SingleBlock eBlPROCESS = eFull.get_block("PROCESS");
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    if (eBlPROCESS.get_bool("DO_STAT"))
      RAW_STATISTICS_ALTIMETER(ListSS[iGrid], eRecSer, ePerm, iGrid);
    if (eBlPROCESS.get_bool("DO_MONTHLY_STAT"))
      RAW_STATISTICS_BREAKDOWN_ALTIMETER(ListSS[iGrid], eRecSer, ePerm, iGrid,
                                         "monthly");
    if (eBlPROCESS.get_bool("DO_HOURLY_STAT"))
      RAW_STATISTICS_BREAKDOWN_ALTIMETER(ListSS[iGrid], eRecSer, ePerm, iGrid,
                                         "hourly");
    if (eBlPROCESS.get_bool("DO_GEOGRAPHIC_STAT"))
      RAW_STATISTICS_BREAKDOWN_ALTIMETER(ListSS[iGrid], eRecSer, ePerm, iGrid,
                                         "geography");
    if (eBlPROCESS.get_bool("DO_BIN_RANGING_STAT"))
      RAW_BIN_RANGING_STATISTICS_BREAKDOWN(ListSS[iGrid], eRecSer, ePerm,
                                           iGrid);
    if (eBlPROCESS.get_bool("DO_SCATTERPLOT"))
      RAW_SCATTER_ALTIMETER(std::cout, ListSS[iGrid], eCall, ePerm, iGrid);
    if (eBlPROCESS.get_bool("WRITE_MEASUREMENT_GEOG_POSITION"))
      WRITE_MEASUREMENT_GEOG_POSITION(ListSS[iGrid], eRecSer, ePerm, iGrid);
    if (eBlPROCESS.get_bool("PLOT_STAT_GEOGRAPHIC_BREAKDOWN"))
      BREAKDOWN_GEOG_POINT(ListSS[iGrid], eRecSer, eCall, ePerm, iGrid);
  }
  if (eBlPROCESS.get_bool("DO_TXTRAW"))
    PRINT_RAW_DATA(ListSS, ListTrackInfo, eRecSer, eCall, ePerm);
  bool PlotAllTrack = eBlPROCESS.get_bool("PLOT_ALL_TRACKS");
  bool PlotIndividualTrack =
    eBlPROCESS.get_bool("PLOT_INDIVIDUAL_TRACKS");
  if (PlotAllTrack || PlotIndividualTrack) {
    RAW_PLOT_ALTIMETER_TRACKS(std::cout, ListTrackInfo, PlotAllTrack,
                              PlotIndividualTrack, eCall, ePerm);
  }
  if (eBlPROCESS.get_bool("PLOT_TRACKS"))
    RAW_PLOT_VALUE_TRACKS(std::cout, ListTrackInfo, ListTotalArr, eCall, ePerm);
}

FullNamelist NAMELIST_Comparison_Altimetry_Source() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["PicPrefix"] = "unset PicPrefix";
  ListStringValues1["Extension"] = "png";
  ListStringValues1["__NaturePlot"] = "ALTIMETRY";
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  ListListStringValues1["ListTypeData"] = {};
  ListListStringValues1["ListDirData"] = {};
  ListIntValues1["NPROC"] = 1;
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues1);
  BlockPROC.setListBoolValues(ListBoolValues1);
  BlockPROC.setListDoubleValues(ListDoubleValues1);
  BlockPROC.setListListDoubleValues(ListListDoubleValues1);
  BlockPROC.setListStringValues(ListStringValues1);
  BlockPROC.setListListStringValues(ListListStringValues1);
  ListBlock["PROC"] = BlockPROC;
  // SELECT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["GEOSELECTION"] = 1;
  ListDoubleValues2["MinLON"] = -7;
  ListDoubleValues2["MaxLON"] = 37;
  ListDoubleValues2["MinLAT"] = 30;
  ListDoubleValues2["MaxLAT"] = 46;
  ListListDoubleValues2["LONPOLY"] = {10, 10, 10};
  ListListDoubleValues2["LATPOLY"] = {10, 10, 10};
  ListDoubleValues2["MaxDistTrackPointKM"] = 20;
  ListDoubleValues2["MaxDeltaTimeTrackPointDay"] = 0.001;
  ListDoubleValues2["MinWIND"] = 0;
  ListDoubleValues2["MaxWIND"] = 998;
  ListListStringValues2["ListMinHS_meas"] = {};
  ListListStringValues2["ListMaxHS_meas"] = {};
  ListListStringValues2["ListRadiusFootprint"] = {};
  ListDoubleValues2["MinHS_model"] = 0;
  ListDoubleValues2["MaxHS_model"] = 998;
  ListStringValues2["BEGTC"] = "20110915.000000";
  ListStringValues2["ENDTC"] = "20110925.000000";
  ListListStringValues2["AllowedSatellites"] =
      GetAllNamesOfSatelliteAltimeter();
  ListIntValues2["MethodTrackSmoothing"] =
      0; // 0 for no smoothing, 1 for smoothing with length given by sat
         // footprint and 2 for smoothing with length given by model resolution
  ListBoolValues2["EliminationShortTrack"] = false;
  ListBoolValues2["EliminationOutsideGrid"] = false;
  ListIntValues2["MinimalTrackSize"] = 30;
  ListBoolValues2["DoMinDistCoast"] = false;
  ListDoubleValues2["MinDistCoastKM"] = 60;
  ListStringValues2["LonLatDiscFile"] = "LonLatDisc.txt";
  SingleBlock BlockSELECT;
  BlockSELECT.setListIntValues(ListIntValues2);
  BlockSELECT.setListBoolValues(ListBoolValues2);
  BlockSELECT.setListDoubleValues(ListDoubleValues2);
  BlockSELECT.setListListDoubleValues(ListListDoubleValues2);
  BlockSELECT.setListStringValues(ListStringValues2);
  BlockSELECT.setListListStringValues(ListListStringValues2);
  ListBlock["SELECT"] = BlockSELECT;
  // PROCESS
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListBoolValues3["USE_CORRECTED"] = true;
  ListDoubleValues3["tolLL"] = 0.01;
  ListDoubleValues3["tolTime"] = 0.0001;
  ListBoolValues3["DO_STAT"] = false;
  ListBoolValues3["DO_TXTRAW"] = false;
  ListBoolValues3["DO_NCOUT"] = false;
  ListBoolValues3["DO_SCATTERPLOT"] = false;
  ListBoolValues3["PLOT_ALL_TRACKS"] = false;
  ListBoolValues3["PLOT_INDIVIDUAL_TRACKS"] = false;
  ListIntValues3["MinEntryTrackPlot"] = 30;
  ListStringValues3["FILE_SAVE_TXT"] = "alldata.txt";
  ListBoolValues3["DO_SAVE_NC"] = false;
  ListStringValues3["FILE_SAVE_NC"] = "NeededFileName.nc";
  SingleBlock BlockPROCESS;
  BlockPROCESS.setListIntValues(ListIntValues3);
  BlockPROCESS.setListBoolValues(ListBoolValues3);
  BlockPROCESS.setListDoubleValues(ListDoubleValues3);
  BlockPROCESS.setListStringValues(ListStringValues3);
  BlockPROCESS.setListListStringValues(ListListStringValues3);
  ListBlock["PROCESS"] = BlockPROCESS;
  // Merging all data
  return FullNamelist(ListBlock);
}

void Process_Comparison_Altimetry_Sources(FullNamelist const &eFull) {
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  SingleBlock BlPLOT = eFull.get_block("PLOT");
  SingleBlock BlPROC = eFull.get_block("PROC");
  SingleBlock BlSELECT = eFull.get_block("SELECT");
  SingleBlock BlPROCESS = eFull.get_block("PROCESS");
  std::vector<std::string> ListTypeData =
    BlPROC.get_list_string("ListTypeData");
  std::vector<std::string> ListDirData =
    BlPROC.get_list_string("ListDirData");
  int nbType = ListTypeData.size();
  if (nbType != 2) {
    std::cerr << "nbType = " << nbType << "\n";
    std::cerr
        << "The comparison of altimetry sources can only be done if nbType=2\n";
    throw TerminalException{1};
  }
  SatelliteSerInfo eRecSer = RetrieveTimeInformation_Begin_End(eFull);
  if (BlSELECT.get_bool("EliminationOutsideGrid")) {
    std::cerr << "EliminationOutsideGrid should be set to false\n";
    std::cerr << "for doing ComparisonAltimetrySource\n";
    throw TerminalException{1};
  }
  //
  FullNamelist eFull0 = eFull;
  std::string eName0 = ListTypeData[0];
  std::string eDir0 = ListDirData[0];
  eFull0.get_block_mut("PROC").get_list_string_mut("ListTypeData") = {eName0};
  eFull0.get_block_mut("PROC").get_list_string_mut("ListDirData") = {eDir0};
  std::vector<SingleEntryMeasurement> PreListSingleEntry0 =
      RETRIEVE_RELEVANT_ALTI_DATA(eRecSer, eFull0);
  std::vector<SingleEntryMeasurement> ListSingleEntry0 =
      FilteringData(PreListSingleEntry0, eFull);
  //
  FullNamelist eFull1 = eFull;
  std::string eName1 = ListTypeData[1];
  std::string eDir1 = ListDirData[1];
  eFull1.get_block_mut("PROC").get_list_string_mut("ListTypeData") = { eName1};
  eFull1.get_block_mut("PROC").get_list_string_mut("ListDirData") = {eDir1};
  std::vector<SingleEntryMeasurement> PreListSingleEntry1 =
      RETRIEVE_RELEVANT_ALTI_DATA(eRecSer, eFull1);
  std::vector<SingleEntryMeasurement> ListSingleEntry1 =
      FilteringData(PreListSingleEntry1, eFull);
  //
  auto Comparison = [&](SingleEntryMeasurement const &x,
                        SingleEntryMeasurement const &y) -> bool {
    if (x.Satellite < y.Satellite)
      return true;
    if (x.Satellite > y.Satellite)
      return false;
    //
    if (x.Time < y.Time)
      return true;
    if (x.Time > y.Time)
      return false;
    //
    if (x.Lon < y.Lon)
      return true;
    if (x.Lon > y.Lon)
      return false;
    //
    if (x.Lat < y.Lat)
      return true;
    if (x.Lat > y.Lat)
      return false;
    //
    return false;
  };
  std::sort(ListSingleEntry0.begin(), ListSingleEntry0.end(), Comparison);
  std::sort(ListSingleEntry1.begin(), ListSingleEntry1.end(), Comparison);
  int nbEntry0 = ListSingleEntry0.size();
  int nbEntry1 = ListSingleEntry1.size();
  std::cerr << "nbEntry0=" << nbEntry0 << " nbEntry1=" << nbEntry1 << "\n";
  //
  SingleBlock BlockPROCESS = eFull.get_block("PROCESS");
  double tolLL = BlockPROCESS.get_double("tolLL");
  double tolTime = BlockPROCESS.get_double("tolTime");
  //
  auto IsAboutIdentical = [&](SingleEntryMeasurement const &x,
                              SingleEntryMeasurement const &y) -> bool {
    if (x.Satellite != y.Satellite)
      return false;
    if (fabs(x.Time - y.Time) > tolTime)
      return false;
    if (fabs(x.Lon - y.Lon) > tolLL)
      return false;
    if (fabs(x.Lat - y.Lat) > tolLL)
      return false;
    return true;
  };
  //
  auto FindPosition =
      [&](SingleEntryMeasurement const &x,
          std::vector<SingleEntryMeasurement> const &ListY) -> int {
    int eFirst = 0;
    int eLast = ListY.size();
    for (int iter = 0; iter < 10; iter++) {
      int eMid = (eFirst + eLast) / 2;
      eMid = std::max(0, eMid);
      eMid = std::min(eMid, eLast - 1);
      if (eLast - eFirst > 100) {
        if (x.Time < ListY[eMid].Time)
          eLast = eMid;
        else
          eFirst = eMid;
      }
    }
    int posRet = -1;
    double deltaTime = 340;
    double deltaLL = 200;
    for (int pos = eFirst; pos < eLast; pos++) {
      if (IsAboutIdentical(x, ListY[pos])) {
        double eDeltaTime = fabs(x.Time - ListY[pos].Time);
        double eDeltaLL =
            fabs(x.Lon - ListY[pos].Lon) + fabs(x.Lat - ListY[pos].Lat);
        if (eDeltaTime < deltaTime && eDeltaLL < deltaLL) {
          deltaTime = eDeltaTime;
          deltaLL = eDeltaLL;
          posRet = pos;
        }
      }
    }
    return posRet;
  };
  //
  std::vector<std::pair<double, double>> ListPairHwave, ListPairWind;
  std::vector<double> ListTimeMatch;
  int nbUnmatched = 0;
  for (auto &eEntry0 : ListSingleEntry0) {
    int pos = FindPosition(eEntry0, ListSingleEntry1);
    if (pos == -1) {
      nbUnmatched++;
    } else {
      SingleEntryMeasurement eEntry1 = ListSingleEntry1[pos];
      ListPairHwave.push_back({eEntry0.Swh_used, eEntry1.Swh_used});
      ListPairWind.push_back({eEntry0.WindSpeed_used, eEntry1.WindSpeed_used});
      ListTimeMatch.push_back(eEntry0.Time);
    }
  }
  std::cerr << "nbUnmatched = " << nbUnmatched << "\n";
  std::cerr << "nbMatched = " << ListTimeMatch.size() << "\n";
  //
  auto PlotDataSource =
      [&](std::string const &VarName,
          std::vector<std::pair<double, double>> const &ListPairSel,
          std::string const &SatName,
          std::vector<double> const &ListTimeSel) -> void {
    int nbEnt = ListPairSel.size();
    MyVector<double> ListMeas0(nbEnt), ListMeas1(nbEnt);
    for (int i = 0; i < nbEnt; i++) {
      ListMeas0(i) = ListPairSel[i].first;
      ListMeas1(i) = ListPairSel[i].second;
    }
    double minTime = VectorMin(ListTimeSel);
    double maxTime = VectorMax(ListTimeSel);
    std::string strPresMin = DATE_ConvertMjd2mystringPres(minTime);
    std::string strPresMax = DATE_ConvertMjd2mystringPres(maxTime);
    double TheMax = std::max(ListMeas0.maxCoeff(), ListMeas1.maxCoeff());
    double TheMin = std::min(ListMeas0.minCoeff(), ListMeas1.minCoeff());
    //
    MyVector<double> ListX(nbEnt);
    for (int i = 0; i < nbEnt; i++)
      ListX(i) = static_cast<double>(i);
    //
    DrawLinesArr eDrawArr;
    eDrawArr.DoTitle = true;
    eDrawArr.TitleStr =
        "Plotting " + VarName + " from " + strPresMin + " to " + strPresMax;
    std::string fVarName = SatName + "_" + VarName;
    eDrawArr.IsTimeSeries = false;
    eDrawArr.PairComparison = true;
    eDrawArr.DoExplicitLabel = false;
    eDrawArr.VarName = fVarName;
    eDrawArr.TheMax = TheMax;
    eDrawArr.TheMin = TheMin;
    eDrawArr.ListX = ListX;
    eDrawArr.ListListVect = {ListMeas0, ListMeas1};
    eDrawArr.ListName_plot = {eName0, eName1};
    std::string FileName = ePerm.eDir + SatName + "_" + VarName + "_comparison";
    LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
  };
  //
  std::cerr << "Before Data output\n";
  if (BlPROCESS.get_bool("PLOT_ALL_TRACKS")) {
    PlotDataSource("Hwave", ListPairHwave, "allsat", ListTimeMatch);
    PlotDataSource("WindSpeed", ListPairWind, "allsat", ListTimeMatch);
  }
  std::cerr << "After PLOT_ALL_TRACKS output\n";
  if (BlPROCESS.get_bool("DO_STAT")) {
    std::ofstream os(ePerm.eDir + "statistics.txt");
    std::cerr << "|ListPairHwave|=" << ListPairHwave.size() << "\n";
    T_stat eStatWave = ComputeStatistics_stdpair(ListPairHwave);
    std::string optStatStr = BlPLOT.get_string("optStatStr");
    T_statString eStatStringWave =
        ComputeStatisticString_from_Statistics(eStatWave, optStatStr);
    std::cerr << "After statistic computations\n";
    os << " ----------------------------------------------\n";
    os << "               H wave                          \n";
    os << " ----------------------------------------------\n";
    Print_Down_Statistics(os, "hwave statistics", eStatWave);
    os << eStatStringWave.strNature << "\n";
    os << eStatStringWave.str << "\n";
    //
    std::cerr << "|ListPairWind|=" << ListPairWind.size() << "\n";
    T_stat eStatWind = ComputeStatistics_stdpair(ListPairWind);
    T_statString eStatStringWind =
        ComputeStatisticString_from_Statistics(eStatWind, optStatStr);
    os << " ----------------------------------------------\n";
    os << "               wind speed                      \n";
    os << " ----------------------------------------------\n";
    Print_Down_Statistics(os, "wind speed statistics", eStatWind);
    os << eStatStringWind.strNature << "\n";
    os << eStatStringWind.str << "\n";
  }
  std::cerr << "After DO_STAT output\n";
  //
  std::cerr << "Before DO_TXTRAW output\n";
  if (BlPROCESS.get_bool("DO_TXTRAW")) {
    std::ofstream os(ePerm.eDir + "RawData.txt");
    if (BlPROCESS.get_bool("DO_WNDMAG")) {
      os << "Data of wind speed (meas/model)\n";
      for (auto &ePair : ListPairWind) {
        os << "meas=" << ePair.first << " model=" << ePair.second << "\n";
      }
    }
    if (BlPROCESS.get_bool("DO_HS")) {
      os << "Data of significant wave height (meas/model)\n";
      for (auto &ePair : ListPairHwave) {
        os << "meas=" << ePair.first << " model=" << ePair.second << "\n";
      }
    }
  }
  std::cerr << "After DO_TXTRAW output\n";
}

//
// SST code
//

FullNamelist NAMELIST_GetStandardSST_COMPARISON() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  ListStringValues1["BEGTC"] = "unset";
  ListStringValues1["ENDTC"] = "unset";
  ListStringValues1["MODELNAME"] = "unset MODELNAME";
  ListStringValues1["GridFile"] = "unset GridFile";
  ListStringValues1["HisPrefix"] = "ROMS_output_";
  ListStringValues1["SST_files_prefix"] = "unset";
  ListDoubleValues1["PreDawnHour"] = 5.0 / 24.0;
  ListStringValues1["Extension"] = "png";
  ListStringValues1["__NaturePlot"] = "SST_COMPARISON";
  ListStringValues1["PicPrefix"] = "Float_Output_";
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Lines_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  ListIntValues1["NPROC"] = 1;
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues1);
  BlockPROC.setListBoolValues(ListBoolValues1);
  BlockPROC.setListDoubleValues(ListDoubleValues1);
  BlockPROC.setListStringValues(ListStringValues1);
  BlockPROC.setListListStringValues(ListListStringValues1);
  BlockPROC.setListListIntValues(ListListIntValues1);
  ListBlock["PROC"] = BlockPROC;
  // STAT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListDoubleValues2["MaxErr_L4"] = true;
  ListIntValues2["GEOSELECTION"] = 1;
  ListDoubleValues2["MinLON"] = -7;
  ListDoubleValues2["MaxLON"] = 37;
  ListDoubleValues2["MinLAT"] = 30;
  ListDoubleValues2["MaxLAT"] = 46;
  ListListDoubleValues2["LONPOLY"] = {10, 10, 10};
  ListListDoubleValues2["LATPOLY"] = {10, 10, 10};
  ListBoolValues2["DoMinDistCoast"] = false;
  ListDoubleValues2["MinDistCoastKM"] = 60;
  ListStringValues2["LonLatDiscFile"] = "Float_Output_";
  SingleBlock BlockSTAT;
  BlockSTAT.setListIntValues(ListIntValues2);
  BlockSTAT.setListBoolValues(ListBoolValues2);
  BlockSTAT.setListDoubleValues(ListDoubleValues2);
  BlockSTAT.setListStringValues(ListStringValues2);
  BlockSTAT.setListListStringValues(ListListStringValues2);
  BlockSTAT.setListListDoubleValues(ListListDoubleValues2);
  ListBlock["STAT"] = BlockSTAT;
  // PLOT
  std::map<std::string, int> ListIntValues4;
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, double> ListDoubleValues4;
  std::map<std::string, std::string> ListStringValues4;
  std::map<std::string, std::vector<std::string>> ListListStringValues4;
  std::map<std::string, std::vector<double>> ListListDoubleValues4;
  ListIntValues4["nbBlock"] = 1;
  ListIntValues4["nbLevelSpa"] = 20;
  ListIntValues4["nbLabelStride"] = 10;
  ListBoolValues4["DoColorBar"] = true;
  ListBoolValues4["DoTitle"] = true;
  ListBoolValues4["FillLand"] = true;
  ListBoolValues4["PrintMMA"] = false;
  ListBoolValues4["DrawRiver"] = false;
  ListBoolValues4["DrawContourBathy"] = false;
  ListBoolValues4["DrawAnnotation"] = false;
  ListBoolValues4["DoTitleString"] = true;
  ListBoolValues4["cnSmoothingOn"] = true;
  ListBoolValues4["UseNativeGrid"] = true;
  ListBoolValues4["UseDynamicRangeInScatter"] = false;
  ListBoolValues4["VariableRange"] = false;
  ListBoolValues4["DoMain"] = true;
  ListListDoubleValues4["ListFrameMinLon"] = {};
  ListListDoubleValues4["ListFrameMinLat"] = {};
  ListListDoubleValues4["ListFrameMaxLon"] = {};
  ListListDoubleValues4["ListFrameMaxLat"] = {};
  ListDoubleValues4["vcRefLengthF"] = 0.02;
  ListListStringValues2["ListSubstitution"] = {};
  ListDoubleValues4["AnnotationLon"] = -400;
  ListDoubleValues4["AnnotationLat"] = -400;
  ListStringValues4["AnnotationText"] = "unset_text";
  ListStringValues4["FileDirectNCLins"] = "irrelevant";
  ListStringValues4["GridResolution"] = "HighRes";
  ListStringValues4["ColorMap"] = "BlAqGrYeOrReVi200";
  ListStringValues4["cnFillMode"] = "RasterFill";
  ListBoolValues4["cnFillOn"] = true;
  ListBoolValues4["cnLinesOn"] = false;
  ListBoolValues4["cnLineLabelsOn"] = false;
  ListStringValues4["LandPortr"] = "Landscape";
  ListStringValues4["optStatStr"] = "double";
  ListBoolValues4["DoPlotScatter"] = true;
  ListBoolValues4["DoPlotDiff"] = true;
  ListBoolValues4["DoPlotTimeAverage"] = true;
  SingleBlock BlockPLOT;
  BlockPLOT.setListIntValues(ListIntValues4);
  BlockPLOT.setListBoolValues(ListBoolValues4);
  BlockPLOT.setListDoubleValues(ListDoubleValues4);
  BlockPLOT.setListStringValues(ListStringValues4);
  BlockPLOT.setListListStringValues(ListListStringValues4);
  BlockPLOT.setListListDoubleValues(ListListDoubleValues4);
  ListBlock["PLOT"] = BlockPLOT;
  // Final part
  return FullNamelist(ListBlock);
}

int GetSmallestIndex(MyVector<double> const &V, double eVal) {
  int pos = -1;
  double MinDist = 100000000;
  for (int idx = 0; idx < static_cast<int>(V.size()); idx++) {
    double eDist = fabs(V(idx) - eVal);
    if (eDist < MinDist) {
      MinDist = eDist;
      pos = idx;
    }
  }
  return pos;
}

void RAW_SCATTER_SST(std::vector<double> const &V_meas,
                     std::vector<double> const &V_model,
                     NCLcaller<GeneralType> &eCall,
                     PermanentInfoDrawing const &ePerm) {
  std::cerr << "Running RAW_SCATTER_SST\n";
  SingleBlock eBlPLOT = ePerm.eFull.get_block("PLOT");
  bool UseDynamicRangeInScatter =
    eBlPLOT.get_bool("UseDynamicRangeInScatter");
  //
  DrawScatterArr eDrw;
  int siz = V_meas.size();
  MyVector<double> eVectA(siz);
  MyVector<double> eVectB(siz);
  for (int i = 0; i < siz; i++) {
    eVectA[i] = V_meas[i];
    eVectB[i] = V_model[i];
  }
  std::vector<double> data_rangeA(2);
  std::vector<double> data_rangeB(2);
  //
  std::string VarType = "SST";
  std::string eUnit = "deg C";
  double TheMin = 10;
  double TheMax = 25;
  if (UseDynamicRangeInScatter && siz > 0) {
    double maxA = eVectA.maxCoeff();
    double maxB = eVectB.maxCoeff();
    TheMax = std::max(maxA, maxB);
    double minA = eVectA.maxCoeff();
    double minB = eVectB.maxCoeff();
    TheMin = std::max(minA, minB);
  }
  data_rangeA[0] = TheMin;
  data_rangeA[1] = TheMax;
  data_rangeB[0] = TheMin;
  data_rangeB[1] = TheMax;
  eDrw.VarNameAB_file = "Scatter_iGrid_" + VarType;
  //
  eDrw.DoTitle = false;
  eDrw.AddStatMeasModel = true;
  eDrw.NameA_plot = "Data (" + eUnit + ")";
  eDrw.NameB_plot = "Model (" + eUnit + ")";
  eDrw.data_rangeA = data_rangeA;
  eDrw.data_rangeB = data_rangeB;
  eDrw.eVectA = eVectA;
  eDrw.eVectB = eVectB;
  eDrw.aSize = 100;
  eDrw.bSize = 100;
  PLOT_SCATTER(eDrw, eCall, ePerm);
}

void Process_sst_Comparison_Request(FullNamelist const &eFull) {
  SingleBlock eBlPROC = eFull.get_block("PROC");
  SingleBlock eBlSTAT = eFull.get_block("STAT");
  SingleBlock eBlPLOT = eFull.get_block("PLOT");
  //
  // Now basic definitions
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm
  //
  // Reading the model
  //
  std::string ModelName = eBlPROC.get_string("MODELNAME");
  std::string GridFile = eBlPROC.get_string("GridFile");
  std::string HisPrefix = eBlPROC.get_string("HisPrefix");
  TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  //
  // Reading the list of files and times.
  //
  std::string SST_files_prefix =
    eBlPROC.get_string("SST_files_prefix");
  std::vector<std::string> ListFile =
      FILE_DirectoryMatchingPrefixExtension(SST_files_prefix, "nc");
  struct SingEnt {
    size_t iFile;
    size_t iTime;
    double eTime;
  };
  std::vector<SingEnt> ListSingEnt;
  for (size_t iFile = 0; iFile < ListFile.size(); iFile++) {
    std::string eFile = ListFile[iFile];
    std::cerr << "iFile=" << iFile << " eFile=" << eFile << "\n";
    std::vector<double> LTime = NC_ReadTimeFromFile(eFile, "time");
    for (size_t iTime = 0; iTime < LTime.size(); iTime++)
      ListSingEnt.push_back({iFile, iTime, LTime[iTime]});
  }
  std::cerr << "|ListSingEnt|=" << ListSingEnt.size() << "\n";
  //
  // Determining the beginning and ending of time for comparison
  //
  std::string strBEGTC = eBlPROC.get_string("BEGTC");
  std::string strENDTC = eBlPROC.get_string("ENDTC");
  double BeginTime = 0, EndTime = 0;
  if (strBEGTC == "earliest") {
    BeginTime = MinimumTimeHistoryArray(TotalArr.eArr);
  } else {
    BeginTime = CT2MJD(strBEGTC);
  }
  if (strENDTC == "latest") {
    EndTime = MaximumTimeHistoryArray(TotalArr.eArr);
  } else {
    EndTime = CT2MJD(strENDTC);
  }
  double PreDawnHour = eBlPROC.get_double("PreDawnHour");
  if (!IsZeroHour(BeginTime)) {
    std::string strPresBegin = DATE_ConvertMjd2mystringPres(BeginTime);
    std::cerr << "The initial date should be a zero hour\n";
    std::cerr << "That is hour=0 , min=0 , sec=0\n";
    std::cerr << "On the other hand we have BeginTime=" << strPresBegin << "\n";
    throw TerminalException{1};
  }
  std::cerr << "We have BeginTime, EndTime, PreDawnTime\n";
  //
  // Reading the SST grid and computing interpolation arrays
  //
  std::string eFileSST = ListFile[0];
  MyVector<double> VectLat = NC_Read1Dvariable(eFileSST, "lat");
  MyVector<double> VectLon = NC_Read1Dvariable(eFileSST, "lon");
  size_t n_lat = VectLat.size();
  size_t n_lon = VectLon.size();
  MyMatrix<double> const &LON = TotalArr.GrdArr.GrdArrRho.LON;
  MyMatrix<double> const &LAT = TotalArr.GrdArr.GrdArrRho.LAT;
  MyMatrix<uint8_t> const &MSK = TotalArr.GrdArr.GrdArrRho.MSK;
  size_t eta_rho = LON.rows();
  size_t xi_rho = LON.cols();
  MyMatrix<int> MatIdxLat(eta_rho, xi_rho);
  MyMatrix<int> MatIdxLon(eta_rho, xi_rho);
  for (size_t iEta = 0; iEta < eta_rho; iEta++)
    for (size_t iXi = 0; iXi < xi_rho; iXi++) {
      double eLon = LON(iEta, iXi);
      double eLat = LAT(iEta, iXi);
      //
      MatIdxLat(iEta, iXi) = GetSmallestIndex(VectLat, eLat);
      MatIdxLon(iEta, iXi) = GetSmallestIndex(VectLon, eLon);
    }
  std::cerr << "We have MatIdxLat, MatIdxLon\n";
  //
  // Now computing the masks at a given time
  //
  int GEOSELECTION = eBlSTAT.get_int("GEOSELECTION");
  double MinLON = eBlSTAT.get_double("MinLON");
  double MaxLON = eBlSTAT.get_double("MaxLON");
  double MinLAT = eBlSTAT.get_double("MinLAT");
  double MaxLAT = eBlSTAT.get_double("MaxLAT");
  std::vector<double> LonPoly = eBlSTAT.get_list_double("LONPOLY");
  std::vector<double> LatPoly = eBlSTAT.get_list_double("LATPOLY");
  MyMatrix<uint8_t> MSK_geog(eta_rho, xi_rho);
  for (size_t iEta = 0; iEta < eta_rho; iEta++)
    for (size_t iXi = 0; iXi < xi_rho; iXi++) {
      double eLon = LON(iEta, iXi);
      double eLat = LAT(iEta, iXi);
      //
      uint8_t eMask_Geog = 1;
      if (GEOSELECTION == 1) {
        if (!(eLon >= MinLON && eLon <= MaxLON && eLat >= MinLAT &&
              eLat <= MaxLAT))
          eMask_Geog = 0;
      }
      if (GEOSELECTION == 2) {
        bool IsInside = IsPointInside(eLon, eLat, LonPoly, LatPoly);
        if (!IsInside)
          eMask_Geog = 0;
      }
      MSK_geog(iEta, iXi) = eMask_Geog;
    }
  std::cerr << "We have MSK_geog\n";
  //
  // Computing the distance mask
  //
  MyMatrix<uint8_t> MSK_dist = MSK;
  bool DoMinDistCoast = eBlSTAT.get_bool("DoMinDistCoast");
  if (DoMinDistCoast) {
    std::string eFileCoast = eBlSTAT.get_string("LonLatDiscFile");
    double MinDistCoastKM = eBlSTAT.get_double("MinDistCoastKM");
    std::cerr << "eFileCoast=" << eFileCoast << "\n";
    std::vector<PairLL> ListPtCoast = ReadLonLatDiscFile(eFileCoast);
    std::vector<PairLL> ListPt;
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++) {
        if (MSK(iEta, iXi) == 1) {
          double eLon = LON(iEta, iXi);
          double eLat = LAT(iEta, iXi);
          PairLL ePt{eLon, eLat};
          ListPt.push_back(ePt);
        }
      }
    std::vector<double> ListMinDist =
        GetListMinimalDistances(ListPtCoast, ListPt);
    int idx = 0;
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++) {
        if (MSK(iEta, iXi) == 1) {
          if (ListMinDist[idx] < MinDistCoastKM)
            MSK_dist(iEta, iXi) = 0;
          idx++;
        }
      }
  }
  std::cerr << "We have MSK_dist\n";
  //
  // Now processing the comparison
  //
  auto GetEntry = [&](double const &eTime) -> std::pair<size_t, size_t> {
    for (auto &eSingEnt : ListSingEnt) {
      if (fabs(eSingEnt.eTime - eTime) < 0.0001) {
        return {eSingEnt.iFile, eSingEnt.iTime};
      }
    }
    return {-1, -1};
  };
  auto ReadSST_entry = [&](std::pair<size_t, size_t> const &ePair,
                           std::string const &eVar) -> MyMatrix<double> {
    std::string eFile = ListFile[ePair.first];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar data = dataFile.getVar(eVar);
    //
    std::vector<size_t> start{ePair.second, 0, 0};
    std::vector<size_t> count{1, n_lat, n_lon};
    MyVector<double> V = NC_ReadVariable_data_start_count(data, start, count);
    MyMatrix<double> M(n_lat, n_lon);
    size_t idx = 0;
    for (size_t i_lat = 0; i_lat < n_lat; i_lat++)
      for (size_t i_lon = 0; i_lon < n_lon; i_lon++) {
        M(i_lat, i_lon) = V(idx);
        idx++;
      }
    return M;
  };
  auto ReadSST_Pair = [&](std::pair<size_t, size_t> const &ePair)
      -> std::pair<MyMatrix<double>, MyMatrix<double>> {
    MyMatrix<double> Mat_SST = ReadSST_entry(ePair, "analysed_sst");
    MyMatrix<double> Mat_ERR = ReadSST_entry(ePair, "analysis_error");
    for (size_t i_lat = 0; i_lat < n_lat; i_lat++)
      for (size_t i_lon = 0; i_lon < n_lon; i_lon++) {
        if (Mat_SST(i_lat, i_lon) < -50)
          Mat_ERR(i_lat, i_lon) = 300;
      }
    return {std::move(Mat_SST), std::move(Mat_ERR)};
  };
  double MaxErr_L4 = eBlSTAT.get_double("MaxErr_L4");
  bool DoPlotScatter = eBlPLOT.get_bool("DoPlotScatter");
  bool DoPlotDiff = eBlPLOT.get_bool("DoPlotDiff");
  bool DoPlotTimeAverage = eBlPLOT.get_bool("DoPlotTimeAverage");
  std::cerr << "MaxErr_L4=" << MaxErr_L4 << "\n";
  std::string PicPrefix = eBlPROC.get_string("PicPrefix");
  std::string FileStatDaily = PicPrefix + "Statistics_Daily.txt";
  std::ofstream os(FileStatDaily);
  std::vector<double> V_meas_total;
  std::vector<double> V_model_total;
  std::vector<QuadDrawInfo> ListQuad =
      GetListQuadArray(eBlPLOT, TotalArr.GrdArr);
  std::cerr << "|ListQuad|=" << ListQuad.size() << "\n";
  int iTime = 0;
  GridArray GrdArr_Plot = TotalArr.GrdArr;
  MyMatrix<uint8_t> &MSK_plot = GrdArr_Plot.GrdArrRho.MSK;
  //
  MyMatrix<double> SumDiff = ZeroMatrix<double>(eta_rho, xi_rho);
  MyMatrix<int> SumAtt = ZeroMatrix<int>(eta_rho, xi_rho);
  for (double eTime = BeginTime; eTime <= EndTime; eTime += 1.0) {
    std::string strPres = DATE_ConvertMjd2mystringPres(eTime);
    std::cerr << "iTime=" << iTime << " date=" << strPres << "\n";
    double eTimeCall = eTime + PreDawnHour;
    RecVar eRecVar =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "TempSurf", eTimeCall);
    //
    std::pair<size_t, size_t> ePair = GetEntry(eTime);
    auto ePairSST_ERR = ReadSST_Pair(ePair);
    MyMatrix<double> const &Mat_SST = ePairSST_ERR.first;
    MyMatrix<double> const &Mat_ERR = ePairSST_ERR.second;
    //
    std::vector<double> V_meas;
    std::vector<double> V_model;
    double CorrKelvin_Celsius = 273.15;
    MyMatrix<double> M_Diff(eta_rho, xi_rho);
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++) {
        //
        // The model and measurement
        //
        double eValModel = eRecVar.F(iEta, iXi);
        int i_lat = MatIdxLat(iEta, iXi);
        int i_lon = MatIdxLon(iEta, iXi);
        double eValMeas = Mat_SST(i_lat, i_lon) - CorrKelvin_Celsius;
        double MeasERR = Mat_ERR(i_lat, i_lon);
        //
        // limitations
        //
        // 1 : mask elimination
        bool IsCorrect = true;
        if (MSK(iEta, iXi) == 0) {
          IsCorrect = false;
        }
        // 2 : elimination via the error
        if (MeasERR > MaxErr_L4) {
          IsCorrect = false;
        }
        // 3 : mask by the distance
        if (MSK_dist(iEta, iXi) == 0) {
          IsCorrect = false;
        }
        // 4: mask by the geographical selection
        if (MSK_geog(iEta, iXi) == 0) {
          IsCorrect = false;
        }
        //
        // Now processing the data
        //
        uint8_t eMSK = 0;
        if (IsCorrect) {
          //          std::cerr << "iEta=" << iEta << " iXi=" << iXi << " meas="
          //          << eValMeas << " err=" << MeasERR << "\n";
          V_meas.push_back(eValMeas);
          V_model.push_back(eValModel);
          V_meas_total.push_back(eValMeas);
          V_model_total.push_back(eValModel);
          M_Diff(iEta, iXi) = eValMeas - eValModel;
          double diff = fabs(eValMeas - eValModel);
          SumDiff(iEta, iXi) += diff;
          SumAtt(iEta, iXi) += 1;
          eMSK = 1;
        }
        MSK_plot(iEta, iXi) = eMSK;
      }
    // Now computing the stats
    T_stat estat = ComputeStatistics_vector(V_meas, V_model);
    T_statString estatstr =
        ComputeStatisticString_from_Statistics(estat, "4dot2f");
    std::cerr << "    date=" << strPres << " stat=" << estatstr.str << "\n";
    std::cerr << "    nbMeas=" << estat.nbMeas << " MeanMeas=" << estat.MeanMeas
              << " MeanModel=" << estat.MeanModel << "\n";
    std::cerr << "    MinMeas=" << estat.MinMeas << " MaxMeas=" << estat.MaxMeas
              << "\n";
    //
    // Plotting the difference
    //
    if (DoPlotDiff) {
      eRecVar.F = M_Diff;
      eRecVar.RecS.minval = -1;
      eRecVar.RecS.maxval = 1;
      for (auto &eQuad : ListQuad) {
        std::string FileName = PicPrefix + "SST_" + eQuad.eFrameName + "_" +
                               StringNumber(iTime, 4);
        eRecVar.RecS.strAll =
            "SST_" + eQuad.eFrameName + "_" + StringNumber(iTime, 4);
        DrawArr eDrw = ePerm.eDrawArr;
        eDrw.eQuadFrame = eQuad.eQuad;
        eDrw.DoTitle = true;
        eDrw.TitleStr = "Difference between SST and model at " + strPres;
        PLOT_PCOLOR(FileName, GrdArr_Plot, eDrw, eRecVar, eCall, ePerm);
      }
    }
    iTime++;
  }
  //
  // Putting the absolute value together
  //
  if (DoPlotTimeAverage) {
    //
    // First the error average at the nodes.
    //
    RecVar eRecVar =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Bathymetry", BeginTime);
    MyMatrix<uint8_t> MSK_final(eta_rho, xi_rho);
    MyMatrix<double> F(eta_rho, xi_rho);
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++) {
        int eval = SumAtt(iEta, iXi);
        double eval_f;
        uint8_t eval_i;
        if (eval > 0) {
          eval_f = SumDiff(iEta, iXi) / eval;
          eval_i = 1;
        } else {
          eval_f = 0;
          eval_i = 0;
        }
        MSK_final(iEta, iXi) = eval_i;
        F(iEta, iXi) = eval_f;
      }
    GrdArr_Plot.GrdArrRho.MSK = MSK_final;
    eRecVar.F = F;
    eRecVar.RecS.minval = 0;
    eRecVar.RecS.maxval = F.maxCoeff();
    for (auto &eQuad : ListQuad) {
      std::string FileName = PicPrefix + "SSTerr_glob_" + eQuad.eFrameName;
      eRecVar.RecS.strAll = "SSTerr_glob_" + eQuad.eFrameName;
      DrawArr eDrw = ePerm.eDrawArr;
      eDrw.eQuadFrame = eQuad.eQuad;
      eDrw.DoTitle = true;
      eDrw.TitleStr = "Average absolute error for the model";
      PLOT_PCOLOR(FileName, GrdArr_Plot, eDrw, eRecVar, eCall, ePerm);
    }
    //
    // The number of attained points.
    //
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++) {
        int eval = SumAtt(iEta, iXi);
        double eval_f;
        uint8_t eval_i;
        if (eval > 0) {
          eval_f = eval;
          eval_i = 1;
        } else {
          eval_f = 0;
          eval_i = 0;
        }
        MSK_final(iEta, iXi) = eval_i;
        F(iEta, iXi) = eval_f;
      }
    GrdArr_Plot.GrdArrRho.MSK = MSK_final;
    eRecVar.F = F;
    eRecVar.RecS.minval = 0;
    eRecVar.RecS.maxval = F.maxCoeff();
    for (auto &eQuad : ListQuad) {
      std::string FileName = PicPrefix + "SSTerr_glob_" + eQuad.eFrameName;
      eRecVar.RecS.strAll = "AttainedNb_glob_" + eQuad.eFrameName;
      DrawArr eDrw = ePerm.eDrawArr;
      eDrw.eQuadFrame = eQuad.eQuad;
      eDrw.DoTitle = true;
      eDrw.TitleStr = "Average absolute error for the model";
      PLOT_PCOLOR(FileName, GrdArr_Plot, eDrw, eRecVar, eCall, ePerm);
    }
  }
  //
  std::string FileStat = PicPrefix + "global_stat.txt";
  std::ofstream os_stat(FileStat);
  T_stat estat = ComputeStatistics_vector(V_meas_total, V_model_total);
  Print_Down_Statistics(os_stat, "global", estat);
  //
  // Scatter plot of measurement and model
  //
  if (DoPlotScatter) {
    std::cerr << "|V_meas_total|=" << V_meas_total.size() << "\n";
    RAW_SCATTER_SST(V_meas_total, V_model_total, eCall, ePerm);
  }
}

//
// Comparison with CTD
//

FullNamelist NAMELIST_GetStandardCTD_COMPARISON() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  ListStringValues1["MODELNAME"] = "unset MODELNAME";
  ListStringValues1["GridFile"] = "unset GridFile";
  ListStringValues1["HisPrefix"] = "ROMS_output_";
  ListStringValues1["File_CTD"] = "unset";
  ListStringValues1["FileOut"] = "unset.out";
  ListStringValues1["FileStat"] = "unset.stat";
  ListStringValues1["FileByStation"] = "unset.station";
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues1);
  BlockPROC.setListBoolValues(ListBoolValues1);
  BlockPROC.setListDoubleValues(ListDoubleValues1);
  BlockPROC.setListStringValues(ListStringValues1);
  BlockPROC.setListListStringValues(ListListStringValues1);
  BlockPROC.setListListIntValues(ListListIntValues1);
  ListBlock["PROC"] = BlockPROC;
  // Final part
  return FullNamelist(ListBlock);
}

void Process_ctd_Comparison_Request(FullNamelist const &eFull) {
  SingleBlock const& eBlPROC = eFull.get_block("PROC");
  //
  // Reading the model
  //
  std::string ModelName = eBlPROC.get_string("MODELNAME");
  std::string GridFile = eBlPROC.get_string("GridFile");
  std::string HisPrefix = eBlPROC.get_string("HisPrefix");
  TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  //
  // Reading the list of files and times.
  //
  std::string File_CTD = eBlPROC.get_string("File_CTD");
  std::vector<std::string> ListLines_CTD = ReadFullFile(File_CTD);
  std::vector<std::string> ListNames;
  std::vector<double> ListLon, ListLat, ListDep, ListDate, ListTempMeas,
      ListSaltMeas;
  std::set<std::string> SetNames;
  for (auto &eLine : ListLines_CTD) {
    std::vector<std::string> LStr = STRING_Split(eLine, ";");
    std::string eName = LStr[0];
    double eLon = std::stod(LStr[1]);
    double eLat = std::stod(LStr[2]);
    double eDep = std::stod(LStr[3]);
    double eDate = CT2MJD(LStr[4]);
    double eTemp = std::stod(LStr[5]);
    double eSalt = std::stod(LStr[6]);
    //
    ListNames.push_back(eName);
    ListLon.push_back(eLon);
    ListLat.push_back(eLat);
    ListDep.push_back(eDep);
    ListDate.push_back(eDate);
    ListTempMeas.push_back(eTemp);
    ListSaltMeas.push_back(eSalt);
    SetNames.insert(eName);
  }
  //
  // Processing the output
  //
  int nbLine = ListLines_CTD.size();
  std::vector<double> ListTempModel(nbLine);
  std::vector<double> ListSaltModel(nbLine);
  std::vector<int> ListStatus(nbLine);
  for (int iLine = 0; iLine < nbLine; iLine++) {
    double eDate = ListDate[iLine];
    double eLon = ListLon[iLine];
    double eLat = ListLat[iLine];
    double eDep = ListDep[iLine];
    double eDep_input = -eDep;
    std::cerr << "iLine=" << iLine << " / " << nbLine << " eLon=" << eLon
              << " eLat=" << eLat << " eDep=" << eDep << "\n";
    //
    RecVar eRecVar_T =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Temp", eDate);
    RecVar eRecVar_S =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "Salt", eDate);
    auto LDim = eRecVar_T.Tens3.dimensions();
    int dim0 = LDim[0];
    int dim1 = LDim[1];
    int dim2 = LDim[2];
    Eigen::Tensor<double, 3> Tens3(dim0, dim1, dim2);
    for (int i0 = 0; i0 < dim0; i0++)
      for (int i1 = 0; i1 < dim1; i1++)
        for (int i2 = 0; i2 < dim2; i2++)
          Tens3(i0, i1, i2) = 1;
    //
    VerticalLevelInfo eVert{2, eDep_input, eDep_input, "irrelevant 1", "irrelevant 2"};
    MyMatrix<double> zeta =
        ModelSpecificVarSpecificTime_Kernel(TotalArr, "ZetaOcean", eDate).F;
    MyMatrix<double> F_horizTemp = ThreeDimensional_to_TwoDimensional(
        eRecVar_T.Tens3, zeta, TotalArr, eVert, eDate);
    MyMatrix<double> F_horizSalt = ThreeDimensional_to_TwoDimensional(
        eRecVar_S.Tens3, zeta, TotalArr, eVert, eDate);
    MyMatrix<double> F_horizOne =
        ThreeDimensional_to_TwoDimensional(Tens3, zeta, TotalArr, eVert, eDate);
    MyMatrix<double> ListXY(2, 1);
    ListXY(0, 0) = eLon;
    ListXY(1, 0) = eLat;
    SingleRecInterp eRec =
        General_FindInterpolationWeight(TotalArr.GrdArr, ListXY, false)[0];
    double eTempModel = 0;
    double eSaltModel = 0;
    int eStatus = 0;
    if (eRec.status) {
      double eStatus_d = 0;
      for (auto &ePart : eRec.LPart) {
        eTempModel += ePart.eCoeff * F_horizTemp(ePart.eEta, ePart.eXi);
        eSaltModel += ePart.eCoeff * F_horizSalt(ePart.eEta, ePart.eXi);
        eStatus_d += ePart.eCoeff * F_horizOne(ePart.eEta, ePart.eXi);
      }
      if (eStatus_d >= 0.99999)
        eStatus = 1;
    }
    ListTempModel[iLine] = eTempModel;
    ListSaltModel[iLine] = eSaltModel;
    ListStatus[iLine] = eStatus;
  }
  //
  // Printing the results.
  //
  std::string FileOut = eBlPROC.get_string("FileOut");
  std::ofstream os(FileOut);
  for (int iLine = 0; iLine < nbLine; iLine++) {
    std::string eName = ListNames[iLine];
    double eDate = ListDate[iLine];
    double eLon = ListLon[iLine];
    double eLat = ListLat[iLine];
    double eDep = ListDep[iLine];
    //
    double eTempMeas = ListTempMeas[iLine];
    double eSaltMeas = ListSaltMeas[iLine];
    double eTempModel = ListTempModel[iLine];
    double eSaltModel = ListSaltModel[iLine];
    int eStatus = ListStatus[iLine];
    std::string strPres = DATE_ConvertMjd2mystringPres(eDate);
    os << strPres << " " << eLon << "/" << eLat << "/" << eDep << "\n";
    os << "       Temp(Meas/Model)=" << eTempMeas << " / " << eTempModel
       << "   Salt(Meas/Model)=" << eSaltMeas << " / " << eSaltModel
       << " status=" << eStatus << "\n";
  }
  std::cerr << "FileOut=" << FileOut << "\n";
  //
  // Keeping the status = 1 entries
  //
  std::vector<std::string> SEL_ListNames;
  std::vector<double> SEL_ListLon, SEL_ListLat, SEL_ListDep, SEL_ListDate,
      SEL_ListTempMeas, SEL_ListSaltMeas, SEL_ListTempModel, SEL_ListSaltModel;
  for (int iLine = 0; iLine < nbLine; iLine++) {
    std::string eName = ListNames[iLine];
    double eDate = ListDate[iLine];
    double eLon = ListLon[iLine];
    double eLat = ListLat[iLine];
    double eDep = ListDep[iLine];
    //
    double eTempMeas = ListTempMeas[iLine];
    double eSaltMeas = ListSaltMeas[iLine];
    double eTempModel = ListTempModel[iLine];
    double eSaltModel = ListSaltModel[iLine];
    if (ListStatus[iLine] == 1) {
      SEL_ListNames.push_back(eName);
      SEL_ListLon.push_back(eLon);
      SEL_ListLat.push_back(eLat);
      SEL_ListDep.push_back(eDep);
      SEL_ListDate.push_back(eDate);
      //
      SEL_ListTempMeas.push_back(eTempMeas);
      SEL_ListSaltMeas.push_back(eSaltMeas);
      SEL_ListTempModel.push_back(eTempModel);
      SEL_ListSaltModel.push_back(eSaltModel);
    }
  }
  int SEL_nbLine = SEL_ListTempMeas.size();
  //
  // Global statistics
  //
  std::string FileStat = eBlPROC.get_string("FileStat");
  std::ofstream os_stat(FileStat);
  T_stat statTemp =
      ComputeStatistics_vector(SEL_ListTempMeas, SEL_ListTempModel);
  T_stat statSalt =
      ComputeStatistics_vector(SEL_ListSaltMeas, SEL_ListSaltModel);
  os_stat << "GLOBAL:\n";
  Print_Down_Statistics(os_stat, "temp", statTemp);
  Print_Down_Statistics(os_stat, "salt", statSalt);
  os_stat << "\n\n";
  //
  // Statistics by station
  //
  std::string FileByStation = eBlPROC.get_string("FileByStation");
  std::ofstream os_bystation(FileByStation);
  os_bystation << " Temp (ME, AE, RMSE)     Salt (ME, AE, RMSE)\n";
  auto DoubleTo5dot2f = [&](double const &x) -> std::string {
    std::stringstream s;
    s << std::setprecision(3);
    s << x;
    std::string s_ret(s.str());
    return s_ret;
  };
  auto fct = [&](double const &x) -> std::string { return DoubleTo5dot2f(x); };
  for (auto &eName : SetNames) {
    os_stat << "Specific to Station=" << eName << "\n";
    std::vector<double> ListTempModel_sel, ListTempMeas_sel;
    std::vector<double> ListSaltModel_sel, ListSaltMeas_sel;
    for (int iLine = 0; iLine < SEL_nbLine; iLine++) {
      if (SEL_ListNames[iLine] == eName) {
        ListTempModel_sel.push_back(SEL_ListTempModel[iLine]);
        ListTempMeas_sel.push_back(SEL_ListTempMeas[iLine]);
        ListSaltModel_sel.push_back(SEL_ListSaltModel[iLine]);
        ListSaltMeas_sel.push_back(SEL_ListSaltMeas[iLine]);
      }
    }
    T_stat statTemp_sel =
        ComputeStatistics_vector(ListTempMeas_sel, ListTempModel_sel);
    T_stat statSalt_sel =
        ComputeStatistics_vector(ListSaltMeas_sel, ListSaltModel_sel);

    Print_Down_Statistics(os_stat, "temp", statTemp_sel);
    Print_Down_Statistics(os_stat, "salt", statSalt_sel);
    os_stat << "\n";
    // By station
    os_bystation << eName << " " << fct(statTemp_sel.MeanError) << " "
                 << fct(statTemp_sel.AbsoluteError) << " "
                 << fct(statTemp_sel.RMSE) << " : "
                 << fct(statSalt_sel.MeanError) << " "
                 << fct(statSalt_sel.AbsoluteError) << " "
                 << fct(statSalt_sel.RMSE) << "\n";
  }
  std::cerr << "FileStat=" << FileStat << "\n";
}

void RadsAscToNetcdf(std::string const& PrefixI, std::string const& PrefixO, int const& method) {
  std::vector<SingleEntryMeasurement> ListEnt;
  std::map<std::string, std::string> MapRadsName_Name;
  MapRadsName_Name["SNTNL-3A"] = "SENTINEL_3A";
  MapRadsName_Name["SNTNL-3B"] = "SENTINEL_3B";
  MapRadsName_Name["SARAL"] = "SARAL";
  MapRadsName_Name["JASON-3"] = "JASON3";
  MapRadsName_Name["CRYOSAT2"] = "CRYOSAT";
  std::vector<std::string> ListMonth = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  std::map<std::string, int> MapMonth;
  for (int iMon=1; iMon<=12; iMon++) {
    MapMonth[ListMonth[iMon-1]] = iMon;
  }
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(PrefixI);
  double minTime = std::numeric_limits<double>::max();
  double maxTime = 0;
  for (auto & eFile : ListFile) {
    std::string eFileF = PrefixI + eFile;
    std::vector<std::string> ListLine = ReadFullFile(eFileF);
    std::string LineSatellite = ListLine[1];
    auto get_satellite=[&]() -> int {
      for (auto & eLine : ListLine) {
        std::vector<std::string> LStrSatellite = STRING_Split(eLine, "Satellite = ");
        if (LStrSatellite.size() == 2) {
          std::string SatelliteNameA = DropSpace(LStrSatellite[1]);
          std::cerr << "SatelliteNameA=\"" << SatelliteNameA << "\"\n";
          std::string SatelliteNameB = MapRadsName_Name.at(SatelliteNameA);
          std::cerr << "SatelliteNameB=\"" << SatelliteNameB << "\"\n";
          int pos = 0;
          for (auto & SatelliteName : GetAllNamesOfSatelliteAltimeter()) {
            pos++;
            if (SatelliteName == SatelliteNameB) {
              return pos;
            }
          }
          std::cerr << "Failed to find the Satellitename\n";
          std::cerr << "SatelliteNameA=" << SatelliteNameA << "\n";
          std::cerr << "SatelliteNameB=" << SatelliteNameB << "\n";
          throw TerminalException{1};
        }
      }
      std::cerr << "Failed to find the satellite\n";
      throw TerminalException{1};
    };
    auto get_equator_time=[&]() -> double {
      for (auto & eLine : ListLine) {
        std::vector<std::string> LStrEquTime = STRING_Split(eLine, "Equ_time ");
        if (LStrEquTime.size() == 2) {
          std::string strA = LStrEquTime[1];

          std::vector<std::string> LStrA = STRING_Split(strA, "(");
          std::string strB = LStrA[1];

          std::vector<std::string> LStrB = STRING_Split(strB, ")");
          std::string strC = LStrB[0];

          int iDay = ParseScalar<int>(strC.substr(0,2));
          int iMonth = MapMonth.at(strC.substr(3,3));
          int iYear = ParseScalar<int>(strC.substr(7,4));
          int iHour = ParseScalar<int>(strC.substr(12,2));
          int iMin = ParseScalar<int>(strC.substr(15,2));
          int iSec = ParseScalar<int>(strC.substr(18,2));
          std::cerr << "iYear=" << iYear << " iMonth=" << iMonth << " iDay=" << iDay << " iHour=" << iHour << " iMin=" << iMin << " iSec=" << iSec << "\n";
          std::vector<int> eDate{iYear, iMonth, iDay, iHour, iMin, iSec};
          double date = DATE_ConvertSix2mjd(eDate);
          return date;
        }
      }
      std::cerr << "Failed to find the equator time\n";
      throw TerminalException{1};
    };
    auto get_columns=[&]() -> std::pair<std::map<std::string, int>, int> {
      std::map<std::string, int> map_column;
      bool HasCOL = false;
      for (int iLine=0; iLine<ListLine.size(); iLine++) {
        std::vector<std::string> LStr = STRING_Split(ListLine[iLine], "Col ");
        if (LStr.size() == 2) {
          HasCOL = true;
          std::string strA = LStr[1];
          std::vector<std::string> LStrA = STRING_Split(strA, "    = ");
          if (LStrA.size() != 2) {
            std::cerr << "LStrA should have length 2\n";
            throw TerminalException{1};
          }
          int iCol = ParseScalar<int>(LStrA[0]);
          std::string desc = LStrA[1];
          map_column[desc] = iCol - 1;
        } else {
          if (HasCOL) {
            return {map_column, iLine};
          }
        }
      }
      std::cerr << "Failed to find the columns\n";
      throw TerminalException{1};
    };
    auto result = get_columns();
    std::map<std::string, int> map_column = result.first;
    int iLineFirst = result.second;
    int satellite = get_satellite();
    int equ_time = get_equator_time();
    std::string strTime = "time rel. to equator passage [s]";
    int colTime = map_column.at(strTime);
    std::string strLon = "longitude [degrees_east]";
    int colLon = map_column.at(strLon);
    std::string strLat = "latitude [degrees_north]";
    int colLat = map_column.at(strLat);
    std::string search_swh = "u-band";
    std::string search_wind = "altimeter wind speed";
    std::vector<int> LIdxSwh, LIdxWind;
    for (auto & kv: map_column) {
      std::string const& name = kv.first;
      if (startswith(name, search_swh)) {
        LIdxSwh.push_back(kv.second);
      }
      if (startswith(name, search_wind)) {
        LIdxWind.push_back(kv.second);
      }
    }
    for (int iLine=iLineFirst; iLine<ListLine.size(); iLine++) {
      std::vector<std::string> LStr = STRING_Split(ListLine[iLine], " ");
      std::vector<double> V;
      for (auto & eStr : LStr) {
        double val = ParseScalar<double>(eStr);
        V.push_back(val);
      }
      double time_rel = V[colTime];
      double time = equ_time + time_rel / 86400.0;
      if (time < minTime) {
        minTime = time;
      }
      if (time > maxTime) {
        maxTime = time;
      }
      double lon = V[colLon];
      double lat = V[colLat];
      double sum_swh = 0;
      for (auto& iCol : LIdxSwh) {
        sum_swh += V[iCol];
      }
      sum_swh /= LIdxSwh.size();
      double sum_wind = 0;
      for (auto& iCol : LIdxWind) {
        sum_wind += V[iCol];
      }
      sum_wind /= LIdxWind.size();
      SingleEntryMeasurement eEnt = GetSingleEntryMeasurement();
      eEnt.Time = time;
      eEnt.Lon = lon;
      eEnt.Lat = lat;
      eEnt.Swh_used = sum_swh;
      eEnt.WindSpeed_used = sum_wind;
      eEnt.Satellite = satellite;
      ListEnt.push_back(eEnt);
    }
  }
  std::cerr << "minTime = " << DATE_ConvertMjd2mystringPres(minTime) << " maxTime=" << DATE_ConvertMjd2mystringPres(maxTime) << "\n";
  std::cerr << "|ListEnt|=" << ListEnt.size() << "\n";
  std::map<int, std::vector<SingleEntryMeasurement>> map;
  for (auto & eEnt : ListEnt) {
    int eKey = DATE_GetDayKey(eEnt.Time);
    map[eKey].push_back(eEnt);
  }
  std::cerr << "|map|=" << map.size() << "\n";
  for (auto & kv : map) {
    std::vector<SingleEntryMeasurement> LocEnt = kv.second;
    std::stable_sort(ListEnt.begin(), ListEnt.end(), [&](SingleEntryMeasurement const& x, SingleEntryMeasurement const& y) -> bool {
      return x.Time < y.Time;
    });
    double time = LocEnt[0].Time;
    std::vector<int> eDate = DATE_ConvertMjd2six(time);
    int n_ent = LocEnt.size();
    int eYear = eDate[0];
    int eMonth = eDate[1];
    int eDay = eDate[2];
    std::string strYear = std::to_string(eYear);
    std::string strMonth = StringNumber(eMonth, 2);
    std::string strDay = StringNumber(eDay, 2);
    std::string eFileO = "wm_" + strYear + strMonth + strDay + ".nc";
    std::string FullFile = PrefixO + strYear + "/" + strMonth + "/" + eFileO;
    std::string command = "mkdir -p " + PrefixO + strYear + "/" + strMonth;
    int iret = system(command.c_str());
    if (iret != 0) {
      std::cerr << "Failed in creation of directory\n";
      throw TerminalException{1};
    }
    std::vector<double> l_time(n_ent);
    std::vector<double> l_lat(n_ent), l_lon(n_ent);
    std::vector<double> l_wind_speed(n_ent);
    std::vector<double> l_blank(n_ent);
    std::vector<double> l_swh(n_ent);
    std::vector<int> l_satellite(n_ent);
    double RefTime = DATE_ConvertSix2mjd({1900, 1, 1, 0, 0, 0});
    for (size_t i=0; i<n_ent; i++) {
      l_time[i] = LocEnt[i].Time - RefTime;
      l_lat[i] = LocEnt[i].Lat;
      l_lon[i] = LocEnt[i].Lon;
      l_wind_speed[i] = LocEnt[i].WindSpeed_used;
      l_swh[i] = LocEnt[i].Swh_used;
      l_satellite[i] = LocEnt[i].Satellite;
    }
    std::cerr << "Writing FullFile=" << FullFile << "\n";
    netCDF::NcFile datafile(FullFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
    datafile.addDim("mes", n_ent);
    std::vector<std::string> LDim{"mes"};
    netCDF::NcVar varTime = datafile.addVar("time", "double", LDim);
    std::vector<double> dataMiss1{-1};
    std::vector<double> dataAddOffset1{0};
    std::vector<double> dataScaleFactor1{1};
    std::vector<double> dataValidRange1{0, 401767};
    varTime.putAtt("long_name", "time");
    varTime.putAtt("units", "days since 1900-1-1 0:0:0");
    varTime.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset1.data());
    varTime.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor1.data());
    varTime.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange1.data());
    varTime.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss1.data());
    varTime.putVar(l_time.data());
    //
    netCDF::NcVar varLat = datafile.addVar("lat", "double", LDim);
    std::vector<double> dataMiss2{-10000};
    std::vector<double> dataAddOffset2{0};
    std::vector<double> dataScaleFactor2{1};
    std::vector<double> dataValidRange2{-180, 180};
    varLat.putAtt("long_name", "latitude in degrees north");
    varLat.putAtt("units", "degrees_north");
    varLat.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset2.data());
    varLat.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor2.data());
    varLat.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange2.data());
    varLat.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss2.data());
    varLat.putVar(l_lat.data());
    //
    netCDF::NcVar varLon = datafile.addVar("lon", "double", LDim);
    std::vector<double> dataMiss3{-10000};
    std::vector<double> dataAddOffset3{0};
    std::vector<double> dataScaleFactor3{1};
    std::vector<double> dataValidRange3{-180, 180};
    varLon.putAtt("long_name", "longitude in degrees east");
    varLon.putAtt("units", "degrees_east");
    varLon.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset3.data());
    varLon.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor3.data());
    varLon.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange3.data());
    varLon.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss3.data());
    varLon.putVar(l_lon.data());
    //
    netCDF::NcVar varWindSpeed = datafile.addVar("wind_speed", "double", LDim);
    std::vector<double> dataMiss4{-10000};
    std::vector<double> dataAddOffset4{0};
    std::vector<double> dataScaleFactor4{1};
    std::vector<double> dataValidRange4{0, 50};
    varWindSpeed.putAtt("long_name", "wind speed");
    varWindSpeed.putAtt("units", "m s-1");
    varWindSpeed.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varWindSpeed.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varWindSpeed.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varWindSpeed.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varWindSpeed.putVar(l_wind_speed.data());
    //
    netCDF::NcVar varWindSpeedCor = datafile.addVar("wind_speed_cor", "double", LDim);
    varWindSpeedCor.putAtt("long_name", "wind speed corrected");
    varWindSpeedCor.putAtt("units", "m s-1");
    varWindSpeedCor.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varWindSpeedCor.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varWindSpeedCor.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varWindSpeedCor.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varWindSpeedCor.putVar(l_wind_speed.data());
    //
    netCDF::NcVar varSigma0 = datafile.addVar("sigma0", "double", LDim);
    varSigma0.putAtt("long_name", "sigma0");
    varSigma0.putAtt("units", "dB");
    varSigma0.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSigma0.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSigma0.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSigma0.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSigma0.putVar(l_blank.data());
    //
    netCDF::NcVar varSigma0cal = datafile.addVar("sigma0_cal", "double", LDim);
    varSigma0cal.putAtt("long_name", "sigma0 calibrated");
    varSigma0cal.putAtt("units", "dB");
    varSigma0cal.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSigma0cal.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSigma0cal.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSigma0cal.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSigma0cal.putVar(l_blank.data());
    //
    netCDF::NcVar varSigma0std = datafile.addVar("sigma0std", "double", LDim);
    varSigma0std.putAtt("long_name", "sigma0 standard deviation");
    varSigma0std.putAtt("units", "dB");
    varSigma0std.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSigma0std.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSigma0std.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSigma0std.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSigma0std.putVar(l_blank.data());
    //
    netCDF::NcVar varSigma0second = datafile.addVar("sigma0second", "double", LDim);
    varSigma0second.putAtt("long_name", "sigma0 second");
    varSigma0second.putAtt("units", "dB");
    varSigma0second.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSigma0second.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSigma0second.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSigma0second.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSigma0second.putVar(l_blank.data());
    //
    netCDF::NcVar varSigma0secondstd = datafile.addVar("sigma0secondstd", "double", LDim);
    varSigma0secondstd.putAtt("long_name", "sigma0 second standard deviation");
    varSigma0secondstd.putAtt("units", "dB");
    varSigma0secondstd.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSigma0secondstd.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSigma0secondstd.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSigma0secondstd.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSigma0secondstd.putVar(l_blank.data());
    //
    netCDF::NcVar varSwh = datafile.addVar("swh", "double", LDim);
    varSwh.putAtt("long_name", "sea_surface_wave_significant_height");
    varSwh.putAtt("units", "m");
    varSwh.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSwh.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSwh.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSwh.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSwh.putVar(l_swh.data());
    //
    netCDF::NcVar varSwhstd = datafile.addVar("swhstd", "double", LDim);
    varSwhstd.putAtt("long_name", "sea_surface_wave_significant_height_standard_deviation");
    varSwhstd.putAtt("units", "m");
    varSwhstd.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSwhstd.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSwhstd.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSwhstd.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSwhstd.putVar(l_blank.data());
    //
    netCDF::NcVar varSwhcor = datafile.addVar("swhcor", "double", LDim);
    varSwhcor.putAtt("long_name", "sea_surface_wave_significant_height_corrected");
    varSwhcor.putAtt("units", "m");
    varSwhcor.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSwhcor.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSwhcor.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSwhcor.putAtt("_FillValue", netCDF::NcType::nc_DOUBLE, 1, dataMiss4.data());
    varSwhcor.putVar(l_swh.data());
    //
    netCDF::NcVar varSatellite = datafile.addVar("satellite", "int", LDim);
    varSatellite.putAtt("long_name", "satellite");
    varSatellite.putAtt("units", "m");
    varSatellite.putAtt("add_offset", netCDF::NcType::nc_DOUBLE, 1, dataAddOffset4.data());
    varSatellite.putAtt("scale_offset", netCDF::NcType::nc_DOUBLE, 1, dataScaleFactor4.data());
    varSatellite.putAtt("valid_range", netCDF::NcType::nc_DOUBLE, 2, dataValidRange4.data());
    varSatellite.putVar(l_satellite.data());
  }
}


// clang-format off
#endif // SRC_OCEAN_SATELLITE_H_
// clang-format on
