// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_COMMONFUNCMODEL_H_
#define SRC_OCEAN_COMMONFUNCMODEL_H_

#include "ArrHistory.h"
#include "Namelist.h"
#include "mjdv2.h"
#include <algorithm>
#include <string>
#include <vector>

std::vector<VarQuery> GetIntervalFLD_query(std::vector<double> const &ListTime,
                                           double const &TimeFrameDay) {
  int nbTime = ListTime.size();
  std::vector<VarQuery> ListQuery(nbTime);
  for (int iTime = 0; iTime < nbTime; iTime++) {
    VarQuery eQuery;
    eQuery.eTimeDay = ListTime[iTime];
    eQuery.iTime = iTime;
    eQuery.TimeFrameDay = TimeFrameDay;
    eQuery.typeQuery = "direct";
    ListQuery[iTime] = eQuery;
  }
  return ListQuery;
}

std::vector<VarQuery> GetIntervalGen_Kernel(SingleBlock const &eBlock,
                                            double const &PropFirstTime,
                                            double const &PropLastTime,
                                            double const &PropDeltaInterval) {
  std::string BEGTC = eBlock.ListStringValues.at("BEGTC");
  std::string ENDTC = eBlock.ListStringValues.at("ENDTC");
  double DELTC = eBlock.ListDoubleValues.at("DELTC");
  std::string UNITC = eBlock.ListStringValues.at("UNITC");
  std::cerr << "We have UNITC\n";
  std::string KindSelect = eBlock.ListStringValues.at("KindSelect");
  std::cerr << "We have KindSelect\n";
  double TimeFrameDay = eBlock.ListDoubleValues.at("TimeFrameDay");
  std::cerr << "We have TimeFrameDay\n";
  if (PositionVect({"direct", "seasonal", "seasonalIvica", "monthly", "yearly",
                    "specific"},
                   KindSelect) == -1) {
    std::cerr << "Allowed values for KindSelect are direct, seasonal, "
                 "seasonalIvica, monthly, specific\n";
    std::cerr << "KindSelect=" << KindSelect << "\n";
    throw TerminalException{1};
  }
  double FirstTime, LastTime, DeltaInterval;
  //
  if (UNITC == "automatic") {
    if (PropDeltaInterval < 0) {
      std::cerr << "PropDeltaInterval has not been set. This should mean that "
                   "no model has been put on input\n";
      std::cerr << "Combined with the choice of \"automatic\" for first time "
                   "this make the system impossible to determine the time\n";
      throw TerminalException{1};
    }
    DeltaInterval = PropDeltaInterval;
  } else {
    DeltaInterval = GetIntervalSize(DELTC, UNITC);
  }
  //
  if (BEGTC == "earliest") {
    if (PropFirstTime < 0) {
      std::cerr << "PropFirstTime has not been set. This should mean that no "
                   "model has been put on input\n";
      std::cerr << "Combined with the choice of \"earliest\" for first time "
                   "this make the system impossible to determine the time\n";
      throw TerminalException{1};
    }
    FirstTime = PropFirstTime;
  } else {
    FirstTime = CT2MJD(BEGTC);
  }
  //
  if (ENDTC == "latest") {
    if (PropLastTime < 0) {
      std::cerr << "PropLastTime has not been set. This should mean that no "
                   "model has been put on input\n";
      std::cerr << "Combined with the choice of \"earliest\" for first time "
                   "this make the system impossible to determine the time\n";
      throw TerminalException{1};
    }
    LastTime = PropLastTime;
  } else {
    LastTime = CT2MJD(ENDTC);
  }
  //
  double tolDay = DeltaInterval / static_cast<double>(1000);
  if (LastTime < FirstTime - tolDay) {
    std::cerr << "Error in GetIntervalGen_Kernel\n";
    std::cerr << "We should have ENDTC >= BEGTC. But instead we have:\n";
    std::cerr << "    BEGTC = " << BEGTC << "\n";
    std::cerr << "    ENDTC = " << ENDTC << "\n";
    std::cerr << "FirstTime = " << FirstTime << "\n";
    std::cerr << " LastTime = " << LastTime << "\n";
    std::cerr << "   tolDay = " << tolDay << "\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  //
  if (KindSelect == "direct") {
    std::cerr << "    FirstTime = " << FirstTime << "\n";
    std::cerr << "     LastTime = " << LastTime << "\n";
    std::cerr << "DeltaInterval = " << DeltaInterval << "\n";
    std::vector<double> ListTime =
        GetIntervalFLD(FirstTime, LastTime, DeltaInterval);
    return GetIntervalFLD_query(ListTime, TimeFrameDay);
  }
  if (KindSelect == "monthly")
    return GetIntervalFLmonthly(FirstTime, LastTime);
  if (KindSelect == "yearly")
    return GetIntervalFLyearly(FirstTime, LastTime);
  if (KindSelect == "seasonal")
    return GetIntervalFLseasonal(FirstTime, LastTime);
  if (KindSelect == "seasonalIvica")
    return GetIntervalFLseasonalIvica(FirstTime, LastTime);
  if (KindSelect == "specific") {
    std::vector<std::string> ListStringTime =
        eBlock.ListListStringValues.at("ListSpecificTimes");
    std::vector<double> ListTime;
    for (auto &eTimeStr : ListStringTime)
      ListTime.push_back(CT2MJD(eTimeStr));
    return GetIntervalFLD_query(ListTime, TimeFrameDay);
  }
  std::cerr << "We should not reach that stage\n";
  throw TerminalException{1};
}

std::vector<VarQuery>
GetIntervalGen_Query(SingleBlock const &eBlock,
                     std::vector<ArrayHistory> const &ListArr) {
  int nbArr = ListArr.size();
  double PropFirstTime;
  double PropLastTime;
  double PropDeltaTime;
  if (nbArr > 0) {
    std::vector<double> ListFirstTime(nbArr);
    std::vector<double> ListLastTime(nbArr);
    std::vector<double> ListDeltaTime(nbArr);
    for (int iArr = 0; iArr < nbArr; iArr++) {
      ListFirstTime[iArr] = MinimumTimeHistoryArray(ListArr[iArr]);
      ListLastTime[iArr] = MaximumTimeHistoryArray(ListArr[iArr]);
      ListDeltaTime[iArr] = ComputeDeltaTimeHistoryArray(ListArr[iArr]);
    }
    PropFirstTime = VectorMax(ListFirstTime);
    PropLastTime = VectorMin(ListLastTime);
    PropDeltaTime = VectorAvg(ListDeltaTime);
  } else {
    PropFirstTime = -1;
    PropLastTime = -1;
    PropDeltaTime = -1;
  }
  std::cerr << "PropFirstTime=" << PropFirstTime
            << " PropLastTime=" << PropLastTime
            << " PropDeltaTime=" << PropDeltaTime << "\n";
  return GetIntervalGen_Kernel(eBlock, PropFirstTime, PropLastTime,
                               PropDeltaTime);
}

std::vector<double> GetIntervalGen(SingleBlock const &eBlock,
                                   std::vector<ArrayHistory> const &ListArr) {
  std::vector<double> ListTime;
  for (auto &eQuery : GetIntervalGen_Query(eBlock, ListArr))
    ListTime.push_back(eQuery.eTimeDay);
  return ListTime;
}

//
// Files ordered sequentially.
// We use ocean_time to order them.
ArrayHistory Sequential_ReadArrayHistory(std::string const &HisPrefix,
                                         std::string const &strTime) {
  std::cerr << "Beginning of Sequential_ReadArrayHistory\n";
  int len = HisPrefix.length();
  std::vector<int> ListPos;
  for (int iChar = 0; iChar < len; iChar++) {
    std::string eChar = HisPrefix.substr(iChar, 1);
    if (eChar == "/")
      ListPos.push_back(iChar);
  }
  if (ListPos.size() == 0) {
    std::cerr << "We should have at least one / in an HisPrefix of "
                 "SequentialReadArrayHistory\n";
    throw TerminalException{1};
  }
  int iCharLast = ListPos[ListPos.size() - 1];
  std::string eDir = HisPrefix.substr(0, iCharLast + 1);
  std::string RawPrefix = HisPrefix.substr(iCharLast + 1, len - iCharLast - 1);
  if (RawPrefix == "") {
    std::cerr << "HisPrefix=" << HisPrefix << "\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "RawPrefix=" << RawPrefix << "\n";
    std::cerr << "We should have RawPrefix being nontrivial\n";
    std::cerr << "For example RepoHistory/WWM_output_\n";
    throw TerminalException{1};
  }
  std::cerr << "Before FILE_GetDirectoryListFile eDir=" << eDir << "\n";
  std::vector<std::string> PreListFile = FILE_GetDirectoryListFile(eDir);
  std::cerr << "|PreListFile|=" << PreListFile.size() << "\n";
  std::vector<std::string> ListFileNames;
  for (auto &eFile : PreListFile) {
    std::string eFileTot = eDir + eFile;
    std::vector<std::string> LStr = STRING_Split(eFileTot, RawPrefix);
    size_t nbBlock = LStr.size();
    std::cerr << "eFileTot=" << eFileTot << " RawPrefix=" << RawPrefix
              << " |LStr|=" << nbBlock << "\n";
    if (nbBlock == 2)
      ListFileNames.push_back(eFileTot);
  }
  std::cerr << "|ListFileNames|=" << ListFileNames.size() << "\n";
  struct FullEntry {
    int iFile;
    int iTime;
    double eTime;
  };
  std::vector<FullEntry> ListFull;
  int nbFile = ListFileNames.size();
  std::cerr << "nbFile=" << nbFile << "\n";
  for (int iFile = 0; iFile < nbFile; iFile++) {
    std::string eFile = ListFileNames[iFile];
    std::vector<double> LTime = NC_ReadTimeFromFile(eFile, strTime);
    int nbTime = LTime.size();
    for (int iTime = 0; iTime < nbTime; iTime++) {
      double eTime = LTime[iTime];
      FullEntry eFull{iFile, iTime, eTime};
      ListFull.push_back(eFull);
    }
  }
  std::cerr << "Now |ListFull|=" << ListFull.size() << "\n";
  sort(ListFull.begin(), ListFull.end(),
       [&](FullEntry const &a, FullEntry const &b) -> bool {
         if (a.eTime < b.eTime)
           return true;
         return false;
       });
  int nbFull = ListFull.size();
  if (nbFull == 0) {
    std::cerr << "nbFull=" << nbFull << "\n";
    std::cerr << "That is we found ZERO files to be usable to us\n";
    std::cerr << "Please correct\n";
    std::cerr << "HisPrefix=" << HisPrefix << "\n";
    throw TerminalException{1};
  }
  std::vector<double> ListTime(nbFull);
  std::vector<int> ListIFile(nbFull);
  std::vector<int> ListIRec(nbFull);
  for (int i = 0; i < nbFull; i++) {
    ListTime[i] = ListFull[i].eTime;
    ListIFile[i] = ListFull[i].iFile;
    ListIRec[i] = ListFull[i].iTime;
  }
  double FirstTime = ListTime[0];
  double LastTime = ListTime[nbFull - 1];
  ArrayHistory eArr;
  eArr.nbFile = nbFile;
  eArr.nbTime = nbFull;
  eArr.FirstTime = FirstTime;
  eArr.LastTime = LastTime;
  eArr.ListIFile = ListIFile;
  eArr.ListIRec = ListIRec;
  eArr.ListFileNames = ListFileNames;
  eArr.ListTime = ListTime;
  eArr.AppendVarName = false;
  eArr.KindArchive = "NETCDF";
  eArr.TimeSteppingInfo = "classic";
  std::cerr << "Sequential array has been completed. Leaving\n";
  return eArr;
}

// Code common to WWM, COSMO, ROMS, WAM
ArrayHistory NC_ReadArrayHistory_Kernel(std::string const &HisPrefix,
                                        std::string const &StringTime,
                                        int const &nbDigit) {
  double FirstTime, LastTime;
  std::vector<std::string> ListFileNames;
  std::vector<int> ListIFile;
  std::vector<int> ListIRec;
  std::vector<double> ListTime;
  ArrayHistory eArr;
  int nbTime = 0;
  if (IsExistingFile(HisPrefix)) {
    std::cerr << "StringTime=" << StringTime << "\n";
    std::vector<double> LTime = NC_ReadTimeFromFile(HisPrefix, StringTime);
    ListFileNames.push_back(HisPrefix);
    int siz = LTime.size();
    for (int i = 0; i < siz; i++) {
      ListIFile.push_back(0);
      ListIRec.push_back(i);
      ListTime.push_back(LTime[i]);
    }
    FirstTime = ListTime[0];
    LastTime = ListTime[siz - 1];
    double TheSep = GetListTimeSeparation(ListTime);
    if (TheSep < 0) {
      eArr.TimeSteppingInfo = "classic";
    } else {
      eArr.TimeSteppingInfo = "singlefile";
      eArr.SeparationTime = TheSep;
    }
  } else {
    int iFileBegin = 0;
    while (true) {
      iFileBegin++;
      std::string TheHisFile =
          HisPrefix + StringNumber(iFileBegin, nbDigit) + ".nc";
      if (IsExistingFile(TheHisFile))
        break;
      if (iFileBegin == 9999) {
        std::cerr << "maybe you specified wrong HisPrefix\n";
        std::cerr << "HisPrefix = " << HisPrefix << "\n";
        std::cerr << " there is no files  HisPrefix????.nc\n";
        std::cerr << "Please correct\n";
        throw TerminalException{1};
      }
    }
    int iFileEnd = iFileBegin;
    while (true) {
      std::string TheHisFile =
          HisPrefix + StringNumber(iFileEnd + 1, nbDigit) + ".nc";
      if (!IsExistingFile(TheHisFile))
        break;
      iFileEnd++;
    }
    std::string TheHisFileBegin =
        HisPrefix + StringNumber(iFileBegin, nbDigit) + ".nc";
    std::vector<double> LTimeBegin =
        NC_ReadTimeFromFile(TheHisFileBegin, StringTime);
    int nbRecBegin = LTimeBegin.size();
    double DeltaTime;
    if (nbRecBegin > 1) {
      DeltaTime = LTimeBegin[1] - LTimeBegin[0];
    } else {
      if (iFileEnd > iFileBegin) {
        std::string TheHisFile =
            HisPrefix + StringNumber(iFileBegin + 1, nbDigit) + ".nc";
        std::vector<double> LTimeBP1 =
            NC_ReadTimeFromFile(TheHisFile, StringTime);
        DeltaTime = LTimeBP1[0] - LTimeBegin[0];
      } else {
        DeltaTime = 0;
      }
    }
    // Determination of number of arrays
    int iFileMiddle, nbRecMiddle;
    if (iFileEnd != iFileBegin) {
      iFileMiddle = iFileBegin + 1;
      std::string TheHisFile =
          HisPrefix + StringNumber(iFileMiddle, nbDigit) + ".nc";
      std::vector<double> LTimeMiddle =
          NC_ReadTimeFromFile(TheHisFile, StringTime);
      nbRecMiddle = LTimeMiddle.size();
    } else {
      iFileMiddle = iFileBegin;
      nbRecMiddle = nbRecBegin;
    }
    int nbRecEnd;
    if (iFileEnd != iFileMiddle) {
      std::string TheHisFile =
          HisPrefix + StringNumber(iFileEnd, nbDigit) + ".nc";
      std::vector<double> LTimeEnd =
          NC_ReadTimeFromFile(TheHisFile, StringTime);
      nbRecEnd = LTimeEnd.size();
    } else {
      nbRecEnd = nbRecMiddle;
    }
    int nbFile = 1 + iFileEnd - iFileBegin;
    std::vector<int> NbPerArray(nbFile);
    for (int iFile = 0; iFile < nbFile; iFile++)
      NbPerArray[iFile] = nbRecMiddle;
    NbPerArray[0] = nbRecBegin;
    NbPerArray[nbFile - 1] = nbRecEnd;
    double currentTime = LTimeBegin[0];
    FirstTime = currentTime;
    for (int iFile = 0; iFile < nbFile; iFile++) {
      int iFileTot = iFile + iFileBegin;
      std::string TheHisFile =
          HisPrefix + StringNumber(iFileTot, nbDigit) + ".nc";
      ListFileNames.push_back(TheHisFile);
      int siz = NbPerArray[iFile];
      for (int i = 0; i < siz; i++) {
        ListIFile.push_back(iFile);
        ListIRec.push_back(i);
        ListTime.push_back(currentTime);
        if (iFile != nbFile - 1 || i != siz - 1)
          currentTime += DeltaTime;
      }
    }
    LastTime = currentTime;
    eArr.TimeSteppingInfo = "multiplenetcdf";
    eArr.SeparationTime = DeltaTime;
    eArr.nbRecBegin = nbRecBegin;
    eArr.nbRecMiddle = nbRecMiddle;
    eArr.HisPrefix = HisPrefix;
  }
  int nbFile = ListFileNames.size();
  std::cerr << "nbFile=" << nbFile << "\n";
  eArr.FirstTime = FirstTime;
  eArr.LastTime = LastTime;
  //  std::string strPresFirst=DATE_ConvertMjd2mystringPres(FirstTime);
  //  std::string strPresLast =DATE_ConvertMjd2mystringPres(LastTime);
  //  std::cerr << "FirstTime=" << FirstTime << " date=" << strPresFirst <<
  //  "\n"; std::cerr << " LastTime=" << LastTime  << " date=" << strPresLast <<
  //  "\n"; exit(1);
  if (eArr.TimeSteppingInfo != "multiplenetcdf") {
    nbTime = ListTime.size();
    eArr.nbFile = nbFile;
    eArr.ListFileNames = ListFileNames;
    eArr.ListIFile = ListIFile;
    eArr.ListIRec = ListIRec;
    eArr.ListTime = ListTime;
  } else {
    eArr.nbFile = 0;
  }
  eArr.nbTime = nbTime;
  eArr.nbDigit = nbDigit;
  eArr.AppendVarName = false;
  eArr.KindArchive = "NETCDF";
  //  std::cerr << "Leaving NC_ReadArrayHistory_Kernel\n";
  //  std::cerr << "TimeSteppingInfo=" << eArr.TimeSteppingInfo << "\n";
  return eArr;
}

ArrayHistory WW3_ReadArrayHistory(std::string const &HisFile,
                                  std::string const &HisPrefix) {
  std::vector<double> LTime = NC_ReadTimeFromFile(HisFile, "time");
  int nbTime = LTime.size();
  std::vector<int> ListIRec, ListIFile;
  for (int iTime = 0; iTime < nbTime; iTime++) {
    ListIFile.push_back(0);
    ListIRec.push_back(iTime);
  }
  double FirstTime = LTime[0];
  double LastTime = LTime[nbTime - 1];
  std::vector<std::string> ListFileNames;
  std::cerr << "NC_ReadArrayHistory, HisPrefix=" << HisPrefix << "\n";
  ListFileNames.push_back(HisPrefix);
  ArrayHistory eArr;
  eArr.nbFile = 1;
  eArr.nbTime = nbTime;
  eArr.FirstTime = FirstTime;
  eArr.LastTime = LastTime;
  eArr.ListFileNames = ListFileNames;
  std::cerr << "|ListFileNames|=" << ListFileNames.size() << "\n";
  eArr.ListIFile = ListIFile;
  eArr.ListIRec = ListIRec;
  eArr.ListTime = LTime;
  eArr.AppendVarName = false;
  eArr.KindArchive = "NETCDF";
  eArr.TimeSteppingInfo = "classic";
  return eArr;
}

bool NETCDF_TestVariableAccessSpecTime(TotalArrGetData const &TotalArr,
                                       std::string const &eVar,
                                       double const &eTimeDay) {
  InterpInfo eInterpInfo =
      GetTimeInterpolationInfoGeneralized(TotalArr.eArr, eTimeDay);
  if (eInterpInfo.UseSingleEntry) {
    int iTime = eInterpInfo.iTimeLow;
    std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTime);
    int iFile = eRecLow[0];
    std::string HisFile = ARR_GetHisFileName(TotalArr.eArr, eVar, iFile);
    bool test = NC_IsVar(HisFile, eVar);
    return test;
  }
  int iTimeLow = eInterpInfo.iTimeLow;
  int iTimeUpp = eInterpInfo.iTimeUpp;
  std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTimeLow);
  std::vector<int> eRecUpp = GetIFileIRec(TotalArr.eArr, iTimeUpp);
  int iFileLow = eRecLow[0];
  int iFileUpp = eRecUpp[0];
  std::string HisFileLow = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileLow);
  std::string HisFileUpp = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileUpp);
  bool testLow = NC_IsVar(HisFileLow, eVar);
  bool testUpp = NC_IsVar(HisFileUpp, eVar);
  return testLow && testUpp;
}


MyMatrix<double> NETCDF_Get2DvariableSpecTime(TotalArrGetData const &TotalArr,
                                              std::string const &eVar,
                                              double const &eTimeDay) {
  InterpInfo eInterpInfo =
      GetTimeInterpolationInfoGeneralized(TotalArr.eArr, eTimeDay);
  if (eInterpInfo.UseSingleEntry) {
    int iTime = eInterpInfo.iTimeLow;
    std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTime);
    int iFile = eRecLow[0];
    int iRec = eRecLow[1];
    std::string HisFile = ARR_GetHisFileName(TotalArr.eArr, eVar, iFile);
    MyMatrix<double> RMat =
        NETCDF_Get2DvariableSpecEntry(HisFile, TotalArr.GrdArr, eVar, iRec);
    return RMat;
  }
  double alphaLow = eInterpInfo.alphaLow;
  double alphaUpp = eInterpInfo.alphaUpp;
  int iTimeLow = eInterpInfo.iTimeLow;
  int iTimeUpp = eInterpInfo.iTimeUpp;
  std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTimeLow);
  std::vector<int> eRecUpp = GetIFileIRec(TotalArr.eArr, iTimeUpp);
  int iFileLow = eRecLow[0];
  int iFileUpp = eRecUpp[0];
  int iRecLow = eRecLow[1];
  int iRecUpp = eRecUpp[1];
  std::string HisFileLow = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileLow);
  std::string HisFileUpp = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileUpp);
  MyMatrix<double> eVarLow =
      NETCDF_Get2DvariableSpecEntry(HisFileLow, TotalArr.GrdArr, eVar, iRecLow);
  MyMatrix<double> eVarUpp =
      NETCDF_Get2DvariableSpecEntry(HisFileUpp, TotalArr.GrdArr, eVar, iRecUpp);
  return alphaLow * eVarLow + alphaUpp * eVarUpp;
}

Eigen::Tensor<double, 3>
NETCDF_Get3DvariableSpecTime(TotalArrGetData const &TotalArr,
                             std::string const &eVar, double const &eTimeDay) {
  InterpInfo eInterpInfo =
      GetTimeInterpolationInfoGeneralized(TotalArr.eArr, eTimeDay);
  //  std::cerr << "eInterpInfo iTimeLow=" << eInterpInfo.iTimeLow
  //            << " iTimeUpp=" << eInterpInfo.iTimeUpp << "\n";
  //  std::cerr << "eInterpInfo alphaLow=" << eInterpInfo.alphaLow
  //            << " alphaUpp=" << eInterpInfo.alphaUpp << "\n";
  if (eInterpInfo.UseSingleEntry) {
    int iTime = eInterpInfo.iTimeLow;
    std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTime);
    int iFile = eRecLow[0];
    int iRec = eRecLow[1];
    std::string HisFile = ARR_GetHisFileName(TotalArr.eArr, eVar, iFile);
    return NETCDF_Get3DvariableSpecEntry(HisFile, TotalArr.GrdArr, eVar, iRec);
  }
  double alphaLow = eInterpInfo.alphaLow;
  double alphaUpp = eInterpInfo.alphaUpp;
  int iTimeLow = eInterpInfo.iTimeLow;
  int iTimeUpp = eInterpInfo.iTimeUpp;
  std::vector<int> eRecLow = GetIFileIRec(TotalArr.eArr, iTimeLow);
  std::vector<int> eRecUpp = GetIFileIRec(TotalArr.eArr, iTimeUpp);
  int iFileLow = eRecLow[0];
  int iFileUpp = eRecUpp[0];
  int iRecLow = eRecLow[1];
  int iRecUpp = eRecUpp[1];
  std::string HisFileLow = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileLow);
  std::string HisFileUpp = ARR_GetHisFileName(TotalArr.eArr, eVar, iFileUpp);
  Eigen::Tensor<double, 3> eVarLow =
      NETCDF_Get3DvariableSpecEntry(HisFileLow, TotalArr.GrdArr, eVar, iRecLow);
  Eigen::Tensor<double, 3> eVarUpp =
      NETCDF_Get3DvariableSpecEntry(HisFileUpp, TotalArr.GrdArr, eVar, iRecUpp);
  auto LDim = eVarLow.dimensions();
  int s_vert = LDim[0];
  int eta = LDim[1];
  int xi = LDim[2];
  //  does not COMPILE:
  //  Eigen::Tensor<double,3> RetVar=alphaLow*eVarLow + alphaUpp*eVarUpp;
  Eigen::Tensor<double, 3> RetVar(s_vert, eta, xi);
  for (int k = 0; k < s_vert; k++)
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++)
        RetVar(k, i, j) =
            alphaLow * eVarLow(k, i, j) + alphaUpp * eVarUpp(k, i, j);
  return RetVar;
}

// clang-format off
#endif  // SRC_OCEAN_COMMONFUNCMODEL_H_
// clang-format on
