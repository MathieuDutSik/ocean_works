#ifndef INCLUDE_ARR_HISTORY
#define INCLUDE_ARR_HISTORY

#include "Temp_common.h"
#include "Basic_Ocean_types.h"


void ARR_PrintHistoryArray(std::ostream & os, ArrayHistory const& eArr)
{
  int nbTime=eArr.nbTime;
  int nbFile=eArr.nbFile;
  os << "------------------- Printing history array -------------------\n";
  os << "TimeSteppingInfo=" << eArr.TimeSteppingInfo << "\n";
  os << "nbTime=" << nbTime << "\n";
  int nbFileReal=eArr.ListFileNames.size();
  os << "nbFile=" << nbFile << " nbFileReal=" << nbFileReal << "\n";
  for (int iFileReal=0; iFileReal<nbFileReal; iFileReal++) {
    os << "iFileReal=" << iFileReal << " eFile=" << eArr.ListFileNames[iFileReal] << "\n";
  }
  std::cerr << "After iFileReal loop\n";
  std::cerr << "|ListIFile| = " << eArr.ListIFile.size() << "\n";
  std::cerr << "|ListIRec|  = " << eArr.ListIRec.size() << "\n";
  std::cerr << "|ListTime|  = " << eArr.ListTime.size() << "\n";
  for (int iTime=0; iTime<nbTime; iTime++) {
    //    std::cerr << "iTime=" << iTime << " / " << nbTime << "\n";
    int iFile=eArr.ListIFile[iTime];
    //    std::cerr << "Step 1\n";
    int iRec=eArr.ListIRec[iTime];
    //    std::cerr << "Step 2\n";
    double eTime=eArr.ListTime[iTime];
    //    std::cerr << "Step 3\n";
    std::string eTimeStr=MJD2CT(eTime);
    os << "iTime=" << iTime << " iFile=" << iFile << " iRec=" << iRec << " time=" << eTimeStr << "\n";
  }
  os << "-----------------------------------------------\n";
  int nbMessage=eArr.ListAllMessage.size();
  std::cerr << "number of grib messages = " << nbMessage << "\n";
  for (int iMessage=0; iMessage<nbMessage; iMessage++) {
    GRIB_MessageInfo eMesg=eArr.ListAllMessage[iMessage];
    std::string strDate=DATE_ConvertMjd2mystringPres(eMesg.time);
    std::cerr << "iMessage=" << iMessage << " / " << nbMessage << " file=" << eMesg.FileName << " date=" << strDate << "\n";
  }

}




std::string ARR_GetHisFileName(ArrayHistory const& eArr, std::string const& eVar, int const& iFile)
{
  if (eArr.TimeSteppingInfo == "multiplenetcdf") {
    return eArr.HisPrefix + StringNumber(iFile+1,eArr.nbDigit) + ".nc";
  }
  int len=eArr.ListFileNames.size();
  if (iFile >= len) {
    ARR_PrintHistoryArray(std::cerr, eArr);
    std::cerr << "WW3 : eVar=" << eVar << "\n";
    std::cerr << "iFile=" << iFile << " len=" << len << "\n";
    std::cerr << "We need iFile < len\n";
    std::cerr << "Error. trying to get eArr.ListFileNames\n";
    std::cerr << "After the last values\n";
    throw TerminalException{1};
  }
  if (eArr.AppendVarName)
    return eArr.ListFileNames[iFile] + eVar + ".nc";
  else
    return eArr.ListFileNames[iFile];
};


std::string ARR_GetFirstFileName(ArrayHistory const& eArr)
{
  std::string eVar="unset_and_bug";
  int iFile=0;
  return ARR_GetHisFileName(eArr, eVar, iFile);
}




double MinimumTimeHistoryArray(ArrayHistory const& eArr)
{
  return eArr.FirstTime;
}

double MaximumTimeHistoryArray(ArrayHistory const& eArr)
{
  return eArr.LastTime;
}


double ComputeDeltaTimeHistoryArray(ArrayHistory const& eArr)
{
  int nbTime=eArr.ListTime.size();
  if (nbTime < 2) {
    return double(-1);
  }
  std::vector<double> ListDeltaTime(nbTime-1);
  for (int i=1; i<nbTime; i++) {
    double eDelta1 = eArr.ListTime[i] - eArr.ListTime[i-1];
    double eDelta2 = std::max(eDelta1, double(0));
    ListDeltaTime[i-1] = eDelta2;
  }
  double eDeltaMin=VectorMin(ListDeltaTime);
  double eDeltaMax=VectorMax(ListDeltaTime);
  std::cerr << "nbTime=" << nbTime << "\n";
  std::cerr << "eDeltaMin=" << eDeltaMin << " eDeltaMax=" << eDeltaMax << "\n";
  double eDelta=(eDeltaMin + eDeltaMax) / double(2);
  std::cerr << "eDelta=" << eDelta << "\n";
  //  throw TerminalException{1};
  return eDelta;
}





struct InterpInfo {
  int iTimeLow, iTimeUpp;
  double alphaLow, alphaUpp;
  bool UseSingleEntry;
};




InterpInfo GetTimeDifferentiationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
  InterpInfo eInterpInfo;
  double tolDay=double(1)/double(1000000);
  int nbTime=LTime.size();
  eInterpInfo.UseSingleEntry=false;
  if (eTimeDay < LTime[0] - tolDay) {
    std::cerr << "Error in GetTimeDifferentiationInfo\n";
    std::cerr << "The asked entry is before the first time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  if (eTimeDay > LTime[nbTime-1] + tolDay) {
    std::cerr << "Error in GetTimeDifferentiationInfo\n";
    std::cerr << "The asked entry is after the last time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  if (nbTime <= 1) {
    std::cerr << "We need at least two entries in order to do the time differential\n";
    std::cerr << "nbTime=" << nbTime << "\n";
    throw TerminalException{1};
  }
  for (int iTimeUpp=1; iTimeUpp<nbTime; iTimeUpp++) {
    int iTimeLow=iTimeUpp-1;
    double eTimeUpp=LTime[iTimeUpp];
    double eTimeLow=LTime[iTimeLow];
    if (eTimeLow - tolDay < eTimeDay && eTimeDay < eTimeUpp + tolDay) {
      double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
      double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
      eInterpInfo.iTimeLow=iTimeLow;
      eInterpInfo.iTimeUpp=iTimeUpp;
      eInterpInfo.alphaLow=alphaLow;
      eInterpInfo.alphaUpp=alphaUpp;
      return eInterpInfo;
    }
  }
  std::cerr << "Failed to find matching record\n";
  std::cerr << "Please debug\n";
  throw TerminalException{1};
}


struct InterpInfoDiff {
  int iTimeLow, iTimeUpp;
  double alphaLow, alphaUpp;
  bool UseSingleEntry;
};



std::vector<InterpInfoDiff> GRIB_GetTimeDifferentiationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
  double tolDay=double(1)/double(1000000);
  int nbTime=LTime.size();
  if (eTimeDay < LTime[0] - tolDay) {
    std::cerr << "Error in GetTimeDifferentiationInfo\n";
    std::cerr << "The asked entry is before the first time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  if (eTimeDay > LTime[nbTime-1] + tolDay) {
    std::cerr << "Error in GetTimeDifferentiationInfo\n";
    std::cerr << "The asked entry is after the last time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  if (nbTime <= 1) {
    std::cerr << "We need at least two entries in order to do the time differential\n";
    std::cerr << "nbTime=" << nbTime << "\n";
    throw TerminalException{1};
  }
  std::vector<InterpInfoDiff> ListReturn;
  for (int iTimeUpp=1; iTimeUpp<nbTime; iTimeUpp++) {
    int iTimeLow=iTimeUpp-1;
    double eTimeUpp=LTime[iTimeUpp];
    double eTimeLow=LTime[iTimeLow];
    if (eTimeLow - tolDay < eTimeDay && eTimeDay < eTimeUpp + tolDay) {
      double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
      double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
      InterpInfoDiff eInterpInfo;
      eInterpInfo.iTimeLow=iTimeLow;
      eInterpInfo.iTimeUpp=iTimeUpp;
      eInterpInfo.alphaLow=alphaLow;
      eInterpInfo.alphaUpp=alphaUpp;
      ListReturn.push_back(eInterpInfo);
    }
  }
  return ListReturn;
}






InterpInfo GetTimeInterpolationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
  InterpInfo eInterpInfo;
  double tolDay=double(1)/double(1000000);
  int nbTime=LTime.size();
  if (nbTime == 0) {
    std::cerr << "Error in GetTimeInterpolationInfo\n";
    std::cerr << "We cannot proceed because nbTime=0\n";
    throw TerminalException{1};
  }
  for (int iTime=0; iTime<nbTime; iTime++) {
    if (fabs(LTime[iTime] - eTimeDay) < tolDay) {
      //      std::cerr << "eDist=" << eDist << "\n";
      //      std::cerr << "eTimeDay=" << eTimeDay << "\n";
      //      std::cerr << "LTime[iTime]=" << LTime[iTime] << "\n";
      //      std::cerr << "iTime=" << iTime << "\n";
      eInterpInfo.UseSingleEntry=true;
      eInterpInfo.iTimeLow=iTime;
      eInterpInfo.iTimeUpp=iTime;
      eInterpInfo.alphaLow=1; // We need to have that in order to read values correctly
      eInterpInfo.alphaUpp=0; // some code do not make distinction between UseSingleEntry=T/F
      //                         and we want it still to work.
      //      std::cerr << "GetTimeInterpolationInfo, leaving case 1\n";
      return eInterpInfo;
    }
  }
  eInterpInfo.UseSingleEntry=false;
  if (eTimeDay < LTime[0] - tolDay) {
    std::cerr << "Error in GetTimeInterpolationInfo\n";
    std::cerr << "The asked entry is before the first time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  if (eTimeDay > LTime[nbTime-1] + tolDay) {
    std::cerr << "Error in GetTimeInterpolationInfo\n";
    std::cerr << "The asked entry is after the last time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    throw TerminalException{1};
  }
  for (int iTimeUpp=1; iTimeUpp<nbTime; iTimeUpp++) {
    int iTimeLow=iTimeUpp-1;
    double eTimeUpp=LTime[iTimeUpp];
    double eTimeLow=LTime[iTimeLow];
    double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
    double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
    if (alphaLow >= 0 && alphaUpp >= 0) {
      eInterpInfo.iTimeLow=iTimeLow;
      eInterpInfo.iTimeUpp=iTimeUpp;
      eInterpInfo.alphaLow=alphaLow;
      eInterpInfo.alphaUpp=alphaUpp;
      //      std::cerr << "GetTimeInterpolationInfo, leaving case 2\n";
      return eInterpInfo;
    }
  }
  std::cerr << "Failed to find matching record\n";
  std::cerr << "Please debug\n";
  throw TerminalException{1};
}

InterpInfo GetTimeInterpolationInfo_infinite(double const& FirstTime, double const& TheSep, double const& eTimeDay)
{
  double tolDay=double(1)/double(1000000);
  if (eTimeDay < FirstTime - tolDay) {
    std::cerr << "Error in GetTimeInterpolationInfo_infinite\n";
    std::cerr << "We have FirstTime = " << FirstTime << "\n";
    std::cerr << "     and eTimeDay = " << eTimeDay << "\n";
    std::cerr << "i.e. eTimeDay < FirstTime\n";
    throw TerminalException{1};
  }
  if (TheSep < 0) {
    std::cerr << "We need TheSep > 0\n";
    std::cerr << "But we have TheSep = " << TheSep << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "  1:  GetTimeInterpolationInfo_infinite\n";
  InterpInfo eInterpInfo;
  int iTime=1;
  while(true) {
    int iTimeUpp = iTime;
    int iTimeLow = iTime-1;
    double eTimeLow = FirstTime + double(iTimeLow)*TheSep;
    double eTimeUpp = FirstTime + double(iTimeUpp)*TheSep;
    if (fabs(eTimeLow - eTimeDay) < tolDay) {
      eInterpInfo.UseSingleEntry=true;
      eInterpInfo.iTimeLow=iTimeLow;
      eInterpInfo.iTimeUpp=iTimeLow;
      eInterpInfo.alphaLow=1; // Needed to avoid some bugs
      eInterpInfo.alphaUpp=0;
      return eInterpInfo;
    }
    if (fabs(eTimeUpp - eTimeDay) < tolDay) {
      eInterpInfo.UseSingleEntry=true;
      eInterpInfo.iTimeLow=iTimeUpp;
      eInterpInfo.iTimeUpp=iTimeUpp;
      eInterpInfo.alphaLow=1; // Needed to avoid some bugs
      eInterpInfo.alphaUpp=0;
      return eInterpInfo;
    }
    double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
    double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
    if (alphaLow >= 0 && alphaUpp >= 0) {
      eInterpInfo.iTimeLow = iTimeLow;
      eInterpInfo.iTimeUpp = iTimeUpp;
      eInterpInfo.alphaLow = alphaLow;
      eInterpInfo.alphaUpp = alphaUpp;
      eInterpInfo.UseSingleEntry=false;
      return eInterpInfo;
    }
    iTime++;
    if (iTime > 1000000) {
      std::cerr << "FirstTime = " << FirstTime << "\n";
      std::cerr << "eTimeFay  = " << eTimeDay << "\n";
      std::cerr << "iTime     = " << iTime << "\n";
      std::cerr << "TheSep    = " << TheSep << "\n";
      std::cerr << "Probably a bug in the infinite loop\n";
      throw TerminalException{1};
    }
  }
}

std::vector<int> GetIFileIRec_infinite(int const& nbRecBegin, int const& nbRecMiddle, int const& iTime)
{
  int iFile=0;
  int iRec;
  while(true) {
    if (iTime < nbRecBegin + iFile*nbRecMiddle) {
      if (iFile == 0) {
	iRec=iTime;
      }
      else {
	iRec = iTime - nbRecBegin - (iFile-1)*nbRecMiddle;
      }
      return {iFile, iRec};
    }
    iFile++;
  }
}

std::vector<int> GetIFileIRec(ArrayHistory const& eArr, int const& iTime)
{
  if (eArr.TimeSteppingInfo == "classic") {
    int iFile=eArr.ListIFile[iTime];
    int iRec=eArr.ListIRec[iTime];
    return {iFile, iRec};
  }
  if (eArr.TimeSteppingInfo == "singlefile") {
    int iFile=0;
    int iRec=iTime;
    return {iFile,iRec};
  }
  return GetIFileIRec_infinite(eArr.nbRecBegin, eArr.nbRecMiddle, iTime);
}


double ARR_GetTime(ArrayHistory const& eArr, int const& iTime)
{
  if (eArr.TimeSteppingInfo == "classic") {
    return eArr.ListTime[iTime];
  }
  return eArr.FirstTime + iTime*eArr.SeparationTime;
}


std::vector<int> GetIntervalListITime(ArrayHistory const& eArr, double const& eTimeDay, double const& TimeFrameDay)
{
  double epsilon=0.0000001;
  std::vector<int> ListRelITime;
  double eTimeLow=eTimeDay - epsilon;
  double eTimeUpp=eTimeDay + TimeFrameDay - epsilon;
  int iTime=0;
  while(true) {
    double eTime=ARR_GetTime(eArr, iTime);
    if (eTime > eTimeLow && eTime < eTimeUpp)
      ListRelITime.push_back(iTime);
    if (eTime > eTimeUpp)
      break;
    iTime++;
  }
  return ListRelITime;
}


double GetListTimeSeparation(std::vector<double> const& ListTime)
{
  int nbTime=ListTime.size();
  std::vector<double> ListVal;
  std::vector<int> ListNb;
  double tolDay = double(1) / double(100000);
  auto InsertDiff=[&](double const& eVal) -> void {
    int len=ListVal.size();
    for (int i=0; i<len; i++) {
      if (fabs(eVal - ListVal[i]) < tolDay) {
	ListNb[i]++;
	return;
      }
    }
    ListVal.push_back(eVal);
    ListNb.push_back(1);
  };
  for (int iTime=1; iTime<nbTime; iTime++) {
    double eDiff=ListTime[iTime] - ListTime[iTime-1];
    InsertDiff(eDiff);
  }
  int siz=ListVal.size();
  int eNb=0;
  double eVal = -1;
  for (int i=0; i<siz; i++) {
    if (ListNb[i] > eNb) {
      eNb = ListNb[i];
      eVal = ListVal[i];
    }
  }
  return eVal;
}


InterpInfo GetTimeInterpolationInfoGeneralized(ArrayHistory const& eArr, double const& eTimeDay)
{
  if (eArr.TimeSteppingInfo == "classic") {
    return GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
  }
  if (eArr.TimeSteppingInfo == "singlefile" || eArr.TimeSteppingInfo == "multiplenetcdf") {
    return GetTimeInterpolationInfo_infinite(eArr.FirstTime, eArr.SeparationTime, eTimeDay);
  }
  std::cerr << "Did not find a matching entry\n";
  throw TerminalException{1};
}



#endif
