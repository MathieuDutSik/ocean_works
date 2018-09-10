#ifndef MJDV2_INCLUDE
#define MJDV2_INCLUDE

#include "Temp_common.h"
#include "Basic_string.h"

struct VarQuery {
  double eTimeDay;
  int iTime;
  std::string NatureQuery; // Can be "instant", "average", "swathMax", "swathMin"
  double TimeFrameDay;
  std::string typeQuery;
};




std::string GetMonthName(int const& iMonth)
{
  if (iMonth <= 0 || iMonth > 12) {
    std::cerr << "We should have iMonth between 1 and 12\n";
    std::cerr << "iMonth=" << iMonth << "\n";
    throw TerminalException{1};
  }
  std::vector<std::string> ListStr{"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  return ListStr[iMonth-1];
}


int MONTH_LEN(int const& year, int const& month)
{
  if (month == 1 || month == 3 || month == 5 || month == 7 || month == 8 || month == 10 || month == 12)
    return 31;
  if (month == 4 || month == 6 || month == 9 || month == 11)
    return 30;
  if (month == 2) {
    int res4=year % 4;
    int res100=year % 100;
    int res400=year % 400;
    if (res4 != 0) {
      return 28;
    }
    else {
      if (res100 != 0) {
	return 29;
      }
      else {
	if (res400 != 0) {
	  return 28;
	}
	else {
	  return 29;
	}
      }
    }
  }
  std::cerr << "Error happened in LEN_MONTH\n";
  throw TerminalException{1};
}


bool TestCorrectnessVectorTime(std::vector<int> const& eDate)
{
  int year=eDate[0];
  int month=eDate[1];
  int day=eDate[2];
  int hour=eDate[3];
  int min=eDate[4];
  int sec=eDate[5];
  if (month < 1 || month > 12) {
    std::cerr << "We have month=" << month << " but possible values are between 1 and 12\n";
    return false;
  }
  int month_len=MONTH_LEN(year, month);
  if (day < 1 || day > month_len) {
    std::cerr << "We have day=" << day << " but possible values are between 1 and " << month_len << "\n";
    return false;
  }
  if (hour < 0 || hour >= 24) {
    std::cerr << "We have hour=" << hour << " but possible values are between 0 and 23\n";
    return false;
  }
  if (min < 0 || min >= 60) {
    std::cerr << "We have min=" << min << " but possible values are between 0 and 59\n";
    return false;
  }
  if (sec < 0 || sec >= 60) {
    std::cerr << "We have sec=" << sec << " but possible values are between 0 and 59\n";
    return false;
  }
  return true;
}

void CheckDateAndDieIfIncorrect(std::vector<int> const& eDate)
{
  bool test=TestCorrectnessVectorTime(eDate);
  if (!test) {
    std::cerr << "We have incoherent date input\n";
    for (auto & eVal : eDate)
      std::cerr << " " << eVal;
    std::cerr << "\n";
    throw TerminalException{1};
  }
}






double DATE2JD(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  double eJDbase, eFracDay, eJD;
  int year, month, day, hour, min, sec;
  year  =eDate[0];
  month =eDate[1];
  day   =eDate[2];
  hour  =eDate[3];
  min   =eDate[4];
  sec   =eDate[5];
  int a, y, m;
  a = int(floor((double(14) - double(month))/double(12)));
  y = year + 4800 - a;
  m = month + 12*a - 3;
  // For a date in the Gregorian calendar:
  eJDbase = double(day) 
    + double(floor((double(153)*double(m) + double(2))/double(5)))
    + double(y)*double(365)                                       
    + double(floor(double(y)/double(4)))                           
    - double(floor(double(y)/double(100)))                          
    + double(floor(double(y)/double(400))) - double(32045);
  eFracDay=(double(sec) +                                         
	    double(60)*double(min) +                               
	    double(3600)*(double(hour) - double(12))               
	    )/double(86400);
  eJD=eJDbase + eFracDay;
  return eJD;
}

double DATE_ConvertSix2mjd(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  double eJD1=DATE2JD(eDate);
  double eJD2=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD=eJD1-eJD2;
  return eMJD;
}


std::vector<int> DATE_ConvertString2six(std::string const& eTimeStr)
{
  if (eTimeStr.size() != 15) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "is not of the right length (should be exactly 16 characters)\n";
    std::cerr << "as in 20160120.000000\n";
    throw TerminalException{1};
  }
  std::string eYear, eMonth, eDay, eHour, eMin, eSec;
  eYear=eTimeStr.substr(0,4);
  eMonth=eTimeStr.substr(4,2);
  eDay=eTimeStr.substr(6,2);
  eHour=eTimeStr.substr(9,2);
  eMin=eTimeStr.substr(11,2);
  eSec=eTimeStr.substr(13,2);

  int year, month, day, hour, min, sec;
  std::istringstream(eYear) >> year;
  std::istringstream(eMonth) >> month;
  std::istringstream(eDay) >> day;
  std::istringstream(eHour) >> hour;
  std::istringstream(eMin) >> min;
  std::istringstream(eSec) >> sec;
  if (month <= 0 || month > 12) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "gives month=" << month << "\n";
    throw TerminalException{1};
  }
  if (day <= 0 || day > MONTH_LEN(year, month)) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "gives day=" << day << "\n";
    throw TerminalException{1};
  }
  if (hour < 0 || hour >= 24) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "gives hour=" << hour << "\n";
    throw TerminalException{1};
  }
  if (min < 0 || min >= 60) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "gives min=" << min << "\n";
    throw TerminalException{1};
  }
  if (sec < 0 || sec > 60) {
    std::cerr << "eTimeStr=" << eTimeStr << "\n";
    std::cerr << "gives sec=" << sec << "\n";
    throw TerminalException{1};
  }
  std::vector<int> Date={year, month, day, hour, min, sec};
  return Date;
}


std::string DATE_ConvertSix2dhmz(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int year, month, day, hour;
  year   =eDate[0];
  month  =eDate[1];
  day    =eDate[2];
  hour   =eDate[3];
  std::string eTimeStr=StringNumber(year, 4) + 
    StringNumber(month, 2) + 
    StringNumber(day, 2) + "_" + 
    StringNumber(hour, 2);
  return eTimeStr;
}




double DATE_ConvertString2mjd(std::string const& eTimeStr)
{
  std::vector<int> eDate = DATE_ConvertString2six(eTimeStr);
  return DATE_ConvertSix2mjd(eDate);
}


std::vector<int> DATE_ConvertSix2tfn(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int year  = eDate[0];
  int month = eDate[1];
  int day   = eDate[2];
  int hour  = eDate[3];
  int min   = eDate[4];
  int sec   = eDate[5];
  int tfn1  = year*10000 + month*100 + day;
  int tfn2  = hour*10000 +   min*100 + sec;
  return {tfn1, tfn2};
}




std::string DATE_ConvertSix2string(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int year, month, day, hour, min, sec;
  year  =eDate[0];
  month =eDate[1];
  day   =eDate[2];
  hour  =eDate[3];
  min   =eDate[4];
  sec   =eDate[5];
  std::string eTimeStr=StringNumber(year, 4) + 
    StringNumber(month, 2) + 
    StringNumber(day, 2) + "." + 
    StringNumber(hour, 2) + 
    StringNumber(min, 2) + 
    StringNumber(sec, 2);
  return eTimeStr;
}


std::string DATE_ConvertSix2string_style1(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int year  = eDate[0];
  int month = eDate[1];
  int day   = eDate[2];
  return IntToString(day) + "." + IntToString(month) + "." + IntToString(year);
}


std::string DATE_ConvertSix2string_style2(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int month = eDate[1];
  int day   = eDate[2];
  return IntToString(day) + " " + GetMonthName(month);
}

std::string DATE_ConvertSix2string_style3(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int year  = eDate[0];
  int month = eDate[1];
  return IntToString(year) + "/" + StringNumber(month,2);
}


std::string DATE_ConvertSix2string_style4(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  int month = eDate[1];
  int day   = eDate[2];
  return StringNumber(month,2) + "." + StringNumber(day,2);
}

std::string DATE_ConvertSix2string_style(std::vector<int> const& eDate, std::string const& style)
{
  CheckDateAndDieIfIncorrect(eDate);
  if (style == "style0")
    return DATE_ConvertSix2string(eDate);
  if (style == "style1")
    return DATE_ConvertSix2string_style1(eDate);
  if (style == "style2")
    return DATE_ConvertSix2string_style2(eDate);
  if (style == "style3")
    return DATE_ConvertSix2string_style3(eDate);
  if (style == "style4")
    return DATE_ConvertSix2string_style4(eDate);
  std::cerr << "Failed to find relevant entry for the time style=" << style << "\n";
  std::cerr << "Allowed styles = style0, style1, style2 and style3\n";
  throw TerminalException{1};
}






std::string DATE_ConvertSix2mystringPres(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  try {
    int year  = eDate[0];
    int month = eDate[1];
    int day   = eDate[2];
    int hour  = eDate[3];
    int min   = eDate[4];
    int sec   = eDate[5];
    std::string eTimeStr=StringNumber(year, 4) + "-" + 
      StringNumber(month, 2) + "-" +
      StringNumber(day, 2) + " " + 
      StringNumber(hour, 2) + ":" + StringNumber(min, 2) + ":" + StringNumber(sec, 2);
    return eTimeStr;
  }
  catch (std::string & eStr) {
    std::stringstream s;
    s << "Error in DATE_ConvertSix2mystringFile\n";
    s << "eDate.size()=" << eDate.size() << "\n";
    s << "eDate=";
    WriteStdVector(s, eDate);
    s << "-----------------------------------------\n";
    s << "exception eStr=\n";
    s << eStr;
    std::string vStr(s.str());
    std::cerr << vStr;
    throw TerminalException{1};
  }
}

std::string DATE_ConvertSix2mystringPresReduced(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  if (eDate[3] != 0 || eDate[4] != 0 || eDate[5] != 0) {
    return DATE_ConvertSix2mystringPres(eDate);
  }
  try {
    int year  = eDate[0];
    int month = eDate[1];
    int day   = eDate[2];
    std::string eTimeStr=StringNumber(year, 4) + "-" + 
      StringNumber(month, 2) + "-" +
      StringNumber(day, 2);
    return eTimeStr;
  }
  catch (std::string & eStr) {
    std::stringstream s;
    s << "Error in DATE_ConvertSix2mystringFile\n";
    s << "eDate.size()=" << eDate.size() << "\n";
    s << "eDate=";
    WriteStdVector(s, eDate);
    s << "-----------------------------------------\n";
    s << "exception eStr=\n";
    s << eStr;
    std::string vStr(s.str());
    std::cerr << vStr;
    throw TerminalException{1};
  }
}






std::string DATE_ConvertSix2mystringFile(std::vector<int> const& eDate)
{
  CheckDateAndDieIfIncorrect(eDate);
  try {
    int year  = eDate[0];
    int month = eDate[1];
    int day   = eDate[2];
    int hour  = eDate[3];
    int min   = eDate[4];
    int sec   = eDate[5];
    std::string eTimeStr=StringNumber(year, 4) + 
      StringNumber(month, 2) + StringNumber(day, 2) + "_" +
      StringNumber(hour, 2) + StringNumber(min, 2) + StringNumber(sec, 2);
    return eTimeStr;
  }
  catch (std::string & eStr) {
    std::stringstream s;
    s << "Error in DATE_ConvertSix2mystringFile\n";
    s << "eDate.size()=" << eDate.size() << "\n";
    s << "eDate=";
    WriteStdVector(s, eDate);
    s << "-----------------------------------------\n";
    s << "exception eStr=\n";
    s << eStr;
    std::string vStr(s.str());
    std::cerr << vStr;
    throw TerminalException{1};
  }
}




// The following algorithm is from the Calendar FAQ.
std::vector<int> JD2DATE(double const& eJD)
{
  int year, month, day, hour, min, sec;
  int ijd, a, b, c, d, e, m;
  double fjd, second;
  ijd = int(floor(eJD + 0.5));
  //
  a = ijd + 32044;
  b = int(floor((double(4)*double(a) + double(3)) / double(146097)));
  c = a - int(floor((double(b) * double(146097)) / double(4)));
  //
  d = int(floor((double(4)*double(c) + double(3)) / double(1461)));
  e = c - int(floor((double(1461)*double(d)) / double(4)));
  m = int(floor((double(5) * double(e) + double(2)) / double(153)));
  //
  day   = e - int(floor((double(153) * double(m) + double(2)) / double(5))) + 1;
  month = m + 3 - 12 * int(floor(double(m) / double(10)));
  year  = b * 100 + d - 4800 + int(floor(double(m) / double(10)));
  //
  fjd    = eJD - double(ijd) + 0.5;
  second = double(86400) * fjd;
  hour   = int(floor(second/double(3600)));
  second = second - double(3600)*double(hour);
  min    = int(floor(second/double(60)));
  sec    = int(floor(second - double(60)*min));
  // Now renormalizing
  int secNear=int(round(second - double(60)*min));
  if (secNear == 60) {
    sec=0;
    min=min+1;
  }
  if (min == 60) {
    min=0;
    hour=hour+1;
  }
  if (hour == 24) {
    hour=0;
    day=day+1;
  }
  int lenmonth=MONTH_LEN(year, month);
  if (day == lenmonth+1) {
    day=1;
    month=month+1;
  }
  if (month == 13) {
    month=1;
    year=year+1;
  }
  std::vector<int> Date={year, month, day, hour, min, sec};
  return Date;
}

double CT2MJD(std::string const& STIME)
{
  double XMJD;
  std::vector<int> eDate=DATE_ConvertString2six(STIME);
  XMJD=DATE_ConvertSix2mjd(eDate);
  return XMJD;
}

std::string MJD2CT(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  return DATE_ConvertSix2string(eDate);
}

std::string DATE_ConvertMjd2mystringPres(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  return DATE_ConvertSix2mystringPres(eDate);
}


std::string DATE_ConvertMjd2mystringPresReduced(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  return DATE_ConvertSix2mystringPresReduced(eDate);
}


std::string DATE_ConvertMjd2mystringPresReducedMilisecond(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  //
  double Time_PresSec_Day = DATE2JD(eDate) - XMJD_1858;
  int Delta_Msec1 = 1000 * 86400 * (XMJD - Time_PresSec_Day);
  int Delta_Msec2 = std::max(std::min(999, Delta_Msec1), 0);
  //
  if (eDate[3] == 0 && eDate[4] == 0 && eDate[5] == 0 && Delta_Msec2 == 0) {
    int year  = eDate[0];
    int month = eDate[1];
    int day   = eDate[2];
    return StringNumber(year, 4) + "-" + 
      StringNumber(month, 2) + "-" +
      StringNumber(day, 2);
  }
  std::string strRet=DATE_ConvertSix2mystringPres(eDate);
  if (Delta_Msec2 > 0) {
    strRet += " " + StringNumber(Delta_Msec2,3);
  }
  return strRet;
}







std::vector<int> DATE_ConvertMjd2six(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  return JD2DATE(eMJD);
}

std::string DATE_ConvertMjd2dhmz(double const& eMJD)
{
  std::vector<int> eDate=DATE_ConvertMjd2six(eMJD);
  return DATE_ConvertSix2dhmz(eDate);
}

std::vector<int> DATE_ConvertMjd2tfn(double const& eTimeDay)
{
  std::vector<int> Date=DATE_ConvertMjd2six(eTimeDay);
  return DATE_ConvertSix2tfn(Date);
}



std::string DATE_ConvertMjd2mystringFile(double const& XMJD)
{
  double XMJD_1858, eMJD;
  XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  return DATE_ConvertSix2mystringFile(eDate);
}

std::string DATE_ConvertMjd2mystringFileMilisecond(double const& XMJD)
{
  double XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  //
  double Time_PresSec_Day = DATE2JD(eDate) - XMJD_1858;
  int Delta_Msec1 = 1000 * 86400 * (XMJD - Time_PresSec_Day);
  int Delta_Msec2 = std::max(std::min(999, Delta_Msec1), 0);
  std::cerr << "Delta_Msec2=" << Delta_Msec2 << " XMJD=" << XMJD << " Time_PresSec_Day=" << Time_PresSec_Day << "\n";
  std::string strRet = DATE_ConvertSix2mystringFile(eDate) + "_" + StringNumber(Delta_Msec2,3);
  return strRet;
}







std::vector<double> GetIntervalFLD(double const& FirstTime, double const& LastTime, double const& DeltaInterval)
{
  double eTime=FirstTime;
  std::vector<double> ListTime;
  double tolDay=DeltaInterval / double(100);
  //  std::cerr.width(15);
  //  std::cerr << "DeltaInterval=" << DeltaInterval << "\n";
  while(true) {
    ListTime.push_back(eTime);
    eTime=eTime + DeltaInterval;
    if (eTime > LastTime + tolDay)
      return ListTime;
  }
}

std::vector<VarQuery> MonSeas_Kernel_GetIntervalFL(double const& FirstTime, double const& LastTime, int const& inc, std::string const& typeQuery)
{
  std::vector<int> eVec1=DATE_ConvertMjd2six(FirstTime);
  int year_first=eVec1[0];
  std::vector<int> eVec2=DATE_ConvertMjd2six(LastTime);
  int year_last=eVec2[0];
  std::vector<VarQuery> ListQuery;
  double eps=0.00001;
  int iTime=0;
  for (int iYear=year_first; iYear<=year_last; iYear++) {
    for (int iMon=1; iMon<=12; iMon += inc) {
      std::vector<int> eVec{iYear, iMon, 1, 0, 0, 0};
      double eMJD=DATE_ConvertSix2mjd(eVec);
      if (eMJD >= FirstTime - eps && eMJD <= LastTime + eps) {
	int iYear_next=iYear, iMon_next=iMon;
	for (int iter=0; iter<inc; iter++) {
	  if (iMon_next == 12) {
	    iYear_next=iYear_next+1;
	    iMon_next=1;
	  }
	  else {
	    //	    iYear_next=iYear_next;
	    iMon_next=iMon_next + 1;
	  }
	}
	std::vector<int> eVecNext{iYear_next, iMon_next, 1, 0, 0, 0};
	double eMJD_next=DATE_ConvertSix2mjd(eVecNext);
	if (eMJD_next <= LastTime + eps) {
	  double TimeFrameDay = eMJD_next - eMJD;
	  VarQuery eQuery;
	  eQuery.eTimeDay=eMJD;
	  eQuery.iTime=iTime;
	  eQuery.TimeFrameDay=TimeFrameDay;
	  eQuery.typeQuery=typeQuery;
	  ListQuery.push_back(eQuery);
	  iTime++;
	}
      }
    }
  }
  return ListQuery;
}


std::vector<VarQuery> MonSeas_Kernel_GetIntervalFLseasonalIvica(double const& FirstTime, double const& LastTime)
{
  std::vector<int> eVec1=DATE_ConvertMjd2six(FirstTime);
  int year_first=eVec1[0];
  std::vector<int> eVec2=DATE_ConvertMjd2six(LastTime);
  int year_last=eVec2[0];
  std::vector<VarQuery> ListQuery;
  double eps=0.00001;
  int iTime=0;
  std::vector<int> ListRelMon{12,3,6,9};
  for (int iYear=year_first; iYear<=year_last; iYear++) {
    for (int i=0; i<4; i++) {
      int iMon=ListRelMon[i];
      std::vector<int> eVec{iYear, iMon, 1, 0, 0, 0};
      double eMJD=DATE_ConvertSix2mjd(eVec);
      if (eMJD >= FirstTime - eps && eMJD <= LastTime + eps) {
	int iYear_next=iYear, iMon_next=iMon;
	for (int iter=0; iter<3; iter++) {
	  if (iMon_next == 12) {
	    iYear_next=iYear_next+1;
	    iMon_next=1;
	  }
	  else {
	    iMon_next=iMon_next + 1;
	  }
	}
	std::vector<int> eVecNext{iYear_next, iMon_next, 1, 0, 0, 0};
	double eMJD_next=DATE_ConvertSix2mjd(eVecNext);
	if (eMJD_next <= LastTime + eps) {
	  double TimeFrameDay = eMJD_next - eMJD;
	  VarQuery eQuery;
	  eQuery.eTimeDay=eMJD;
	  eQuery.iTime=iTime;
	  eQuery.TimeFrameDay=TimeFrameDay;
	  eQuery.typeQuery="seasonalIvica";
	  ListQuery.push_back(eQuery);
	  iTime++;
	}
      }
    }
  }
  return ListQuery;
}




std::vector<VarQuery> GetIntervalFLmonthly(double const& FirstTime, double const& LastTime)
{
  return MonSeas_Kernel_GetIntervalFL(FirstTime, LastTime, 1, "monthly");
}



std::vector<VarQuery> GetIntervalFLseasonal(double const& FirstTime, double const& LastTime)
{
  return MonSeas_Kernel_GetIntervalFL(FirstTime, LastTime, 3, "seasonal");
}


std::vector<VarQuery> GetIntervalFLseasonalIvica(double const& FirstTime, double const& LastTime)
{
  return MonSeas_Kernel_GetIntervalFLseasonalIvica(FirstTime, LastTime);
}



double GetIntervalSize(double const& DELTC, std::string const& UNITC)
{
  int IsDone=0;
  double eMult = 0;
  if (UNITC == "DAY") {
    IsDone=1;
    eMult=double(1);
  }
  if (UNITC == "HOUR") {
    IsDone=1;
    eMult=double(1)/double(24);
  }
  if (UNITC == "HR") {
    IsDone=1;
    eMult=double(1)/double(24);
  }
  if (UNITC == "MIN") {
    IsDone=1;
    eMult=double(1)/double(1440);
  }
  if (UNITC == "SEC") {
    IsDone=1;
    eMult=double(1)/double(86400);
  }
  if (IsDone == 0) {
    std::cerr << "UNITC has not been found\n";
    std::cerr << "Allowed: DAY, HR, HOUR, MIN, SEC\n";
    std::cerr << "UNITC=" << UNITC << "\n";
    throw TerminalException{1};
  }
  double DeltaInterval = DELTC * eMult;
  return DeltaInterval;
}



std::vector<double> GetInterval(std::string const& BEGTC, std::string const& ENDTC, double const& DELTC, std::string const& UNITC)
{
  double DeltaInterval=GetIntervalSize(DELTC, UNITC);
  double FirstTime=CT2MJD(BEGTC);
  double LastTime=CT2MJD(ENDTC);
  double tolDay= DeltaInterval / double(10000);
  if (LastTime < FirstTime - tolDay) {
    std::cerr << "We should have ENDTC >= BEGTC. But instead we have:\n";
    std::cerr << "BEGTC = " << BEGTC << "\n";
    std::cerr << "ENDTC = " << ENDTC << "\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  return GetIntervalFLD(FirstTime, LastTime, DeltaInterval);
}




#endif
