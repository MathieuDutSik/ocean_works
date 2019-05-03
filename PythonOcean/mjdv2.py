import os
import math
import decimal

def DATE2JD(Date):
    import math
    year=Date[0]
    month=Date[1]
    day=Date[2]
    hour=Date[3]
    minD=Date[4]
    sec=Date[5]
    a = math.floor((14. - month)/12.)
    y = year + 4800 - a
    m = month + 12*a - 3
    eJDbase = day + math.floor((153.*m + 2.)/5.) + y*365. + math.floor(y/4.) - math.floor(y/100.) + math.floor(y/400.) - 32045.
    eFracDay=(sec + 60 * minD + 3600 * (hour - 12))/86400.
    eJD=eJDbase + eFracDay
    return eJD



def DATE_ConvertSix2mjd(eDate):
    eJD1=DATE2JD(eDate)
    eJD2=DATE2JD([1858, 11, 17, 0, 0, 0])
    eMJD=eJD1-eJD2
    return eMJD


def DATE_ConvertSix2mjdROMS(eDate):
    eJD1=DATE2JD(eDate)
    eJD2=DATE2JD([1968, 5, 23, 0, 0, 0])
    eMJD=eJD1-eJD2
    return eMJD


def DATE_ConvertString2six(eTimeStr):
    eYear=int(eTimeStr[0:3])
    eMonth=int(eTimeStr[4:5])
    eDay=int(eTimeStr[6:7])
    eHour=int(eTimeStr[9:10])
    eMin=int(eTimeStr[11:12])
    eSec=int(eTimeStr[13:14])
    return [eYear, eMonth, eDay, eHour, eMin, eSec]

def DATE_ConvertString2mjd(eTimeStr):
    eDate = DATE_ConvertString2six(eTimeStr)
    return DATE_ConvertSix2mjd(eDate)



def StringNumber(eNb, nbDigit):
    if (eNb < 10):
        return "0" + str(eNb)
    if (eNb < 100):
        return str(eNb)
    print("eNb=", eNb)
    print("Error in StringNumber\n")
    os.sys.exit()
                            
    
    





def DATE_ConvertSix2string(Date):
    year=Date[0]
    month=Date[1]
    day=Date[2]
    hour=Date[3]
    minD=Date[4]
    sec=Date[5]
    str1 = str(year) + StringNumber(month,2) + StringNumber(day,2)
    str2 = StringNumber(hour,2) + StringNumber(minD,2) + StringNumber(sec,2)
    return str1 + "." + str2

def MONTH_LEN(year, month):
    if (month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12):
        return 31;
    if (month == 4 or month == 6 or month == 9 or month == 11):
        return 30;
    if (month == 2):
        res4=year % 4;
        res100=year % 100;
        res400=year % 400;
        if (res4 != 0):
            return 28;
        if (res100 != 0):
            return 29;
        if (res400 != 0):
            return 28;
        return 29;
    print("year=", year, " month=", month)
    print("Error in MONTH_LEN\n");
    os.sys.exit()

def NEXT_DAY(year, month, day):
    lenMonth=MONTH_LEN(year,month);
    if (day < lenMonth):
        return [year, month, day+1];
    if (month < 12):
        return [year, month+1, 1];
    return [year+1, 1, 1];


def NEXT_MONTH(year, month):
    if (month < 12):
        return [year, month+1];
    return [year+1, 1];

def PREV_MONTH(year, month):
    if (month > 1):
        return [year, month-1];
    return [year-1, 12];



def StringOfDay(year, month, day):
    if (month < 10):
        strMonth = "0" + str(month);
    else:
        strMonth = str(month);
    if (day < 10):
        strDay = "0" + str(day);
    else:
        strDay = str(day);
    return str(year) + "-" + strMonth + "-" + strDay;



    
def JD2DATE(eJD):
    ijd = math.floor(eJD + 0.5)
    a = ijd + 32044
    b = math.floor((4.*a + 3.) / 146097.)
    c = a - math.floor((b * 146097.) / 4.)
    d = math.floor((4.*c + 3.) / 1461.)
    e = c - math.floor((1461.*d) / 4.)
    m = math.floor((5. * e + 2.) / 153.)
    day   = e - math.floor((153. * m + 2.) / 5.) + 1;
    month = m + 3 - 12 * math.floor(m / 10.);
    year  = b * 100 + d - 4800 + math.floor(m / 10.);
    #
    fjd    = eJD - ijd + 0.5;
    second = 86400. * fjd;
    hour   = math.floor(second/3600.);
    second = second - 3600.*hour;
    minD   = math.floor(second/60.);
    sec    = math.floor(second - 60.*minD);
    #
#    secNear = math.round(second - 60.*minD)
    secNear = int(round(second - 60.*minD))
    if (secNear == 0):
        sec=0
        minD=minD+1
    if (minD == 60):
        minD=0
        hour=hour+1
    if (hour == 24):
        hour=0
        day=day+1
    lenmonth=MONTH_LEN(year,month)
    if (day == lenmonth+1):
        day=1
        month=month+1
    if (month == 13):
        month=1
        year=year+1
    return [year, month, day, hour, minD, sec]

def CT2MJD(STIME):
    eDate=DATE_ConvertString2six(STIME)
    return DATE_ConvertSix2mjd(eDate)

def MJD2CT(XMJD):
    XMJD_1858=DATE2JD([1858, 11, 17, 0, 0, 0])
    eMJD = XMJD + XMJD_1858
    eDate=JD2DATE(eMJD)
    return DATE_ConvertSix2string(eDate);

def DATE_ConvertSix2mystringPres(eDate):
    eYear=eDate[0]
    eMonth=eDate[1]
    eDay=eDate[2]
    eHour=eDate[3]
    eMin=eDate[4]
    eSec=eDate[5]
    strYear = str(eYear)
    strMonth = StringNumber(eMonth, 2)
    strDay = StringNumber(eDay, 2)
    strHour = StringNumber(eHour, 2)
    strMin = StringNumber(eMin, 2)
    strSec = StringNumber(eSec, 2)
    return strYear + "-" + strMonth + "-" + strDay + " " + strHour + ":" + strMin + ":" + strSec
    

def DATE_ConvertMjd2mystringPres(XMJD):
    XMJD_1858=DATE2JD([1858, 11, 17, 0, 0, 0])
    eMJD = XMJD + XMJD_1858
    eDate=JD2DATE(eMJD)
    return DATE_ConvertSix2mystringPres(eDate)


def DATE_ConvertMjd2six(XMJD):
    XMJD_1858=DATE2JD([1858, 11, 17, 0, 0, 0]);
    eMJD = XMJD + XMJD_1858;
    return JD2DATE(eMJD);



                
    
