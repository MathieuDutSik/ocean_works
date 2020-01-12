#ifndef BASIC_NETCDF_INCLUDE
#define BASIC_NETCDF_INCLUDE

#include "MAT_Matrix.h"
#include "MAT_Tensor.h"
#include "Temp_common.h"
#include "Basic_Ocean_types.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "mjdv2.h"
#include <netcdf>

struct CFinformation {
  std::string LongName;
  std::string Units;
  std::string StdName;
};



double GetStandardMissingValue()
{
  return double(-100000000);
}


CFinformation GetCFnames(std::string const& var)
{
  if (var == "lon")
    return {"longitude", "degrees_east", "longitude"};
  if (var == "lat")
    return {"latitude", "degrees_north", "latitude"};
  if (var == "ang")
    return {"grid_Angle", "degrees", "grid_angle"};
  if (var == "mask")
    return {"mask", "integer", "grid_mask"};
  std::cerr << "Did not find matching name\n";
  std::cerr << "var = " << var << "\n";
  throw TerminalException{1};
}







bool NC_IsVar(std::string const& eFile, std::string const& eVar)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_IsVar\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << " eVar=" << eVar << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "eFile = " << eFile << "\n";
  try {
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    //    std::cerr << "After the dataFile creation\n";
    netCDF::NcVar data=dataFile.getVar(eVar);
    if(data.isNull()) {
      return false;
    }
    return true;
  }
  catch (...) {
    return false;
  }
}

/*
  This appears impossible in practice to write.
  The netcdf functionality does not allow to download the full list
  of variable names.
  Design error. */
/*
std::vector<std::string> NC_ListVar(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ListVar\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << "\n";
    throw TerminalException{1};
  }
  std::vector<std::string> ListVarName;
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  std::multimap<std::string,netCDF::NcVar> TotalList = dataFile.getVars(netCDF::NcGroup::ChildrenGrps);

  
  return {};
  }*/


struct RecTime {
  netCDF::NcDim timeDim;
  netCDF::NcVar timeVarSec;
  netCDF::NcVar timeVarDay;
  netCDF::NcVar timeVarStr;
  std::string strTime;
};

RecTime AddTimeArray(netCDF::NcFile & dataFile, std::string const& strTime, double const& RefTime)
{
  netCDF::NcDim timeDim=dataFile.addDim(strTime);
  std::vector<std::string> LDim1{strTime};
  std::vector<std::string> LDim2{strTime, "dateString"};
  netCDF::NcVar timeVarSec=dataFile.addVar(strTime, "double", LDim1);
  netCDF::NcVar timeVarDay=dataFile.addVar(strTime + "_day", "double", LDim1);
  netCDF::NcVar timeVarStr=dataFile.addVar(strTime + "_str", "char", LDim2);
  std::string dateStr=DATE_ConvertMjd2mystringPres(RefTime);
  std::string attStr1="seconds since " + dateStr;
  std::string attStr2="days since " + dateStr;
  timeVarSec.putAtt("long_name", "time");
  timeVarDay.putAtt("long_name", "time");
  timeVarSec.putAtt("units", attStr1);
  timeVarDay.putAtt("units", attStr2);
  timeVarSec.putAtt("calendar", "gregorian");
  timeVarDay.putAtt("calendar", "gregorian");
  return {timeDim, timeVarSec, timeVarDay, timeVarStr, strTime};
}



void AddTimeArrayRomsBound(netCDF::NcFile & dataFile, std::string const& strTime, double const& RefTime)
{
  netCDF::NcDim timeDim=dataFile.addDim(strTime);
  std::vector<std::string> LDim1{strTime};
  std::vector<std::string> LDim2{strTime, "dateString"};
  netCDF::NcVar timeVarDay=dataFile.addVar(strTime, "double", LDim1);
  netCDF::NcVar timeVarStr=dataFile.addVar(strTime + "_str", "char", LDim2);
  std::string dateStr=DATE_ConvertMjd2mystringPres(RefTime);
  std::string attStr="days since " + dateStr;
  timeVarDay.putAtt("long_name", "time");
  timeVarDay.putAtt("units", attStr);
  timeVarDay.putAtt("calendar", "gregorian");
}


void AddTimeArrayROMS(netCDF::NcFile & dataFile, std::string const& strTime, double const& RefTime)
{
  netCDF::NcDim timeDim=dataFile.addDim(strTime);
  std::vector<std::string> LDim1{strTime};
  std::vector<std::string> LDim2{strTime, "dateString"};
  netCDF::NcVar timeVarDay=dataFile.addVar(strTime, "double", LDim1);
  netCDF::NcVar timeVarStr=dataFile.addVar(strTime + "_str", "char", LDim2);
  std::string dateStr=DATE_ConvertMjd2mystringPres(RefTime);
  std::string attStr1="days since " + dateStr;
  timeVarDay.putAtt("long_name", "time");
  timeVarDay.putAtt("units", attStr1);
  timeVarDay.putAtt("calendar", "gregorian");
  netCDF::NcVar timeVarDayOT1=dataFile.getVar("ocean_time");
  if (timeVarDayOT1.isNull() ) {
    netCDF::NcVar timeVarDayOT2=dataFile.addVar("ocean_time", "double", LDim1);
    timeVarDayOT2.putAtt("long_name", "time");
    timeVarDayOT2.putAtt("units", attStr1);
    timeVarDayOT2.putAtt("calendar", "gregorian");
  }
}





void PutTimeDay(RecTime & eRec, size_t const& pos, double const& eTimeDay)
{
  double eTimeSec=eTimeDay * double(86400);
  std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
  //  std::cerr << "strPres=" << strPres << "\n";
  //  std::cerr << "pos=" << pos << "\n";
  std::vector<size_t> start2{size_t(pos)};
  std::vector<size_t> count2{1};
  eRec.timeVarSec.putVar(start2, count2, &eTimeSec);
  eRec.timeVarDay.putVar(start2, count2, &eTimeDay);
  std::vector<size_t> start3{size_t(pos),0};
  std::vector<size_t> count3{1,19};
  eRec.timeVarStr.putVar(start3, count3, strPres.c_str());
}




MyVector<int> NC_ReadVariable_StatusFill_data(netCDF::NcVar const& data)
{
  if (data.isNull()) {
    std::cerr << "NC_ReadVariable_StatusFill_data : The array data is null\n";
    throw TerminalException{1};
  }
  netCDF::NcType eType=data.getType();
  if (eType.isNull()) {
    std::cerr << "NC_ReadVariable_StatusFill_data : eType is null is an error\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  int nbTot=1;
  for (int iDim=0; iDim<nbDim; iDim++) {
    netCDF::NcDim eDim=data.getDim(iDim);
    int valDim=eDim.getSize();
    nbTot *= valDim;
  }
  //  std::cerr << "nbTot=" << nbTot << "\n";
  bool IsMatch=false;
  MyVector<int> StatusFill=ZeroVector<int>(nbTot);
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double eFillValue;
    netCDF::NcVarAtt eAtt=data.getAtt("_FillValue");
    eAtt.getValues(&eFillValue);
    //
    double *eVal;
    eVal=new double[nbTot];
    data.getVar(eVal);
    for (int i=0; i<nbTot; i++)
      if (eVal[i] == eFillValue)
	StatusFill(i)=1;
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float eFillValue;
    netCDF::NcVarAtt eAtt=data.getAtt("_FillValue");
    eAtt.getValues(&eFillValue);
    //
    float *eVal;
    eVal=new float[nbTot];
    data.getVar(eVal);
    for (int i=0; i<nbTot; i++)
      if (eVal[i] == eFillValue)
	StatusFill(i)=1;
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_INT) {
    int eFillValue;
    netCDF::NcVarAtt eAtt=data.getAtt("_FillValue");
    eAtt.getValues(&eFillValue);
    //
    int *eVal;
    eVal=new int[nbTot];
    data.getVar(eVal);
    for (int i=0; i<nbTot; i++)
      if (eVal[i] == eFillValue)
	StatusFill(i)=1;
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_SHORT) {
    signed short int eFillValue;
    netCDF::NcVarAtt eAtt=data.getAtt("_FillValue");
    eAtt.getValues(&eFillValue);
    //
    signed short int *eVal;
    eVal=new signed short int[nbTot];
    data.getVar(eVal);
    for (int i=0; i<nbTot; i++)
      if (eVal[i] == eFillValue)
	StatusFill(i)=1;
    delete [] eVal;
    IsMatch=true;
  }
  if (!IsMatch) {
    std::cerr << "NC_ReadVariable_StatusFill_data : Did not find the right number type\n";
    throw TerminalException{1};
  }
  return StatusFill;
}

std::vector<size_t> NC_ReadVariable_listdim(netCDF::NcVar const& data)
{
  if (data.isNull()) {
    std::cerr << "Error in NC_ReadVariable_listdim : The array data is null\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  std::vector<size_t> ListDim(nbDim);
  for (int i=0; i<nbDim; i++) {
    netCDF::NcDim eDim=data.getDim(i);
    size_t eDim_i=eDim.getSize();
    ListDim[i]=eDim_i;
  }
  return ListDim;
}


int NC_ReadVariable_NbFillValue_data(netCDF::NcVar const& data)
{
  MyVector<int> StatusFill=NC_ReadVariable_StatusFill_data(data);
  int len=StatusFill.size();
  int nbFillValue=0;
  for (int i=0; i<len; i++)
    nbFillValue += StatusFill(i);
  return nbFillValue;
}

// We cannot return a netCDF::NcVar out of scope
// because it apparently depends on netCDF::NcFile which would
// go out of scope.
void CheckNetcdfDataArray(std::string const& CallFct, std::string const& eFile, std::string const& eVar)
{
  //  std::cerr << "Beginning of CheckNetcdfDataArray\n";
  try {
    //    std::cerr << "CallFct=" << CallFct << " eFile = " << eFile << " eVar=" << eVar << " step 1\n";
    if (!IsExistingFile(eFile)) {
      std::cerr << "Error in CheckNetcdfDataArray\n";
      std::cerr << "Called from CallFct=" << CallFct << "\n";
      std::cerr << "Trying to open non-existing file\n";
      std::cerr << "eFile = " << eFile << "\n";
      throw TerminalException{1};
    }
    //    std::cerr << "CallFct=" << CallFct << " eFile = " << eFile << " eVar=" << eVar << " step 2\n";
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    //    std::cerr << "CallFct=" << CallFct << " eFile = " << eFile << " eVar=" << eVar << " step 3\n";
    if (dataFile.isNull()) {
      std::cerr << "Error in CheckNetcdfDataArray : we found dataFile to be null\n";
      std::cerr << "Called from CallFct=" << CallFct << "\n";
      throw TerminalException{1};
    }
    //    std::cerr << "CallFct=" << CallFct << " eFile = " << eFile << " eVar=" << eVar << " step 4\n";
    netCDF::NcVar data=dataFile.getVar(eVar);
    //    std::cerr << "CallFct=" << CallFct << " eFile=" << eFile << " eVar=" << eVar << " step 5\n";
    if (data.isNull()) {
      std::cerr << "Error in CheckNetcdfDataArray. Variable data is null\n";
      std::cerr << "Likely the variable is absent from the netcdf file\n";
      std::cerr << "Called from CallFct = " << CallFct << "\n";
      std::cerr << "eFile = " << eFile << "\n";
      std::cerr << "eVar  = " << eVar << "\n";
      throw TerminalException{1};
    }
    //    std::cerr << "CallFct=" << CallFct << " eFile=" << eFile << " eVar=" << eVar << " step 6\n";
  }
  catch(...) {
    std::cerr << "Catch an exception in trying to read file\n";
    throw TerminalException{1};
  }
}



MyMatrix<int> NC_Read2Dvariable_Mask_data(netCDF::NcVar const& data)
{
  if (data.isNull()) {
    std::cerr << "NC_Read2Dvariable_Mask_data : The array data is null\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2) {
    std::cerr << "NC_Read2Dvariable_Mask_data : The number of dimensions is not correct\n";
    std::cerr << "nbDim=" << nbDim << " instead of 2\n";
    throw TerminalException{1};
  }
  netCDF::NcDim eDim=data.getDim(0);
  int eta=eDim.getSize();
  netCDF::NcDim fDim=data.getDim(1);
  int xi=fDim.getSize();
  MyVector<int> StatusFill=NC_ReadVariable_StatusFill_data(data);
  MyMatrix<int> eArr(eta, xi);
  int idx=0;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      eArr(i,j)=StatusFill(idx);
      idx++;
    }
  return eArr;
}




Eigen::Tensor<int,3> NC_Read3Dvariable_Mask_data(netCDF::NcVar const& data)
{
  if (data.isNull()) {
    std::cerr << "NC_Read2Dvariable_Mask_data : The array data is null\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 3) {
    std::cerr << "NC_Read2Dvariable_Mask_data : The number of dimensions is not correct\n";
    std::cerr << "nbDim=" << nbDim << " instead of 3\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int dim0=ListDim[0];
  int dim1=ListDim[1];
  int dim2=ListDim[2];
  MyVector<int> StatusFill=NC_ReadVariable_StatusFill_data(data);
  Eigen::Tensor<int,3> eTens(dim0, dim1, dim2);
  int idx=0;
  for (int i0=0; i0<dim0; i0++)
    for (int i1=0; i1<dim1; i1++)
      for (int i2=0; i2<dim2; i2++) {
	eTens(i0, i1, i2)=StatusFill(idx);
	idx++;
      }
  return eTens;
}

Eigen::Tensor<int,3> NC_Read3Dvariable_Mask_file(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read3Dvariable_Mask_file", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read3Dvariable_Mask_data(data);
}



MyMatrix<int> StrictProjectionMask(Eigen::Tensor<int,3> const& eTens, int eDim)
{
  if (eDim == 0) {
    auto LDim=eTens.dimensions();
    int dim0=LDim[0];
    MyMatrix<int> eProj(LDim[1], LDim[2]);
    for (int i1=0; i1<LDim[1]; i1++)
      for (int i2=0; i2<LDim[2]; i2++) {
	int sum=0;
	for (int i0=0; i0<dim0; i0++)
	  sum += eTens(i0, i1, i2);
	if (sum != 0 && sum != dim0) {
	  std::cerr << "We should have sum=0 or dim0\n";
	  throw TerminalException{1};
	}
	eProj(i1, i2) = sum / dim0;
      }
    return eProj;
  }
  std::cerr << "Code missing or incorrect value. eDim=" << eDim << "\n";
  throw TerminalException{1};
}




int NC_ReadDimension(netCDF::NcFile const& dataFile, std::string const& dimName)
{
  netCDF::NcDim eDim = dataFile.getDim(dimName);
  return eDim.getSize();
}







MyVector<double> NC_ReadVariable_data_start_count(netCDF::NcVar const& data, std::vector<size_t> const& start, std::vector<size_t> const& count)
{
  if (data.isNull()) {
    std::cerr << "NC_ReadVariable_data_start_count : the variable data satisfies isNull\n";
    throw TerminalException{1};
  }
  netCDF::NcType eType=data.getType();
  if (eType.isNull()) {
    std::cerr << "NC_ReadVariable_data_start_count : eType is null is an error\n";
    throw TerminalException{1};
  }
  int eDimTot=1;
  for (auto & eVal : count)
    eDimTot *= eVal;
  MyVector<double> eArr(eDimTot);
  bool IsMatch=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[eDimTot];
    data.getVar(start, count, eVal);
    for (int i=0; i<eDimTot; i++)
      eArr(i)=eVal[i];
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eValFLOAT;
    eValFLOAT=new float[eDimTot];
    //    std::cerr << "Before eValFLOAT read\n";
    data.getVar(start, count, eValFLOAT);
    //    std::cerr << " After eValFLOAT read\n";
    for (int i=0; i<eDimTot; i++) {
      float eValF=eValFLOAT[i];
      double eValD=double(eValF);
      eArr(i)=eValD;
    }
    delete [] eValFLOAT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[eDimTot];
    data.getVar(start, count, eValINT);
    for (int i=0; i<eDimTot; i++) {
      int eValI=eValINT[i];
      double eValD=double(eValI);
      eArr(i)=eValD;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_SHORT) {
    signed short int *eValINT;
    eValINT=new signed short int[eDimTot];
    data.getVar(start, count, eValINT);
    for (int i=0; i<eDimTot; i++) {
      double eValD=double(eValINT[i]);
      eArr(i)=eValD;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (!IsMatch) {
    std::cerr << "NC_ReadVariable_data_start_count : Did not find any matching number type\n";
    throw TerminalException{1};
  }
  // Now reading the offset and scaling_factor
  double eScal, eOff;
  try {
    //    std::cerr << "Before reading scale_factor\n";
    netCDF::NcVarAtt eScalAtt=data.getAtt("scale_factor");
    //    std::cerr << "After reading scale_factor\n";
    if (eScalAtt.isNull()) {
      eScal=1;
    }
    else {
      eScalAtt.getValues(&eScal);
    }
  }
  catch (...) {
    eScal=1;
  }
  try {
    netCDF::NcVarAtt eOffAtt=data.getAtt("add_offset");
    if (eOffAtt.isNull()) {
      eOff=0;
    }
    else {
      eOffAtt.getValues(&eOff);
    }
  }
  catch (...) {
    eOff=0;
  }
  for (int i=0; i<eDimTot; i++)
    eArr(i) = eOff + eScal*eArr(i);
  return eArr;
}


MyVector<double> NC_ReadVariable_data(netCDF::NcVar const& data)
{
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbDim=ListDim.size();
  std::vector<size_t> start(nbDim, 0);
  return NC_ReadVariable_data_start_count(data, start, ListDim);
}





std::vector<size_t> NC_ReadVariable_listdim_file(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_ReadVariable_listdim(data);
}


MyMatrix<int> NC_Read2Dvariable_Mask_file(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read2Dvariable_Mask_data(data);
}






MyMatrix<double> NC_Read2Dvariable_data(netCDF::NcVar const& data)
{
  MyVector<double> VecData = NC_ReadVariable_data(data);
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int eta=ListDim[0];
  int xi =ListDim[1];
  MyMatrix<double> MatData(eta, xi);
  int idx=0;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      MatData(i,j)=VecData(idx);
      idx++;
    }
  return MatData;
}








MyMatrix<double> NC_Read2Dvariable(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read2Dvariable_data(data);
}






MyMatrix<int> NC_Read2Dvariable_int_data(netCDF::NcVar const& data)
{
  if(data.isNull()) {
    std::cerr << "NC_Read2Dvariable_int_data : data is null\n";
    throw TerminalException{1};
  }
  netCDF::NcType eType=data.getType();
  if (eType.isNull()) {
    std::cerr << "NC_Read2Dvariable_int_data : eType is null\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2) {
    std::cerr << "NC_Read2Dvariable_int_data : The number of dimensions is not correct\n";
    throw TerminalException{1};
  }
  netCDF::NcDim eDim=data.getDim(0);
  int eta=eDim.getSize();
  netCDF::NcDim fDim=data.getDim(1);
  int xi=fDim.getSize();
  MyMatrix<int> eArr(eta, xi);
  bool IsMatch=false;
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[eta*xi];
    data.getVar(eValINT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	int eValI=eValINT[idx];
	double eValD=double(eValI);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_SHORT) {
    signed short int *eValINT;
    eValINT=new signed short int[eta*xi];
    data.getVar(eValINT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	signed short int eValI=eValINT[idx];
	double eValD=double(eValI);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_BYTE) {
    signed char *eValINT;
    eValINT=new signed char[eta*xi];
    data.getVar(eValINT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	signed short int eValI=eValINT[idx];
	double eValD=double(eValI);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValINT;
    IsMatch=true;
  }
  if (!IsMatch) {
    std::cerr << "NC_Read2Dvariable_int_data : Error in the call\n";
    throw TerminalException{1};
  }
  return eArr;
}


MyMatrix<int> NC_Read2Dvariable_int(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read2Dvariable_int_data(data);
}








MyVector<double> NC_Read1Dvariable_data(netCDF::NcVar const& data)
{
  return NC_ReadVariable_data(data);
}

MyVector<double> NC_Read1Dvariable(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read1Dvariable_data(data);
}











MyVector<int> NC_Read1Dvariable_int_data(netCDF::NcVar const& data)
{
  if (data.isNull()) {
    std::cerr << "NC_Read1Dvariable_int_data : Error in accessing to data\n";
    throw TerminalException{1};
  }
  netCDF::NcType eType=data.getType();
  if (eType.isNull()) {
    std::cerr << "NC_Read1Dvariable_int_data : Error in accessing to type information\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 1) {
    std::cerr << "NC_Read1Dvariable_int_data : The number of dimensions is not correct\n";
    throw TerminalException{1};
  }
  netCDF::NcDim eDim=data.getDim(0);
  int dim=eDim.getSize();
  MyVector<int> eArr(dim);
  bool IsMatch=false;
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[dim];
    data.getVar(eValINT);
    for (int i=0; i<dim; i++) {
      int eValI=eValINT[i];
      //    std::cerr << "i=" << i << " eValI=" << eValI << "\n";
      eArr(i)=eValI;
      //  std::cerr << "After the write\n";
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_SHORT) {
    signed short int *eValINT;
    eValINT=new signed short int[dim];
    data.getVar(eValINT);
    for (int i=0; i<dim; i++) {
      int eValI=int(eValINT[i]);
      eArr(i)=eValI;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_BYTE) {
    signed char *eVal;
    eVal=new signed char[dim];
    data.getVar(eVal);
    for (int i=0; i<dim; i++) {
      int eValI=int(eVal[i]);
      eArr(i)=eValI;
    }
    delete [] eVal;
    IsMatch=true;
  }
  if (!IsMatch) {
    std::cerr << "NC_Read1Dvariable_int_data : Did not find any matching number type\n";
    throw TerminalException{1};
  }
  return eArr;
}


MyVector<int> NC_Read1Dvariable_int(std::string const& eFile, std::string const& eVar)
{
  CheckNetcdfDataArray("NC_Read1Dvariable", eFile, eVar);
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  return NC_Read1Dvariable_int_data(data);
}











void CF_EXTRACT_TIME(std::string const& eStrUnitTime, double & ConvertToDay, double & eTimeStart)
{
  std::string YnameYear, YnameMonth, YnameDay;
  std::string YnameHour, YnameMin, YnameSec;
  std::string YnameD, YnameE;
  std::string YnameDate, YnameTime, YnameTimeP;
  std::string eStrTime;
  std::string strSpace=" ";
  int posBlank=STRING_GetCharPositionInString(eStrUnitTime, strSpace);
  std::string Xname=eStrUnitTime.substr(0, posBlank);
  int IsDone=0;
  if (Xname == "days") {
    IsDone=1;
    ConvertToDay=double(1);
  }
  if (Xname == "hours") {
    IsDone=1;
    ConvertToDay=double(1)/double(24);
  }
  if (Xname == "seconds") {
    IsDone=1;
    ConvertToDay=double(1)/double(86400);
  }
  if (IsDone == 0) {
    std::cerr << "We did not find a match for the time unit\n";
    std::cerr << "eStrUnitTime=" << eStrUnitTime << "\n";
    std::cerr << "Xname=" << Xname << "\n";
    std::cerr << "allowed Xname=days/hours/seconds\n";
    throw TerminalException{1};
  }
  //
  int alen=eStrUnitTime.length();
  std::string Yname=eStrUnitTime.substr(posBlank+1, alen - 1 - posBlank);
  int alenB=Yname.length();
  int posBlankB=STRING_GetCharPositionInString(Yname, strSpace);
  std::string YnameB=Yname.substr(posBlankB+1, alenB - 1 - posBlankB);
  // YnameB should be 1990-01-01 0:0:0 or something like that

  std::string strT="T";
  std::string strZ="Z";
  std::vector<std::string> LStrDateT=STRING_Split(YnameB, strT);
  int sizStrDateT=LStrDateT.size();
  if (sizStrDateT > 1) {
    // case of WW3 that has dates with the file format
    // "days since 1990-01-01T00:00:00Z"
    YnameDate=LStrDateT[0]; // should be 1990-01-01
    std::string eStrB=LStrDateT[1]; // 00:00:00Z
    int alenC=eStrUnitTime.length();
    YnameTime=eStrB.substr(0,alenC-2);
    //    std::cerr << "Case of WW3\n";
    //    std::cerr << "YnameDate=" << YnameDate << "\n";
    //    std::cerr << "YnameTime=" << YnameTime << "\n";
  }
  else {
    std::vector<std::string> LStrDate=STRING_Split(YnameB, strSpace);
    int sizStrDate=LStrDate.size();
    if (sizStrDate > 1) {
      YnameDate=LStrDate[0]; // should be 1990-01-01
      YnameTime=LStrDate[1]; // should be 0:0:0
      if (sizStrDate > 2) {
	std::string StrShift=LStrDate[2];
	if (StrShift != "GMT") {
	  std::cerr << "Yname=" << Yname << "\n";
	  std::cerr << "YnameB=" << YnameB << "\n";
	  std::cerr << "YnameDate=" << YnameDate << "\n";
	  std::cerr << "YnameTimeP=" << YnameTimeP << "\n";
	  std::cerr << "Need to program that case\n";
	  std::cerr << "Basically time can be of the form 0:0:0 GMT\n";
	  std::cerr << "or other stuff like that\n";
	  throw TerminalException{1};
	}
      }
    }
    else {
      YnameDate=LStrDate[0];
      YnameTime="00:00:00";
    }
  }
  //
  std::vector<std::string> eVectDate=STRING_Split(YnameDate,"-");
  std::string eStrYear, eStrMonth, eStrDay;
  eStrYear=eVectDate[0];
  eStrMonth=eVectDate[1];
  eStrDay=eVectDate[2];
  int year, month, day;
  std::istringstream(eStrYear) >> year;
  std::istringstream(eStrMonth) >> month;
  std::istringstream(eStrDay) >> day;
  // 
  std::vector<std::string> eVectTime=STRING_Split(YnameTime,":");
  std::string eStrHour, eStrMin, eStrSec;
  eStrHour=eVectTime[0];
  eStrMin=eVectTime[1];
  eStrSec=eVectTime[2];
  int hour, min, sec;
  std::istringstream(eStrHour) >> hour;
  std::istringstream(eStrMin) >> min;
  std::istringstream(eStrSec) >> sec;
  // Now collating
  eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, sec});
}




std::vector<double> NC_ReadTimeFromFile(std::string const& eFile, std::string const& StringTime)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ReadTimeFromFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "eFile=" << eFile << "\n";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(StringTime);
  if (data.isNull()) {
    std::cerr << "Error in accessing to the file (Case 5)\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "StringTime = " << StringTime << "\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 1) {
    std::cerr << "The number of dimensions is not correct\n";
    throw TerminalException{1};
  }
  netCDF::NcDim eDim=data.getDim(0);
  int siz=eDim.getSize();
  //  std::cerr << "siz=" << siz << "\n";
  double *eVal;
  eVal=new double[siz];
  data.getVar(eVal);
  netCDF::NcVarAtt eTimeAtt=data.getAtt("units");
  char eString[1024]="";
  eTimeAtt.getValues(eString);
  std::string eStrUnitTime=eString;
  double ConvertToDay, eTimeStart;
  CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart);
  std::vector<double> LTime(siz);
  if (siz == 0) {
    std::cerr << "We found siz=0\n";
    std::cerr << "This means that we find zero times in the file\n";
    std::cerr << "which is not what we expected\n";
    throw TerminalException{1};
  }
  double minTime=0, maxTime=0;
  for (int i=0; i<siz; i++) {
    //    std::cerr << "i=" << i << " eVal=" << eVal[i] << "\n";
    //    std::cerr << "ConvertToDay=" << ConvertToDay << "\n";
    double eTimeDay = eVal[i]*ConvertToDay + eTimeStart;
    if (i == 0) {
      minTime=eTimeDay;
      maxTime=eTimeDay;
    }
    else {
      if (eTimeDay > maxTime)
	maxTime=eTimeDay;
      if (eTimeDay < minTime)
	minTime=eTimeDay;
    }
    LTime[i]=eTimeDay;
  }
  delete [] eVal;
  bool ShowMinMax=false;
  if (ShowMinMax) {
    std::cerr << "minTime=" << minTime << " maxTime=" << maxTime << "\n";
    std::string strPresMin=DATE_ConvertMjd2mystringPres(minTime);
    std::string strPresMax=DATE_ConvertMjd2mystringPres(maxTime);
    std::cerr << "strPresMin=" << strPresMin << "\n";
    std::cerr << "strPresMax=" << strPresMax << "\n";
  }
  return LTime;
}

MyMatrix<double> NETCDF_Get2DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get2DvariableSpecEntry_FD\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "NC_Read2DvariableSpecEntry_FD : dataFile is null\n";
    throw TerminalException{1};
  }
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "NC_Read2DvariableSpecEntry_FD : data is null\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 2)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  int nbDim=ListDim.size();
  if (nbDim == 3) {
    std::vector<size_t> start{size_t(iRec), 0, 0};
    size_t eta=ListDim[1];
    size_t xi=ListDim[2];
    std::vector<size_t> count{1, eta, xi};
    MyVector<double> eVal=NC_ReadVariable_data_start_count(data, start, count);
    MyMatrix<double> eArr(eta, xi);
    int idx=0;
    for (size_t i=0; i<eta; i++)
      for (size_t j=0; j<xi; j++) {
	eArr(i, j)=eVal[idx];
	idx++;
      }
    return eArr;  
  }
  if (nbDim != 2) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "We should have nbDim = 2 or 3. nbDim=" << nbDim << "\n";
    throw TerminalException{1};
  }
  int nbWet=ListDim[1];
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, size_t(nbWet)};
  MyVector<double> eVal=NC_ReadVariable_data_start_count(data, start, count);
  if (nbWet == GrdArr.GrdArrRho.nbWet) {
    int eta=GrdArr.GrdArrRho.LON.rows();
    int xi=GrdArr.GrdArrRho.LON.cols();
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrRho.Idx[iWet];
      int j=GrdArr.GrdArrRho.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }
  if (nbWet == GrdArr.GrdArrU.nbWet) {
    int eta=GrdArr.GrdArrU.LON.rows();
    int xi=GrdArr.GrdArrU.LON.cols();
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrU.Idx[iWet];
      int j=GrdArr.GrdArrU.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }
  if (nbWet == GrdArr.GrdArrV.nbWet) {
    int eta=GrdArr.GrdArrV.LON.rows();
    int xi=GrdArr.GrdArrV.LON.cols();
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrV.Idx[iWet];
      int j=GrdArr.GrdArrV.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }  
  std::cerr << "Routine is NETCDF_Get2DvariableSpecEntry_FD\n";
  std::cerr << "eVar=" << eVar << "\n";
  std::cerr << "We did not find the size\n";
  throw TerminalException{1};
}



MyMatrix<double> NETCDF_Get2DvariableSpecEntry_FE(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get2DvariableSpecEntry_FE\n";
    std::cerr << "eVar=" << eVar << " iRec=" << iRec << "\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NETCDF_Get2DvariableSpecEntry_FE\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbDim=ListDim.size();
  if (nbDim != 2) {
    std::cerr << "nbDim=" << nbDim << " but it should be 2\n";
    throw TerminalException{1};
  }
  int nbRec=ListDim[0];
  int mnp=ListDim[1];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 3)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, size_t(mnp)};
  MyVector<double> eVal=NC_ReadVariable_data_start_count(data, start, count);
  MyMatrix<double> eArr(mnp, 1);
  for (int i=0; i<mnp; i++)
    eArr(i,0) = eVal(i);
  //  std::cerr << "eVar=" << eVar << " iRec=" << iRec << " eArr(min/max)=" << eArr.minCoeff() << " / " << eArr.maxCoeff() << "\n";
  if (GrdArr.L_IndexSelect) {
    int siz=GrdArr.I_IndexSelect.size();
    MyMatrix<double> eArrRet(siz, 1);
    for (int i=0; i<siz; i++) {
      int iGlob=GrdArr.I_IndexSelect[i];
      eArrRet(i,0)=eArr(iGlob,0);
    }
    return eArrRet;
  }
  return eArr;
}



Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_FE(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get2DvariableSpecEntry_FE\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NETCDF_Get3DvariableSpecEntry_FE\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbDim=ListDim.size();
  if (nbDim != 3) {
    std::cerr << "This command will certainly not work\n";
    std::cerr << "Dimensions are not correct\n";
    throw TerminalException{1};
  }
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 4)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  size_t mnp=ListDim[1];
  size_t NTR=ListDim[2];
  //
  std::vector<size_t> start{size_t(iRec), 0, 0};
  std::vector<size_t> count{1, mnp, NTR};
  MyVector<double> eVal = NC_ReadVariable_data_start_count(data, start, count);
  Eigen::Tensor<double,3> eArr(int(NTR), int(mnp), 1);
  int idx=0;
  for (size_t i=0; i<mnp; i++) {
    for (size_t iTr=0; iTr<NTR; iTr++) {
      eArr(iTr,i,0)=eVal(idx);
      idx++;
    }
  }
  if (GrdArr.L_IndexSelect) {
    int siz=GrdArr.I_IndexSelect.size();
    Eigen::Tensor<double,3> eArrRet(siz, int(NTR), 1);
    for (int i=0; i<siz; i++) {
      int iGlob=GrdArr.I_IndexSelect[i];
      for (size_t iTr=0; iTr<NTR; iTr++)
	eArrRet(iTr,i,0)=eArr(iTr,iGlob,0);
    }
    return eArrRet;
  }
  return eArr;
}

void NETCDF_Write2DvariableSpecEntry(std::string const& eFile, std::string const& eVar, int const& iRec, MyMatrix<double> const& M)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Write2DvariableSpecEntry\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::write);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in accessing to the file (Case 9)\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    std::cerr << "iRec  = " << iRec << "\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2 && nbDim != 3) {
    std::cerr << "Inconsistency in the number of dimension. Allowed is 2 or 3\n";
    std::cerr << "nbDim=" << nbDim << "\n";
    throw TerminalException{1};
  }
  std::cerr << "nbDim=" << nbDim << "\n";
  netCDF::NcDim eDim = data.getDim(0);
  int nbRec=eDim.getSize();
  if (iRec > nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 5)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec <= nbRec\n";
    std::cerr << "Only increases by 1 at most\n";
    throw TerminalException{1};
  }
  size_t eProd=1;
  std::vector<size_t> ListDim;
  for (int iDim=1; iDim<nbDim; iDim++) {
    eDim=data.getDim(1);
    size_t eSize=eDim.getSize();
    ListDim.push_back(eSize);
    eProd *= eSize;
  }
  std::vector<size_t> start;
  std::vector<size_t> count;
  if (nbDim == 2) {
    start={size_t(iRec), 0};
    count={1, ListDim[0]};
  }
  if (nbDim == 3) {
    start={size_t(iRec), 0, 0};
    count={1, ListDim[0], ListDim[1]};
  }
  netCDF::NcType eType=data.getType();
  bool IsDone=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[eProd];
    for (int i=0; i<M.size(); i++)
      eVal[i]=M(i);
    data.putVar(start, count, eVal);
    delete [] eVal;
    IsDone=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eVal;
    eVal=new float[eProd];
    for (int i=0; i<M.size(); i++)
      eVal[i]=float(M(i));
    data.putVar(start, count, eVal);
    delete [] eVal;
    IsDone=true;
  }
  if (!IsDone) {
    std::cerr << "Data wriding failed for the function\n";
    throw TerminalException{1};
  }
}




void NETCDF_Write3DvariableSpecEntry(std::string const& eFile, std::string const& eVar, int const& iRec, Eigen::Tensor<double,3> const& Tens)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Write2DvariableSpecEntry\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::write);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in accessing to the file (Case 9)\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    std::cerr << "iRec  = " << iRec << "\n";
    throw TerminalException{1};
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2 && nbDim != 3) {
    std::cerr << "Inconsistency in the number of dimension. Allowed is 2 or 3\n";
    std::cerr << "nbDim=" << nbDim << "\n";
    throw TerminalException{1};
  }
  std::cerr << "nbDim=" << nbDim << "\n";
  netCDF::NcDim eDim = data.getDim(0);
  int nbRec=eDim.getSize();
  if (iRec > nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 5)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec <= nbRec\n";
    std::cerr << "Only increases by 1 at most\n";
    throw TerminalException{1};
  }
  size_t eProd=1;
  std::vector<size_t> ListDim;
  for (int iDim=1; iDim<nbDim; iDim++) {
    eDim=data.getDim(1);
    size_t eSize=eDim.getSize();
    ListDim.push_back(eSize);
    eProd *= eSize;
  }
  std::vector<size_t> start{size_t(iRec)};
  std::vector<size_t> count{1};
  for (int iDim=1; iDim<nbDim; iDim++) {
    start.push_back(0);
    count.push_back(ListDim[iDim-1]);
  }
  netCDF::NcType eType=data.getType();
  bool IsDone=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[eProd];
    for (int i=0; i<Tens.size(); i++)
      eVal[i]=Tens(i);
    data.putVar(start, count, eVal);
    delete [] eVal;
    IsDone=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eVal;
    eVal=new float[eProd];
    for (int i=0; i<Tens.size(); i++)
      eVal[i]=float(Tens(i));
    data.putVar(start, count, eVal);
    delete [] eVal;
    IsDone=true;
  }
  if (!IsDone) {
    std::cerr << "Data wriding failed for the function\n";
    throw TerminalException{1};
  }
}







MyMatrix<double> NETCDF_Get2DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  //  std::cerr << "NETCDF_Get2DvariableSpecEntry, eFile=" << eFile << "\n";
  if (GrdArr.IsFE == 1)
    return NETCDF_Get2DvariableSpecEntry_FE(eFile, GrdArr, eVar, iRec);
  return NETCDF_Get2DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}







Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_ROMS_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get3DvariableSpecEntry_ROMS_FD\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  int s_rho=NC_ReadDimension(dataFile, "s_rho");
  int s_w=NC_ReadDimension(dataFile, "s_w");
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NETCDF_Get3DvariableSpecEntry_ROMS_FD\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 6)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  int nbDim=ListDim.size();
  if (nbDim == 4) {
    int s_vert=ListDim[1];
    int eta=ListDim[2];
    int xi=ListDim[3];
    std::vector<size_t> start{size_t(iRec), 0, 0, 0};
    std::vector<size_t> count{1, size_t(s_vert), size_t(eta), size_t(xi)};
    MyVector<double> eVal = NC_ReadVariable_data_start_count(data, start, count);
    Eigen::Tensor<double,3> eArr(s_vert, eta, xi);
    int idx=0;
    for (int k=0; k<s_vert; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++) {
	  eArr(k, i, j)=eVal[idx];
	  idx++;
	}
    return eArr;  
  }
  // The file format uses only WET points.
  int nbWet=ListDim[1];
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, size_t(nbWet)};
  MyVector<double> eVal = NC_ReadVariable_data_start_count(data, start, count);
  if (nbWet == s_rho*GrdArr.GrdArrRho.nbWet) {
    int eta=GrdArr.GrdArrRho.LON.rows();
    int xi=GrdArr.GrdArrRho.LON.cols();
    Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr(k, i, j)=0;
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrRho.nbWet; iWet++) {
	int i=GrdArr.GrdArrRho.Idx[iWet];
	int j=GrdArr.GrdArrRho.Jdx[iWet];
	eArr(k, i, j)=eVal(idx);
	idx++;
      }
    }
    return eArr;
  }
  if (nbWet == s_rho*GrdArr.GrdArrU.nbWet) {
    int eta=GrdArr.GrdArrU.LON.rows();
    int xi=GrdArr.GrdArrU.LON.cols();
    Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr(k, i, j)=0;
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrU.nbWet; iWet++) {
	int i=GrdArr.GrdArrU.Idx[iWet];
	int j=GrdArr.GrdArrU.Jdx[iWet];
	eArr(k, i, j)=eVal(idx);
	idx++;
      }
    }
    return eArr;
  }
  if (nbWet == s_rho*GrdArr.GrdArrV.nbWet) {
    int eta=GrdArr.GrdArrV.LON.rows();
    int xi=GrdArr.GrdArrV.LON.cols();
    Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr(k, i, j)=0;
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrV.nbWet; iWet++) {
	int i=GrdArr.GrdArrV.Idx[iWet];
	int j=GrdArr.GrdArrV.Jdx[iWet];
	eArr(k, i, j)=eVal(idx);
	idx++;
      }
    }
    return eArr;
  }
  std::cerr << "Routine is NETCDF_Get3DvariableSpecEntry_ROMS_FD\n";
  std::cerr << "eVar = " << eVar << "\n";
  std::cerr << "s_rho = " << s_rho << " s_w=" << s_w << "\n";
  std::cerr << "nbWet = " << nbWet << "\n";
  std::cerr << "nbWetRho = " << s_rho * GrdArr.GrdArrRho.nbWet << "\n";
  std::cerr << "  nbWetU = " << s_rho * GrdArr.GrdArrU.nbWet << "\n";
  std::cerr << "  nbWetV = " << s_rho * GrdArr.GrdArrV.nbWet << "\n";
  std::cerr << "We did not find the size\n";
  throw TerminalException{1};
}





Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_Direct_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec, std::string const& dimVert)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get3DvariableSpecEntry_FD\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  int s_vert_read=NC_ReadDimension(dataFile, dimVert);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NETCDF_Get3DvariableSpecEntry_FD\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 6)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  int nbDim=ListDim.size();
  if (nbDim != 4) {
    std::cerr << "nbDim = " << nbDim << " but it should be equal to 4\n";
    throw TerminalException{1};
  }
  int s_vert=ListDim[1];
  if (s_vert != s_vert_read) {
    std::cerr << "s_vert=" << s_vert << " s_vert_read=" << s_vert_read << " but should be equal\n";
    throw TerminalException{1};
  }
  int eta=ListDim[2];
  int xi=ListDim[3];
  std::vector<size_t> start{size_t(iRec), 0, 0, 0};
  std::vector<size_t> count{1, size_t(s_vert), size_t(eta), size_t(xi)};
  MyVector<double> eVal = NC_ReadVariable_data_start_count(data, start, count);
  Eigen::Tensor<double,3> eArr(s_vert, eta, xi);
  int idx=0;
  for (int k=0; k<s_vert; k++)
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	eArr(k, i, j)=eVal[idx];
	idx++;
      }
  return eArr;  
}


Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (GrdArr.ModelName == "ROMS")
    return NETCDF_Get3DvariableSpecEntry_ROMS_FD(eFile, GrdArr, eVar, iRec);
  if (GrdArr.ModelName == "HYCOM")
    return NETCDF_Get3DvariableSpecEntry_Direct_FD(eFile, GrdArr, eVar, iRec, "depth");
  if (GrdArr.ModelName == "NEMO")
    return NETCDF_Get3DvariableSpecEntry_Direct_FD(eFile, GrdArr, eVar, iRec, "depth");
  std::cerr << "The InfoVertical is not matched. GrdArr.ModelName=" << GrdArr.ModelName << "\n";
  std::cerr << "Maybe missing code\n";
  throw TerminalException{1};
}


Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (GrdArr.IsFE == 1)
    return NETCDF_Get3DvariableSpecEntry_FE(eFile, GrdArr, eVar, iRec);
  return NETCDF_Get3DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}


struct TimeNEMO {
  std::vector<int> ListIdxTime1;
  std::vector<int> ListIdxTime2;
  double mjdDay1, mjdDay2;
  double alpha1, alpha2;
};

TimeNEMO GetPairListTime(double const& mjdDay)
{
  std::vector<int> vectTime=DATE_ConvertMjd2six(mjdDay);
  int eYear1=vectTime[0];
  int eMonth1=vectTime[1];
  int eDay1=vectTime[2];
  int eHour1=vectTime[3];
  double mjdDay1=DATE_ConvertSix2mjd({eYear1, eMonth1, eDay1, eHour1, 0, 0});
  double mjdDay2=mjdDay1 + double(1)/double(24);
  std::vector<int> vectTime2=DATE_ConvertMjd2six(mjdDay2);
  int eYear2=vectTime2[0];
  int eMonth2=vectTime2[1];
  int eDay2=vectTime2[2];
  int eHour2=vectTime2[3];
  double alpha1=(mjdDay2 - mjdDay ) / (mjdDay2 - mjdDay1);
  double alpha2=(mjdDay  - mjdDay1) / (mjdDay2 - mjdDay1);
  return {{eYear1, eMonth1, eDay1, eHour1}, {eYear2, eMonth2, eDay2, eHour2}, mjdDay1, mjdDay2, alpha1, alpha2};
}



Eigen::Tensor<double,3> NEMO_Get3DvariableSpecEntry_Kernel(std::string const& eFile, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NEMO_Get3DvariableSpecEntry_Kernele\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NEMO_Get3DvariableSpecEntry_Kernel\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbDim=ListDim.size();
  if (nbDim != 4) {
    std::cerr << "nbDim=" << nbDim << "\n";
    std::cerr << "We should have dimension equal to 4\n";
    throw TerminalException{1};
  }
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 6)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  int s_vert=ListDim[1];
  int eta=ListDim[2];
  int xi=ListDim[3];
  std::vector<size_t> start{size_t(iRec), 0, 0, 0};
  std::vector<size_t> count{1, size_t(s_vert), size_t(eta), size_t(xi)};
  MyVector<double> eVal = NC_ReadVariable_data_start_count(data, start, count);
  Eigen::Tensor<double,3> eArr(s_vert, eta, xi);
  int idx=0;
  for (int k=0; k<s_vert; k++)
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	eArr(k, i, j)=eVal(idx);
	idx++;
      }
  return eArr;
}


Eigen::Tensor<double,3> NEMO_Get3DvariableSpecEntry(std::string const& HisPrefix, std::vector<int> const& ListIdxTime, std::string const& eVarName1, std::string const& eVarName2)
{
  int eYear=ListIdxTime[0];
  int eMonth=ListIdxTime[1];
  int eDay=ListIdxTime[2];
  int eHour=ListIdxTime[3];
  std::string strYear = StringNumber(eYear,4);
  std::string strMon  = StringNumber(eMonth,2);
  std::string strDay  = StringNumber(eDay,2);
  std::string strHour = StringNumber(eHour,2);
  std::string eFile = HisPrefix + strYear + "/" + strMon + "/NEMO_" + strYear + "_" + strMon + "_" + strDay + "_" + eVarName1 + ".nc";
  //  std::cerr << "eFile = " << eFile << "\n";
  //  std::cerr << "eVarName1=" << eVarName1 << " eVarName2=" << eVarName2 << "\n";
  int iRec = eHour;
  //  std::cerr << "iRec=" << iRec << "\n";
  return NEMO_Get3DvariableSpecEntry_Kernel(eFile, eVarName2, iRec);
}


Eigen::Tensor<double,3> NEMO_Get3DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVarName1, std::string const& eVarName2, double const& mjdDay)
{
  std::string HisPrefix=TotalArr.eArr.HisPrefix;
  TimeNEMO recNEMO = GetPairListTime(mjdDay);
  //  std::cerr << "Before getting Tens1\n";
  Eigen::Tensor<double,3> Tens1 = NEMO_Get3DvariableSpecEntry(HisPrefix, recNEMO.ListIdxTime1,  eVarName1, eVarName2);
  //  std::cerr << "Before getting Tens2\n";
  Eigen::Tensor<double,3> Tens2 = NEMO_Get3DvariableSpecEntry(HisPrefix, recNEMO.ListIdxTime2,  eVarName1, eVarName2);
  //  std::cerr << " After getting Tens2\n";
  auto LDim=Tens1.dimensions();
  int s_vert=LDim[0];
  int eta=LDim[1];
  int xi=LDim[2];
  //  does not COMPILE:
  //  Eigen::Tensor<double,3> RetVar=alphaLow*eVarLow + alphaUpp*eVarUpp;
  Eigen::Tensor<double,3> RetVar(s_vert, eta, xi);
  for (int k=0; k<s_vert; k++)
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	RetVar(s_vert - 1 - k, i, j)=recNEMO.alpha1*Tens1(k, i, j) + recNEMO.alpha2*Tens2(k, i, j);
  return RetVar;
}



MyMatrix<double> NEMO_Get2DvariableSpecEntry_Kernel(std::string const& eFile, std::string const& eVar, int const& iRec)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "NETCDF_Get2DvariableSpecEntry_FD\n";
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if (data.isNull()) {
    std::cerr << "Error in NEMO_Get2DvariableSpecEntry_Kernel\n";
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "eVar  = " << eVar << "\n";
    throw TerminalException{1};
  }
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbDim=ListDim.size();
  if (nbDim != 3) {
    std::cerr << "nbDim=" << nbDim << " but should be 3\n";
    throw TerminalException{1};
  }
  int nbRec=ListDim[0];
  if (iRec < 0 || iRec >= nbRec) {
    std::cerr << "eFile = " << eFile << "\n";
    std::cerr << "Error, iRec is too large (Case 2)\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    std::cerr << "We need C-convention iRec < nbRec\n";
    throw TerminalException{1};
  }
  std::vector<size_t> start{size_t(iRec), 0, 0};
  size_t eta=ListDim[1];
  size_t xi=ListDim[2];
  std::vector<size_t> count{1, eta, xi};
  MyVector<double> eVal=NC_ReadVariable_data_start_count(data, start, count);
  MyMatrix<double> eArr(eta, xi);
  int idx=0;
  for (size_t i=0; i<eta; i++)
    for (size_t j=0; j<xi; j++) {
      eArr(i, j)=eVal(idx);
      idx++;
    }
  return eArr;  
}




MyMatrix<double> NEMO_Get2DvariableSpecEntry(std::string const& HisPrefix, std::vector<int> const& ListIdxTime, std::string const& eVarName1, std::string const& eVarName2)
{  
  int eYear=ListIdxTime[0];
  int eMonth=ListIdxTime[1];
  int eDay=ListIdxTime[2];
  int eHour=ListIdxTime[3];
  std::string strYear = StringNumber(eYear,4);
  std::string strMon  = StringNumber(eMonth,2);
  std::string strDay  = StringNumber(eDay,2);
  std::string strHour = StringNumber(eHour,2);
  std::string eFile = HisPrefix + strYear + "/" + strMon + "/NEMO_" + strYear + "_" + strMon + "_" + strDay + "_" + eVarName1 + ".nc";
  int iRec = eHour;
  return NEMO_Get2DvariableSpecEntry_Kernel(eFile, eVarName2, iRec);
}


MyMatrix<double> NEMO_Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVarName1, std::string const& eVarName2, double const& mjdDay)
{
  std::string HisPrefix=TotalArr.eArr.HisPrefix;
  TimeNEMO recNEMO = GetPairListTime(mjdDay);
  MyMatrix<double> Mat1 = NEMO_Get2DvariableSpecEntry(HisPrefix, recNEMO.ListIdxTime1, eVarName1, eVarName2);
  MyMatrix<double> Mat2 = NEMO_Get2DvariableSpecEntry(HisPrefix, recNEMO.ListIdxTime2, eVarName1, eVarName2);
  return recNEMO.alpha1*Mat1 + recNEMO.alpha2*Mat2;
}






#endif
