#ifndef INCLUDE_RIVER_FUNCTION
#define INCLUDE_RIVER_FUNCTION

#include "Basic_file.h"
#include "Basic_plot.h"
#include "Model_grids.h"
#include "ROMSfunctionality.h"
#include "SphericalGeom.h"
#include "Namelist.h"
#include "Basic_netcdf.h"
#include "mjdv2.h"
#include "SVGfunctions.h"


FullNamelist Individual_River_File()
{
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["name"]="";
  ListStringValues1["TypeVaryingTransport"]="YearlyClimatology"; // other possibility is Po
  ListStringValues1["TypeVaryingTemperature"]="unset"; // possibilities: Constant
  ListStringValues1["TypeVaryingSalinity"]="unset"; // possibilities: Constant
  ListStringValues1["verticalShapeOption"]="UpperLayer"; // possibilities: Constant
  ListListDoubleValues1["ListMonthlyFlux"] = {};
  ListListDoubleValues1["ListMonthlyTemp"] = {};
  ListDoubleValues1["ConstantFlux"] = -1;
  ListDoubleValues1["ConstantFactorFlux"] = 1.0;
  ListStringValues1["PoPrefixData"] = "unset";
  ListStringValues1["WScase"] = "River";
  ListStringValues1["PrefixPoData"] = "unset";
  ListBoolValues1["SetRiverTemperature"]=true;
  ListBoolValues1["SetRiverSalinity"]=true;
  ListDoubleValues1["ConstantRiverTemperature"]=14;
  ListDoubleValues1["ConstantRiverSalinity"]=0;
  ListDoubleValues1["lon"] = -400;
  ListDoubleValues1["lat"] = -400;
  ListDoubleValues1["direction"] = -400;
  ListDoubleValues1["MaxDepth"] = 2;
  ListDoubleValues1["targetDepth"] = 2;
  ListDoubleValues1["FrequencyDay"] = 7;
  ListDoubleValues1["DurationHour"] = 2;
  ListDoubleValues1["TotalFlux"] = 15300;
  ListIntValues1["iSelect"]=-1;
  ListIntValues1["jSelect"]=-1;
  ListIntValues1["SignSelect"]=-400;
  ListIntValues1["DirSelect"]=-1;
  ListIntValues1["ChoiceSelect"]=-1;
  SingleBlock BlockDESC;
  BlockDESC.ListIntValues=ListIntValues1;
  BlockDESC.ListBoolValues=ListBoolValues1;
  BlockDESC.ListDoubleValues=ListDoubleValues1;
  BlockDESC.ListListDoubleValues=ListListDoubleValues1;
  BlockDESC.ListListIntValues=ListListIntValues1;
  BlockDESC.ListStringValues=ListStringValues1;
  BlockDESC.ListListStringValues=ListListStringValues1;
  ListBlock["DESCRIPTION"]=BlockDESC;
  //
  return {ListBlock, "undefined"};
}

struct PairTimeMeas {
  double time;
  double meas;
};

double InterpolateMeasurement(std::vector<PairTimeMeas> const& ListPairTimeMeas, double const& eTime)
{
  int siz=ListPairTimeMeas.size();
  if (siz == 0) {
    std::cerr << "We have |ListPairTimeMeas| = 0\n";
    std::cerr << "InterpolateMeasurement cannot be run correctly\n";
    throw TerminalException{1};
  }
  for (int i=1; i<siz; i++) {
    double time0=ListPairTimeMeas[i-1].time;
    double time1=ListPairTimeMeas[i].time;
    if (time0 <= eTime && eTime < time1) {
      double alpha0 = (time1 - eTime) / (time1 - time0);
      double alpha1 = (eTime - time0) / (time1 - time0);
      double meas0=ListPairTimeMeas[i-1].meas;
      double meas1=ListPairTimeMeas[i].meas;
      double eValInterp = alpha0 * meas0 + alpha1 * meas1;
      return eValInterp;
    }
  }
  std::cerr << "siz=" << siz << "\n";
  double minTime=ListPairTimeMeas[0].time;
  double maxTime=ListPairTimeMeas[siz-1].time;
  std::cerr << "Failed to find the right entry in the list of values\n";
  std::cerr << "eTime=" << eTime << " strPres=" << DATE_ConvertMjd2mystringPres(eTime) << "\n";
  std::cerr << "minTime=" << minTime << " strPres=" << DATE_ConvertMjd2mystringPres(minTime) << "\n";
  std::cerr << "maxTime=" << maxTime << " strPres=" << DATE_ConvertMjd2mystringPres(maxTime) << "\n";
  throw TerminalException{1};
}



std::vector<PairTimeMeas> ReadListPairTimeMeas_PoStyle(std::string const& PrefixData, std::string const& CharSel)
{
  std::cerr << "ReadListPairTimeMeas_PoStyle, step 1 CharSel=" << CharSel << "\n";
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(PrefixData);
  std::cerr << "ReadListPairTimeMeas_PoStyle, step 2 |ListFile|=" << ListFile.size() << "\n";
  std::vector<PairTimeMeas> ListPairTimeMeas;
  auto IsCorrectLine=[](std::string const& eLine) -> bool {
    std::vector<std::string> LStrA = STRING_Split(eLine, ",,,");
    if (LStrA.size() != 1)
      return false;
    return true;
  };
  auto GetPairTimeMeas=[](std::string const& eLine) -> PairTimeMeas {
    std::vector<std::string> LStrA = STRING_Split(eLine, ",,");
    if (LStrA.size() != 2) {
      std::cerr << "|LStrA|=" << LStrA.size() << " but should be 2\n";
      throw TerminalException{1};
    }
    std::vector<std::string> LStrB = STRING_Split(LStrA[0], ",");
    if (LStrB.size() != 4) {
      std::cerr << "eLine=" << eLine << "\n";
      std::cerr << "LStrA[0]=" << LStrA[0] << "\n";
      std::cerr << "|LStrB|=" << LStrB.size() << " but should be 4\n";
      throw TerminalException{1};
    }
    std::string strTime=LStrB[0];
    std::string strMeas=LStrB[3];
    double eMeas;
    std::istringstream(strMeas) >> eMeas;
    std::vector<std::string> LStrC = STRING_Split(strTime, ":");
    if (LStrC.size() != 3) {
      std::cerr << "|LStrC|=" << LStrC.size() << " but should be 3\n";
      throw TerminalException{1};
    }
    std::string strHour=LStrC[1];
    std::string strMin =LStrC[2];
    std::string strDate=LStrC[0];
    std::vector<std::string> LStrD = STRING_Split(strDate, "/");
    if (LStrD.size() != 3) {
      std::cerr << "|LStrD|=" << LStrD.size() << " but should be 3\n";
      throw TerminalException{1};
    }
    std::string strYear=LStrD[0];
    std::string strMon =LStrD[1];
    std::string strDay =LStrD[2];
    int eYear, eMon, eDay, eHour, eMin, eSec=0;
    std::istringstream(strYear) >> eYear;
    std::istringstream(strMon) >> eMon;
    std::istringstream(strDay) >> eDay;
    std::istringstream(strHour) >> eHour;
    std::istringstream(strMin) >> eMin;
    double eDate=DATE_ConvertSix2mjd({eYear,eMon,eDay,eHour,eMin,eSec});
    return {eDate, eMeas};
  };
  std::cerr << "ReadListPairTimeMeas_PoStyle, step 3\n";
  for (auto & eFile : ListFile) {
    std::string FirstChar = eFile.substr(0,1);
    if (FirstChar == CharSel) {
      std::string FullFile = PrefixData + eFile;
      std::vector<std::string> ListLine = ReadFullFile(FullFile);
      for (auto & eLine : ListLine)
	if (IsCorrectLine(eLine))
	  ListPairTimeMeas.push_back(GetPairTimeMeas(eLine));
    }
  }
  std::cerr << "ReadListPairTimeMeas_PoStyle, step 4 |ListPairTimeMeas|=" << ListPairTimeMeas.size() << "\n";
  std::sort(ListPairTimeMeas.begin(), ListPairTimeMeas.end(),
	    [](PairTimeMeas const& x1, PairTimeMeas const& x2) -> bool {
	      return x1.time < x2.time;
	    });
  std::cerr << "ReadListPairTimeMeas_PoStyle, step 5 |ListPairTimeMeas|=" << ListPairTimeMeas.size() << "\n";
  return ListPairTimeMeas;
}


struct DescriptionRiver {
  double lon;
  double lat;
  double direction;
  std::string verticalShapeOption;
  double MaxDepth;
  double targetDepth;
  double ConstantRiverTemperature;
  double ConstantRiverSalinity;
  bool SetRiverTemperature;
  bool SetRiverSalinity;
  std::vector<double> ListMonthlyFlux;
  std::vector<double> ListMonthlyTemp;
  double ConstantFlux;
  double ConstantFactorFlux;
  std::vector<PairTimeMeas> ListPairTimeFlux;
  std::vector<PairTimeMeas> ListPairTimeTemp;
  std::string PoPrefixData;
  std::string TypeVaryingTransport;
  std::string TypeVaryingTemperature;
  std::string TypeVaryingSalinity;
  std::string name;
  double FrequencyDay;
  double DurationHour;
  double TotalFlux;
  std::string WScase;
  int iSelect;
  int jSelect;
  int SignSelect;
  int DirSelect;
  int ChoiceSelect;
};



DescriptionRiver ReadRiverDescription(std::string const& RiverDescriptionFile)
{
  //  std::cerr << "ReadRiverDescription, step 1\n";
  FullNamelist eFull = Individual_River_File();
  //  std::cerr << "ReadRiverDescription, step 2\n";
  NAMELIST_ReadNamelistFile(RiverDescriptionFile, eFull);
  //  std::cerr << "ReadRiverDescription, step 3\n";
  SingleBlock eBlDESC=eFull.ListBlock.at("DESCRIPTION");
  //  std::cerr << "ReadRiverDescription, step 4\n";
  DescriptionRiver eDesc;
  //  std::cerr << "ReadRiverDescription, step 5\n";
  eDesc.lon = eBlDESC.ListDoubleValues.at("lon");
  eDesc.lat = eBlDESC.ListDoubleValues.at("lat");
  eDesc.direction = eBlDESC.ListDoubleValues.at("direction");
  eDesc.MaxDepth = eBlDESC.ListDoubleValues.at("MaxDepth");
  eDesc.targetDepth = eBlDESC.ListDoubleValues.at("targetDepth");
  eDesc.SetRiverTemperature = eBlDESC.ListBoolValues.at("SetRiverTemperature");
  eDesc.SetRiverSalinity    = eBlDESC.ListBoolValues.at("SetRiverSalinity");
  eDesc.ConstantRiverTemperature = eBlDESC.ListDoubleValues.at("ConstantRiverTemperature");
  eDesc.ConstantRiverSalinity    = eBlDESC.ListDoubleValues.at("ConstantRiverSalinity");
  eDesc.TypeVaryingTransport   = eBlDESC.ListStringValues.at("TypeVaryingTransport");
  eDesc.TypeVaryingTemperature = eBlDESC.ListStringValues.at("TypeVaryingTemperature");
  eDesc.TypeVaryingSalinity    = eBlDESC.ListStringValues.at("TypeVaryingSalinity");
  eDesc.verticalShapeOption    = eBlDESC.ListStringValues.at("verticalShapeOption");
  eDesc.FrequencyDay = eBlDESC.ListDoubleValues.at("FrequencyDay");
  eDesc.DurationHour = eBlDESC.ListDoubleValues.at("DurationHour");
  eDesc.TotalFlux = eBlDESC.ListDoubleValues.at("TotalFlux");
  eDesc.WScase    = eBlDESC.ListStringValues.at("WScase");
  eDesc.PoPrefixData = eBlDESC.ListStringValues.at("PoPrefixData");
  eDesc.name = eBlDESC.ListStringValues.at("name");
  eDesc.ListMonthlyFlux = eBlDESC.ListListDoubleValues.at("ListMonthlyFlux");
  eDesc.ListMonthlyTemp = eBlDESC.ListListDoubleValues.at("ListMonthlyTemp");
  eDesc.ConstantFlux = eBlDESC.ListDoubleValues.at("ConstantFlux");
  eDesc.ConstantFactorFlux = eBlDESC.ListDoubleValues.at("ConstantFactorFlux");
  //  std::cerr << "ReadRiverDescription, step 6\n";
  auto CheckListMonth=[&](std::vector<double> const& ListMon) -> void {
    int nbMonth=ListMon.size();
    std::cerr << "nbMonth=" << nbMonth << "\n";
    if (nbMonth != 12 && nbMonth != 0) {
      std::cerr << "RiverDescriptionFile = " << RiverDescriptionFile << "\n";
      std::cerr << "We have nbMonth=" << nbMonth << "\n";
      std::cerr << "It should be 12 or not set\n";
      throw TerminalException{1};
    }
  };
  CheckListMonth(eDesc.ListMonthlyFlux);
  CheckListMonth(eDesc.ListMonthlyTemp);
  if (eDesc.TypeVaryingTransport == "PoFlux") {
    std::string PrefixPoData  = eBlDESC.ListStringValues.at("PrefixPoData");
    eDesc.ListPairTimeFlux = ReadListPairTimeMeas_PoStyle(PrefixPoData, "Q");
  }
  if (eDesc.TypeVaryingTemperature == "PoTemp") {
    std::string PrefixPoData  = eBlDESC.ListStringValues.at("PrefixPoData");
    eDesc.ListPairTimeTemp = ReadListPairTimeMeas_PoStyle(PrefixPoData, "T");
  }
  eDesc.iSelect = eBlDESC.ListIntValues.at("iSelect");
  eDesc.jSelect = eBlDESC.ListIntValues.at("iSelect");
  eDesc.SignSelect = eBlDESC.ListIntValues.at("SignSelect");
  eDesc.DirSelect = eBlDESC.ListIntValues.at("DirSelect");
  eDesc.ChoiceSelect = eBlDESC.ListIntValues.at("ChoiceSelect");
  //  std::cerr << "ReadRiverDescription, step 7\n";
  return eDesc;
}



FullNamelist NAMELIST_PLOT_River()
{
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["GridFile"]="unset";
  ListStringValues1["RiverFile"]="unset";
  ListStringValues1["SVGfile"]="unset";
  ListStringValues1["PicPrefix"]="unset";
  ListStringValues1["Extension"]="png";
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110915.000000";
  ListStringValues1["__NaturePlot"]="RIVER";
  ListBoolValues1["FirstCleanDirectory"]=true;
  ListBoolValues1["KeepNC_NCL"]=false;
  ListBoolValues1["InPlaceRun"]=false;
  ListBoolValues1["PrintDebugInfo"]=false;
  ListBoolValues1["OnlyCreateFiles"]=false;
  ListIntValues1["NPROC"]=1;
  ListStringValues1["Pcolor_method"]="ncl";
  ListStringValues1["Quiver_method"]="ncl";
  ListStringValues1["Lines_method"]="ncl";
  ListStringValues1["Scatter_method"]="ncl";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  ListBlock["PROC"]=BlockPROC;
  //
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<int>> ListListIntValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["PlotFlux"]=false;
  ListStringValues2["BEGTC"]="20090915.000000";
  ListStringValues2["ENDTC"]="20110915.000000";
  ListDoubleValues2["WidthLine"]=0.4;
  ListDoubleValues2["WidthLineRiver"]=0.4;
  ListIntValues2["SizeText"]=2;
  ListIntValues2["nbLabel"]=5;
  ListStringValues2["StyleDate"]="unset";
  ListStringValues2["StyleTitle"]="unset";
  ListStringValues2["TitleString"]="unset";
  ListStringValues2["StyleTitle"]="unset";
  ListDoubleValues2["ShiftLonText"]=0.05;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  BlockPLOT.ListListIntValues=ListListIntValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  ListBlock["PLOT"]=BlockPLOT;
  //
  return {ListBlock, "undefined"};
}


struct ijSeaLand {
  int iSea;
  int jSea;
  int iLand;
  int jLand;
};


struct ijdsInfo {
  int eStatus;
  int iSelect;
  int jSelect;
  int DirSelect;
  int SignSelect;
};


ijSeaLand GetArrayIjSeaLand(ijdsInfo const& RecInfo)
{
  int iSelect=RecInfo.iSelect;
  int jSelect=RecInfo.jSelect;
  int DirSelect=RecInfo.DirSelect;
  int SignSelect=RecInfo.SignSelect;
  int iSea=-1, jSea=-1, iLand=-1, jLand=-1;
  if (DirSelect == 0 && SignSelect == 1) {
    // Right case in map_rivers.m
    iSea=iSelect+1;
    jSea=jSelect+1;
    iLand=iSea;
    jLand=jSea-1;
  }
  if (DirSelect == 1 && SignSelect == 1) {
    // Up case in map_rivers.m
    iSea=iSelect+1;
    jSea=jSelect+1;
    iLand=iSea-1;
    jLand=jSea;
  }
  if (DirSelect == 0 && SignSelect == -1) {
    // Left case in map_rivers.m
    iSea=iSelect+1;
    jSea=jSelect;
    iLand=iSea;
    jLand=jSea+1;
  }
  if (DirSelect == 1 && SignSelect == -1) {
    // Down case in map_rivers.m
    iSea=iSelect;
    jSea=jSelect+1;
    iLand=iSea+1;
    jLand=jSea;
  }
  return {iSea-1, jSea-1, iLand-1, jLand-1};
}






/*
  Choice=0 means right
  Choice=1 means up
  Choice=2 means left
  Choice=3 means down
 */
ijdsInfo RetrieveIJDSarray(int const& eEtaSea, int const& eXiSea, int const& iChoice)
{
  int iSelect=-1, jSelect=-1, DirSelect=-1, SignSelect=-1;
  if (iChoice == 0) {
    iSelect=eEtaSea-1;
    jSelect=eXiSea-1;
    DirSelect=0;
    SignSelect=1;
  }
  if (iChoice == 1) {
    iSelect=eEtaSea-1;
    jSelect=eXiSea-1;
    DirSelect=1;
    SignSelect=1;
  }
  if (iChoice == 2) {
    iSelect=eEtaSea-1;
    jSelect=eXiSea;
    DirSelect=0;
    SignSelect=-1;
  }
  if (iChoice == 3) {
    iSelect=eEtaSea;
    jSelect=eXiSea-1;
    DirSelect=1;
    SignSelect=-1;
  }
  if (iChoice == -1) {
    std::cerr << "RetrieveIJDSarray INVALID ICHOICE iChoice=-1\n";
    return {0, -1, -1, -1, -1};
  }
  iSelect++; // increment conversion between matlab and C++ indexing.
  jSelect++;
  return {1, iSelect, jSelect, DirSelect, SignSelect};
}









void PlotRiverInformation(FullNamelist const& eFull)
{
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  //
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  std::string RiverFile=eBlPROC.ListStringValues.at("RiverFile");
  std::string SVGfile=eBlPROC.ListStringValues.at("SVGfile");
  //
  double SizeLine=eBlPLOT.ListDoubleValues.at("WidthLine");
  double SizeLineRiver=eBlPLOT.ListDoubleValues.at("WidthLineRiver");
  int SizeText=eBlPLOT.ListIntValues.at("SizeText");
  double ShiftLonText=eBlPLOT.ListDoubleValues.at("ShiftLonText");
  //
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  MyVector<double> ListETA_v = NC_Read1Dvariable(RiverFile, "river_Eposition");
  MyVector<double> ListXI_v  = NC_Read1Dvariable(RiverFile, "river_Xposition");
  MyVector<double> ListDir_v = NC_Read1Dvariable(RiverFile, "river_direction");
  MyMatrix<double> MatTransport = NC_Read2Dvariable(RiverFile, "river_transport");
  std::vector<double> ListRiverTime=NC_ReadTimeFromFile(RiverFile, "river_time");
  
  int nbRiver=ListETA_v.size();
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  std::vector<std::string> ListRiverName;
  try {
    netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::read);
    netCDF::NcGroupAtt RiversAtt=dataFile.getAtt("rivers");
    std::string strRiverName;
    RiversAtt.getValues(strRiverName);
    ListRiverName = STRING_Split(strRiverName, ",");
  }
  catch (...) {
    for (int i=0; i<nbRiver; i++)
      ListRiverName.push_back(IntToString(i));
  }
  //  std::cerr << "strRiverName = " << strRiverName << "\n";
  //
  int nbTime=MatTransport.rows();
  std::cerr << "MatTransport Rows=" << MatTransport.rows() << " Cols=" << MatTransport.cols() << "\n";
  std::vector<int> ListETA(nbRiver), ListXI(nbRiver), ListDir(nbRiver), ListSign(nbRiver);
  std::vector<double> ListAvgFlux(nbRiver);
  for (int iRiver=0; iRiver<nbRiver; iRiver++) {
    ListETA[iRiver] = ListETA_v(iRiver);
    ListXI[iRiver] = ListXI_v(iRiver);
    ListDir[iRiver] = ListDir_v(iRiver);
    double TheSum=0;
    for (int iTime=0; iTime<nbTime; iTime++)
      TheSum += MatTransport(iTime,iRiver);
    double TheAvgFlux=TheSum / double(nbTime);
    int eSign;
    if (TheSum > 0)
      eSign=1;
    else
      eSign=-1;
    ListSign[iRiver] = eSign;
    ListAvgFlux[iRiver] = TheAvgFlux;
  }
  //
  // Some standard definitions
  //
  MyMatrix<int> MatDir(4,2);
  MatDir(0,0)=0;
  MatDir(0,1)=0;
  MatDir(1,0)=1;
  MatDir(1,1)=0;
  MatDir(2,0)=1;
  MatDir(2,1)=1;
  MatDir(3,0)=0;
  MatDir(3,1)=1;
  CoordGridArrayFD RecPsi2=GRID_ExtendedPsiThi(GrdArr.GrdArrRho, GrdArr.GrdArrU, GrdArr.GrdArrV, GrdArr.GrdArrPsi);
  int eta_psi2=RecPsi2.LON.rows();
  int xi_psi2 =RecPsi2.LON.cols();
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> MatRadius = GetMatrixRadiusROMS(GrdArr);
  //
  // The list of entries to be printed.
  //
  std::vector<SVGgeneral> ListGeneral;
  //
  // Computing the coloring of the points
  //
  bool PolylineColoring=true;
  if (PolylineColoring) {
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	std::vector<coor> ListCoor(4);
	for (int u=0; u<4; u++) {
	  int i2=i + MatDir(u,0);
	  int j2=j + MatDir(u,1);
	  coor c{RecPsi2.LON(i2,j2), RecPsi2.LAT(i2,j2)};
	  ListCoor[u] = c;
	}
	std::string MarkerEnd="";
	std::string clip="full";
	double Size=0;
	std::vector<int> colorstroke{0,0,0};
	std::vector<int> colorfill;
	if (GrdArr.GrdArrRho.MSK(i,j) == 0)
	  colorfill={255,0,0};
	else
	  colorfill={0,0,255};
	SVGqualInfoPolyline eQual{colorfill, colorstroke, Size, MarkerEnd, clip};
	SVGpolyline ePolyline{ListCoor, eQual};
	ListGeneral.push_back(SVGgeneral(ePolyline));
      }
    std::cerr << "ListPolyline inserted\n";
  }
  //
  // Inserting the lines of the grid
  //
  auto InsertGridLine=[&](int const& i1, int const& j1, int const& i2, int const& j2) -> void {
    //    std::cerr << "InsertGridLine,  step 1, i1=" << i1 << " j1=" << j1 << " i2=" << i2 << " j2=" << j2 << "\n";
    coor ePt{RecPsi2.LON(i1, j1), RecPsi2.LAT(i1, j1)};
    //    std::cerr << "InsertGridLine,  step 2\n";
    coor fPt{RecPsi2.LON(i2, j2), RecPsi2.LAT(i2, j2)};
    //    std::cerr << "InsertGridLine,  step 3\n";
    std::string MarkerEnd="";
    std::string clip="";
    double Size=SizeLine;
    std::vector<int> color{0,0,0};
    SVGqualInfo eQual{color, Size, MarkerEnd, clip};
    SVGline eLine{ePt, fPt, eQual};
    ListGeneral.push_back(SVGgeneral(eLine));
  };
  for (int i=0; i<eta_psi2-1; i++)
    for (int j=0; j<xi_psi2; j++)
      InsertGridLine(i, j, i+1, j);
  std::cerr << "Block 1 inserted\n";
  for (int i=0; i<eta_psi2; i++)
    for (int j=0; j<xi_psi2-1; j++)
      InsertGridLine(i, j, i, j+1);
  std::cerr << "Block 2 inserted\n";
  //
  // Computing the river information
  //
  bool PlotRiver=true;
  if (PlotRiver) {
    for (int iRiver=0; iRiver<nbRiver; iRiver++) {
      int iSelect=ListETA[iRiver];
      int jSelect=ListXI[iRiver];
      int DirSelect=ListDir[iRiver];
      int SignSelect=ListSign[iRiver];
      ijSeaLand recIJSL = GetArrayIjSeaLand({1, iSelect, jSelect, DirSelect, SignSelect});
      for (int u=0; u<2; u++) {
	std::vector<int> color;
	int i, j;
	if (u == 0) {
	  color={255,0,0};
	  i=recIJSL.iSea;
	  j=recIJSL.jSea;
	}
	else {
	  color={0,0,255};
	  i=recIJSL.iLand;
	  j=recIJSL.jLand;
	}
	double coef=0.3;
	double eRad=coef * MatRadius(i,j);
	coor r{eRad, eRad};
	double lon=GrdArr.GrdArrRho.LON(i,j);
	double lat=GrdArr.GrdArrRho.LAT(i,j);
	coor c{lon, lat};
	std::string MarkerEnd = "";
	std::string clip = "";
	double Size=SizeLineRiver;
	SVGqualInfo eQual{color, Size, MarkerEnd, clip};
	SVGellipse eEll{c, r, eQual};
	ListGeneral.push_back(SVGgeneral(eEll));
      }
      int iSea=recIJSL.iSea;
      int jSea=recIJSL.jSea;
      int iLand=recIJSL.iLand;
      int jLand=recIJSL.jLand;
      coor ePt{GrdArr.GrdArrRho.LON(iSea, jSea), GrdArr.GrdArrRho.LAT(iSea, jSea)};
      coor fPt{GrdArr.GrdArrRho.LON(iLand, jLand), GrdArr.GrdArrRho.LAT(iLand, jLand)};
      int eMSKsea =GrdArr.GrdArrRho.MSK(iSea, jSea);
      int eMSKland=GrdArr.GrdArrRho.MSK(iLand, jLand);
      if (eMSKsea != 1 || eMSKland != 0) {
	std::cerr << "iRiver=" << iRiver << " has eMSKsea=" << eMSKsea << " eMSKland=" << eMSKland << "\n";
      }
      double lonSea=GrdArr.GrdArrRho.LON(iSea, jSea);
      double latSea=GrdArr.GrdArrRho.LAT(iSea, jSea);
      int eDEPsea =GrdArr.GrdArrRho.DEP(iSea, jSea);
      std::cerr << "iRiver=" << iRiver << " Sea(lon/lat/dep)=" << lonSea << " / " << latSea << " / " << eDEPsea << "\n";
      //
      {
	std::string MarkerEnd="";
	std::string clip="";
	double Size=SizeLineRiver;
	std::vector<int> color{255,255,0};
	SVGqualInfo eQual{color, Size, MarkerEnd, clip};
	SVGline eLine{ePt, fPt, eQual};
	ListGeneral.push_back(SVGgeneral(eLine));
      }
      //
      {
	double lonLnd=GrdArr.GrdArrRho.LON(iLand, jLand);
	double latLnd=GrdArr.GrdArrRho.LAT(iLand, jLand);
	double eLon=(lonSea + lonLnd)/2;
	double eLat=(latSea + latLnd)/2;
	coor point{eLon + ShiftLonText, eLat};
	std::vector<int> color{0,255,255};
	std::string str = ListRiverName[iRiver] + " : " + DoubleToString(ListAvgFlux[iRiver]);
	SVGtext eText{point, str, color, SizeText};
	ListGeneral.push_back(SVGgeneral(eText));
      }
    }
    std::cerr << "River lines inserted\n";
  }
  //
  // Plot the flux
  //
  bool PlotFlux=eBlPLOT.ListBoolValues.at("PlotFlux");
  std::cerr << "PlotFlux=" << PlotFlux << "\n";
  if (PlotFlux) {
    std::string strBEGTC=eBlPROC.ListStringValues.at("BEGTC");
    std::string strENDTC=eBlPROC.ListStringValues.at("ENDTC");
    double BeginTime=CT2MJD(strBEGTC);
    double EndTime  =CT2MJD(strENDTC);
    std::cerr << "BeginTime=" << BeginTime << " EndTime=" << EndTime << "\n";
    std::vector<int> ListIdx;
    int nbTime=ListRiverTime.size();
    for (int iTime=0; iTime<nbTime; iTime++) {
      double eTime = ListRiverTime[iTime];
      if (BeginTime <= eTime && eTime <= EndTime)
	ListIdx.push_back(iTime);
    }
    int nbCorr=ListIdx.size();
    std::cerr << "nbTime=" << nbTime << " nbCorr=" << nbCorr << "\n";
    MyVector<double> ListRiverTime_sel(nbCorr);
    for (int iCorr=0; iCorr<nbCorr; iCorr++) {
      int idx=ListIdx[iCorr];
      ListRiverTime_sel(iCorr) = ListRiverTime[idx];
    }
    std::cerr << "nbRiver=" << nbRiver << "\n";
    int nbLabel=eBlPLOT.ListIntValues.at("nbLabel");
    std::string StyleDate=eBlPLOT.ListStringValues.at("StyleDate");
    std::string StyleTitle=eBlPLOT.ListStringValues.at("StyleTitle");
    for (int iRiver=0; iRiver<nbRiver; iRiver++) {
      MyVector<double> ListTransport_sel(nbCorr);
      for (int iCorr=0; iCorr<nbCorr; iCorr++) {
	int idx=ListIdx[iCorr];
	ListTransport_sel(iCorr) = MatTransport(idx, iRiver);
      }
      double TheMax=ListTransport_sel.maxCoeff();
      DrawLinesArr eDrawArr;
      eDrawArr.DoTitle=true;
      if (StyleTitle == "style0")
	eDrawArr.TitleStr="River " + IntToString(iRiver) + " name=" + ListRiverName[iRiver];
      if (StyleTitle == "style1")
	eDrawArr.TitleStr=ListRiverName[iRiver];
      if (StyleTitle == "style2")
	eDrawArr.TitleStr=eBlPLOT.ListStringValues.at("TitleString");
      eDrawArr.IsTimeSeries=true;
      eDrawArr.PairComparison=false;
      eDrawArr.DoExplicitLabel=true;
      eDrawArr.nbLabel=nbLabel;
      eDrawArr.StyleDate=StyleDate;
      eDrawArr.VarName="River_flux";
      eDrawArr.TheMax=TheMax;
      eDrawArr.TheMin=0;
      eDrawArr.ListX=ListRiverTime_sel;
      eDrawArr.ListListVect={ListTransport_sel};
      eDrawArr.YAxisString = "River flux [m3/s]";
      //      eDrawArr.ListName_plot={"flux"};
      std::string FileName=ePerm.eDir + "FLUX_" + IntToString(iRiver);
      LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
    }
  }
  //
  // Computing the river information
  //
  SVGplotDescription SVGplot;
  SVGplot.ListGeneral = ListGeneral;
  SVGplot.FrameOption=1;
  SVGplot.height=400;
  SVGplot.width=400;
  SVGplot.scale_factor=1;
  SVGplot.add_offsetX=0;
  SVGplot.add_offsetY=0;
  SVGplot.RoundMethod=2;
  GeneralWriteSVGfile(SVGfile, SVGplot);
}











FullNamelist NAMELIST_GetStandard_ComputeRiverForcing_ROMS()
{
  std::map<std::string, SingleBlock> ListBlock;
  // INPUT
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListListStringValues1["ListRiverName"]={};
  ListStringValues1["RiverPrefix"]="unset";
  ListStringValues1["RiverSuffix"]="unset";
  ListStringValues1["GridFile"]="UNK";
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListStringValues1["RefTime"]="19680523.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListIntValues1["ARVD_N"]=-1;
  ListIntValues1["ARVD_Vtransform"]=-1;
  ListIntValues1["ARVD_Vstretching"]=-1;
  ListDoubleValues1["ARVD_Tcline"]=-1;
  ListDoubleValues1["ARVD_hc"]=-1;
  ListDoubleValues1["ARVD_theta_s"]=-1;
  ListDoubleValues1["ARVD_theta_b"]=-1;
  ListStringValues1["RiverFile"]="unset.nc";
  SingleBlock BlockINPUT;
  BlockINPUT.ListIntValues=ListIntValues1;
  BlockINPUT.ListBoolValues=ListBoolValues1;
  BlockINPUT.ListDoubleValues=ListDoubleValues1;
  BlockINPUT.ListListDoubleValues=ListListDoubleValues1;
  BlockINPUT.ListListIntValues=ListListIntValues1;
  BlockINPUT.ListStringValues=ListStringValues1;
  BlockINPUT.ListListStringValues=ListListStringValues1;
  ListBlock["INPUT"]=BlockINPUT;
  //
  return {ListBlock, "undefined"};
}


FullNamelist NAMELIST_RetrieveData()
{
  std::map<std::string, SingleBlock> ListBlock;
  // INPUT
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["RiverDescriptionFile"]="unset";
  ListListStringValues1["ListTimes"]={};
  ListIntValues1["StylePrint"]=1;
  SingleBlock BlockINPUT;
  BlockINPUT.ListIntValues=ListIntValues1;
  BlockINPUT.ListBoolValues=ListBoolValues1;
  BlockINPUT.ListDoubleValues=ListDoubleValues1;
  BlockINPUT.ListListDoubleValues=ListListDoubleValues1;
  BlockINPUT.ListListIntValues=ListListIntValues1;
  BlockINPUT.ListStringValues=ListStringValues1;
  BlockINPUT.ListListStringValues=ListListStringValues1;
  ListBlock["INPUT"]=BlockINPUT;
  //
  return {ListBlock, "undefined"};
}








double MonthlyInterpolation(std::vector<double> const& ListMonthly, double const& eTimeDay)
{
  std::vector<int> LDate = DATE_ConvertMjd2six(eTimeDay);
  int eYear=LDate[0];
  int eMonth=LDate[1];
  double eTimeDayMid = DATE_ConvertSix2mjd({eYear, eMonth, 15, 0, 0, 0});
  double eTimeDayMidSEQapprox;
  if (eTimeDay < eTimeDayMid) {
    eTimeDayMidSEQapprox = eTimeDayMid - 30;
  }
  else {
    eTimeDayMidSEQapprox = eTimeDayMid + 30;
  }
  std::vector<int> LDate2 = DATE_ConvertMjd2six(eTimeDayMidSEQapprox);
  int eYear2=LDate2[0];
  int eMonth2=LDate2[1];
  double eTimeDayMidSEQ = DATE_ConvertSix2mjd({eYear2, eMonth2, 15, 0, 0, 0});
  //
  double alpha2=(eTimeDay - eTimeDayMidSEQ) / (eTimeDayMid - eTimeDayMidSEQ);
  double alpha1=(eTimeDayMid - eTimeDay) / (eTimeDayMid - eTimeDayMidSEQ);
  double epsilon = 0.0001;
  if (alpha1 < -epsilon || alpha2 < -epsilon) {
    std::cerr << "--------------------------------\n";
    std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "eTimeDay=" << eTimeDay << " date=" << strPres << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMid=DATE_ConvertMjd2mystringPres(eTimeDayMid);
    std::cerr << "eTimeDayMid=" << eTimeDayMid << " date=" << strPresMid << "\n";
    std::cerr << "eYear=" << eYear << " eMonth=" << eMonth << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMidSEQ=DATE_ConvertMjd2mystringPres(eTimeDayMidSEQ);
    std::cerr << "eTimeDayMidSEQ=" << eTimeDayMidSEQ << " date=" << strPresMidSEQ << "\n";
    std::cerr << "eYear2=" << eYear2 << " eMonth2=" << eMonth2 << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMidSEQapprox=DATE_ConvertMjd2mystringPres(eTimeDayMidSEQapprox);
    std::cerr << "eTimeDayMidSEQapprox=" << eTimeDayMidSEQapprox << " date=" << strPresMidSEQapprox << "\n";
    std::cerr << "--------------------------------\n";
    std::cerr << "eTimeDayMid=" << eTimeDayMid << " eTimeDayMidSEQ=" << eTimeDayMidSEQ << "\n";
    std::cerr << "alpha1=" << alpha1 << " alpha2=" << alpha2 << "\n";
    std::cerr << "Error on alpha1 / alpha2\n";
    throw TerminalException{1};
  }
  double Flux1 = ListMonthly[eMonth2-1];
  double Flux2 = ListMonthly[eMonth-1];
  return alpha1 * Flux1 + alpha2 * Flux2;
}



struct TransTempSalt {
  double eTransport;
  double eTemp;
  double eSalt;
};



TransTempSalt RetrieveTTS(DescriptionRiver const& eDescRiv, double const& eTimeDay)
{
  double eTransport = 0, eTemp = 0, eSalt = 0;
  bool HasTransport=false, HasTemp=false, HasSalt=false;
  if (eDescRiv.TypeVaryingTransport == "PoFlux") {
    eTransport=InterpolateMeasurement(eDescRiv.ListPairTimeFlux, eTimeDay);
    HasTransport=true;
  }
  if (eDescRiv.TypeVaryingTransport == "ConstantFlux") {
    eTransport=eDescRiv.ConstantFlux;
    HasTransport=true;
  }
  if (eDescRiv.TypeVaryingTransport == "YearlyClimatology") {
    eTransport=MonthlyInterpolation(eDescRiv.ListMonthlyFlux, eTimeDay);
    HasTransport=true;
  }
  if (eDescRiv.TypeVaryingTransport == "RegularBurstOutflow") {
    double FrequencyDay=eDescRiv.FrequencyDay;
    double DurationHour = eDescRiv.DurationHour;
    double TotalFlux = eDescRiv.TotalFlux;  // In M3
    //
    double DurationSec = DurationHour * 3600;
    double DurationDay = DurationHour / 24;
    double InstantFlux = TotalFlux / DurationSec;
    //
    double eQuot=round(eTimeDay / FrequencyDay);
    double eDiff = eTimeDay - eQuot * FrequencyDay;
    if (fabs(eDiff) <= DurationDay / 2) {
      eTransport=InstantFlux;
    }
    else {
      eTransport=0;
    }
    HasTransport=true;
  }
  if (eDescRiv.TypeVaryingTemperature == "PoTemp") {
    eTemp=InterpolateMeasurement(eDescRiv.ListPairTimeTemp, eTimeDay);
    HasTemp=true;
  }
  if (eDescRiv.TypeVaryingTemperature == "YearlyClimatology") {
    eTemp=MonthlyInterpolation(eDescRiv.ListMonthlyTemp, eTimeDay);
    HasTemp=true;
  }
  if (eDescRiv.TypeVaryingTemperature == "Constant") {
    eTemp = eDescRiv.ConstantRiverTemperature;
    HasTemp=true;
  }
  if (eDescRiv.TypeVaryingSalinity == "Constant") {
    eSalt = eDescRiv.ConstantRiverSalinity;
    HasSalt=true;
  }
  //
  if (!HasTransport) {
    std::cerr << "We have HasTransport = " << HasTransport << "\n";
    std::cerr << "TypeVaryingTransport = " << eDescRiv.TypeVaryingTransport << "\n";
    std::cerr << "AllowedValues are : YearlyClimatology, RegularBurstOutflow\n";
    throw TerminalException{1};
  }
  if (!HasTemp) {
    std::cerr << "eDescRiv.TypeVaryingTemperature = " << eDescRiv.TypeVaryingTemperature << "\n";
    std::cerr << "We have HasTemp = " << HasTemp << "\n";
    throw TerminalException{1};
  }
  if (!HasSalt) {
    std::cerr << "eDescRiv.TypeVaryingSalinity = " << eDescRiv.TypeVaryingSalinity << "\n";
    std::cerr << "We have HasSalt = " << HasSalt << "\n";
    throw TerminalException{1};
  }
  eTransport *= eDescRiv.ConstantFactorFlux;
  return {eTransport, eTemp, eSalt};
}


int RIVER_ExtendedMask(MyMatrix<int> const& MSK_rho, int const& iEta, int const& iXi)
{
  int eta_rho=MSK_rho.rows();
  int xi_rho =MSK_rho.cols();
  if (iEta >= 0 && iEta < eta_rho && iXi >= 0 && iXi < xi_rho)
    return MSK_rho(iEta, iXi);
  else
    return 0;
}

double RIVER_AngleReduction(double const& TheAng)
{
  double pi = GetPI();
  double RetAng = TheAng;
  while(1) {
    std::cerr << "RetAng=" << RetAng << "\n";
    if (RetAng > pi)
      RetAng -= 2*pi;
    else {
      if (RetAng < -pi)
	RetAng += 2*pi;
      else
	break;
    }
  }
  return RetAng;
}



struct RecordAngleStatusRiver {
  Eigen::Tensor<int, 3> TotalListStatus;
  Eigen::Tensor<double,3> TotalListAngle;
  MyMatrix<int> StatusRULD;
};


/*
  StatusRULD is a matrix of size (eta_rho, xi_rho)
  It is 1 if a point is wet but has one of Right, Up, Left or Down point land.
  
 */
RecordAngleStatusRiver DetermineRiverPossibleCandidates(MyMatrix<int> const& MSK_rho, MyMatrix<double> const& ANG_rho)
{
  double pi = 3.141592653589793;
  int eta_rho=ANG_rho.rows();
  int xi_rho =ANG_rho.cols();
  Eigen::Tensor<int,3>    TotalListStatus(4, eta_rho, xi_rho);
  Eigen::Tensor<double,3> TotalListAngle(4, eta_rho, xi_rho);
  MyMatrix<int> StatusRULD(eta_rho, xi_rho);
  for (int iEta=0; iEta<eta_rho; iEta++)
    for (int iXi=0; iXi<xi_rho; iXi++) {
      int iRight=1, iLeft=1, iUp=1, iDown=1;
      double dRight=0, dLeft=0, dUp=0, dDown=0;
      if (MSK_rho(iEta, iXi) == 1) {
	int iEtaRight=iEta;
	int iXiRight=iXi-1;
	iRight = RIVER_ExtendedMask(MSK_rho, iEtaRight, iXiRight);
	dRight = ANG_rho(iEta, iXi);
	//
	int iEtaUp=iEta-1;
	int iXiUp=iXi;
	iUp = RIVER_ExtendedMask(MSK_rho, iEtaUp, iXiUp);
	dUp = ANG_rho(iEta, iXi) + pi/2;
	//
	int iEtaLeft=iEta;
	int iXiLeft=iXi+1;
	iLeft = RIVER_ExtendedMask(MSK_rho, iEtaLeft, iXiLeft);
	dLeft = ANG_rho(iEta, iXi) + pi;
	//
	int iEtaDown=iEta+1;
	int iXiDown=iXi;
	iDown = RIVER_ExtendedMask(MSK_rho, iEtaDown, iXiDown);
	dDown = ANG_rho(iEta, iXi) + 3*(pi/2);
      }
      TotalListStatus(0,iEta,iXi)=iRight;
      TotalListStatus(1,iEta,iXi)=iUp;
      TotalListStatus(2,iEta,iXi)=iLeft;
      TotalListStatus(3,iEta,iXi)=iDown;
      TotalListAngle(0,iEta,iXi)=dRight;
      TotalListAngle(1,iEta,iXi)=dUp;
      TotalListAngle(2,iEta,iXi)=dLeft;
      TotalListAngle(3,iEta,iXi)=dDown;
      StatusRULD(iEta,iXi) = 1 - iRight * iUp * iLeft * iDown;
    }
  int eProd=eta_rho * xi_rho;
  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << " eProd=" << eProd << "\n";
  std::cerr << "sum(MSK_rho)=" << MSK_rho.sum() << "\n";
  std::cerr << "sum(StatusRULD)=" << StatusRULD.sum() << " min/max=" << StatusRULD.minCoeff() << " / " << StatusRULD.maxCoeff() << "\n";
  for (int i=0; i<4; i++) {
    MyMatrix<int> ArrStatus = DimensionExtraction(TotalListStatus, 0, i);
    MyMatrix<double> ArrAngle = DimensionExtraction(TotalListAngle, 0, i);
    std::cerr << "i=" << i << " sum(ArrStatus)=" << ArrStatus.sum() << " sum(ArrAngle)=" << ArrAngle.sum() << "\n";
  }
  return {TotalListStatus, TotalListAngle, StatusRULD};
}


MyVector<double> RetrieveListOfWeight(MyVector<double> const& Zr, MyVector<double> const& Zw, DescriptionRiver const& eDescRiv)
{
  int N = Zr.size();
  std::cerr << "N = " << N << "\n";
  for (int i=0; i<N; i++) {
    std::cerr << "i=" << i << " Zr=" << Zr(i) << "\n";
  }
  int Np1 = Zw.size();
  std::cerr << "Np1 = " << Np1 << "\n";
  for (int i=0; i<Np1; i++) {
    std::cerr << "i=" << i << " Zw=" << Zw(i) << "\n";
  }
  MyVector<double> PreListWeight(N);
  if (eDescRiv.verticalShapeOption == "UniformUpperLayer") {
    for (int iS=0; iS<N; iS++) {
      double eDepLevel=Zr(iS);
      double eHeight=Zw(iS+1) - Zw(iS);
      double eWeight;
      if (eDepLevel > - eDescRiv.MaxDepth) {
	eWeight=eHeight;
      }
      else {
	eWeight=0;
      }
      PreListWeight(iS) = eWeight;
    }
    double eSum = PreListWeight.sum();
    return PreListWeight / eSum;
  }
  if (eDescRiv.verticalShapeOption == "FixedDepth") {
    for (int iS=0; iS<N; iS++)
      PreListWeight(iS) = 0;
    int iSfound = -1;
    for (int iS=0; iS<N; iS++) {
      double eDep=- eDescRiv.targetDepth;
      if (Zw(iS+1) > eDep && eDep > Zw(iS))
	iSfound = iS;
    }
    if (iSfound == -1) {
      std::cerr << "We did not find wanted depth\n";
      throw TerminalException{1};
    }
    PreListWeight(iSfound) = 1;
    return PreListWeight;
  }
  std::cerr << "No option chosen for the vertical shape\n";
  throw TerminalException{1};
}



void CreateRiverFile(FullNamelist const& eFull)
{
  SingleBlock eBlINPUT=eFull.ListBlock.at("INPUT");
  //
  // List of river names
  //
  std::vector<std::string> ListRiverName = eBlINPUT.ListListStringValues.at("ListRiverName");
  int nbRiver=ListRiverName.size();
  std::string RiverPrefix=eBlINPUT.ListStringValues.at("RiverPrefix");
  std::string RiverSuffix=eBlINPUT.ListStringValues.at("RiverSuffix");
  std::string RiverFile=eBlINPUT.ListStringValues.at("RiverFile");
  int N=eBlINPUT.ListIntValues.at("ARVD_N");
  //
  // Timings
  //
  std::string RefTimeStr=eBlINPUT.ListStringValues.at("RefTime");
  double RefTime=CT2MJD(RefTimeStr);
  std::string strBeginTime=eBlINPUT.ListStringValues.at("BEGTC");
  double BeginTime=CT2MJD(strBeginTime);
  std::string strEndTime=eBlINPUT.ListStringValues.at("ENDTC");
  double EndTime=CT2MJD(strEndTime);
  std::string UNITC=eBlINPUT.ListStringValues.at("UNITC");
  double DELTC=eBlINPUT.ListDoubleValues.at("DELTC");
  double DeltaTime = GetIntervalSize(DELTC, UNITC);
  //
  // Now reading the vertical stratification
  //
  int Vtransform=eBlINPUT.ListIntValues.at("ARVD_Vtransform");
  int Vstretching=eBlINPUT.ListIntValues.at("ARVD_Vstretching");
  double Tcline=eBlINPUT.ListDoubleValues.at("ARVD_Tcline");
  double hc=eBlINPUT.ListDoubleValues.at("ARVD_hc");
  double theta_s=eBlINPUT.ListDoubleValues.at("ARVD_theta_s");
  double theta_b=eBlINPUT.ListDoubleValues.at("ARVD_theta_b");
  ARVDtyp ARVD = ROMSgetARrayVerticalDescription(N, Vtransform, Vstretching, Tcline, hc, theta_s, theta_b);
  //
  // Now reading the input file for the rivers
  //
  std::vector<DescriptionRiver> ListDescriptionRiver(nbRiver);
  for (int iRiver=0; iRiver<nbRiver; iRiver++) {
    std::string eRiverName = ListRiverName[iRiver];
    std::string eRiverNameFull = RiverPrefix + eRiverName + RiverSuffix;
    std::cerr << "Before ReadRiverDescription iRiver=" << iRiver << "\n";
    ListDescriptionRiver[iRiver] = ReadRiverDescription(eRiverNameFull);
    std::cerr << " After ReadRiverDescription\n";
  }
  //
  // Now reading the grid arrays and related stuff
  //
  std::string GridFile = eBlINPUT.ListStringValues.at("GridFile");
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho=GrdArr.GrdArrRho.LON.cols();
  struct ijLL {
    int i;
    int j;
    double lon;
    double lat;
  };
  RecordAngleStatusRiver RecAngStatRiv = DetermineRiverPossibleCandidates(GrdArr.GrdArrRho.MSK, GrdArr.GrdArrRho.ANG);
  std::vector<ijLL> ListIJLLruld;
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++)
      if (RecAngStatRiv.StatusRULD(i,j) == 1) {
	double lon=GrdArr.GrdArrRho.LON(i,j);
	double lat=GrdArr.GrdArrRho.LAT(i,j);
	ListIJLLruld.push_back({i,j,lon,lat});
      }
  std::vector<ijLL> ListIJLLwet;
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++)
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	double lon=GrdArr.GrdArrRho.LON(i,j);
	double lat=GrdArr.GrdArrRho.LAT(i,j);
	ListIJLLwet.push_back({i,j,lon,lat});
      }  
  std::cerr << "sum(RecAngStatRiv.StatusRULD)=" << RecAngStatRiv.StatusRULD.sum() << "\n";
  auto GetNearest=[&](double const& lon, double const& lat, std::vector<ijLL> const& ListIJLL) -> ijLL {
    bool IsFirst=true;
    ijLL TheSel;
    double MinNorm=-1;
    for (auto & eEnt : ListIJLL) {
      double dist=GeodesicDistanceKM(eEnt.lon, eEnt.lat, lon, lat);
      if (IsFirst) {
	TheSel = eEnt;
	IsFirst=false;
	MinNorm=dist;
      }
      else {
	if (dist < MinNorm) {
	  MinNorm = dist;
	  TheSel = eEnt;
	}
      }
    }
    return TheSel;
  };
  MyMatrix<double> InflMatrixSph_rho=GRID_GetMaxRadiusInfluence_kernel(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT, "spherical");
  //
  // Now creating the list of signs, ETA, XI and direction
  //
  std::vector<int> ListSign, ListETA, ListXI, ListDirection, ListIRiver;
  std::vector<int> ListETAland, ListXIland, ListETAsea, ListXIsea;
  std::vector<double> ListDepArrival;
  std::vector<std::string> ListStringName;
  double pi=GetPI();
  auto GetIJDSriverCase=[&](double const& lon, double const& lat, double const& direction) -> ijdsInfo {
    double TheAngOrig = direction * (pi / 180);
    ijLL ijllNear = GetNearest(lon, lat, ListIJLLruld);
    double distance=GeodesicDistanceKM(ijllNear.lon, ijllNear.lat, lon, lat);
    int eEtaSea=ijllNear.i;
    int eXiSea =ijllNear.j;
    double dist1=InflMatrixSph_rho(eEtaSea, eXiSea);
    if (distance > dist1) {
      std::cerr << "GetIJDSriverCase distance=" << distance << " dist1=" << dist1 << " INVALID distance.\n";
      return {0, -1, -1, -1, -1};
    }
    int iChoice=-1;
    double minDeltaAng = 4545;
    for (int i=0; i<4; i++) {
      if (RecAngStatRiv.TotalListStatus(i,eEtaSea,eXiSea) == 0) {
	double TheAngNew=RecAngStatRiv.TotalListAngle(i,eEtaSea,eXiSea);
	double deltaAng=TheAngNew - TheAngOrig;
	double deltaAngRed=abs(RIVER_AngleReduction(deltaAng));
	if (deltaAngRed < minDeltaAng) {
	  iChoice=i;
	  minDeltaAng = deltaAngRed;
	}
      }
    }
    return RetrieveIJDSarray(eEtaSea, eXiSea, iChoice);
  };
  auto GetIJDSwetCase=[&](double const& lon, double const& lat, int const& ChoiceSelect) -> ijdsInfo {
    ijLL ijllNear = GetNearest(lon, lat, ListIJLLwet);
    double MinDist=GeodesicDistanceKM(ijllNear.lon, ijllNear.lat, lon, lat);
    int eEtaSea=ijllNear.i;
    int eXiSea =ijllNear.j;
    double distInfl=InflMatrixSph_rho(eEtaSea, eXiSea);
    if (MinDist > distInfl) {
      std::cerr << "GetIJDSwetCase : distInfl=" << distInfl << " MinDist=" << MinDist << " INVALID distance\n";
      return {0, -1, -1, -1, -1};
    }
    return RetrieveIJDSarray(eEtaSea, eXiSea, ChoiceSelect);
  };
  for (int iRiver=0; iRiver<nbRiver; iRiver++) {
    DescriptionRiver eDescRiv = ListDescriptionRiver[iRiver];
    //
    ijdsInfo RecordInfo;
    bool IsDone = false;
    if (eDescRiv.WScase == "River") {
      RecordInfo = GetIJDSriverCase(eDescRiv.lon, eDescRiv.lat, eDescRiv.direction);
      IsDone = true;
    }
    if (eDescRiv.WScase == "Wet") {
      RecordInfo = GetIJDSwetCase(eDescRiv.lon, eDescRiv.lat, eDescRiv.ChoiceSelect);
      IsDone = true;
    }
    if (eDescRiv.WScase == "Direct") {
      RecordInfo = {1, eDescRiv.iSelect, eDescRiv.jSelect, eDescRiv.DirSelect, eDescRiv.SignSelect};
      IsDone = true;
    }
    if (!IsDone) {
      std::cerr << "We failed to find matching algorithm\n";
      throw TerminalException{1};
    }
    std::cerr << "iRiver=" << iRiver << " name=" << eDescRiv.name << "\n";
    if (RecordInfo.eStatus == 0) {
      std::cerr << "  No valid choice found WScase = " << eDescRiv.WScase << "\n";
    }
    else {
      if (RecordInfo.DirSelect != 0 && RecordInfo.DirSelect != 1) {
	std::cerr << "DirSelect=" << RecordInfo.DirSelect << "\n";
	std::cerr << "DirSelect should be 0 or 1\n";
	throw TerminalException{1};
      }
      if (RecordInfo.SignSelect != -1 && RecordInfo.SignSelect != 1) {
	std::cerr << "SignSelect=" << RecordInfo.SignSelect << "\n";
	std::cerr << "SignSelect should be -1 or 1\n";
	throw TerminalException{1};
      }
      ijSeaLand recIJSL = GetArrayIjSeaLand(RecordInfo);
      ListETAsea.push_back(recIJSL.iSea);
      ListXIsea.push_back(recIJSL.jSea);
      ListETAland.push_back(recIJSL.iLand);
      ListXIland.push_back(recIJSL.jLand);
      double eDEP = GrdArr.GrdArrRho.DEP(recIJSL.iSea, recIJSL.jSea);
      ListDepArrival.push_back(eDEP);
      //
      ListETA.push_back(RecordInfo.iSelect);
      ListXI.push_back(RecordInfo.jSelect);
      ListDirection.push_back(RecordInfo.DirSelect);
      ListSign.push_back(RecordInfo.SignSelect);
      ListIRiver.push_back(iRiver);
      ListStringName.push_back(eDescRiv.name);
      std::cerr << "  Found river iS=" << RecordInfo.iSelect << " jS=" << RecordInfo.jSelect << " iD=" << RecordInfo.DirSelect << " iN=" << RecordInfo.SignSelect << "\n";
      std::cerr << "       Land(i,j)=" << recIJSL.iLand << " / " << recIJSL.jLand << " Sea(i,j)=" << recIJSL.iSea << " / " << recIJSL.jSea << " dep=" << eDEP << "\n";
    }
  }
  int nbRiverReal = ListIRiver.size();
  std::cerr << "nbRiver=" << nbRiver << " nbRiverReal=" << nbRiverReal << "\n";
  MyMatrix<int> MatrixIdx(N, nbRiverReal);
  int idx=0;
  for (int iN=0; iN<N; iN++)
    for (int iRiverReal=0; iRiverReal<nbRiverReal; iRiverReal++) {
      MatrixIdx(iN, iRiverReal) = idx;
      idx++;
    }
  //
  // Basic definition of the river file
  //
  netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  // Now dimensions
  netCDF::NcDim eDimRiver=dataFile.addDim("river", nbRiverReal);
  netCDF::NcDim eDimSvert=dataFile.addDim("s_rho", N);
  netCDF::NcDim eDimRiverTime=dataFile.addDim("river_time");
  netCDF::NcDim eDim19=dataFile.addDim("dateString", 19);
  // Now variables
  std::vector<std::string> LDim1{"river"};
  std::vector<std::string> LDim2{"s_rho", "river"};
  std::vector<std::string> LDim3{"river_time", "s_rho", "river"};
  std::vector<std::string> LDim4{"river_time"};
  std::vector<std::string> LDim5{"river_time", "dateString"};
  std::vector<std::string> LDim6{"river_time", "river"};
  netCDF::NcVar varRiver = dataFile.addVar("river", "double", LDim1);
  varRiver.putAtt("long_name", "river runoff identification number");
  varRiver.putAtt("units", "nondimensional");
  varRiver.putAtt("field", "river, scalar");
  netCDF::NcVar varEposition = dataFile.addVar("river_Eposition", "double", LDim1);
  varEposition.putAtt("long_name", "river ETA-position at RHO-points");
  varEposition.putAtt("units", "nondimensional");
  varEposition.putAtt("field", "river_Eposition, scalar");
  netCDF::NcVar varXposition = dataFile.addVar("river_Xposition", "double", LDim1);
  varXposition.putAtt("long_name", "river XI-position at RHO-points");
  varXposition.putAtt("units", "nondimensional");
  varXposition.putAtt("field", "river_Xposition, scalar");
  netCDF::NcVar varVshape = dataFile.addVar("river_Vshape", "double", LDim2);
  varVshape.putAtt("long_name", "river runoff mass transport vertical profile");
  varVshape.putAtt("units", "nondimensional");
  varVshape.putAtt("field", "river_Vshape, scalar");
  netCDF::NcVar varDirection = dataFile.addVar("river_direction", "double", LDim1);
  varDirection.putAtt("long_name", "river runoff direction");
  varDirection.putAtt("units", "nondimensional");
  varDirection.putAtt("field", "river_direction, scalar");
  netCDF::NcVar varFlag = dataFile.addVar("river_flag", "double", LDim1);
  varFlag.putAtt("long_name", "river runoff tracer flag");
  varFlag.putAtt("units", "nondimensional");
  varFlag.putAtt("option_0", "all tracers are off");
  varFlag.putAtt("option_1", "only temperature is on");
  varFlag.putAtt("option_2", "only salinity is on");
  varFlag.putAtt("option_3", "only both are on");
  varFlag.putAtt("field", "river_flag, scalar");
  netCDF::NcVar varSalt = dataFile.addVar("river_salt", "double", LDim3);
  varSalt.putAtt("long_name", "river runoff salinity");
  varSalt.putAtt("units", "PSU");
  varSalt.putAtt("field", "river_salt, scalar, series");
  netCDF::NcVar varTemp = dataFile.addVar("river_temp", "double", LDim3);
  varTemp.putAtt("long_name", "river runoff potential temperature");
  varTemp.putAtt("units", "Celsius");
  varTemp.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTime = dataFile.addVar("river_time", "double", LDim4);
  varRiverTime.putAtt("long_name", "river runoff time");
  std::string dateStr=DATE_ConvertMjd2mystringPres(RefTime);
  std::string attTime="days since " + dateStr;
  varRiverTime.putAtt("units", attTime);
  varRiverTime.putAtt("calendar", "gregorian");
  varRiverTime.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTimeStr = dataFile.addVar("river_time_str", "char", LDim5);
  netCDF::NcVar varTransport = dataFile.addVar("river_transport", "double", LDim6);
  varTransport.putAtt("long_name", "river runoff vertically integrated mass transport");
  varTransport.putAtt("units", "meter3 second-1");
  varTransport.putAtt("field", "river_transport, scalar, series");
  //
  // List of all names
  //
  std::string stringListNameRiver;
  for (int iRiverReal=0; iRiverReal<nbRiverReal; iRiverReal++) {
    if (iRiverReal > 0)
      stringListNameRiver += ",";
    stringListNameRiver += ListStringName[iRiverReal];
  }
  dataFile.putAtt("rivers", stringListNameRiver);
  //
  // Function for write downs
  //
  auto WriteDownNbRiver=[&](netCDF::NcVar & eVAR, std::vector<int> const& ListVal) -> void {
    double *A;
    A = new double[nbRiverReal];
    for (int i=0; i<nbRiverReal; i++)
      A[i] = double(ListVal[i]);
    eVAR.putVar(A);
    delete [] A;
  };
  //
  // Now easy definitions
  //
  std::vector<int> ListIdxRiver(nbRiverReal);
  for (int i=0; i<nbRiverReal; i++)
    ListIdxRiver[i] = i+1;
  WriteDownNbRiver(varRiver, ListIdxRiver);
  WriteDownNbRiver(varEposition, ListETA);
  WriteDownNbRiver(varXposition, ListXI);
  WriteDownNbRiver(varDirection, ListDirection);
  //
  // Now creating whether we assign salinity and/or temperature
  //
  std::vector<int> ListFlag(nbRiver);
  for (int iRiverReal=0; iRiverReal<nbRiverReal; iRiverReal++) {
    int iRiver=ListIRiver[iRiverReal];
    int eFlag=0;
    bool SetTemp=ListDescriptionRiver[iRiver].SetRiverTemperature;
    bool SetSalt   =ListDescriptionRiver[iRiver].SetRiverSalinity;
    if (SetTemp && SetSalt)
      eFlag=3;
    if (!SetTemp && SetSalt)
      eFlag=2;
    if (SetTemp && !SetSalt)
      eFlag=1;
    if (!SetTemp && !SetSalt)
      eFlag=0;
    std::cerr << "iRiverReal=" << iRiverReal << " iRiver=" << iRiver << " name=" << ListDescriptionRiver[iRiver].name << " SetTemp=" << SetTemp << " SetSalt=" << SetSalt << " eFlag=" << eFlag << "\n";
    ListFlag[iRiverReal] = eFlag;
  }
  WriteDownNbRiver(varFlag, ListFlag);
  //
  // Now we write down the vertical shape
  //
  double *Ashape;
  Ashape = new double[N * nbRiverReal];
  for (int iRiverReal=0; iRiverReal<nbRiverReal; iRiverReal++) {
    int iRiver=ListIRiver[iRiverReal];
    DescriptionRiver eDescRiv = ListDescriptionRiver[iRiver];
    double eDep=ListDepArrival[iRiverReal];
    double eZeta=0;
    MyVector<double> Zr_out=GetVertCoord_R(ARVD, eDep, eZeta);
    MyVector<double> Zw_out=GetVertCoord_W(ARVD, eDep, eZeta);
    MyVector<double> ListWeight = RetrieveListOfWeight(Zr_out, Zw_out, eDescRiv);
    for (int iS=0; iS<N; iS++) {
      int idx=MatrixIdx(iS, iRiverReal);
      Ashape[idx] = ListWeight(iS);
    }
  }
  varVshape.putVar(Ashape);
  delete [] Ashape;
  //
  // The time loop of creating the data
  //
  double CurrentTime=BeginTime;
  double epsilon = 0.00001;
  int pos=0;
  double *Asalt, *Atemp, *Atransport;
  Asalt = new double[N * nbRiverReal];
  Atemp = new double[N * nbRiverReal];
  Atransport = new double[nbRiverReal];
  while(1) {
    std::cerr << "pos=" << pos << "\n";
    // First writing the time
    std::string strPres=DATE_ConvertMjd2mystringPres(CurrentTime);
    std::vector<size_t> start1{size_t(pos)};
    std::vector<size_t> count1{1};
    double DiffTime=CurrentTime - RefTime;
    varRiverTime.putVar(start1, count1, &DiffTime);
    std::vector<size_t> start2{size_t(pos),0};
    std::vector<size_t> count2{1,19};
    varRiverTimeStr.putVar(start2, count2, strPres.c_str());
    //
    for (int iRiverReal=0; iRiverReal<nbRiverReal; iRiverReal++) {
      int iRiver = ListIRiver[iRiverReal];
      TransTempSalt eTTS = RetrieveTTS(ListDescriptionRiver[iRiver], CurrentTime);
      for (int iS=0; iS<N; iS++) {
	int idx=MatrixIdx(iS, iRiverReal);
	Atemp[idx] = eTTS.eTemp;
	Asalt[idx] = eTTS.eSalt;
      }
      Atransport[iRiverReal] = ListSign[iRiverReal] * eTTS.eTransport;
    }
    // Now writing your data
    std::vector<size_t> startTS{size_t(pos), 0, 0};
    std::vector<size_t> countTS{1, size_t(N), size_t(nbRiverReal)};
    varSalt.putVar(startTS, countTS, Asalt);
    varTemp.putVar(startTS, countTS, Atemp);
    std::vector<size_t> startTrans{size_t(pos), 0};
    std::vector<size_t> countTrans{1, size_t(nbRiverReal)};
    varTransport.putVar(startTrans, countTrans, Atransport);
    // Now maybe leaving
    CurrentTime += DeltaTime;
    pos++;
    if (CurrentTime > EndTime + epsilon)
      break;
  }
  delete [] Asalt;
  delete [] Atemp;
  delete [] Atransport;  
}







void MergeRiverFile(std::string const& RiverFile, std::vector<std::string> const& ListRiverFile)
{
  int nbFile = ListRiverFile.size();
  if (nbFile == 0) {
    std::cerr << "We have |ListRiverFile| = 0\n";
    throw TerminalException{1};
  }
  std::vector<int> ListNbRiver(nbFile);
  std::vector<int> ListNbTime(nbFile);
  std::vector<int> ListSrho(nbFile);
  for (int iFile=0; iFile<nbFile; iFile++) {
    std::string eFile=ListRiverFile[iFile];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    std::string strRiver="river";
    std::string strRiverTime="river_time";
    std::string strSrho="s_rho";
    ListNbRiver[iFile] = NC_ReadDimension(dataFile, strRiver);
    ListNbTime[iFile] = NC_ReadDimension(dataFile, strRiverTime);
    ListSrho[iFile] = NC_ReadDimension(dataFile, strSrho);
  }
  auto CheckIdentityDimension=[&](std::vector<int> const& ListVal) -> void {
    int eMin=VectorMin(ListVal);
    int eMax=VectorMax(ListVal);
    if (eMin != eMax) {
      std::cerr << "Wrong identity of the dimensions\n";
      std::cerr << "eMin=" << eMin << " eMax=" << eMax << "\n";
      throw TerminalException{1};
    }
  };
  CheckIdentityDimension(ListNbTime);
  CheckIdentityDimension(ListSrho);
  int nbRiver=VectorSum(ListNbRiver);
  int nbTime=ListNbTime[0];
  int s_rho=ListSrho[0];
  //
  // Reading time from Single river file
  //
  MyVector<double> ListTime;
  auto GetRefTime=[&]() -> double {
    std::string SingleRiverFile = ListRiverFile[0];
    netCDF::NcFile dataFile(SingleRiverFile, netCDF::NcFile::read);
    netCDF::NcVar data=dataFile.getVar("river_time");
    ListTime = NC_ReadVariable_data(data);
    return CT2MJD("19680523.000000");
  };
  double RefTime=GetRefTime();
  //
  // Basic definition of the river file
  //
  MyMatrix<int> MatrixIdx(s_rho, nbRiver);
  int idxM=0;
  for (int iS=0; iS<s_rho; iS++)
    for (int iRiver=0; iRiver<nbRiver; iRiver++) {
      MatrixIdx(iS, iRiver) = idxM;
      idxM++;
    }
  
  //
  // Basic definition of the river file
  //
  netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  // Now dimensions
  netCDF::NcDim eDimRiver=dataFile.addDim("river", nbRiver);
  netCDF::NcDim eDimSvert=dataFile.addDim("s_rho", s_rho);
  netCDF::NcDim eDimRiverTime=dataFile.addDim("river_time");
  netCDF::NcDim eDim19=dataFile.addDim("dateString", 19);
  // Now variables
  std::vector<std::string> LDim1{"river"};
  std::vector<std::string> LDim2{"s_rho", "river"};
  std::vector<std::string> LDim3{"river_time", "s_rho", "river"};
  std::vector<std::string> LDim4{"river_time"};
  std::vector<std::string> LDim5{"river_time", "dateString"};
  std::vector<std::string> LDim6{"river_time", "river"};
  netCDF::NcVar varRiver = dataFile.addVar("river", "double", LDim1);
  varRiver.putAtt("long_name", "river runoff identification number");
  varRiver.putAtt("units", "nondimensional");
  varRiver.putAtt("field", "river, scalar");
  netCDF::NcVar varEposition = dataFile.addVar("river_Eposition", "double", LDim1);
  varEposition.putAtt("long_name", "river ETA-position at RHO-points");
  varEposition.putAtt("units", "nondimensional");
  varEposition.putAtt("field", "river_Eposition, scalar");
  netCDF::NcVar varXposition = dataFile.addVar("river_Xposition", "double", LDim1);
  varXposition.putAtt("long_name", "river XI-position at RHO-points");
  varXposition.putAtt("units", "nondimensional");
  varXposition.putAtt("field", "river_Xposition, scalar");
  netCDF::NcVar varVshape = dataFile.addVar("river_Vshape", "double", LDim2);
  varVshape.putAtt("long_name", "river runoff mass transport vertical profile");
  varVshape.putAtt("units", "nondimensional");
  varVshape.putAtt("field", "river_Vshape, scalar");
  netCDF::NcVar varDirection = dataFile.addVar("river_direction", "double", LDim1);
  varDirection.putAtt("long_name", "river runoff direction");
  varDirection.putAtt("units", "nondimensional");
  varDirection.putAtt("field", "river_direction, scalar");
  netCDF::NcVar varFlag = dataFile.addVar("river_flag", "double", LDim1);
  varFlag.putAtt("long_name", "river runoff tracer flag");
  varFlag.putAtt("units", "nondimensional");
  varFlag.putAtt("option_0", "all tracers are off");
  varFlag.putAtt("option_1", "only temperature is on");
  varFlag.putAtt("option_2", "only salinity is on");
  varFlag.putAtt("option_3", "only both are on");
  varFlag.putAtt("field", "river_flag, scalar");
  netCDF::NcVar varSalt = dataFile.addVar("river_salt", "double", LDim3);
  varSalt.putAtt("long_name", "river runoff salinity");
  varSalt.putAtt("units", "PSU");
  varSalt.putAtt("field", "river_salt, scalar, series");
  netCDF::NcVar varTemp = dataFile.addVar("river_temp", "double", LDim3);
  varTemp.putAtt("long_name", "river runoff potential temperature");
  varTemp.putAtt("units", "Celsius");
  varTemp.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTime = dataFile.addVar("river_time", "double", LDim4);
  varRiverTime.putAtt("long_name", "river runoff time");
  std::string dateStr=DATE_ConvertMjd2mystringPres(RefTime);
  std::string attTime="days since " + dateStr;
  varRiverTime.putAtt("units", attTime);
  varRiverTime.putAtt("calendar", "gregorian");
  varRiverTime.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTimeStr = dataFile.addVar("river_time_str", "char", LDim5);
  netCDF::NcVar varTransport = dataFile.addVar("river_transport", "double", LDim6);
  varTransport.putAtt("long_name", "river runoff vertically integrated mass transport");
  varTransport.putAtt("units", "meter3 second-1");
  varTransport.putAtt("field", "river_transport, scalar, series");
  //
  // List of all names
  //
  std::string stringListNameRiver;
  for (int iFile=0; iFile<nbFile; iFile++) {
    std::string SingleRiverFile = ListRiverFile[0];
    netCDF::NcFile dataFile(SingleRiverFile, netCDF::NcFile::read);
    netCDF::NcGroupAtt RiversAtt=dataFile.getAtt("rivers");
    std::string strRiverName;
    RiversAtt.getValues(strRiverName);
    if (iFile>0)
      stringListNameRiver += ",";
    stringListNameRiver += strRiverName;
  }
  dataFile.putAtt("rivers", stringListNameRiver);
  //
  // Function for write downs
  //
  auto ExtractGlobalField=[&](std::string const& VarName) -> std::vector<double> {
    std::vector<double> retVal;
    int idx=0;
    for (int iFile=0; iFile<nbFile; iFile++) {
      std::string eFile=ListRiverFile[iFile];
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar data=dataFile.getVar(VarName);
      MyVector<double> ListVal = NC_ReadVariable_data(data);
      for (int iRiver=0; iRiver<ListNbRiver[iFile]; iRiver++) {
	retVal.push_back(ListVal(iRiver));
	idx++;
      }
    }
    return retVal;
  };
  auto WriteDownEntry=[&](netCDF::NcVar & eVAR, std::vector<double> ListVal) -> void {
    double *A;
    A = new double[nbRiver];
    for (int iRiver=0; iRiver<nbRiver; iRiver++)
      A[iRiver] = ListVal[iRiver];
    eVAR.putVar(A);
    delete [] A;
  };
  auto MergeAndWriteDownEntry=[&](netCDF::NcVar & eVAR,std::string const& VarName) -> void {
    std::vector<double> ListVal = ExtractGlobalField(VarName);
    WriteDownEntry(eVAR, ListVal);
  };
  //
  // Now easy definitions
  //
  std::vector<double> ListIdxRiver(nbRiver);
  for (int i=0; i<nbRiver; i++)
    ListIdxRiver[i] = i+1;
  WriteDownEntry(varRiver, ListIdxRiver);
  MergeAndWriteDownEntry(varEposition, "river_Eposition");
  MergeAndWriteDownEntry(varXposition, "river_Xposition");
  MergeAndWriteDownEntry(varDirection, "river_direction");
  MergeAndWriteDownEntry(varFlag, "river_flag");
  //
  // Now we write down the vertical shape
  //
  double *Ashape;
  Ashape = new double[s_rho * nbRiver];
  int idx2=0;
  for (int iFile=0; iFile<nbFile; iFile++) {
    std::string eFile=ListRiverFile[iFile];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar data=dataFile.getVar("river_Vshape");
    MyMatrix<double> eMat = NC_Read2Dvariable_data(data);
    int nbRiverLoc=ListNbRiver[iFile];
    for (int iRiverLoc=0; iRiverLoc<nbRiverLoc; iRiverLoc++) {
      for (int iS=0; iS<s_rho; iS++) {
	int idxB=MatrixIdx(iS, idx2);
	Ashape[idxB] = eMat(iS, iRiverLoc);
      }
      idx2++;
    }
  }
  varVshape.putVar(Ashape);
  delete [] Ashape;
  //
  // The time loop of creating the data
  //
  double *Asalt, *Atemp, *Atransport;
  Asalt = new double[s_rho * nbRiver];
  Atemp = new double[s_rho * nbRiver];
  Atransport = new double[nbRiver];
  for (int iTime=0; iTime<nbTime; iTime++) {
    // First writing the time
    double CurrentTime=ListTime(iTime);
    std::string strPres=DATE_ConvertMjd2mystringPres(RefTime + CurrentTime);
    std::vector<size_t> start1{size_t(iTime)};
    std::vector<size_t> count1{1};
    varRiverTime.putVar(start1, count1, &CurrentTime);
    std::vector<size_t> start2{size_t(iTime),0};
    std::vector<size_t> count2{1,19};
    varRiverTimeStr.putVar(start2, count2, strPres.c_str());
    //
    int idx3=0;
    for (int iFile=0; iFile<nbFile; iFile++) {
      int nbRiverLoc=ListNbRiver[iFile];
      std::vector<size_t> startTS{size_t(iTime), 0, 0};
      std::vector<size_t> countTS{1, size_t(s_rho), size_t(nbRiverLoc)};
      std::vector<size_t> startTrans{size_t(iTime), 0};
      std::vector<size_t> countTrans{1, size_t(nbRiverLoc)};
      //
      std::string eFile=ListRiverFile[iFile];
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar dataSalt=dataFile.getVar("river_salt");
      netCDF::NcVar dataTemp=dataFile.getVar("river_temp");
      netCDF::NcVar dataTrans=dataFile.getVar("river_transport");
      MyVector<double> ArrSalt = NC_ReadVariable_data_start_count(dataSalt, startTS, countTS);
      MyVector<double> ArrTemp = NC_ReadVariable_data_start_count(dataTemp, startTS, countTS);
      MyVector<double> ArrTrans = NC_ReadVariable_data_start_count(dataTrans, startTrans, countTrans);
      MyMatrix<int> MatrixIdxLoc(s_rho, nbRiverLoc);
      int idxM2=0;
      for (int iS=0; iS<s_rho; iS++)
	for (int iRiver=0; iRiver<nbRiverLoc; iRiver++) {
	  MatrixIdxLoc(iS, iRiver) = idxM2;
	  idxM2++;
	}
      for (int iRiverLoc=0; iRiverLoc<nbRiverLoc; iRiverLoc++) {
	for (int iS=0; iS<s_rho; iS++) {
	  int idx4=MatrixIdxLoc(iS, iRiverLoc);
	  Asalt[s_rho*idx3 + iS] = ArrSalt(idx4);
	  Atemp[s_rho*idx3 + iS] = ArrTemp(idx4);
	}
	Atransport[idx3] = ArrTrans(iRiverLoc);
      }
    }
    // Now writing your data
    std::vector<size_t> startTS{size_t(iTime), 0, 0};
    std::vector<size_t> countTS{1, size_t(s_rho), size_t(nbRiver)};
    varSalt.putVar(startTS, countTS, Asalt);
    varTemp.putVar(startTS, countTS, Atemp);
    std::vector<size_t> startTrans{size_t(iTime), 0};
    std::vector<size_t> countTrans{1, size_t(nbRiver)};
    varTransport.putVar(startTrans, countTrans, Atransport);
  }
  delete [] Asalt;
  delete [] Atemp;
  delete [] Atransport;  
}



void PrintRiverInformation(FullNamelist const& eFull)
{
  SingleBlock eBlINPUT=eFull.ListBlock.at("INPUT");
  std::string RiverDescriptionFile=eBlINPUT.ListStringValues.at("RiverDescriptionFile");
  std::vector<std::string> ListTimes=eBlINPUT.ListListStringValues.at("ListTimes");
  int StylePrint=eBlINPUT.ListIntValues.at("StylePrint");
  
  DescriptionRiver eDescRiv = ReadRiverDescription(RiverDescriptionFile);
  int nbTime=ListTimes.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  for (int iTime=0; iTime<nbTime; iTime++) {
    std::string eTimeStr = ListTimes[iTime];
    double eTime = CT2MJD(eTimeStr);
    std::string strPres=DATE_ConvertMjd2mystringPres(eTime);
    TransTempSalt eTTS = RetrieveTTS(eDescRiv, eTime);
    if (StylePrint == 1)
      std::cerr << "iTime=" << iTime << " date=" << strPres << " transport=" << eTTS.eTransport << " temp=" << eTTS.eTemp << " salt=" << eTTS.eSalt << "\n";
    if (StylePrint == 2)
      std::cerr << eTTS.eTransport << " " << eTTS.eTemp << " " << eTTS.eSalt << "\n";
    if (StylePrint == 3)
      std::cerr << eTTS.eTransport << "\n";
    if (StylePrint == 4)
      std::cerr << eTTS.eTemp << "\n";
  }
  
}




#endif
