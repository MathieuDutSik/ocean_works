#ifndef INCLUDE_TRANSECT_PLOT
#define INCLUDE_TRANSECT_PLOT

#include "NamelistExampleOcean.h"
#include "SphericalGeom.h"
#include "Interpolation.h"
#include "Model_interpolation.h"
#include "Model_grids.h"
#include "CommonFuncModel.h"
#include "Kernel_Transect.h"
#include "Plotting_fct.h"

FullNamelist NAMELIST_GetStandard_PlotBuoy()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListStringValues1["KindSelect"]="direct"; // possible values: direct, monthly, seasonal, yearly, specific
  ListDoubleValues1["TimeFrameDay"]=1;
  ListIntValues1["nbBlock"]=1;
  ListListStringValues1["ListMODELNAME"]={"unset MODELNAME in ListMODELNAME"};
  ListListStringValues1["ListGridFile"]={"unset GridFile in ListGridFile"};
  ListListStringValues1["ListHisPrefix"]={"ROMS_output_"};
  ListListStringValues1["ListRunName"]={"Explicit scheme"};
  ListStringValues1["PrefixBuoyFileName"]="./";
  ListListStringValues1["ListBuoyFileName"]={};
  ListListStringValues1["ListBuoyFileNamePlot"]={};
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListStringValues1["__NaturePlot"]="TIMESERIES";
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
  //  ListListStringValues1["ListSpecificTimes"]={};
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  ListBlock["PROC"]=BlockPROC;
  // PLOT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["nbLabel"]=5;
  ListIntValues2["nbLevelSpa"]=20;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["DoTitle"]=true;
  ListBoolValues2["DoColorBar"]=true;
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["VariableMin"]=false;
  ListBoolValues2["VariableMax"]=false;
  ListBoolValues2["PrintMMA"]=false;
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["DrawHorizVertLines"]=false;
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListBoolValues2["cnSmoothingOn"]=false;
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["FillLand"]=true;
  ListBoolValues2["DoTitleString"]=false;
  ListBoolValues2["DoExplicitLabel"]=false;
  ListBoolValues2["DoScatterPlot"]=false;
  ListBoolValues2["DoLinePlot"]=false;
  ListStringValues2["StyleDate"]="unsetstyledate";
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListStringValues2["LandPortr"]="Landscape";
  ListStringValues2["optStatStr"]="double";
  ListStringValues2["AnnotationText"]="unset";
  ListStringValues2["FileDirectNCLins"]="irrelevant";
  ListStringValues2["GridResolution"]="HighRes";
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnFillOn"]=true;
  ListBoolValues2["cnLinesOn"]=false;
  ListBoolValues2["cnLineLabelsOn"]=false;
  ListDoubleValues2["vcRefLengthF"]=0.2;
  ListDoubleValues2["AnnotationLon"]=-400;
  ListDoubleValues2["AnnotationLat"]=-400;
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // FILTER
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::vector<double>> ListListDoubleValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListDoubleValues3["MaximumLengthInterpolationIntervalSeconds"]=3600;
  ListDoubleValues3["MinInvalidatingValue"] = 0.02;
  ListDoubleValues3["MaxInvalidatingValue"] = 2.0;
  ListDoubleValues3["MaxAllowedDecrease"] = 1.3;
  ListListStringValues3["BoundSingle_var"]={};
  ListListDoubleValues3["BoundSingle_min"]={};
  ListListDoubleValues3["BoundSingle_max"]={};
  SingleBlock BlockFILTER;
  BlockFILTER.ListBoolValues=ListBoolValues3;
  BlockFILTER.ListDoubleValues=ListDoubleValues3;
  BlockFILTER.ListListDoubleValues=ListListDoubleValues3;
  BlockFILTER.ListListStringValues=ListListStringValues3;
  ListBlock["FILTER"]=BlockFILTER;
  // VARS
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, std::vector<std::string>> ListListStringValues4;
  ListListStringValues4["ListVarSynonymIn"]={};
  ListListStringValues4["ListVarSynonymOut"]={};
  std::vector<std::string> ListVarOut=GetAllPossibleVariables();
  for (auto& eVal : ListVarOut)
    ListBoolValues4[eVal]=false;
  SingleBlock BlockVARS;
  BlockVARS.ListBoolValues=ListBoolValues4;
  BlockVARS.ListListStringValues=ListListStringValues4;
  ListBlock["VARS"]=BlockVARS;
  // Final part
  return {std::move(ListBlock), "undefined"};
}



FullNamelist NAMELIST_GetStandard_MultipleVarPlot()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListStringValues1["KindSelect"]="direct"; // possible values: direct, monthly, seasonal, yearly, specific
  ListDoubleValues1["TimeFrameDay"]=1;
  ListIntValues1["nbBlock"]=1;
  ListListStringValues1["ListMODELNAME"]={"unset MODELNAME in ListMODELNAME"};
  ListListStringValues1["ListGridFile"]={"unset GridFile in ListGridFile"};
  ListListStringValues1["ListHisPrefix"]={"ROMS_output_"};
  ListListStringValues1["ListRunName"]={};
  ListListStringValues1["ListVarName"]={};
  //
  ListListDoubleValues1["ListPointLongitude"]={};
  ListListDoubleValues1["ListPointLatitude"]={};
  ListListStringValues1["ListPointName"]={};
  ListListStringValues1["ListAverageRegionLongitude"]={};
  ListListStringValues1["ListAverageRegionLatitude"]={};
  ListListStringValues1["ListAverageRegionName"]={};
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListStringValues1["__NaturePlot"]="TIMESERIES";
  ListBoolValues1["FirstCleanDirectory"]=true;
  ListBoolValues1["KeepNC_NCL"]=false;
  ListBoolValues1["InPlaceRun"]=false;
  ListBoolValues1["PrintDebugInfo"]=false;
  ListBoolValues1["OnlyCreateFiles"]=false;
  ListBoolValues1["DoLinePlot"]=true;
  ListBoolValues1["DoCsvFile"]=true;
  ListBoolValues1["DoSeasonMonthlyAverages"]=true;
  ListIntValues1["NPROC"]=1;
  ListStringValues1["Pcolor_method"]="ncl";
  ListStringValues1["Quiver_method"]="ncl";
  ListStringValues1["Lines_method"]="ncl";
  ListStringValues1["Scatter_method"]="ncl";
  //  ListListStringValues1["ListSpecificTimes"]={};
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  ListBlock["PROC"]=BlockPROC;
  // PLOT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["nbLabel"]=5;
  ListIntValues2["nbLevelSpa"]=20;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["DoTitle"]=true;
  ListBoolValues2["DoColorBar"]=true;
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["VariableMin"]=false;
  ListBoolValues2["VariableMax"]=false;
  ListDoubleValues2["SpecifiedMin"]=double(0);
  ListDoubleValues2["SpecifiedMax"]=double(10);
  ListBoolValues2["PrintMMA"]=false;
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListBoolValues2["cnSmoothingOn"]=false;
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["FillLand"]=true;
  ListBoolValues2["DoTitleString"]=false;
  ListBoolValues2["DoExplicitLabel"]=false;
  ListBoolValues2["DrawHorizVertLines"]=false;
  ListStringValues2["StyleDate"]="unsetstyledate";
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListStringValues2["LandPortr"]="Landscape";
  ListStringValues2["optStatStr"]="double";
  ListStringValues2["AnnotationText"]="unset";
  ListStringValues2["FileDirectNCLins"]="irrelevant";
  ListStringValues2["GridResolution"]="HighRes";
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnFillOn"]=true;
  ListBoolValues2["cnLinesOn"]=false;
  ListBoolValues2["cnLineLabelsOn"]=false;
  ListDoubleValues2["vcRefLengthF"]=0.2;
  ListDoubleValues2["AnnotationLon"]=-400;
  ListDoubleValues2["AnnotationLat"]=-400;
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // Final
  return {std::move(ListBlock), "undefined"};
}





PairLL ReadStationCoordinate(std::string const& eFile)
{
  MyVector<double> ListLon = NC_Read1Dvariable(eFile, "LONGITUDE");
  MyVector<double> ListLat = NC_Read1Dvariable(eFile, "LATITUDE");
  double eLon = ListLon(0);
  double eLat = ListLat(0);
  return {eLon, eLat};
}



void CheckingGridPointStatus(std::vector<GridArray> const& ListGrdArr, MyMatrix<double> const& ListXY)
{
  int nbGrid=ListGrdArr.size();
  int nbBuoy=ListXY.cols();
  std::vector<std::vector<bool>> ListListBelong(nbGrid);
  for (int iGrid=0; iGrid<nbGrid; iGrid++)
    ListListBelong[iGrid] = DetermineBelonging_ListXY(ListGrdArr[iGrid], ListXY);
  int nbWrongBuoy=0;
  for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
    int nbError=0;
    for (int iGrid=0; iGrid<nbGrid; iGrid++) {
      bool val = ListListBelong[iGrid][iBuoy];
      if (!val)
	nbError++;
    }
    if (nbError > 0) {
      std::cerr << "Error for buoy iBuoy=" << iBuoy << "\n";
      nbWrongBuoy++;
    }
  }
  if (nbWrongBuoy > 0) {
    std::cerr << "nbWrongBuoy = " << nbWrongBuoy << "\n";
    std::cerr << "Several buoys are not contained in grids. Cannot work\n";
    throw TerminalException{1};
  }
}




void BUOY_Plot(FullNamelist const& eFull)
{
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  //
  // Reading grid arrays and the like
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::vector<std::string> ListModelName=eBlPROC.ListListStringValues.at("ListMODELNAME");
  std::vector<std::string> ListGridFile=eBlPROC.ListListStringValues.at("ListGridFile");
  std::vector<std::string> ListHisPrefix=eBlPROC.ListListStringValues.at("ListHisPrefix");
  std::vector<std::string> ListRunName=eBlPROC.ListListStringValues.at("ListRunName");
  std::string PicPrefix = eBlPROC.ListStringValues.at("PicPrefix");
  std::vector<std::string> ListRunNameExt{"buoy"};
  for (auto & eVal : ListRunName)
    ListRunNameExt.push_back(eVal);
  int nbGrid=ListGridFile.size();
  size_t nbGrid_t=ListGridFile.size();
  std::cerr << "nbGrid=" << nbGrid << "\n";
  if (nbGrid_t != ListHisPrefix.size() || nbGrid_t != ListModelName.size()) {
    std::cerr << "Error for the length of\n";
    std::cerr << "ListMODELNAME, ListGridFile, ListHisPrefix\n";
    std::cerr << "which should be of the same length\n";
    throw TerminalException{1};
  }
  std::vector<GridArray> ListGrdArr;
  std::vector<TotalArrGetData> ListTotalArr;
  for (int iGrid=0; iGrid<nbGrid; iGrid++) {
    std::cerr << "iGrid=" << iGrid << " / " << nbGrid << "\n";
    std::string eModelName=ListModelName[iGrid];
    std::string GridFile=ListGridFile[iGrid];
    std::string HisPrefix=ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    ListGrdArr.push_back(GrdArr);
    TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
    ListTotalArr.push_back(TotalArr);
  }
  //
  // Reading time information
  //
  std::string strBEGTC=eBlPROC.ListStringValues.at("BEGTC");
  std::string strENDTC=eBlPROC.ListStringValues.at("ENDTC");
  double BeginTime=CT2MJD(strBEGTC);
  double EndTime  =CT2MJD(strENDTC);
  double DELTC=eBlPROC.ListDoubleValues.at("DELTC");
  std::string UNITC=eBlPROC.ListStringValues.at("UNITC");
  double DeltaInterval=GetIntervalSize(DELTC, UNITC);
  std::vector<double> ListTime = GetIntervalFLD(BeginTime, EndTime, DeltaInterval);
  int nbTime = ListTime.size();
  MyVector<double> ListTime_vect(nbTime);
  for (int iTime=0; iTime<nbTime; iTime++)
    ListTime_vect(iTime) = ListTime[iTime];
  //
  // Reading filter information.
  //
  SingleBlock eBlockFILTER=eFull.ListBlock.at("FILTER");
  double MaximumLengthInterpolationIntervalSeconds=eBlockFILTER.ListDoubleValues.at("MaximumLengthInterpolationIntervalSeconds");
  double MaximumLengthInterpolationIntervalDay = MaximumLengthInterpolationIntervalSeconds / double(86400);
  double MinInvalidatingValue = eBlockFILTER.ListDoubleValues.at("MinInvalidatingValue");
  double MaxInvalidatingValue = eBlockFILTER.ListDoubleValues.at("MaxInvalidatingValue");
  std::cerr << "MinInvalidatingValue = " << MinInvalidatingValue << "\n";
  std::cerr << "MaxInvalidatingValue = " << MaxInvalidatingValue << "\n";
  std::vector<std::string> BoundSingle_var = eBlockFILTER.ListListStringValues.at("BoundSingle_var");
  std::vector<double> BoundSingle_min = eBlockFILTER.ListListDoubleValues.at("BoundSingle_min");
  std::vector<double> BoundSingle_max = eBlockFILTER.ListListDoubleValues.at("BoundSingle_max");
  if (BoundSingle_var.size() != BoundSingle_min.size() || BoundSingle_var.size() != BoundSingle_max.size()) {
    std::cerr << "Wrong sizes for BoundSingle\n";
    throw TerminalException{1};
  }
  double MaxAllowedDecrease = eBlockFILTER.ListDoubleValues.at("MaxAllowedDecrease");

  //
  // Reading variable information.
  //
  SingleBlock eBlockVAR=eFull.ListBlock.at("VARS");
  std::vector<std::string> ListVarOut=ExtractMatchingBool(eBlockVAR);
  std::vector<std::string> ListVarSynonymIn  = eBlockVAR.ListListStringValues.at("ListVarSynonymIn");
  std::vector<std::string> ListVarSynonymOut = eBlockVAR.ListListStringValues.at("ListVarSynonymOut");
  auto MapSynonym=[&](std::string const& eFile, std::string const& eVARin) -> std::string {
    int len=ListVarSynonymIn.size();
    std::vector<std::string> ListMatch;
    for (int i=0; i<len; i++) {
      if (ListVarSynonymIn[i] == eVARin) {
	std::string eVARout = ListVarSynonymOut[i];
	if (NC_IsVar(eFile, eVARout))
	  ListMatch.push_back(eVARout);
      }
    }
    if (ListMatch.size() != 1) {
      std::cerr << "|ListMatch|=" << ListMatch.size() << "\n";
      std::cerr << "should be exactly 1\n";
      std::cerr << "eFile = " << eFile << "\n";
      std::cerr << "Failed to find matching synonym for eVARin = " << eVARin << "\n";
      throw TerminalException{1};
    }
    return ListMatch[0];
  };
  //
  // Reading buoy location and computing interpolation arrays
  //
  std::vector<std::string> ListBuoyFileName=eBlPROC.ListListStringValues.at("ListBuoyFileName");
  std::vector<std::string> ListBuoyFileNamePlot=eBlPROC.ListListStringValues.at("ListBuoyFileNamePlot");
  if (ListBuoyFileName.size() != ListBuoyFileNamePlot.size()) {
    std::cerr << "ListBuoyFileName does not have the same size as ListBuoyFileNamePlot\n";
    throw TerminalException{1};
  }
  std::string PrefixBuoyFileName=eBlPROC.ListStringValues.at("PrefixBuoyFileName");
  int nbBuoy = ListBuoyFileName.size();
  MyMatrix<double> ListXY(2,nbBuoy);
  for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
    std::string eBuoyFileName = PrefixBuoyFileName + ListBuoyFileName[iBuoy];
    PairLL ePair = ReadStationCoordinate(eBuoyFileName);
    ListXY(0,iBuoy) = ePair.eLon;
    ListXY(1,iBuoy) = ePair.eLat;
  }
  CheckingGridPointStatus(ListGrdArr, ListXY);
  std::vector<SingleArrayInterpolation> ListRec(nbGrid);
  for (int iGrid=0; iGrid<nbGrid; iGrid++)
    ListRec[iGrid] = ComputeArrayInterpolation_ListXY(ListGrdArr[iGrid], ListXY);
  //
  // Reading and showing times
  //
  std::vector<int> NbCharBuoyName(nbBuoy);
  for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
    std::string eNameRed = ListBuoyFileName[iBuoy];
    int nbChar=eNameRed.size();
    NbCharBuoyName[iBuoy] = nbChar;
  }
  int eMaxNbChar = VectorMax(NbCharBuoyName);
  for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
    std::string eNameRed = ListBuoyFileName[iBuoy];
    std::string eBuoyFileName = PrefixBuoyFileName + eNameRed;
    std::vector<double> ListTimeBuoy = NC_ReadTimeFromFile(eBuoyFileName, "TIME");
    double minTime=VectorMin(ListTimeBuoy);
    double maxTime=VectorMax(ListTimeBuoy);
    int diffNb=eMaxNbChar - NbCharBuoyName[iBuoy];
    std::string minTimePres=DATE_ConvertMjd2mystringPres(minTime);
    std::string maxTimePres=DATE_ConvertMjd2mystringPres(maxTime);
    std::cerr << "iBuoy=" << iBuoy << " filename=" << eNameRed;
    for (int u=0; u<diffNb; u++)
      std::cerr << " ";
    std::cerr << " minTimePres=" << minTimePres << " maxTimePres=" << maxTimePres << "\n";
  }
  //
  // Reading full data set
  //
  bool VariableMin=eBlPLOT.ListBoolValues.at("VariableMin");
  bool VariableMax=eBlPLOT.ListBoolValues.at("VariableMax");
  bool DoScatterPlot=eBlPLOT.ListBoolValues.at("DoScatterPlot");
  bool DoLinePlot=eBlPLOT.ListBoolValues.at("DoLinePlot");
  bool DrawHorizVertLines=eBlPLOT.ListBoolValues.at("DrawHorizVertLines");
  bool DoExplicitLabel=eBlPLOT.ListBoolValues.at("DoExplicitLabel");
  int nbLabel=eBlPLOT.ListIntValues.at("nbLabel");
  std::string StyleDate = eBlPLOT.ListStringValues.at("StyleDate");
  PlotBound ePlotBound=ReadPlotBound(eFull);
  double MissingValueFinal=-9999;
  std::string StatFile = PicPrefix + "Statistics.txt";
  std::ofstream os(StatFile);
  for (auto & eVarName : ListVarOut) {
    std::cerr << "eVarName = " << eVarName << "\n";
    std::vector<std::vector<MyVector<double>>> ListListVect(nbBuoy);
    for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
      std::string eBuoyFileName = PrefixBuoyFileName + ListBuoyFileName[iBuoy];
      std::string eVarCall = MapSynonym(eBuoyFileName, eVarName);
      std::cerr << "eBuoyFileName = " << eBuoyFileName << "    eVarCall = " << eVarCall << "\n";
      std::vector<double> ListTimeBuoy = NC_ReadTimeFromFile(eBuoyFileName, "TIME");
      std::cerr << "ListTimeBuoy(min/max)=" << VectorMin(ListTimeBuoy) << " / " << VectorMax(ListTimeBuoy) << "\n";
      MyVector<double> ListVal;
      std::vector<size_t> LDim = NC_ReadVariable_listdim_file(eBuoyFileName, eVarCall);
      if (LDim.size() == 2) {
	MyMatrix<double> eMat = NC_Read2Dvariable(eBuoyFileName, eVarCall);
	int nbRow=eMat.rows();
	int nbCol=eMat.cols();
	std::cerr << "eMat rows = " << eMat.rows() << "    cols = " << eMat.cols() << "\n";
	int iColFind=-1;
	int nbFind=0;
	for (int iCol=0; iCol<nbCol; iCol++) {
	  MyVector<double> V(nbRow);
	  for (int iRow=0; iRow<nbRow; iRow++)
	    V(iRow) = eMat(iRow, iCol);
	  if (V.minCoeff() != V.maxCoeff()) {
	    iColFind = iCol;
	    nbFind++;
	  }
	}
	if (nbFind != 1) {
	  std::cerr << "iColFind=" << iColFind << " nbFind=" << nbFind << "\n";
	  std::cerr << "Our heuristic failed\n";
	  throw TerminalException{1};
	}
	MyVector<double> V(nbRow);
	for (int iRow=0; iRow<nbRow; iRow++)
	  V(iRow) = eMat(iRow, iColFind);
	ListVal = V;
      } else {
	if (LDim.size() == 1) {
	  ListVal = NC_Read1Dvariable(eBuoyFileName, eVarCall);
	} else {
	  std::cerr << "Wrong dimension in LDim\n";
	  throw TerminalException{1};
	}
      }
      std::cerr << "|ListTimeBuoy| = " << ListTimeBuoy.size() << "     |ListVal| = " << ListVal.size() << "\n";
      if (int(ListTimeBuoy.size()) != int(ListVal.size())) {
	std::cerr << "ListTimeBuoy and ListVal are of different size\n";
	throw TerminalException{1};
      }
      std::cerr << "ListVal(min/max)=" << ListVal.minCoeff() << " / " << ListVal.maxCoeff() << "\n";
      MyVector<double> ListInterpVal(nbTime);
      std::cerr << "Before computing ListInterpVal\n";
      double MissingValueDetect=-999;
      std::vector<int> ListITime;
      for (int iTime=0; iTime<nbTime; iTime++) {
	double eTimeDay=ListTime[iTime];
	//	std::cerr << "iTime=" << iTime << " / " << nbTime << " eTimeDay=" << eTimeDay << "\n";
	InterpInfo eInterp = GetTimeInterpolationInfo(ListTimeBuoy, eTimeDay);
	//	std::cerr << "eInterp.iTimeLow=" << eInterp.iTimeLow << " iTimeUpp=" << eInterp.iTimeUpp << "\n";
	double deltaTime=ListTimeBuoy[eInterp.iTimeUpp] - ListTimeBuoy[eInterp.iTimeLow];
	std::cerr << "iTime=" << iTime << " deltaTime=" << deltaTime << "\n";
	double val = MissingValueFinal;
	if (deltaTime < MaximumLengthInterpolationIntervalDay) {
	  std::cerr << "Matching in time\n";
	  if (ListVal(eInterp.iTimeLow) > MissingValueDetect && ListVal(eInterp.iTimeUpp) > MissingValueDetect) {
	    std::cerr << "Pass the MissingValue test\n";
	    double eVal = eInterp.alphaLow * ListVal(eInterp.iTimeLow) + eInterp.alphaUpp * ListVal(eInterp.iTimeUpp);
	    if (MinInvalidatingValue < eVal && eVal < MaxInvalidatingValue) {
	      std::cerr << "Matching in MinMax invalidating\n";
	      val = eVal;
	      ListITime.push_back(iTime);
	    }
	  }
	}
	ListInterpVal(iTime) = val;
      }
      int nbTimeSel=ListITime.size();
      for (int uTimeSel=0; uTimeSel<nbTimeSel; uTimeSel++) {
	int eTimeSel=ListITime[uTimeSel];
	bool IsCorrect=true;
	if (uTimeSel > 0) {
	  int eTimePrev=ListITime[uTimeSel-1];
	  if (ListInterpVal(eTimePrev) < ListInterpVal(eTimeSel) - MaxAllowedDecrease)
	    IsCorrect=false;
	}
	if (uTimeSel < nbTimeSel-1) {
	  int eTimeNext=ListITime[uTimeSel+1];
	  if (ListInterpVal(eTimeNext) < ListInterpVal(eTimeSel) - MaxAllowedDecrease)
	    IsCorrect=false;
	}
	if (!IsCorrect)
	  ListInterpVal(eTimeSel) = MissingValueFinal;
      }
      std::cerr << " After computing ListInterpVal\n";
      std::vector<MyVector<double>> ListVect(nbGrid+1);
      ListVect[0] = ListInterpVal;
      MyVector<double> eVect = ZeroVector<double>(nbTime);
      for (int iGrid=0; iGrid<nbGrid; iGrid++)
	ListVect[iGrid+1] = eVect;
      ListListVect[iBuoy] = ListVect;
    }
    std::cerr << "After reading of the buoy data\n";
    RecSymbolic RecS;
    for (int iGrid=0; iGrid<nbGrid; iGrid++)
      for (int iTime=0; iTime<nbTime; iTime++) {
	double eTimeDay = ListTime[iTime];
	VarQuery eQuery;
	eQuery.eTimeDay = eTimeDay;
	eQuery.iTime = -1;
	eQuery.NatureQuery = "instant";
	eQuery.TimeFrameDay = -1;
	eQuery.typeQuery = "unset";
	std::cerr << "iTime=" << iTime << " / " << nbTime << "\n";
	RecVar eRecVar = ModelSpecificVarSpecificTimeGeneral(ListTotalArr[iGrid], eVarName, eQuery, ePlotBound);
	RecS = eRecVar.RecS;
	RecVar fRecVar=INTERPOL_SingleRecVarInterpolation(ListRec[iGrid], eRecVar);
	for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
	  double eVal = fRecVar.F(iBuoy,0);
	  ListListVect[iBuoy][iGrid+1](iTime) = eVal;
	}
      }
    int len=BoundSingle_var.size();
    for (int u=0; u<len; u++)
      if (BoundSingle_var[u] == eVarName) {
	RecS.minval = BoundSingle_min[u];
	RecS.maxval = BoundSingle_max[u];
      }
    std::string eUnit = RecS.Unit;
    std::cerr << "After reading of the model data\n";
    for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
      std::string eBuoyFileName = ListBuoyFileName[iBuoy];
      std::string eBuoyFileNamePlot = ListBuoyFileNamePlot[iBuoy];
      int nbBlock = eBlPROC.ListIntValues.at("nbBlock");
      std::vector<int> ListPos = DivideListPosition(nbTime, nbBlock);
      std::cerr << "-------- LIST VAL iBuoy = " << iBuoy << " ----------\n";
      double minDeltaTime=100000, maxDeltaTime=0;
      for (int iTime=0; iTime<nbTime; iTime++) {
	std::cerr << "iTime=" << iTime << " / " << nbTime;
	//
	double deltaTime=0;
	if (iTime == 0) {
	  deltaTime = 0;
	} else {
	  deltaTime = ListTime_vect(iTime) - ListTime_vect(iTime-1);
	  if (deltaTime > maxDeltaTime)
	    maxDeltaTime = deltaTime;
	  if (deltaTime < minDeltaTime)
	    minDeltaTime = deltaTime;
	}
	std::cerr << " deltaT=" << deltaTime;
	//
	std::string strPres=DATE_ConvertMjd2mystringPres(ListTime_vect(iTime));
	std::cerr << " date=" << strPres;
	//
	std::cerr << " V =";
	for (int iGrid=0; iGrid<nbGrid; iGrid++)
	  std::cerr << " " << ListListVect[iBuoy][iGrid](iTime);
	std::cerr << "\n";
      }
      std::cerr << "deltaTime min=" << minDeltaTime << " max=" << maxDeltaTime << "\n";
      if (DoLinePlot) {
	for (int iBlock=0; iBlock<nbBlock; iBlock++) {
	  DrawLinesArr eDrawArr;
	  eDrawArr.DoTitle=true;
	  eDrawArr.TitleStr="Time series of " + eVarName + " for buoy " + eBuoyFileNamePlot;
	  eDrawArr.IsTimeSeries=true;
	  eDrawArr.PairComparison=false;
	  eDrawArr.nbLabel=nbLabel;
	  eDrawArr.DoExplicitLabel=DoExplicitLabel;
          eDrawArr.DrawHorizVertLines=DrawHorizVertLines;
	  eDrawArr.StyleDate=StyleDate;
	  eDrawArr.VarName=eVarName + "_" + std::to_string(iBuoy+1) + "_" + std::to_string(iBlock);
	  eDrawArr.ListName_plot=ListRunNameExt;
	  eDrawArr.YAxisString=RecS.VarName2 + "(" + RecS.Unit + ")";
	  //
	  // Determination of dimension variable
	  //
	  int pos1 = ListPos[iBlock];
	  int pos2 = ListPos[iBlock+1];
	  int len = pos2 - pos1;
	  MyVector<double> ListTime_SEL(len);
	  for (int pos=pos1; pos<pos2; pos++)
	    ListTime_SEL(pos - pos1) = ListTime_vect(pos);
	  std::vector<MyVector<double>> ListVect_SEL(nbGrid+1);
	  for (int iGrid=0; iGrid<=nbGrid; iGrid++) {
	    MyVector<double> eVect(len);
	    for (int pos=pos1; pos<pos2; pos++) {
	      double eVal = ListListVect[iBuoy][iGrid](pos);
	      eVect(pos - pos1) = eVal;
	    }
	    ListVect_SEL[iGrid] = eVect;
	    double maxVal = eVect.maxCoeff();
	    std::cerr << "iGrid=" << iGrid << " maxVal=" << maxVal << "\n";
	  }
	  std::cerr << "iBuoy=" << iBuoy << " iBlock=" << iBlock << "\n";
	  for (int u=0; u<len; u++) {
	    std::cerr << "u=" << u << " / " << len << "   V =";
	    for (int iGrid=0; iGrid<=nbGrid; iGrid++)
	      std::cerr << " " << ListVect_SEL[iGrid](u);
	    std::cerr << "\n";
	  }
	  eDrawArr.ListX = ListTime_SEL;
	  eDrawArr.ListListVect = ListVect_SEL;
	  //
	  // Determinantion of min/max ranges
	  //
	  double TheMax, TheMin;
	  if (VariableMax) {
	    TheMax = - 1000000;
	    for (auto & eVect : eDrawArr.ListListVect)
	      TheMax = std::max(TheMax, eVect.maxCoeff());
	  } else {
	    TheMax=RecS.maxval;
	  }
	  if (VariableMin) {
	    TheMin = 1000000;
	    for (auto & eVect : eDrawArr.ListListVect)
	      TheMin = std::min(TheMin, eVect.minCoeff());
	  } else {
	    TheMin=RecS.minval;
	  }
	  std::cerr << "TheMin=" << TheMin << " TheMax=" << TheMax << "\n";
	  eDrawArr.TheMax=TheMax;
	  eDrawArr.TheMin=TheMin;
	  //
	  std::string FileName=ePerm.eDir + "TimeSeries_iBuoy" + std::to_string(iBuoy+1) + "_iBlock" + std::to_string(iBlock) + "_" + eVarName;
	  LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
	}
      }
      //
      double TheMax_sc, TheMin_sc;
      if (VariableMax) {
	TheMax_sc = - 1000000;
	for (auto & eVect : ListListVect[iBuoy])
	  TheMax_sc = std::max(TheMax_sc, eVect.maxCoeff());
      } else {
	TheMax_sc = RecS.maxval;
      }
      if (VariableMin) {
	TheMin_sc = 1000000;
	for (auto & eVect : ListListVect[iBuoy])
	  TheMin_sc = std::min(TheMin_sc, eVect.minCoeff());
      } else {
	TheMin_sc = RecS.minval;
      }




      std::vector<double> data_rangeA(2);
      std::vector<double> data_rangeB(2);
      data_rangeA[0] = TheMin_sc;
      data_rangeA[1] = TheMax_sc;
      data_rangeB[0] = TheMin_sc;
      data_rangeB[1] = TheMax_sc;
      //
      std::vector<MyVector<double>> ListVect = ListListVect[iBuoy];
      std::vector<std::vector<double>> ListVectSel(nbGrid+1);
      for (int iTime=0; iTime<nbTime; iTime++) {
	bool IsCorrect=true;
	for (int iGrid=0; iGrid<=nbGrid; iGrid++) {
	  double eVal = ListVect[iGrid](iTime);
	  if (eVal < TheMin_sc)
	    IsCorrect=false;
	  if (eVal > TheMax_sc)
	    IsCorrect=false;
	}
	if (IsCorrect) {
	  for (int iGrid=0; iGrid<=nbGrid; iGrid++) {
	    double eVal = ListVect[iGrid](iTime);
	    ListVectSel[iGrid].push_back(eVal);
	  }
	}
      }
      for (int iGrid=0; iGrid<nbGrid; iGrid++) {
	MyVector<double> eVectA = VectorFromStdVector(ListVectSel[0]);
	MyVector<double> eVectB = VectorFromStdVector(ListVectSel[iGrid+1]);
	int len=eVectA.size();
	for (int u=0; u<len; u++) {
	  std::cerr << "u=" << u << " / " << len << " eVectA=" << eVectA(u) << " eVectB=" << eVectB(u) << "\n";
	}
	T_stat eStat = ComputeStatistics_MyVector(eVectA, eVectB);
	std::string eName = "iBuoy" + std::to_string(iBuoy) + "_iGrid" + std::to_string(iGrid) + "_Global";
	Print_Down_Statistics(std::cerr, eName, eStat);
	os << "eVarName=" << eVarName << " iBuoy=" << iBuoy << " iGrid=" << iGrid << "\n";
	Print_Down_Statistics(os, eName, eStat);
	if (DoScatterPlot) {
	  DrawScatterArr eDrw;
	  eDrw.VarNameAB_file="Scatter_" + std::to_string(iBuoy+1) + "_" + eVarName;
	  eDrw.DoTitle=true;
	  eDrw.AddStatMeasModel=true;
	  eDrw.NameA_plot="Data (" + eUnit + ")";
	  eDrw.NameB_plot="Model (" + eUnit + ")";
	  eDrw.data_rangeA=data_rangeA;
	  eDrw.data_rangeB=data_rangeB;
	  eDrw.eVectA=eVectA;
	  eDrw.eVectB=eVectB;
	  eDrw.aSize=100;
	  eDrw.bSize=100;
	  PLOT_SCATTER(eDrw, eCall, ePerm);
	}
      }
    }
    std::cerr << "After plotting of the model data + Buoy\n";
  }
}


template<typename T>
std::vector<T> ReplicateInformation(std::vector<T> const& V, size_t const& n_elt, std::string const& entry)
{
  if (n_elt == V.size()) {
    return V;
  }
  if (V.size() == 1) {
    return std::vector<T>(n_elt, V[0]);
  }
  std::cerr << "ReplicateInformation operation failed entry=" << entry << "\n";
  std::cerr << "The number of element n_elt=" << n_elt << "\n";
  std::cerr << "|V|=" << V.size() << "\n";
  std::cerr << "The allowed sizes are 1 or n_elt\n";
  throw TerminalException{1};
}


void PointOutputPlot(FullNamelist const& eFull)
{
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  //
  // Reading grid arrays and the like
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::vector<std::string> ListModelName=eBlPROC.ListListStringValues.at("ListMODELNAME");
  std::vector<std::string> ListGridFile=eBlPROC.ListListStringValues.at("ListGridFile");
  std::vector<std::string> ListHisPrefix=eBlPROC.ListListStringValues.at("ListHisPrefix");
  std::vector<std::string> ListRunName=eBlPROC.ListListStringValues.at("ListRunName");
  std::vector<std::string> ListVarName=eBlPROC.ListListStringValues.at("ListVarName");
  bool DoLinePlot = eBlPROC.ListBoolValues.at("DoLinePlot");
  bool DoCsvFile = eBlPROC.ListBoolValues.at("DoCsvFile");
  bool DoSeasonMonthlyAverages = eBlPROC.ListBoolValues.at("DoSeasonMonthlyAverages");
  std::string PicPrefix = eBlPROC.ListStringValues.at("PicPrefix");
  std::vector<size_t> LSiz = {ListModelName.size(), ListGridFile.size(), ListHisPrefix.size(), ListRunName.size(), ListVarName.size()};
  size_t nbGridVar_t = *std::max_element(LSiz.begin(), LSiz.end());
  int nbGridVar = nbGridVar_t;
  //
  size_t nbGrid_t = ListModelName.size();
  int nbGrid = nbGrid_t;
  if (nbGrid_t != ListGridFile.size() || nbGrid_t != ListHisPrefix.size()) {
    std::cerr << "inconsistent input\n";
    std::cerr << "|ListModelName| = " << ListModelName.size() << "\n";
    std::cerr << "|ListGridFile|  = " << ListGridFile.size()  << "\n";
    std::cerr << "|ListHisPrefix| = " << ListHisPrefix.size() << "\n";
    throw TerminalException{1};
  }
  //
  std::vector<GridArray> ListGrdArr;
  std::vector<TotalArrGetData> ListTotalArr;
  for (int iGrid=0; iGrid<nbGrid; iGrid++) {
    std::cerr << "PointOutputPlot : iGrid=" << iGrid << " / " << nbGrid << "\n";
    std::string eModelName=ListModelName[iGrid];
    std::string GridFile=ListGridFile[iGrid];
    std::string HisPrefix=ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    ListGrdArr.push_back(GrdArr);
    TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
    ListTotalArr.push_back(TotalArr);
  }
  std::cerr << "Before replication step\n";
  ListGrdArr = ReplicateInformation(ListGrdArr, nbGridVar_t, "ListGrdArr");
  ListTotalArr = ReplicateInformation(ListTotalArr, nbGridVar_t, "ListTotalArr");
  ListRunName   = ReplicateInformation(ListRunName  , nbGridVar_t, "ListRunName");
  ListVarName   = ReplicateInformation(ListVarName  , nbGridVar_t, "ListVarName");
  std::cerr << "After replication step\n";
  //
  // Reading time information
  //
  std::string strBEGTC=eBlPROC.ListStringValues.at("BEGTC");
  std::string strENDTC=eBlPROC.ListStringValues.at("ENDTC");
  double BeginTime=CT2MJD(strBEGTC);
  double EndTime  =CT2MJD(strENDTC);
  double DELTC=eBlPROC.ListDoubleValues.at("DELTC");
  std::string UNITC=eBlPROC.ListStringValues.at("UNITC");
  double DeltaInterval=GetIntervalSize(DELTC, UNITC);
  std::vector<double> ListTime = GetIntervalFLD(BeginTime, EndTime, DeltaInterval);
  int nbTime = ListTime.size();
  //
  // Reading buoy location and computing interpolation arrays
  //
  std::vector<double> ListPointLon=eBlPROC.ListListDoubleValues.at("ListPointLongitude");
  std::vector<double> ListPointLat=eBlPROC.ListListDoubleValues.at("ListPointLatitude");
  std::vector<std::string> ListPointName=eBlPROC.ListListStringValues.at("ListPointName");
  int nbBuoy = ListPointLon.size();
  MyMatrix<double> ListXY(2,nbBuoy);
  for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
    ListXY(0,iBuoy) = ListPointLon[iBuoy];
    ListXY(1,iBuoy) = ListPointLat[iBuoy];
  }
  CheckingGridPointStatus(ListGrdArr, ListXY);
  std::vector<SingleArrayInterpolation> ListRec(nbGridVar);
  for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
    ListRec[iGridVar] = ComputeArrayInterpolation_ListXY(ListGrdArr[iGridVar], ListXY);
  //
  // Reading the regions of the averaging.
  //
  std::vector<std::string> ListAverageRegionLon=eBlPROC.ListListStringValues.at("ListAverageRegionLongitude");
  std::vector<std::string> ListAverageRegionLat=eBlPROC.ListListStringValues.at("ListAverageRegionLatitude");
  std::vector<std::string> ListAverageRegionName=eBlPROC.ListListStringValues.at("ListAverageRegionName");
  int nbRegion = ListAverageRegionLon.size();
  std::vector<std::pair<std::vector<double>, std::vector<double>>> ListRegions;
  for (int iRegion=0; iRegion<nbRegion; iRegion++) {
    std::vector<double> LonPoly, LatPoly;
    std::vector<std::string> LStr_lon = STRING_Split(ListAverageRegionLon[iRegion], " ");
    std::vector<std::string> LStr_lat = STRING_Split(ListAverageRegionLat[iRegion], " ");
    size_t len = LStr_lon.size();
    for (size_t i=0; i<len; i++) {
      LonPoly.push_back(ParseScalar<double>(LStr_lon[i]));
      LatPoly.push_back(ParseScalar<double>(LStr_lat[i]));
    }
    ListRegions.push_back({LonPoly, LatPoly});
  }
  std::vector<SingleArrayRegionAveraging> ListRecRegAve(nbGridVar);
  for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
    ListRecRegAve[iGridVar] = ComputeArrayRegionAveraging_ListPolygon(ListGrdArr[iGridVar], ListRegions);
  //
  // Reading full data set
  //
  bool VariableMin=eBlPLOT.ListBoolValues.at("VariableMin");
  bool VariableMax=eBlPLOT.ListBoolValues.at("VariableMax");
  double SpecifiedMin = eBlPLOT.ListDoubleValues.at("SpecifiedMin");
  double SpecifiedMax = eBlPLOT.ListDoubleValues.at("SpecifiedMax");
  bool DoExplicitLabel=eBlPLOT.ListBoolValues.at("DoExplicitLabel");
  bool DrawHorizVertLines = eBlPLOT.ListBoolValues.at("DrawHorizVertLines");
  int nbLabel=eBlPLOT.ListIntValues.at("nbLabel");
  std::string StyleDate = eBlPLOT.ListStringValues.at("StyleDate");
  PlotBound ePlotBound=ReadPlotBound(eFull);
  int nbBlock = eBlPROC.ListIntValues.at("nbBlock");
  //
  // Loading data for the plots.
  // This structure is needed for efficient loading.
  //
  std::vector<int> ListPos = DivideListPosition(nbTime, nbBlock);
  std::vector<std::vector<MyVector<double>>> ListListVect(nbBuoy+nbRegion, std::vector<MyVector<double>>(nbGridVar, MyVector<double>(nbTime)));
  RecSymbolic RecS;
  for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
    std::string eVarName = ListVarName[iGridVar];
    for (int iTime=0; iTime<nbTime; iTime++) {
      double eTimeDay = ListTime[iTime];
      VarQuery eQuery;
      eQuery.eTimeDay = eTimeDay;
      eQuery.iTime = -1;
      eQuery.NatureQuery = "instant";
      eQuery.TimeFrameDay = -1;
      eQuery.typeQuery = "unset";
      std::cerr << "iTime=" << iTime << " / " << nbTime << "\n";
      RecVar eRecVar = ModelSpecificVarSpecificTimeGeneral(ListTotalArr[iGridVar], eVarName, eQuery, ePlotBound);
      RecS = eRecVar.RecS;
      std::string eUnit = RecS.Unit;
      RecVar fRecVar = INTERPOL_SingleRecVarInterpolation(ListRec[iGridVar], eRecVar);
      for (int iBuoy=0; iBuoy<nbBuoy; iBuoy++) {
        double eVal = fRecVar.F(iBuoy,0);
        ListListVect[iBuoy][iGridVar](iTime) = eVal;
      }
      RecVar gRecVar = REGAVE_SingleRecVarAveraging(ListRecRegAve[iGridVar], eRecVar);
      for (int iRegion=0; iRegion<nbRegion; iRegion++) {
        double eVal = gRecVar.F(iRegion,0);
        ListListVect[nbBuoy+iRegion][iGridVar](iTime) = eVal;
      }
    }
  }
  //
  // Creating CSV files
  //
  if (DoCsvFile) {
    for (int iBuoyReg=0; iBuoyReg<nbBuoy+nbRegion; iBuoyReg++) {
      std::cerr << "iBuoyReg=" << iBuoyReg << " nbBuoy=" << nbBuoy << " nbRegion=" << nbRegion << "\n";
      std::string OutputFile = PicPrefix + "/csv_file" + std::to_string(iBuoyReg+1) + ".csv";
      std::ofstream os(OutputFile);
      //
      os << "Short Name";
      for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
        os << ";" << ListVarName[iGridVar];
      os << "\n";
      //
      std::vector<RecSymbolic> ListRecSymb(nbGridVar);
      for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
        ListRecSymb[iGridVar] = RetrieveTrivialRecVar(ListVarName[iGridVar]).RecS;
      //
      os << "Long name";
      for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
        os << ";" << ListRecSymb[iGridVar].VarName2;
      os << "\n";
      //
      os << "unit";
      for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
        os << ";" << ListRecSymb[iGridVar].Unit;
      os << "\n";
      //
      for (int iTime=0; iTime<nbTime; iTime++) {
        double eTime = ListTime[iTime];
        std::string dateTimeStr = DATE_ConvertMjd2mystringPres(eTime);
        os << dateTimeStr;
        for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
          os << ";" << ListListVect[iBuoyReg][iGridVar](iTime);
        }
        os << "\n";
      }
    }
  }
  //
  // Text Output
  //
  if (DoSeasonMonthlyAverages) {
    for (int iBuoyReg=0; iBuoyReg<nbBuoy+nbRegion; iBuoyReg++) {
      MyMatrix<int>    SumMonth_I  = ZeroMatrix<int>(nbGridVar,12);
      MyMatrix<int>    SumSeason_I = ZeroMatrix<int>(nbGridVar,4);
      MyMatrix<double> SumMonth_D  = ZeroMatrix<double>(nbGridVar,12);
      MyMatrix<double> SumSeason_D = ZeroMatrix<double>(nbGridVar,4);
      std::map<std::pair<int,int>,std::vector<double>> MapMonth_D;
      std::map<std::pair<int,int>,std::vector<int>> MapMonth_I;
      std::map<std::pair<int,int>,std::vector<double>> MapSeason_D;
      std::map<std::pair<int,int>,std::vector<int>> MapSeason_I;
      //
      std::string OutputFile = PicPrefix + "/interpolated_results" + std::to_string(iBuoyReg+1) + ".txt";
      std::ofstream os(OutputFile);
      os << "nbTime=" << nbTime << "\n";
      os << "VARS";
      for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
        os << "," << ListVarName[iGridVar];
      }
      os << "\n";
      for (int iTime=0; iTime<nbTime; iTime++) {
        double eTime = ListTime[iTime];
        int eYear = DATE_GetYear(eTime);
        int eMonth = DATE_GetMonth(eTime);
        int eSeason = DATE_GetSeason(eTime);
        std::pair<int,int> eYearMonth{eYear,eMonth};
        std::pair<int,int> eYearSeason{eYear,eSeason};
        std::string eTimeStr = DATE_ConvertMjd2mystringPres(eTime);
        os << eTimeStr;
        if (MapMonth_D.find(eYearMonth) == MapMonth_D.end()) {
          MapMonth_D[eYearMonth] = std::vector<double>(nbGridVar,0);
          MapMonth_I[eYearMonth] = std::vector<int>(nbGridVar,0);
        }
        if (MapSeason_D.find(eYearSeason) == MapSeason_D.end()) {
          MapSeason_D[eYearSeason] = std::vector<double>(nbGridVar,0);
          MapSeason_I[eYearSeason] = std::vector<int>(nbGridVar,0);
        }
        for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
          double eVal = ListListVect[iBuoyReg][iGridVar](iTime);
          os << " " << eVal;
          SumMonth_I(iGridVar, eMonth)++;
          SumMonth_D(iGridVar, eMonth) += eVal;
          SumSeason_I(iGridVar, eSeason)++;
          SumSeason_D(iGridVar, eSeason) += eVal;
          //
          MapMonth_D[eYearMonth][iGridVar] += eVal;
          MapMonth_I[eYearMonth][iGridVar] ++;
          MapSeason_D[eYearSeason][iGridVar] += eVal;
          MapSeason_I[eYearSeason][iGridVar] ++;
        }
        os << "\n";
      }
      //
      {
        std::string OutputFile = PicPrefix + "/interpolated_results" + std::to_string(iBuoyReg+1) + "_month_season.txt";
        std::ofstream os(OutputFile);
        os << "VARS";
        for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
          os << "," << ListVarName[iGridVar];
        }
        os << "\n";
        for (int iMonth=0; iMonth<12; iMonth++) {
          os << GetMonthName(iMonth+1);
          for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
            double avgVal = SumMonth_D(iGridVar, iMonth) / double(SumMonth_I(iGridVar, iMonth));
            os << "," << avgVal;
          }
          os << "\n";
        }
        for (int iSeason=0; iSeason<4; iSeason++) {
          os << GetSeasonName(iSeason+1);
          for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
            std::cerr << "sumseason_D=" << SumSeason_D(iGridVar, iSeason) << "\n";
            std::cerr << "sumseason_i=" << SumSeason_I(iGridVar, iSeason) << "\n";
            double avgVal = SumSeason_D(iGridVar, iSeason) / double(SumSeason_I(iGridVar, iSeason));
            os << "," << avgVal;
          }
          os << "\n";
        }
        os << "-----------------------------------------\n";
        for (auto & ePair : MapMonth_D) {
          std::vector<double> LMonth_D = MapMonth_D[ePair.first];
          std::vector<int> LMonth_I = MapMonth_I[ePair.first];
          os << ePair.first.first << " " << GetMonthName(ePair.first.second+1);
          for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
            double avgVal = LMonth_D[iGridVar] / double(LMonth_I[iGridVar]);
            os << "," << avgVal;
          }
          os << "\n";
        }
        for (auto & ePair : MapSeason_D) {
          std::vector<double> LSeason_D = MapSeason_D[ePair.first];
          std::vector<int> LSeason_I = MapSeason_I[ePair.first];
          os << ePair.first.first << " " << GetSeasonName(ePair.first.second+1);
          for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
            double avgVal = LSeason_D[iGridVar] / double(LSeason_I[iGridVar]);
            os << "," << avgVal;
          }
          os << "\n";
        }
      }
    }
  }
  //
  // line plots
  //
  if (DoLinePlot) {
    for (int iBuoyReg=0; iBuoyReg<nbBuoy+nbRegion; iBuoyReg++) {
      std::cerr << "After reading of the model data\n";
      for (int iBlock=0; iBlock<nbBlock; iBlock++) {
        DrawLinesArr eDrawArr;
        eDrawArr.DoTitle=true;
        eDrawArr.TitleStr="Time series for location " + ListPointName[iBuoyReg];
        eDrawArr.IsTimeSeries=true;
        eDrawArr.PairComparison=false;
        eDrawArr.nbLabel=nbLabel;
        eDrawArr.DoExplicitLabel=DoExplicitLabel;
        eDrawArr.DrawHorizVertLines=DrawHorizVertLines;
        eDrawArr.StyleDate=StyleDate;
        eDrawArr.VarName=std::to_string(iBuoyReg+1) + "_" + std::to_string(iBlock);
        eDrawArr.ListName_plot=ListRunName;
        eDrawArr.YAxisString="";
        //
        // Determination of dimension variable
        //
        int pos1 = ListPos[iBlock];
        int pos2 = ListPos[iBlock+1];
        int len = pos2 - pos1;
        MyVector<double> ListTime_SEL(len);
        for (int pos=pos1; pos<pos2; pos++)
          ListTime_SEL(pos - pos1) = ListTime[pos];
        std::vector<MyVector<double>> ListVect_SEL(nbGridVar);
        for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++) {
          MyVector<double> eVect(len);
          for (int pos=pos1; pos<pos2; pos++) {
            double eVal = ListListVect[iBuoyReg][iGridVar](pos);
            eVect(pos - pos1) = eVal;
          }
          ListVect_SEL[iGridVar] = eVect;
          double maxVal = eVect.maxCoeff();
          std::cerr << "iGridVar=" << iGridVar << " maxVal=" << maxVal << "\n";
        }
        std::cerr << "iBuoyReg=" << iBuoyReg << " iBlock=" << iBlock << "\n";
        for (int u=0; u<len; u++) {
          std::cerr << "u=" << u << " / " << len << "   V =";
          for (int iGridVar=0; iGridVar<nbGridVar; iGridVar++)
            std::cerr << " " << ListVect_SEL[iGridVar](u);
          std::cerr << "\n";
        }
        eDrawArr.ListX = ListTime_SEL;
        eDrawArr.ListListVect = ListVect_SEL;
        //
        // Determinantion of min/max ranges
        //
        double TheMax, TheMin;
        if (VariableMax) {
          TheMax = - 1000000;
          for (auto & eVect : eDrawArr.ListListVect)
            TheMax = std::max(TheMax, eVect.maxCoeff());
        } else {
          TheMax = SpecifiedMax;
        }
        if (VariableMin) {
          TheMin = 1000000;
          for (auto & eVect : eDrawArr.ListListVect)
            TheMin = std::min(TheMin, eVect.minCoeff());
        } else {
          TheMin = SpecifiedMin;
        }
        std::cerr << "TheMin=" << TheMin << " TheMax=" << TheMax << "\n";
        eDrawArr.TheMax=TheMax;
        eDrawArr.TheMin=TheMin;
        eDrawArr.YAxisString="(" + RecS.Unit + ")";
        //
        std::string FileName=ePerm.eDir + "TimeSeries_iBuoyReg" + std::to_string(iBuoyReg+1) + "_iBlock" + std::to_string(iBlock);
        LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
      }
    }
  }
  std::cerr << "After plotting of the interpolated model data\n";
}



#endif
