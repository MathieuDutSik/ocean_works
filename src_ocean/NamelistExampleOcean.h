#ifndef INCLUDE_NAMELIST_EXAMPLE_OCEAN
#define INCLUDE_NAMELIST_EXAMPLE_OCEAN

#include "Namelist.h"
#include "Model_data_loading.h"
#include "Model_grids.h"






FullNamelist NAMELIST_ComparisonSequentialRuns()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  ListIntValues1["GEOSELECTION"]=1;
  ListDoubleValues1["MinLON"]=-7;
  ListDoubleValues1["MaxLON"]=37;
  ListDoubleValues1["MinLAT"]=30;
  ListDoubleValues1["MaxLAT"]=46;
  ListListDoubleValues1["LONPOLY"]={10, 10, 10};
  ListListDoubleValues1["LATPOLY"]={10, 10, 10};
  ListStringValues1["HisPrefix"]="unset";
  ListStringValues1["ModelName"]="unset";
  ListStringValues1["shortName"]="unset";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  ListBlock["PROC"]=BlockPROC;
  // Final part
  return {std::move(ListBlock), "undefined"};
}


FullNamelist NAMELIST_GetStandard_CREATE_TracerSourceTerm()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["MODELNAME"]="unset";
  ListStringValues1["GridFile"]="unset";
  ListStringValues1["OutFile"]="Float_Output_";
  ListListDoubleValues1["ListLON"]={};
  ListListDoubleValues1["ListLAT"]={};
  ListListDoubleValues1["ListDistKM"]={};
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  ListBlock["PROC"]=BlockPROC;
  // Final part
  return {std::move(ListBlock), "undefined"};
}






FullNamelist NAMELIST_GetStandard_ComputeFloatTrajectories()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListListStringValues1["ListMODELNAME"]={"unset MODELNAME in ListMODELNAME"};
  ListListStringValues1["ListGridFile"]={"unset GridFile in ListGridFile"};
  ListListStringValues1["ListHisPrefix"]={"ROMS_output_"};
  ListListIntValues1["ListFatherGrid"]={-1};
  ListStringValues1["OutFile"]="Float_Output_";
  ListDoubleValues1["DEFINETC"]=86400;
  ListDoubleValues1["HISTIME"]=1800;
  ListIntValues1["NPROC"]=1;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  ListBlock["PROC"]=BlockPROC;
  // FLOAT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["DoTitle"]=true;
  ListBoolValues2["VariableRange"]=false;
  SingleBlock BlockFLOAT;
  BlockFLOAT.ListIntValues=ListIntValues2;
  BlockFLOAT.ListBoolValues=ListBoolValues2;
  BlockFLOAT.ListDoubleValues=ListDoubleValues2;
  BlockFLOAT.ListStringValues=ListStringValues2;
  BlockFLOAT.ListListStringValues=ListListStringValues2;
  BlockFLOAT.ListListDoubleValues=ListListDoubleValues2;
  ListBlock["FLOAT"]=BlockFLOAT;
  // Final part
  return {std::move(ListBlock), "undefined"};
}



FullNamelist NAMELIST_GetStandard_PlotRomsFloats()
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
  ListStringValues1["MODELNAME"]="unset MODELNAME in ListMODELNAME";
  ListStringValues1["GridFile"]="unset GridFile in ListGridFile";
  ListStringValues1["HisPrefix"]="ROMS_output_";
  ListStringValues1["BoundFile"]="bound.nc";
  ListStringValues1["FloatFile"]="irrelevant.nc";
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListStringValues1["__NaturePlot"]="ROMS_DRIFTER";
  ListStringValues1["FileDescFloat"]="unset";
  ListStringValues1["FileListBlocks"]="unset";
  ListStringValues1["FileListBlockNames"]="unset";
  ListStringValues1["FileDrifterStartEnd"]="unset";
  ListListStringValues1["ListNatureQuery"]={"instant"}; // By default instantaneous values
  ListStringValues1["Sphericity"]="unset";
  ListBoolValues1["CutWorldMap"]=false;
  ListBoolValues1["HigherLatitudeCut"]=false;
  ListBoolValues1["SplittingAt180"]=false;
  ListDoubleValues1["MinLatCut"]=-80;
  ListDoubleValues1["MaxLatCut"]=80;
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
  ListBoolValues2["DoTitle"]=true;
  ListBoolValues2["DoMain"]=true;
  ListListDoubleValues2["SpatialResolutionTransectKM"]={1};
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  ListStringValues2["TypeListPoint"]="empty"; // possible values empty, fileRovinj, namelist
  ListStringValues2["ListPointFile"]="";
  ListListDoubleValues2["ListFrameMinLon"]={};
  ListListDoubleValues2["ListFrameMinLat"]={};
  ListListDoubleValues2["ListFrameMaxLon"]={};
  ListListDoubleValues2["ListFrameMaxLat"]={};
  ListListDoubleValues2["ListPointLon"]={};
  ListListDoubleValues2["ListPointLat"]={};
  ListListDoubleValues2["ListPointDepth"]={};
  ListListStringValues2["ListPointName"]={};
  ListDoubleValues2["VertResolM"]=0.2;
  ListBoolValues2["PlotDensity"]=false;
  ListBoolValues2["PlotTrajectory"]=false;
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["PrintMMA"]=false;
  ListIntValues2["nbLevelSpa"]=50;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListDoubleValues2["AnnotationLon"]=0;
  ListDoubleValues2["AnnotationLat"]=0;
  ListStringValues2["AnnotationText"]="something to write";
  ListBoolValues2["FillLand"]=true;
  ListStringValues2["GridResolution"]="HighRes";
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["DoTitle"]=true;
  ListDoubleValues2["vcRefLengthF"]=0.02;
  ListBoolValues2["DoColorBar"]=true;
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnFillOn"]=true;
  ListBoolValues2["cnLinesOn"]=false;
  ListBoolValues2["cnLineLabelsOn"]=false;
  ListBoolValues2["cnSmoothingOn"]=true;
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DoTitleString"]=true;
  ListStringValues2["LandPortr"]="Landscape";
  ListStringValues2["optStatStr"]="double";
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // Final part
  return {std::move(ListBlock), "undefined"};
}






FullNamelist NAMELIST_GetStandard_PlotTransect()
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
  ListStringValues1["KindSelect"]="direct"; // possible values: direct, monthly, seasonal, specific
  ListListStringValues1["ListSpecificTimes"]={};
  ListDoubleValues1["TimeFrameDay"]=1;
  ListListStringValues1["ListMODELNAME"]={"unset MODELNAME in ListMODELNAME"};
  ListListStringValues1["ListGridFile"]={"unset GridFile in ListGridFile"};
  ListListStringValues1["ListHisPrefix"]={"ROMS_output_"};
  ListListStringValues1["ListRunName"]={"Explicit scheme"};
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListStringValues1["__NaturePlot"]="TRANSECT";
  ListListStringValues1["ListNatureQuery"]={"instant"}; // By default instantaneous values
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
  ListBoolValues2["DoTitle"]=true;
  ListListDoubleValues2["TransectLonStart"]={0};
  ListListDoubleValues2["TransectLatStart"]={0};
  ListListDoubleValues2["TransectLonEnd"]={0};
  ListListDoubleValues2["TransectLatEnd"]={0};
  ListListDoubleValues2["SpatialResolutionTransectKM"]={1};
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  ListStringValues2["TypeListPoint"]="empty"; // possible values empty, fileRovinj, namelist
  ListStringValues2["ListPointFile"]="";
  ListListDoubleValues2["ListPointLon"]={};
  ListListDoubleValues2["ListPointLat"]={};
  ListListDoubleValues2["ListPointDepth"]={};
  ListListStringValues2["ListPointName"]={};
  ListDoubleValues2["VertResolM"]=0.2;
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["PrintMMA"]=false;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // VARS
  std::map<std::string, bool> ListBoolValues3;
  std::vector<std::string> ListVarOut=GetAllPossibleVariables();
  for (auto& eVal : ListVarOut)
    ListBoolValues3[eVal]=false;
  SingleBlock BlockVARS;
  BlockVARS.ListBoolValues=ListBoolValues3;
  ListBlock["VARS"]=BlockVARS;
  // Final part
  return {std::move(ListBlock), "undefined"};
}





FullNamelist NAMELIST_GetStandard_PlotGrid()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::string LPoss="Possibilities:";
  bool IsFirst=true;
  for (auto & eStr : GetAllPossibleModels()) {
    if (!IsFirst)
      LPoss += ",";
    IsFirst=false;
    LPoss += " " + eStr;
  }
  ListStringValues1["MODELNAME"]=LPoss;
  ListStringValues1["GridFile"]="unset GridFile";
  ListStringValues1["BoundFile"]="unset";
  ListStringValues1["HisPrefix"]="irrelevant_should_stay_that_way";
  ListStringValues1["Sphericity"]="unset";
  ListBoolValues1["CutWorldMap"]=false;
  ListBoolValues1["HigherLatitudeCut"]=false;
  ListBoolValues1["SplittingAt180"]=false;
  ListDoubleValues1["MinLatCut"]=-80;
  ListDoubleValues1["MaxLatCut"]=80;
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListBoolValues1["FirstCleanDirectory"]=true;
  ListBoolValues1["KeepNC_NCL"]=false;
  ListBoolValues1["InPlaceRun"]=false;
  ListBoolValues1["PrintDebugInfo"]=false;
  ListBoolValues1["OnlyCreateFiles"]=false;
  ListBoolValues1["OverwritePrevious"]=false;
  ListBoolValues1["WriteITimeInFileName"]=true;
  ListIntValues1["NPROC"]=1;
  ListIntValues1["MaxNbTime"]=-1;
  ListDoubleValues1["ThresholdApplyAverage"]=10000;
  ListBoolValues1["ApplyThresholdAveraging"]=false;
  ListStringValues1["Pcolor_method"]="ncl";
  ListStringValues1["Quiver_method"]="ncl";
  ListStringValues1["Lines_method"]="ncl";
  ListStringValues1["Scatter_method"]="ncl";
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
  std::map<std::string, std::vector<int>> ListListIntValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListStringValues2["ColorMapDiff"]="BlWhRe";
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnFillOn"]=true;
  ListBoolValues2["cnLinesOn"]=false;
  ListBoolValues2["cnLineLabelsOn"]=false;
  ListBoolValues2["cnSmoothingOn"]=true;
  ListBoolValues2["DoColorBar"]=true;
  ListIntValues2["nbLevelSpa"]=50;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["DoTitle"]=true;
  ListStringValues2["GridResolution"]="HighRes";
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["PrintMMA"]=false;
  ListBoolValues2["LocateMM"]=false;
  ListBoolValues2["DoMain"]=true;
  ListBoolValues2["PlotDepth"]=true;
  ListBoolValues2["PlotMesh"]=false;
  ListBoolValues2["PlotIOBP"]=false;
  ListBoolValues2["DoTitleString"]=true;
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListBoolValues2["ExcludeLargeValues"]=false;
  ListDoubleValues2["ThresholdExclusionPlot"]=100000;
  ListDoubleValues2["MultiplierResolutionRegrid"]=1;
  ListListIntValues2["Tens3ListLevel"]={};
  ListDoubleValues2["AnnotationLon"]=0;
  ListDoubleValues2["AnnotationLat"]=0;
  ListStringValues2["AnnotationText"]="something to write";
  ListDoubleValues2["vcRefLengthF"]=0.02;
  ListStringValues2["LandPortr"]="Landscape";
  ListStringValues2["optStatStr"]="double";
  ListListStringValues2["RenameVariable_VarName1"]={};
  ListListStringValues2["RenameVariable_VarName2"]={};
  ListListStringValues2["RenameVariable_VarNameUF"]={};
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  ListListStringValues2["BoundDiff_var"]={};
  ListListDoubleValues2["BoundDiff_min"]={};
  ListListDoubleValues2["BoundDiff_max"]={};
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["FillLand"]=true;
  ListListDoubleValues2["ListFrameMinLon"]={};
  ListListDoubleValues2["ListFrameMinLat"]={};
  ListListDoubleValues2["ListFrameMaxLon"]={};
  ListListDoubleValues2["ListFrameMaxLat"]={};
  ListListDoubleValues2["TransectLonStart"]={};
  ListListDoubleValues2["TransectLatStart"]={};
  ListListDoubleValues2["TransectLonEnd"]={};
  ListListDoubleValues2["TransectLatEnd"]={};
  ListListDoubleValues2["TransectSpatialResolutionKM"]={};
  ListListDoubleValues2["TransectVerticalResolutionM"]={};
  ListDoubleValues2["FrameLonLat"]=0.5;
  ListBoolValues2["UseRegridArray"]=false;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  BlockPLOT.ListListIntValues=ListListIntValues2;
  ListBlock["PLOT"]=BlockPLOT;
  //
  return {std::move(ListBlock), "undefined"};
}





FullNamelist NAMELIST_GetStandard_PlotRoutine_common()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::string LPoss="Possibilities:";
  bool IsFirst=true;
  for (auto & eStr : GetAllPossibleModels()) {
    if (!IsFirst)
      LPoss += ",";
    IsFirst=false;
    LPoss += " " + eStr;
  }
  ListStringValues1["MODELNAME"]=LPoss;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListStringValues1["KindSelect"]="direct"; // possible values: direct, monthly, seasonal, specific
  ListListStringValues1["ListSpecificTimes"]={};
  ListStringValues1["GridFile"]="unset GridFile";
  ListStringValues1["BoundFile"]="unset";
  ListStringValues1["Sphericity"]="unset";
  ListBoolValues1["CutWorldMap"]=false;
  ListBoolValues1["HigherLatitudeCut"]=false;
  ListBoolValues1["SplittingAt180"]=false;
  ListDoubleValues1["MinLatCut"]=-80;
  ListDoubleValues1["MaxLatCut"]=80;
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListListStringValues1["ListNatureQuery"]={"instant"}; // By default instantaneous values
  ListDoubleValues1["TimeFrameDay"]=1;
  ListBoolValues1["FirstCleanDirectory"]=true;
  ListBoolValues1["KeepNC_NCL"]=false;
  ListBoolValues1["InPlaceRun"]=false;
  ListBoolValues1["PrintDebugInfo"]=false;
  ListBoolValues1["OnlyCreateFiles"]=false;
  ListBoolValues1["OverwritePrevious"]=false;
  ListBoolValues1["WriteITimeInFileName"]=true;
  ListIntValues1["NPROC"]=1;
  ListIntValues1["MaxNbTime"]=-1;
  ListDoubleValues1["ThresholdApplyAverage"]=10000;
  ListBoolValues1["ApplyThresholdAveraging"]=false;
  ListStringValues1["Pcolor_method"]="ncl";
  ListStringValues1["Quiver_method"]="ncl";
  ListStringValues1["Lines_method"]="ncl";
  ListStringValues1["Scatter_method"]="ncl";
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
  std::map<std::string, std::vector<int>> ListListIntValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListStringValues2["ColorMapDiff"]="BlWhRe";
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnFillOn"]=true;
  ListBoolValues2["cnLinesOn"]=false;
  ListBoolValues2["cnLineLabelsOn"]=false;
  ListBoolValues2["DoColorBar"]=true;
  ListBoolValues2["cnSmoothingOn"]=true;
  ListIntValues2["nbLevelSpa"]=50;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["DoTitle"]=true;
  ListStringValues2["GridResolution"]="HighRes";
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["PrintMMA"]=false;
  ListBoolValues2["LocateMM"]=false;
  ListBoolValues2["DoMain"]=true;
  ListBoolValues2["PlotDepth"]=true;
  ListBoolValues2["PlotMesh"]=false;
  ListBoolValues2["PlotIOBP"]=false;
  ListBoolValues2["DoTitleString"]=true;
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListBoolValues2["ExcludeLargeValues"]=false;
  ListDoubleValues2["ThresholdExclusionPlot"]=100000;
  ListDoubleValues2["MultiplierResolutionRegrid"]=1;
  ListListIntValues2["Tens3ListLevel"]={};
  ListStringValues2["AnnotationText"]="something to write";
  ListDoubleValues2["vcRefLengthF"]=0.02;
  ListDoubleValues2["AnnotationLon"]=0;
  ListDoubleValues2["AnnotationLat"]=0;
  ListStringValues2["LandPortr"]="Landscape";
  ListStringValues2["optStatStr"]="double";
  ListListStringValues2["RenameVariable_VarName1"]={};
  ListListStringValues2["RenameVariable_VarName2"]={};
  ListListStringValues2["RenameVariable_VarNameUF"]={};
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  ListListStringValues2["BoundDiff_var"]={};
  ListListDoubleValues2["BoundDiff_min"]={};
  ListListDoubleValues2["BoundDiff_max"]={};
  ListBoolValues2["VariableRange"]=false;
  ListBoolValues2["FillLand"]=true;
  ListListDoubleValues2["ListFrameMinLon"]={};
  ListListDoubleValues2["ListFrameMinLat"]={};
  ListListDoubleValues2["ListFrameMaxLon"]={};
  ListListDoubleValues2["ListFrameMaxLat"]={};
  ListListDoubleValues2["TransectLonStart"]={};
  ListListDoubleValues2["TransectLatStart"]={};
  ListListDoubleValues2["TransectLonEnd"]={};
  ListListDoubleValues2["TransectLatEnd"]={};
  ListListDoubleValues2["TransectSpatialResolutionKM"]={};
  ListListDoubleValues2["TransectVerticalResolutionM"]={};
  ListDoubleValues2["FrameLonLat"]=0.5;
  ListBoolValues2["UseRegridArray"]=false;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  BlockPLOT.ListListIntValues=ListListIntValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // VARS
  std::map<std::string, bool> ListBoolValues3;
  std::vector<std::string> ListVarOut=GetAllPossibleVariables_with_pairs();
  for (auto& eVal : ListVarOut)
    ListBoolValues3[eVal]=false;
  SingleBlock BlockVARS;
  BlockVARS.ListBoolValues=ListBoolValues3;
  ListBlock["VARS"]=BlockVARS;
  // Final part
  //  std::cerr << "NAMELIST_GetStandard_PlotRoutine_common CutWorldMap=" << ListBlock.at("PROC").ListBoolValues.at("CutWorldMap") << "\n";
  return {std::move(ListBlock), "undefined"};
}





FullNamelist NAMELIST_GetStandard_PlotRoutine_single()
{
  FullNamelist eFull=NAMELIST_GetStandard_PlotRoutine_common();
  eFull.ListBlock["PROC"].ListStringValues["__NaturePlot"]="SINGLE";
  eFull.ListBlock["PROC"].ListStringValues["HisPrefix"]="relevant model entry";
  return eFull;
}

FullNamelist NAMELIST_GetStandard_PlotRoutine_pair()
{
  FullNamelist eFull=NAMELIST_GetStandard_PlotRoutine_common();
  eFull.ListBlock["PROC"].ListStringValues["__NaturePlot"]="PAIR";
  eFull.ListBlock["PROC"].ListStringValues["HisPrefix1"]="relevant model entry";
  eFull.ListBlock["PROC"].ListStringValues["HisPrefix2"]="relevant model entry";
  eFull.ListBlock["PROC"].ListStringValues["Name1"]="name to be written1";
  eFull.ListBlock["PROC"].ListStringValues["Name2"]="name to be written2";
  return eFull;
}




FullNamelist NAMELIST_GetStandardPLOT_DRIFTER_TRACK()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["DrifterFile"]="UNK";
  ListStringValues1["PicPrefix"]="UNK";
  ListStringValues1["Extension"]="png";
  ListStringValues1["__NaturePlot"]="drifters";
  ListBoolValues1["KeepNC_NCL"]=false;
  ListBoolValues1["InPlaceRun"]=false;
  ListBoolValues1["PrintDebugInfo"]=false;
  ListBoolValues1["OnlyCreateFiles"]=false;
  ListBoolValues1["FirstCleanDirectory"]=true;
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
  // PLOT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["IndividualTrajectories"]=false;
  ListBoolValues2["FixedFrame"]=false;
  ListDoubleValues2["MinimalDistanceTrajectoriesKM"]=0;
  ListDoubleValues2["MinimalTimeTrajectories"]=0;
  ListDoubleValues2["FrameLonLat"]=0;
  ListDoubleValues2["MinLon"]=12.2; // Those are for the Adriatic.
  ListDoubleValues2["MaxLon"]=19.5;
  ListDoubleValues2["MinLat"]=40.3;
  ListDoubleValues2["MaxLat"]=45.7;
  ListDoubleValues2["DeltaLonDensPlot"]=0.05;
  ListDoubleValues2["DeltaLatDensPlot"]=0.05;
  ListBoolValues2["DensityPassingPlot"]=false;
  ListBoolValues2["DensityFinalPlot"]=false;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  ListBlock["PLOT"]=BlockPLOT;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}








FullNamelist NAMELIST_GetStandard_CREATE_LTransInput()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  ListStringValues1["MODELNAME"]="unset MODELNAME";
  ListStringValues1["GridFile"]="unset GridFile";
  ListStringValues1["BoundFile"]="unset";
  ListStringValues1["HisPrefix"]="unset HisPrefix";
  ListDoubleValues1["MaxDistKM"]=10;
  ListListDoubleValues1["ListLonFloat"]={double(-40), double(-60)};
  ListListDoubleValues1["ListLatFloat"]={double(-40), double(-60)};
  ListListIntValues1["ListTime"]={0, 3600, 7200, 10800};
  ListListDoubleValues1["ListDepth"]={double(-2), double(-4), double(-6)};
  ListStringValues1["FloatFile"]="floats.in";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  ListBlock["PROC"]=BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}




FullNamelist NAMELIST_GetStandardWWM()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["PROCNAME"]="limon";
  ListIntValues1["DIMMODE"]=2;
  ListBoolValues1["LSTEA"]=false;
  ListBoolValues1["LQSTEA"]=false;
  ListBoolValues1["LSPHE"]=false;
  ListBoolValues1["LNAUTIN"]=true;
  ListBoolValues1["LMONO_IN"]=false;
  ListBoolValues1["LMONO_OUT"]=false;
  ListBoolValues1["LNAUTOUT"]=true;
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListDoubleValues1["DMIN"]=0.0001;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  ListBlock["PROC"]=BlockPROC;
  // COUPL
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["LCPL"]=false;
  ListBoolValues2["LROMS"]=false;
  ListBoolValues2["LTIMOR"]=false;
  ListBoolValues2["LSHYFEM"]=false;
  ListStringValues2["RADFLAG"]="LON";
  ListBoolValues2["LETOT"]=false;
  ListIntValues2["NLVT"]=10;
  ListDoubleValues2["DTCOUP"]=600;
  SingleBlock BlockCOUPL;
  BlockCOUPL.ListIntValues=ListIntValues2;
  BlockCOUPL.ListBoolValues=ListBoolValues2;
  BlockCOUPL.ListDoubleValues=ListDoubleValues2;
  BlockCOUPL.ListStringValues=ListStringValues2;
  BlockCOUPL.ListListStringValues=ListListStringValues2;
  ListBlock["COUPL"]=BlockCOUPL;
  // GRID
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListBoolValues3["LCIRD"]=true;
  ListBoolValues3["LSTAG"]=true;
  ListDoubleValues3["MINDIR"]=0;
  ListDoubleValues3["MAXDIR"]=360;
  ListIntValues3["MDC"]=36;
  ListIntValues3["MSC"]=36;
  ListDoubleValues3["FRLOW"]=0.04;
  ListDoubleValues3["FRHIGH"]=1;
  ListIntValues3["IGRIDTYPE"]=3;
  ListStringValues3["FILEGRID"]="hgrid.gr3";
  ListBoolValues3["LSLOP"]=false;
  ListDoubleValues3["SLMAX"]=0.2;
  ListBoolValues3["LVAR1D"]=false;
  SingleBlock BlockGRID;
  BlockGRID.ListIntValues=ListIntValues3;
  BlockGRID.ListBoolValues=ListBoolValues3;
  BlockGRID.ListDoubleValues=ListDoubleValues3;
  BlockGRID.ListStringValues=ListStringValues3;
  BlockGRID.ListListStringValues=ListListStringValues3;
  ListBlock["GRID"]=BlockGRID;
  // INIT
  std::map<std::string, int> ListIntValues4;
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, double> ListDoubleValues4;
  std::map<std::string, std::string> ListStringValues4;
  std::map<std::string, std::vector<std::string>> ListListStringValues4;
  ListBoolValues4["LHOTR"]=false;
  ListBoolValues4["LINID"]=true;
  ListIntValues4["INITSTYLE"]=1;
  SingleBlock BlockINIT;
  BlockINIT.ListIntValues=ListIntValues4;
  BlockINIT.ListBoolValues=ListBoolValues4;
  BlockINIT.ListDoubleValues=ListDoubleValues4;
  BlockINIT.ListStringValues=ListStringValues4;
  BlockINIT.ListListStringValues=ListListStringValues4;
  ListBlock["INIT"]=BlockINIT;
  // BOUC
  std::map<std::string, int> ListIntValues5;
  std::map<std::string, bool> ListBoolValues5;
  std::map<std::string, double> ListDoubleValues5;
  std::map<std::string, std::string> ListStringValues5;
  std::map<std::string, std::vector<std::string>> ListListStringValues5;
  ListBoolValues5["LBCSE"]=false;
  ListBoolValues5["LBINTER"]=false;
  ListBoolValues5["LBCWA"]=true;
  ListBoolValues5["LINHOM"]=false;
  ListBoolValues5["LBCSP"]=false;
  ListBoolValues5["LINDSPRDEG"]=false;
  ListBoolValues5["LPARMDIR"]=false;
  ListStringValues5["FILEWAVE"]="wwmbnd.gr3";
  ListBoolValues5["LBSP1D"]=false;
  ListBoolValues5["LBSP2D"]=false;
  ListDoubleValues5["DELTC"]=600;
  ListStringValues5["UNITC"]="SEC";
  ListStringValues5["FILEBOUND"]="wwmbnd.gr3";
  ListIntValues5["IBOUNDFORMAT"]=1;
  ListDoubleValues5["WBHS"]=4.;
  ListDoubleValues5["WBSS"]=2.;
  ListDoubleValues5["WBTP"]=8.;
  ListDoubleValues5["WBDM"]=90.0;
  ListDoubleValues5["WBDSMS"]=1.;
  ListDoubleValues5["WBDS"]=20.;
  ListDoubleValues5["WBGAUSS"]=0.1;
  ListDoubleValues5["WBPKEN"]=3.3;
  ListStringValues5["NCDF_HS_NAME"]="hs";
  ListStringValues5["NCDF_DIR_NAME"]="dir";
  ListStringValues5["NCDF_SPR_NAME"]="spr";
  ListStringValues5["NCDF_FP_NAME"]="fp";
  ListStringValues5["NCDF_F02_NAME"]="t02";
  SingleBlock BlockBOUC;
  BlockBOUC.ListIntValues=ListIntValues5;
  BlockBOUC.ListBoolValues=ListBoolValues5;
  BlockBOUC.ListDoubleValues=ListDoubleValues5;
  BlockBOUC.ListStringValues=ListStringValues5;
  BlockBOUC.ListListStringValues=ListListStringValues5;
  ListBlock["BOUC"]=BlockBOUC;
  // WIND
  std::map<std::string, int> ListIntValues6;
  std::map<std::string, bool> ListBoolValues6;
  std::map<std::string, double> ListDoubleValues6;
  std::map<std::string, std::string> ListStringValues6;
  std::map<std::string, std::vector<std::string>> ListListStringValues6;
  ListBoolValues6["LSEWD"]=false;
  ListStringValues6["BEGTC"]="20121218.120000";
  ListStringValues6["ENDTC"]="20121223.000000";
  ListDoubleValues6["DELTC"]=3;
  ListStringValues6["UNITC"]="HR";
  ListBoolValues6["LINTERWD"]=true;
  ListBoolValues6["LSTWD"]=false;
  ListBoolValues6["LCWIN"]=false;
  ListBoolValues6["LWDIR"]=true;
  ListDoubleValues6["WDIR"]=120.0;
  ListDoubleValues6["WVEL"]=30.0;
  ListDoubleValues6["CWINDX"]=0.0;
  ListDoubleValues6["CWINDY"]=0.0;
  ListStringValues6["FILEWIND"]="wind1.dat";
  SingleBlock BlockWIND;
  BlockWIND.ListIntValues=ListIntValues6;
  BlockWIND.ListBoolValues=ListBoolValues6;
  BlockWIND.ListDoubleValues=ListDoubleValues6;
  BlockWIND.ListStringValues=ListStringValues6;
  BlockWIND.ListListStringValues=ListListStringValues6;
  ListBlock["WIND"]=BlockWIND;
  // CURR
  std::map<std::string, int> ListIntValues7;
  std::map<std::string, bool> ListBoolValues7;
  std::map<std::string, double> ListDoubleValues7;
  std::map<std::string, std::string> ListStringValues7;
  std::map<std::string, std::vector<std::string>> ListListStringValues7;
  ListBoolValues7["LSECU"]=false;
  ListStringValues7["BEGTC"]="20121218.120000";
  ListDoubleValues7["DELTC"]=600;
  ListStringValues7["UNITC"]="SEC";
  ListStringValues7["ENDTC"]="20121223.000000";
  ListBoolValues7["LINTERCU"]=false;
  ListBoolValues7["LSTCU"]=false;
  ListBoolValues7["LCCUR"]=false;
  ListDoubleValues7["CCURTX"]=0.0;
  ListDoubleValues7["CCURTY"]=0.0;
  ListStringValues7["FILECUR"]="current.dat";
  ListBoolValues7["LERGINP"]=false;
  ListDoubleValues7["CURFAC"]=1.000000;
  ListIntValues7["ICURRFORMAT"]=1;
  SingleBlock BlockCURR;
  BlockCURR.ListIntValues=ListIntValues7;
  BlockCURR.ListBoolValues=ListBoolValues7;
  BlockCURR.ListDoubleValues=ListDoubleValues7;
  BlockCURR.ListStringValues=ListStringValues7;
  BlockCURR.ListListStringValues=ListListStringValues7;
  ListBlock["CURR"]=BlockCURR;
  // WALV
  std::map<std::string, int> ListIntValues8;
  std::map<std::string, bool> ListBoolValues8;
  std::map<std::string, double> ListDoubleValues8;
  std::map<std::string, std::string> ListStringValues8;
  std::map<std::string, std::vector<std::string>> ListListStringValues8;
  ListBoolValues8["LSEWL"]=false;
  ListStringValues8["BEGTC"]=" ";
  ListStringValues8["ENDTC"]=" ";
  ListDoubleValues8["DELTC"]=1;
  ListStringValues8["UNITC"]="HR";
  ListBoolValues8["LINTERWL"]=false;
  ListBoolValues8["LSTWL"]=false;
  ListBoolValues8["LCWLV"]=false;
  ListDoubleValues8["CWATLV"]=0.0;
  ListStringValues8["FILEWATL"]=" ";
  SingleBlock BlockWALV;
  BlockWALV.ListIntValues=ListIntValues8;
  BlockWALV.ListBoolValues=ListBoolValues8;
  BlockWALV.ListDoubleValues=ListDoubleValues8;
  BlockWALV.ListStringValues=ListStringValues8;
  BlockWALV.ListListStringValues=ListListStringValues8;
  ListBlock["WALV"]=BlockWALV;
  // ENGS
  std::map<std::string, int> ListIntValues9;
  std::map<std::string, bool> ListBoolValues9;
  std::map<std::string, double> ListDoubleValues9;
  std::map<std::string, std::string> ListStringValues9;
  std::map<std::string, std::vector<std::string>> ListListStringValues9;
  ListIntValues9["MESNL"]=0;
  ListIntValues9["MESIN"]=1;
  ListIntValues9["IFRIC"]=1;
  ListIntValues9["MESBF"]=1;
  ListDoubleValues9["FRICC"]=0.067;
  ListIntValues9["MESBR"]=1;
  ListIntValues9["ICRIT"]=1;
  ListDoubleValues9["ALPBJ"]=1.;
  ListDoubleValues9["BRHD"]=0.78;
  ListBoolValues9["LMAXETOT"]=true;
  ListIntValues9["MESDS"]=1;
  ListIntValues9["MESTR"]=0;
  ListDoubleValues9["TRICO"]=1.;
  ListDoubleValues9["TRIRA"]=2.5;
  ListDoubleValues9["TRIURS"]=0.1;
  SingleBlock BlockENGS;
  BlockENGS.ListIntValues=ListIntValues9;
  BlockENGS.ListBoolValues=ListBoolValues9;
  BlockENGS.ListDoubleValues=ListDoubleValues9;
  BlockENGS.ListStringValues=ListStringValues9;
  BlockENGS.ListListStringValues=ListListStringValues9;
  ListBlock["ENGS"]=BlockENGS;
  // NUMS
  std::map<std::string, int> ListIntValues10;
  std::map<std::string, bool> ListBoolValues10;
  std::map<std::string, double> ListDoubleValues10;
  std::map<std::string, std::string> ListStringValues10;
  std::map<std::string, std::vector<std::string>> ListListStringValues10;
  ListIntValues10["ICOMP"]=0;
  ListIntValues10["AMETHOD"]=0;
  ListIntValues10["SMETHOD"]=2;
  ListIntValues10["DMETHOD"]=2;
  ListDoubleValues10["RTHETA"]=0.5;
  ListBoolValues10["LITERSPLIT"]=false;
  ListBoolValues10["LFILTERTH"]=false;
  ListDoubleValues10["MAXCFLTH"]=1.0;
  ListIntValues10["FMETHOD"]=1;
  ListBoolValues10["LFILTERSIG"]=false;
  ListDoubleValues10["MAXCFLSIG"]=1.0;
  ListBoolValues10["LSIGBOUND"]=false;
  ListBoolValues10["LTHBOUND"]=false;
  ListBoolValues10["LSOUBOUND"]=false;
  ListBoolValues10["LLIMT"]=true;
  ListIntValues10["MELIM"]=1;
  ListDoubleValues10["LIMFAK"]=0.1;
  ListBoolValues10["LDIFR"]=false;
  ListIntValues10["IDIFFR"]=1;
  ListBoolValues10["LCONV"]=false;
  ListBoolValues10["LCFL"]=false;
  ListIntValues10["NQSITER"]=10;
  ListDoubleValues10["QSCONV1"]=0.98;
  ListDoubleValues10["QSCONV2"]=0.98;
  ListDoubleValues10["QSCONV3"]=0.98;
  ListDoubleValues10["QSCONV4"]=0.98;
  ListDoubleValues10["QSCONV5"]=0.98;
  ListBoolValues10["LEXPIMP"]=false;
  ListDoubleValues10["FREQEXP"]=0.1;
  ListDoubleValues10["EPSH1"]=0.01;
  ListDoubleValues10["EPSH2"]=0.01;
  ListDoubleValues10["EPSH3"]=0.01;
  ListDoubleValues10["EPSH4"]=0.01;
  ListDoubleValues10["EPSH5"]=0.01;
  ListBoolValues10["LVECTOR"]=false;
  ListIntValues10["IVECTOR"]=2;
  ListBoolValues10["LADVTEST"]=false;
  ListBoolValues10["LCHKCONV"]=false;
  SingleBlock BlockNUMS;
  BlockNUMS.ListIntValues=ListIntValues10;
  BlockNUMS.ListBoolValues=ListBoolValues10;
  BlockNUMS.ListDoubleValues=ListDoubleValues10;
  BlockNUMS.ListStringValues=ListStringValues10;
  BlockNUMS.ListListStringValues=ListListStringValues10;
  ListBlock["NUMS"]=BlockNUMS;
  // HISTORY
  std::map<std::string, int> ListIntValues11;
  std::map<std::string, bool> ListBoolValues11;
  std::map<std::string, double> ListDoubleValues11;
  std::map<std::string, std::string> ListStringValues11;
  std::map<std::string, std::vector<std::string>> ListListStringValues11;
  ListStringValues11["BEGTC"]="20110915.000000";
  ListDoubleValues11["DELTC"]=300;
  ListStringValues11["UNITC"]="SEC";
  ListStringValues11["ENDTC"]="20110925.000000";
  ListIntValues11["DEFINETC"]=0;
  ListStringValues11["OUTSTYLE"]="NC";
  ListIntValues11["MULTIPLEOUT"]=0;
  ListBoolValues11["USE_SINGLE_OUT"]=true;
  ListBoolValues11["PARAMWRITE"]=true;
  ListBoolValues11["GRIDWRITE"]=true;
  ListBoolValues11["PRINTMMA"]=false;
  ListStringValues11["FILEOUT"]="field.dat";
  ListBoolValues11["LWXFN"]=false;
  ListBoolValues11["HS"]=false;
  ListBoolValues11["TM01"]=false;
  ListBoolValues11["TM10"]=false;
  ListBoolValues11["TM02"]=false;
  ListBoolValues11["KLM"]=false;
  ListBoolValues11["WLM"]=false;
  ListBoolValues11["ETOTC"]=false;
  ListBoolValues11["ETOTS"]=false;
  ListBoolValues11["DM"]=false;
  ListBoolValues11["DSPR"]=false;
  ListBoolValues11["TPPD"]=false;
  ListBoolValues11["TPP"]=false;
  ListBoolValues11["CPP"]=false;
  ListBoolValues11["WNPP"]=false;
  ListBoolValues11["CGPD"]=false;
  ListBoolValues11["CGPP"]=false;
  ListBoolValues11["CPPD"]=false;
  ListBoolValues11["KPPD"]=false;
  ListBoolValues11["KPP"]=false;
  ListBoolValues11["LPP"]=false;
  ListBoolValues11["PEAKD"]=false;
  ListBoolValues11["PEAKDSPR"]=false;
  ListBoolValues11["DPEAK"]=false;
  ListBoolValues11["UBOT"]=false;
  ListBoolValues11["ORBITAL"]=false;
  ListBoolValues11["BOTEXPER"]=false;
  ListBoolValues11["TMBOT"]=false;
  ListBoolValues11["URSELL"]=false;
  ListBoolValues11["UFRIC"]=false;
  ListBoolValues11["Z0"]=false;
  ListBoolValues11["ALPHA_CH"]=false;
  ListBoolValues11["WINDX"]=false;
  ListBoolValues11["WINDY"]=false;
  ListBoolValues11["WINDMAG"]=false;
  ListBoolValues11["CD"]=false;
  ListBoolValues11["CURRTX"]=false;
  ListBoolValues11["CURRTY"]=false;
  ListBoolValues11["WATLEV"]=false;
  ListBoolValues11["WATLEVOLD"]=false;
  ListBoolValues11["DEP"]=false;
  ListBoolValues11["DEPDT"]=false;
  ListBoolValues11["TAUW"]=false;
  ListBoolValues11["TAUWX"]=false;
  ListBoolValues11["TAUWY"]=false;
  ListBoolValues11["TAUHF"]=false;
  ListBoolValues11["TAUTOT"]=false;
  ListBoolValues11["STOKESSURFX"]=false;
  ListBoolValues11["STOKESSURFY"]=false;
  ListBoolValues11["STOKESBAROX"]=false;
  ListBoolValues11["STOKESBAROY"]=false;
  ListBoolValues11["RSXX"]=false;
  ListBoolValues11["RSXY"]=false;
  ListBoolValues11["RSYY"]=false;
  ListBoolValues11["CFL1"]=false;
  ListBoolValues11["CFL2"]=false;
  ListBoolValues11["CFL3"]=false;
  SingleBlock BlockHISTORY;
  BlockHISTORY.ListIntValues=ListIntValues11;
  BlockHISTORY.ListBoolValues=ListBoolValues11;
  BlockHISTORY.ListDoubleValues=ListDoubleValues11;
  BlockHISTORY.ListStringValues=ListStringValues11;
  BlockHISTORY.ListListStringValues=ListListStringValues11;
  ListBlock["HISTORY"]=BlockHISTORY;
  // STATION
  std::map<std::string, int> ListIntValues12;
  std::map<std::string, bool> ListBoolValues12;
  std::map<std::string, double> ListDoubleValues12;
  std::map<std::string, std::vector<double>> ListListDoubleValues12;
  std::map<std::string, std::string> ListStringValues12;
  std::map<std::string, std::vector<std::string>> ListListStringValues12;
  ListDoubleValues12["DELTC"]=600;
  ListStringValues12["UNITC"]="SEC";
  ListStringValues12["OUTSTYLE"]="NO";
  ListStringValues12["FILEOUT"]="station_single.dat";
  ListBoolValues12["LOUTITER"]=false;
  ListBoolValues12["LLOUTS"]=false;
  ListIntValues12["ILOUTS"]=1;
  std::vector<std::string> eListStr={"P-1"};
  std::vector<double> eListX={1950};
  std::vector<double> eListY={304};
  std::vector<double> eListCutoff={0.44};
  ListListStringValues12["NLOUTS"]=eListStr;
  ListListStringValues12["NOUTS"]=eListStr;
  ListIntValues12["IOUTS"]=1;
  ListListDoubleValues12["XOUTS"]=eListX;
  ListListDoubleValues12["YOUTS"]=eListY;
  ListListDoubleValues12["CUTOFF"]=eListCutoff;
  ListBoolValues12["LSP1D"]=false;
  ListBoolValues12["LSP2D"]=false;
  ListBoolValues12["LSIGMAX"]=true;
  ListBoolValues12["AC"]=false;
  ListBoolValues12["WK"]=false;
  ListBoolValues12["ACOUT_1D"]=false;
  ListBoolValues12["ACOUT_2D"]=false;
  ListBoolValues12["HS"]=false;
  ListBoolValues12["TM01"]=false;
  ListBoolValues12["TM10"]=false;
  ListBoolValues12["TM02"]=false;
  ListBoolValues12["KLM"]=false;
  ListBoolValues12["WLM"]=false;
  ListBoolValues12["ETOTC"]=false;
  ListBoolValues12["ETOTS"]=false;
  ListBoolValues12["DM"]=false;
  ListBoolValues12["DSPR"]=false;
  ListBoolValues12["TPPD"]=false;
  ListBoolValues12["TPP"]=false;
  ListBoolValues12["CPP"]=false;
  ListBoolValues12["WNPP"]=false;
  ListBoolValues12["CGPD"]=false;
  ListBoolValues12["CGPP"]=false;
  ListBoolValues12["CPPD"]=false;
  ListBoolValues12["KPPD"]=false;
  ListBoolValues12["KPP"]=false;
  ListBoolValues12["LPP"]=false;
  ListBoolValues12["PEAKD"]=false;
  ListBoolValues12["PEAKDSPR"]=false;
  ListBoolValues12["DPEAK"]=false;
  ListBoolValues12["UBOT"]=false;
  ListBoolValues12["ORBITAL"]=false;
  ListBoolValues12["BOTEXPER"]=false;
  ListBoolValues12["TMBOT"]=false;
  ListBoolValues12["URSELL"]=false;
  ListBoolValues12["UFRIC"]=false;
  ListBoolValues12["Z0"]=false;
  ListBoolValues12["ALPHA_CH"]=false;
  ListBoolValues12["WINDX"]=false;
  ListBoolValues12["WINDY"]=false;
  ListBoolValues12["WINDMAG"]=false;
  ListBoolValues12["CD"]=false;
  ListBoolValues12["CURRTX"]=false;
  ListBoolValues12["CURRTY"]=false;
  ListBoolValues12["WATLEV"]=false;
  ListBoolValues12["WATLEVOLD"]=false;
  ListBoolValues12["DEP"]=false;
  ListBoolValues12["DEPDT"]=false;
  ListBoolValues12["TAUW"]=false;
  ListBoolValues12["TAUWX"]=false;
  ListBoolValues12["TAUWY"]=false;
  ListBoolValues12["TAUHF"]=false;
  ListBoolValues12["TAUTOT"]=false;
  ListBoolValues12["STOKESSURFX"]=false;
  ListBoolValues12["STOKESSURFY"]=false;
  ListBoolValues12["STOKESBAROX"]=false;
  ListBoolValues12["STOKESBAROY"]=false;
  ListBoolValues12["RSXX"]=false;
  ListBoolValues12["RSXY"]=false;
  ListBoolValues12["RSYY"]=false;
  ListBoolValues12["CFL1"]=false;
  ListBoolValues12["CFL2"]=false;
  ListBoolValues12["CFL3"]=false;
  SingleBlock BlockSTATION;
  BlockSTATION.ListIntValues=ListIntValues12;
  BlockSTATION.ListBoolValues=ListBoolValues12;
  BlockSTATION.ListDoubleValues=ListDoubleValues12;
  BlockSTATION.ListListDoubleValues=ListListDoubleValues12;
  BlockSTATION.ListStringValues=ListStringValues12;
  BlockSTATION.ListListStringValues=ListListStringValues12;
  ListBlock["STATION"]=BlockSTATION;
  // HOTFILE
  std::map<std::string, int> ListIntValues13;
  std::map<std::string, bool> ListBoolValues13;
  std::map<std::string, double> ListDoubleValues13;
  std::map<std::string, std::string> ListStringValues13;
  std::map<std::string, std::vector<std::string>> ListListStringValues13;
  ListBoolValues13["LHOTF"]=false;
  ListStringValues13["BEGTC"]="19980901.000000";
  ListDoubleValues13["DELTC"]=5;
  ListStringValues13["UNITC"]="SEC";
  ListStringValues13["ENDTC"]="19980901.010000";
  ListBoolValues13["LCYCLEHOT"]=true;
  ListIntValues13["HOTSTYLE_OUT"]=2;
  ListIntValues13["MULTIPLEOUT"]=0;
  ListStringValues13["FILEHOT_OUT"]="hotfile_out.dat";
  ListIntValues13["HOTSTYLE_IN"]=2;
  ListIntValues13["IHOTPOS_IN"]=1;
  ListIntValues13["MULTIPLEIN"]=0;
  ListStringValues13["FILEHOT_IN"]="hotfile_in.dat";
  SingleBlock BlockHOTFILE;
  BlockHOTFILE.ListIntValues=ListIntValues13;
  BlockHOTFILE.ListBoolValues=ListBoolValues13;
  BlockHOTFILE.ListDoubleValues=ListDoubleValues13;
  BlockHOTFILE.ListStringValues=ListStringValues13;
  BlockHOTFILE.ListListStringValues=ListListStringValues13;
  ListBlock["HOTFILE"]=BlockHOTFILE;
  // PETScOptions
  std::map<std::string, int> ListIntValues14;
  std::map<std::string, bool> ListBoolValues14;
  std::map<std::string, double> ListDoubleValues14;
  std::map<std::string, std::string> ListStringValues14;
  std::map<std::string, std::vector<std::string>> ListListStringValues14;
  ListStringValues14["KSPTYPE"]="bcgsl";
  ListDoubleValues14["RTOL"]=1.E-10;
  ListDoubleValues14["ABSTOL"]=1.E-14;
  ListDoubleValues14["DTOL"]=10000.;
  ListIntValues14["MAXITS"]=1000;
  ListBoolValues14["INITIALGUESSNONZERO"]=false;
  ListBoolValues14["GMRESPREALLOCATE"]=true;
  ListStringValues14["PCTYPE"]="sor";
  SingleBlock BlockPETSC;
  BlockPETSC.ListIntValues=ListIntValues14;
  BlockPETSC.ListBoolValues=ListBoolValues14;
  BlockPETSC.ListDoubleValues=ListDoubleValues14;
  BlockPETSC.ListStringValues=ListStringValues14;
  BlockPETSC.ListListStringValues=ListListStringValues14;
  ListBlock["PETScOptions"]=BlockPETSC;
  return {std::move(ListBlock), "undefined"};
}


#endif
