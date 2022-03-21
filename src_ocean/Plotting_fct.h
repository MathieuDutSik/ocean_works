#ifndef PLOTTING_FCT_INCLUDE
#define PLOTTING_FCT_INCLUDE

#include "Interpolation.h"
#include "Namelist.h"
#include "Model_grids.h"
#include "Model_data_loading.h"
#include "Basic_plot.h"
#include "CommonFuncModel.h"
#include "Statistics.h"
#include "Kernel_Transect.h"

void LocateMM_FCT(MyMatrix<double> const& F, GridArray const& GrdArr, std::string const& VarName)
{
  int eta=F.rows();
  int xi=F.cols();
  double minval=0, maxval=0;
  int nb=0;
  bool IsFirst=true;
  int iMax=-1;
  int jMax=-1;
  int iMin=-1;
  int jMin=-1;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++)
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
        double eVal=F(i,j);
        nb++;
        if (IsFirst) {
          IsFirst=false;
          maxval=eVal;
          minval=eVal;
          iMax=i;
          jMax=j;
          iMin=i;
          jMin=j;
        } else {
          if (eVal > maxval) {
            maxval=eVal;
            iMax=i;
            jMax=j;
          }
          if (eVal < minval) {
            minval=eVal;
            iMin=i;
            jMin=j;
          }
        }
      }
  int OptChoice=3;
  if (OptChoice == 1) {
    int nbChar=VarName.size();
    std::cerr << "  " << VarName << " min(p,i,j)=(" << minval << "," << iMin << "," << jMin << ") lon=" << GrdArr.GrdArrRho.LON(iMin,jMin) << " lat=" << GrdArr.GrdArrRho.LAT(iMin,jMin) << "\n";
    std::cerr << "  ";
    for (int iChar=0; iChar<nbChar; iChar++)
      std::cerr << " ";
    std::cerr << " max(p,i,j)=(" << maxval << "," << iMax << "," << jMax << ") lon=" << GrdArr.GrdArrRho.LON(iMax,jMax) << " lat=" << GrdArr.GrdArrRho.LAT(iMax,jMax) << "\n";
  }
  if (OptChoice == 2) {
    std::cerr << " mm " << VarName << " min(p,i,j)=(" << minval << "," << iMin << "," << jMin << ") lon=" << GrdArr.GrdArrRho.LON(iMin,jMin) << " lat=" << GrdArr.GrdArrRho.LAT(iMin,jMin) << "\n";
    std::cerr << "    " << VarName << " max(p,i,j)=(" << maxval << "," << iMax << "," << jMax << ") lon=" << GrdArr.GrdArrRho.LON(iMax,jMax) << " lat=" << GrdArr.GrdArrRho.LAT(iMax,jMax) << "\n";
  }
  if (OptChoice == 3) {
    std::cerr << " min(p,i,j)=(" << minval << "," << iMin << "," << jMin << ") lon=" << GrdArr.GrdArrRho.LON(iMin,jMin) << " lat=" << GrdArr.GrdArrRho.LAT(iMin,jMin);
    std::cerr << "    max(p,i,j)=(" << maxval << "," << iMax << "," << jMax << ") lon=" << GrdArr.GrdArrRho.LON(iMax,jMax) << " lat=" << GrdArr.GrdArrRho.LAT(iMax,jMax) << "\n";
  }
}




std::string RetrieveVarName2(RecSymbolic const& RecS, FullNamelist const& eFull)
{
  std::string VarName1=RecS.VarName1;
  SingleBlock BlPLOT=eFull.ListBlock.at("PLOT");
  std::vector<std::string> RenameVariable_VarName1=BlPLOT.ListListStringValues.at("RenameVariable_VarName1");
  std::vector<std::string> RenameVariable_VarName2=BlPLOT.ListListStringValues.at("RenameVariable_VarName2");
  int len1=RenameVariable_VarName1.size();
  int len2=RenameVariable_VarName2.size();
  if (len1 != len2) {
    std::cerr << "RenameVariable_VarName1/2 are not of the same length. Illegal\n";
    throw TerminalException{1};
  }
  for (int i=0; i<len1; i++)
    if (RenameVariable_VarName1[i] == VarName1)
      return RenameVariable_VarName2[i];
  return RecS.VarName2;
}

std::string RetrieveVarNameUF_kernel(RecSymbolic const& RecS, FullNamelist const& eFull)
{
  std::string VarName1=RecS.VarName1;
  SingleBlock BlPLOT=eFull.ListBlock.at("PLOT");
  std::vector<std::string> RenameVariable_VarName1 =BlPLOT.ListListStringValues.at("RenameVariable_VarName1");
  std::vector<std::string> RenameVariable_VarNameUF=BlPLOT.ListListStringValues.at("RenameVariable_VarNameUF");
  int len1=RenameVariable_VarName1.size();
  int len2=RenameVariable_VarNameUF.size();
  if (len1 != len2) {
    std::cerr << "RenameVariable_VarName1/2 are not of the same length. Illegal\n";
    throw TerminalException{1};
  }
  for (int i=0; i<len1; i++)
    if (RenameVariable_VarName1[i] == VarName1)
      return RenameVariable_VarNameUF[i];
  return VarName1;
}


std::vector<QuadDrawInfo> GetListQuadArray(SingleBlock const& eBlPLOT, GridArray const& GrdArr)
{
  std::vector<double> ListFrameMinLon=eBlPLOT.ListListDoubleValues.at("ListFrameMinLon");
  std::vector<double> ListFrameMinLat=eBlPLOT.ListListDoubleValues.at("ListFrameMinLat");
  std::vector<double> ListFrameMaxLon=eBlPLOT.ListListDoubleValues.at("ListFrameMaxLon");
  std::vector<double> ListFrameMaxLat=eBlPLOT.ListListDoubleValues.at("ListFrameMaxLat");
  size_t nbFrame=ListFrameMinLon.size();
  if (nbFrame != ListFrameMinLat.size() || nbFrame != ListFrameMaxLon.size() || nbFrame != ListFrameMaxLat.size()) {
    std::cerr << "Error the ListFrame(Min/Max)(Lon/Lat) arrays should all be of the same size\n";
    throw TerminalException{1};
  }
  std::vector<QuadDrawInfo> RetList;
  int iFrameIdx=0;
  if (eBlPLOT.ListBoolValues.at("DoMain")) {
    RetList.push_back({"Main", iFrameIdx, GetQuadArray(GrdArr)});
    iFrameIdx++;
  }
  for (int iFrame=0; iFrame<int(nbFrame); iFrame++) {
    double MinLon=ListFrameMinLon[iFrame];
    double MinLat=ListFrameMinLat[iFrame];
    double MaxLon=ListFrameMaxLon[iFrame];
    double MaxLat=ListFrameMaxLat[iFrame];
    QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};
    std::string eFrameName = "fr" + std::to_string(iFrame);
    RetList.push_back({eFrameName, iFrameIdx, eQuad});
    iFrameIdx++;
  }
  int nbFrameTot=RetList.size();
  if (eBlPLOT.ListListDoubleValues.count("TransectLonStart") > 0) {
    int nbTrans=eBlPLOT.ListListDoubleValues.at("TransectLonStart").size();
    if (nbFrameTot == 0 && nbTrans == 0) {
      std::cerr << "It might not be clever to call the plotting software\n";
      std::cerr << "with zero frames selected and zero transect selected\n";
      std::cerr << "In section PLOT\n";
      std::cerr << "edit the variables DoMain, \n";
      std::cerr << "and/or ListFrameMinLon, ListFrameMinLat, ListFrameMaxLon, ListFrameMaxLat\n";
      throw TerminalException{1};
    }
    if (nbFrameTot > 1) {
      for (int iFrameTot=0; iFrameTot<nbFrameTot; iFrameTot++) {
        RetList[iFrameTot].iFrame = iFrameTot;
        RetList[iFrameTot].eFrameName = "_fr" + StringNumber(iFrameTot, 2);
      }
    }
  }
  return RetList;
}





DrawArr CommonAssignation_DrawArr(FullNamelist const& eFull)
{
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  DrawArr eDrawArr;
  eDrawArr.nbLevelSpa=eBlPLOT.ListIntValues.at("nbLevelSpa");
  eDrawArr.nbLabelStride=eBlPLOT.ListIntValues.at("nbLabelStride");
  if (eBlPLOT.ListBoolValues.count("DrawRiver") == 0) {
    eDrawArr.DrawRiver=false;
  } else {
    eDrawArr.DrawRiver=eBlPLOT.ListBoolValues.at("DrawRiver");
  }
  eDrawArr.TheAnnot.DrawAnnotation=eBlPLOT.ListBoolValues.at("DrawAnnotation");
  eDrawArr.TheAnnot.AnnotationLon=eBlPLOT.ListDoubleValues.at("AnnotationLon");
  eDrawArr.TheAnnot.AnnotationLat=eBlPLOT.ListDoubleValues.at("AnnotationLat");
  eDrawArr.TheAnnot.AnnotationText=eBlPLOT.ListStringValues.at("AnnotationText");
  std::string FileDirectNCLins = eBlPLOT.ListStringValues.at("FileDirectNCLins");
  if (FileDirectNCLins != "irrelevant") {
    eDrawArr.ListInsertLines = ReadFullFile(FileDirectNCLins);
  }
  if (eBlPLOT.ListBoolValues.count("FillLand") > 0) {
    eDrawArr.FillLand=eBlPLOT.ListBoolValues.at("FillLand");
  } else {
    eDrawArr.FillLand=false;
  }
  if (eBlPLOT.ListStringValues.count("GridResolution") > 0) {
    eDrawArr.GridResolution=eBlPLOT.ListStringValues.at("GridResolution");
  } else {
    eDrawArr.GridResolution="unset";
  }
  if (eBlPLOT.ListBoolValues.count("UseNativeGrid") > 0) {
    eDrawArr.UseNativeGrid=eBlPLOT.ListBoolValues.at("UseNativeGrid");
  } else {
    eDrawArr.UseNativeGrid=true;
  }
  eDrawArr.DoTitle=eBlPLOT.ListBoolValues.at("DoTitle");
  eDrawArr.vcRefLengthF=eBlPLOT.ListDoubleValues.at("vcRefLengthF");
  eDrawArr.DoColorBar=eBlPLOT.ListBoolValues.at("DoColorBar");
  eDrawArr.cnFillMode=eBlPLOT.ListStringValues.at("cnFillMode");
  eDrawArr.cnFillOn=eBlPLOT.ListBoolValues.at("cnFillOn");
  eDrawArr.cnLinesOn=eBlPLOT.ListBoolValues.at("cnLinesOn");
  eDrawArr.cnLineLabelsOn=eBlPLOT.ListBoolValues.at("cnLineLabelsOn");
  eDrawArr.cnSmoothingOn=eBlPLOT.ListBoolValues.at("cnSmoothingOn");
  eDrawArr.ColorMap=eBlPLOT.ListStringValues.at("ColorMap");
  if (eBlPLOT.ListBoolValues.count("DrawContourBathy") > 0) {
    eDrawArr.DrawContourBathy=eBlPLOT.ListBoolValues.at("DrawContourBathy");
  } else {
    eDrawArr.DrawContourBathy=false;
  }
  eDrawArr.PrintMMA=eBlPLOT.ListBoolValues.at("PrintMMA");
  eDrawArr.DoTitleString=eBlPLOT.ListBoolValues.at("DoTitleString");
  eDrawArr.LandPortr=eBlPLOT.ListStringValues.at("LandPortr");
  eDrawArr.optStatStr=eBlPLOT.ListStringValues.at("optStatStr");
  return eDrawArr;
}




void PLOT_DIFF_FD_RHO_PCOLOR(GridArray const& GrdArr,
			     RecVar const& eRecVar1,
			     RecVar const& eRecVar2,
			     NCLcaller<GeneralType> & eCall,
			     PermanentInfoDrawing const& ePerm)
{
  SingleBlock eBlPLOT=ePerm.eFull.ListBlock.at("PLOT");
  SingleBlock eBlPROC=ePerm.eFull.ListBlock.at("PROC");
  DrawArr eDrawArr=ePerm.eDrawArr;
  std::string Name1=eBlPROC.ListStringValues.at("Name1");
  std::string Name2=eBlPROC.ListStringValues.at("Name2");
  MyMatrix<double> eDiff = eRecVar1.F - eRecVar2.F;
  bool PrintMMA=eBlPLOT.ListBoolValues.at("PrintMMA");
  if (PrintMMA) {
    std::cerr << "diff :";
    PrintMMA_FCT(eDiff, GrdArr.GrdArrRho.MSK, eRecVar1.RecS.VarName1, eRecVar1.RecS.Unit);
    std::cerr << "  F1 :";
    PrintMMA_FCT(eRecVar1.F, GrdArr.GrdArrRho.MSK, eRecVar1.RecS.VarName1, eRecVar1.RecS.Unit);
    std::cerr << "  F2 :";
    PrintMMA_FCT(eRecVar2.F, GrdArr.GrdArrRho.MSK, eRecVar1.RecS.VarName1, eRecVar1.RecS.Unit);
  }
  bool LocateMM=eBlPLOT.ListBoolValues.at("LocateMM");
  if (LocateMM) {
    std::cerr << "diff :";
    LocateMM_FCT(eDiff, GrdArr, eRecVar1.RecS.VarName1);
    std::cerr << "  F1 :";
    LocateMM_FCT(eRecVar1.F, GrdArr, eRecVar1.RecS.VarName1);
    std::cerr << "  F2 :";
    LocateMM_FCT(eRecVar2.F, GrdArr, eRecVar1.RecS.VarName1);
  }
  std::string strColorMapDiff=eBlPLOT.ListStringValues.at("ColorMapDiff");
  eDrawArr.ColorMap=strColorMapDiff;
  RecVar eRecVar=eRecVar1;
  eRecVar.F=eDiff;
  eRecVar.RecS.minval=eRecVar1.RecS.mindiff;
  eRecVar.RecS.maxval=eRecVar1.RecS.maxdiff;
  for (auto & eQuadInfo : ePerm.ListQuadInfo) {
    std::string VarName2=RetrieveVarName2(eRecVar1.RecS, ePerm.eFull);
    std::string TitleStr="diff. " + VarName2 + " (" + Name1 + "-" + Name2 + ")";
    TitleStr += " " + eRecVar1.RecS.strPres;
    eDrawArr.TitleStr=TitleStr;
    std::string VarNameUF = RetrieveVarNameUF_kernel(eRecVar1.RecS, ePerm.eFull);
    //    std::cerr << "VarNameUF=" << VarNameUF << "\n";
    eDrawArr.VarNameUF = VarNameUF + eQuadInfo.eFrameName;
    std::string FileName = ePerm.eDir + eDrawArr.VarNameUF + "_" + eRecVar1.RecS.strAll;
    //    bool testFile=IsExistingFile(FileName + "." + ePerm.Extension);
    //    std::cerr << "testFile=" << testFile << "\n";
    if (!IsExistingFile(FileName) || eBlPROC.ListBoolValues.at("OverwritePrevious")) {
      eDrawArr.eQuadFrame = eQuadInfo.eQuad;
      PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    }
  }
}


void SINGLE_PLOT_QUIVER(GridArray const& GrdArr,
			RecVar const& eRecVar,
			NCLcaller<GeneralType> & eCall,
			PermanentInfoDrawing const& ePerm)
{
  DrawArr eDrawArr=ePerm.eDrawArr;
  bool test=ePerm.eFull.ListBlock.at("PLOT").ListBoolValues.at("ExcludeLargeValues");
  if (test) {
    double eMax=eRecVar.F.maxCoeff();
    double Thr=ePerm.eFull.ListBlock.at("PLOT").ListDoubleValues.at("ThresholdExclusionPlot");
    if (eMax > Thr)
      return;
  }
  if (eDrawArr.PrintMMA) {
    PrintMMA_FCT(eRecVar.F, GrdArr.GrdArrRho.MSK, eRecVar.RecS.VarName1, eRecVar.RecS.Unit);
  }
  bool OverwritePrevious=ePerm.eFull.ListBlock.at("PROC").ListBoolValues.at("OverwritePrevious");
  for (auto & eQuadInfo : ePerm.ListQuadInfo) {
    std::string VarNameUF = RetrieveVarNameUF_kernel(eRecVar.RecS, ePerm.eFull);
    eDrawArr.VarNameUF = VarNameUF + eQuadInfo.eFrameName;
    std::string FileName=ePerm.eDir + eDrawArr.VarNameUF + "_" + eRecVar.RecS.strAll;
    bool testFile=IsExistingFile(FileName + "." + ePerm.Extension);
    if (!testFile || OverwritePrevious) {
      std::string VarName2=RetrieveVarName2(eRecVar.RecS, ePerm.eFull);
      eDrawArr.TitleStr=VarName2 + " " + eRecVar.RecS.strPres;
      eDrawArr.eQuadFrame = eQuadInfo.eQuad;
      if (GrdArr.IsFE == 0) {
	//	std::cerr << "Call to PLOT_QUIVER, case 1\n";
	PLOT_QUIVER(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
      } else {
	int iFrame=eQuadInfo.iFrame;
	MyMatrix<double> U=SingleInterpolationOfField_2D(ePerm.ListInterpol[iFrame].InterpArr, eRecVar.U);
	MyMatrix<double> V=SingleInterpolationOfField_2D(ePerm.ListInterpol[iFrame].InterpArr, eRecVar.V);
	MyMatrix<double> F=SingleInterpolationOfField_2D(ePerm.ListInterpol[iFrame].InterpArr, eRecVar.F);
	RecVar NewRecVar;
	NewRecVar.RecS = eRecVar.RecS;
	NewRecVar.U=U;
	NewRecVar.V=V;
	NewRecVar.F=F;
	//	std::cerr << "Call to PLOT_QUIVER, case 2\n";
	PLOT_QUIVER(FileName, ePerm.ListInterpol[iFrame].GrdArr, eDrawArr, NewRecVar, eCall, ePerm);
      }
    }
  }
}



void SINGLE_PLOT_PCOLOR(const GridArray& GrdArr, const RecVar& eRecVar, NCLcaller<GeneralType> & eCall, const PermanentInfoDrawing& ePerm)
{
  DrawArr eDrawArr=ePerm.eDrawArr;
  bool test=ePerm.eFull.ListBlock.at("PLOT").ListBoolValues.at("ExcludeLargeValues");
  if (test) {
    double eMax=eRecVar.F.maxCoeff();
    double Thr=ePerm.eFull.ListBlock.at("PLOT").ListDoubleValues.at("ThresholdExclusionPlot");
    if (eMax > Thr)
      return;
  }
  if (eDrawArr.PrintMMA) {
    PrintMMA_FCT(eRecVar.F, GrdArr.GrdArrRho.MSK, eRecVar.RecS.VarName1, eRecVar.RecS.Unit);
  }
  bool OverwritePrevious=ePerm.eFull.ListBlock.at("PROC").ListBoolValues.at("OverwritePrevious");
  for (auto & eQuadInfo : ePerm.ListQuadInfo) {
    std::string VarNameUF = RetrieveVarNameUF_kernel(eRecVar.RecS, ePerm.eFull);
    eDrawArr.VarNameUF = VarNameUF + eQuadInfo.eFrameName;
    std::string FileName=ePerm.eDir + eDrawArr.VarNameUF + "_" + eRecVar.RecS.strAll;
    bool testFile=IsExistingFile(FileName + "." + ePerm.Extension);
    if (!testFile || OverwritePrevious) {
      std::string VarName2=RetrieveVarName2(eRecVar.RecS, ePerm.eFull);
      eDrawArr.TitleStr=VarName2 + " " + eRecVar.RecS.strPres;
      eDrawArr.eQuadFrame = eQuadInfo.eQuad;
      bool UseRegridArray=ePerm.eFull.ListBlock.at("PLOT").ListBoolValues.at("UseRegridArray");
      //      std::cerr << "UseRegridArray = " << UseRegridArray << "\n";
      if (GrdArr.IsFE == 0 && !UseRegridArray) {
        //        std::cerr << "Call to PLOT_PCOLOR, case 1\n";
	PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
      } else {
	int iFrame=eQuadInfo.iFrame;
        if (size_t(iFrame) >= ePerm.ListInterpol.size()) {
          std::cerr << "iFrame=" << iFrame << "\n";
          std::cerr << "|ePerm.ListInterpol|=" << ePerm.ListInterpol.size() << "\n";
          std::cerr << "This is a recurring problem in our code\n";
          throw TerminalException{1};
        }
	MyMatrix<double> F=SingleInterpolationOfField_2D(ePerm.ListInterpol[iFrame].InterpArr, eRecVar.F);
	RecVar NewRecVar;
	NewRecVar.RecS = eRecVar.RecS;
	NewRecVar.F=F;
        // std::cerr << "Call to PLOT_PCOLOR, case 2\n";
	PLOT_PCOLOR(FileName, ePerm.ListInterpol[iFrame].GrdArr, eDrawArr, NewRecVar, eCall, ePerm);
      }
    }
  }
}



void TRANSECT_PLOT_PCOLOR(TransectInformation_3D const& eTrans3,
			  RecVar const& eRecVar,
			  NCLcaller<GeneralType> & eCall,
			  PermanentInfoDrawing const& ePerm)
{
  RecVar NewRecVar;
  NewRecVar.RecS=eRecVar.RecS;
  if (eRecVar.RecS.VarNature == "3Drho") {
    NewRecVar.F=TransectInterpolation_3D(eTrans3, eRecVar.Tens3);
  }
  if (eRecVar.RecS.VarNature == "3Duv") {
    MyMatrix<double> U=TransectInterpolation_3D(eTrans3, eRecVar.Uthree);
    MyMatrix<double> V=TransectInterpolation_3D(eTrans3, eRecVar.Vthree);
    NewRecVar.F = eTrans3.normU * U + eTrans3.normV * V;
  }
  GridArray GrdArr=GetGridArrayFromTransect3(eTrans3);
  bool VariableRange=ePerm.eFull.ListBlock.at("PLOT").ListBoolValues.at("VariableRange");
  if (VariableRange) {
    PairMinMax ePair=ComputeMinMax(GrdArr, NewRecVar.F);
    NewRecVar.RecS.mindiff=ePair.TheMin;
    NewRecVar.RecS.maxdiff=ePair.TheMax;
    NewRecVar.RecS.minval=ePair.TheMin;
    NewRecVar.RecS.maxval=ePair.TheMax;
  }
  //
  DrawArr eDrawArr=ePerm.eDrawArr;
  eDrawArr.DrawContourBathy=false;
  eDrawArr.ListLineSegment.clear();
  eDrawArr.ListMarker.clear();
  std::string VarName2=RetrieveVarName2(eRecVar.RecS, ePerm.eFull);
  eDrawArr.TitleStr=VarName2 + " " + eRecVar.RecS.strPres;
  //
  std::string VarNameUF = RetrieveVarNameUF_kernel(eRecVar.RecS, ePerm.eFull);
  eDrawArr.VarNameUF = VarNameUF;
  std::string FileName=ePerm.eDir + eDrawArr.VarNameUF + "_" + eRecVar.RecS.strAll;
  //
  PLOT_PCOLOR(FileName, GrdArr, eDrawArr, NewRecVar, eCall, ePerm);
}



void GENERAL_PLOT_SINGLE(GridArray const& GrdArr,
			 RecVar const& eRecVar,
			 NCLcaller<GeneralType> & eCall,
			 PermanentInfoDrawing const& ePerm)
{
  SingleBlock eBlPLOT=ePerm.eFull.ListBlock.at("PLOT");
  bool LocateMM=eBlPLOT.ListBoolValues.at("LocateMM");
  if (eRecVar.RecS.VarNature == "uv") {
    SINGLE_PLOT_QUIVER(GrdArr, eRecVar, eCall, ePerm);
    return;
  }
  if (eRecVar.RecS.VarNature == "rho") {
    if (LocateMM) {
      LocateMM_FCT(eRecVar.F, GrdArr, eRecVar.RecS.VarName1);
    }
    SINGLE_PLOT_PCOLOR(GrdArr, eRecVar, eCall, ePerm);
    return;
  }
  if (eRecVar.RecS.VarNature == "3Drho") {
    auto LDim=eRecVar.Tens3.dimensions();
    int nbPlane=LDim[0];
    RecVar NewRecVar;
    NewRecVar.RecS=eRecVar.RecS;
    std::vector<int> ListPos=ePerm.eFull.ListBlock.at("PLOT").ListListIntValues.at("Tens3ListLevel");
    int siz=ListPos.size();
    if (siz == 0) {
      for (int iPlane=0; iPlane<nbPlane; iPlane++)
	ListPos.push_back(iPlane);
    }
    int nbDigit = GetNumberDigit(1 + VectorMax(ListPos));
    for (int& ePlane : ListPos) {
      NewRecVar.F=DimensionExtraction(eRecVar.Tens3, 0, ePlane);
      std::string ePlaneStr=StringNumber(ePlane+1, nbDigit);
      NewRecVar.RecS.VarName1=eRecVar.RecS.VarName1 + "_lev" + ePlaneStr;
      NewRecVar.RecS.VarName2=eRecVar.RecS.VarName2 + " (level " + ePlaneStr + ")";
      SINGLE_PLOT_PCOLOR(GrdArr, NewRecVar, eCall, ePerm);
    }
    int nbTrans=ePerm.ListTransect.size();
    for (int iTrans=0; iTrans<nbTrans; iTrans++) {
      RecVar NewRecVar=eRecVar;
      std::string eTransStr=StringNumber(iTrans+1, 2);
      std::string strAddi;
      if (nbTrans > 1) {
	strAddi=" (trans " + eTransStr + ")";
      } else {
	strAddi="";
      }
      NewRecVar.RecS.VarName1=eRecVar.RecS.VarName1;
      NewRecVar.RecS.VarName2=eRecVar.RecS.VarName2 + strAddi;
      TRANSECT_PLOT_PCOLOR(ePerm.ListTransect[iTrans], NewRecVar, eCall, ePerm);
    }
    return;
  }
  if (eRecVar.RecS.VarNature == "3Duv") {
    auto LDim=eRecVar.Tens3.dimensions();
    int nbPlane=LDim[0];
    RecVar NewRecVar;
    NewRecVar.RecS=eRecVar.RecS;
    std::vector<int> ListPos=ePerm.eFull.ListBlock.at("PLOT").ListListIntValues.at("Tens3ListLevel");
    int siz=ListPos.size();
    if (siz == 0) {
      for (int iPlane=0; iPlane<nbPlane; iPlane++)
	ListPos.push_back(iPlane);
    }
    int nbDigit = GetNumberDigit(1 + VectorMax(ListPos));
    for (int& ePlane : ListPos) {
      NewRecVar.U=DimensionExtraction(eRecVar.Uthree, 0, ePlane);
      NewRecVar.V=DimensionExtraction(eRecVar.Vthree, 0, ePlane);
      NewRecVar.F=DimensionExtraction(eRecVar.Tens3 , 0, ePlane);
      std::string ePlaneStr=StringNumber(ePlane+1, nbDigit);
      NewRecVar.RecS.VarName1=eRecVar.RecS.VarName1 + "_lev" + ePlaneStr;
      NewRecVar.RecS.VarName2=eRecVar.RecS.VarName2 + " (level " + ePlaneStr + ")";
      SINGLE_PLOT_QUIVER(GrdArr, NewRecVar, eCall, ePerm);
    }
    int nbTrans=ePerm.ListTransect.size();
    for (int iTrans=0; iTrans<nbTrans; iTrans++) {
      RecVar NewRecVar=eRecVar;
      std::string eTransStr=StringNumber(iTrans+1, 2);
      NewRecVar.RecS.VarName1=eRecVar.RecS.VarName1 + "_trans" + eTransStr;
      NewRecVar.RecS.VarName2=eRecVar.RecS.VarName2 + " (trans " + eTransStr + ")";
      TRANSECT_PLOT_PCOLOR(ePerm.ListTransect[iTrans], NewRecVar, eCall, ePerm);
    }
    return;
  }
  std::cerr << "No matching VarNature\n";
  std::cerr << "eRecVar.RecS.VarNature = " << eRecVar.RecS.VarNature << "\n";
  throw TerminalException{1};
}




void GRID_PLOTTING(GridArray const& GrdArr, std::string const& GridFile,
		   NCLcaller<GeneralType> & eCall,
		   PermanentInfoDrawing const& ePerm)
{
  SingleBlock eBlPLOT=ePerm.eFull.ListBlock.at("PLOT");
  DrawArr eDrawArr=ePerm.eDrawArr;
  eDrawArr.eQuadFrame=GetQuadArray(GrdArr);
  if (IsExistingFile(GridFile)) {
    if (NC_IsVar(GridFile, "MSK_att_ocn")) {
      MyMatrix<double> MSK_att_ocn=NC_Read2Dvariable(GridFile, "MSK_att_ocn");
      RecVar eRecVar;
      eRecVar.RecS.strAll="MASKattocn";
      eRecVar.RecS.VarName1="MSK_att_ocn";
      eRecVar.RecS.VarName2="Mask attainment oceanic model";
      eRecVar.RecS.minval=0;
      eRecVar.RecS.maxval=1;
      eRecVar.RecS.Unit="nondim.";
      eRecVar.F=MSK_att_ocn;
      std::string FileName=ePerm.eDir + "MSK_att_ocn";
      eDrawArr.TitleStr=eRecVar.RecS.VarName2;
      eDrawArr.VarNameUF="MSK_att_ocn";
      PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    }
    if (NC_IsVar(GridFile, "MSK_att_wav")) {
      MyMatrix<double> MSK_att_wav=NC_Read2Dvariable(GridFile, "MSK_att_wav");
      RecVar eRecVar;
      eRecVar.RecS.strAll="MSKattwav";
      eRecVar.RecS.VarName1="MSK_att_wav";
      eRecVar.RecS.VarName2="Mask attainment wav model";
      eRecVar.RecS.minval=0;
      eRecVar.RecS.maxval=1;
      eRecVar.RecS.Unit="nondim.";
      eRecVar.F=MSK_att_wav;
      std::string FileName=ePerm.eDir + "MSK_att_wav";
      eDrawArr.TitleStr=eRecVar.RecS.VarName2;
      eDrawArr.VarNameUF="MSK_att_wav";
      PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
    }
  }
  bool PlotDepth=eBlPLOT.ListBoolValues.at("PlotDepth");
  int SizeLON=GrdArr.GrdArrRho.LON.size();
  int SizeDEP=0;
  if (GrdArr.GrdArrRho.DEP)
    SizeDEP = GetDEP(GrdArr.GrdArrRho).size();
  if (PlotDepth && SizeLON == SizeDEP) {
    const MyMatrix<double> & DEP = GetDEP(GrdArr.GrdArrRho);
    PairMinMax ePair=ComputeMinMax(GrdArr, DEP);
    RecVar eRecVar;
    eRecVar.RecS.strAll="bathymetry";
    eRecVar.RecS.VarName1="Bathymetry";
    eRecVar.RecS.VarName2="Bathymetry of model";
    eRecVar.RecS.minval=ePair.TheMin;
    eRecVar.RecS.maxval=ePair.TheMax;
    eRecVar.RecS.Unit="m";
    eRecVar.F=DEP;
    std::string FileName=ePerm.eDir + "Bathymetry";
    eDrawArr.TitleStr=eRecVar.RecS.VarName2;
    eDrawArr.VarNameUF="Bathymetry";
    PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
  }
  bool PlotMesh=eBlPLOT.ListBoolValues.at("PlotMesh");
  if (PlotMesh && GrdArr.IsFE == 1) {
    PLOT_MESH(eDrawArr, GrdArr, eCall, ePerm);
    std::string eFileSVG=ePerm.eDir + "mesh.svg";
    DEFINE_MESH_SVG(eFileSVG, GrdArr);
  }
  bool PlotIOBP=eBlPLOT.ListBoolValues.at("PlotIOBP");
  if (PlotIOBP && GrdArr.IsFE == 1) {
    PairMinMax ePair{double(0), double(2)};
    RecVar eRecVar;
    eRecVar.RecS.strAll="IOBP";
    eRecVar.RecS.VarName1="IOBP";
    eRecVar.RecS.VarName2="IOBP";
    eRecVar.RecS.minval=ePair.TheMin;
    eRecVar.RecS.maxval=ePair.TheMax;
    eRecVar.RecS.Unit="m";
    int mnp=GrdArr.IOBP.size();
    MyMatrix<double> IOBP_f(mnp,1);
    for (int i=0; i<mnp; i++)
      IOBP_f(i,0) = double(GrdArr.IOBP(i));
    eRecVar.F=IOBP_f;
    std::string FileName=ePerm.eDir + "IOBP";
    eDrawArr.TitleStr=eRecVar.RecS.VarName2;
    eDrawArr.VarNameUF="IOBP";
    PLOT_PCOLOR(FileName, GrdArr, eDrawArr, eRecVar, eCall, ePerm);
  }
  int nbTrans=eBlPLOT.ListListDoubleValues.at("TransectLonStart").size();
  std::cerr << "GRID_PLOTTING nbTrans=" << nbTrans << "\n";
  if (nbTrans > 0) {
    std::vector<double> ListLonStart=eBlPLOT.ListListDoubleValues.at("TransectLonStart");
    std::vector<double> ListLonEnd  =eBlPLOT.ListListDoubleValues.at("TransectLonEnd");
    std::vector<double> ListLatStart=eBlPLOT.ListListDoubleValues.at("TransectLatStart");
    std::vector<double> ListLatEnd  =eBlPLOT.ListListDoubleValues.at("TransectLatEnd");
    double FrameLonLat=eBlPLOT.ListDoubleValues.at("FrameLonLat");
    //
    const MyMatrix<double> & DEP = GetDEP(GrdArr.GrdArrRho);
    PairMinMax ePair=ComputeMinMax(GrdArr, DEP);
    RecVar eRecVar;
    eRecVar.RecS.strAll="bathymetry2";
    eRecVar.RecS.VarName1="Bathymetry";
    eRecVar.RecS.VarName2="Bathymetry of model";
    eRecVar.RecS.minval=ePair.TheMin;
    eRecVar.RecS.maxval=ePair.TheMax;
    eRecVar.RecS.Unit="m";
    //
    for (int iTrans=0; iTrans<nbTrans; iTrans++) {
      double eLonStart=ListLonStart[iTrans];
      double eLonEnd  =ListLonEnd[iTrans];
      double eLatStart=ListLatStart[iTrans];
      double eLatEnd  =ListLatEnd[iTrans];
      std::vector<PairLL> ListCoord{{eLonStart, eLatStart}, {eLonEnd, eLatEnd}};
      double MinLon=std::min(eLonStart, eLonEnd) - FrameLonLat;
      double MaxLon=std::max(eLonStart, eLonEnd) + FrameLonLat;
      double MinLat=std::min(eLatStart, eLatEnd) - FrameLonLat;
      double MaxLat=std::max(eLatStart, eLatEnd) + FrameLonLat;
      QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};

      double distLON=GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat, eQuad.MaxLon, eQuad.MinLat);
      double distLAT=GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat, eQuad.MinLon, eQuad.MaxLat);
      double avgDistKM = 0.5;
      //      std::cerr << "avgDistKM=" << avgDistKM << " eMult=" << eMult << "\n";
      double nbLON=int(distLON / avgDistKM);
      double nbLAT=int(distLAT / avgDistKM);
      std::cerr << "nbLON=" << nbLON << " nbLAT=" << nbLAT << "\n";
      //
      // Now computing the GridArray necessary.
      //
      GridArray GrdArrOut=RECTANGULAR_GRID_ARRAY(eQuad, nbLON, nbLAT);
      SingleArrayInterpolation eInterp=GetSingleArrayInterpolationTrivialCase(GrdArrOut, GrdArr);
      eRecVar.F=SingleInterpolationOfField_2D(eInterp, GetDEP(GrdArr.GrdArrRho));
      std::cerr << "F(min/max)=" << eRecVar.F.minCoeff() << " / " << eRecVar.F.maxCoeff() << "\n";
      MyMatrix<uint8_t> MSK=ComputeInsideMask(eInterp);
      GrdArrOut.GrdArrRho.MSK=MSK;
      //
      SeqLineSegment eSeq{ListCoord, false};
      std::vector<SeqLineSegment> TheList{eSeq};
      //
      std::string eStrNumber=StringNumber(iTrans, 2);
      std::string TitleStr="Transect nr " + eStrNumber;
      std::string FileName=ePerm.eDir + "Transect_position_" + eStrNumber;
      double deltaLL=0.1;
      eQuad.MinLon -= deltaLL;
      eQuad.MinLat -= deltaLL;
      eQuad.MaxLon += deltaLL;
      eQuad.MaxLat += deltaLL;
      DrawArr eDrw=eDrawArr;
      eDrw.ListLineSegment=TheList;
      eDrw.eQuadFrame=eQuad;
      eDrw.DoTitle=false;
      eDrw.TitleStr=TitleStr;
      eDrw.VarNameUF="Transect_" + eStrNumber;
      //
      PLOT_PCOLOR(FileName, GrdArrOut, eDrw, eRecVar, eCall, ePerm);
    }
  }
}


struct PairDiscrete {
  MyMatrix<int> MatPoint;
  GraphSparseImmutable GR;
};

PairDiscrete ComputePairDiscrete(GridArray const& GrdArr)
{
  std::cerr << "ComputePairDiscrete IsFE=" << GrdArr.IsFE << "\n";
  if (GrdArr.IsFE == 1) {
    std::cerr << "Unstructured case\n";
    int nbNode=GrdArr.GrdArrRho.LON.size();
    std::cerr << "nbNode=" << nbNode << "\n";
    MyMatrix<int> MatPoint(nbNode,2);
    for (int iNode=0; iNode<nbNode; iNode++) {
      MatPoint(iNode,0) = iNode;
      MatPoint(iNode,1) = 0;
    }
    GraphSparseImmutable GR=GetUnstructuredVertexAdjInfo(GrdArr.INE, nbNode);
    return {MatPoint, GR};
  }
  //
  std::cerr << "Structured case\n";
  int nbRow=GrdArr.GrdArrRho.MSK.rows();
  int nbCol=GrdArr.GrdArrRho.MSK.cols();
  size_t nbWet=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      if (GrdArr.GrdArrRho.MSK(iRow,iCol) == 1)
        nbWet++;
  std::cerr << "nbWet=" << nbWet << "\n";
  MyMatrix<int> MatPoint(nbWet,2);
  int idx=0;
  MyMatrix<int> MatIdx(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      if (GrdArr.GrdArrRho.MSK(iRow,iCol) == 1) {
        MatPoint(idx,0) = iRow;
        MatPoint(idx,1) = iCol;
        MatIdx(iRow, iCol) = idx;
        idx++;
      }
    }
  std::cerr << "MatPoint / MatIdx built\n";
  std::vector<int> LVal = {1,0 , -1,0 , 0,1 , 0,-1};
  std::vector<int> ListNbAdj(nbWet,0);
  for (size_t iWet=0; iWet<nbWet; iWet++) {
    int iRow = MatPoint(iWet, 0);
    int iCol = MatPoint(iWet, 1);
    int nbAdj=0;
    for (int iAdj=0; iAdj<4; iAdj++) {
      int iRowAdj=iRow + LVal[2*iAdj];
      int iColAdj=iCol + LVal[2*iAdj+1];
      if (iRowAdj >= 0 && iRowAdj < nbRow && iColAdj >= 0 && iColAdj < nbCol) {
        if (GrdArr.GrdArrRho.MSK(iRowAdj,iColAdj) == 1)
          nbAdj++;
      }
    }
    ListNbAdj[iWet] = nbAdj;
  }
  std::cerr << "ListNbAdj built\n";
  std::vector<size_t> ListStart(nbWet+1,0);
  for (size_t iWet=0; iWet<nbWet; iWet++)
    ListStart[iWet+1] = ListStart[iWet] + ListNbAdj[iWet];
  int TotalSum = ListStart[nbWet];
  std::vector<size_t> ListListAdj(TotalSum);
  std::cerr << "Other assignation\n";
  idx=0;
  for (size_t iWet=0; iWet<nbWet; iWet++) {
    int iRow = MatPoint(iWet, 0);
    int iCol = MatPoint(iWet, 1);
    for (size_t iAdj=0; iAdj<4; iAdj++) {
      int iRowAdj=iRow + LVal[2*iAdj];
      int iColAdj=iCol + LVal[2*iAdj+1];
      if (iRowAdj >= 0 && iRowAdj < nbRow && iColAdj >= 0 && iColAdj < nbCol) {
        if (GrdArr.GrdArrRho.MSK(iRowAdj,iColAdj) == 1) {
          int iWetAdj = MatIdx(iRowAdj, iColAdj);
          ListListAdj[idx] = iWetAdj;
          idx++;
        }
      }
    }
  }
  std::cerr << "ListListAdj built\n";
  GraphSparseImmutable GR(nbWet, ListStart, ListListAdj);
  return {MatPoint, GR};
}



double ComputeAvgDist(PairDiscrete const& ePair, QuadArray const& eQuad, GridArray const& GrdArr)
{
  int nbNode = ePair.MatPoint.rows();
  std::vector<int> ListNodeStatus(nbNode,0);
  for (int iNode=0; iNode<nbNode; iNode++) {
    int iRow = ePair.MatPoint(iNode,0);
    int iCol = ePair.MatPoint(iNode,1);
    double eLon=GrdArr.GrdArrRho.LON(iRow,iCol);
    double eLat=GrdArr.GrdArrRho.LAT(iRow,iCol);
    int eStatus=0;
    if (eLon >= eQuad.MinLon && eLon <= eQuad.MaxLon && eLat >= eQuad.MinLat && eLat <= eQuad.MaxLat)
      eStatus=1;
    ListNodeStatus[iNode]=eStatus;
  }
  double sumDist=0;
  int nbDist=0;
  bool IsSpherical=GrdArr.IsSpherical;
  for (int iNode=0; iNode<nbNode; iNode++)
    if (ListNodeStatus[iNode] == 1) {
      int iRow = ePair.MatPoint(iNode,0);
      int iCol = ePair.MatPoint(iNode,1);
      double eLon1=GrdArr.GrdArrRho.LON(iRow,iCol);
      double eLat1=GrdArr.GrdArrRho.LAT(iRow,iCol);
      std::vector<size_t> ListAdj=ePair.GR.Adjacency(iNode);
      for (size_t & jNode : ListAdj) {
	  if (ListNodeStatus[jNode] == 1) {
            int jRow = ePair.MatPoint(jNode,0);
            int jCol = ePair.MatPoint(jNode,1);
	    double eLon2=GrdArr.GrdArrRho.LON(jRow, jCol);
	    double eLat2=GrdArr.GrdArrRho.LAT(jRow, jCol);
	    double eDist=GeodesicDistanceM_General(eLon1, eLat1, eLon2, eLat2, IsSpherical);
	    sumDist += eDist;
	    nbDist++;
	  }
	}
      }
  if (nbDist == 0) {
    std::cerr << "We cannot work out the subgrid:\n";
    std::cerr << "sumDist=" << sumDist << " nbDist=" << nbDist << "\n";
    std::cerr << "eQuad, lon(min/max)=" << eQuad.MinLon << " / " << eQuad.MaxLon << "\n";
    std::cerr << "eQuad, lat(min/max)=" << eQuad.MinLat << " / " << eQuad.MaxLat << "\n";
    std::cerr << "GrdArr.GrdArrRho.LON(min/max)=" << GrdArr.GrdArrRho.LON.minCoeff() << " / " << GrdArr.GrdArrRho.LON.maxCoeff() << "\n";
    std::cerr << "GrdArr.GrdArrRho.LAT(min/max)=" << GrdArr.GrdArrRho.LAT.minCoeff() << " / " << GrdArr.GrdArrRho.LAT.maxCoeff() << "\n";
    throw TerminalException{1};
  }
  double avgDist=sumDist / double(nbDist);
  return avgDist;
}



std::vector<InterpolToUVpoints> ComputeSpecificGrdArrInterpol(GridArray const& GrdArr, std::vector<QuadDrawInfo> const& ListQuadInfo, double const& eMult)
{
  std::cerr << "IsFE=" << GrdArr.IsFE << "\n";
  PairDiscrete ePair = ComputePairDiscrete(GrdArr);
  int nbQuad=ListQuadInfo.size();
  std::vector<InterpolToUVpoints> ListInterpol(nbQuad);
  bool IsSpherical=GrdArr.IsSpherical;
  for (int iQuad=0; iQuad<nbQuad; iQuad++) {
    QuadArray eQuad=ListQuadInfo[iQuad].eQuad;
    double avgDist = ComputeAvgDist(ePair, eQuad, GrdArr);
    double midLon=(eQuad.MinLon + eQuad.MaxLon)/double(2);
    double midLat=(eQuad.MinLat + eQuad.MaxLat)/double(2);
    double distLON=GeodesicDistanceM_General(eQuad.MinLon, midLat, eQuad.MaxLon, midLat, IsSpherical);
    double distLAT=GeodesicDistanceM_General(midLon, eQuad.MinLat, midLon, eQuad.MaxLat, IsSpherical);
    std::cerr << "distLON=" << distLON << " distLAT=" << distLAT << "\n";
    std::cerr << "avgDist=" << avgDist << " eMult=" << eMult << "\n";
    double nbLON=int(distLON / (avgDist*eMult));
    double nbLAT=int(distLAT / (avgDist*eMult));
    std::cerr << "nbLON=" << nbLON << " nbLAT=" << nbLAT << "\n";
    //
    // Now computing the GridArray necessary.
    //
    GridArray GrdArrOut=RECTANGULAR_GRID_ARRAY(eQuad, nbLON, nbLAT);
    GrdArrOut.IsSpherical=IsSpherical;
    SingleArrayInterpolation eInterp=GetSingleArrayInterpolationTrivialCase(GrdArrOut, GrdArr);
    MyMatrix<uint8_t> MSK=ComputeInsideMask(eInterp);
    GrdArrOut.GrdArrRho.MSK=MSK;
    //
    ListInterpol[iQuad] = {std::move(GrdArrOut), std::move(eInterp)};
  }
  return ListInterpol;
}






//
// We compute here several arrays that are needed for more advanced plots
//
void Compute_Additional_array(PermanentInfoDrawing & ePerm, const TotalArrGetData& TotalArr)
{
  SingleBlock eBlPLOT=ePerm.eFull.ListBlock.at("PLOT");
  SingleBlock eBlVAR=ePerm.eFull.ListBlock.at("VARS");
  std::vector<std::string> ListVar=ExtractMatchingBool(eBlVAR);
  auto NeedFDarrayPlot=[&](std::string const& eVar, std::vector<std::string> const&ListNatSought) -> bool {
    RecVar eRec=RetrieveTrivialRecVar(eVar);
    if (PositionVect(ListNatSought, eRec.RecS.VarNature) != -1)
      return true;
    return false;
  };
  auto ComputeNeedVar=[&](std::vector<std::string> const& ListNatSought) -> bool {
    for (auto & eVar : ListVar)
      if (NeedFDarrayPlot(eVar, ListNatSought))
	return true;
    return false;
  };
  //
  // The drawing array assignation
  //
  ePerm.eDrawArr=CommonAssignation_DrawArr(ePerm.eFull);
  //
  // The list of quad arrays
  //
  ePerm.ListQuadInfo=GetListQuadArray(eBlPLOT, TotalArr.GrdArr);
  //
  // The interpolation to a finite difference grid
  //
  int IsFE = TotalArr.GrdArr.IsFE;
  bool UseRegridArray=eBlPLOT.ListBoolValues.at("UseRegridArray");
  //  std::cerr << "|ListVar|=" << ListVar.size() << "\n";
  //  std::cerr << "UseRegridArray=" << UseRegridArray << "\n";
  bool NeedFDarray;
  if (IsFE == 0 && !UseRegridArray) {
    std::vector<std::string> ListNatUV={"uv", "3Duv"};
    NeedFDarray=ComputeNeedVar(ListNatUV);
  } else {
    NeedFDarray=true;
  }
  std::cerr << "NeedFDarray=" << NeedFDarray << "\n";
  if (NeedFDarray) {
    double eMult=ePerm.eFull.ListBlock.at("PLOT").ListDoubleValues.at("MultiplierResolutionRegrid");
    std::cerr << "eMult=" << eMult << "\n";
    ePerm.ListInterpol = ComputeSpecificGrdArrInterpol(TotalArr.GrdArr, ePerm.ListQuadInfo, eMult);
  }
  //
  // The computation of the arrays for transects
  //
  std::vector<std::string> ListNat3D={"3Drho", "3Duv"};
  bool NeedTransArray=ComputeNeedVar(ListNat3D);
  int nbTrans=0;
  if (NeedTransArray) {
    std::vector<double> ListLonStart=eBlPLOT.ListListDoubleValues.at("TransectLonStart");
    std::vector<double> ListLatStart=eBlPLOT.ListListDoubleValues.at("TransectLatStart");
    std::vector<double> ListLonEnd=eBlPLOT.ListListDoubleValues.at("TransectLonEnd");
    std::vector<double> ListLatEnd=eBlPLOT.ListListDoubleValues.at("TransectLatEnd");
    std::vector<double> ListResolKM=eBlPLOT.ListListDoubleValues.at("TransectSpatialResolutionKM");
    std::vector<double> ListVertResolM=eBlPLOT.ListListDoubleValues.at("TransectVerticalResolutionM");
    nbTrans=ListResolKM.size();
    size_t nbTrans_t=ListResolKM.size();
    if (nbTrans_t != ListLonStart.size() || nbTrans_t != ListLatStart.size() || nbTrans_t != ListLonEnd.size() || nbTrans_t != ListLonEnd.size() || nbTrans_t != ListVertResolM.size()) {
      std::cerr << "Error in the transect information input\n";
      std::cerr << "The number should all be the same for:\n";
      std::cerr << "|TransectLonStart|=" << ListLonStart.size() << "\n";
      std::cerr << "|TransectLatStart|=" << ListLatStart.size() << "\n";
      std::cerr << "  |TransectLonEnd|=" << ListLonEnd.size() << "\n";
      std::cerr << "  |TransectLatEnd|=" << ListLatEnd.size() << "\n";
      std::cerr << "|TransectSpatialResolutionKM|=" << ListResolKM.size() << "\n";
      std::cerr << "|TransectVerticalResolutionM|=" << ListVertResolM.size() << "\n";
      std::cerr << "TransectLonStart, TransectLatStart, TransectLonEnd, TransectLatEnd, TransectSpatialResolutionKM and TransectVerticalResolutionM\n";
      throw TerminalException{1};
    }
    if (nbTrans > 0) {
      std::vector<TransectInformation_3D> ListTransect(nbTrans);
      Eigen::Tensor<double,3> VertCoord=RetrieveStandardVerticalCoordinate(TotalArr);
      for (int iTrans=0; iTrans<nbTrans; iTrans++) {
	double eLonStart=ListLonStart[iTrans];
	double eLonEnd  =ListLonEnd[iTrans];
	double eLatStart=ListLatStart[iTrans];
	double eLatEnd  =ListLatEnd[iTrans];
	double eResolKM=ListResolKM[iTrans];
	double eResolM=ListVertResolM[iTrans];
	TransectInformation eTrans=GetTransectInformation({TotalArr.GrdArr},
							  eLonStart, eLatStart,
							  eLonEnd, eLatEnd,
							  eResolKM);
	ListTransect[iTrans]=GetTransectInformation_3D(eTrans, TotalArr.GrdArr, VertCoord, eResolM);
      }
      ePerm.ListTransect = std::move(ListTransect);
    }
  }
}



TripleModelDesc Retrieve_triple_from_array(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  std::string BoundFile=eBlPROC.ListStringValues.at("BoundFile");
  std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
  std::string Sphericity=eBlPROC.ListStringValues.at("Sphericity");
  bool CutWorldMap=eBlPROC.ListBoolValues.at("CutWorldMap");
  //  std::cerr << "1: CutWorldMap=" << CutWorldMap << "\n";
  bool HigherLatitudeCut=eBlPROC.ListBoolValues.at("HigherLatitudeCut");
  double MinLatCut=eBlPROC.ListDoubleValues.at("MinLatCut");
  double MaxLatCut=eBlPROC.ListDoubleValues.at("MaxLatCut");
  GridSymbolic RecGridSymb(Sphericity, CutWorldMap, HigherLatitudeCut, MinLatCut, MaxLatCut, 0, 0, 0, 0, 0);
  //  std::cerr << "2: CutWorldMap=" << RecGridSymb.CutWorldMap << "\n";
  TripleModelDesc eTriple{eModelName, GridFile, BoundFile, HisPrefix, RecGridSymb};
  return eTriple;
}






void SINGLE_Plotting_Function(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  bool WriteITimeInFileName=eBlPROC.ListBoolValues.at("WriteITimeInFileName");
  //
  // Creating the triple for data reading (grid and history)
  //
  //  std::cerr << "SINGLE_Plotting_Function, step 0\n";
  //  std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
  TripleModelDesc eTriple = Retrieve_triple_from_array(eFull);
  //
  // Retrieving the grid array
  //
  GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
  std::cerr << "Sphericity=" << GrdArr.IsSpherical << "\n";
  //
  // Setting up the timings.
  //
  ArrayHistory eArr=ReadArrayHistory(eTriple);
  std::vector<VarQuery> ListQuery=GetIntervalGen_Query(eBlPROC, {eArr});
  int nbTime=ListQuery.size();
  int MaxNbTime=eBlPROC.ListIntValues.at("MaxNbTime");
  if (MaxNbTime != -1)
    nbTime=MaxNbTime;
  std::cerr << "nbTime=" << nbTime << "\n";
  //
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  //
  // Setting up DrawArr
  // It contains the permanent feature of plots that will not change from variable
  // to variable and time to time
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // It has to be put there so that it is destroyed before ePerm.PrefixTemp
  Compute_Additional_array(ePerm, TotalArr);
  //
  // Preliminary drawings
  //
  GRID_PLOTTING(GrdArr, eTriple.GridFile, eCall, ePerm);
  //
  // Now the major time loop
  //
  SingleBlock eBlockVAR=eFull.ListBlock.at("VARS");
  std::vector<std::string> ListVarOut=ExtractMatchingBool(eBlockVAR);
  std::vector<std::string> ListNatureQuery=eBlPROC.ListListStringValues.at("ListNatureQuery");
  PlotBound ePlotBound=ReadPlotBound(eFull);
  std::cerr << "|ListVarOut|=" << ListVarOut.size() << "\n";
  for (auto &eVar : ListVarOut)
    std::cerr << "eVar=" << eVar << "\n";
  for (int iTime=0; iTime<nbTime; iTime++) {
    VarQuery eQuery=ListQuery[iTime];
    if (!WriteITimeInFileName)
      eQuery.iTime=-1;
    std::string strPres=DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << " date=" << strPres << "\n";
    for (auto & eVarName : ListVarOut)
      for (auto & eNatureQuery : ListNatureQuery) {
	eQuery.NatureQuery=eNatureQuery;
	RecVar eRecVar=ModelSpecificVarSpecificTimeGeneral(TotalArr, eVarName, eQuery, ePlotBound);
	GENERAL_PLOT_SINGLE(GrdArr, eRecVar, eCall, ePerm);
      }
  }
}

void GENERAL_PLOT_PAIR(GridArray const& GrdArr,
		       RecVar const& eRecVar1,
		       RecVar const& eRecVar2,
		       NCLcaller<GeneralType> & eCall,
		       PermanentInfoDrawing const& ePerm)
{
  if (eRecVar1.RecS.VarNature == "uv") {
    std::cerr << "For the plotting of pair and the field of the type UV\n";
    std::cerr << "we do not yet have relevant code\n";
    std::cerr << "VarName1=" << eRecVar1.RecS.VarName1 << "\n";
    return;
  }
  if (eRecVar1.RecS.VarNature == "rho") {
    PLOT_DIFF_FD_RHO_PCOLOR(GrdArr, eRecVar1, eRecVar2, eCall, ePerm);
    return;
  }
  std::cerr << "No matching VarNature\n";
  std::cerr << "eRecVar1.RecS.VarNature = " << eRecVar1.RecS.VarNature << "\n";
  throw TerminalException{1};
}





void PAIR_Plotting_Function(FullNamelist const& eFull)
{
  //
  // Retrieving the grid array and other arrays
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  std::string BoundFile=eBlPROC.ListStringValues.at("BoundFile");
  std::string HisPrefix1=eBlPROC.ListStringValues.at("HisPrefix1");
  std::string HisPrefix2=eBlPROC.ListStringValues.at("HisPrefix2");
  bool WriteITimeInFileName=eBlPROC.ListBoolValues.at("WriteITimeInFileName");
  std::string Sphericity=eBlPROC.ListStringValues.at("Sphericity");
  bool CutWorldMap=eBlPROC.ListBoolValues.at("CutWorldMap");
  bool HigherLatitudeCut=eBlPROC.ListBoolValues.at("HigherLatitudeCut");
  double MinLatCut=eBlPROC.ListDoubleValues.at("MinLatCut");
  double MaxLatCut=eBlPROC.ListDoubleValues.at("MaxLatCut");
  GridSymbolic RecGridSymb(Sphericity, CutWorldMap, HigherLatitudeCut, MinLatCut, MaxLatCut, 0, 0, 0, 0, 0);
  TripleModelDesc eTriple1{eModelName, GridFile, BoundFile, HisPrefix1, RecGridSymb};
  TripleModelDesc eTriple2{eModelName, GridFile, BoundFile, HisPrefix2, RecGridSymb};
  GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple1); // Right now we have same grid for both
  ArrayHistory eArr1=ReadArrayHistory(eTriple1);
  ArrayHistory eArr2=ReadArrayHistory(eTriple2);
  TotalArrGetData TotalArr1 = RetrieveTotalArr(eTriple1);
  TotalArrGetData TotalArr2 = RetrieveTotalArr(eTriple2);
  std::cerr << "We have eArr. Printout of eArr\n";
  std::vector<VarQuery> ListQuery=GetIntervalGen_Query(eBlPROC, {eArr1, eArr2});
  int nbTime=ListQuery.size();
  int MaxNbTime=eBlPROC.ListIntValues.at("MaxNbTime");
  if (MaxNbTime != -1)
    nbTime=MaxNbTime;
  std::cerr << "nbTime=" << nbTime << "\n";
  SingleBlock eBlockVAR=eFull.ListBlock.at("VARS");
  //
  // Setting up DrawArr
  // It contains the permanent feature of plots that will not change from variable
  // to variable and time to time
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // It has to be put there so that it is destroyed before
                                      // ePerm.PrefixTemp
  Compute_Additional_array(ePerm, TotalArr1);
  //
  // Now the major time loop
  //
  std::vector<std::string> ListVarOut=ExtractMatchingBool(eBlockVAR);
  std::vector<std::string> ListNatureQuery=eBlPROC.ListListStringValues.at("ListNatureQuery");
  PlotBound ePlotBound=ReadPlotBound(eFull);
  for (int iTime=0; iTime<nbTime; iTime++) {
    VarQuery eQuery=ListQuery[iTime];
    if (!WriteITimeInFileName)
      eQuery.iTime=-1;
    std::string strPres=DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << " date=" << strPres << "\n";
    for (auto & eVarName : ListVarOut)
      for (auto & eNatureQuery : ListNatureQuery) {
	eQuery.NatureQuery=eNatureQuery;
	PairRecVar ePairRecVar=ModelPairSpecificVarSpecificTimeGeneral(TotalArr1, TotalArr2, eVarName, eQuery, ePlotBound);
	GENERAL_PLOT_PAIR(GrdArr, ePairRecVar.RecVar1, ePairRecVar.RecVar2, eCall, ePerm);
      }
  }
}




#endif
