#ifndef FLOATS_INCLUDE
#define FLOATS_INCLUDE

#include "Plotting_fct.h"
#include "Basic_netcdf.h"
#include "Basic_string.h"
#include "Basic_plot.h"
#include "NamelistExampleOcean.h"
#include "SphericalGeom.h"
#include "Model_grids.h"





void PLOT_ROMS_float(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  std::string BoundFile=eBlPROC.ListStringValues.at("BoundFile");
  std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
  std::string PicPrefix=eBlPROC.ListStringValues.at("PicPrefix");
  std::string FloatFile=eBlPROC.ListStringValues.at("FloatFile");
  std::string FileDescFloat=eBlPROC.ListStringValues.at("FileDescFloat");
  //
  std::string strBEGTC=eBlPROC.ListStringValues.at("BEGTC");
  std::string strENDTC=eBlPROC.ListStringValues.at("ENDTC");
  //
  std::string Sphericity=eBlPROC.ListStringValues.at("Sphericity");
  bool CutWorldMap=eBlPROC.ListBoolValues.at("CutWorldMap");
  bool HigherLatitudeCut=eBlPROC.ListBoolValues.at("HigherLatitudeCut");
  double MinLatCut=eBlPROC.ListDoubleValues.at("MinLatCut");
  double MaxLatCut=eBlPROC.ListDoubleValues.at("MaxLatCut");
  GridSymbolic RecGridSymb(Sphericity, CutWorldMap, HigherLatitudeCut, MinLatCut, MaxLatCut, 0, 0, 0, 0, 0);
  TripleModelDesc eTriple{eModelName, GridFile, BoundFile, HisPrefix, RecGridSymb};
  //
  ArrayHistory eArr=ReadArrayHistory(eTriple);
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  std::vector<QuadDrawInfo> ListQuad = GetListQuadArray(eBlPLOT, TotalArr.GrdArr);
  bool PlotDensity = eBlPLOT.ListBoolValues.at("PlotDensity");
  bool PlotTrajectory = eBlPLOT.ListBoolValues.at("PlotTrajectory");
  std::cerr << "PlotDensity=" << PlotDensity << " PlotTrajectory=" << PlotTrajectory << "\n";
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm
  //
  std::vector<double> LTime=NC_ReadTimeFromFile(FloatFile, "ocean_time");
  int nbTime = LTime.size();
  int idx_first = 0;
  int idx_last  = nbTime;
  if (strBEGTC != "earliest") {
    double BeginTime=CT2MJD(strBEGTC);
    for (int idx=0; idx<nbTime; idx++) {
      if (LTime[idx] <= BeginTime)
        idx_first = idx;
    }
  }
  if (strENDTC != "latest") {
    double EndTime  =CT2MJD(strENDTC);
    for (int idx=0; idx<nbTime; idx++) {
      int jdx = nbTime - 1 - idx;
      if (LTime[jdx] >= EndTime)
        idx_last = jdx;
    }
  }
  std::cerr << "idx_first=" << idx_first << " idx_last=" << idx_last << " nbTime=" << nbTime << "\n";
  std::cerr << "Begin=" << DATE_ConvertMjd2mystringPres(LTime[idx_first]) << " End=" << DATE_ConvertMjd2mystringPres(LTime[idx_last-1]) << "\n";
  RecVar eRecVar = ModelSpecificVarSpecificTime_Kernel(TotalArr, "Bathymetry", LTime[0]);
  int idx_len = idx_last - idx_first;
  //
  netCDF::NcFile dataFile(FloatFile, netCDF::NcFile::read);
  netCDF::NcVar lon_var = dataFile.getVar("lon");
  netCDF::NcVar lat_var = dataFile.getVar("lat");
  netCDF::NcVar depth_var = dataFile.getVar("depth");
  netCDF::NcVar temp_var = dataFile.getVar("temp");
  netCDF::NcVar salt_var = dataFile.getVar("salt");
  netCDF::NcDim dimDrifter = lon_var.getDim(1);
  MyMatrix<double> LON_mat = NC_Read2Dvariable_data(lon_var);
  MyMatrix<double> LAT_mat = NC_Read2Dvariable_data(lat_var);
  int nb_drifter = dimDrifter.getSize();
  auto const& GrdAR = TotalArr.GrdArr.GrdArrRho;
  double LONmax = GrdAR.LON.maxCoeff();
  double LONmin = GrdAR.LON.minCoeff();
  double LATmax = GrdAR.LAT.maxCoeff();
  double LATmin = GrdAR.LAT.minCoeff();

  // Assigning the values
  std::vector<std::string> ListFloatDesc;
  if (FileDescFloat == "unset") {
    for (int i=0; i<nb_drifter; i++)
      ListFloatDesc.push_back(std::to_string(i+1));
  } else {
    ListFloatDesc = ReadFullFile(FileDescFloat);
  }
  if (nb_drifter != int(ListFloatDesc.size())) {
    std::cerr << "The number of entries ni ListFloatDesc does not match the number of drifters\n";
    std::cerr << "nb_drifter=" << nb_drifter << " |ListFloatDesc|=" << ListFloatDesc.size() << "\n";
    std::cerr << "FileDescFloat = " << FileDescFloat << "\n";
    throw TerminalException{1};
  }
  // plotting the drifters themselves
  size_t TotalNbPoint = 0;
  for (int i_drifter=0; i_drifter<nb_drifter; i_drifter++) {
    std::cerr << "i_drifter=" << i_drifter << "/" << nb_drifter << "\n";

    std::vector<PairLL> ListPairLL;
    double deltaLL = 1;
    for (int idx=0; idx<idx_len; idx++) {
      double eLon = LON_mat(idx + idx_first, i_drifter);
      double eLat = LAT_mat(idx + idx_first, i_drifter);
      if (eLon < LONmax + deltaLL && eLon > LONmin - deltaLL &&
          eLat < LATmax + deltaLL && eLat > LATmin - deltaLL) {
        PairLL eP{eLon, eLat};
        ListPairLL.push_back(eP);
      }
    }
    std::cerr << "idx_len=" << idx_len << " |ListPairLL|=" << ListPairLL.size() << "\n";
    TotalNbPoint += ListPairLL.size();
    if (PlotTrajectory) {
      SeqLineSegment eSeq = {ListPairLL, false};
      for (auto & eQuad : ListQuad) {
        std::cerr << "iFrame=" << eQuad.iFrame << " eFrameName=" << eQuad.eFrameName << "\n";
        std::string FileName = PicPrefix + "Drifter" + StringNumber(i_drifter+1, 4) + "_" + eQuad.eFrameName;
        std::cerr << "FileName = " << FileName << "\n";
        DrawArr eDrw = BasicArrayDraw(eQuad.eQuad);
        eDrw.DoTitle = true;
        eDrw.TitleStr = "Drifter " + ListFloatDesc[i_drifter];
        eRecVar.RecS.strAll = std::to_string(i_drifter) + "_" + eQuad.eFrameName;
        eDrw.ListLineSegment = {eSeq};
        std::cerr << "Before PLOT_PCOLOR 1\n";
        PLOT_PCOLOR(FileName, TotalArr.GrdArr, eDrw, eRecVar, eCall, ePerm);
        std::cerr << "After PLOT_PCOLOR 1\n";
      }
    }
  }
  //
  // Now plotting the data
  //
  std::cerr << "TotalNbPoint=" << TotalNbPoint << "\n";
  if (PlotDensity) {
    MyMatrix<double> ListXY(2,TotalNbPoint);
    size_t pos=0;
    for (int i_drifter=0; i_drifter<nb_drifter; i_drifter++) {
      double deltaLL = 1;
      for (int idx=0; idx<idx_len; idx++) {
        double eLon = LON_mat(idx + idx_first, i_drifter);
        double eLat = LAT_mat(idx + idx_first, i_drifter);
        if (eLon < LONmax + deltaLL && eLon > LONmin - deltaLL &&
            eLat < LATmax + deltaLL && eLat > LATmin - deltaLL) {
          ListXY(0,pos) = eLon;
          ListXY(1,pos) = eLat;
          pos++;
        }
      }
    }
    std::vector<SingleRecInterp> LRec = General_FindInterpolationWeight(TotalArr.GrdArr, ListXY);
    int eta_rho = GrdAR.LON.rows();
    int xi_rho = GrdAR.LON.cols();
    MyMatrix<double> DensMat = ZeroMatrix<double>(eta_rho, xi_rho);
    for (auto & eRec : LRec) {
      if (eRec.status) {
        for (auto & ePart : eRec.LPart)
          DensMat(ePart.eEta, ePart.eXi) += ePart.eCoeff;
      }
    }
    double eMax = DensMat.maxCoeff();
    DensMat /= eMax;
    eRecVar.F = DensMat;
    eRecVar.RecS.minval = 0;
    eRecVar.RecS.maxval = 1;
    std::cerr << "eMax = " << eMax << "\n";
    for (auto & eQuad : ListQuad) {
      std::string FileName = PicPrefix + "Density_Plot_"  + eQuad.eFrameName;
      eRecVar.RecS.strAll = "Density_plot_" + eQuad.eFrameName;
      DrawArr eDrw = ePerm.eDrawArr;
      eDrw.eQuadFrame = eQuad.eQuad;
      eDrw.DoTitle = true;
      eDrw.TitleStr = "Density plot";
      std::cerr << "Before PLOT_PCOLOR 2\n";
      PLOT_PCOLOR(FileName, TotalArr.GrdArr, eDrw, eRecVar, eCall, ePerm);
      std::cerr << "After PLOT_PCOLOR 2\n";
    }
  }

}





void ROMS_CreateLTransFile(std::string const& eFile, GridArray const& GrdArr, double const& DistanceKM, std::vector<double> const& ListLONpt, std::vector<double> const& ListLATpt, std::vector<double> const& ListDepth, std::vector<int> const& ListTime)
{
  int eta_psi=GrdArr.GrdArrPsi.LON.rows();
  int xi_psi=GrdArr.GrdArrPsi.LON.cols();
  std::vector<double> ListLon;
  std::vector<double> ListLat;
  int nbPt=ListLONpt.size();
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++) {
      double eLon=GrdArr.GrdArrPsi.LON(i,j);
      double eLat=GrdArr.GrdArrPsi.LAT(i,j);
      double MinDistKM=10000;
      for (int iPt=0; iPt<nbPt; iPt++) {
	double eDistKM=GeodesicDistanceKM(eLon, eLat, ListLONpt[iPt], ListLATpt[iPt]);
	if (eDistKM < MinDistKM)
	  MinDistKM=eDistKM;
      }
      if (MinDistKM < DistanceKM) {
	ListLon.push_back(eLon);
	ListLat.push_back(eLat);
      }
    }
  //
  std::ofstream OUT(eFile);
  int nbPart=ListLon.size();
  int nbDep=ListDepth.size();
  int nbTime=ListTime.size();
  int numpar=0;
  for (int iPart=0; iPart<nbPart; iPart++)
    for (int iDep=0; iDep<nbDep; iDep++)
      for (int iTime=0; iTime<nbTime; iTime++) {
	OUT << ListLon[iPart] << "," << ListLat[iPart] << "," << ListDepth[iDep] << "," << ListTime[iTime] << ",101001\n";
	numpar++;
      }
  std::cerr << "numpar=" << numpar << "\n";
}




void ROMS_CreateRomsOfflineFile(std::string const& eFile, GridArray const& GrdArr, double const& DistanceKM, std::vector<double> const& ListLONpt, std::vector<double> const& ListLATpt, std::vector<double> const& ListDepth, std::vector<int> const& ListTime)
{
  int eta_psi=GrdArr.GrdArrPsi.LON.rows();
  int xi_psi=GrdArr.GrdArrPsi.LON.cols();
  std::vector<double> ListLon;
  std::vector<double> ListLat;
  int nbPt=ListLONpt.size();
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++) {
      double eLon=GrdArr.GrdArrPsi.LON(i,j);
      double eLat=GrdArr.GrdArrPsi.LAT(i,j);
      double MinDistKM=10000;
      for (int iPt=0; iPt<nbPt; iPt++) {
	double eDistKM=GeodesicDistanceKM(eLon, eLat, ListLONpt[iPt], ListLATpt[iPt]);
	if (eDistKM < MinDistKM)
	  MinDistKM=eDistKM;
      }
      if (MinDistKM < DistanceKM) {
	ListLon.push_back(eLon);
	ListLat.push_back(eLat);
      }
    }
  //
  std::ofstream OUT(eFile);
  int nbPart=ListLon.size();
  int nbDep=ListDepth.size();
  int nbTime=ListTime.size();
  int numpar=0;
  OUT << "1  Ftitle (a80)\n";
  OUT << "ROMS 1.0 - Initial Drifters Locations -\n";
  OUT << "2  Ft0,   Fx0,   Fy0,     Fz0, Fcoor,Ftype,Fcount, Fdt,   Fdx,   Fdy,   Fdz\n";
  for (int iPart=0; iPart<nbPart; iPart++)
    for (int iDep=0; iDep<nbDep; iDep++)
      for (int iTime=0; iTime<nbTime; iTime++) {
	OUT << ListTime[iTime] << "   " << ListLon[iPart] << "   " << ListLat[iPart] << "    " << ListDepth[iDep] << "    1    0    1      0.00   0.000   0.000   0.0\n";
	numpar++;
      }
  OUT << "99 END OF FILE\n";
  std::cerr << "numpar=" << numpar << "\n";
}

/*
std::vector<FloatState> CREATE_ListFloatStates(GridArray const& GrdArr, double const& DistanceKM, std::vector<double> const& ListLONpt, std::vector<double> const& ListLATpt, std::vector<double> const& ListDepth, std::vector<int> const& ListTime)
{
  int eta_psi=GrdArr.GrdArrPsi.eta;
  int xi_psi=GrdArr.GrdArrPsi.xi;
  std::vector<double> ListLon;
  std::vector<double> ListLat;
  int nbPt=ListLONpt.size();
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++) {
      double eLon=GrdArr.GrdArrPsi.LON(i,j);
      double eLat=GrdArr.GrdArrPsi.LAT(i,j);
      double MinDistKM=10000;
      for (int iPt=0; iPt<nbPt; iPt++) {
	double eDistKM=GeodesicDistanceKM(eLon, eLat, ListLONpt[iPt], ListLATpt[iPt]);
	if (eDistKM < MinDistKM)
	  MinDistKM=eDistKM;
      }
      if (MinDistKM < DistanceKM) {
	ListLon.push_back(eLon);
	ListLat.push_back(eLat);
      }
    }
  //
  std::ofstream OUT(eFile);
  int nbPart=ListLon.size();
  int nbDep=ListDepth.size();
  int nbTime=ListTime.size();
  int numpar=0;
  OUT << "1  Ftitle (a80)\n";
  OUT << "ROMS 1.0 - Initial Drifters Locations -\n";
  OUT << "2  Ft0,   Fx0,   Fy0,     Fz0, Fcoor,Ftype,Fcount, Fdt,   Fdx,   Fdy,   Fdz\n";
  for (int iPart=0; iPart<nbPart; iPart++)
    for (int iDep=0; iDep<nbDep; iDep++)
      for (int iTime=0; iTime<nbTime; iTime++) {
	OUT << ListTime[iTime] << "   " << ListLon[iPart] << "   " << ListLat[iPart] << "    " << ListDepth[iDep] << "    1    0    1      0.00   0.000   0.000   0.0\n";
	numpar++;
      }
  OUT << "99 END OF FILE\n";
  std::cerr << "numpar=" << numpar << "\n";
}
*/







std::vector<double> ICHTHYOP_ReadListTime(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ReadTimeFromFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile=" << eFile << "\n";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  std::string StringTime="time";
  netCDF::NcVar data=dataFile.getVar(StringTime);
  if (data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "StringTime=" << StringTime << "\n";
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
  netCDF::NcVarAtt eTimeAtt=data.getAtt("origin");
  char eString[1024]="";
  //  fprintf(stderr, "Before call to getValues\n");
  eTimeAtt.getValues(eString);
  std::string eStringStr=eString;
  std::vector<std::string> ListStr=STRING_Split(eStringStr, " ");
  std::string eStrYear=ListStr[1];
  std::string eStrMonth=ListStr[3];
  std::string eStrDay=ListStr[5];
  std::string eStrTime=ListStr[7];
  int year, month, day;
  std::istringstream(eStrYear) >> year;
  std::istringstream(eStrMonth) >> month;
  std::istringstream(eStrDay) >> day;
  //
  std::vector<std::string> eVectTime=STRING_Split(eStrTime,":");
  std::string eStrHour=eVectTime[0];
  std::string eStrMin=eVectTime[1];
  int hour, min, sec;
  std::istringstream(eStrHour) >> hour;
  std::istringstream(eStrMin) >> min;
  sec=0;
  double eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, sec});
  //
  std::vector<double> ListTime(siz);
  for (int iTime=0; iTime<siz; iTime++) {
    double eValSec=eVal[iTime];
    double eTime=eTimeStart + eValSec / double(86400);
    ListTime[iTime]=eTime;
  }
  return ListTime;
}


void ICHTHYOP_PlotTrajectories(FullNamelist const& eFull)
{
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm for destructor ordering
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  std::string DrifterFile=eBlPROC.ListStringValues.at("DrifterFile");
  std::vector<double> ListTime=ICHTHYOP_ReadListTime(DrifterFile);
  int nbTime=ListTime.size();
  double MinDistanceTrajKM=eBlPLOT.ListDoubleValues.at("MinimalDistanceTrajectoriesKM");
  double FrameLonLat=eBlPLOT.ListDoubleValues.at("FrameLonLat");
  double MinLon_frame=eBlPLOT.ListDoubleValues.at("MinLon");
  double MaxLon_frame=eBlPLOT.ListDoubleValues.at("MaxLon");
  double MinLat_frame=eBlPLOT.ListDoubleValues.at("MinLat");
  double MaxLat_frame=eBlPLOT.ListDoubleValues.at("MaxLat");
  bool FixedFrame=eBlPLOT.ListBoolValues.at("FixedFrame");
  bool IndividualTrajectories=eBlPLOT.ListBoolValues.at("IndividualTrajectories");
  bool DensityPassingPlot=eBlPLOT.ListBoolValues.at("DensityPassingPlot");
  bool DensityFinalPlot=eBlPLOT.ListBoolValues.at("DensityFinalPlot");
  //
  std::string VarLon="lon";
  std::string VarLat="lat";
  MyMatrix<double> ListLon=NC_Read2Dvariable(DrifterFile, VarLon);
  MyMatrix<double> ListLat=NC_Read2Dvariable(DrifterFile, VarLat);
  int nbDrifter=ListLon.cols();
  std::cerr << "nbDrifter=" << nbDrifter << "\n";
  if (IndividualTrajectories) {
    for (int iDrifter=0; iDrifter<nbDrifter; iDrifter++) {
      std::cerr << "iDrifter=" << iDrifter << "\n";
      std::vector<PairLL> ListCoord;
      double MinLon=0, MinLat=0, MaxLon=0, MaxLat=0;
      bool IsFirst=true;
      for (int iTime=0; iTime<nbTime; iTime++) {
	double eLon=ListLon(iTime, iDrifter);
	double eLat=ListLat(iTime, iDrifter);
	PairLL ePairLL{eLon, eLat};
	ListCoord.push_back(ePairLL);
	if (IsFirst) {
	  MinLon=eLon;
	  MinLat=eLat;
	  MaxLon=eLon;
	  MaxLat=eLat;
	  IsFirst=false;
	}
	else {
	  if (eLon < MinLon)
	    MinLon=eLon;
	  if (eLon > MaxLon)
	    MaxLon=eLon;
	  if (eLat < MinLat)
	    MinLat=eLat;
	  if (eLat > MaxLat)
	    MaxLat=eLat;
	}
      }
      MinLon -= FrameLonLat;
      MinLat -= FrameLonLat;
      MaxLon += FrameLonLat;
      MaxLat += FrameLonLat;
      if (FixedFrame) {
	MinLon=MinLon_frame;
	MaxLon=MaxLon_frame;
	MinLat=MinLat_frame;
	MaxLat=MaxLat_frame;
      }
      double TotalDistKM=0;
      for (int iTime=1; iTime<nbTime; iTime++) {
	double eLon1=ListLon(iTime-1, iDrifter);
	double eLat1=ListLat(iTime-1, iDrifter);
	double eLon2=ListLon(iTime, iDrifter);
	double eLat2=ListLat(iTime, iDrifter);
	double eDistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	TotalDistKM += eDistKM;
      }
      int nbSplitLon=100;
      int nbSplitLat=100;
      QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};
      GridArray GrdArr=RECTANGULAR_GRID_ARRAY(eQuad, nbSplitLon, nbSplitLat);
      RecVar eRecVar=GetTrivialArrayPlot(GrdArr);
      //
      SeqLineSegment eSeq{ListCoord, false};
      std::vector<SeqLineSegment> TheList{eSeq};
      //
      std::cerr << "Lon(max - min)=" << MaxLon - MinLon << "\n";
      std::cerr << "Lat(max - min)=" << MaxLat - MinLat << "\n";
      std::cerr << "TotalDistKM=" << TotalDistKM << "\n";
      if (TotalDistKM > MinDistanceTrajKM) {
	std::string eStrNumber=StringNumber(iDrifter, 6);
	std::string TitleStr="Drifter trajectory nr " + eStrNumber;
	std::string FileName=ePerm.eDir + "Drifter_traj_" + eStrNumber;
	DrawArr eDrw;
	eDrw.DoTitle=true;
	eDrw.TitleStr=TitleStr;
	eDrw.DoColorBar=false;
	eDrw.ColorMap="WhBlGrYeRe";
	eDrw.eQuadFrame=eQuad;
	eDrw.DoTitle=false;
	eDrw.ListLineSegment=TheList;
	//
	PLOT_PCOLOR(FileName, GrdArr, eDrw, eRecVar, eCall, ePerm);
      }
    }
  }
  double MinLonPl=MinLon_frame;
  double MaxLonPl=MaxLon_frame;
  double MinLatPl=MinLat_frame;
  double MaxLatPl=MaxLat_frame;
  QuadArray eQuad{MinLonPl, MaxLonPl, MinLatPl, MaxLatPl};
  double DeltaLonDensPlot=eBlPLOT.ListDoubleValues.at("DeltaLonDensPlot");
  double DeltaLatDensPlot=eBlPLOT.ListDoubleValues.at("DeltaLatDensPlot");
  int nbSplitLon=round((MaxLonPl - MinLonPl)/DeltaLonDensPlot);
  int nbSplitLat=round((MaxLatPl - MinLatPl)/DeltaLatDensPlot);
  GridArray GrdArr=RECTANGULAR_GRID_ARRAY(eQuad, nbSplitLon, nbSplitLat);
  double DeltaLonWork=(MaxLonPl - MinLonPl)/double(nbSplitLon);
  double DeltaLatWork=(MaxLatPl - MinLatPl)/double(nbSplitLat);
  if (DensityPassingPlot) {
    MyMatrix<double> TheDens=ZeroMatrix<double>(nbSplitLon, nbSplitLat);
    for (int iDrifter=0; iDrifter<nbDrifter; iDrifter++) {
      for (int iTime=1; iTime<nbTime; iTime++) {
	double eLon1=ListLon(iTime-1, iDrifter);
	double eLat1=ListLat(iTime-1, iDrifter);
	double eLon2=ListLon(iTime, iDrifter);
	double eLat2=ListLat(iTime, iDrifter);
	double eDistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	if (eDistKM > 0) {
	  int i=floor((eLon2 - MinLonPl)/DeltaLonWork);
	  int j=floor((eLat2 - MinLatPl)/DeltaLatWork);
	  if (i>= 0 && i<nbSplitLon && j>=0 && j<nbSplitLat) {
	    TheDens(i,j)++;
	  }
	}
      }
    }
    double pNorm=6;
    double maxVal=TheDens.maxCoeff();
    double eSum=0;
    for (int i=0; i<nbSplitLon; i++)
      for (int j=0; j<nbSplitLat; j++)
	eSum += pow(TheDens(i,j), pNorm);
    double eLpNorm=pow(eSum/double(nbSplitLon*nbSplitLat), 1/pNorm);
    std::cerr << "maxVal=" << maxVal << " eLpNorm=" << eLpNorm << "\n";
    std::string TitleStr="Density of passing points";
    std::string FileName=ePerm.eDir + "Density_passing_point";
    DrawArr eDrw;
    eDrw.DoTitle=true;
    eDrw.TitleStr=TitleStr;
    eDrw.DoColorBar=true;
    eDrw.eQuadFrame=eQuad;
    eDrw.DoTitle=false;
    //
    RecVar RecVarPassDens;
    RecVarPassDens.RecS.VarName1="DensPass";
    RecVarPassDens.RecS.VarName2="Density of passing points";
    RecVarPassDens.RecS.minval=0;
    RecVarPassDens.RecS.maxval=eLpNorm;
    RecVarPassDens.RecS.Unit="nondim.";
    RecVarPassDens.F=TheDens;
    //
    PLOT_PCOLOR(FileName, GrdArr, eDrw, RecVarPassDens, eCall, ePerm);
  }
  if (DensityFinalPlot) {
    MyMatrix<double> TheDens=ZeroMatrix<double>(nbSplitLon, nbSplitLat);
    for (int iDrifter=0; iDrifter<nbDrifter; iDrifter++) {
      double eLon=ListLon(nbTime-1, iDrifter);
      double eLat=ListLat(nbTime-1, iDrifter);
      int i=floor((eLon - MinLonPl)/DeltaLonWork);
      int j=floor((eLat - MinLatPl)/DeltaLatWork);
      if (i>= 0 && i<nbSplitLon && j>=0 && j<nbSplitLat) {
	TheDens(i,j)++;
      }
    }
    double maxVal=TheDens.maxCoeff();
    std::cerr << "maxVal=" << maxVal << "\n";
    std::string TitleStr="Density of final points";
    std::string FileName=ePerm.eDir + "Density_final_point";
    DrawArr eDrw;
    eDrw.DoTitle=true;
    eDrw.TitleStr=TitleStr;
    eDrw.DoColorBar=true;
    eDrw.eQuadFrame=eQuad;
    eDrw.DoTitle=false;
    //
    RecVar RecVarF;
    RecVarF.RecS.VarName1="DensFinal";
    RecVarF.RecS.VarName2="Density of final points";
    RecVarF.RecS.minval=0;
    RecVarF.RecS.maxval=maxVal;
    RecVarF.RecS.Unit="nondim.";
    RecVarF.F=TheDens;
    //
    PLOT_PCOLOR(FileName, GrdArr, eDrw, RecVarF, eCall, ePerm);
  }
}









struct TotalGridStruct {
  int nbGrid;
  std::vector<ROMSgridArray> ListGrdArr;
  std::vector<std::vector<PairLL>> ListBoundaries;
  std::vector<int> ListFatherGrid;
  std::vector<MyMatrix<int>> ListNestedStatus;
  std::vector<std::vector<int>> ListListChild;
  int iMainGrid;
  std::vector<MyMatrix<int>> ListChildBelonging;
};











struct SingleStateOcean {
  MyMatrix<double> Zeta;
  Eigen::Tensor<double,3> U;
  Eigen::Tensor<double,3> V;
  Eigen::Tensor<double,3> W;
};

// We have to adopt a shifting index strategy for rolling over
// between one state and the next. Otherwise, we need to do
// huge useles copies.

struct IntervalResolutionState {
  std::vector<int> ListIndex;
  std::vector<int> ListITime;
  std::vector<double> ListTimeDay;
  std::vector<SingleStateOcean> ListSingleState;
};


struct FullStateOcean {
  double minTimeDay;
  double maxTimeDay;
  std::vector<IntervalResolutionState> ListState;
};



struct SingleFloatEntry {
  double eTimeDay;
  double lon, lat, dep;
  double xgrd, ygrd, zGrid;
  int iGrid;
};

struct LocalFieldVar {
  double zeta;
  double u, v, w;
};


enum class RomsVarType {
  p2dvar = 1, r2dvar = 2,
    u2dvar = 3,
    v2dvar = 4,
    p3dvar = 5,
    r3dvar = 6,
    u3dvar = 7,
    v3dvar = 8,
    w3dvar = 9,
  b3dvar = 10
};



/*

LocalFieldVar GetVariables(SingleFloatEntry const & eEnt, FullStateOcean const& OceanState, TotalGridStruct const& eTotalGrid)
{
  int iGrid=eEnt.iGrid;
  double eTimeDay=eEnt.eTimeDay;
  int nbState=OceanState.ListState[iGrid].ListIndex.size();
  for (int jState=1; jState<nbState; jState++) {
    int iState=jState-1;
    int iIdx=OceanState.ListState[iGrid].ListIndex[iState];
    int jIdx=OceanState.ListState[iGrid].ListIndex[jState];
    double eTimeLow=OceanState.ListState[iGrid].ListTimeDay[iIdx];
    double eTimeUpp=OceanState.ListState[iGrid].ListTimeDay[jIdx];
    if (eTimeDay >= eTimeLow && eTimeUpp >= eTimeDay) {
      double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
      double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);

    }
  }
  std::cerr << "Failed to find entry. Please correct\n";
  throw TerminalException{1};
}
*/



std::vector<PairLL> GetRomsGridBoundary(ROMSgridArray const& eGrdArr)
{
  int eta_rho=eGrdArr.LON_rho.rows();
  int xi_rho =eGrdArr.LON_rho.cols();
  int iStart, iEnd, jStart, jEnd;
  iStart=0;
  iEnd=eta_rho-1;
  jStart=0;
  jEnd=xi_rho-1;
  return GetGridBoundary(eGrdArr.LON_rho, eGrdArr.LAT_rho, iStart, iEnd, jStart, jEnd);
}






TotalGridStruct ComputeTotalGridStruct(std::vector<std::string> const& ListGridFile, std::vector<int> const& ListFatherGrid)
{
  int nbGrid=ListGridFile.size();
  std::vector<ROMSgridArray> ListGrdArr(nbGrid);
  std::vector<std::vector<PairLL>> ListBoundaries(nbGrid);
  std::vector<std::vector<int>> ListListChild(nbGrid);
  int nbTopGrid=0;
  int iMainGrid=-1;
  for (int iGrid=0; iGrid<nbGrid; iGrid++) {
    ListGrdArr[iGrid]=ReadFullROMSgridArray(ListGridFile[iGrid]);
    ListBoundaries[iGrid]=GetRomsGridBoundary(ListGrdArr[iGrid]);
    int jGridFather=ListFatherGrid[iGrid];
    if (jGridFather != -1) {
      ListListChild[jGridFather].push_back(iGrid);
    }
    else {
      nbTopGrid++;
      iMainGrid=iGrid;
    }
  }
  if (nbTopGrid != 1) {
    std::cerr << "The number of top grids should be exactly 1\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<int>> ListChildBelonging(nbGrid);
  for (int iGrid=0; iGrid<nbGrid; iGrid++) {
    int eta_rho=ListGrdArr[iGrid].LON_rho.rows();
    int xi_rho=ListGrdArr[iGrid].LON_rho.cols();
    MyMatrix<int> MatrixBelong(eta_rho, xi_rho);
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	int eBelong=-1;
	PairLL ePt{ListGrdArr[iGrid].LON_rho(i,j), ListGrdArr[iGrid].LAT_rho(i,j)};
	for (auto &iGridChild : ListListChild[iGrid]) {
	  bool test=IsPointInside_Point(ePt, ListBoundaries[iGridChild]);
	  if (test && eBelong != -1) {
	    std::cerr << "Some children grids are overlapping, which is not allowed";
	    throw TerminalException{1};
	  }
	  if (test) {
	    if (eBelong != -1) {
	      std::cerr << "Error in eBelong. A master grid point should belog only to\n";
	      std::cerr << "one nested grid\n";
	      std::cerr << "We have iGridChild=" << iGridChild << "\n";
	      std::cerr << "eBelong=" << eBelong << "\n";
	      throw TerminalException{1};
	    }
	    eBelong=iGridChild;
	  }
	  MatrixBelong(i,j)=eBelong;
	}
      }
    ListChildBelonging[iGrid]=MatrixBelong;
  }
  TotalGridStruct eTotal;
  eTotal.nbGrid=nbGrid;
  eTotal.ListGrdArr=ListGrdArr;
  eTotal.ListBoundaries=ListBoundaries;
  eTotal.ListFatherGrid=ListFatherGrid;
  //  eTotal.ListNestedStatus=ListNestedStatus;
  eTotal.ListListChild=ListListChild;
  eTotal.iMainGrid=iMainGrid;
  eTotal.ListChildBelonging=ListChildBelonging;
  return eTotal;
};




int FindContainingGrid(TotalGridStruct const& eTotalGrid, PairLL const& ePt)
{
  int iGrid=eTotalGrid.iMainGrid;
  while(true) {
    bool test1=IsPointInside_Point(ePt, eTotalGrid.ListBoundaries[iGrid]);
    if (!test1) {
      return -1;
    }
    bool DoSomething=false;
    for (auto & iGridChild : eTotalGrid.ListListChild[iGrid]) {
      if (!DoSomething) {
	bool test2=IsPointInside_Point(ePt, eTotalGrid.ListBoundaries[iGridChild]);
	if (test2) {
	  DoSomething=true;
	  iGrid=iGridChild;
	}
      }
    }
    if (!DoSomething) {
      break;
    }
  }
  return iGrid;
}

struct PairCoordCurv {
  double dx;
  double dy;
};


/* Adapted from Utility/interpolate */
PairCoordCurv GetPosition(PairLL const& ePt, MyMatrix<double> const& LON, MyMatrix<double> const& LAT, MyMatrix<double> const& angler, PairCoord const& ePairCoord)
{
  double eLon=ePt.eLon;
  double eLat=ePt.eLat;
  int Imin=ePairCoord.i;
  int Jmin=ePairCoord.j;
  //
  double Eradius=6371315.0;
  double deg2rad=3.1415926535 / double(180);
  double yfac=Eradius*deg2rad;
  double xfac=yfac*cos(eLon*deg2rad);
  double xpp=(eLon - LON(Imin,Jmin))*xfac;
  double ypp=(eLat - LAT(Imin,Jmin))*yfac;
  auto sqr=[](double const& x) -> double {
    return x*x;
  };
  double diag2=sqr((LON(Imin+1,Jmin) - LON(Imin,Jmin+1))*xfac) +
    sqr((LAT(Imin+1,Jmin)-LAT(Imin,Jmin+1))*yfac);
  double aa2=sqr((LON(Imin,Jmin)-LON(Imin+1,Jmin))*xfac) +
    sqr((LAT(Imin,Jmin) - LAT(Imin+1,Jmin))*yfac);
  double bb2=sqr((LON(Imin,Jmin)-LON(Imin,Jmin+1))*xfac)+
    sqr((LAT(Imin,Jmin) - LAT(Imin,Jmin+1))*yfac);
  double phi=asin((diag2-aa2-bb2)/(2*sqrt(aa2*bb2)));
  double ang=angler(Imin,Jmin);
  double dx=xpp*cos(ang) + ypp*sin(ang);
  double dy=ypp*cos(ang) - xpp*sin(ang);
  dx += dy*tan(phi);
  dy /= cos(phi);
  dx /= sqrt(aa2);
  dy /= sqrt(bb2);
  auto fRenorm=[&](double const& x) -> double {
    if (x<0) return 0;
    if (x>1) return 1;
    return x;
  };
  dx=fRenorm(dx);
  dy=fRenorm(dy);
  return {dx, dy};
}






SingleFloatEntry GetFloatEntry(double const& eLon, double const& eLat, double const& eDep, TotalGridStruct const& eTotalGrid, double const& eTimeDay)
{
  PairLL ePt{eLon, eLat};
  int iGrid=FindContainingGrid(eTotalGrid, ePt);
  SingleFloatEntry eEnt;
  eEnt.iGrid=-1;
  eEnt.lon = eLon;
  eEnt.lat = eLat;
  eEnt.dep = eDep;
  if (iGrid == -1)
    return eEnt;
  PairCoord ePair=FindContaining(ePt, eTotalGrid.ListGrdArr[iGrid].LON_rho, eTotalGrid.ListGrdArr[iGrid].LAT_rho);
  PairCoordCurv ePairCoordCurv=GetPosition(ePt, eTotalGrid.ListGrdArr[iGrid].LON_rho, eTotalGrid.ListGrdArr[iGrid].LAT_rho, eTotalGrid.ListGrdArr[iGrid].angle, ePair);
  double dx=ePairCoordCurv.dx;
  double dy=ePairCoordCurv.dy;
  int i=ePair.i;
  int j=ePair.j;
  double xgrd=double(i) + dx;
  double ygrd=double(j) + dy;
  eEnt.xgrd=xgrd;
  eEnt.ygrd=ygrd;
  double eDepSea=(1-dx)*(1-dy)*eTotalGrid.ListGrdArr[iGrid].h(i,j) +
    dx*(1-dy)*eTotalGrid.ListGrdArr[iGrid].h(i+1,j) +
    (1-dx)*dy*eTotalGrid.ListGrdArr[iGrid].h(i,j+1) +
    dx*dy*eTotalGrid.ListGrdArr[iGrid].h(i+1,j+1);
  int N=eTotalGrid.ListGrdArr[iGrid].ARVD.N;
  VerticalInfo eVert=GetVerticalInfo(N);
  int eta_rho=eTotalGrid.ListGrdArr[iGrid].LON_rho.rows();
  int xi_rho=eTotalGrid.ListGrdArr[iGrid].LON_rho.cols();
  MyMatrix<double> zeta;
  zeta.setZero(eta_rho, xi_rho);
  ComputeHz(eTotalGrid.ListGrdArr[iGrid].ARVD, eTotalGrid.ListGrdArr[iGrid].h(i,j), zeta(i,j), eVert);
  auto GetZgrid=[&](double const& zfloat) -> double {
    for (int k=0; k<N; k++)
      if (eVert.z_w(k) <= zfloat && zfloat <= eVert.z_w(k+1))
	return double(k) + (zfloat-eVert.z_w(k))/eVert.Hz(k);
    std::cerr << "Failed to find position\n";
    throw TerminalException{1};
  };
  eEnt.zGrid = GetZgrid(eDepSea);
  return eEnt;
}










struct FrameTime {
  double BaseTime;
  double DeltaTime;
};


// We select a frame of times together and then do the computation
// on the interval [A, A + T]
// a frame is correct if for A and A+T all models contain a state exactly
// on that point.
FrameTime ComputeFrameTime(std::vector<TotalArrGetData> const& ListRec)
{
  int nbArr=ListRec.size();
  ArrayHistory eArrRef=ListRec[0].eArr;
  int nbTime=eArrRef.nbTime;
  struct Answer {
    double eTime;
    int iTime;
  };
  auto fAnswer=[&](int const& timeStart) -> Answer {
    for (int iTime=timeStart; iTime<nbTime; iTime++) {
      double eTimeDay=eArrRef.ListTime[iTime];
      bool IsCorrect=true;
      for (int iArr=1; iArr<nbArr; iArr++)
	if (IsCorrect) {
	  InterpInfo eInterpInfo=GetTimeInterpolationInfo(ListRec[iArr].eArr.ListTime, eTimeDay);
	  if (!eInterpInfo.UseSingleEntry)
	    IsCorrect=false;
	}
      if (IsCorrect)
	return {eTimeDay, iTime};
    }
    return {double(-1), -1};
  };
  Answer eAns1=fAnswer(0);
  Answer eAns2=fAnswer(eAns1.iTime + 1);
  return {eAns1.eTime, eAns2.eTime - eAns1.eTime};
}



FullStateOcean RetrieveFullStateOcean(double const& minTimeDay, double const& maxTimeDay, std::vector<TotalArrGetData> const& ListRec)
{
  int nbArr=ListRec.size();
  std::vector<IntervalResolutionState> ListState(nbArr);
  double tolDay=double(1)/double(100000);
  for (int iArr=0; iArr<nbArr; iArr++) {
    std::vector<int> ListIndex;
    std::vector<int> ListITime;
    std::vector<double> ListTimeDay;
    std::vector<SingleStateOcean> ListSingleState;
    int nbTime=ListRec[iArr].eArr.ListTime.size();
    int idx=0;
    for (int iTime=0; iTime<nbTime; iTime++) {
      double eTimeDay=ListRec[iArr].eArr.ListTime[iTime];
      if (eTimeDay > minTimeDay - tolDay && eTimeDay < maxTimeDay + tolDay) {
	ListITime.push_back(iTime);
	ListIndex.push_back(idx);
	idx++;
	ListTimeDay.push_back(eTimeDay);
	int iFile=ListRec[iArr].eArr.ListIFile[iTime];
	int iRec=ListRec[iArr].eArr.ListIRec[iTime];
	std::string HisFile=ListRec[iArr].eArr.ListFileNames[iFile];
	MyMatrix<double> zeta=NETCDF_Get2DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "zeta", iRec);
	Eigen::Tensor<double,3> U=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "u", iRec);
	Eigen::Tensor<double,3> V=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "v", iRec);
	Eigen::Tensor<double,3> W=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "w", iRec);
	ListSingleState.push_back({zeta, U, V, W});
      }
    }

    ListState.push_back({ListIndex, ListITime, ListTimeDay, ListSingleState});
  }
  return {minTimeDay, maxTimeDay, ListState};
}


void ShiftFullStateOcean(double const& minTimeDay, double const& maxTimeDay, std::vector<TotalArrGetData> const& ListRec, FullStateOcean & TheState)
{
  TheState.minTimeDay=minTimeDay;
  TheState.maxTimeDay=maxTimeDay;
  int nbArr=ListRec.size();
  double tolDay=double(1)/double(100000);
  for (int iArr=0; iArr<nbArr; iArr++) {
    std::vector<int> ListITime;
    std::vector<double> ListTimeDay;
    std::vector<SingleStateOcean> ListSingleState;
    int nbTime=ListRec[iArr].eArr.ListTime.size();
    for (int iTime=0; iTime<nbTime; iTime++) {
      double eTimeDay=ListRec[iArr].eArr.ListTime[iTime];
      if (eTimeDay > minTimeDay - tolDay && eTimeDay < maxTimeDay + tolDay) {
	ListITime.push_back(iTime);
	ListTimeDay.push_back(eTimeDay);
      }
    }
    int nbTimeRel=ListITime.size();
    int nbTimeOld=TheState.ListState[iArr].ListITime.size();
    if (nbTimeRel != nbTimeOld) {
      std::cerr << "We have a logical problem in the database system. Please correct\n";
      throw TerminalException{1};
    }
    std::vector<int> ListStatusRel(nbTimeRel,-1);
    std::vector<int> ListStatusOld(nbTimeRel,-1);
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++)
      for (int iTimeOld=0; iTimeOld<nbTimeOld; iTimeOld++) {
	if (ListITime[iTimeRel] == TheState.ListState[iArr].ListITime[iTimeOld]) {
	  TheState.ListState[iArr].ListIndex[iTimeRel]=iTimeOld;
	  ListStatusRel[iTimeRel]=0;
	  ListStatusOld[iTimeOld]=0;
	}
      }
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++)
      if (ListStatusRel[iTimeRel] == -1) {
	bool WeAreDone=false;
	for (int iTimeOld=0; iTimeOld<nbTimeOld; iTimeOld++)
	  if (!WeAreDone && ListStatusOld[iTimeOld] == -1) {
	    WeAreDone=true;
	    TheState.ListState[iArr].ListIndex[iTimeRel]=iTimeOld;
	    TheState.ListState[iArr].ListTimeDay[iTimeOld]=ListTimeDay[iTimeRel];
	    int iTimeFound=ListITime[iTimeRel];
	    TheState.ListState[iArr].ListITime[iTimeOld]=iTimeFound;
	    int iFile=ListRec[iArr].eArr.ListIFile[iTimeFound];
	    int iRec=ListRec[iArr].eArr.ListIRec[iTimeFound];
	    std::string HisFile=ListRec[iArr].eArr.ListFileNames[iFile];
	    TheState.ListState[iArr].ListSingleState[iTimeOld].Zeta=NETCDF_Get2DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "zeta", iRec);
	    TheState.ListState[iArr].ListSingleState[iTimeOld].U=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "u", iRec);
	    TheState.ListState[iArr].ListSingleState[iTimeOld].V=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "v", iRec);
	    TheState.ListState[iArr].ListSingleState[iTimeOld].W=NETCDF_Get3DvariableSpecEntry(HisFile, ListRec[iArr].GrdArr, "w", iRec);
	  }
      }
  }
}




void FLOAT_CreateNetcdfFile(std::string const& eFileNC, int const& nbFloat)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  netCDF::NcDim eDimOcean=dataFile.addDim("ocean_time");
  netCDF::NcDim eDimDateString=dataFile.addDim("dateString", 19);
  netCDF::NcDim eDimFloat=dataFile.addDim("float", size_t(nbFloat));
  std::vector<std::string> LDim1{"ocean_time"};
  std::vector<std::string> LDim2{"ocean_time", "dateString"};
  std::vector<std::string> LDim3{"ocean_time", "float"};
  netCDF::NcVar eVAR_time    = dataFile.addVar("ocean_time"    , "double", LDim1);
  netCDF::NcVar eVAR_time_day= dataFile.addVar("ocean_time_day", "double", LDim1);
  netCDF::NcVar eVAR_time_str= dataFile.addVar("ocean_time_str", "char"  , LDim2);
  netCDF::NcVar eVAR_lon     = dataFile.addVar("lon", "double", LDim3);
  netCDF::NcVar eVAR_lat     = dataFile.addVar("lat", "double", LDim3);
  netCDF::NcVar eVAR_dep     = dataFile.addVar("dep", "double", LDim3);
}

void FLOAT_AppendFloatsNetcdfFile(std::string const& eFileNC, std::vector<SingleFloatEntry> const& ListFloatEntry, double const& eTimeDay)
{
  //  int nbEnt=ListTimeDay.size();
  int nbFloat=ListFloatEntry.size();
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::write);
  std::string eTime="ocean_time";
  std::string eTimeStr="ocean_time_str";
  netCDF::NcVar eVAR_d=dataFile.getVar(eTime);
  netCDF::NcVar eVAR_char=dataFile.getVar(eTimeStr);
  std::multimap<std::string,netCDF::NcDim> MapDims=dataFile.getDims();
  std::multimap<std::string,netCDF::NcDim>::iterator iter=MapDims.begin();
  netCDF::NcDim timeDim;
  while (iter != MapDims.end()) {
    if (iter->first == "ocean_time")
      timeDim=iter->second;
    iter++;
  }
  //  netCDF::NcDim timeDim=MapDims.at("ocean_time");
  if (!timeDim.isUnlimited()) {
    std::cerr << "Error the dimension should be unlimited\n";
    throw TerminalException{1};
  }
  size_t siz=timeDim.getSize();
  double *Alon, *Alat, *Adep;
  Alon=new double[nbFloat];
  Alat=new double[nbFloat];
  Adep=new double[nbFloat];
  netCDF::NcVar eVar_lon=dataFile.getVar("lon");
  netCDF::NcVar eVar_lat=dataFile.getVar("lat");
  netCDF::NcVar eVar_dep=dataFile.getVar("dep");
  //
  std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
  std::vector<size_t> start2{siz};
  std::vector<size_t> count2{1};
  eVAR_d.putVar(start2, count2, &eTimeDay);
  std::vector<size_t> start3{siz,0};
  std::vector<size_t> count3{1,19};
  eVAR_char.putVar(start3, count3, strPres.c_str());
  //
  for (int iFloat=0; iFloat<nbFloat; iFloat++) {
    Alon[iFloat]=ListFloatEntry[iFloat].lon;
    Alat[iFloat]=ListFloatEntry[iFloat].lat;
    Adep[iFloat]=ListFloatEntry[iFloat].dep;
  }
  std::vector<size_t> start{siz,0};
  std::vector<size_t> count{1, size_t(nbFloat)};
  eVar_lon.putVar(start, count, Alon);
  eVar_lat.putVar(start, count, Alat);
  eVar_dep.putVar(start, count, Adep);
  delete [] Alon;
  delete [] Alat;
  delete [] Adep;
}




std::vector<PairLL> GetListPairLL(double const& eLon, double const& eLat, double const& DistKM, int const& NbAsk)
{
  double EarthRadiusKM=6370;
  std::vector<PairLL> ListSolution(NbAsk);
  double deltaLL=(180/GetPI())*(DistKM/EarthRadiusKM);
  int idx=0;
  while(true) {
    double x1 = ((double) rand() / (RAND_MAX));
    double x2 = ((double) rand() / (RAND_MAX));
    double NewLon = eLon + (2*x1-1)*deltaLL;
    double NewLat = eLat + (2*x2-1)*deltaLL;
    double eDistPartKM=GeodesicDistanceKM(eLon, eLat, NewLon, NewLat);
    if (eDistPartKM <= DistKM) {
      ListSolution[idx]={NewLon, NewLat};
      idx++;
    }
    if (idx == NbAsk)
      break;
  }
  return ListSolution;
}




std::vector<SingleFloatEntry> GetListSingleFloatEntry(TotalGridStruct const& eTotalGrid, FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::vector<double> ListLONpt=eBlPROC.ListListDoubleValues.at("ListLonFloat");
  std::vector<double> ListLATpt=eBlPROC.ListListDoubleValues.at("ListLatFloat");
  std::vector<double> ListDepth=eBlPROC.ListListDoubleValues.at("ListDepth");
  std::vector<double> ListTime=eBlPROC.ListListDoubleValues.at("ListTime");
  int NbAsk=eBlPROC.ListIntValues.at("NbAsk");
  double DistKM=eBlPROC.ListDoubleValues.at("DistKM");
  int nbPart=ListLONpt.size();
  //  int nbDep=ListDepth.size();
  //  int nbTime=ListTime.size();
  std::vector<SingleFloatEntry> ListEnt;
  for (int iPart=0; iPart<nbPart; iPart++) {
    double eLon=ListLONpt[iPart];
    double eLat=ListLATpt[iPart];
    std::vector<PairLL> ListPairLL=GetListPairLL(eLon, eLat, DistKM, NbAsk);
    for (auto & ePair : ListPairLL) {
      double eLon1=ePair.eLon;
      double eLat1=ePair.eLat;
      for (auto & eDep : ListDepth)
	for (auto & eTimeDay : ListTime) {
	  SingleFloatEntry eEnt=GetFloatEntry(eLon1, eLat1, eDep, eTotalGrid, eTimeDay);
	  ListEnt.push_back(eEnt);
	}
    }
  }
  return ListEnt;
}

struct FloatInfo {
  double dt;
  int nbTime;
  int nbOutSeq;
};



// We take the FloatEntry as entries and in return it is the final position.
// eEnt is at the beginning the float at the starting position.
// The TargetFinal is built over the computation.
/*
void FloatEvolution(SingleFloatEntry & eEnt, std::vector<SingleFloatEntry> & TargetFinal, std::vector<TotalArrGetData> const& ListRec, TotalGridStruct const& eTotalGrid, FloatInfo const& eFlInfo)
{
  double dt=eFlInfo.dt;
  int nbTime=eFlInfo.nbTime;
  int nbOutSeq=eFlInfo.nbOutSeq;
  int idx=0;
  int iGrid=eEnt.iGrid;


  for (int iTime=0; iTime<nbTime; iTime++) {
    int res=iTime%nbOutSeq;
    if (res == 0) {
      TargetFinal[idx]=eEnt;
      idx++;
    }
    int Ir=int(eEnt.xgrd);
    int Jr=int(eEnt.xgrd);
  }

}
*/


 /*
void ProcessFloatComputation(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::vector<std::string> ListGridFile=eBlPROC.ListListStringValues.at("ListGridFile");
  std::vector<std::string> ListHisPrefix=eBlPROC.ListListStringValues.at("ListHisPrefix");
  std::vector<int> ListFatherGrid=eBlPROC.ListListIntValues.at("ListFatherGrid");
  std::string OutFile=eBlPROC.ListStringValues.at("OutFile");
  double DEFINETC=eBlPROC.ListDoubleValues.at("DEFINETC");
  double HISTIME=eBlPROC.ListDoubleValues.at("HISTIME");
  double dt=eBlPROC.ListDoubleValues.at("dt");
  TotalGridStruct eTotalGrid=ComputeTotalGridStruct(ListGridFile, ListFatherGrid);
  int nbGrid=ListGridFile.size();
  std::vector<TotalArrGetData> ListRec(nbGrid);
  std::vector<SingleFloatEntry> ListFloatEntry=GetListSingleFloatEntry(eTotalGrid, eFull);
  int nbFloat=ListFloatEntry.size();
  for (int iGrid=0; iGrid<nbGrid; iGrid++) {
    std::string eModelName="ROMS";
    std::string GridFile=ListGridFile[iGrid];
    std::string BoundFile="unset";
    std::string HisPrefix=ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, BoundFile, HisPrefix, {}};
    GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
    ArrayHistory eArr=ReadArrayHistory(eTriple);
    TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
    ListRec[iGrid]=TotalArr;
  }
  FrameTime eFrame=ComputeFrameTime(ListRec);
  // The data is stored in an interval [eTime, eTime + eFrame.DeltaTime]
  std::cerr << " BaseTime=" << eFrame.BaseTime << "\n";
  std::cerr << "DeltaTime=" << eFrame.DeltaTime << "\n";
  double DeltaTime=eFrame.DeltaTime;
  std::vector<double> ListTime=GetInterval(eBlPROC.ListStringValues.at("BEGTC"),
					   eBlPROC.ListStringValues.at("ENDTC"),
					   eBlPROC.ListDoubleValues.at("DELTC"),
					   eBlPROC.ListStringValues.at("UNITC"));
  int nbTimeTot=ListTime.size();
  double eTimeLast=ListTime[nbTimeTot-1];
  double eTimeFirst=ListTime[0];
  double DeltaTimeTot=eTimeLast - eTimeFirst;
  double DeltaTimeSing=ListTime[1] - eTimeFirst;
  double minTimeDay=eFrame.BaseTime;
  double maxTimeDay=eFrame.BaseTime + DeltaTime;
  double tolDay=double(1)/double(100000);
  while(true) {
    if (eTimeFirst > minTimeDay - tolDay && eTimeFirst < maxTimeDay + tolDay) {
      break;
    }
    if (eTimeFirst < minTimeDay) {
      minTimeDay -= DeltaTime;
      maxTimeDay -= DeltaTime;
    }
    if (eTimeFirst > maxTimeDay) {
      minTimeDay += DeltaTime;
      maxTimeDay += DeltaTime;
    }
  }
  FullStateOcean TheState=RetrieveFullStateOcean(minTimeDay, maxTimeDay, ListRec);
  int nbBlock=int(eFrame.DeltaTime / HISTIME);
  std::vector<std::vector<SingleFloatEntry>> ListListFloatEntry(nbBlock);
  for (int iBlock=0; iBlock<nbBlock; iBlock++)
    ListListFloatEntry[iBlock].resize(nbFloat);
  int nbTime=int(eFrame.DeltaTime / DeltaTimeSing);
  int nbTimeFrame=int(DeltaTimeTot / eFrame.DeltaTime);
  int nbOutSeq=int(HISTIME / DeltaTimeSing);
  //
  int sizNetcdf = int(DEFINETC / HISTIME);
  int posNetcdf = 0;
  int iNetcdfFile=0;
  std::string eFileNC;
  std::vector<double> ListTimeDay;
  auto CreateNetcdfFile=[&]() -> void {
    iNetcdfFile++;
    eFileNC = OutFile + StringNumber(iNetcdfFile, 4) + ".nc";
    FLOAT_CreateNetcdfFile(eFileNC, nbFloat);
  };
  auto SingleWriteNetcdf=[&](std::vector<SingleFloatEntry> const& eListFloatEnt, double const& eTimeDay) -> void {
    if (posNetcdf == 0)
      CreateNetcdfFile();
    FLOAT_AppendFloatsNetcdfFile(eFileNC, eListFloatEnt, eTimeDay);
    posNetcdf++;
    if (posNetcdf == sizNetcdf)
      posNetcdf=0;
  };
  auto WriteNetcdfFile=[&]() -> void {
    for (int iBlock=0; iBlock<nbBlock; iBlock++)
      SingleWriteNetcdf(ListListFloatEntry[iBlock], ListTimeDay[iBlock]);
  };
  FloatInfo eFlInfo;
  eFlInfo.dt = dt;
  eFlInfo.nbTime=nbTime;
  eFlInfo.nbOutSeq=nbOutSeq;
  std::vector<SingleFloatEntry> SingleEvolution(nbBlock);
  for (int iTimeFrame=0; iTimeFrame<nbTimeFrame; iTimeFrame++) {
    double eTimeDay = ListTime[iTimeFrame];
    //
    for (int iFloat=0; iFloat<nbFloat; iFloat++) {
      FloatEvolution(ListFloatEntry[iFloat], SingleEvolution, ListRec, eTotalGrid, eFlInfo);
      for (int iBlock=0; iBlock<nbBlock; iBlock++)
	ListListFloatEntry[iBlock][iFloat]=SingleEvolution[iBlock];
    }
    WriteNetcdfFile();
    ShiftFullStateOcean(minTimeDay, maxTimeDay, ListRec, TheState);
  }
  SingleWriteNetcdf(ListFloatEntry, eTimeLast);
}
 */





#endif
