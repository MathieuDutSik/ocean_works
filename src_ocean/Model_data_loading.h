#ifndef MODEL_DATA_LOADING_INCLUDE
#define MODEL_DATA_LOADING_INCLUDE

#include "Model_grids.h"
#include "SmoothingBathyBasic.h"
#include "ROMSfunctionality.h"


struct PairMinMax {
  double TheMin;
  double TheMax;
};


PairMinMax ComputeMinMax(GridArray const& GrdArr, MyMatrix<double> const& F)
{
  //  std::cerr << "Begin of ComputeMinMax\n";
  bool IsFirst=true;
  int eta_rho=F.rows();
  int xi_rho=F.cols();
  int eta_rho_msk=GrdArr.GrdArrRho.MSK.rows();
  int xi_rho_msk=GrdArr.GrdArrRho.MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho     << " xi_rho=" << xi_rho     << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk << "\n";
    throw TerminalException{1};
  }
  double TheMin=0;
  double TheMax=0;
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eMsk;
      if (GrdArr.IsFE == 1) {
	eMsk=1;
      }
      else {
	eMsk=GrdArr.GrdArrRho.MSK(i,j);
      }
      if (eMsk == 1) {
	double eVal=F(i,j);
	if (IsFirst) {
	  TheMin=eVal;
	  TheMax=eVal;
	}
	else {
	  if (eVal < TheMin)
	    TheMin=eVal;
	  if (eVal > TheMax)
	    TheMax=eVal;
	}
	IsFirst=false;
      }
    }
  //  std::cerr << "  End of ComputeMinMax\n";
  return {TheMin, TheMax};
}




PairMinMax ComputeMinMaxMask(MyMatrix<int> const& MSK, MyMatrix<double> const& F)
{
  bool IsFirst=true;
  int eta_rho=F.rows();
  int xi_rho=F.cols();
  int eta_rho_msk=MSK.rows();
  int xi_rho_msk=MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho     << " xi_rho=" << xi_rho     << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk << "\n";
    throw TerminalException{1};
  }
  double TheMin=0;
  double TheMax=0;
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eMsk=MSK(i,j);
      if (eMsk == 1) {
	double eVal=F(i,j);
	if (IsFirst) {
	  TheMin=eVal;
	  TheMax=eVal;
	}
	else {
	  if (eVal < TheMin)
	    TheMin=eVal;
	  if (eVal > TheMax)
	    TheMax=eVal;
	}
	IsFirst=false;
      }
    }
  return {TheMin, TheMax};
}



PairMinMax ComputeMinMax_3D(GridArray const& GrdArr, Eigen::Tensor<double,3> const& F)
{
  //  std::cerr << "Begin of ComputeMinMax\n";
  auto LDim=F.dimensions();
  int nVert=LDim[0];
  int eta_rho=LDim[1];
  int xi_rho=LDim[2];
  int eta_rho_msk=GrdArr.GrdArrRho.MSK.rows();
  int xi_rho_msk=GrdArr.GrdArrRho.MSK.cols();
  if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
    std::cerr << "ComputeMinMax_3D error : Inconsistency in dimension\n";
    std::cerr << "  F: eta_rho=" << eta_rho     << " xi_rho=" << xi_rho     << "\n";
    std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk << "\n";
    throw TerminalException{1};
  }
  bool IsFirst=true;
  double TheMin=0;
  double TheMax=0;
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eMsk;
      if (GrdArr.IsFE == 1)
	eMsk=1;
      else
	eMsk=GrdArr.GrdArrRho.MSK(i,j);
      if (eMsk == 1) {
	for (int iVert=0; iVert<nVert; iVert++) {
	  double eVal=F(iVert, i,j);
	  if (IsFirst) {
	    TheMin=eVal;
	    TheMax=eVal;
	  }
	  else {
	    if (eVal < TheMin)
	      TheMin=eVal;
	    if (eVal > TheMax)
	      TheMax=eVal;
	  }
	  IsFirst=false;
	}
      }
    }
  //  std::cerr << "  End of ComputeMinMax\n";
  return {TheMin, TheMax};
}






MyMatrix<double> COMPUTE_NORM(MyMatrix<double> const& U, MyMatrix<double> const& V)
{
  int eta=U.rows();
  int xi=U.cols();
  MyMatrix<double> WINDMAG(eta, xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eU=U(i,j);
      double eV=V(i,j);
      WINDMAG(i,j)=sqrt(eU*eU + eV*eV);
    }
  return WINDMAG;
}

MyMatrix<double> FreqPeriodChange(MyMatrix<double> const& F)
{
  int eta=F.rows();
  int xi=F.cols();
  MyMatrix<double> Fret(eta, xi);
  double pi=3.1415926535;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal=F(i,j);
      double NewVal=2*pi/eVal;
      Fret(i,j)=NewVal;
    }
  return Fret;
}



bool TOTALARR_IsVar(TotalArrGetData const& TotalArr, std::string const& eVar)
{
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    std::string HisFile=ARR_GetFirstFileName(TotalArr.eArr);
    return NC_IsVar(HisFile, eVar);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    for (auto & eVarName : TotalArr.eArr.RawVarNames)
      if (eVarName == eVar)
	return true;
    return false;
  }
  std::cerr << "Error in TOTALARR_IsVar\n";
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}

std::string TOTALARR_GetUnits(TotalArrGetData const& TotalArr, std::string const& eVar)
{
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    std::string HisFile=ARR_GetFirstFileName(TotalArr.eArr);
    netCDF::NcFile dataFile(HisFile, netCDF::NcFile::read);
    netCDF::NcVar data=dataFile.getVar(eVar);
    netCDF::NcVarAtt eAtt=data.getAtt("units");
    char eString[1024]="";
    eAtt.getValues(eString);
    return std::string(eString);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    for (auto & eRec : TotalArr.eArr.ListListMessages[0])
      if (eRec.shortName == eVar)
	return eRec.units;
    std::cerr << "Failed to find the entry\n";
    throw TerminalException{1};
  }
  std::cerr << "Error in TOTALARR_GetUnits\n";
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}







MyMatrix<double> Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& VarName, double const& eTimeDay)
{
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    //    std::cerr << "Before call to NETCDF_Get2DvariableSpecTime\n";
    return NETCDF_Get2DvariableSpecTime(TotalArr, VarName, eTimeDay);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    return GRIB_Get2DvariableSpecTime(TotalArr, "shortName", VarName, eTimeDay);
  }
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  throw TerminalException{1};
}



std::string DecltypeString(std::string const& FullVarName)
{
  std::vector<std::string> ListStr=STRING_Split(FullVarName, ":");
  return ListStr[0];
}

// From Ivica code.
MyMatrix<double> mixing_ratio2relative_humidity(MyMatrix<double> const& Q2, MyMatrix<double> const& PSFC, MyMatrix<double> const& T2K)
{
  int eta=Q2.rows();
  int xi=Q2.cols();
  MyMatrix<double> retField(eta,xi);
  double t_sat=373.16;
  double c1=7.90298;
  double c2=5.02808;
  double c3=1.3816e-7;
  double c4=11.344;
  double c5=8.1328e-3;
  double c6=3.49149;
  double water_air_mass_ratio = 0.621970585;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double mix_ratio=Q2(i,j)*1000;
      double pres=PSFC(i,j)/100;
      double temp=T2K(i,j);
      double t_ratio=t_sat/temp;
      double rt_ratio=1/t_ratio;
      double sl_pressure=1013.246;
      double svp1 =  -1.0 * c1 * ( t_ratio - 1.0 ) 
	+  c2 * log10( t_ratio ) 
	-  c3 * ( pow(10.0, c4 * ( 1.0 - rt_ratio ) ) - 1.0 ) 
	+  c5 * ( pow(10.0, -1.0 * c6 * ( t_ratio - 1.0 ) ) - 1.0 ) 
	+  log10( sl_pressure );
      double svp=pow(10.0 , svp1);
      double vapor_pressure_ratio = svp / ( pres - svp );
      double smr = 1000.0 * water_air_mass_ratio * vapor_pressure_ratio;
      double rh=(mix_ratio / smr) * 100;
      retField(i,j)=rh;
    }
  return retField;
}





MyMatrix<double> ConvertSpecHumid_to_RelativeHumidity(TotalArrGetData const& TotalArr, std::string const& VarName2r, std::string const& VarName2t, std::string const& VarNameMSL, std::string const& VarNameQ, double const& eTimeDay)
{
  MyMatrix<double> F;
  if (TOTALARR_IsVar(TotalArr, VarName2r)) {
    F=Get2DvariableSpecTime(TotalArr, VarName2r, eTimeDay);
  }
  else {
    int MethodRH = 2;
    if (MethodRH == 1) {
      double airDens=1.225; // kg/m3 : Density of air, see https://en.wikipedia.org/wiki/Density_of_air
      double waterDens=0.804; // g/L = kg/m3 : Density of water vapor https://en.wikipedia.org/wiki/Water_vapor
      // We use formula from http://www.engineeringtoolbox.com/specific-relative-humidity-air-d_688.html
      // formula is then phi = 100 x /(0.622 * rho_ws/(rho - rho_ws) )
      // It seems not to work correctly.
      double quot=waterDens / (airDens - waterDens);
      MyMatrix<double> Fspecific=Get2DvariableSpecTime(TotalArr, VarNameQ, eTimeDay);
      double TheMult=100 / ( 0.622 * quot);
      F = Fspecific * TheMult;
    }
    if (MethodRH == 2) {
      // We use formula from
      // http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
      MyMatrix<double> F_q=Get2DvariableSpecTime(TotalArr, VarNameQ, eTimeDay);
      MyMatrix<double> F_p;
      if (TOTALARR_IsVar(TotalArr, VarNameMSL)) {
	F_p=Get2DvariableSpecTime(TotalArr, VarNameMSL, eTimeDay);
      }
      else {
	int eta=F_q.rows();
	int xi=F_q.cols();
	F_p.setConstant(eta, xi, 103000);
      }
      MyMatrix<double> F_TK=Get2DvariableSpecTime(TotalArr, VarName2t, eTimeDay);
      int eta=F_q.rows();
      int xi=F_q.cols();
      F=MyMatrix<double>(eta,xi);
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++) {
	  double eT=F_TK(i,j);
	  double eQ=F_q(i,j);
	  double eP=F_p(i,j);
	  double eT0=double(273.15);
	  double TheQuot=double(17.67) * (eT - eT0)/(eT - double(29.65));
	  double eRH=0.263 * eP *eQ /(exp(TheQuot));
	  F(i,j)=std::min(eRH, double(100));
	}
    }
  }
  return F;
}



std::string GetBasicModelName(std::string const& eModelName)
{
  if (eModelName == "ROMS_IVICA")
    return "ROMS";
  if (eModelName == "WWM_DAILY")
    return "WWM";
  return eModelName;
}


ARVDtyp TOTALARR_GetARVD(TotalArrGetData const& TotalArr)
{
  int iFile=0;
  std::string eFile=ARR_GetHisFileName(TotalArr.eArr, "irrelevant", iFile);
  return ReadROMSverticalStratification(eFile);
}


std::vector<std::string> GetAllPossibleVariables()
{
  std::vector<std::string> ListVarOut{
    "IOBP", "MAPSTA", "FieldOut1", "CFL1", "CFL2", "CFL3", "ThreeDfield1", "NbIterSolv", 
    "WIND10", "Uwind", "Vwind", "WINDMAG",
    "ChlorophyllConcOCI",
    "ChlorophyllConcOCX",
    "CalciteConc",
    "ParticleOrganicCarbon",
    "InstPhotosyntheticallyAvailableRad",
    "PhotosyntheticallyAvailableRad",
    "BottomCurr", "SurfStress", "SurfCurr", "UsurfCurr", "VsurfCurr", "SurfCurrMag",
    "Curr", "CurrMag", "HorizCurr",
    "CurrBaro", "CurrBaroMag", "ChlorophylA", 
    "Temp", "Salt", "HorizTemp", "TempSurf", "TempBottom", "SaltSurf",
    "AIRT2", "AIRT2K", "Rh2", "Rh2frac", "AIRD", "SurfPres", 
    "ZetaOcean", "ZetaOceanDerivative", "DynBathy", "Bathymetry", "RoughnessFactor", "ZetaSetup", 
    "CdWave", "AlphaWave", "AirZ0", "AirFricVel", "CGwave", 
    "shflux", "ssflux", "evaporation", "CloudFraction", 
    "Hwave", "BreakingFraction", 
    "rain", "swrad", "lwrad", "latent", "sensible", 
    "MeanWaveFreq", "PeakWaveFreq", "TM02", 
    "MeanWavePer", "PeakWavePer", 
    "MeanWaveDirSpread", "PeakWaveDirSpread", 
    "MeanWaveDir", "PeakWaveDir", "MeanWaveDirVect", "PeakWaveDirVect", 
    "DiscPeakWaveDir",
    "MeanWaveLength", "PeakWaveLength", "MeanWaveNumber", "PeakWaveNumber", 
    "TotSurfStr", "WaveSurfStr", "SurfStrHF",
    "oxygen", "PO4", "NO3", "NH4", "NitrogenSink", "SiOH4", "ReductionEquivalent",
    "bacteriaC", "bacteriaN", "bacteriaP",
    "diatomsC", "diatomsN", "diatomsP", "diatomsL", "diatomsS",
    "flagellatesC", "flagellatesN", "flagellatesP", "flagellatesL",
    "picophytoplanktonC", "picophytoplanktonN", "picophytoplanktonP", "picophytoplanktonL",
    "largephytoplanktonC", "largephytoplanktonN", "largephytoplanktonP", "largephytoplanktonL",
    "CarnPlanktonC", "CarnPlanktonN", "CarnPlanktonP",
    "OmniPlanktonC", "OmniPlanktonN", "OmniPlanktonP",
    "MicroPlanktonC", "MicroPlanktonN", "MicroPlanktonP", 
    "HeteroNanoflagelattesC", "HeteroNanoflagelattesN", "HeteroNanoflagelattesP",
    "LabileDOM1c", "LabileDOM1n", "LabileDOM1p",
    "LabileDOM2c", "RefractoryDOMc",
    "ParticleOMc", "ParticleOMn", "ParticleOMp", "ParticleOMs",
    "DissolvedICc", "DissolvedICh",
    "Irradiance", "DIC",
    "chlorophyll",
    "NetProductionP1", "NetProductionP2", "NetProductionP3", "NetProductionP4",
    "RegFactorP1", "RegFactorP2", "RegFactorP3", "RegFactorP4",
    "GrossPP", "SecondPP", "ExtinctionCoeff"};
  return ListVarOut;
}



Eigen::Tensor<double,3> RetrieveStandardVerticalCoordinate(TotalArrGetData const& TotalArr)
{
  std::string eModelName=GetBasicModelName(TotalArr.GrdArr.ModelName);
  if (eModelName == "SCHISM_NETCDF_OUT") {
    double eTimeDay=MinimumTimeHistoryArray(TotalArr.eArr);
    Eigen::Tensor<double,3> znl=NETCDF_Get3DvariableSpecTime(TotalArr, "znl", eTimeDay);
    return znl;
  }
  if (eModelName == "ROMS") {
    int eta_rho=TotalArr.GrdArr.GrdArrRho.LON.rows();
    int xi_rho =TotalArr.GrdArr.GrdArrRho.LON.cols();
    MyMatrix<double> zeta=ZeroMatrix<double>(eta_rho, xi_rho);
    return ROMS_ComputeVerticalGlobalCoordinate(TotalArr.GrdArr, zeta);
  }
  std::cerr << "Error in RetrieveStandardVerticalCoordinate\n";
  std::cerr << "eModelName=" << eModelName << "\n";
  std::cerr << "Programmed models : SCHISM_NETCDF_OUT and ROMS\n";
  std::cerr << "We cannot find a matching model\n";
  throw TerminalException{1};
}

struct VerticalLevelInfo {
  int Choice;
  double dep;
  std::string strDepth;
  std::string type;
};


VerticalLevelInfo RetrieveVerticalInformation(std::string const& FullVarName, std::string const& eModelName)
{
  int Choice=-1;
  std::string strDepth="unset", type="unset";
  std::vector<std::string> ListStr=STRING_Split(FullVarName, ":");
  std::string eVarName=ListStr[0];
  if (ListStr.size() != 2) {
    std::cerr << "Error in the variable name. ListStr should be of length 2\n";
    std::cerr << "FullVarName=" << FullVarName << "\n";
    std::cerr << "eVarName=" << eVarName << "\n";
    std::cerr << "eModelName=" << eModelName << "\n";
    std::cerr << "|ListStr|=" << ListStr.size() << "\n";
    throw TerminalException{1};
  }
  std::string str=ListStr[1];
  std::vector<std::string> ListStrB=STRING_SplitCharNb(str);
  if (ListStrB.size() != 3) {
    std::cerr << "Error in the variable name should have 3 blocks\n";
    std::cerr << "|ListStrB|=" << ListStrB.size() << "\n";
    std::cerr << "str=" << str << "\n";
    int siz=ListStrB.size();
    for (int j=0; j<siz; j++)
      std::cerr << "j=" << j << " ListStrB[j]=" << ListStrB[j] << "\n";
    throw TerminalException{1};
  }
  if (ListStrB[0] == "VA")
    Choice=1;
  if (ListStrB[0] == "VR")
    Choice=2;
  if (Choice == -1) {
    std::cerr << "We should have VA or VR as possible choice\n";
    throw TerminalException{1};
  }
  double eVal;
  std::istringstream(ListStrB[1]) >> eVal;
  double dep=GetUnitInMeter(eVal,ListStrB[2]);
  strDepth=" at " + ListStr[1];
  //
  if (eModelName == "WWM") {
    std::string str=ListStr[1];
    std::vector<std::string> ListStrB=STRING_SplitCharNb(str);
    std::string type=ListStrB[0];
    strDepth=ListStrB[1];
  }
  //
  return {Choice, dep, strDepth, type};
}


std::string TransformVarName(std::string const& str)
{
  std::string RetStr;
  int len=str.size();
  for (int i=0; i<len; i++) {
    std::string eChar=str.substr(i,1);
    if (eChar == ":")
      eChar="_";
    RetStr += eChar;
  }
  return RetStr;
}


MyMatrix<double> GetNormMatrix(MyMatrix<double> const& U, MyMatrix<double> const& V)
{
  int nbRow=U.rows();
  int nbCol=U.cols();
  MyMatrix<double> Fwr(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      double eU=U(iRow,iCol);
      double eV=V(iRow,iCol);
      //	  std::cerr << "iRow=" << iRow << " iCol=" << iCol << "\n";
      //	  std::cerr << "eU=" << eU << " eV=" << eV << "\n";
      double eNorm=sqrt(eU*eU + eV*eV);
      //	  std::cerr << "eNorm=" << eNorm << "\n";
      Fwr(iRow,iCol) = eNorm;
    }
  return Fwr;
}





RecVar ModelSpecificVarSpecificTime_Kernel(TotalArrGetData const& TotalArr, std::string const& FullVarName, double const& eTimeDay)
{
  //  std::cerr << "TotalArr.GrdArr.ModelName = " << TotalArr.GrdArr.ModelName << "\n";
  std::string eModelName=GetBasicModelName(TotalArr.GrdArr.ModelName);
  //  std::cerr << "eModelName=" << eModelName << "\n";
  std::vector<std::string> ListStr=STRING_Split(FullVarName, ":");
  std::string eVarName=ListStr[0];
  //  std::cerr << "   ModelSpecificVarSpecificTime_Kernel, FullVarName=" << FullVarName << "\n";
  //  std::cerr << "   ModelSpecificVarSpecificTime_Kernel, eVarName=" << eVarName << " eTimeDay=" << eTimeDay << " eModelName=" << eModelName << "\n";
  int eta_rho=TotalArr.GrdArr.GrdArrRho.LON.rows();
  int xi_rho=TotalArr.GrdArr.GrdArrRho.LON.cols();
  std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
  std::string strFile=DATE_ConvertMjd2mystringFile(eTimeDay);
  RecSymbolic RecS;
  RecS.eTimeDay=eTimeDay;
  RecS.strPres="at " + strPres;
  RecS.strFile=strFile;
  RecS.VarNature="rho";
  RecS.FullVarName=FullVarName;
  RecS.VarName1=TransformVarName(FullVarName);
  RecS.VarName2="unset";
  RecS.strTime_ROMS="unset";
  RecS.varName_GRIB="unset";
  MyMatrix<double> F;
  MyMatrix<double> U;
  MyMatrix<double> V;
  Eigen::Tensor<double,3> Tens3;
  Eigen::Tensor<double,3> Uthree;
  Eigen::Tensor<double,3> Vthree;
  //
  // Generic model kind of variables
  //
  if (eVarName == "NbIterSolv") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "NB_ITER_SOLV", eTimeDay);
    RecS.VarName2="nb Iteration Solver";
    RecS.minval=0;
    RecS.maxval=50;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="nondim.";
  }
  if (eVarName == "CFL1") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "CFL1", eTimeDay);
    RecS.VarName2="CFL1";
    RecS.minval=0;
    RecS.maxval=5;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "CFL2") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "CFL2", eTimeDay);
    RecS.VarName2="CFL2";
    RecS.minval=0;
    RecS.maxval=5;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "CFL3") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "CFL3", eTimeDay);
    RecS.VarName2="CFL3";
    RecS.minval=0;
    RecS.maxval=5;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "ThreeDfield1") {
    if (eModelName == "WWM" || eModelName == "ROMS")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ThreeDfield1", eTimeDay);
    RecS.VarName2="Generic three dim. field 1";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.VarNature="3Drho";
    RecS.Unit="unspecified";
  }
  if (eVarName == "CGwave") {
    if (eModelName == "WWM")
      Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "CG", eTimeDay);
    RecS.VarName2="group velocity";
    RecS.minval=0;
    RecS.maxval=30;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.VarNature="3Drho";
    RecS.Unit="m/s";
  }
  if (eVarName == "FieldOut1") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "FieldOut1", eTimeDay);
    RecS.VarName2="Generic Field Out 1";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="unspecified";
  }
  if (eVarName == "ChlorophyllConcOCI") {
    RecS.VarName2="Chlorophyll Concentration, OCI Algorithm";
    RecS.CFshortName="chlor_a";
    RecS.minval=0;
    RecS.maxval=100;
    RecS.mindiff=-50;
    RecS.maxdiff=50;
    RecS.Unit="mg m-3";
  }
  if (eVarName == "ChlorophyllConcOCX") {
    RecS.VarName2="Chlorophyll Concentration OCX Algorithm";
    RecS.CFshortName="chl_ocx";
    RecS.minval=0;
    RecS.maxval=100;
    RecS.mindiff=-50;
    RecS.maxdiff=50;
    RecS.Unit="mg m-3";
  }
  if (eVarName == "CalciteConc") {
    RecS.VarName2="Calcite Concentration";
    RecS.CFshortName="pic";
    RecS.minval=0;
    RecS.maxval=0.125;
    RecS.mindiff=-0.03;
    RecS.maxdiff=0.03;
    RecS.Unit="mol m-3";
  }
  if (eVarName == "ParticleOrganicCarbon") {
    RecS.VarName2="Particle Organic Carbon";
    RecS.CFshortName="poc";
    RecS.minval=0;
    RecS.maxval=1000;
    RecS.mindiff=-300;
    RecS.maxdiff=300;
    RecS.Unit="mg m-3";
  }
  if (eVarName == "InstPhotosyntheticallyAvailableRad") {
    RecS.VarName2="Instantaneous Photosynthetically Available Radiation";
    RecS.CFshortName="ipar";
    RecS.minval=0;
    RecS.maxval=0.0032;
    RecS.mindiff=-0.001;
    RecS.maxdiff= 0.001;
    RecS.Unit="einstein m-2 s-1";
  }
  if (eVarName == "PhotosyntheticallyAvailableRad") {
    RecS.VarName2="Photosynthetically Available Radiation";
    RecS.CFshortName="par";
    RecS.minval=0;
    RecS.maxval=130;
    RecS.mindiff=-50;
    RecS.maxdiff=50;
    RecS.Unit="einstein m-2 day-1";
  }
  if (eVarName == "ChlorophylA") {
    if (eModelName == "CFCONVENTION")
      F=Get2DvariableSpecTime(TotalArr, "FieldOut1", eTimeDay);
    RecS.VarName2="Chlorophyl A concentration";
    RecS.minval=0;
    RecS.maxval=300;
    RecS.mindiff=0;
    RecS.maxdiff=300;
    RecS.Unit="mg m-3";
  }
  if (eVarName == "IOBP") {
    if (eModelName == "WWM") // we should have a NETCDF_WW3
      F=Get2DvariableSpecTime(TotalArr, "IOBP_WW3", eTimeDay);
    if (eModelName == "UNRUNOFF") {
      int mnp=TotalArr.GrdArr.IOBP.size();
      if (mnp == 0) {
	std::cerr << "The IOBP has not been assigned originally. So cannot plot IOBP\n";
	std::cerr << "Most likely the BoundFile is set to \"unset\" and that is wrong\n";
	throw TerminalException{1};
      }
      MyVector<double> Fret(mnp);
      for (int i=0; i<mnp; i++)
	Fret(i) = double(TotalArr.GrdArr.IOBP(i));
      F=Fret;
    }
    RecS.VarName2="IOBP of the unstructured model";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "RoughnessFactor") {
    if (eModelName == "WWM") {
      F=GetRoughnessFactor(TotalArr.GrdArr.GrdArrRho.DEP, TotalArr.GrdArr.INE);
    }
    RecS.VarName2="Roughness factor";
    RecS.minval=0;
    RecS.maxval=0.2;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "MAPSTA") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "MAPSTA", eTimeDay);
    RecS.VarName2="MAPSTA of wavewatchIII";
    RecS.minval=-2;
    RecS.maxval=2;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="nondim.";
  }
  //
  // Atmospheric variables
  //
  if (eVarName == "Uwind") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
    F=RecVarWork.U;
    RecS.VarName2="Eastward wind";
    RecS.minval=-10;
    RecS.maxval=10;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="m/s";
  }
  if (eVarName == "Vwind") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
    F=RecVarWork.V;
    RecS.VarName2="Northward wind";
    RecS.minval=-10;
    RecS.maxval=10;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="m/s";
  }
  if (eVarName == "WIND10") {
    if (eModelName == "UNRUNOFF") {
      U=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "SCHISM_SFLUX") {
      U=Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      U=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "ROMS" || eModelName == "WWM") {
      U=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U=Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
    }
    if (eModelName == "WRF") {
      U=Get2DvariableSpecTime(TotalArr, "U10", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "V10", eTimeDay);
    }
    if (eModelName == "WW3") {
      Eigen::Tensor<double,3> Utens = NETCDF_Get3DvariableSpecTime(TotalArr, "u10m", eTimeDay);
      Eigen::Tensor<double,3> Vtens = NETCDF_Get3DvariableSpecTime(TotalArr, "v10m", eTimeDay);
      U = DimensionExtraction(Utens, 0, 0);
      V = DimensionExtraction(Vtens, 0, 0);
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
      U=Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2="10m wind";
    RecS.minval=0;
    RecS.maxval=13;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="m/s";
    RecS.VarNature="uv";
    RecS.nameU="Uwind";
    RecS.nameV="Vwind";
    RecS.strTime_ROMS="wind_time";
    RecS.varName_ROMS_U="Uwind";
    RecS.varName_ROMS_V="Vwind";
  }
  if (eVarName == "WINDMAG") {
    if (eModelName == "ROMS" || eModelName == "WWM") {
      if (TOTALARR_IsVar(TotalArr, "Uwind") && TOTALARR_IsVar(TotalArr, "Vwind") ) {
	MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
	MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
	F=COMPUTE_NORM(Us, Vs);
      }
      else {
	if (eModelName == "WWM")
	  F=Get2DvariableSpecTime(TotalArr, "WINDMAG", eTimeDay);
	else
	  F=Get2DvariableSpecTime(TotalArr, "WNDMAG", eTimeDay);
      }
    }
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
      F=COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "SCHISM_SFLUX") {
      MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
      MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
      F=COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
      F=COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
      F=COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
      MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
      F=COMPUTE_NORM(Us, Vs);
    }
    if (eModelName == "GRIB_WAM_FORT30")
      F=Get2DvariableSpecTime(TotalArr, "wind", eTimeDay);
    RecS.VarName2="10m wind speed";
    RecS.minval=0;
    RecS.maxval=13;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="m/s";
  }
  if (eVarName == "AIRD") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "AIRD", eTimeDay);
    RecS.VarName2="air density";
    RecS.minval=1.12;
    RecS.maxval=1.20;
    RecS.mindiff=-0.02;
    RecS.maxdiff=0.02;
    RecS.Unit="kg/m3";
  }
  if (eVarName == "rain") {
    if (eModelName == "UNRUNOFF")
      F=Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
    if (eModelName == "SCHISM_SFLUX")
      F=Get2DvariableSpecTime(TotalArr, "prate", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
    //    if (eModelName == "GRIB_DWD")
    //      F=Get2DvariableSpecTime(TotalArr, "tp", eTimeDay);
    if (eModelName == "WRF")
      F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "RAINNC", eTimeDay);
    std::vector<std::string> ListGRIBmodel{"GRIB_COSMO", "GRIB_DWD", "GRIB_ECMWF", "GRIB_ALADIN"};
    if (PositionVect(ListGRIBmodel, eModelName) != -1)
      F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "tp", eTimeDay);
    int siz=F.size();
    for (int u=0; u<siz; u++)
      F(u) = std::max(F(u), double(0));
    RecS.VarName2="rainfall rate";
    RecS.minval=0;
    RecS.maxval=0.001;
    RecS.mindiff=-0.001;
    RecS.maxdiff=0.001;
    RecS.Unit="kg/m^2/s";
    RecS.strTime_ROMS="rain_time";
    RecS.varName_ROMS="rain";
  }
  if (eVarName == "swrad") {
    if (eModelName == "SCHISM_SFLUX")
      F=Get2DvariableSpecTime(TotalArr, "dswrf", eTimeDay);
    if (eModelName == "WRF")
      F=Get2DvariableSpecTime(TotalArr, "SWDOWN", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "swrad", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "srad", eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F=Get2DvariableSpecTime(TotalArr, "sobs_rad", eTimeDay);
    if (eModelName == "GRIB_ALADIN")
      F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "nswrs", eTimeDay);
    if (eModelName == "GRIB_ECMWF")
      F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "ssrd", eTimeDay);
    //      F=Get2DvariableSpecTime(TotalArr, "ssrd", eTimeDay);
    RecS.VarName2="Shortwave flux";
    RecS.minval=100;
    RecS.maxval=1000;
    RecS.mindiff=-100;
    RecS.maxdiff=100;
    RecS.Unit="W/m2";
    RecS.strTime_ROMS="srf_time";
    RecS.varName_ROMS="swrad";
  }
  if (eVarName == "lwrad") {
    if (eModelName == "SCHISM_SFLUX")
      F=Get2DvariableSpecTime(TotalArr, "dlwrf", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT") {
      MyMatrix<double> F1=Get2DvariableSpecTime(TotalArr, "hradu", eTimeDay);
      MyMatrix<double> F2=Get2DvariableSpecTime(TotalArr, "hradd", eTimeDay);
      F = F1 - F2;
    }
    if (eModelName == "WRF")
      F=Get2DvariableSpecTime(TotalArr, "GLW", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "lwrad", eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F=Get2DvariableSpecTime(TotalArr, "thbs_rad", eTimeDay);
    if (eModelName == "GRIB_ECMWF")
      F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "strd", eTimeDay);
    //      F=Get2DvariableSpecTime(TotalArr, "thbs_rad", eTimeDay);
    RecS.VarName2="Longwave flux";
    RecS.minval=200;
    RecS.maxval=500;
    RecS.mindiff=-50;
    RecS.maxdiff= 50;
    RecS.Unit="W/m2";
    RecS.strTime_ROMS="lrf_time";
    RecS.varName_ROMS="lwrad_down";
  }
  if (eVarName == "latent") {
    if (eModelName == "UNRUNOFF")
      F=Get2DvariableSpecTime(TotalArr, "lat_flux", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "latent", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "fluxlu", eTimeDay);
    RecS.VarName2="Latent flux";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="W/m2";
  }
  if (eVarName == "SurfPres") {
    if (eModelName == "ROMS") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
      F=100*Fin;
    }
    if (eModelName == "UNRUNOFF")
      F=Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
    if (eModelName == "WRF")
      F=Get2DvariableSpecTime(TotalArr, "PSFC", eTimeDay);
    if (eModelName == "SCHISM_SFLUX")
      F=Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "pr", eTimeDay);
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS")
      F=Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
    if (eModelName == "GRIB_ECMWF" || eModelName == "GRIB_ALADIN")
      F=Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
    if (eModelName == "GRIB_COSMO") {
      if (TOTALARR_IsVar(TotalArr, "pmsl"))
	F=Get2DvariableSpecTime(TotalArr, "pmsl", eTimeDay);
      if (TOTALARR_IsVar(TotalArr, "msl"))
	F=Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
    }
    RecS.VarName2="mean sea level pressure";
    RecS.minval=100000;
    RecS.maxval=103000;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="Pa";
    RecS.strTime_ROMS="pair_time";
    RecS.varName_ROMS="Pair";
  }
  if (eVarName == "sensible") {
    if (eModelName == "UNRUNOFF")
      F=Get2DvariableSpecTime(TotalArr, "sen_flux", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "sensible", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "fluxsu", eTimeDay);
    RecS.VarName2="Sensible heat flux";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="W/m2";
  }
  if (eVarName == "shflux") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "shflux", eTimeDay);
    RecS.VarName2="Surface heat flux";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="W/m2";
  }
  if (eVarName == "ssflux") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ssflux", eTimeDay);
    RecS.VarName2="Surface salinity flux";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="PSU/m2s";
  }
  if (eVarName == "evaporation") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "evaporation", eTimeDay);
    RecS.VarName2="Evaporation rate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="kg/m2s";
  }
  if (eVarName == "CloudFraction") {
    if (eModelName == "ROMS")
      F = Get2DvariableSpecTime(TotalArr, "cloud", eTimeDay);
    if (eModelName == "GRIB_ALADIN")
      F = Get2DvariableSpecTime(TotalArr, "tcc", eTimeDay);
    RecS.VarName2="Cloud fraction";
    RecS.minval=0;
    RecS.maxval=100;
    RecS.mindiff=-20;
    RecS.maxdiff=20;
    RecS.Unit="nondimensional";
    RecS.strTime_ROMS="cloud_time";
    RecS.varName_ROMS="cloud";
  }
  if (eVarName == "AIRT2K") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "AIRT2", eTimeDay);
    F=RecVarWork.F;
    int siz=F.size();
    for (int i=0; i<siz; i++)
      F(i) += double(273.15);
    RecS.VarName2="2m air temperature (K)";
    RecS.minval=273.15 + 10;
    RecS.maxval=273.15 + 20;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="deg";
  }
  if (eVarName == "AIRT2") {
    if (eModelName == "SCHISM_SFLUX") {
      F=Get2DvariableSpecTime(TotalArr, "stmp", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      F=Get2DvariableSpecTime(TotalArr, "airt1", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    if (eModelName == "ROMS") {
      F=Get2DvariableSpecTime(TotalArr, "Tair", eTimeDay);
    }
    if (eModelName == "COSMO") {
      F=Get2DvariableSpecTime(TotalArr, "t_2m", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    if (eModelName == "WRF") {
      F=Get2DvariableSpecTime(TotalArr, "T2", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    std::vector<std::string> ListModel{"GRIB_DWD", "GRIB_ECMWF", "GRIB_GFS", "GRIB_COSMO", "GRIB_ALADIN"};
    if (PositionVect(ListModel, eModelName) != -1) {
      F=Get2DvariableSpecTime(TotalArr, "2t", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    RecS.VarName2="2m air temperature";
    RecS.minval=10;
    RecS.maxval=20;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="deg";
    RecS.strTime_ROMS="tair_time";
    RecS.varName_ROMS="Tair";
  }
  if (eVarName == "Rh2frac") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "Rh2", eTimeDay);
    F=RecVarWork.F / double(100);
    RecS.VarName2="2m relative humidity";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-0.2;
    RecS.maxdiff=0.2;
    RecS.Unit="nondim.";
  }
  if (eVarName == "Rh2") {
    if (eModelName == "COSMO")
      F=Get2DvariableSpecTime(TotalArr, "rh_2m", eTimeDay);
    if (eModelName == "GRIB_ALADIN")
      F=Get2DvariableSpecTime(TotalArr, "r", eTimeDay);
    if (eModelName == "GRIB_DWD")
      F=Get2DvariableSpecTime(TotalArr, "RELHUM_2M", eTimeDay);
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "shum1", eTimeDay);
    if (eModelName == "SCHISM_SFLUX") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "spfh", eTimeDay);
      F=100 * Fin;
    }
    if (eModelName == "ROMS") {
      F=Get2DvariableSpecTime(TotalArr, "Qair", eTimeDay);
    }
    if (eModelName == "GRIB_ECMWF") {
      F=ConvertSpecHumid_to_RelativeHumidity(TotalArr, "2r", "2t", "msl", "q", eTimeDay);
    }
    if (eModelName == "GRIB_COSMO") {
      F=ConvertSpecHumid_to_RelativeHumidity(TotalArr, "2r", "2t", "msl", "QV_S", eTimeDay);
    }
    if (eModelName == "WRF") {
      MyMatrix<double> Q2=Get2DvariableSpecTime(TotalArr, "Q2", eTimeDay);
      MyMatrix<double> PSFC=Get2DvariableSpecTime(TotalArr, "PSFC", eTimeDay);
      MyMatrix<double> T2K=Get2DvariableSpecTime(TotalArr, "T2", eTimeDay);
      F=mixing_ratio2relative_humidity(Q2, PSFC, T2K);
    }
    RecS.VarName2="2m relative humidity";
    RecS.minval=0;
    RecS.maxval=100;
    RecS.mindiff=-20;
    RecS.maxdiff=20;
    RecS.Unit="nondim.";
    RecS.strTime_ROMS="qair_time";
    RecS.varName_ROMS="Qair";
  }
  //
  // Oceanic variables
  //
  if (eVarName == "HorizCurr") {
    VerticalLevelInfo VertInfo;
    if (eModelName != "TRIVIAL")
      VertInfo = RetrieveVerticalInformation(FullVarName, eModelName);
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      Uthree=My_u2rho_3D(Utot, TotalArr.GrdArr.GrdArrU.MSK);
      Vthree=My_v2rho_3D(Vtot, TotalArr.GrdArr.GrdArrV.MSK);
      MyMatrix<double> zeta=Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
      ARVDtyp ARVD=TOTALARR_GetARVD(TotalArr);
      U=VerticalInterpolation_P2_R(ARVD, TotalArr.GrdArr.GrdArrRho.DEP, zeta, TotalArr.GrdArr.GrdArrRho.MSK, VertInfo.dep, Uthree, VertInfo.Choice);
      V=VerticalInterpolation_P2_R(ARVD, TotalArr.GrdArr.GrdArrRho.DEP, zeta, TotalArr.GrdArr.GrdArrRho.MSK, VertInfo.dep, Vthree, VertInfo.Choice);
      AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
      Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
      Eigen::Tensor<double,3> znl=NETCDF_Get3DvariableSpecTime(TotalArr, "znl", eTimeDay);
      MyMatrix<double> zeta=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
      U=VerticalInterpolation_SCHISM_ZNL(znl, zeta, VertInfo.dep, Utot, VertInfo.Choice);
      V=VerticalInterpolation_SCHISM_ZNL(znl, zeta, VertInfo.dep, Vtot, VertInfo.Choice);
    }
    if (eModelName == "WWM") {
      std::string strCallU="UhorizCurr_" + VertInfo.type + VertInfo.strDepth + "m";
      std::string strCallV="VhorizCurr_" + VertInfo.type + VertInfo.strDepth + "m";
      U=Get2DvariableSpecTime(TotalArr, strCallU, eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, strCallV, eTimeDay);
    }
    RecS.VarName2="horizontal current" + VertInfo.strDepth;
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
    RecS.VarNature="uv";
    RecS.nameU="UhorizCurr";
    RecS.nameV="VhorizCurr";
  }
  if (eVarName == "HorizTemp") {
    VerticalLevelInfo VertInfo;
    if (eModelName != "TRIVIAL")
      VertInfo = RetrieveVerticalInformation(FullVarName, eModelName);
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> TEMPtot=NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
      MyMatrix<double> zeta=Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
      ARVDtyp ARVD=TOTALARR_GetARVD(TotalArr);
      F=VerticalInterpolation_P2_R(ARVD, TotalArr.GrdArr.GrdArrRho.DEP, zeta, TotalArr.GrdArr.GrdArrRho.MSK, VertInfo.dep, TEMPtot, VertInfo.Choice);
    }
    RecS.VarName2="horizontal temperature" + VertInfo.strDepth;
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
  }
  if (eVarName == "UsurfCurr") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
    F=RecVarWork.U;
    RecS.VarName2="Eastward current";
    RecS.minval=-0.3;
    RecS.maxval=0.3;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
  }
  if (eVarName == "VsurfCurr") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
    F=RecVarWork.V;
    RecS.VarName2="Northward current";
    RecS.minval=-0.3;
    RecS.maxval=0.3;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
  }
  if (eVarName == "Curr") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      Uthree=My_u2rho_3D(Utot, TotalArr.GrdArr.GrdArrU.MSK);
      Vthree=My_v2rho_3D(Vtot, TotalArr.GrdArr.GrdArrV.MSK);
    }
    if (eModelName == "WWM") {
      Uthree=NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
      Vthree=NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
    }
    if (eModelName == "NEMO") {
      Uthree=NEMO_Get3DvariableSpecTime(TotalArr, "cur", "uo", eTimeDay);
      Vthree=NEMO_Get3DvariableSpecTime(TotalArr, "cur", "vo", eTimeDay);
    }
    RecS.VarName2="baroclinic current";
    RecS.minval=0;
    RecS.maxval=0.2;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.VarNature="3Duv";
    RecS.Unit="m/s";
  }
  if (eVarName == "Temp") {
    if (eModelName == "ROMS")
      Tens3=NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
    if (eModelName == "NEMO")
      Tens3=NEMO_Get3DvariableSpecTime(TotalArr, "tem", "thetao", eTimeDay);
    RecS.VarName2="temperature";
    RecS.minval=17;
    RecS.maxval=25;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.VarNature="3Drho";
    RecS.Unit="deg C";
  }
  if (eVarName == "Salt") {
    if (eModelName == "ROMS")
      Tens3=NETCDF_Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
    if (eModelName == "NEMO")
      Tens3=NEMO_Get3DvariableSpecTime(TotalArr, "sal", "so", eTimeDay);
    RecS.VarName2="salinity";
    RecS.minval=35;
    RecS.maxval=38;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.VarNature="3Drho";
    RecS.Unit="PSU";
  }
  if (eVarName == "CurrMag") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "Curr", eTimeDay);
    auto LDim=RecVarWork.Uthree.dimensions();
    int dim0=LDim[0];
    int dim1=LDim[1];
    int dim2=LDim[2];
    Eigen::Tensor<double,3> Tens3(dim0, dim1, dim2);
    for (int i0=0; i0<dim0; i0++)
      for (int i1=0; i1<dim1; i1++)
	for (int i2=0; i2<dim2; i2++) {
	  double eU=RecVarWork.Uthree(i0, i1, i2);
	  double eV=RecVarWork.Vthree(i0, i1, i2);
	  double eNorm=sqrt(eU*eU + eV*eV);
	  Tens3(i0, i1, i2) = eNorm;
	}
    RecS.VarName2="baroclinic current magnitude";
    RecS.minval=0;
    RecS.maxval=0.2;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.VarNature="3Drho";
    RecS.Unit="m/s";
  }
  if (eVarName == "BottomCurr") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      MyMatrix<double> Usurf=DimensionExtraction(Utot, 0, 0);
      MyMatrix<double> Vsurf=DimensionExtraction(Vtot, 0, 0);
      U=My_u2rho(Usurf, TotalArr.GrdArr.GrdArrU.MSK);
      V=My_v2rho(Vsurf, TotalArr.GrdArr.GrdArrV.MSK);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2="bottom current";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
    RecS.VarNature="uv";
    RecS.nameU="UbottomCurr";
    RecS.nameV="VbottomCurr";
  }
  if (eVarName == "CurrBaro") {
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> Ucurr=Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
      MyMatrix<double> Vcurr=Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
      MyMatrix<double> Helev = Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
      U = Helev.cwiseProduct(Ucurr);
      V = Helev.cwiseProduct(Vcurr);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2="barotropic current";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m^2/s";
    RecS.VarNature="uv";
    RecS.nameU="UsurfCurr";
    RecS.nameV="VsurfCurr";
  }
  if (eVarName == "CurrBaroMag") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "CurrBaro", eTimeDay);
    F=COMPUTE_NORM(RecVarWork.U, RecVarWork.V);
    RecS.VarName2="barotropic current magnitude";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m^2/s";
  }
  if (eVarName == "SurfCurr") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      int s_rho=Utot.dimension(0);
      MyMatrix<double> Usurf=DimensionExtraction(Utot, 0, s_rho-1);
      MyMatrix<double> Vsurf=DimensionExtraction(Vtot, 0, s_rho-1);
      U=My_u2rho(Usurf, TotalArr.GrdArr.GrdArrU.MSK);
      V=My_v2rho(Vsurf, TotalArr.GrdArr.GrdArrV.MSK);
    }
    if (eModelName == "UNRUNOFF") {
      U=Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
    }
    if (eModelName == "WWM") {
      if (TOTALARR_IsVar(TotalArr, "CURTX")  && TOTALARR_IsVar(TotalArr, "CURTX") ) {
	U=Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
	V=Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
      }
      else {
	if (TOTALARR_IsVar(TotalArr, "UsurfCurr")  && TOTALARR_IsVar(TotalArr, "VsurfCurr") ) {
	  U=Get2DvariableSpecTime(TotalArr, "UsurfCurr", eTimeDay);
	  V=Get2DvariableSpecTime(TotalArr, "VsurfCurr", eTimeDay);
	}
	else {
	  Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
	  Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
	  int s_rho=Utot.dimension(0);
	  U=DimensionExtraction(Utot, 0, s_rho-1);
	  V=DimensionExtraction(Vtot, 0, s_rho-1);
	}
      }
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U=Get2DvariableSpecTime(TotalArr, "ucurr", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "vcurr", eTimeDay);
    }
    if (eModelName == "WW3") {
      U=Get2DvariableSpecTime(TotalArr, "U", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "V", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2="surface current";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
    RecS.VarNature="uv";
    RecS.nameU="UsurfCurr";
    RecS.nameV="VsurfCurr";
  }
  if (eVarName == "SurfStress") {
    if (eModelName == "UNRUNOFF") {
      U=Get2DvariableSpecTime(TotalArr, "sustr", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "svstr", eTimeDay);
    }
    if (eModelName == "ROMS") {
      U=Get2DvariableSpecTime(TotalArr, "sustr", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "svstr", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr.GrdArrRho.ANG);
    RecS.VarName2="surface stress";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
    RecS.VarNature="uv";
    RecS.nameU="UsurfStress";
    RecS.nameV="VsurfStress";
  }
  if (eVarName == "SurfCurrMag") {
    RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
    F=COMPUTE_NORM(RecVarWork.U, RecVarWork.V);
    RecS.VarName2="surface current magnitude";
    RecS.minval=0;
    RecS.maxval=0.5;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m/s";
  }
  if (eVarName == "TempSurf") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> TheTemp=NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
      int s_rho=TheTemp.dimension(0);
      F=DimensionExtraction(TheTemp, 0, s_rho-1);
    }
    if (eModelName == "COSMO") {
      F=Get2DvariableSpecTime(TotalArr, "t_s", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      Eigen::Tensor<double,3> TheTemp=NETCDF_Get3DvariableSpecTime(TotalArr, "tr_nd1", eTimeDay);
      int s_rho=TheTemp.dimension(0);
      F=DimensionExtraction(TheTemp, 0, s_rho-1);
    }
    RecS.VarName2="sea surface temperature";
    RecS.minval=10;
    RecS.maxval=20;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="deg";
  }
  if (eVarName == "TempBottom") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> TheTemp=NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
      F=DimensionExtraction(TheTemp, 0, 0);
    }
    RecS.VarName2="sea bottom temperature";
    RecS.minval=10;
    RecS.maxval=20;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="deg";
  }
  if (eVarName == "SaltSurf") {
    if (eModelName == "ROMS") {
      Eigen::Tensor<double,3> TheSalt=NETCDF_Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
      int s_rho=TheSalt.dimension(0);
      F=DimensionExtraction(TheSalt, 0, s_rho-1);
    }
    if (eModelName == "SCHISM_NETCDF_OUT") {
      Eigen::Tensor<double,3> TheTemp=NETCDF_Get3DvariableSpecTime(TotalArr, "tr_nd2", eTimeDay);
      int s_rho=TheTemp.dimension(0);
      F=DimensionExtraction(TheTemp, 0, s_rho-1);
    }
    RecS.VarName2="sea surface salinity";
    RecS.minval=30;
    RecS.maxval=40;
    RecS.mindiff=-2;
    RecS.maxdiff=2;
    RecS.Unit="PSU";
  }
  if (eVarName == "ZetaOcean") {
    if (eModelName == "COSMO")
      F=Get2DvariableSpecTime(TotalArr, "ZetaOcean", eTimeDay);
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
    if (eModelName == "WW3")
      F=Get2DvariableSpecTime(TotalArr, "XE", eTimeDay);
    if (eModelName == "UNRUNOFF") {
      MyMatrix<double> TotalElev = Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
      MyMatrix<double> DEP=TotalArr.GrdArr.GrdArrRho.DEP;
      if (!IsEqualSizeMatrices(TotalElev, DEP)) {
	std::cerr << "The matrices TotalElev and DEP have different sizes\n";
	std::cerr << "Most likely the grid does not match the history used\n";
	throw TerminalException{1};
      }
      F=TotalElev + DEP;
    }
    if (eModelName == "SCHISM_NETCDF_OUT")
      F=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
    if (eModelName == "NEMO")
      F=NEMO_Get2DvariableSpecTime(TotalArr, "ssh", "zos", eTimeDay);
    RecS.VarName2="free surface elevation";
    RecS.minval=-0.2;
    RecS.maxval=0.2;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m";
  }
  if (eVarName == "ZetaOceanDerivative") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "DEPDT", eTimeDay);
    RecS.VarName2="free surface elevation derivative";
    RecS.minval=-0.01;
    RecS.maxval=0.01;
    RecS.mindiff=-0.001;
    RecS.maxdiff=0.001;
    RecS.Unit="m/s";
  }
  //
  // Wave variables
  //
  if (eVarName == "MeanWaveLength") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "WLM", eTimeDay);
    RecS.VarName2="mean wave length";
    RecS.minval=2;
    RecS.maxval=30;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="m";
  }
  if (eVarName == "PeakWaveLength") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "LPP", eTimeDay);
    RecS.VarName2="peak wave length";
    RecS.minval=2;
    RecS.maxval=30;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="m";
  }
  if (eVarName == "MeanWaveNumber") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "KLM", eTimeDay);
    RecS.VarName2="mean wave number";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m-1";
  }
  if (eVarName == "PeakWaveNumber") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "KPP", eTimeDay);
    RecS.VarName2="peak wave number";
    RecS.minval=0;
    RecS.maxval=1;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m-1";
  }
  if (eVarName == "MeanWaveDir") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
    RecS.VarName2="mean wave direction";
    RecS.minval=0;
    RecS.maxval=360;
    RecS.mindiff=-30;
    RecS.maxdiff=30;
    RecS.Unit="deg";
  }
  if (eVarName == "PeakWaveDir") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
    RecS.VarName2="peak wave direction";
    RecS.minval=0;
    RecS.maxval=360;
    RecS.mindiff=-30;
    RecS.maxdiff=30;
    RecS.Unit="deg";
  }
  if (eVarName == "MeanWaveDirVect") {
    if (eModelName == "WWM") {
      F=Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
      int nbRow=F.rows();
      int nbCol=F.cols();
      double deg2rad=3.1415926535 / double(180);
      U=MyMatrix<double>(nbRow,nbCol);
      V=MyMatrix<double>(nbRow,nbCol);
      for (int iRow=0; iRow<nbRow; iRow++)
	for (int iCol=0; iCol<nbCol; iCol++) {
	  double eAngRad=deg2rad*F(iRow,iCol);
	  U(iRow,iCol)=cos(eAngRad);
	  V(iRow,iCol)=sin(eAngRad);
	}
    }
    RecS.VarName2="mean wave direction";
    RecS.minval=0;
    RecS.maxval=360;
    RecS.mindiff=-30;
    RecS.maxdiff=30;
    RecS.VarNature="uv";
    RecS.Unit="deg";
  }
  if (eVarName == "PeakWaveDirVect") {
    if (eModelName == "WWM") {
      F=Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
      int nbRow=F.rows();
      int nbCol=F.cols();
      double deg2rad=3.1415926535 / double(180);
      U=MyMatrix<double>(nbRow,nbCol);
      V=MyMatrix<double>(nbRow,nbCol);
      for (int iRow=0; iRow<nbRow; iRow++)
	for (int iCol=0; iCol<nbCol; iCol++) {
	  double eAngRad=deg2rad*F(iRow,iCol);
	  U(iRow,iCol)=cos(eAngRad);
	  V(iRow,iCol)=sin(eAngRad);
	}
    }
    RecS.VarName2="peak wave direction";
    RecS.minval=0;
    RecS.maxval=360;
    RecS.mindiff=-30;
    RecS.maxdiff=30;
    RecS.VarNature="uv";
    RecS.Unit="deg";
  }
  if (eVarName == "DiscPeakWaveDir") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "DPEAK", eTimeDay);
    RecS.VarName2="discrete peak wave direction";
    RecS.minval=0;
    RecS.maxval=360;
    RecS.mindiff=-30;
    RecS.maxdiff=30;
    RecS.Unit="deg";
  }
  if (eVarName == "ZetaSetup") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "ZETA_SETUP", eTimeDay);
    RecS.VarName2="free surface setup";
    RecS.minval=0;
    RecS.maxval=0.76;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="m";
  }
  if (eVarName == "BreakingFraction") {
    MyMatrix<double> Fhs, Fzeta;
    if (eModelName == "WWM")
      Fhs=Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
    if (eModelName == "WWM")
      Fzeta=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
    if (!IsEqualSizeMatrices(Fhs, TotalArr.GrdArr.GrdArrRho.DEP)) {
      std::cerr << "The matrices Fhs and DEP have different sizes\n";
      std::cerr << "Most likely the grid does not match the history used\n";
      throw TerminalException{1};
    }
    F=MyMatrix<double>(eta_rho, xi_rho);
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++)
	F(i,j)=Fhs(i,j) / (Fzeta(i,j) + TotalArr.GrdArr.GrdArrRho.DEP(i,j));
    RecS.VarName2="Breaking fraction";
    RecS.minval=0;
    RecS.maxval=0.76;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "Hwave") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
    if (eModelName == "WW3")
      F=Get2DvariableSpecTime(TotalArr, "hs", eTimeDay);
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "Hwave", eTimeDay);
    if (eModelName == "GRIB_WAM_FORT30" || eModelName == "GRIB_IFS")
      F=Get2DvariableSpecTime(TotalArr, "swh", eTimeDay);
    RecS.VarName2="Significant wave height";
    RecS.varName_GRIB="swh";
    RecS.minval=0;
    RecS.maxval=4.5;
    RecS.mindiff=-0.5;
    RecS.maxdiff=0.5;
    RecS.Unit="m";
  }
  if (eVarName == "MeanWaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
    if (eModelName == "WWM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    RecS.VarName2="mean wave frequency";
    RecS.minval=0;
    RecS.maxval=0.9;
    RecS.mindiff=-0.2;
    RecS.maxdiff=0.2;
    RecS.Unit="Hz";
  }
  if (eVarName == "PeakWaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
    if (eModelName == "WWM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    RecS.VarName2="peak wave frequency";
    RecS.minval=0;
    RecS.maxval=0.9;
    RecS.mindiff=-0.2;
    RecS.maxdiff=0.2;
    RecS.Unit="Hz";
  }
  if (eVarName == "MeanWavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
    RecS.VarName2="mean wave period";
    RecS.minval=2;
    RecS.maxval=10;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="s";
  }
  if (eVarName == "PeakWavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
    RecS.VarName2="peak wave period";
    RecS.minval=2;
    RecS.maxval=10;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="s";
  }
  if (eVarName == "TM02") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TM02", eTimeDay);
    RecS.VarName2="zero crossing wave period";
    RecS.minval=2;
    RecS.maxval=10;
    RecS.mindiff=-1;
    RecS.maxdiff=1;
    RecS.Unit="s";
  }
  if (eVarName == "DynBathy") {
    if (eModelName == "UNRUNOFF")
      F = Get2DvariableSpecTime(TotalArr, "H", eTimeDay);
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "DW", eTimeDay);
    RecS.VarName2="dynamic bathymetry";
    RecS.minval=0;
    RecS.maxval=30;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="deg";
  }
  if (eVarName == "Bathymetry") {
    if (eModelName == "WWM")
      F=TotalArr.GrdArr.GrdArrRho.DEP;
    PairMinMax ePair=ComputeMinMax(TotalArr.GrdArr, TotalArr.GrdArr.GrdArrRho.DEP);
    RecS.VarName2="bathymetry";
    RecS.minval=0;
    RecS.maxval=ePair.TheMax;
    RecS.mindiff=-50;
    RecS.maxdiff=50;
    RecS.Unit="m";
  }
  if (eVarName == "MeanWaveDirSpread") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "DSPR", eTimeDay);
    RecS.VarName2="directional spreading";
    RecS.minval=0;
    RecS.maxval=30;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="deg";
  }
  if (eVarName == "PeakWaveDirSpread") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "PEAKDSPR", eTimeDay);
    RecS.VarName2="peak directional spreading";
    RecS.minval=0;
    RecS.maxval=30;
    RecS.mindiff=-5;
    RecS.maxdiff=5;
    RecS.Unit="deg";
  }
  if (eVarName == "AirZ0") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "Z0", eTimeDay);
    RecS.VarName2="air roughness length";
    RecS.minval=0;
    RecS.maxval=  0.0002;
    RecS.mindiff=-0.00005;
    RecS.maxdiff= 0.00005;
    RecS.Unit="m";
  }
  if (eVarName == "AirFricVel") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "UFRIC", eTimeDay);
    RecS.VarName2="air roughness length";
    RecS.minval=0;
    RecS.maxval=0.3;
    RecS.mindiff=-0.05;
    RecS.maxdiff= 0.05;
    RecS.Unit="m";
  }
  if (eVarName == "CdWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "CdWave", eTimeDay);
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "CD", eTimeDay);
    RecS.VarName2="drag coefficient from the wave model";
    RecS.minval=0;
    RecS.maxval=0.20;
    RecS.mindiff=-0.05;
    RecS.maxdiff=0.05;
    RecS.Unit="nondim.";
  }
  if (eVarName == "AlphaWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "AlphaWave", eTimeDay);
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "ALPHA_CH", eTimeDay);
    RecS.VarName2="Charnock coefficient from the wave model";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "TotSurfStr") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TAUTOT", eTimeDay);
    RecS.VarName2="Total Surface stress";
    RecS.minval=0;
    RecS.maxval=0.06;
    RecS.mindiff=-0.01;
    RecS.maxdiff= 0.01;
    RecS.Unit="unknown";
  }
  if (eVarName == "WaveSurfStr") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TAUW", eTimeDay);
    RecS.VarName2="wave supported Surface stress";
    RecS.minval=0;
    RecS.maxval=0.06;
    RecS.mindiff=-0.01;
    RecS.maxdiff= 0.01;
    RecS.Unit="unknown";
  }
  if (eVarName == "SurfStrHF") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "TAUHF", eTimeDay);
    RecS.VarName2="high frequency Surface stress";
    RecS.minval=0;
    RecS.maxval=0.06;
    RecS.mindiff=-0.01;
    RecS.maxdiff= 0.01;
    RecS.Unit="unknown";
  }
  //
  // BFM model variables.
  //
  if (eVarName == "oxygen") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "oxygen", eTimeDay);
    RecS.VarName2="dissolved oxygen concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }

  if (eVarName == "PO4") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "PO4", eTimeDay);
    RecS.VarName2="phosphate concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "NO3") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NO3", eTimeDay);
    RecS.VarName2="nitrate concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "NH4") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NH4", eTimeDay);
    RecS.VarName2="ammonium concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "NitrogenSink") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NitrogenSink", eTimeDay);
    RecS.VarName2="aerobic and anaerobic bacteria";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "SiOH4") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "SiOH4", eTimeDay);
    RecS.VarName2="silicate concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "ReductionEquivalent") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ReductionEquivalent", eTimeDay);
    RecS.VarName2="reduction equivalent";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "bacteriaC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "bateria_c", eTimeDay);
    RecS.VarName2="aerobic and anaerobic bacteria(C)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "bacteriaN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "bateria_n", eTimeDay);
    RecS.VarName2="aerobic and anaerobic bacteria(N)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "bacteriaP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "bateria_p", eTimeDay);
    RecS.VarName2="aerobic and anaerobic bacteria(P)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "diatomsC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "diatoms_c", eTimeDay);
    RecS.VarName2="diatoms carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "diatomsN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "diatoms_n", eTimeDay);
    RecS.VarName2="diatoms nitrogen";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "diatomsP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "diatoms_p", eTimeDay);
    RecS.VarName2="diatoms phosphate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "diatomsL") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "diatoms_l", eTimeDay);
    RecS.VarName2="diatoms chlorophyl";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "diatomsS") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "diatoms_s", eTimeDay);
    RecS.VarName2="diatoms silicate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "flagellatesC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "flagellates_c", eTimeDay);
    RecS.VarName2="flagellates carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "flagellatesN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "flagellates_n", eTimeDay);
    RecS.VarName2="flagellates nitrogen";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "flagellatesP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "flagellates_p", eTimeDay);
    RecS.VarName2="flagellates phosphate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "flagellatesL") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "flagellates_l", eTimeDay);
    RecS.VarName2="flagellates chlorophyl";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "picophytoplanktonC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "picophytoplankton_c", eTimeDay);
    RecS.VarName2="picophytoplankton carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "picophytoplanktonN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "picophytoplankton_n", eTimeDay);
    RecS.VarName2="picophytoplankton nitrogen";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "picophytoplanktonP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "picophytoplankton_p", eTimeDay);
    RecS.VarName2="picophytoplankton phosphate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "picophytoplanktonL") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "picophytoplankton_l", eTimeDay);
    RecS.VarName2="picophytoplankton chlorophyl";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "largephytoplanktonC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "largephytoplankton_c", eTimeDay);
    RecS.VarName2="largephytoplankton carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "largephytoplanktonN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "largephytoplankton_n", eTimeDay);
    RecS.VarName2="largephytoplankton nitrogen";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "largephytoplanktonP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "largephytoplankton_p", eTimeDay);
    RecS.VarName2="largephytoplankton phosphate";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "largephytoplanktonL") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "largephytoplankton_l", eTimeDay);
    RecS.VarName2="largephytoplankton chlorophyl";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "CarnPlanktonC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "CarnPlankton_c", eTimeDay);
    RecS.VarName2="Carnivorous Mesozooplankton(C)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "CarnPlanktonN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "CarnPlankton_n", eTimeDay);
    RecS.VarName2="Carnivorous Mesozooplankton(N)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "CarnPlanktonP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "CarnPlankton_p", eTimeDay);
    RecS.VarName2="Carnivorous Mesozooplankton(P)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "OmniPlanktonC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "OmniPlankton_c", eTimeDay);
    RecS.VarName2="Omnivorous Mesozooplankton(C)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "OmniPlanktonN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "OmniPlankton_n", eTimeDay);
    RecS.VarName2="Omnivorous Mesozooplankton(N)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "OmniPlanktonP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "OmniPlankton_p", eTimeDay);
    RecS.VarName2="Omnivorous Mesozooplankton(P)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "MicroPlanktonC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "MicroPlankton_c", eTimeDay);
    RecS.VarName2="Microzooplankton(C)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "MicroPlanktonN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "MicroPlankton_n", eTimeDay);
    RecS.VarName2="Microzooplankton(N)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "MicroPlanktonP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "MicroPlankton_p", eTimeDay);
    RecS.VarName2="Microzooplankton(P)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "HeteroNanoflagelattesC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_c", eTimeDay);
    RecS.VarName2="Heterotrophic Nanoflagellates (HNAN) C";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "HeteroNanoflagelattesN") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_n", eTimeDay);
    RecS.VarName2="Heterotrophic Nanoflagellates (HNAN) N";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "HeteroNanoflagelattesP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "HeteroNanoflagelattes_p", eTimeDay);
    RecS.VarName2="Heterotrophic Nanoflagellates (HNAN) P";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "LabileDOM1c") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "LabileDOM1_c", eTimeDay);
    RecS.VarName2="Labile Dissolved Organic Matter (Carbon)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "LabileDOM1n") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "LabileDOM1_n", eTimeDay);
    RecS.VarName2="Labile Dissolved Organic Matter (Nitrogen)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "LabileDOM1p") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "LabileDOM1_p", eTimeDay);
    RecS.VarName2="Labile Dissolved Organic Matter (Phosphate)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "LabileDOM2c") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "LabileDOM2_c", eTimeDay);
    RecS.VarName2="Semi-labile Dissolved Organic Carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "RefractoryDOMc") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "RefractoryDOM_c", eTimeDay);
    RecS.VarName2="Semi-refractory Dissolved Organic Carbon";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "ParticleOMc") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ParticleOM_c", eTimeDay);
    RecS.VarName2="Particle Organic Matter (Carbon)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "ParticleOMn") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ParticleOM_n", eTimeDay);
    RecS.VarName2="Particle Organic Matter (Nitrogen)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "ParticleOMp") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ParticleOM_p", eTimeDay);
    RecS.VarName2="Particle Organic Matter (Phosphate)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "ParticleOMs") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ParticleOM_s", eTimeDay);
    RecS.VarName2="Particle Organic Matter (Silicate)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "DissolvedICc") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "DissolvedIC_c", eTimeDay);
    RecS.VarName2="Dissolved Inorganic Carbon C";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "DissolvedICh") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "DissolvedIC_h", eTimeDay);
    RecS.VarName2="Dissolved Inorganic Carbon H";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "Irradiance") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "Irradiance", eTimeDay);
    RecS.VarName2="Irradiance";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="uE/m2s";
  }
  if (eVarName == "DIC") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "DIC", eTimeDay);
    RecS.VarName2="DIC concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mm/m3";
  }
  if (eVarName == "chlorophyll") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "chlorophyll", eTimeDay);
    RecS.VarName2="chlorophyll concentration";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3";
  }
  if (eVarName == "NetProductionP1") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NetProductionP1", eTimeDay);
    RecS.VarName2="Specific Net Production of P1(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="1/d";
  }
  if (eVarName == "NetProductionP2") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NetProductionP2", eTimeDay);
    RecS.VarName2="Specific Net Production of P2(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="1/d";
  }
  if (eVarName == "NetProductionP3") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NetProductionP3", eTimeDay);
    RecS.VarName2="Specific Net Production of P3(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="1/d";
  }
  if (eVarName == "NetProductionP4") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "NetProductionP4", eTimeDay);
    RecS.VarName2="Specific Net Production of P4(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="1/d";
  }
  if (eVarName == "RegFactorP1") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "RegFactorP1", eTimeDay);
    RecS.VarName2="Regular Factor for Light in P1(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "RegFactorP2") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "RegFactorP2", eTimeDay);
    RecS.VarName2="Regular Factor for Light in P2(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "RegFactorP3") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "RegFactorP3", eTimeDay);
    RecS.VarName2="Regular Factor for Light in P3(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "RegFactorP4") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "RegFactorP4", eTimeDay);
    RecS.VarName2="Regular Factor for Light in P4(Phytoplankton)";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="nondim.";
  }
  if (eVarName == "GrossPP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "GrossPP", eTimeDay);
    RecS.VarName2="Gross Primary Production";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3d";
  }
  if (eVarName == "SecondPP") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "SecondPP", eTimeDay);
    RecS.VarName2="Second Primary Production";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="mg/m3d";
  }
  if (eVarName == "ExtinctionCoeff") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ExtinctionCoeff", eTimeDay);
    RecS.VarName2="Total Extinction Coefficient";
    RecS.minval=0;
    RecS.maxval=0.033;
    RecS.mindiff=-0.1;
    RecS.maxdiff=0.1;
    RecS.Unit="1/m";
  }
  

  
  //
  // Now error parsing and assignations
  //
  RecVar eRecVar;
  eRecVar.RecS=RecS;
  if (RecS.VarName2 == "unset") {
    std::cerr << "We did not find the variable\n";
    std::cerr << "eVarName = " << eVarName << "\n";
    std::cerr << "in the list of allowed ones\n";
    std::cerr << "possibly missspelling or lack of relevant code\n";
    throw TerminalException{1};
  }
  if (eModelName != "TRIVIAL") {
    if (RecS.VarNature == "rho") {
      if (F.size() == 0) {
	std::cerr << "VarNature = " << RecS.VarNature << "\n";
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program as authorized variable\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	throw TerminalException{1};
      }
      eRecVar.F=F;
    }
    if (RecS.VarNature == "uv") {
      if (U.size() == 0 || V.size() == 0) {
	std::cerr << "VarNature = " << RecS.VarNature << "\n";
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	throw TerminalException{1};
      }
      eRecVar.U=U;
      eRecVar.V=V;
      if (!IsEqualSizeMatrices(U, V)) {
	std::cerr << "Matrices U and V are not equal sized\n";
	std::cerr << "|U|=" << StringSizeMatrix(U) << "\n";
	std::cerr << "|V|=" << StringSizeMatrix(V) << "\n";
	throw TerminalException{1};
      }
      eRecVar.F=GetNormMatrix(U, V);
    }
    if (RecS.VarNature == "3Drho") {
      auto LDim=Tens3.dimensions();
      if (LDim[0] == 0) {
	std::cerr << "VarNature = " << RecS.VarNature << "\n";
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	throw TerminalException{1};
      }
      eRecVar.Tens3=Tens3;
    }
    if (RecS.VarNature == "3Duv") {
      auto LDim=Uthree.dimensions();
      if (LDim[0] == 0) {
	std::cerr << "VarNature = " << RecS.VarNature << "\n";
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	throw TerminalException{1};
      }
      int dim0=LDim[0];
      int dim1=LDim[1];
      int dim2=LDim[2];
      eRecVar.Uthree=Uthree;
      eRecVar.Vthree=Vthree;
      Eigen::Tensor<double,3> Fwr(dim0, dim1, dim2);
      for (int i0=0; i0<dim0; i0++)
	for (int i1=0; i1<dim1; i1++)
	  for (int i2=0; i2<dim2; i2++) {
	    double eU=Uthree(i0, i1, i2);
	    double eV=Vthree(i0, i1, i2);
	    double eNorm=sqrt(eU*eU + eV*eV);
	    Fwr(i0, i1, i2) = eNorm;
	  }
      eRecVar.Tens3=Fwr;
    }
  }
  return eRecVar;
}


RecVar ModelSpecificVarSpecificTime(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay)
{
  std::string eSep="_";
  std::vector<std::string> ListStr=STRING_Split(eVarName, eSep);
  int len=ListStr.size();
  //  std::cerr << "len=" << len << "\n";
  if (len == 1)
    return ModelSpecificVarSpecificTime_Kernel(TotalArr, eVarName, eTimeDay);
  if (len != 2) {
    std::cerr << "We should clearly rethink len. Expected value is 1 or 2\n";
    std::cerr << "But len=" << len << "\n";
    throw TerminalException{1};
  }
  std::string eVar_rho=ListStr[0];
  std::string eVar_uv=ListStr[1];
  RecVar RecVar_rho=ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_rho, eTimeDay);
  RecVar RecVar_uv =ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_uv , eTimeDay);
  std::string VarNat_rho=RecVar_rho.RecS.VarNature;
  std::string VarNat_uv =RecVar_uv.RecS.VarNature;
  if (VarNat_rho != "rho" && VarNat_rho != "3Drho") {
    std::cerr << "The RecVar_rho is not a rho type variable. Error!\n";
    std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
    std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and SurfCurr as uv variable\n";
    std::cerr << "The call was with eVarName=" << eVarName << "\n";
    throw TerminalException{1};
  }
  if (VarNat_uv != "uv" && VarNat_uv != "3Duv") {
    std::cerr << "The RecVar_uv is not a uv type variable. Error!\n";
    std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
    std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and SurfCurr as uv variable\n";
    std::cerr << "The call was with eVarName=" << eVarName << "\n";
    throw TerminalException{1};
  }
  if ((VarNat_rho == "3Drho" && VarNat_uv == "uv") || (VarNat_rho == "rho" && VarNat_uv == "3Duv") ) {
    std::cerr << "Error. variables do not have the same dimensionality\n";
    std::cerr << "It should be both 3D or both 2D\n";
    std::cerr << "Right now, we have VarNat_rho=" << VarNat_rho << "\n";
    std::cerr << "Right now, we have VarNat_uv =" << VarNat_uv << "\n";
    throw TerminalException{1};
  }
  if (VarNat_rho == "rho") {
    RecVar_rho.U = RecVar_uv.U;
    RecVar_rho.V = RecVar_uv.V;
  }
  if (VarNat_rho == "3Drho") {
    RecVar_rho.Uthree = RecVar_uv.Uthree;
    RecVar_rho.Vthree = RecVar_uv.Vthree;
  }
  RecVar_rho.RecS.VarName1 += "_" + RecVar_uv.RecS.VarName1;
  RecVar_rho.RecS.VarName2 += " + " + RecVar_uv.RecS.VarName2;
  RecVar_rho.RecS.VarNature="uv";
  return RecVar_rho;
}



RecVar RetrieveTrivialRecVar(std::string const& eVarName)
{
  TotalArrGetData TotalArrTrivial;
  TotalArrTrivial.GrdArr.ModelName="TRIVIAL";
  double eTimeDayTrivial=0;
  return ModelSpecificVarSpecificTime(TotalArrTrivial, eVarName, eTimeDayTrivial);
}

std::vector<std::string> GetAllPossibleVariables_with_pairs()
{
  std::vector<std::string> ListVar=GetAllPossibleVariables();
  std::vector<std::string> ListVar_rho;
  std::vector<std::string> ListVar_uv;
  std::vector<std::string> ListVar_3Drho;
  std::vector<std::string> ListVar_3Duv;
  std::vector<std::string> ListVar_Ret;
  for (auto & eVar : ListVar) {
    RecVar eRec=RetrieveTrivialRecVar(eVar);
    if (eRec.RecS.VarNature == "rho")
      ListVar_rho.push_back(eVar);
    if (eRec.RecS.VarNature == "uv")
      ListVar_uv.push_back(eVar);
    if (eRec.RecS.VarNature == "3Drho")
      ListVar_3Drho.push_back(eVar);
    if (eRec.RecS.VarNature == "3Duv")
      ListVar_3Duv.push_back(eVar);
    ListVar_Ret.push_back(eVar);
  }
  for (auto & eVar_uv : ListVar_uv) {
    for (auto & eVar_rho : ListVar_rho) {
      std::string eVarTot = eVar_rho + "_" + eVar_uv;
      ListVar_Ret.push_back(eVarTot);
    }
  }
  for (auto & eVar_uv : ListVar_3Duv) {
    for (auto & eVar_rho : ListVar_3Drho) {
      std::string eVarTot = eVar_rho + "_" + eVar_uv;
      ListVar_Ret.push_back(eVarTot);
    }
  }
  return ListVar_Ret;
}






void ApplyPlotBound(TotalArrGetData const& TotalArr, RecVar & eRecVar, std::string const& eVarName, PlotBound const& ePlotBound)
{
  //
  // Setting up bounds for the plots.
  //
  int nbSingle=ePlotBound.BoundSingle_var.size();
  int nbSingleMin=ePlotBound.BoundSingle_min.size();
  int nbSingleMax=ePlotBound.BoundSingle_max.size();
  //  std::cerr << "nbSingle=" << nbSingle << "\n";
  if (nbSingle != nbSingleMin || nbSingle != nbSingleMax) {
    std::cerr << "Number of entries in BoundSingle_var, BoundSingle_min, BoundSingle_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbSingle    = " << nbSingle << "\n";
    std::cerr << "nbSingleMin = " << nbSingleMin << "\n";
    std::cerr << "nbSingleMax = " << nbSingleMax << "\n";
  }
  for (int iS=0; iS<nbSingle; iS++)
    if (ePlotBound.BoundSingle_var[iS] == eVarName) {
      eRecVar.RecS.minval=ePlotBound.BoundSingle_min[iS];
      eRecVar.RecS.maxval=ePlotBound.BoundSingle_max[iS];
    }
  int nbDiff=ePlotBound.BoundDiff_var.size();
  int nbDiffMin=ePlotBound.BoundDiff_min.size();
  int nbDiffMax=ePlotBound.BoundDiff_max.size();
  //  std::cerr << "nbDiff=" << nbDiff << "\n";
  if (nbDiff != nbDiffMin || nbDiff != nbDiffMax) {
    std::cerr << "Number of entries in BoundDiff_var, BoundDiff_min, BoundDiff_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbDiff    = " << nbDiff << "\n";
    std::cerr << "nbDiffMin = " << nbDiffMin << "\n";
    std::cerr << "nbDiffMax = " << nbDiffMax << "\n";
  }
  for (int iD=0; iD<nbDiff; iD++)
    if (ePlotBound.BoundDiff_var[iD] == eVarName) {
      eRecVar.RecS.mindiff=ePlotBound.BoundDiff_min[iD];
      eRecVar.RecS.maxdiff=ePlotBound.BoundDiff_max[iD];
    }
  int eSize=eRecVar.F.size();
  if (ePlotBound.VariableRange && eSize > 0) {
    PairMinMax ePair=ComputeMinMax(TotalArr.GrdArr, eRecVar.F);
    eRecVar.RecS.minval=ePair.TheMin;
    eRecVar.RecS.maxval=ePair.TheMax;
  }
}


void ApplyPlotBoundPair(TotalArrGetData const& TotalArr1, TotalArrGetData const& TotalArr2, RecVar & eRecVar1, RecVar & eRecVar2, std::string const& eVarName, PlotBound const& ePlotBound)
{

  int eSize=eRecVar1.F.size();
  std::cerr << "ApplyPlotBoundPair : eSize = " << eSize << "\n";
  if (ePlotBound.VariableRange && eSize > 0) {
    MyMatrix<double> eDiff12=eRecVar1.F - eRecVar2.F;
    PairMinMax ePair=ComputeMinMax(TotalArr1.GrdArr, eDiff12);
    double MaxChange=std::max(ePair.TheMax, -ePair.TheMin);
    std::cerr << "ApplyPlotBoundPair : min/max = " << ePair.TheMin << " / " << ePair.TheMax << "\n";
    eRecVar1.RecS.mindiff=-MaxChange;
    eRecVar1.RecS.maxdiff= MaxChange;
    eRecVar2.RecS.mindiff=-MaxChange;
    eRecVar2.RecS.maxdiff= MaxChange;
  }
}




RecVar ModelSpecificVarSpecificTimeBound(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay, PlotBound const& ePlotBound)
{
  RecVar eRecVar=ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
  ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
  return eRecVar;
}




std::string GetStrAllOfPlot(VarQuery const& eQuery)
{
  int iTime=eQuery.iTime;
  std::string strAll;
  std::string strFile;
  if (eQuery.typeQuery == "direct") {
    strFile=DATE_ConvertMjd2mystringFile(eQuery.eTimeDay);
  }
  else {
    std::vector<int> eDate=DATE_ConvertMjd2six(eQuery.eTimeDay);
    int iYear=eDate[0];
    int iMon=eDate[1];
    if (eQuery.typeQuery == "monthly") {
      std::vector<std::string> ListMon{"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
      strFile=ListMon[iMon-1] + IntToString(iYear);
    }
    if (eQuery.typeQuery == "seasonal") {
      std::vector<std::string> ListSeas{"Winter", "Spring", "Summer", "Autumn"};
      int iSeas=(iMon - 1)/3;
      strFile=ListSeas[iSeas] + IntToString(iYear);
    }
  }
  if (iTime == -1) {
    strAll=strFile;
  }
  else {
    if (iTime > 10000) {
      std::cerr << "Error in the code\n";
      std::cerr << "iTime is too large\n";
      throw TerminalException{1};
    }
    strAll=StringNumber(iTime, 4) + "_" + strFile;
  }
  if (eQuery.NatureQuery != "instant")
    strAll += "_" + eQuery.NatureQuery;
  return strAll;
}



std::string GetStrPresOfPlot(VarQuery const& eQuery)
{
  std::string strPres1=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay);
  if (eQuery.NatureQuery == "instant") {
    return "at " + strPres1;
  }
  if (eQuery.NatureQuery == "swathMax") {
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
    return "max from " + strPres1 + " to " + strPres2;
  }
  if (eQuery.NatureQuery == "swathMin") {
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
    return "min from " + strPres1 + " to " + strPres2;
  }
  if (eQuery.NatureQuery == "average") {
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
    return "avg. from " + strPres1 + " to " + strPres2;
  }
  std::cerr << "Failed to find NatureQuery in list of available options\n";
  std::cerr << "eQuery.NatureQuery=" << eQuery.NatureQuery << "\n";
  throw TerminalException{1};
}



RecVar ModelSpecificVarSpecificTimeGeneral(TotalArrGetData const& TotalArr, std::string const& eVarName, VarQuery const& eQuery, PlotBound const& ePlotBound)
{
  //
  // Check correctness 
  //
  std::vector<std::string> ListAllow{"instant", "average", "swathMax", "swathMin"};
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) == ListAllow.end()) {
    std::cerr << "We failed to find NatureQuery=" << eQuery.NatureQuery << "\n";
    std::cerr << "List of allowed queries:\n";
    for (auto & eStr : ListAllow)
      std::cerr << "  eStr=" << eStr << "\n";
    throw TerminalException{1};
  }
  std::string strPres=DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
  std::cerr << "Query ModelSpecificVarSpecificTimeGeneral NatureQuery=" << eQuery.NatureQuery << " date=" << strPres << " VarName=" << eVarName << "\n";
  //
  // Reading array and doing operations "instant", "average", etc.
  //
  RecVar eRecVar;
  if (eQuery.NatureQuery == "instant") {
    eRecVar=ModelSpecificVarSpecificTimeBound(TotalArr, eVarName, eQuery.eTimeDay, ePlotBound);
  }
  else {
    double eTimeDay=eQuery.eTimeDay;
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::vector<int> ListRelITime=GetIntervalListITime(TotalArr.eArr, eTimeDay, TimeFrameDay);
    int nbTimeRel=ListRelITime.size();
    if (nbTimeRel == 0) {
      std::cerr << "  NatureQuery = " << eQuery.NatureQuery << "\n";
      std::cerr << "  We have nbTimeRel = 0 therefore we cannot process rhe query\n";
      std::cerr << "  eTimeDay=" << eTimeDay << "\n";
      std::cerr << "  TimeFrameDay=" << TimeFrameDay << "\n";
      throw TerminalException{1};
    }
    RecVar RecVarTrivial=RetrieveTrivialRecVar(eVarName);
    MyMatrix<double> F, U, V;
    Eigen::Tensor<double,3> Tens3;
    Eigen::Tensor<double,3> Uthree;
    Eigen::Tensor<double,3> Vthree;
    std::cerr << "nbTimeRel = " << nbTimeRel << "\n";
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++) {
      int iTime=ListRelITime[iTimeRel];
      double eTimeDayB=ARR_GetTime(TotalArr.eArr, iTime);
      std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDayB);
      eRecVar=ModelSpecificVarSpecificTimeBound(TotalArr, eVarName, eTimeDayB, ePlotBound);
      std::cerr << "iTimeRel=" << iTimeRel << "/" << nbTimeRel <<  " date=" << strPres << "\n";
      //      std::cerr << "iTimeRel=" << iTimeRel << " iTime=" << iTime << " VarNature=" << eRecVar.RecS.VarNature << "\n";
      if (iTimeRel == 0) {
	if (RecVarTrivial.RecS.VarNature == "rho")
	  F=eRecVar.F;
	if (RecVarTrivial.RecS.VarNature == "uv") {
	  U=eRecVar.U;
	  V=eRecVar.V;
	  F=eRecVar.F;
	}
	if (RecVarTrivial.RecS.VarNature == "3Drho")
	  Tens3=eRecVar.Tens3;
	if (RecVarTrivial.RecS.VarNature == "3Duv") {
	  Uthree=eRecVar.Uthree;
	  Vthree=eRecVar.Vthree;
	  Tens3=eRecVar.Tens3;
	}
      }
      else {
	if (eQuery.NatureQuery == "average") {
	  if (RecVarTrivial.RecS.VarNature == "rho")
	    F += eRecVar.F;
	  if (RecVarTrivial.RecS.VarNature == "uv") {
	    U += eRecVar.U;
	    V += eRecVar.V;
	    F += eRecVar.F;
	  }
	  if (RecVarTrivial.RecS.VarNature == "3Drho")
	    Tens3 += eRecVar.Tens3;
	  if (RecVarTrivial.RecS.VarNature == "3Duv") {
	    Uthree += eRecVar.Uthree;
	    Vthree += eRecVar.Vthree;
	    Tens3 += eRecVar.Tens3;
	  }
	}
	if (eQuery.NatureQuery == "swathMax") {
	  if (RecVarTrivial.RecS.VarNature == "rho")
	    F=F.cwiseMax(eRecVar.F);
	  if (RecVarTrivial.RecS.VarNature == "3Drho")
	    Tens3=Tens3.cwiseMax(eRecVar.Tens3);
	  if (RecVarTrivial.RecS.VarNature == "uv" || RecVarTrivial.RecS.VarNature == "3Duv") {
	    std::cerr << "swathMax for uv does not have any sense\n";
	    throw TerminalException{1};
	  }
	}
	if (eQuery.NatureQuery == "swathMin") {
	  if (RecVarTrivial.RecS.VarNature == "rho")
	    F=F.cwiseMin(eRecVar.F);
	  if (RecVarTrivial.RecS.VarNature == "3Drho")
	    Tens3=Tens3.cwiseMin(eRecVar.Tens3);
	  if (RecVarTrivial.RecS.VarNature == "uv" || RecVarTrivial.RecS.VarNature == "3Duv") {
	    std::cerr << "swathMin for uv does not have any sense\n";
	    throw TerminalException{1};
	  }
	}
      }
    }
    if (eQuery.NatureQuery == "average") {
      if (RecVarTrivial.RecS.VarNature == "rho") {
	F /= double(nbTimeRel);
      }
      if (RecVarTrivial.RecS.VarNature == "uv") {
	U /= double(nbTimeRel);
	V /= double(nbTimeRel);
	F /= double(nbTimeRel);
      }
      if (RecVarTrivial.RecS.VarNature == "3Drho") {
	// does not compile
	// Tens3 /= double(nbTimeRel);
	auto LDim=Tens3.dimensions();
	int dim0=LDim[0];
	int dim1=LDim[1];
	int dim2=LDim[2];
	for (int i0=0; i0<dim0; i0++)
	  for (int i1=0; i1<dim1; i1++)
	    for (int i2=0; i2<dim2; i2++)
	      Tens3(i0, i1, i2) /= double(nbTimeRel);
      }
      if (RecVarTrivial.RecS.VarNature == "3Duv") {
	// does not compile
	// Tens3 /= double(nbTimeRel);
	auto LDim=Tens3.dimensions();
	int dim0=LDim[0];
	int dim1=LDim[1];
	int dim2=LDim[2];
	for (int i0=0; i0<dim0; i0++)
	  for (int i1=0; i1<dim1; i1++)
	    for (int i2=0; i2<dim2; i2++) {
	      Uthree(i0, i1, i2) /= double(nbTimeRel);
	      Vthree(i0, i1, i2) /= double(nbTimeRel);
	      Tens3 (i0, i1, i2) /= double(nbTimeRel);
	    }
      }
    }
    if (RecVarTrivial.RecS.VarNature == "rho") {
      eRecVar.F=F;
    }
    if (RecVarTrivial.RecS.VarNature == "3Drho") {
      eRecVar.Tens3=Tens3;
    }
    if (RecVarTrivial.RecS.VarNature == "uv") {
      eRecVar.U=U;
      eRecVar.V=V;
      eRecVar.F=F;
    }
    if (RecVarTrivial.RecS.VarNature == "3Duv") {
      eRecVar.Uthree=Uthree;
      eRecVar.Vthree=Vthree;
      eRecVar.Tens3=Tens3;
    }
  }
  //  std::cerr << "End of ModelSpecificVarSpecificTimeGeneral\n";
  std::string strAll=GetStrAllOfPlot(eQuery);
  std::string strPresPlot=GetStrPresOfPlot(eQuery);
  //  std::cerr << "strAll=" << strAll << "\n";
  eRecVar.RecS.strAll=strAll;
  eRecVar.RecS.strPres=strPresPlot;
  ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
  //  std::cerr << "After ApplyPlotBound\n";
  return eRecVar;
}



PairRecVar ModelPairSpecificVarSpecificTimeGeneral(TotalArrGetData const& TotalArr1, TotalArrGetData const& TotalArr2, std::string const& eVarName, VarQuery const& eQuery, PlotBound const& ePlotBound)
{
  //  std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
  std::vector<std::string> ListAllow{"instant", "average", "swathMax", "swathMin"};
  RecVar eRecVar1, eRecVar2;
  //  std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral\n";
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) != ListAllow.end()) {
    eRecVar1=ModelSpecificVarSpecificTimeGeneral(TotalArr1, eVarName, eQuery, ePlotBound);
    eRecVar2=ModelSpecificVarSpecificTimeGeneral(TotalArr2, eVarName, eQuery, ePlotBound);
    ApplyPlotBoundPair(TotalArr1, TotalArr2, eRecVar1, eRecVar2, eVarName, ePlotBound);
    //    std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral, case 1\n";
  }
  else {
    //    std::cerr << "Begin ModelPairSpecificVarSpecificTimeGeneral, case 2\n";
    std::vector<std::string> ListAllowB{"MaxDiff", "MinDiff"};
    if (std::find(ListAllowB.begin(), ListAllowB.end(), eQuery.NatureQuery) == ListAllowB.end()) {
      std::cerr << "Only two difference operator are allowed: MaxDiff and MinDiff\n";
      std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
    }
    RecVar RecVarTrivial=RetrieveTrivialRecVar(eVarName);
    //
    double eTimeDay=eQuery.eTimeDay;
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::vector<int> ListRelITime1=GetIntervalListITime(TotalArr1.eArr, eTimeDay, TimeFrameDay);
    std::vector<int> ListRelITime2=GetIntervalListITime(TotalArr2.eArr, eTimeDay, TimeFrameDay);
    int nbTimeRel1=ListRelITime1.size();
    int nbTimeRel2=ListRelITime2.size();
    if (nbTimeRel1 == 0 || nbTimeRel2 == 0) {
      std::cerr << "We have nbTimeRel1 = " << nbTimeRel1 << "\n";
      std::cerr << "We have nbTimeRel2 = " << nbTimeRel2 << "\n";
      std::cerr << "We need both to be non-zero\n";
      throw TerminalException{1};
    }
    double DeltaTime1=(ARR_GetTime(TotalArr1.eArr, ListRelITime1[nbTimeRel1-1]) - ARR_GetTime(TotalArr1.eArr, ListRelITime1[0]))/double(nbTimeRel1-1);
    double DeltaTime2=(ARR_GetTime(TotalArr2.eArr, ListRelITime2[nbTimeRel2-1]) - ARR_GetTime(TotalArr2.eArr, ListRelITime2[0]))/double(nbTimeRel2-1);
    double DeltaTime=std::min(DeltaTime1, DeltaTime2);
    double FirstTime=eTimeDay;
    double LastTime=eTimeDay + TimeFrameDay;
    std::vector<double> ListRelTimeDay=GetIntervalFLD(FirstTime, LastTime, DeltaTime);
    int nbTimeRel=ListRelTimeDay.size();
    MyMatrix<double> diffF;
    if (RecVarTrivial.RecS.VarNature == "uv") {
      std::cerr << "Cannot do MaxDiff or MinDiff for variables of type uv\n";
      throw TerminalException{1};
    }
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++) {
      double eTimeDayB=ListRelTimeDay[iTimeRel];
      eRecVar1=ModelSpecificVarSpecificTimeBound(TotalArr1, eVarName, eTimeDayB, ePlotBound);
      eRecVar2=ModelSpecificVarSpecificTimeBound(TotalArr2, eVarName, eTimeDayB, ePlotBound);
      if (iTimeRel == 0) {
	if (RecVarTrivial.RecS.VarNature == "rho") {
	  diffF = eRecVar1.F - eRecVar2.F;
	}
      }
      else {
	if (eQuery.NatureQuery == "MaxDiff") {
	  if (RecVarTrivial.RecS.VarNature == "rho")
	    diffF=diffF.cwiseMax(eRecVar1.F - eRecVar2.F);
	}
	if (eQuery.NatureQuery == "MinDiff") {
	  if (RecVarTrivial.RecS.VarNature == "rho")
	    diffF=diffF.cwiseMin(eRecVar1.F - eRecVar2.F);
	}
      }
    }
    if (RecVarTrivial.RecS.VarNature == "rho") {
      eRecVar1.F=diffF;
      eRecVar2.F.fill(double(0));
    }
    ApplyPlotBound(TotalArr1, eRecVar1, eVarName, ePlotBound);
    ApplyPlotBound(TotalArr2, eRecVar2, eVarName, ePlotBound);
    ApplyPlotBoundPair(TotalArr1, TotalArr2, eRecVar1, eRecVar2, eVarName, ePlotBound);
  }
  std::string strAll=GetStrAllOfPlot(eQuery);
  eRecVar1.RecS.strAll=strAll;
  eRecVar2.RecS.strAll=strAll;
  return {eRecVar1, eRecVar2};
}


RecVar GetTrivialArrayPlot(GridArray const& GrdArr)
{
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> F(eta_rho, xi_rho);
  double MaxLon=GrdArr.GrdArrRho.LON.maxCoeff();
  double MinLon=GrdArr.GrdArrRho.LON.minCoeff();
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eLon=GrdArr.GrdArrRho.LON(i,j);
      double eVal=(eLon - MinLon)/(MaxLon - MinLon);
      F(i, j)=eVal;
    }
  RecVar eRecVar;
  eRecVar.RecS.strAll="trivialarrayplot";
  eRecVar.RecS.VarName1="Track";
  eRecVar.RecS.VarName2="Track";
  eRecVar.RecS.minval=2;
  eRecVar.RecS.maxval=3;
  eRecVar.RecS.Unit="nondim.";
  eRecVar.F=F;
  return eRecVar;
}


#endif
