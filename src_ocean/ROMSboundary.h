#ifndef ROMS_BOUNDARY
#define ROMS_BOUNDARY


#include "ROMSfunctionality.h"
#include "NamelistExampleOcean.h"
#include "NCL_Kernel.h"
#include "Basic_plot.h"


void BOUND_Plotting_Function(FullNamelist const& eFull)
{
  struct ArrSide {
    std::string InputName;
    std::string NcName;
    int eta_ts;
    int eta_u;
    int eta_v;
    MyVector<double> DEP_rho;
    MyVector<double> DEP_u;
    MyVector<double> DEP_v;
  };
  struct TypeVar {
    std::string VarName;
    std::string SystemName;
    std::string Nature;
  };
  auto GetMatrixSide=[](MyMatrix<double> const& M, std::string const& eSide) -> MyVector<double> {
    int eta=M.rows();
    int xi=M.cols();
    if (eSide == "South") {
      MyVector<double> V(xi);
      for (int i=0; i<xi; i++)
	V(i) = M(0,i);
      return V;
    }
    if (eSide == "North") {
      MyVector<double> V(xi);
      for (int i=0; i<xi; i++)
	V(i) = M(eta-1,i);
      return V;
    }
    if (eSide == "East") {
      MyVector<double> V(eta);
      for (int i=0; i<eta; i++)
	V(i) = M(i,xi-1);
      return V;
    }
    if (eSide == "West") {
      MyVector<double> V(eta);
      for (int i=0; i<eta; i++)
	V(i) = M(i,0);
      return V;
    }
    std::cerr << "Failed to find Matching entry in GetMatrixSide\n";
    throw TerminalException{1};
  };
  auto GetVectorDEP=[&](std::string const& typeName, ArrSide const& eSide) -> MyVector<double> {
    if (typeName == "rho")
      return eSide.DEP_rho;
    if (typeName == "u")
      return eSide.DEP_u;
    if (typeName == "v")
      return eSide.DEP_v;
    std::cerr << "Failed to find Matching entry in GetVectorDEP\n";
    throw TerminalException{1};
  };
  std::string strSRho="s_rho";
  //
  // PROC entries
  //
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string BoundaryFile=eBlPROC.ListStringValues.at("BoundaryFile");
  netCDF::NcFile dataFile(BoundaryFile, netCDF::NcFile::read);
  int s_rho  =NC_ReadDimension(dataFile, strSRho);
  std::vector<double> ListTime=NC_ReadTimeFromFile(BoundaryFile, "zeta_time");
  int nbTime =ListTime.size();
  //  bool WriteITimeInFileName=eBlPROC.ListBoolValues.at("WriteITimeInFileName");
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  GridArray GrdArr=NC_ReadRomsGridFile(GridFile);
  ARVDtyp ARVD = ReadROMSverticalStratification(BoundaryFile);
  std::cerr << "PROC entries read\n";
  //
  // PLOT parameter
  //
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  bool PlotTemp = eBlPLOT.ListBoolValues.at("PlotTemp");
  bool PlotSalt = eBlPLOT.ListBoolValues.at("PlotSalt");
  bool PlotU = eBlPLOT.ListBoolValues.at("PlotU");
  bool PlotV = eBlPLOT.ListBoolValues.at("PlotV");
  bool VariableRange=eBlPLOT.ListBoolValues.at("VariableRange");
  std::vector<std::string> ListSides = eBlPLOT.ListListStringValues.at("ListSides");
  std::vector<std::string> ListSidesTot = {"South", "North", "West", "East"};
  std::vector<ArrSide> ListArrSide;
  for (auto & eStr : ListSides) {
    bool IsPresent=false;
    for (auto & eStrTot : ListSidesTot)
      if (eStr == eStrTot)
	IsPresent=true;
    if (!IsPresent) {
      std::cerr << "The variable eStr = " << eStr << " is not allowed\n";
      std::cerr << "Allowed variables are South, North, West, East\n";
      throw TerminalException{1};
    }
    ArrSide eArrSide;
    eArrSide.InputName = eStr;
    eArrSide.NcName = UpperCaseToLowerCase(eStr);
    eArrSide.DEP_rho = GetMatrixSide(GrdArr.GrdArrRho.DEP, eStr);
    eArrSide.DEP_u = GetMatrixSide(GrdArr.GrdArrU.DEP, eStr);
    eArrSide.DEP_v = GetMatrixSide(GrdArr.GrdArrV.DEP, eStr);
    ListArrSide.push_back(eArrSide);
  }
  std::vector<TypeVar> ListTypeVar;
  if (PlotTemp) {
    ListTypeVar.push_back({"temp", "Temp", "rho"});
  }
  if (PlotSalt) {
    ListTypeVar.push_back({"salt", "Salt", "rho"});
  }
  if (PlotU) {
    ListTypeVar.push_back({"u", "Curr", "u"});
  }
  if (PlotV) {
    ListTypeVar.push_back({"v", "Curr", "v"});
  }
  int nbLevelSpa = eBlPLOT.ListIntValues.at("nbLevelSpa");
  std::cerr << "PLOT entries read\n";

  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  std::cerr << "ePerm obtained\n";
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  std::cerr << "eCall obtained\n";

  for (int iTime=0; iTime<nbTime; iTime++) {
    double eTimeDay=ListTime[iTime];
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << " date=" << strPres << "\n";
    for (auto& eArrSide : ListArrSide) {
      for (auto& eTypeVar : ListTypeVar) {
	std::string varName = eTypeVar.VarName + "_" + eArrSide.NcName;
	netCDF::NcVar data=dataFile.getVar(varName);
	//
	// Reading the bathymetry
	//
	MyVector<double> DEP = GetVectorDEP(eTypeVar.Nature, eArrSide);
	int siz=DEP.size();
	double maxDep = DEP.maxCoeff();
	int NbVert=100;
	double DeltaZ = maxDep / double(NbVert);
	MyVector<double> ListVertPos(NbVert+1);
	for (int i=0; i<=NbVert; i++) {
	  double eVertPos = -maxDep + i*DeltaZ;
	  ListVertPos(i)=eVertPos;
	}
	//
	// Reading the data sets
	//
	std::vector<size_t> start{size_t(iTime), 0, 0};
	std::vector<size_t> count{1, size_t(s_rho), size_t(siz)};
	MyVector<double> eVal=NC_ReadVariable_data_start_count(data, start, count);
	MyMatrix<double> M(s_rho, siz);
	int idx=0;
	for (int i=0; i<s_rho; i++)
	  for (int j=0; j<siz; j++) {
	    M(i,j)=eVal[idx];
	    idx++;
	  }
        //        std::cerr << "M assigned\n";
	//
	// The vertical coordinate
	//
	MyMatrix<double> F(NbVert+1, siz);
	MyMatrix<uint8_t> MSK(NbVert+1, siz);
	for (int i=0; i<siz; i++) {
	  double eZeta= 0;
	  double eDep = DEP(i);
	  MyVector<double> Zr_out = GetVertCoord_R(ARVD, eDep, eZeta);
	  double eps=0.00001;
	  for (int iV=0; iV<=NbVert; iV++) {
	    double eVert = ListVertPos(iV);
	    int eMSK=0;
	    double eF=0;
	    if (eVert >= Zr_out(0)) {
	      eMSK=1;
	      if (eVert > Zr_out(s_rho-1) - eps) {
		eF = M(s_rho-1, i);
	      }
	      else {
		bool IsMatch=false;
		for (int iS=0; iS<s_rho-1; iS++) {
		  double dep1 = Zr_out(iS);
		  double dep2 = Zr_out(iS+1);
		  if (eVert >= dep1 - eps && eVert <= dep2 + eps) {
		    IsMatch=true;
		    double alpha1=(dep2 - eVert)/(dep2 - dep1);
		    double alpha2=(eVert - dep1)/(dep2 - dep1);
		    eF = M(iS,i) * alpha1 + M(iS+1,i) * alpha2;
		  }
		}
		if (!IsMatch) {
		  std::cerr << "Failed to find matching depth\n";
		  throw TerminalException{1};
		}
	      }
	    }
	    MSK(iV,i) = eMSK;
	    F(iV,i) = eF;
	  }
	}
        //        std::cerr << "We have F and MSK\n";
	//
	// The eta/xi coordinate
	//
	MyMatrix<double> LON(NbVert+1, siz);
	MyMatrix<double> LAT(NbVert+1, siz);
	for (int i=0; i<siz; i++)
	  for (int iV=0; iV<=NbVert; iV++) {
	    LON(iV,i) = double(i);
	    LAT(iV,i) = ListVertPos(iV);
	  }
        //        std::cerr << "We have LON/LAT\n";
	//
	// The grid array
	//
	GridArray GrdArr;
	GrdArr.IsFE = 0;
	GrdArr.IsSpherical=false;
	GrdArr.GrdArrRho.LON=LON;
	GrdArr.GrdArrRho.LAT=LAT;
	GrdArr.GrdArrRho.MSK=MSK;
        //        std::cerr << "We have GrdArr\n";
	//
	// The plotting function
	//
	TotalArrGetData TotalArrTrivial;
	TotalArrTrivial.GrdArr.ModelName="TRIVIAL";
	RecVar eRecVarTriv = ModelSpecificVarSpecificTime(TotalArrTrivial, eTypeVar.SystemName, eTimeDay);
	RecVar NewRecVar;
	NewRecVar.RecS=eRecVarTriv.RecS;
	if (VariableRange) {
	  PairMinMax ePair=ComputeMinMaxMask(MSK, F);
	  NewRecVar.RecS.mindiff=ePair.TheMin;
	  NewRecVar.RecS.maxdiff=ePair.TheMax;
	  NewRecVar.RecS.minval=ePair.TheMin;
	  NewRecVar.RecS.maxval=ePair.TheMax;
	}
	NewRecVar.F = F;
        //        std::cerr << "We have NewRecVar\n";
	//
	DrawArr eDrawArr=ePerm.eDrawArr;
	eDrawArr.DrawContourBathy=false;
        eDrawArr.DoTitle = true;
	eDrawArr.TitleStr=varName + " " + NewRecVar.RecS.strPres;
        eDrawArr.nbLevelSpa = nbLevelSpa;
        //        std::cerr << "We have eDrawArr\n";
	//
	eDrawArr.VarNameUF = varName;
	std::string FileName=ePerm.eDir + varName + "_" + NewRecVar.RecS.strAll;
	//
	PLOT_PCOLOR(FileName, GrdArr, eDrawArr, NewRecVar, eCall, ePerm);
      }
    }
  }
}








#endif
