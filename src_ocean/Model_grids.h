#ifndef MODEL_GRIDS_INCLUDE
#define MODEL_GRIDS_INCLUDE

#include "Basic_netcdf.h"
#include "Basic_grib.h"
#include "SphericalGeom.h"
#include "Triangulations.h"
#include "CommonFuncModel.h"


QuadArray GetQuadArray(GridArray const& GrdArr)
{
  double MinLon=0, MaxLon=0, MinLat=0, MaxLat=0;
  if (GrdArr.IsFE == 1) {
    int siz1=GrdArr.GrdArrRho.LON.rows();
    int siz2=GrdArr.GrdArrRho.LON.cols();
    if (siz1 == 0 || siz2 == 0) {
      std::cerr << "We need to have a nontrivial matrix\n";
      std::cerr << "siz1=" << siz1 << " siz2=" << siz2 << "\n";
      throw TerminalException{1};
    }
    MinLon=GrdArr.GrdArrRho.LON.minCoeff();
    MaxLon=GrdArr.GrdArrRho.LON.maxCoeff();
    MinLat=GrdArr.GrdArrRho.LAT.minCoeff();
    MaxLat=GrdArr.GrdArrRho.LAT.maxCoeff();
  }
  else {
    bool IsFirst=true;
    int eta_rho=GrdArr.GrdArrRho.LON.rows();
    int xi_rho =GrdArr.GrdArrRho.LON.cols();
    int eta_rho_msk=GrdArr.GrdArrRho.MSK.rows();
    int xi_rho_msk =GrdArr.GrdArrRho.MSK.cols();
    if (eta_rho_msk != eta_rho || xi_rho_msk != xi_rho) {
      std::cerr << "Dimension error in the arrays\n";
      throw TerminalException{1};
    }
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++)
	if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	  double eLon=GrdArr.GrdArrRho.LON(i,j);
	  double eLat=GrdArr.GrdArrRho.LAT(i,j);
	  if (IsFirst) {
	    MinLon=eLon;
	    MaxLon=eLon;
	    MinLat=eLat;
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
  }
  return {MinLon, MaxLon, MinLat, MaxLat};
}


QuadArray GetQuadArray(MyMatrix<double> const& LON, MyMatrix<double> const& LAT, double const& deltaLL)
{
  double MinLon=LON.minCoeff() - deltaLL;
  double MaxLon=LON.maxCoeff() + deltaLL;
  double MinLat=LAT.minCoeff() - deltaLL;
  double MaxLat=LAT.maxCoeff() + deltaLL;
  return {MinLon, MaxLon, MinLat, MaxLat};
}


std::ostream& operator<<(std::ostream& os, QuadArray const& eQ)
{
  os << "(Lon(min/max)=" << eQ.MinLon << " / " << eQ.MaxLon << " Lat(min/max)=" << eQ.MinLat << " / " << eQ.MaxLat << ")";
  return os;
}




std::vector<std::string> GetAllPossibleModels()
{
  std::vector<std::string> vec{"COSMO", "WAM", "ROMS", "ROMS_IVICA", "WWM", "WWM_DAILY", "WW3", "GRIB_DWD", "GRIB_ALADIN", "GRIB_ECMWF", "GRIB_GFS", "GRIB_IFS", "GRIB_COSMO", "GRIB_WAM_FORT30", "SCHISM_SFLUX", "SCHISM_NETCDF_OUT", "RECTANGULAR", "WRF", "UNRUNOFF", "IVICA_UVP", "NEMO"};
  return vec;
}


std::string GetKernelModelName(std::string const& eModelName)
{
  std::vector<std::string> ListStr=STRING_Split(eModelName, ":");
  return ListStr[0];
}





void CHECK_Model_Allowedness(std::string const& PreModelName)
{
  std::string eModelName = GetKernelModelName(PreModelName);
  std::vector<std::string> vec=GetAllPossibleModels();
  bool isPresent = (std::find(vec.begin(), vec.end(), eModelName) != vec.end());
  if (!isPresent) {
    std::cerr << "We did not find the MODEL NAME\n";
    std::cerr << "eModelName = " << eModelName << "     PreModelName = " << PreModelName << "\n";
    std::cerr << "List of allowed models\n";
    for (size_t iModel=0; iModel<vec.size(); iModel++)
      std::cerr << "iModel=" << iModel << " eModel=" << vec[iModel] << "\n";
    throw TerminalException{1};
  }
}


void InitializeIdxJdxWet(CoordGridArrayFD & eCoordGrdArr)
{
  int eta=eCoordGrdArr.LON.rows();
  int xi=eCoordGrdArr.LON.cols();
  eCoordGrdArr.Idx.clear();
  eCoordGrdArr.Jdx.clear();
  int nbWet=0;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++)
      if (eCoordGrdArr.MSK(i, j) == 1) {
        nbWet++;
        eCoordGrdArr.Idx.push_back(i);
        eCoordGrdArr.Jdx.push_back(j);
      }
  eCoordGrdArr.nbWet=nbWet;
}


bool TestEqualityGridArray(GridArray const& GrdArr1, GridArray const& GrdArr2)
{
  if (GrdArr1.IsFE != GrdArr2.IsFE)
    return false;
  if (GrdArr1.IsFE) {
    if (GrdArr1.INE.rows() != GrdArr2.INE.rows())
      return false;
    int nbTrig=GrdArr1.INE.rows();
    for (int iTrig=0; iTrig<nbTrig; iTrig++)
      for (int i=0; i<3; i++)
	if (GrdArr1.INE(iTrig,i) != GrdArr2.INE(iTrig,i))
	  return false;
  }
  int eta1=GrdArr1.GrdArrRho.LON.rows();
  int eta2=GrdArr2.GrdArrRho.LON.rows();
  if (eta1 != eta2)
    return false;
  int eta=eta1;
  int xi1=GrdArr1.GrdArrRho.LON.cols();
  int xi2=GrdArr2.GrdArrRho.LON.cols();
  if (xi1 != xi2)
    return false;
  if (GrdArr1.GrdArrRho.HaveDEP != GrdArr2.GrdArrRho.HaveDEP)
    return false;
  bool HaveDEP=GrdArr2.GrdArrRho.HaveDEP;
  //  std::cerr << "GrdArr1.ModelName = " << GrdArr1.ModelName << "\n";
  //  std::cerr << "GrdArr2.ModelName = " << GrdArr2.ModelName << "\n";
  int xi=xi1;
  double err=0;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      //      std::cerr << "i=" << i << " j=" << j << "\n";
      double lon1=GrdArr1.GrdArrRho.LON(i,j);
      //      std::cerr << "step 1\n";
      double lon2=GrdArr2.GrdArrRho.LON(i,j);
      //      std::cerr << "step 2\n";
      err += fabs(lon1 - lon2);
      double lat1=GrdArr1.GrdArrRho.LAT(i,j);
      //      std::cerr << "step 3\n";
      double lat2=GrdArr2.GrdArrRho.LAT(i,j);
      //      std::cerr << "step 4\n";
      err += fabs(lat1 - lat2);
      if (HaveDEP) {
	double dep1=GrdArr1.GrdArrRho.DEP(i,j);
	//	std::cerr << "step 5\n";
	double dep2=GrdArr2.GrdArrRho.DEP(i,j);
	//	std::cerr << "step 6\n";
	err += fabs(dep1 - dep2);
      }
    }
  if (err > double(1))
    return false;
  return true;
}





GridArray NC_ReadRomsGridFile(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Missing roms grid file eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  std::cerr << "eFile = " << eFile << "\n";
  std::function<int(double const&)> fConv=[](double const& x) -> int {
    return int(x);
  };
  std::string xName, yName;
  if (NC_IsVar(eFile, "lon_rho")) {
    xName = "lon";
    yName = "lat";
  }
  else {
    xName = "x";
    yName = "y";
  }
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="ROMS";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  // Rho part of the arrays
  GrdArr.GrdArrRho.LON=NC_Read2Dvariable(eFile, xName + "_rho");
  GrdArr.GrdArrRho.LAT=NC_Read2Dvariable(eFile, yName + "_rho");
  GrdArr.GrdArrRho.DEP=NC_Read2Dvariable(eFile, "h");
  GrdArr.GrdArrRho.HaveDEP=true;
  GrdArr.GrdArrRho.ANG=NC_Read2Dvariable(eFile, "angle");
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho=GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> eMSK_rho_double=NC_Read2Dvariable(eFile, "mask_rho");
  GrdArr.GrdArrRho.MSK=ConvertMatrix(eMSK_rho_double, fConv);
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  // U
  GrdArr.GrdArrU.LON=NC_Read2Dvariable(eFile, xName + "_u");
  GrdArr.GrdArrU.LAT=NC_Read2Dvariable(eFile, yName + "_u");
  int eta_u=GrdArr.GrdArrU.LON.rows();
  int xi_u=GrdArr.GrdArrU.LON.cols();
  MyMatrix<int> MSKu(eta_u, xi_u);
  MyMatrix<double> DEPu(eta_u, xi_u);
  MyMatrix<double> ANGu(eta_u, xi_u);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u; j++) {
      DEPu(i, j)=
	(GrdArr.GrdArrRho.DEP(i, j)+
	 GrdArr.GrdArrRho.DEP(i, j+1) )/double(2);
      ANGu(i, j)=
	(GrdArr.GrdArrRho.ANG(i, j)+
	 GrdArr.GrdArrRho.ANG(i, j+1) )/double(2);
      MSKu(i, j)=
	GrdArr.GrdArrRho.MSK(i, j)*
	GrdArr.GrdArrRho.MSK(i, j+1);
    }
  GrdArr.GrdArrU.MSK=MSKu;
  GrdArr.GrdArrU.DEP=DEPu;
  GrdArr.GrdArrU.HaveDEP=true;
  GrdArr.GrdArrU.ANG=ANGu;
  InitializeIdxJdxWet(GrdArr.GrdArrU);
  // V
  GrdArr.GrdArrV.LON=NC_Read2Dvariable(eFile, xName + "_v");
  GrdArr.GrdArrV.LAT=NC_Read2Dvariable(eFile, yName + "_v");
  int eta_v=GrdArr.GrdArrV.LON.rows();
  int xi_v=GrdArr.GrdArrV.LON.cols();
  MyMatrix<int> MSKv(eta_v, xi_v);
  MyMatrix<double> DEPv(eta_v, xi_v);
  MyMatrix<double> ANGv(eta_v, xi_v);
  for (int i=0; i<eta_v; i++)
    for (int j=0; j<xi_v; j++) {
      DEPv(i, j)=
	(GrdArr.GrdArrRho.DEP(i, j)+

	 GrdArr.GrdArrRho.DEP(i+1,j) )/double(2);
      ANGv(i, j)=
	(GrdArr.GrdArrRho.ANG(i, j)+
	 GrdArr.GrdArrRho.ANG(i+1, j) )/double(2);
      MSKv(i, j)=
	GrdArr.GrdArrRho.MSK(i, j)*
	GrdArr.GrdArrRho.MSK(i+1, j);
    }
  GrdArr.GrdArrV.MSK=MSKv;
  GrdArr.GrdArrV.DEP=DEPv;
  GrdArr.GrdArrV.HaveDEP=true;
  GrdArr.GrdArrV.ANG=ANGv;
  InitializeIdxJdxWet(GrdArr.GrdArrV);
  // PSI
  GrdArr.GrdArrPsi.LON=NC_Read2Dvariable(eFile, xName + "_psi");
  GrdArr.GrdArrPsi.LAT=NC_Read2Dvariable(eFile, yName + "_psi");
  int eta_psi=GrdArr.GrdArrPsi.LON.rows();
  int xi_psi=GrdArr.GrdArrPsi.LON.cols();
  MyMatrix<int> MSKp(eta_psi, xi_psi);
  MyMatrix<double> DEPp(eta_psi, xi_psi);
  MyMatrix<double> ANGp(eta_psi, xi_psi);
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++) {
      DEPp(i, j)=
	(GrdArr.GrdArrRho.DEP(i  ,j+1)+
	 GrdArr.GrdArrRho.DEP(i+1,j+1)+
	 GrdArr.GrdArrRho.DEP(i  ,j  )+
	 GrdArr.GrdArrRho.DEP(i+1,j  ))/double(4);
      ANGp(i, j)=
	(GrdArr.GrdArrRho.ANG(i  ,j+1)+
	 GrdArr.GrdArrRho.ANG(i+1,j+1)+
	 GrdArr.GrdArrRho.ANG(i  ,j  )+
	 GrdArr.GrdArrRho.ANG(i+1,j  ))/double(4);
      MSKp(i, j)=
	 GrdArr.GrdArrRho.MSK(i  ,j+1)*
	 GrdArr.GrdArrRho.MSK(i+1,j+1)*
	 GrdArr.GrdArrRho.MSK(i  ,j  )*
	 GrdArr.GrdArrRho.MSK(i+1,j  );
    }
  GrdArr.GrdArrPsi.MSK=MSKp;
  GrdArr.GrdArrPsi.DEP=DEPp;
  GrdArr.GrdArrPsi.HaveDEP=true;
  GrdArr.GrdArrPsi.ANG=ANGp;
  InitializeIdxJdxWet(GrdArr.GrdArrPsi);
  std::cerr << "The ROMS grid has been read\n";
  bool PrintKeyInformation=true;
  if (PrintKeyInformation) {
    auto ThePrint=[&](int const& i1, int const& i2) -> void {
      std::cerr << "(i,j) = (" << i1 << "," << i2 << ") lon=" << GrdArr.GrdArrRho.LON(i1,i2) << " lat=" << GrdArr.GrdArrRho.LAT(i1,i2) << "\n";
    };
    ThePrint(0,0);
    ThePrint(0,xi_rho-1);
    ThePrint(eta_rho-1,xi_rho-1);
    ThePrint(eta_rho-1,0);
  }
  return GrdArr;
}

CoordGridArrayFD GRID_ExtendedPsiThi(CoordGridArrayFD const& RecRho, CoordGridArrayFD const& RecU, CoordGridArrayFD const& RecV, CoordGridArrayFD const& RecPsi)
{
  int eta_rho=RecRho.LON.rows();
  int xi_rho=RecRho.LON.cols();
  int eta_psi = RecPsi.LON.rows();
  int xi_psi = RecPsi.LON.cols();
  int eta_u = RecU.LON.rows();
  //  int xi_u = RecU.LON.cols();
  //  int eta_v = RecV.LON.rows();
  int xi_v = RecV.LON.cols();
  int eta_psi2 = eta_rho + 1;
  int xi_psi2 = xi_rho + 1;
  // We have eta_u = eta_rho    ;  eta_v = eta_rho-1  ;  eta_psi = eta_rho-1
  // We have xi_u  = xi_rho -1  ;  xi_v  = xi_rho     ;  xi_psi  = xi_rho-1
  auto CompField2=[&](MyMatrix<double> const& FieldPsi, MyMatrix<double> const& FieldU, MyMatrix<double> const& FieldV) -> MyMatrix<double> {
    MyMatrix<double> FieldPsi2(eta_psi2, xi_psi2);
    for (int i=0; i<eta_psi; i++)
      for (int j=0; j<xi_psi; j++)
	FieldPsi2(i+1,j+1) = FieldPsi(i,j);
    for (int i=0; i<eta_psi; i++) {
      FieldPsi2(i+1,0)         = 2*FieldV(i,0)      - FieldPsi(i,0);
      FieldPsi2(i+1,xi_psi2-1) = 2*FieldV(i,xi_v-1) - FieldPsi(i,xi_psi-1);
    }
    for (int j=0; j<xi_psi; j++) {
      FieldPsi2(0,j+1)          = 2*FieldU(0,j)       - FieldPsi(0,j);
      FieldPsi2(eta_psi2-1,j+1) = 2*FieldU(eta_u-1,j) - FieldPsi(eta_psi-1,j);
    }
    FieldPsi2(0,0)         = 2*FieldPsi2(1,0)         - FieldPsi2(2,0);
    FieldPsi2(0,xi_psi2-1) = 2*FieldPsi2(1,xi_psi2-1) - FieldPsi2(2,xi_psi2-1);
    //
    FieldPsi2(eta_psi2-1,0)         = 2*FieldPsi2(eta_psi2-2,0)         - FieldPsi2(eta_psi2-3,0);
    FieldPsi2(eta_psi2-1,xi_psi2-1) = 2*FieldPsi2(eta_psi2-2,xi_psi2-1) - FieldPsi2(eta_psi2-3,xi_psi2-1);
    return FieldPsi2;
  };
  MyMatrix<double> LON_psi2 = CompField2(RecPsi.LON, RecU.LON, RecV.LON);
  MyMatrix<double> LAT_psi2 = CompField2(RecPsi.LAT, RecU.LAT, RecV.LAT);
  MyMatrix<double> DEP_psi2 = CompField2(RecPsi.DEP, RecU.DEP, RecV.DEP);
  MyMatrix<double> ANG_psi2 = CompField2(RecPsi.ANG, RecU.ANG, RecV.ANG);
  //
  MyMatrix<int> MSK_psi2=ZeroMatrix<int>(eta_psi2, xi_psi2);
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++)
      MSK_psi2(i+1,j+1) = RecPsi.MSK(i,j);
  for (int i=0; i<eta_psi; i++) {
    MSK_psi2(i+1,0) = MSK_psi2(i+1,1);
    MSK_psi2(i+1,xi_psi2-1) = MSK_psi2(i+1, xi_psi2-2);
  }
  for (int j=0; j<xi_psi; j++) {
    MSK_psi2(0,j+1) = MSK_psi2(1,j+1);
    MSK_psi2(eta_psi2-1,j+1) = MSK_psi2(eta_psi2-2, j+1);
  }
  if (MSK_psi2(1,0) && MSK_psi2(0,1))
    MSK_psi2(0,0) = 1;
  if (MSK_psi2(eta_psi2-2,0) && MSK_psi2(eta_psi2-1,1))
    MSK_psi2(eta_psi2-1,0) = 1;
  if (MSK_psi2(1,xi_psi2-1) && MSK_psi2(0,xi_psi2-2))
    MSK_psi2(0,xi_psi2-1) = 1;
  if (MSK_psi2(eta_psi2-2,xi_psi2-1) && MSK_psi2(eta_psi2-1,xi_psi2-2))
    MSK_psi2(eta_psi2-1,xi_psi2-1) = 1;
  //
  CoordGridArrayFD RecPsi2;
  RecPsi2.MSK = MSK_psi2;
  RecPsi2.LON = LON_psi2;
  RecPsi2.LAT = LAT_psi2;
  RecPsi2.DEP = DEP_psi2;
  RecPsi2.ANG = ANG_psi2;
  RecPsi2.HaveDEP=true;
  InitializeIdxJdxWet(RecPsi2);
  return RecPsi2;
}


MyMatrix<double> GetMatrixRadiusROMS(GridArray const& GrdArr)
{
  CoordGridArrayFD RecPsi2=GRID_ExtendedPsiThi(GrdArr.GrdArrRho, GrdArr.GrdArrU, GrdArr.GrdArrV, GrdArr.GrdArrPsi);
  MyMatrix<int> MatDir(4,2);
  MatDir(0,0)=0;
  MatDir(0,1)=0;
  MatDir(1,0)=1;
  MatDir(1,1)=0;
  MatDir(2,0)=1;
  MatDir(2,1)=1;
  MatDir(3,0)=0;
  MatDir(3,1)=1;
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> MatRadius(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double lon1=GrdArr.GrdArrRho.LON(i,j);
      double lat1=GrdArr.GrdArrRho.LON(i,j);
      double MinDist=45555555;
      for (int u=0; u<4; u++) {
	int i2=i + MatDir(u,0);
	int j2=j + MatDir(u,1);
	double lon2=RecPsi2.LON(i2,j2);
	double lat2=RecPsi2.LON(i2,j2);
	double dx=lon1 - lon2;
	double dy=lat1 - lat2;
	double dist=sqrt(dx*dx + dy*dy);
	if (dist < MinDist)
	  MinDist=dist;
      }
      MatRadius(i,j) = MinDist;
    }
  return MatRadius;
}




GridArray NC_ReadWrfGridFile(std::string const& eFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="WRF";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  // Rho part of the arrays
  MyMatrix<double> LON=NC_Read2Dvariable(eFile, "XLONG");
  MyMatrix<double> LAT=NC_Read2Dvariable(eFile, "XLAT");
  GrdArr.GrdArrRho.LON=LON;
  GrdArr.GrdArrRho.LAT=LAT;
  GrdArr.GrdArrRho.HaveDEP=false;
  GrdArr.GrdArrRho.ANG=CreateAngleMatrix(LON, LAT);
  int eta_rho=LON.rows();
  int xi_rho=LON.cols();
  MyMatrix<int> MSK;
  MSK.setConstant(eta_rho, xi_rho, 1);
  GrdArr.GrdArrRho.MSK=MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  return GrdArr;
}



// It must be a TEM file
GridArray NC_ReadNemoGridFile(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in NC_ReadNemoGridFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile = " << eFile << "\n";
    throw TerminalException{1};
  }
  GridArray GrdArr;
  GrdArr.ModelName="NEMO";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  // Rho part of the arrays
  MyVector<double> lon1d=NC_Read1Dvariable(eFile, "lon");
  MyVector<double> lat1d=NC_Read1Dvariable(eFile, "lat");
  MyVector<double> dep1d_pre=NC_Read1Dvariable(eFile, "depth");
  int nbLon=lon1d.size();
  int nbLat=lat1d.size();
  int nbDep=dep1d_pre.size();
  MyVector<double> dep1d(nbDep);
  // We want index 0 to be deepest and index nbDep-1 to be near surface
  for (int iDep=0; iDep<nbDep; iDep++)
    dep1d(nbDep-1-iDep) = dep1d_pre(iDep);
  for (int iDep=0; iDep<nbDep; iDep++)
    std::cerr << "iDep=" << iDep << " dep=" << dep1d(iDep) << "\n";
  MyMatrix<double> LON(nbLat, nbLon);
  MyMatrix<double> LAT(nbLat, nbLon);
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      LON(i,j) = lon1d(j);
      LAT(i,j) = lat1d(i);
    }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  if (dataFile.isNull()) {
    std::cerr << "Error while opening dataFile\n";
    throw TerminalException{1};
  }
  netCDF::NcVar data=dataFile.getVar("votemper");
  if (dataFile.isNull()) {
    std::cerr << "Error while reading sossheig\n";
    throw TerminalException{1};
  }
  MyVector<int> StatusFill=NC_ReadVariable_StatusFill_data(data);
  MyVector<double> VarFill=NC_ReadVariable_data(data);
  std::cerr << "|StatusFill|=" << StatusFill.size() << " min/max=" << StatusFill.minCoeff() << " / " << StatusFill.maxCoeff() << " sum=" << StatusFill.sum() << "\n";
  std::vector<size_t> ListDim = NC_ReadVariable_listdim(data);
  int nbTime=ListDim[0];
  int s_vert=ListDim[1];
  int eta=ListDim[2];
  int xi=ListDim[3];
  std::cerr << "nbTime=" << nbTime << " nbDep=" << nbDep << " eta=" << eta << " xi=" << xi << "\n";
  if (eta != nbLat || xi != nbLon || s_vert != nbDep) {
    std::cerr << "eta=" << eta << " nbLat=" << nbLat << "\n";
    std::cerr << "xi=" << xi << " nbLon=" << nbLon << "\n";
    std::cerr << "s_vert=" << s_vert << " nbDep=" << nbDep << "\n";
    std::cerr << "Inconsistency between array sizes. Logical error\n";
    throw TerminalException{1};
  }
  //
  // Computing MSK and DEP
  //
  MyMatrix<int> MSK(nbLat, nbLon);
  MyMatrix<double> DEP(nbLat, nbLon);
  Eigen::Tensor<int,4> StatusTens(nbTime, nbDep, nbLat, nbLon);
  Eigen::Tensor<double,4> VarTens(nbTime, nbDep, nbLat, nbLon);
  MyMatrix<int> StatusSum=ZeroMatrix<int>(nbLat, nbLon);
  int idx=0;
  for (int iTime=0; iTime<nbTime; iTime++)
    for (int iDep=0; iDep<nbDep; iDep++)
      for (int i=0; i<nbLat; i++)
	for (int j=0; j<nbLon; j++) {
	  StatusTens(iTime, nbDep -1 - iDep, i, j) = StatusFill(idx);
	  VarTens(iTime, nbDep - 1 - iDep, i, j) = VarFill(idx);
	  StatusSum(i,j) += StatusFill(idx);
	  idx++;
	}
  bool CoherencyCheck=true;
  if (CoherencyCheck) {
    for (int iTime=0; iTime<nbTime; iTime++)
      for (int i=0; i<nbLat; i++)
	for (int j=0; j<nbLon; j++) {
	  for (int iDep=1; iDep<nbDep; iDep++) {
	    if (StatusTens(iTime, iDep-1,i,j) == 0 && StatusTens(iTime, iDep,i,j) == 1) {
	      std::cerr << "Found inconsistency at i=" << i << " j=" << j << " iTime = " << iTime << " iDep=" << iDep << "\n";
	    }
	  }
	}
    std::cerr << "After coherency checks\n";
  }
  int ValLand = nbTime * nbDep;
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      int eVal=1;
      if (StatusSum(i,j) == ValLand)
	eVal=0;
      MSK(i,j)=eVal;
    }
  int iTimeRef=0;
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      if (StatusSum(i,j) == 48) {
	double lat=lat1d(i);
	double lon=lon1d(j);
	std::cerr << "i=" << i << " j=" << j << " lon/lat=" << lon << " / " << lat << "\n";
	for (int iDep=0; iDep<nbDep; iDep++)
	  std::cerr << "iDep=" << iDep << " dep=" << dep1d(iDep) << " status=" << StatusTens(iTimeRef,iDep,i,j) << "\n";
      }
    }
  std::cerr << "StatusFill min/max=" << StatusFill.minCoeff() << " / " << StatusFill.maxCoeff() << "\n";
  std::cerr << "StatusSum  min/max=" << StatusSum.minCoeff() << " / " << StatusSum.maxCoeff() << "\n";
  std::cerr << "NEMO MSK min / max / sum=" << MSK.minCoeff() << " / " << MSK.maxCoeff() << " / " << MSK.sum() << "\n";
  std::cerr << "nbLat=" << nbLat << " nbLon=" << nbLon << "\n";
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      double eDep=0;
      if (MSK(i,j) == 1) {
	for (int iDep=1; iDep<nbDep; iDep++)  {
	  if (StatusTens(iTimeRef, iDep-1, i, j) == 1 && StatusTens(iTimeRef, iDep, i, j) == 0) {
	    double dep1=dep1d(iDep-1); // rock
	    double dep2=dep1d(iDep);  // sea
	    eDep = (dep1 + dep2)/double(2);
	  }
	}
	if (eDep == 0) {
	  int sumStatus=0;
	  for (int iDep=0; iDep<nbDep; iDep++)
	    sumStatus += StatusTens(iTimeRef, iDep, i, j);
	  if (sumStatus != nbDep) {
	    double lon=lon1d(j);
	    double lat=lat1d(i);
	    std::cerr << "i=" << i << " j=" << j << " lon=" << lon << " lat=" << lat << "\n";
	    for (int iDep=0; iDep<nbDep; iDep++)
	      std::cerr << "iDep=" << iDep << " status=" << StatusTens(iTimeRef, iDep, i, j) << " dep=" << dep1d(iDep) << "\n";
	    std::cerr << "sumStatus=" << sumStatus << " nbDep=" << nbDep << "\n";
	    std::cerr << "sumStatus is not equal to nbDep\n";
	    throw TerminalException{1};
	  }
	  eDep = dep1d(nbDep-1);
	}
      }
      /*
      if (eDep > 1500) {
	double lon=lon1d(j);
	double lat=lat1d(i);
	std::cerr << "i=" << i << "j=" << j << " lon=" << lon << " lat=" << lat << " eDep=" << eDep << "\n";
	} */
      DEP(i,j) = eDep;
    }
  std::cerr << "DEP min/max=" << DEP.minCoeff() << " / " << DEP.maxCoeff() << "\n";
  double sumTEM=0;
  for (int i=0; i<nbLat; i++)
    for (int j=0; j<nbLon; j++) {
      if (MSK(i,j) == 1) {
	sumTEM += VarTens(iTimeRef, nbDep-1, i, j);
      }
    }
  std::cerr << "sumTEM = " << sumTEM << "\n";
  GrdArr.ARVD.N=nbDep;
  GrdArr.ARVD.IsAssigned=true;
  GrdArr.ARVD.ListZ_r = -dep1d;
  GrdArr.ARVD.Zcoordinate=true;
  GrdArr.ARVD.ModelName="NEMO";
  //
  // More standard assignation
  //
  GrdArr.GrdArrRho.LON=LON;
  GrdArr.GrdArrRho.LAT=LAT;
  GrdArr.GrdArrRho.DEP=DEP;
  GrdArr.GrdArrRho.HaveDEP=true;
  GrdArr.GrdArrRho.ANG=CreateAngleMatrix(LON, LAT);
  GrdArr.GrdArrRho.MSK=MSK;
  InitializeIdxJdxWet(GrdArr.GrdArrRho);
  return GrdArr;
}




GridArray NC_ReadCosmoWamStructGridFile(std::string const& eFile, std::string const& postfix)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  if (postfix == "atm")
    GrdArr.ModelName="COSMO";
  else
    GrdArr.ModelName="WAM";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  // Rho part of the arrays
  std::string LONstr="LON_" + postfix;
  //  std::cerr << "Before reading " << LONstr << "\n";
  GrdArr.GrdArrRho.LON=NC_Read2Dvariable(eFile, LONstr);
  std::string LATstr = "LAT_" + postfix;
  //  std::cerr << "Before reading " << LATstr << "\n";
  GrdArr.GrdArrRho.LAT=NC_Read2Dvariable(eFile, LATstr);
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho=GrdArr.GrdArrRho.LON.cols();
  //  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
  // The bathymetry if available
  std::string DEPstr = "DEP_" + postfix;
  if (NC_IsVar(eFile, DEPstr) ) {
    GrdArr.GrdArrRho.DEP=NC_Read2Dvariable(eFile, DEPstr);
    GrdArr.GrdArrRho.HaveDEP=true;
  }
  else {
    GrdArr.GrdArrRho.DEP=ZeroMatrix<double>(eta_rho, xi_rho);
    GrdArr.GrdArrRho.HaveDEP=false;
  }
  // The mask if available
  std::string MSKstr="MSK_" + postfix;
  MyMatrix<double> MSK_double(eta_rho, xi_rho);
  if (NC_IsVar(eFile, DEPstr) ) {
    MSK_double=NC_Read2Dvariable(eFile, MSKstr);
  }
  else {
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	MSK_double(i,j)=double(1);
      }
  }
  MyMatrix<int> MSK_int(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      MSK_int(i,j)=int(MSK_double(i,j));
    }
  GrdArr.GrdArrRho.MSK=MSK_int;
  // The angle if available
  //  std::cerr << "Before reading angle\n";
  if (postfix == "atm")
    GrdArr.GrdArrRho.ANG=NC_Read2Dvariable(eFile, "ANG_atm");
  else
    GrdArr.GrdArrRho.ANG=CreateAngleMatrix(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
  // U / V / PSI we do not need a priori
  return GrdArr;
}



GridArray NC_ReadSCHISM_sflux_grid(std::string const& eFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="SCHISM_SFLUX";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  // Rho part of the arrays
  std::string LONstr="lon";
  GrdArr.GrdArrRho.LON=NC_Read2Dvariable(eFile, LONstr);
  std::string LATstr="lat";
  GrdArr.GrdArrRho.LAT=NC_Read2Dvariable(eFile, LATstr);
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho=GrdArr.GrdArrRho.LON.cols();
  GrdArr.GrdArrRho.HaveDEP=false;
  //
  MyMatrix<int> MSK_int(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++)
      MSK_int(i,j)=1;
  GrdArr.GrdArrRho.MSK=MSK_int;
  GrdArr.GrdArrRho.ANG=CreateAngleMatrix(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
  // U / V / PSI we do not need a priori
  return GrdArr;
}



MyMatrix<int> NC_ReadElements(std::string const& eFile, std::string const& eStr)
{
  MyMatrix<int> INE=NC_Read2Dvariable_int(eFile, eStr);
  int nbRow=INE.rows();
  int nbCol=INE.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      int IP=INE(iRow, iCol);
      INE(iRow, iCol)=IP-1;
    }
  return INE;
}


GridArray NC_ReadWamGridFile(std::string const& eFile)
{
  MyVector<int> eVectLLUNSTR=NC_Read1Dvariable_int(eFile, "LLUNSTR");
  int LLUNSTR=eVectLLUNSTR(0);
  //  std::cerr << "LLUNSTR=" << LLUNSTR << "\n";
  if (LLUNSTR == 0)
    return NC_ReadCosmoWamStructGridFile(eFile, "wav");
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="WAM";
  GrdArr.IsFE=1;
  GrdArr.IsSpherical=true;
  GrdArr.L_IndexSelect=false;
  //
  GrdArr.INE=NC_ReadElements(eFile, "ele");
  //  std::cerr << "NC_ReadWam, step 1\n";
  MyVector<double> LON=NC_Read1Dvariable(eFile, "LON_wav");
  //  std::cerr << "NC_ReadWam, step 2\n";
  MyVector<double> LAT=NC_Read1Dvariable(eFile, "LAT_wav");
  //  std::cerr << "NC_ReadWam, step 3\n";
  MyVector<double> DEP=NC_Read1Dvariable(eFile, "DEP_wav");
  //  std::cerr << "NC_ReadWam, step 4\n";
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  //  std::cerr << "NC_ReadWam, step 5\n";
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=DEP(iPoint);
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
  }
  GrdArr.GrdArrRho.LON=LONarr;
  GrdArr.GrdArrRho.LAT=LATarr;
  GrdArr.GrdArrRho.DEP=DEParr;
  GrdArr.GrdArrRho.ANG=ANGarr;
  GrdArr.GrdArrRho.MSK=MSKarr;
  return GrdArr;
}





GridArray WWM_ReadGridFile_netcdf(std::string const& GridFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.IsFE=1;
  GrdArr.L_IndexSelect=false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_netcdf\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  GrdArr.INE=NC_ReadElements(GridFile, "ele");
  bool IsVar_X=NC_IsVar(GridFile, "x");
  bool IsVar_Y=NC_IsVar(GridFile, "y");
  bool IsVar_lon=NC_IsVar(GridFile, "lon");
  bool IsVar_lat=NC_IsVar(GridFile, "lat");
  if (IsVar_X != IsVar_Y) {
    std::cerr << "IsVar_X/Y=" << IsVar_X << " / " << IsVar_Y << "\n";
    std::cerr << "They should be identical\n";
    throw TerminalException{1};
  }
  if (IsVar_lon != IsVar_lat) {
    std::cerr << "IsVar_lon/lat=" << IsVar_lon << " / " << IsVar_lat << "\n";
    std::cerr << "They should be identical\n";
    throw TerminalException{1};
  }
  if ((IsVar_X && IsVar_lon) || (!IsVar_X && !IsVar_lon)) {
    std::cerr << "IsVar_X/lon=" << IsVar_X << " / " << IsVar_lon << "\n";
    std::cerr << "One and exactly one option should be selected\n";
    throw TerminalException{1};
  }
  bool LSPHE=IsVar_lon;
  std::string Xname, Yname;
  if (LSPHE) {
    Xname="lon";
    Yname="lat";
    GrdArr.IsSpherical=true;
  }
  else {
    Xname="x";
    Yname="y";
    GrdArr.IsSpherical=false;
  }
  MyVector<double> LON=NC_Read1Dvariable(GridFile, Xname);
  MyVector<double> LAT=NC_Read1Dvariable(GridFile, Yname);
  MyVector<double> DEP=NC_Read1Dvariable(GridFile, "depth");
  MyVector<double> IOBP=NC_Read1Dvariable(GridFile, "IOBP");
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  MyVector<int> IOBParr(nbPoint);
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=DEP(iPoint);
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
    IOBParr(iPoint)=int(IOBP(iPoint));
  }
  GrdArr.GrdArrRho.LON=LONarr;
  GrdArr.GrdArrRho.LAT=LATarr;
  GrdArr.GrdArrRho.DEP=DEParr;
  GrdArr.GrdArrRho.ANG=ANGarr;
  GrdArr.GrdArrRho.MSK=MSKarr;
  GrdArr.IOBP=IOBParr;
  return GrdArr;
}



// IOBP should be the WWM type of IOBP
// Possible values are:
// --- 0: normal point
// --- 1: Island
// --- 2: Dirichlet
// --- 3: Neumann
// --- 4: New mixed kind
void WriteGridFile_msh(std::string const& GridFile, GridArray const& GrdArr)
{
  if (GrdArr.IsFE == 0) {
    std::cerr << "We need the grid to be finite element\n";
    throw TerminalException{1};
  }

  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  int np_total=GrdArr.GrdArrRho.LON.rows();
  std::vector<int> IPbound;
  std::vector<int> IPisland;
  if (GrdArr.IOBP.size() != np_total) {
    std::cerr << "We have |GrdArr.IOBP|=" << GrdArr.IOBP.size() << "\n";
    std::cerr << "when it should be np_total=" << np_total << "\n";
    std::cerr << "Most likely you forgot to provide for boundary file\n";
    throw TerminalException{1};
  }
  for (int i=0; i<np_total; i++) {
    int eVal=GrdArr.IOBP(i);
    if (eVal == 2) {
      IPbound.push_back(i);
    }
    if (eVal == 1 || eVal == 3 || eVal == 4) {
      IPisland.push_back(i);
    }
  }
  int nbDirichlet=IPbound.size();
  int nbIsland=IPisland.size();
  std::vector<int> ACTIVE(nbDirichlet,1);
  os << "$MeshFormat\n";
  os << "2 0 8\n";
  os << "$EndMeshFormat\n";
  os << "$Nodes\n";
  os << np_total << "\n";
  for (int i=0; i<np_total; i++) {
    int IP=i+1;
    double eXP=GrdArr.GrdArrRho.LON(i,0);
    double eYP=GrdArr.GrdArrRho.LAT(i,0);
    double eDEP=GrdArr.GrdArrRho.DEP(i,0);
    os << IP << " " << eXP << " " << eYP << " " << eDEP << "\n";
  }
  //
  // Writing the elements
  //
  int ne_total=GrdArr.INE.rows();
  int nbEle=ne_total + nbDirichlet + nbIsland;
  os << "$EndNodes\n";
  os << "$Elements\n";
  os << nbEle << "\n";
  //
  // Write boundary
  //
  int ie2=0;
  for (int i=0; i<nbDirichlet; i++) {
    ie2++;
    int eAct=ACTIVE[i];
    os << ie2 << " 15 2 " << eAct << " 0 " << IPbound[i] << "\n";
  }
  //
  // Write island
  //
  for (int i=0; i<nbIsland; i++) {
    ie2++;
    int ip=i+1;
    int eIPisland=IPisland[i] + 1;
    os << ie2 << " 15 2 0 " << ip << " " << eIPisland << "\n";
  }
  //
  // Write gmsh elements
  //
  for (int ie=0; ie<ne_total; ie++) {
    ie2++;
    int ieRel=ie+1;
    int ip1=GrdArr.INE(ie,0) + 1;
    int ip2=GrdArr.INE(ie,1) + 1;
    int ip3=GrdArr.INE(ie,2) + 1;
    os << ie2 << " 2 3 0 " << ieRel << " 0 " << ip1 << " " << ip2 << " " << ip3 << "\n";
  }
  os << "$EndElements\n";
}






void WriteWWMboundaryGR3(std::string const& BndFile, GridArray const& GrdArr)
{
  int nbVert=GrdArr.IOBP.size();
  int nbNode=GrdArr.GrdArrRho.LON.size();
  int nbTrig=GrdArr.INE.rows();
  if (nbVert != nbNode) {
    std::cerr << "We have nbVert not equal to nbNode\n";
    std::cerr << "nbNode = " << nbNode << "\n";
    std::cerr << "nbVert = " << nbVert << "\n";
    throw TerminalException{1};
  }
  std::ofstream os(BndFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << "boundary file\n";
  os << nbTrig << " " << nbNode << "\n";
  for (int iVert=0; iVert<nbVert; iVert++) {
    double eLon=GrdArr.GrdArrRho.LON(iVert,0);
    double eLat=GrdArr.GrdArrRho.LAT(iVert,0);
    int iPos=iVert+1;
    int eIOBP=GrdArr.IOBP(iVert);
    os << iPos << " " << eLon << " " << eLat << " " << eIOBP << "\n";
  }
}


MyVector<int> WWM_ReadBoundFile_gr3(std::string const& BoundFile)
{
  if (!IsExistingFile(BoundFile)) {
    std::cerr << "Error in WWM_ReadBoundFile_gr3\n";
    std::cerr << "Missing BoundFile=" << BoundFile << "\n";
    throw TerminalException{1};
  }
  std::ifstream IN(BoundFile);
  std::string line;
  std::getline(IN, line);
  int mne, mnp;
  IN >> mne;
  IN >> mnp;
  MyVector<int> eVect(mnp);
  for (int i=0; i<mnp; i++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != i+1) {
      std::cerr << "Inconsistency at this level\n";
      throw TerminalException{1};
    }
    int eIOBP=int(ZPDTMP);
    eVect(i)=eIOBP;
  }
  return eVect;
}


/* This strategy fails if the domain has a width of less than 90m */
bool GuessIsSpherical(GridArray const& GrdArr)
{
  bool IsSpherical=true;
  double LONmax = GrdArr.GrdArrRho.LON.maxCoeff();
  double LONmin = GrdArr.GrdArrRho.LON.minCoeff();
  double deltaLON=LONmax - LONmin;
  if (deltaLON > 360)
    IsSpherical=false;
  double LATmax = GrdArr.GrdArrRho.LAT.maxCoeff();
  double LATmin = GrdArr.GrdArrRho.LAT.minCoeff();
  if (LATmax > 90 || LATmin < -90)
    IsSpherical=false;
  return IsSpherical;
}




GridArray WWM_ReadGridFile_msh(std::string const& GridFile)
{
  GridArray GrdArr;
  std::ifstream is(GridFile);
  std::string line;
  std::getline(is, line);
  if (line != "$MeshFormat") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "First line should be $MeshFormat\n";
    throw TerminalException{1};
  }
  std::getline(is, line);
  std::cerr << "GMSH version number / fileType / dataSize = " << line << "\n";
  //
  std::getline(is, line);
  if (line != "$EndMeshFormat") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "Line should be $EndMeshFormat\n";
    throw TerminalException{1};
  }
  //
  auto WaitForString=[&](std::string const& strSearch) -> void {
    while(true) {
      try {
	std::string strRead;
	std::getline(is, strRead);
	if (strRead == strSearch)
	  return;
      }
      catch (...) {
	std::cerr << "Error in data reading\n";
	throw TerminalException{1};
      }
    }
  };
  //
  WaitForString("$Nodes");
  //
  std::getline(is, line);
  int mnp;
  std::istringstream(line) >> mnp;
  std::cerr << "mnp=" << mnp << "\n";
  MyMatrix<double> LON(mnp,1), LAT(mnp,1), DEP(mnp,1);
  for (int ip=0; ip<mnp; ip++) {
    int idx;
    double eLon, eLat, eDep;
    std::getline(is, line);
    std::vector<std::string> LStr=STRING_Split(line, " ");
    if (LStr.size() < 4) {
      std::cerr << "|LStr|=" << LStr.size() << "\n";
      std::cerr << "line=" << line << "\n";
      std::cerr << "Unfortunately LStr is too small\n";
      throw TerminalException{1};
    }
    std::istringstream(LStr[0]) >> idx;
    std::istringstream(LStr[1]) >> eLon;
    std::istringstream(LStr[2]) >> eLat;
    std::istringstream(LStr[3]) >> eDep;
    if (idx != ip+1) {
      std::cerr << "idx=" << idx << " eLon=" << eLon << " eLat=" << eLat << " eDep=" << eDep << "\n";
      std::cerr << "idx=" << idx << " ip=" << ip << "\n";
      std::cerr << "Error in the indices\n";
      throw TerminalException{1};
    }
    LON(ip,0)=eLon;
    LAT(ip,0)=eLat;
    DEP(ip,0)=eDep;
  }
  GrdArr.GrdArrRho.LON=LON;
  GrdArr.GrdArrRho.LAT=LAT;
  GrdArr.GrdArrRho.DEP=DEP;
  //
  std::getline(is,line);
  if (line != "$EndNodes") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "Line should be $EndNodes\n";
    throw TerminalException{1};
  }
  //
  WaitForString("$Elements");
  //
  std::getline(is, line);
  int nbElemTot;
  std::istringstream(line) >> nbElemTot;
  int nb15=0;
  int nb1=0;
  int nb2=0;
  int nb3=0;
  std::cerr << "nbElemTot=" << nbElemTot << "\n";
  std::vector<int> List15, List1, List2, List3;
  for (int ie=0; ie<nbElemTot; ie++) {
    std::getline(is, line);
    std::vector<int> LInt=STRING_Split_Int(line, " ");
    int len=LInt.size();
    if (len < 2) {
      std::cerr << "len is not large enough. len=" << len << "\n";
      throw TerminalException{1};
    }
    int idx=LInt[0];
    int elm_type=LInt[1];
    if (idx != ie+1) {
      std::cerr << "idx=" << idx << " ie=" << ie << "\n";
      std::cerr << "Error in the indices\n";
      throw TerminalException{1};
    }
    if (elm_type != 15 && elm_type != 1 && elm_type != 2 && elm_type != 3) {
      std::cerr << "elm_type=" << elm_type << "\n";
      std::cerr << "Only values allowed are 15, 1, 2 and 3\n";
      throw TerminalException{1};
    }
    if (elm_type == 15) {
      int a4=LInt[len-1];
      List15.push_back(a4-1);
      nb15++;
    }
    if (elm_type == 1) {
      int a1=LInt[len-2];
      int a2=LInt[len-1];
      List1.push_back(a1-1);
      List1.push_back(a2-1);
      nb1++;
    }
    if (elm_type == 2) {
      int a1=LInt[len-3];
      int a2=LInt[len-2];
      int a3=LInt[len-1];
      List2.push_back(a1-1);
      List2.push_back(a2-1);
      List2.push_back(a3-1);
      nb2++;
    }
    if (elm_type == 3) {
      int a1=LInt[len-4];
      int a2=LInt[len-3];
      int a3=LInt[len-2];
      int a4=LInt[len-1];
      List3.push_back(a1-1);
      List3.push_back(a2-1);
      List3.push_back(a3-1);
      List3.push_back(a4-1);
      nb3++;
    }
  }
  if (nb2 > 0 && nb3 > 0) {
    std::cerr << "nb2=" << nb2 << " nb3=" << nb3 << "\n";
    std::cerr << "Right now we cannot mix triangles and quadrangles\n";
    throw TerminalException{1};
  }
  std::getline(is, line);
  if (line != "$EndElements") {
    std::cerr << "line=" << line << "\n";
    std::cerr << "We reach an error here. Should be $EndElements\n";
    throw TerminalException{1};
  }
  if (nb2 > 0) {
    MyMatrix<int> INE(nb2,3);
    int idx=0;
    for (int ie=0; ie<nb2; ie++) {
      INE(ie,0)=List2[idx];
      idx++;
      INE(ie,1)=List2[idx];
      idx++;
      INE(ie,2)=List2[idx];
      idx++;
    }
    GrdArr.INE=INE;
  }
  if (nb3 > 0) {
    MyMatrix<int> INE(nb3,4);
    int idx=0;
    for (int ie=0; ie<nb3; ie++) {
      INE(ie,0)=List3[idx];
      idx++;
      INE(ie,1)=List3[idx];
      idx++;
      INE(ie,2)=List3[idx];
      idx++;
      INE(ie,3)=List3[idx];
      idx++;
    }
    GrdArr.INE=INE;
  }
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ARVD.Zcoordinate=false;
  GrdArr.GrdArrU.HaveDEP=false;
  GrdArr.GrdArrV.HaveDEP=false;
  GrdArr.GrdArrPsi.HaveDEP=false;
  return GrdArr;
}





GridArray WWM_ReadGridFile_gr3(std::string const& GridFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.IsFE=1;
  GrdArr.L_IndexSelect=false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_gr3\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::ifstream IN(GridFile);
  // read first line
  std::string line;
  std::getline(IN, line);
  std::cerr << "line=" << line << "\n";
  //
  int mne, mnp;
  IN >> mne;
  IN >> mnp;
  std::cerr << "mne=" << mne << " mnp=" << mnp << "\n";
  GrdArr.INE=MyMatrix<int>(mne,3);
  GrdArr.GrdArrRho.LON=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
  for (int iP=0; iP<mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    //    std::cerr << "iP=" << iP << " XYZ=" << XPDTMP << " " << YPDTMP << " " << ZPDTMP << "\n";
    GrdArr.GrdArrRho.LON(iP,0)=XPDTMP;
    GrdArr.GrdArrRho.LAT(iP,0)=YPDTMP;
    GrdArr.GrdArrRho.DEP(iP,0)=ZPDTMP;
    GrdArr.GrdArrRho.ANG(iP,0)=0;
    GrdArr.GrdArrRho.MSK(iP,0)=1;
  }
  for (int iE=0; iE<mne; iE++) {
    int KTMP, LTMP, ip1, ip2, ip3;
    IN >> KTMP >> LTMP >> ip1 >> ip2 >> ip3;
    //    std::cerr << "iE=" << iE << " IP123=" << ip1 << " " << ip2 << " " << ip3 << "\n";
    GrdArr.INE(iE,0)=ip1 - 1;
    GrdArr.INE(iE,1)=ip2 - 1;
    GrdArr.INE(iE,2)=ip3 - 1;
  }
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}


MyVector<int> WWM_ReadBoundFile_DAT(std::string const& BoundFile)
{
  if (!IsExistingFile(BoundFile)) {
    std::cerr << "Error in WWM_ReadBoundFile_DAT\n";
    std::cerr << "Missing BoundFile=" << BoundFile << "\n";
    throw TerminalException{1};
  }
  std::ifstream IN(BoundFile);
  std::string line;
  for (int i=0; i<2; i++)
    std::getline(IN, line);
  int ITMP, JTMP;
  IN >> ITMP;
  std::getline(IN, line);
  IN >> JTMP;
  int mnp=ITMP + JTMP;
  for (int i=0; i<7; i++)
    std::getline(IN, line);
  MyVector<int> eVect(mnp);
  for (int i=0; i<mnp; i++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != i+1) {
      std::cerr << "Inconsistency error\n";
      throw TerminalException{1};
    }
    int eIOBP=int(ZPDTMP);
    eVect(i)=eIOBP;
  }
  return eVect;
}


GridArray WWM_ReadGridFile_obj(std::string const& GridFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.IsFE=1;
  GrdArr.L_IndexSelect=false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_obj\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::string line;
  std::ifstream IN(GridFile);
  std::getline(IN, line);
  //  std::cerr << "line=" << line << "\n";
  std::vector<std::vector<double>> ListVect;
  std::vector<std::vector<int>> ListTrig;
  while(true) {
    if (IN.eof() == 1)
      break;
    std::getline(IN, line);
    if (line.size() == 0)
      break;
    //    std::cerr << "line=" << line << " |line|=" << line.size() << "\n";
    std::istringstream is(line);
    std::string eChar;
    is >> eChar;
    bool IsDone=false;
    if (eChar == "v") {
      IsDone=true;
      double lon, lat, dep;
      is >> lon >> lat >> dep;
      ListVect.push_back({lon, lat, dep});
    }
    if (eChar == "f") {
      IsDone=true;
      int i1, i2, i3;
      is >> i1 >> i2 >> i3;
      ListTrig.push_back({i1-1,i2-1,i3-1});
    }
    if (!IsDone) {
      std::cerr << "Error while reading data\n";
      std::cerr << "eChar=" << eChar << "\n";
      throw TerminalException{1};
    }
  }
  int mnp=ListVect.size();
  GrdArr.GrdArrRho.LON=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
  for (int i=0; i<mnp; i++) {
    GrdArr.GrdArrRho.LON(i,0)=ListVect[i][0];
    GrdArr.GrdArrRho.LAT(i,0)=ListVect[i][1];
    GrdArr.GrdArrRho.DEP(i,0)=ListVect[i][2];
    GrdArr.GrdArrRho.ANG(i,0)=0;
    GrdArr.GrdArrRho.MSK(i,0)=1;
  }
  int mne=ListTrig.size();
  GrdArr.INE=MyMatrix<int>(mne,3);
  for (int i=0; i<mne; i++)
    for (int j=0; j<3; j++)
      GrdArr.INE(i,j)=ListTrig[i][j];
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}








GridArray WWM_ReadGridFile_DAT(std::string const& GridFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.IsFE=1;
  GrdArr.L_IndexSelect=false;
  //
  if (!IsExistingFile(GridFile)) {
    std::cerr << "Error in WWM_ReadGridFile_DAT\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    throw TerminalException{1};
  }
  std::string line;
  int ITMP, JTMP;
  std::ifstream IN(GridFile);
  for (int i=0; i<2; i++)
    std::getline(IN, line);
  std::getline(IN, line);
  std::istringstream(line) >> ITMP;
  std::cerr << "ITMP=" << ITMP << "\n";
  std::getline(IN, line);
  std::getline(IN, line);
  std::istringstream(line) >> JTMP;
  std::cerr << "JTMP=" << JTMP << "\n";
  int mnp=ITMP + JTMP;
  std::cerr << "mnp=" << mnp << "\n";
  for (int i=0; i<7; i++)
    std::getline(IN, line);
  GrdArr.GrdArrRho.LON=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
  GrdArr.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
  for (int iP=0; iP<mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    std::getline(IN, line);
    std::istringstream(line) >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    if (KTMP != iP) {
      std::cerr << "KTMP=" << KTMP << " iP=" << iP << "\n";
      std::cerr << "Inconsistency in the values\n";
      throw TerminalException{1};
    }
    GrdArr.GrdArrRho.LON(iP,0)=XPDTMP;
    GrdArr.GrdArrRho.LAT(iP,0)=YPDTMP;
    GrdArr.GrdArrRho.DEP(iP,0)=ZPDTMP;
    GrdArr.GrdArrRho.ANG(iP,0)=0;
    GrdArr.GrdArrRho.MSK(iP,0)=1;
  }
  for (int i=0; i<2; i++)
    std::getline(IN, line);
  std::getline(IN, line);
  int mne;
  std::istringstream(line) >> mne;
  GrdArr.INE=MyMatrix<int>(mne,3);
  for (int i=0; i<3; i++)
    std::getline(IN, line);
  for (int iE=0; iE<mne; iE++) {
    int KTMP, LTMP, ip1, ip2, ip3;
    std::getline(IN, line);
    std::istringstream(line) >> ip1 >> ip2 >> ip3 >> KTMP >> LTMP;
    if (LTMP != iE || KTMP != 0) {
      std::cerr << "LTMP=" << LTMP << " iE=" << iE << " (should be equal)\n";
      std::cerr << "KTMP=" << KTMP << " (should be 0)\n";
      std::cerr << "Inconsistent values in the .dat file\n";
      throw TerminalException{1};
    }
    GrdArr.INE(iE,0)=ip1;
    GrdArr.INE(iE,1)=ip2;
    GrdArr.INE(iE,2)=ip3;
  }
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}























GridArray NC_ReadWW3_GridFile(std::string const& eFile)
{
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="WWM";
  GrdArr.IsFE=1;
  GrdArr.L_IndexSelect=false;
  std::cerr << "NC_ReadWW3_GRidFile\n";
  //
  GrdArr.INE=NC_ReadElements(eFile, "tri");
  MyVector<double> LON=NC_Read1Dvariable(eFile, "longitude");
  MyVector<double> LAT=NC_Read1Dvariable(eFile, "latitude");
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=0;
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
  }
  GrdArr.GrdArrRho.LON=LONarr;
  GrdArr.GrdArrRho.LAT=LATarr;
  GrdArr.GrdArrRho.DEP=DEParr;
  GrdArr.GrdArrRho.ANG=ANGarr;
  GrdArr.GrdArrRho.MSK=MSKarr;
  GrdArr.ModelName="WW3";
  GrdArr.IsSpherical = GuessIsSpherical(GrdArr);
  return GrdArr;
}




int TheSignFct(double const& eVal)
{
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}



MyMatrix<double> get_angle_corr_rho(MyMatrix<double> const& LON_rho, MyMatrix<double> const& LAT_rho)
{
  std::string spheroid="wgs84";
  double A=-1, B=-1, E=-1;
  if (spheroid == "sph") {
    A = 6371000.0;
    B = A;
    E = sqrt(A*A-B*B)/A;
  }
  if (spheroid == "cla") {
    A = 6378206.4E0;
    B = 6356583.8E0;
    E= sqrt(A*A-B*B)/A;
  }
  if (spheroid == "iau") {
    A = 6378160.e0;
    B = 6356774.516E0;
    E = sqrt(A*A-B*B)/A;
  }
  if (spheroid == "wgs84") {
    A = 6378137.;
    E = 0.081819191;
    double alpha=A*E;
    B = sqrt(A*A - alpha*alpha);
  }
  double eps= E*E/(1-E*E);
  int eta_rho=LON_rho.rows();
  int xi_rho=LON_rho.cols();
  int eta_u=eta_rho;
  int xi_u=xi_rho-1;
  MyMatrix<double> LONrad_u(eta_u, xi_u);
  MyMatrix<double> LATrad_u(eta_u, xi_u);
  double pi=3.1415926535;
  double eFact=pi/double(360);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u; j++) {
      double eLON=(LON_rho(i,j) + LON_rho(i,j+1))*eFact;
      double eLAT=(LAT_rho(i,j) + LAT_rho(i,j+1))*eFact;
      if (eLAT == 0)
	eLAT=eps;
      LONrad_u(i,j)=eLON;
      LATrad_u(i,j)=eLAT;
    }
  MyMatrix<double> azim(eta_u, xi_u-1);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u-1; j++) {
      double PHI1=LATrad_u(i, j);
      double XLAM1=LONrad_u(i, j);
      double PHI2=LATrad_u(i, j+1);
      double XLAM2=LONrad_u(i, j+1);
      if (PHI1 == PHI2)
	PHI2=PHI2 + 1e-14;
      if (XLAM1 == XLAM2)
	XLAM2=XLAM2 + 1e-14;
      //
      double EsPHI1=E*sin(PHI1);
      double EsPHI2=E*sin(PHI2);
      double xnu1=A/sqrt(1 - EsPHI1*EsPHI1);
      double xnu2=A/sqrt(1 - EsPHI2*EsPHI2);
      double TPSI2=(1-E*E)*tan(PHI2) + E*E*xnu1*sin(PHI1)/(xnu2*cos(PHI2));
      double DLAM=XLAM2 - XLAM1;
      double CTA12=(cos(PHI1)*TPSI2 - sin(PHI1)*cos(DLAM))/sin(DLAM);
      double DLAM2=DLAM;
      if (DLAM2 >= pi)
	DLAM2=DLAM2 - 2*pi;
      if (DLAM2 <= -pi)
	DLAM2=DLAM2 + 2*pi;
      double eAzim=atan(1/CTA12);
      if (eAzim < -pi)
	eAzim += 2*pi;
      if (eAzim > pi)
	eAzim += - 2*pi;
      if (TheSignFct(eAzim) != TheSignFct(DLAM2)) {
	eAzim += pi*double(TheSignFct(-eAzim));
      }
      azim(i,j)=eAzim;
    }
  MyMatrix<double> angle(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=1; j<xi_u; j++) {
      double eAzim=azim(i, j-1);
      double eAngle=(pi/double(2)) - eAzim;
      angle(i, j)=eAngle;
    }
  for (int i=0; i<eta_rho; i++) {
    double eAngle=angle(i,1);
    angle(i,0)=eAngle;
    eAngle=angle(i, xi_u-1);
    angle(i, xi_u)=eAngle;
  }
  return angle;
}


void DifferenceLonRenormalize(double & Lon)
{
  if (Lon > 180)
    Lon=Lon - 360;
  if (Lon < -180)
    Lon=Lon + 360;
}



GridArray CFONE_GRID_ARRAY(std::string const& GridFile)
{
  MyVector<double> LonArr=NC_Read1Dvariable(GridFile, "lon");
  MyVector<double> LatArr=NC_Read1Dvariable(GridFile, "lat");
  int nbLon=LonArr.size();
  int nbLat=LatArr.size();
  MyMatrix<double> LON(nbLon, nbLat);
  MyMatrix<double> LAT(nbLon, nbLat);
  MyMatrix<int> MSK(nbLon, nbLat);
  MyMatrix<double> DEP(nbLon, nbLat);
  MyMatrix<double> ANG(nbLon, nbLat);
  for (int iLon=0; iLon<nbLon; iLon++)
    for (int iLat=0; iLat<nbLat; iLat++) {
      LON(iLon, iLat)=LonArr(iLon);
      LAT(iLon, iLat)=LatArr(iLat);
      MSK(iLon, iLat)=1;
      DEP(iLon, iLat)=0;
      ANG(iLon, iLat)=0;
    }
  CoordGridArrayFD GrdArrRho;
  GrdArrRho.LON=LON;
  GrdArrRho.LAT=LAT;
  GrdArrRho.MSK=MSK;
  GrdArrRho.DEP=DEP;
  GrdArrRho.ANG=ANG;
  //
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="CFconvention";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  GrdArr.GrdArrRho=GrdArrRho;
  return GrdArr;
}

MyMatrix<double> MatrixSubsample(MyMatrix<double> const& F, int const& splitRow, int const& splitCol)
{
  int nbRow=F.rows();
  int nbCol=F.cols();
  auto GetListIdx=[](int const& nbPos, int const& splitPos) -> std::vector<int> {
    int nbPosRed=nbPos / splitPos;
    double multCoef=double(nbPos-1) / double(nbPosRed - 1);
    std::vector<int> ListIdx(nbPosRed);
    for (int i=0; i<nbPosRed; i++) {
      double xPos=double(i) * multCoef;
      int iPos=int(round(xPos));
      int ePos=std::max(0, std::min(nbPos-1, iPos));
      ListIdx[i]=ePos;
    }
    ListIdx[0]=0;
    ListIdx[nbPosRed-1]=nbPos-1;
    return ListIdx;
  };
  std::vector<int> ListIdxRow=GetListIdx(nbRow, splitRow);
  std::vector<int> ListIdxCol=GetListIdx(nbCol, splitCol);
  int nbRowRed=ListIdxRow.size();
  int nbColRed=ListIdxCol.size();
  MyMatrix<double> Fred(nbRowRed, nbColRed);
  for (int iRowRed=0; iRowRed<nbRowRed; iRowRed++)
    for (int iColRed=0; iColRed<nbColRed; iColRed++) {
      int iRow=ListIdxRow[iRowRed];
      int iCol=ListIdxCol[iColRed];
      Fred(iRowRed,iColRed) = F(iRow,iCol);
    }
  return Fred;
}







GridArray CURVILINEAR_GRID_ARRAY(MyMatrix<double> const& LON, MyMatrix<double> const& LAT)
{
  if (LON.rows() != LAT.rows() || LON.cols() != LAT.cols()) {
    std::cerr << "LON and LAT should have same size\n";
    throw TerminalException{1};
  }
  int nbRow=LON.rows();
  int nbCol=LON.cols();
  MyMatrix<int> MSK(nbRow, nbCol);
  MyMatrix<double> DEP(nbRow, nbCol);
  MyMatrix<double> ANG(nbRow, nbCol);
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      MSK(iRow, iCol)=1;
      DEP(iRow, iCol)=0;
      ANG(iRow, iCol)=0;
    }
  CoordGridArrayFD GrdArrRho;
  GrdArrRho.LON=LON;
  GrdArrRho.LAT=LAT;
  GrdArrRho.MSK=MSK;
  GrdArrRho.DEP=DEP;
  GrdArrRho.ANG=ANG;
  //
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="CURVILINEAR";
  GrdArr.IsFE=0;
  GrdArr.IsSpherical=true;
  GrdArr.GrdArrRho=GrdArrRho;
  return GrdArr;
}



GridArray RECTANGULAR_GRID_ARRAY(QuadArray const& eQuad, int const& nbSplitLon, int const& nbSplitLat)
{
  double MinLon = eQuad.MinLon;
  double MinLat = eQuad.MinLat;
  double MaxLon = eQuad.MaxLon;
  double MaxLat = eQuad.MaxLat;
  double deltaLon=(MaxLon - MinLon)/double(nbSplitLon-1);
  double deltaLat=(MaxLat - MinLat)/double(nbSplitLat-1);
  std::cerr << "nbSplitLon=" << nbSplitLon << " nbSplitLat=" << nbSplitLat << "\n";
  MyMatrix<double> LON(nbSplitLon, nbSplitLat);
  MyMatrix<double> LAT(nbSplitLon, nbSplitLat);
  for (int iLon=0; iLon<nbSplitLon; iLon++)
    for (int iLat=0; iLat<nbSplitLat; iLat++) {
      LON(iLon, iLat)=MinLon + deltaLon*iLon;
      LAT(iLon, iLat)=MinLat + deltaLat*iLat;
    }
  GridArray GrdArr=CURVILINEAR_GRID_ARRAY(LON, LAT);
  GrdArr.ARVD.IsAssigned=false;
  GrdArr.ModelName="RECTANGULAR";
  return GrdArr;
}






void CutWorldMap(GridArray & GrdArr)
{
  // We cut at -180 - eps
  double eps=1e-8;
  int nbPoint=GrdArr.GrdArrRho.LON.rows();
  double LonSplit=0;
  while(true) {
    double MinDist=2400;
    for (int iPoint=0; iPoint<nbPoint; iPoint++) {
      double eLon=GrdArr.GrdArrRho.LON(iPoint,0);
      if (eLon > 0)
	eLon -= 360;
      LonSplit= - 180 - eps;
      double dist=fabs(eLon - LonSplit);
      if (dist < MinDist)
	MinDist=dist;
    }
    std::cerr << "eps=" << eps << " MinDist=" << MinDist << "\n";
    if (MinDist > eps/2)
      break;
    eps *= 2;
  }
  int nbTrig=GrdArr.INE.rows();
  std::vector<int> ListStatus(nbTrig);
  int SumStatus=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    double eLon1=GrdArr.GrdArrRho.LON(i1,0);
    double eLon2=GrdArr.GrdArrRho.LON(i2,0);
    double eLon3=GrdArr.GrdArrRho.LON(i3,0);
    eLon1 -= LonSplit;
    eLon2 -= LonSplit;
    eLon3 -= LonSplit;
    DifferenceLonRenormalize(eLon1);
    DifferenceLonRenormalize(eLon2);
    DifferenceLonRenormalize(eLon3);
    int eStatus=1;
    double UpperLimit=90;
    if (fabs(eLon1) < UpperLimit &&
	fabs(eLon2) < UpperLimit &&
	fabs(eLon3) < UpperLimit) {
      double eProd12=eLon1*eLon2;
      double eProd23=eLon2*eLon3;
      double eProd31=eLon3*eLon1;
      if (eProd12 < 0 || eProd23 < 0 || eProd31 < 0)
	eStatus=0;
    }
    ListStatus[iTrig]=eStatus;
    SumStatus += eStatus;
  }
  std::cerr << "SumStatus = " << SumStatus << "   nbTrig = " << nbTrig << "\n";
  int nbTrigNew=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1)
      nbTrigNew++;
  MyMatrix<int> INEnew(nbTrigNew,3);
  int iTrigNew=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1) {
      int i1=GrdArr.INE(iTrig,0);
      int i2=GrdArr.INE(iTrig,1);
      int i3=GrdArr.INE(iTrig,2);
      INEnew(iTrigNew,0)=i1;
      INEnew(iTrigNew,1)=i2;
      INEnew(iTrigNew,2)=i3;
      iTrigNew++;
    }
  GrdArr.INE=INEnew;
}









void CUT_HigherLatitude(GridArray & GrdArr, double MinLatCut, double MaxLatCut)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  std::vector<int> ListStatus(mnp);
  std::vector<int> Index(mnp);
  std::vector<int> RevIndex(mnp);
  int iNodeNew=0;
  std::vector<int> I_IndexSelectOld;
  if (GrdArr.L_IndexSelect) {
    I_IndexSelectOld=GrdArr.I_IndexSelect;
  }
  else {
    for (int i=0; i<mnp; i++)
      I_IndexSelectOld.push_back(i);
  }
  std::vector<int> I_IndexSelect;
  for (int iNode=0; iNode<mnp; iNode++) {
    double eLat=GrdArr.GrdArrRho.LAT(iNode,0);
    if (MinLatCut < eLat && eLat < MaxLatCut) {
      ListStatus[iNode]=1;
      Index[iNode]=iNodeNew;
      RevIndex[iNodeNew]=iNode;
      int iNodeMain=I_IndexSelectOld[iNode];
      I_IndexSelect.push_back(iNodeMain);
      iNodeNew++;
    }
    else {
      ListStatus[iNode]=0;
    }
  }
  int nbNodeNew=iNodeNew;
  MyMatrix<double> LONnew(nbNodeNew,1);
  MyMatrix<double> LATnew(nbNodeNew,1);
  MyMatrix<double> DEPnew(nbNodeNew,1);
  MyMatrix<double> ANGnew(nbNodeNew,1);
  MyMatrix<int>    MSKnew(nbNodeNew,1);
  for (int iNodeNewB=0; iNodeNewB<nbNodeNew; iNodeNewB++) {
    int iNode=RevIndex[iNodeNewB];
    double eLon=GrdArr.GrdArrRho.LON(iNode,0);
    double eLat=GrdArr.GrdArrRho.LAT(iNode,0);
    double eDep=GrdArr.GrdArrRho.DEP(iNode,0);
    double eAng=GrdArr.GrdArrRho.ANG(iNode,0);
    int eMsk=GrdArr.GrdArrRho.MSK(iNode,0);
    LONnew(iNodeNewB,0)=eLon;
    LATnew(iNodeNewB,0)=eLat;
    DEPnew(iNodeNewB,0)=eDep;
    ANGnew(iNodeNewB,0)=eAng;
    MSKnew(iNodeNewB,0)=eMsk;
  }
  GrdArr.GrdArrRho.LON=LONnew;
  GrdArr.GrdArrRho.LAT=LATnew;
  GrdArr.GrdArrRho.DEP=DEPnew;
  GrdArr.GrdArrRho.ANG=ANGnew;
  GrdArr.GrdArrRho.MSK=MSKnew;
  int nbTrigNew=0;
  for (int iTrig=0; iTrig<mne; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1)
      nbTrigNew++;
  }
  MyMatrix<int> INEnew(nbTrigNew, 3);
  int iTrigNew=0;
  for (int iTrig=0; iTrig<mne; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1) {
      INEnew(iTrigNew,0)=Index[i1];
      INEnew(iTrigNew,1)=Index[i2];
      INEnew(iTrigNew,2)=Index[i3];
      iTrigNew++;
    }
  }
  GrdArr.INE=INEnew;
  GrdArr.L_IndexSelect = true;
  GrdArr.I_IndexSelect = I_IndexSelect;
}





double GetGridSpacing(GridArray const& GrdArr)
{
  int IsFE=GrdArr.IsFE;
  double SumDistKM=0;
  int SumNb=0;
  if (IsFE == 1) {
    int nbEle=GrdArr.INE.rows();
    for (int iEle=0; iEle<nbEle; iEle++)
      for (int i=0; i<3; i++) {
	int j=NextIdx(3,i);
	int iNode1=GrdArr.INE(iEle,i);
	int iNode2=GrdArr.INE(iEle,j);
	double eLon1=GrdArr.GrdArrRho.LON(iNode1,0);
	double eLat1=GrdArr.GrdArrRho.LAT(iNode1,0);
	double eLon2=GrdArr.GrdArrRho.LON(iNode2,0);
	double eLat2=GrdArr.GrdArrRho.LAT(iNode2,0);
	double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	SumDistKM += DistKM;
	SumNb += 1;
      }
  }
  else {
    int nbRow=GrdArr.GrdArrRho.LON.rows();
    int nbCol=GrdArr.GrdArrRho.LON.cols();
    for (int iRow=0; iRow<nbRow; iRow++) {
      for (int iCol=1; iCol<nbCol; iCol++) {
	double eLon1=GrdArr.GrdArrRho.LON(iRow, iCol);
	double eLat1=GrdArr.GrdArrRho.LAT(iRow, iCol);
	double eLon2=GrdArr.GrdArrRho.LON(iRow, iCol-1);
	double eLat2=GrdArr.GrdArrRho.LAT(iRow, iCol-1);
	double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	SumDistKM += DistKM;
	SumNb += 1;
      }
    }
    for (int iRow=1; iRow<nbRow; iRow++) {
      for (int iCol=0; iCol<nbCol; iCol++) {
	double eLon1=GrdArr.GrdArrRho.LON(iRow, iCol);
	double eLat1=GrdArr.GrdArrRho.LAT(iRow, iCol);
	double eLon2=GrdArr.GrdArrRho.LON(iRow-1, iCol);
	double eLat2=GrdArr.GrdArrRho.LAT(iRow-1, iCol);
	double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
	SumDistKM += DistKM;
	SumNb++;
      }
    }
  }
  double avgDistKM=SumDistKM / double(SumNb);
  return avgDistKM;
}





struct GridSymbolic {
public:  
  GridSymbolic() {
    Sphericity="unset";
    CutWorldMap=false;
    HigherLatitudeCut=false;
    MinLatCut=0;
    MaxLatCut=0;
    MinLat=0;
    MaxLat=0;
    MinLon=0;
    MaxLon=0;
    deltaKM=0;
  }
  GridSymbolic(std::string _Sphericity, bool _CutWorldMap, bool _HigherLatitudeCut, double _MinLatCut, double _MaxLatCut, double _MinLat, double _MaxLat, double _MinLon, double _MaxLon, double _deltaKM)
  {
    Sphericity=_Sphericity;
    CutWorldMap=_CutWorldMap;
    HigherLatitudeCut=_HigherLatitudeCut;
    MinLatCut = _MinLatCut;
    MaxLatCut = _MaxLatCut;
    MinLat    = _MinLat;
    MaxLat    = _MaxLat;
    MinLon    = _MinLon;
    MaxLon    = _MaxLon;
    deltaKM   = _deltaKM;
  }
  std::string Sphericity;
  bool CutWorldMap;
  bool HigherLatitudeCut;
  double MinLatCut;
  double MaxLatCut;
  double MinLat;
  double MaxLat;
  double MinLon;
  double MaxLon;
  double deltaKM;
};


struct TripleModelDesc {
  std::string ModelName;
  std::string GridFile;
  std::string BoundFile;
  std::string HisPrefix;
  GridSymbolic RecGridSymb;
};




std::string GET_GRID_FILE(TripleModelDesc const& eTriple)
{
  std::string PreModelName = eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  std::string HisPrefix=eTriple.HisPrefix;
  if (eModelName == "RECTANGULAR")
    return "irrelevant";
  if (eModelName == "COSMO")
    return HisPrefix + "0001.nc";
  if (eModelName == "WAM")
    return HisPrefix + "0001.nc";
  if (eModelName == "ROMS" || eModelName == "ROMS_IVICA")
    return eTriple.GridFile;
  if (eModelName == "WWM" || eModelName == "WWM_DAILY" || eModelName == "UNRUNOFF")
    return eTriple.GridFile;
  if (eModelName == "SCHISM_SFLUX")
    return eTriple.GridFile;
  if (eModelName == "WRF")
    return HisPrefix + "0001.nc"; // but maybe should be something else
  if (eModelName == "SCHISM_NETCDF_OUT")
    return eTriple.GridFile;
  if (eModelName == "WW3") {
    std::string ThePrefix=HisPrefix + "*";
    std::vector<std::string> ListFile=ls_operation(ThePrefix);
    return ListFile[0];
  }
  if (eModelName == "NEMO") {
    std::vector<std::string> ListFile=FILE_DirectoryFilesSpecificExtension(HisPrefix, "nc");
    for (auto & eFile : ListFile) {
      std::vector<std::string> ListStr=STRING_Split(eFile, "tem");
      if (ListStr.size() == 2)
	return eFile;
    }
    std::cerr << "We failed to find the matching file with a tem in the title\n";
    throw TerminalException{1};
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
    std::vector<std::string> ListFile=FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      throw TerminalException{1};
    }
    return ListFile[0];
  }
  if (eModelName == "GRIB_WAM_FORT30") {
    std::string eFile=HisPrefix;
    if (!IsExistingFile(eFile)) {
      std::cerr << "The file eFile = " << eFile << " is missing\n";
      std::cerr << "It serves as grid and should be put in HisPrefix\n";
      throw TerminalException{1};
    }
    return eFile;
  }
  std::cerr << "Error in GET_GRID_FILE\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  throw TerminalException{1};
}





void WriteUnstructuredGrid_GR3(std::string const& GridFile, GridArray const& GrdArr)
{
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << GridFile << "\n";
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  os << mne << " " << mnp << "\n";
  for (int ip=0; ip<mnp; ip++) {
    double lon=GrdArr.GrdArrRho.LON(ip);
    double lat=GrdArr.GrdArrRho.LAT(ip);
    double dep=GrdArr.GrdArrRho.DEP(ip);
    int idx=ip+1;
    os << idx << " " << lon << " " << lat << " " << dep << "\n";
  }
  for (int ie=0; ie<mne; ie++) {
    int ip1=GrdArr.INE(ie,0)+1;
    int ip2=GrdArr.INE(ie,1)+1;
    int ip3=GrdArr.INE(ie,2)+1;
    int idx=ie+1;
    os << idx << " 3 " << ip1 << " " << ip2 << " " << ip3 << "\n";
  }
}



GridArray WWM_ReadGridFile_Ricchiuto_grd(std::string const& GridFile)
{
  std::ifstream IN(GridFile);
  GridArray GrdArr;
  GrdArr.ARVD.IsAssigned=false;
  //
  std::string FirstLine;
  std::getline(IN, FirstLine);
  //
  std::string lineDATA;
  std::getline(IN, lineDATA);
  std::vector<std::string> LStr=STRING_Split(lineDATA, " ");
  int mne, mnp;
  std::istringstream(LStr[1]) >> mne;
  std::istringstream(LStr[2]) >> mnp;
  //
  MyMatrix<int> INE(mne,3);
  for (int ie=0; ie<mne; ie++) {
    std::string lineINE;
    std::getline(IN, lineINE);
    std::vector<std::string> LStrB=STRING_Split(lineINE, " ");
    for (int i=0; i<3; i++) {
      int IP;
      std::istringstream(LStrB[i]) >> IP;
      INE(ie,i)=IP;
    }
  }
  GrdArr.INE=INE;
  GrdArr.IsFE=true;
  //
  MyMatrix<double> LON(mnp,1), LAT(mnp,1), DEP(mnp,1);
  for (int ip=0; ip<mnp; ip++) {
    std::string lineXYD;
    std::getline(IN, lineXYD);
    std::vector<std::string> LStrC=STRING_Split(lineXYD, " ");
    double eLon, eLat, eDep;
    std::istringstream(LStrC[0]) >> eLon;
    std::istringstream(LStrC[1]) >> eLat;
    std::istringstream(LStrC[2]) >> eDep;
    LON(ip,0)=eLon;
    LAT(ip,0)=eLat;
    DEP(ip,0)=eDep;
  }
  GrdArr.GrdArrRho.LON=LON;
  GrdArr.GrdArrRho.LAT=LAT;
  GrdArr.GrdArrRho.DEP=DEP;
  return GrdArr;
}

void WriteUnstructuredGrid_Ricchiuto_GRD(std::string const& GridFile, GridArray const& GrdArr)
{
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  os << GridFile << "\n";
  int NN=GrdArr.GrdArrRho.LON.rows();
  int NE=GrdArr.INE.rows();
  int NBF=0;
  os << " 2 " << NE << " " << NN << " " << NBF << "\n";
  for (int ie=0; ie<NE; ie++) {
    int ip1=GrdArr.INE(ie,0);
    int ip2=GrdArr.INE(ie,1);
    int ip3=GrdArr.INE(ie,2);
    os << ip1 << " " << ip2 << " " << ip3 << "\n";
  }
  for (int ip=0; ip<NN; ip++) {
    double lon=GrdArr.GrdArrRho.LON(ip);
    double lat=GrdArr.GrdArrRho.LAT(ip);
    double dep=GrdArr.GrdArrRho.DEP(ip);
    os << lon << " " << lat << " " << dep << "\n";
  }
  for (int ib=0; ib<NBF; ib++) {
    /* Write something here */
  }
}






void WriteUnstructuredGrid_DAT(std::string const& GridFile, GridArray const& GrdArr)
{
  std::ofstream os(GridFile);
  os << std::fixed;
  os << std::setprecision(9);
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  os << "C system.dat, erzeugt von xf am sometime\n";
  os << "C Anzahl der Randknoten:\n";
  os << 0 << "\n";
  os << "C Anzahl der Gebietsknoten:\n";
  os << mnp << "\n";
  os << "C Koordinaten und Skalarwerte der Knoten\n";
  os << "C --------------------------------------\n";
  os << "C Zuerst die Randknoten  (Anzahl s.o.),\n";
  os << "C dann die Gebietsknoten (Anzahl s.o.).\n";
  os << "C ------------+-------------+-------------+---------------\n";
  os << "C     Nr.     |  x-Koord.   |   y-Koord.  |   Skalarwert\n";
  os << "C ------------+-------------+-------------+---------------\n";
  for (int ip=0; ip<mnp; ip++) {
    double lon=GrdArr.GrdArrRho.LON(ip);
    double lat=GrdArr.GrdArrRho.LAT(ip);
    double dep=GrdArr.GrdArrRho.DEP(ip);
    int idx=ip;
    os << idx << " " << lon << " " << lat << " " << dep << "\n";
  }
  os << "C ------------------------------------------------------------\n";
  os << "C Anzahl der Elemente:\n";
  os << mne << "\n";
  os << "C Elementverzeichnis\n";
  os << "C ------------------------------------------------------------\n";
  os << "C    Knoten i  Knoten j  Knoten k   Kennung     Nr.\n";
  for (int ie=0; ie<mne; ie++) {
    int ip1=GrdArr.INE(ie,0);
    int ip2=GrdArr.INE(ie,1);
    int ip3=GrdArr.INE(ie,2);
    int idx=ie;
    os << ip1 << " " << ip2 << " " << ip3 << " 0 " << idx << "\n";
  }
  os << "C ------------------------------------------------------------\n";
}

void WriteUnstructuredGrid_NC(std::string const& GridFile, GridArray const& GrdArr)
{
  netCDF::NcFile dataFile(GridFile, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  netCDF::NcDim eDimOne   = dataFile.addDim("one", 1);
  netCDF::NcDim eDimTwo   = dataFile.addDim("two", 2);
  netCDF::NcDim eDimThree = dataFile.addDim("three", 3);
  netCDF::NcDim eDimMnp   = dataFile.addDim("mnp", mnp);
  netCDF::NcDim eDimMne   = dataFile.addDim("mne", mne);
  std::vector<std::string> ListDimLL{"mnp"};
  std::vector<std::string> ListDimINE{"mne", "three"};
  //
  std::string strLON, strLAT;
  if (GrdArr.IsSpherical) {
    strLON="lon";
    strLAT="lat";
  }
  else {
    strLON="x";
    strLAT="y";
  }
  auto putVARmnp=[&](std::string const& name, MyMatrix<double> const& VAR) -> void {
    netCDF::NcVar eVar = dataFile.addVar(name, "double", ListDimLL);
    double *A;
    A=new double[mnp];
    for (int ip=0; ip<mnp; ip++)
      A[ip] = VAR(ip,0);
    eVar.putVar(A);
    delete [] A;
  };
  putVARmnp(strLON, GrdArr.GrdArrRho.LON);
  putVARmnp(strLAT, GrdArr.GrdArrRho.LAT);
  putVARmnp("depth", GrdArr.GrdArrRho.DEP);
  //
  /*
  bool *Abool;
  std::vector<std::string> ListDimSPHE{"one"};
  netCDF::NcVar eVarSphe = dataFile.addVar("LSPHE", "bool", ListDimSPHE);
  Abool=new bool[1];
  Abool[0]=GrdArr.IsSpherical;
  eVarSphe.putVar(Abool);
  delete [] Abool;*/
  //
  if (GrdArr.IOBP.size() != mnp) {
    std::cerr << "We have |GrdArr.IOBP.size()|=" << GrdArr.IOBP.size() << "\n";
    std::cerr << "Which is distinct from mnp=" << mnp << "\n";
    std::cerr << "Most likely you forgot to put the input file\n";
    throw TerminalException{1};
  }
  int *A;
  netCDF::NcVar eVar = dataFile.addVar("IOBP", "int", ListDimLL);
  A=new int[mnp];
  for (int ip=0; ip<mnp; ip++)
    A[ip] = GrdArr.IOBP(ip,0);
  eVar.putVar(A);
  delete [] A;
  //
  int *Aine;
  netCDF::NcVar eVarIne = dataFile.addVar("ele", "int", ListDimINE);
  Aine=new int[mne*3];
  int idx=0;
  for (int ie=0; ie<mne; ie++)
    for (int j=0; j<3; j++) {
      Aine[idx] = GrdArr.INE(ie,j) + 1;
      idx++;
    }
  eVarIne.putVar(Aine);
  delete [] Aine;
  //
  MyMatrix<int> LEdge=GetEdgeSet(GrdArr.INE, mnp);
  int nbEdge=LEdge.rows();
  netCDF::NcDim eDimNbEdge = dataFile.addDim("nbedge", nbEdge);
  std::vector<std::string> ListDimEDGE{"nbedge", "two"};
  int *Aedge;
  Aedge=new int[nbEdge * 2];
  netCDF::NcVar eVarEdge = dataFile.addVar("edges", "int", ListDimEDGE);
  int idx2=0;
  for (int iedge=0; iedge<nbEdge; iedge++)
    for (int j=0; j<2; j++) {
      Aedge[idx2] = LEdge(iedge,j) + 1;
      idx2++;
    }
  eVarEdge.putVar(Aedge);
  delete [] Aedge;
}




void WriteUnstructuredGrid(std::string const& GridFile, GridArray const& GrdArr)
{
  std::string eExtension=FILE_GetExtension(GridFile);
  std::cerr << "WriteUnstructuredGrid  with  eExtension=" << eExtension << "\n";
  if (eExtension == "gr3") {
    WriteUnstructuredGrid_GR3(GridFile, GrdArr);
    return;
  }
  if (eExtension == "dat") {
    WriteUnstructuredGrid_DAT(GridFile, GrdArr);
    return;
  }
  if (eExtension == "grd") {
    WriteUnstructuredGrid_Ricchiuto_GRD(GridFile, GrdArr);
    return;
  }
  if (eExtension == "nc") {
    WriteUnstructuredGrid_NC(GridFile, GrdArr);
    return;
  }
  if (eExtension == "msh") {
    std::cerr << "Before call to WriteGridFile_msh\n";
    WriteGridFile_msh(GridFile, GrdArr);
    std::cerr << " After call to WriteGridFile_msh\n";
    return;
  }
  std::cerr << "eExtension = " << eExtension << "\n";
  std::cerr << "Programmed extensions are (gr3 , dat , grd , nc , msh)\n";
  throw TerminalException{1};
}



void WriteUnstructuredGridTot(std::string const& GridFile, std::string const& BoundFile, GridArray const& GrdArr)
{
  std::string eExtension=FILE_GetExtension(GridFile);
  GridArray GrdArrCopy=GrdArr;
  int nbNode=GrdArr.IOBP.size();
  MyMatrix<double> IOBPmat(nbNode,1);
  for (int ip=0; ip<nbNode; ip++) {
    double eVal=double(GrdArr.IOBP(ip));
    IOBPmat(ip,0)=eVal;
  }
  GrdArrCopy.GrdArrRho.DEP=IOBPmat;
  if (eExtension == "gr3") {
    WriteUnstructuredGrid_GR3(GridFile, GrdArr);
    WriteUnstructuredGrid_GR3(BoundFile, GrdArrCopy);
    return;
  }
  if (eExtension == "dat") {
    WriteUnstructuredGrid_DAT(GridFile, GrdArr);
    WriteUnstructuredGrid_DAT(BoundFile, GrdArrCopy);
    return;
  }
  if (eExtension == "nc") {
    WriteUnstructuredGrid_NC(GridFile, GrdArr);
    return;
  }
  if (eExtension == "msh") {
    WriteGridFile_msh(GridFile, GrdArr);
    return;
  }
  std::cerr << "eExtension = " << eExtension << "\n";
  std::cerr << "Programmed extensions are .gr3 , .dat and .nc only\n";
  throw TerminalException{1};
}


GridArray ReadUnstructuredGrid(std::string const& GridFile, std::string const& BoundFile)
{
  std::string eExtension=FILE_GetExtension(GridFile);
  if (!IsExistingFile(GridFile)) {
    std::cerr << "The file GridFile = " << GridFile << " is not existent\n";
    throw TerminalException{1};
  }
  //    std::cerr << "eExtension=" << eExtension << "\n";
  if (eExtension == "gr3" || eExtension == "ll") {
    GridArray GrdArr=WWM_ReadGridFile_gr3(GridFile);
    if (BoundFile != "unset") {
      std::string fExtension=FILE_GetExtension(BoundFile);
      if (fExtension != "gr3" && fExtension != "ll") {
	std::cerr << "Error in ReasUnstructuredGrid\n";
	std::cerr << "GridFile = " << GridFile << " so we are in gr3 case\n";
	std::cerr << "But BoundFile = " << BoundFile << "\n";
	std::cerr << "which is not ending by gr3 or ll\n";
	throw TerminalException{1};
      }
      std::cerr << "BoundFile=" << BoundFile << "\n";
      MyVector<int> eVect=WWM_ReadBoundFile_gr3(BoundFile);
      if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
	std::cerr << "not same number of vertices between grid file and boundary file\n";
	std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
	std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
	throw TerminalException{1};
      }
      GrdArr.IOBP=eVect;
    }
    return GrdArr;
  }
  if (eExtension == "dat") {
    std::cerr << "Before WWM_ReadGridFile_DAT\n";
    GridArray GrdArr=WWM_ReadGridFile_DAT(GridFile);
    std::cerr << "After  WWM_ReadGridFile_DAT\n";
    if (BoundFile != "unset") {
      std::string fExtension=FILE_GetExtension(BoundFile);
      if (fExtension != "dat") {
	std::cerr << "Error in ReasUnstructuredGrid\n";
	std::cerr << "GridFile = " << GridFile << " so we are in dat/xfn case\n";
	std::cerr << "But BoundFile = " << BoundFile << "\n";
	std::cerr << "which is not ending by dat\n";
	throw TerminalException{1};
      }
      MyVector<int> eVect=WWM_ReadBoundFile_DAT(BoundFile);
      if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
	std::cerr << "not same number of vertices between grid file and boundary file\n";
	std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
	std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
	throw TerminalException{1};
      }
      GrdArr.IOBP=eVect;
    }
    return GrdArr;
  }
  if (eExtension == "grd") {
    GridArray GrdArr=WWM_ReadGridFile_Ricchiuto_grd(GridFile);
    return GrdArr;
  }
  if (eExtension == "obj") {
    std::cerr << "Before WWM_ReadGridFile_obj\n";
    GridArray GrdArr=WWM_ReadGridFile_obj(GridFile);
    std::cerr << "After  WWM_ReadGridFile_obj\n";
    return GrdArr;
  }
  if (eExtension == "nc")
    return WWM_ReadGridFile_netcdf(GridFile);
  if (eExtension == "msh")
    return WWM_ReadGridFile_msh(GridFile);
  std::cerr << "Error in reading grid for WWM\n";
  std::cerr << "We did not find the right kind\n";
  throw TerminalException{1};
}




GridArray PRE_RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
  std::string PreModelName=eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  CHECK_Model_Allowedness(eModelName);
  std::string GridFile=GET_GRID_FILE(eTriple);
  std::cerr << "eModelName=" << eModelName << "\n";
  if (eModelName == "RECTANGULAR") {
    QuadArray eQuad;
    eQuad.MaxLat=eTriple.RecGridSymb.MaxLat;
    eQuad.MinLat=eTriple.RecGridSymb.MinLat;
    eQuad.MaxLon=eTriple.RecGridSymb.MaxLon;
    eQuad.MinLon=eTriple.RecGridSymb.MinLon;
    std::cerr << "Lat (min/max)=" << eQuad.MinLat << " / " << eQuad.MaxLat << "\n";
    std::cerr << "Lon (min/max)=" << eQuad.MinLon << " / " << eQuad.MaxLon << "\n";
    double deltaKM=eTriple.RecGridSymb.deltaKM;
    std::cerr << "deltaKM=" << deltaKM << "\n";
    double distLON=GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat, eQuad.MaxLon, eQuad.MinLat);
    double distLAT=GeodesicDistanceKM(eQuad.MinLon, eQuad.MinLat, eQuad.MinLon, eQuad.MaxLat);
    double nbLON=int(distLON / deltaKM);
    double nbLAT=int(distLAT / deltaKM);
    return RECTANGULAR_GRID_ARRAY(eQuad, nbLON, nbLAT);
  }
  if (eModelName == "COSMO")
    return NC_ReadCosmoWamStructGridFile(GridFile, "atm");
  if (eModelName == "WAM")
    return NC_ReadWamGridFile(GridFile);
  if (eModelName == "ROMS" || eModelName == "ROMS_IVICA")
    return NC_ReadRomsGridFile(GridFile);
  if (eModelName == "WRF")
    return NC_ReadWrfGridFile(GridFile);
  if (eModelName == "NEMO")
    return NC_ReadNemoGridFile(GridFile);
  if (eModelName == "WWM" || eModelName == "WWM_DAILY") {
    std::string BoundFile=eTriple.BoundFile;
    std::cerr << "GridFile=" << GridFile << "\n";
    std::cerr << "BoundFile=" << BoundFile << "\n";
    GridArray GrdArr=ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName="WWM";
    return GrdArr;
  }
  if (eModelName == "UNRUNOFF") {
    std::string BoundFile=eTriple.BoundFile;
    GridArray GrdArr=ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName="UNRUNOFF";
    return GrdArr;
  }
  if (eModelName == "WW3")
    return NC_ReadWW3_GridFile(GridFile);
  if (eModelName == "SCHISM_SFLUX")
    return NC_ReadSCHISM_sflux_grid(GridFile);
  if (eModelName == "SCHISM_NETCDF_OUT") {
    std::string BoundFile=eTriple.BoundFile;
    GridArray GrdArr=ReadUnstructuredGrid(GridFile, BoundFile);
    GrdArr.ModelName="SCHISM_NETCDF_OUT";
    return GrdArr;
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ALADIN" || eModelName == "GRIB_IFS") {
    std::string HisPrefix=eTriple.HisPrefix;
    std::vector<std::string> ListFile=FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      throw TerminalException{1};
    }
    std::string eFileName=ListFile[0];
    GridArray GrdArr=GRIB_ReadGridArray(eFileName, eModelName);
    std::cerr << "We have GrdArr in that case\n";
    return GrdArr;
  }
  if (eModelName == "GRIB_WAM_FORT30") {
    std::string eFileName=eTriple.HisPrefix;
    if (!IsExistingFile(eFileName)) {
      std::cerr << "The file eFileName = " << eFileName << " is missing\n";
      std::cerr << "This is set by HisPRefix and serves for the data storage\n";
      throw TerminalException{1};
    }
    return GRIB_ReadGridArray(eFileName, eModelName);
  }
  std::cerr << "Error in PRE_RETRIEVE_GRID_ARRAY\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  throw TerminalException{1};
}


void PrintGridArray(std::ostream & os, GridArray const& GrdArr)
{
  os << "IsFE=" << GrdArr.IsFE << "\n";
  QuadArray eArr=GetQuadArray(GrdArr);
  os << "Lon (min/max)=" << eArr.MinLon << " / " << eArr.MaxLon << "\n";
  os << "Lat (min/max)=" << eArr.MinLat << " / " << eArr.MaxLat << "\n";
}





GridArray RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
  //  std::cerr << "Before PRE_RETRIEVE_GRID_ARRAY\n";
  GridArray GrdArr=PRE_RETRIEVE_GRID_ARRAY(eTriple);
  std::string strSphericity=eTriple.RecGridSymb.Sphericity;
  if (strSphericity != "unset") {
    if (strSphericity != "Spherical" && strSphericity != "Cartesian") {
      std::cerr << "Error, we need the Sphericity option to be set to either\n";
      std::cerr << "unset: then the plotting software is guessing the right value or reading directly from netcdf file\n";
      std::cerr << "Spherical: if the grid is spherical (needed only if the guess is wrong)\n";
      std::cerr << "Cartesian: if the grid is cartesian (needed only if the guess is wrong)\n";
      throw TerminalException{1};
    }
    if (strSphericity == "Spherical")
      GrdArr.IsSpherical=true;
    if (strSphericity == "Cartesian")
      GrdArr.IsSpherical=false;
  }
  std::cerr << "IsSpherical=" << GrdArr.IsSpherical << "\n";
  //  std::cerr << "After PRE_RETRIEVE_GRID_ARRAY\n";
  if (GrdArr.IsFE == 0)
    return GrdArr;
  if (eTriple.RecGridSymb.CutWorldMap) {
    CutWorldMap(GrdArr);
    std::cerr << "After CutWorldMap\n";
  }
  if (eTriple.RecGridSymb.HigherLatitudeCut) {
    double MinLatCut=eTriple.RecGridSymb.MinLatCut;
    double MaxLatCut=eTriple.RecGridSymb.MaxLatCut;
    CUT_HigherLatitude(GrdArr, MinLatCut, MaxLatCut);
    std::cerr << "After CUT_HigherLatitude\n";
  }
  CHECK_UnstructuredGrid(GrdArr);
  CHECK_CombinatorialGrid(GrdArr);
  CHECK_COORDINATE_ORIENTATION(GrdArr);
  std::cerr << "Before returning GrdArr in RETRIEVE_GRID_ARRAY\n";
  return GrdArr;
}




ArrayHistory NC_ReadArrayHistory_NEMO(std::string const& HisPrefix)
{
  ArrayHistory eArr;
  eArr.KindArchive="NETCDF";
  eArr.HisPrefix=HisPrefix;
  return eArr;
}



ArrayHistory NC_ReadArrayHistory(TripleModelDesc const& eTriple)
{
  std::string StringTime="ocean_time";
  std::string PreModelName=eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  std::string HisPrefix=eTriple.HisPrefix;
  std::cerr << "Debug: NC_ReadArrayHistory\n";
  std::cerr << "eModelName=" << eModelName << "\n";
  // special models first
  if (eModelName == "WW3") {
    std::string HisFile=GET_GRID_FILE(eTriple);
    return WW3_ReadArrayHistory(HisFile, HisPrefix);
  }  
  if (eModelName == "ROMS_IVICA" || eModelName == "WWM_DAILY")
    return Sequential_ReadArrayHistory(HisPrefix);
  if (eModelName == "SCHISM_SFLUX")
    return NC_ReadArrayHistory_Kernel(HisPrefix, "time", 3);
  if (eModelName == "NEMO")
    return NC_ReadArrayHistory_NEMO(HisPrefix);
  // generic cases of well behaved models
  return NC_ReadArrayHistory_Kernel(HisPrefix, StringTime, 4);
}



/* 
   When we have several runs (typical in operational runs)
   then we need to select the best one.
 */
double GetOptimalTimeShiftLength(std::string const& eModelName)
{
  std::vector<std::string> LStr=STRING_Split(eModelName, ":");
  for (int i=1; i<int(LStr.size()); i++) {
    std::vector<std::string> LStrB=STRING_Split(LStr[i], "_");
    if (LStrB[0] == "optimaltime") {
      double eValHour;
      std::istringstream(LStrB[1]) >> eValHour;
      double eValDay = eValHour / double(24);
      return eValDay;
    }
  }
  return 0;
}


bool RetrieveAllStates(std::string const& eModelName)
{
  std::vector<std::string> LStr=STRING_Split(eModelName, ":");
  for (int i=1; i<int(LStr.size()); i++)
    if (LStr[i] == "retrieveallstates")
      return true;
  return false;
}





ArrayHistory GRIB_ReadArrayHistory_Kernel(std::vector<std::string> const& ListFile, std::string const& eModelName)
{
  int nbFile=ListFile.size();
  //
  // Determining parameters of the search
  //
  GRIB_CheckAllowedKeywords(eModelName);
  double OptimalShiftDay = GetOptimalTimeShiftLength(eModelName);
  std::cerr << "OptimalShiftDay = " << OptimalShiftDay << "\n";
  double MaxErrorTime=0.01;
  bool RetAllStates = RetrieveAllStates(eModelName);
  //
  // Determining the list of messages
  //
  std::vector<GRIB_MessageInfo> ListAllMessage;
  std::set<std::string> TotalListShortName;
  for (int iFile=0; iFile<nbFile; iFile++) {
    std::string eFile=ListFile[iFile];
    std::cerr << "iFile=" << iFile << " / " << nbFile << " eFile=" << eFile << "\n";
    //    std::cerr << "Before getting ListMessage eFile=" << eFile << "\n";
    std::vector<GRIB_MessageInfo> ListMessage=GRIB_GetAllMessagesFromFile(eFile, eModelName);
    //    std::cerr << "After getting ListMessage\n";
    int nbMessage=ListMessage.size();
    //    std::cerr << "eFile = " << eFile << " nbMessage=" << nbMessage << "\n";
    //PrintVectorGRIBmessageInfo(std::cerr, ListMessage);
    //    for (auto & eMesg : ListMessage)
    //      std::cerr << "  shortName=" << eMesg.shortName << "\n";
    //    std::cerr << "nbMessage = " << nbMessage << "\n";
    if (nbMessage == 0) {
      std::cerr << "Remark: eFile = " << eFile << " has zero messages\n";
    }
    for (auto & eMesg : ListMessage) {
      ListAllMessage.push_back(eMesg);
      TotalListShortName.insert(eMesg.shortName);
    }
  }
  int TotalNbMessage=ListAllMessage.size();
  if (TotalNbMessage == 0) {
    std::cerr << "TotalNbMessage=" << TotalNbMessage << "\n";
    std::cerr << "|ListFile|=" << ListFile.size() << "\n";
    std::cerr << "We have zero messages. No work can be done\n";
    throw TerminalException{1};
  }
  //  std::cerr << "1: |ListAllMessage|=" << ListAllMessage.size() << "\n";
  //
  // Now reordering the messages
  //
  sort(ListAllMessage.begin(), ListAllMessage.end(), 
       [&](GRIB_MessageInfo const& a, GRIB_MessageInfo const& b) -> bool {
	 if (a.time < b.time)
	   return true;
	 return false;
       });

  //  std::cerr << "2: |ListAllMessage|=" << ListAllMessage.size() << "\n";
  //
  // Determining the list of times
  //
  double TimePrev=ListAllMessage[0].time;
  std::vector<double> ListTime{TimePrev};
  std::vector<int> ListITime(TotalNbMessage);
  std::vector<int> ListIFile(TotalNbMessage,-1);
  std::vector<int> ListIRec(TotalNbMessage,-1);
  int posTime=0;
  for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
    GRIB_MessageInfo eMesg=ListAllMessage[iMesg];
    //    std::cerr << "shortName=" << eMesg.shortName << "\n";
    double eTime=eMesg.time;
    double TimeDiff = fabs(eTime - TimePrev);
    //    std::cerr << "iMesg=" << iMesg << " posTime=" << posTime << "\n";
    if (TimeDiff > MaxErrorTime) {
      ListTime.push_back(eTime);
      TimePrev=eTime;
      posTime++;
      //      std::cerr << "  TimeDiff=" << TimeDiff << " TimePrev=" << TimePrev << "\n";
    }
    ListITime[iMesg]=posTime;
  }
  //
  // Determination of List of starttime from the list of messages.
  //
  std::vector<double> ListStartTime;
  std::vector<int> ListIStartTime(TotalNbMessage);
  double tolDay=double(1)/double(10000);
  auto GetIStartTime=[&](double const& timeStart) -> int {
    int nbTimeStart=ListStartTime.size();
    for (int iTimeStart=0; iTimeStart<nbTimeStart; iTimeStart++)
      if (fabs(timeStart - ListStartTime[iTimeStart]) < tolDay)
	return iTimeStart;
    ListStartTime.push_back(timeStart);
    return nbTimeStart;
  };
  for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
    double eStartTime=ListAllMessage[iMesg].timeStart;
    int iTimeStart=GetIStartTime(eStartTime);
    ListIStartTime[iMesg]=iTimeStart;
  }
  bool PrintEssentialInfo=true;
  if (PrintEssentialInfo) {
    for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
      GRIB_MessageInfo eMesg=ListAllMessage[iMesg];
      double time=eMesg.time;
      double timestart=eMesg.timeStart;
      double deltaTimeTimeStart=timestart - time;
      std::cerr << "iMesg=" << iMesg << " shortName=" << eMesg.shortName << " time=" << eMesg.time << " deltaTT=" << deltaTimeTimeStart << " file=" << eMesg.FileName << " p=" << ListITime[iMesg] << "\n";
    }
  }
  int nbTime=ListTime.size();
  //  std::cerr << "nbTime=" << nbTime << "\n";
  //
  // Determining the list of messages according to time
  //
  std::vector<std::vector<GRIB_MessageInfo>> ListListMessages(nbTime);
  for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
    GRIB_MessageInfo eMesg=ListAllMessage[iMesg];
    //    std::cerr << "iMesg=" << iMesg << " eMesg.shortName=" << eMesg.shortName << "\n";
    int iTime=ListITime[iMesg];
    //    std::cerr << "iTime=" << iTime << "\n";
    if (iTime < 0 || iTime >= nbTime) {
      std::cerr << "iTime=" << iTime << " but nbTime=" << nbTime << "\n";
      throw TerminalException{1};
    }
    ListListMessages[iTime].push_back(eMesg);
  }
  //
  // Now cleaning the entries by prefering the entries with timeStart + OptimalTimeDay as small as possible.
  //
  std::map<std::string, std::vector<std::pair<double, std::vector<GRIB_MessageInfo>>>> FullOrganizedInfo;
  if (RetAllStates) {
    for (auto & eShortName : TotalListShortName) {
      FullOrganizedInfo[eShortName] = {};
    }
  }
  std::set<std::string> SetRawNames;
  for (int iTime=0; iTime<nbTime; iTime++) {
    double eTime=ListTime[iTime];
    std::string strPres=DATE_ConvertMjd2mystringPres(eTime);
    std::cerr << "iTime=" << iTime << " / " << nbTime << " eTime=" << eTime << " date=" << strPres << "\n";
    std::vector<GRIB_MessageInfo> ListMessages=ListListMessages[iTime];
    std::cerr << "1: |ListMessages|=" << ListMessages.size() << "\n";
    std::set<std::string> ListShortName;
    for (auto & eMesg : ListMessages) {
      ListShortName.insert(eMesg.shortName);
      SetRawNames.insert(eMesg.shortName);
    }
    //    std::cerr << "ListMessages:\n";
    //    for (auto & eMesg : ListMessages)
    //      std::cerr << "  shortName=" << eMesg.shortName << " file=" << eMesg.FileName << "\n";
    //    std::cerr << "2: |ListShortName|=" << ListShortName.size() << "\n";
    //
    std::vector<GRIB_MessageInfo> NewListMessages;
    std::vector<GRIB_MessageInfo> ListMatchingMessages;
    for (auto & eShortName : ListShortName) {
      std::cerr << "1: eShortName=" << eShortName << "\n";
      GRIB_MessageInfo NewMesg;
      double minPenaltyFct=10^(30);
      int nbMatch=0;
      for (auto & eMesg : ListMessages) {
	double ePenaltyFct = fabs(eMesg.timeStart + OptimalShiftDay - eMesg.time);
	if (eMesg.shortName == eShortName) {
	  nbMatch++;
	  if (RetAllStates) {
	    ListMatchingMessages.push_back(eMesg);
	  }
	  std::cerr << " " << eMesg.FileName << "\n";
	  if (ePenaltyFct < minPenaltyFct) {
	    NewMesg = eMesg;
	    minPenaltyFct = ePenaltyFct;
	  }
	}
      }
      std::cerr << "nbMatch=" << nbMatch << " Selected=" << NewMesg.FileName << "\n";
      NewListMessages.push_back(NewMesg);
      if (RetAllStates) {
	std::pair<double, std::vector<GRIB_MessageInfo>> ePair{eTime, ListMatchingMessages};
	FullOrganizedInfo[eShortName].push_back(ePair);
      }
      //      std::cerr << "2: eShortName=" << eShortName << "\n";
    }
    //    std::cerr << "NewListMessages:\n";
    //    for (auto & eMesg : NewListMessages)
    //      std::cerr << "  shortName=" << eMesg.shortName << " file=" << eMesg.FileName << "\n";
    ListListMessages[iTime] = NewListMessages;
  }
  std::vector<std::string> RawVarNames;
  std::map<std::string, std::vector<int>> MatchingByVariable;
  for (auto& eName : SetRawNames) {
    RawVarNames.push_back(eName);
    MatchingByVariable[eName]={};
  }
  for (int iTime=0; iTime<nbTime; iTime++) {
    for (auto & eMesg : ListListMessages[iTime]) {
      std::string eName=eMesg.shortName;
      MatchingByVariable[eName].push_back(iTime);
    }
  }
  double FirstTime=ListTime[0];
  double LastTime=ListTime[nbTime-1];
  ArrayHistory eArr;
  eArr.nbFile = nbFile;
  eArr.nbTime = nbTime;
  eArr.ListListMessages = ListListMessages;
  eArr.ListAllMessage = ListAllMessage;
  eArr.ListStartTime = ListStartTime;
  eArr.ListIStartTime = ListIStartTime;
  eArr.RawVarNames = RawVarNames;
  eArr.MatchingByVariable = MatchingByVariable;
  eArr.FullOrganizedInfo = FullOrganizedInfo;
  eArr.ListITime = ListITime;
  eArr.ListIFile = ListIFile;
  eArr.ListIRec = ListIRec;
  eArr.ListTime = ListTime;
  eArr.FirstTime = FirstTime;
  eArr.LastTime = LastTime;
  eArr.KindArchive = "GRIB";
  eArr.TimeSteppingInfo = "classic";
  eArr.AppendVarName = false;
  //  std::cerr << "Ending of GRIB_ReadArrayHistory\n";
  return eArr;
}







ArrayHistory GRIB_ReadArrayHistory(std::string const& HisPrefix, std::string const& eModelName)
{
  //
  // Determining the list of files.
  //
  std::cerr << "Begin of GRIB_ReadArrayHistory\n";
  std::vector<std::string> ListFile;
  if (IsExistingFile(HisPrefix) && FILE_IsRegularFile(HisPrefix)) {
    ListFile = {HisPrefix};
  }
  else {
    ListFile = FILE_DirectoryFilesSpecificExtension(HisPrefix, "grb");
  }
  return GRIB_ReadArrayHistory_Kernel(ListFile, eModelName);
}




ArrayHistory ReadArrayHistory(TripleModelDesc const& eTriple)
{
  ArrayHistory eArr;
  std::string HisPrefix=eTriple.HisPrefix;
  std::string PreModelName=eTriple.ModelName;
  std::string eModelName = GetKernelModelName(PreModelName);
  CHECK_Model_Allowedness(eModelName);
  std::vector<std::string> ListModelGrib{"GRIB_DWD", "GRIB_GFS", "GRIB_COSMO", "GRIB_ECMWF", "GRIB_ALADIN", "GRIB_IFS", "GRIB_WAM_FORT30"};
  if (PositionVect(ListModelGrib, eModelName) != -1) {
    return GRIB_ReadArrayHistory(HisPrefix, PreModelName);
  }
  std::vector<std::string> ListModelNetcdf{"COSMO", "WAM", "ROMS", "ROMS_IVICA", "WWM", "WWM_DAILY", "WW3", "SCHISM_SFLUX", "SCHISM_NETCDF_OUT", "RECTANGULAR", "WRF", "UNRUNOFF", "IVICA_UVP", "NEMO"};
  if (PositionVect(ListModelNetcdf, eModelName) != -1) {
    return NC_ReadArrayHistory(eTriple);
  }
  std::cerr << "Error in ReadArrayHistory. Could not find a matching method for creating the array history\n";
  throw TerminalException{1};
}


VerticalInfo GetVerticalInfo(int const& N)
{
  MyVector<double> Hz(N);
  MyVector<double> z_w(N+1);
  MyVector<double> z_r(N);
  return {Hz, z_w, z_r};
}




void ComputeHz(ARVDtyp const& ARVD, double const& hwater, double const& eZeta, VerticalInfo & eVert)
{
  int Vtrans=ARVD.Vtransform;
  int N=ARVD.N;
  if (Vtrans == 1) {
    double hinv=1/hwater;
    eVert.z_w(0)=-hwater;
    for (int k=1; k<=N; k++) {
      double cff_r=ARVD.hc * (ARVD.sc_r(k-1) - ARVD.Cs_r(k-1));
      double cff_w=ARVD.hc * (ARVD.sc_w(k) - ARVD.Cs_w(k));
      double cff1_r=ARVD.Cs_r(k-1);
      double cff1_w=ARVD.Cs_w(k);
      double z_w0=cff_w + cff1_w*hwater;
      eVert.z_w(k)=z_w0 + eZeta*(1+z_w0*hinv);
      double z_r0=cff_r + cff1_r*hwater;
      eVert.z_r(k-1)=z_r0 + eZeta*(1+z_r0*hinv);
      eVert.Hz(k-1)=eVert.z_w(k) - eVert.z_w(k-1);
    }
    return;
  }
  if (Vtrans == 2) {
    eVert.z_w(0)=-hwater;
    double hinv=1/(ARVD.hc + hwater);
    for (int k=1; k<=N; k++) {
      double cff_r=ARVD.hc * ARVD.sc_r(k-1);
      double cff_w=ARVD.hc * ARVD.sc_w(k);
      double cff1_r=ARVD.Cs_r(k-1);
      double cff1_w=ARVD.Cs_w(k);
      double cff2_r=(cff_r+cff1_r*hwater)*hinv;
      double cff2_w=(cff_w+cff1_w*hwater)*hinv;
      eVert.z_w(k)   = eZeta + (eZeta + hwater)*cff2_w;
      eVert.z_r(k-1) = eZeta + (eZeta + hwater)*cff2_r;
      eVert.Hz(k-1)  = eVert.z_w(k) - eVert.z_w(k-1);
    }
    return;
  }
  std::cerr << "Failed to find matching entry of Vtrans Vtrans=" << Vtrans << "\n";
  throw TerminalException{1};
}

Eigen::Tensor<double,3> ROMS_ComputeVerticalGlobalCoordinate(GridArray const& GrdArr, MyMatrix<double> const& zeta)
{
  int N=GrdArr.ARVD.N;
  int eta_rho=zeta.rows();
  int xi_rho=zeta.cols();
  VerticalInfo eVert=GetVerticalInfo(N);
  Eigen::Tensor<double,3> Zmat(N, eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      for (int k=0; k<N; k++)
	eVert.z_r(k)=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1)
	ComputeHz(GrdArr.ARVD, GrdArr.GrdArrRho.DEP(i,j), zeta(i,j), eVert);
      for (int k=0; k<N; k++)
	Zmat(k, i, j)=eVert.z_r(k);
    }
  return Zmat;
}



struct PairMSKfield {
  MyMatrix<int> MSK;
  MyMatrix<double> field;
};

PairMSKfield VerticalInterpolation_P1_W(ARVDtyp const& ARVD, MyMatrix<double> const& h, MyMatrix<double> const& zeta, MyMatrix<int> const& MSK, double const& dep, Eigen::Tensor<double,3> const& VertField_W)
{
  int eta=h.rows();
  int xi=h.cols();
  MyMatrix<int> MSKret(eta,xi);
  MyMatrix<double> FieldRet(eta,xi);
  int N=ARVD.N;
  VerticalInfo eVert=GetVerticalInfo(N);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      int eMSK=MSK(i,j);
      double eField=0;
      if (eMSK == 1) {
	if (dep < -h(i,j)) {
	  eMSK=0;
	}
	else {
	  ComputeHz(ARVD, h(i,j), zeta(i,j), eVert);
	  for (int iVert=0; iVert<N; iVert++) {
	    double dep1=eVert.z_w(iVert);
	    double dep2=eVert.z_w(iVert+1);
	    if (dep1 >= dep && dep <= dep2) {
	      double alpha1=(dep2 - dep)/(dep2 - dep1);
	      double alpha2=(dep - dep1)/(dep2 - dep1);
	      eField=alpha1*VertField_W(iVert,i,j) + alpha2*VertField_W(iVert+1,i,j);
	    }
	  }
	}
      }
      MSKret(i,j)=eMSK;
      FieldRet(i,j)=eField;
    }
  return {MSKret, FieldRet};
}





MyMatrix<double> VerticalInterpolation_P2_W(ARVDtyp const& ARVD, MyMatrix<double> const& h, MyMatrix<double> const& zeta, MyMatrix<int> const& MSK, double const& dep, Eigen::Tensor<double,3> const& VertField_W)
{
  PairMSKfield ePair=VerticalInterpolation_P1_W(ARVD, h, zeta, MSK, dep, VertField_W);
  int eta=h.rows();
  int xi=h.cols();
  MyMatrix<double> FieldRet(eta,xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal;
      if (ePair.MSK(i,j) == 0) {
	eVal=0;
      }
      else {
	eVal=ePair.field(i,j);
      }
      FieldRet(i,j)=eVal;
    }
  return FieldRet;
}


MyVector<double> GetVertCoord_R(ARVDtyp const& ARVD, double const& eDep, double const& eZeta)
{
  if (ARVD.Zcoordinate)
    return ARVD.ListZ_r;
  if (ARVD.ModelName == "ROMS") {
    VerticalInfo eVert=GetVerticalInfo(ARVD.N);
    ComputeHz(ARVD, eDep, eZeta, eVert);
    return eVert.z_r;
  }
  std::cerr << "Did not find any matching entry\n";
  throw TerminalException{1};
}


MyVector<double> GetVertCoord_W(ARVDtyp const& ARVD, double const& eDep, double const& eZeta)
{
  if (ARVD.Zcoordinate)
    return ARVD.ListZ_w;
  if (ARVD.ModelName == "ROMS") {
    VerticalInfo eVert=GetVerticalInfo(ARVD.N);
    ComputeHz(ARVD, eDep, eZeta, eVert);
    return eVert.z_w;
  }
  std::cerr << "Did not find any matching entry\n";
  throw TerminalException{1};
}



Eigen::Tensor<double,3> VerticalInterpolationTensor_R(GridArray const& GrdArrOut, ARVDtyp const& ARVDin, Eigen::Tensor<double,3> const& TensIn)
{
  int eta_rho=GrdArrOut.GrdArrRho.LON.rows();
  int xi_rho =GrdArrOut.GrdArrRho.LON.cols();
  int NvertOut=GrdArrOut.ARVD.N;
  int NvertIn =ARVDin.N;
  Eigen::Tensor<double,3> Fvert(NvertOut, eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eDep=GrdArrOut.GrdArrRho.DEP(i,j);
      double eZeta=0;
      MyVector<double> Zr_out = GetVertCoord_R(GrdArrOut.ARVD, eDep, eZeta);
      MyVector<double> Zr_in  = GetVertCoord_R(ARVDin, eDep, eZeta);
      //      std::cerr << "|Zr_out|=" << Zr_out.size() << " |Zr_in|=" << Zr_in.size() << "\n";
      //      std::cerr << "NvertOut=" << NvertOut << " NvertIn=" << NvertIn << "\n";
      for (int k=0; k<NvertOut; k++) {
	//	std::cerr << "k=" << k << "\n";
	double depW=Zr_out(k);
	//	std::cerr << "depW=" << depW << "\n";
	double eValOut=0;
	if (depW < Zr_in(0)) {
	  eValOut=TensIn(0,i,j);
	}
	else {
	  if (depW > Zr_in(NvertIn-1)) {
	    eValOut=TensIn(NvertIn-1,i,j);
	  }
	  else {
	    for (int u=1; u<NvertIn; u++) {
	      double dep1=Zr_in(u-1);
	      double dep2=Zr_in(u);
	      double alpha1=(dep2 - depW)/(dep2 - dep1);
	      double alpha2=(depW - dep1)/(dep2 - dep1);
	      eValOut = alpha1*TensIn(u-1, i, j) + alpha2*TensIn(u, i, j);
	    }
	  }
	}
	Fvert(k,i,j)=eValOut;
      }
    }
  return Fvert;
}




PairMSKfield VerticalInterpolation_P1_R(ARVDtyp const& ARVD, MyMatrix<double> const& h, MyMatrix<double> const& zeta, MyMatrix<int> const& MSK, double const& dep, Eigen::Tensor<double,3> const& VertField_R, int const& Choice)
{
  int eta=h.rows();
  int xi=h.cols();
  MyMatrix<int> MSKret(eta,xi);
  MyMatrix<double> FieldRet(eta,xi);
  int N=ARVD.N;
  VerticalInfo eVert=GetVerticalInfo(N);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      int eMSK=MSK(i,j);
      double depW;
      if (Choice == 1)
	depW=dep;
      else
	depW=dep + zeta(i,j);
      double eField=0;
      if (eMSK == 1) {
	if (depW < -h(i,j)) {
	  eMSK=1;
	}
	else {
	  ComputeHz(ARVD, h(i,j), zeta(i,j), eVert);
	  bool WeMatch=false;
	  if (depW <= eVert.z_r(0)) {
	    eField=VertField_R(0,i,j);
	    WeMatch=true;
	  }
	  if (eVert.z_r(N-1) <= depW) {
	    eField=VertField_R(N-1,i,j);
	    WeMatch=true;
	  }
	  else {
	    for (int iVert=0; iVert<N-1; iVert++) {
	      double dep1=eVert.z_r(iVert);
	      double dep2=eVert.z_r(iVert+1);
	      //	      std::cerr << "iVert=" << iVert << " dep1=" << dep1 << " dep2=" << dep2 << "\n";
	      if (dep1 <= depW && depW <= dep2) {
		double alpha1=(dep2 - depW)/(dep2 - dep1);
		double alpha2=(depW - dep1)/(dep2 - dep1);
		eField=alpha1*VertField_R(iVert,i,j) + alpha2*VertField_R(iVert+1,i,j);
		WeMatch=true;
	      }
	    }
	  }
	  if (!WeMatch) {
	    std::cerr << "No assignation for the vertical interpolation\n";
	    std::cerr << "i=" << i << " j=" << j << "\n";
	    std::cerr << "zeta=" << zeta(i,j) << " h=" << h(i,j) << "\n";
	    std::cerr << "dep=" << dep << " N=" << N << "\n";
	    std::cerr << "depW=" << depW << "\n";
	    for (int iVert=0; iVert<N; iVert++)
	      std::cerr << "  iVert=" << iVert << " z_z=" << eVert.z_r(iVert) << "\n";
	    throw TerminalException{1};
	  }
	}
      }
      MSKret(i,j)=eMSK;
      FieldRet(i,j)=eField;
    }
  return {MSKret, FieldRet};
}



MyMatrix<double> ConvertBaroclinic_to_Barotropic(Eigen::Tensor<double,3> const& F3, MyMatrix<double> const& zeta, GridArray const& GrdArr)
{
  auto LDim=F3.dimensions();
  int N=LDim[0];
  int eta_rho=LDim[1];
  int xi_rho=LDim[2];
  int eta_rho_grid=GrdArr.GrdArrRho.LON.rows();
  int xi_rho_grid =GrdArr.GrdArrRho.LON.cols();
  if (N != GrdArr.ARVD.N) {
    std::cerr << "First dimension of F1 is N=" << N << "\n";
    std::cerr << "But GrdArr.ARVD.N=" << GrdArr.ARVD.N << "\n";
    throw TerminalException{1};
  }
  if (eta_rho != eta_rho_grid || xi_rho != xi_rho_grid) {
    std::cerr << "eta_rho      = " << eta_rho      << " xi_rho      = " << xi_rho << "\n";
    std::cerr << "eta_rho_grid = " << eta_rho_grid << " xi_rho_grid = " << "\n";
    std::cerr << "They should be equal\n";
    throw TerminalException{1};
  }
  MyMatrix<double> F(eta_rho, xi_rho);
  VerticalInfo eVert=GetVerticalInfo(N);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eVal=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	ComputeHz(GrdArr.ARVD, GrdArr.GrdArrRho.DEP(i,j), zeta(i,j), eVert);
	double VertInt=0;
	for (int k=0; k<N; k++) {
	  double dep=eVert.z_w(k+1) - eVert.z_w(k);
	  VertInt += dep*F3(k,i,j);
	}
	double DeltaDep=zeta(i,j) + GrdArr.GrdArrRho.DEP(i,j);
	eVal=VertInt/DeltaDep;
      }
      F(i,j)=eVal;
    }
  return F;
}





MyMatrix<double> VerticalInterpolation_SCHISM_ZNL(Eigen::Tensor<double,3> const& znl, MyMatrix<double> const& zeta, double const& dep, Eigen::Tensor<double,3> const& VertField_R, int const& Choice)
{
  int eta=zeta.rows();
  int xi =zeta.cols();
  auto LDim=znl.dimensions();
  int N=LDim[0];
  MyMatrix<double> FieldRet(eta,xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double depSearch;
      if (Choice == 1)
	depSearch=dep;
      else
	depSearch=dep + zeta(i,j);
      double eField=0;
      if (depSearch <= znl(0,i,j)) {
	eField=VertField_R(0,i,j);
      }
      if (znl(N-1,i,j) <= depSearch) {
	eField=VertField_R(N-1,i,j);
      }
      else {
	for (int iVert=0; iVert<N-1; iVert++) {
	  double dep1=znl(iVert  ,i,j);
	  double dep2=znl(iVert+1,i,j);
	  //	      std::cerr << "iVert=" << iVert << " dep1=" << dep1 << " dep2=" << dep2 << "\n";
	  if (dep1 <= depSearch && depSearch <= dep2) {
	    double alpha1=(dep2 - depSearch)/(dep2 - dep1);
	    double alpha2=(depSearch - dep1)/(dep2 - dep1);
	    eField=alpha1*VertField_R(iVert,i,j) + alpha2*VertField_R(iVert+1,i,j);
	  }
	}
      }
      FieldRet(i,j)=eField;
    }
  return FieldRet;
}



MyMatrix<double> VerticalInterpolation_P2_R(ARVDtyp const& ARVD, MyMatrix<double> const& h, MyMatrix<double> const& zeta, MyMatrix<int> const& MSK, double const& dep, Eigen::Tensor<double,3> const& VertField_R, int const& Choice)
{
  PairMSKfield ePair=VerticalInterpolation_P1_R(ARVD, h, zeta, MSK, dep, VertField_R, Choice);
  int eta=h.rows();
  int xi=h.cols();
  MyMatrix<double> FieldRet(eta,xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal;
      if (ePair.MSK(i,j) == 0) {
	eVal=0;
      }
      else {
	eVal=ePair.field(i,j);
      }
      FieldRet(i,j)=eVal;
    }
  return FieldRet;  
}


double GetUnitInMeter(double const& eVal, std::string const& unit)
{
  if (unit == "m")
    return eVal;
  if (unit == "dm")
    return eVal*double(0.1);
  if (unit == "cm")
    return eVal*double(0.01);
  if (unit == "mm")
    return eVal*double(0.001);
  if (unit == "km")
    return eVal*double(1000);
  std::cerr << "We should never reach that stage in GetUnitInMeter\n";
  throw TerminalException{1};
}







#endif
