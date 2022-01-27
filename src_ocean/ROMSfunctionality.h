#ifndef ROMS_FUNCTIONALITY_INCLUDE
#define ROMS_FUNCTIONALITY_INCLUDE


#include "Basic_netcdf.h"
#include "OneDimInterpolation.h"
#include "Namelist.h"
#include "Model_grids.h"
#include "Statistics.h"


FullNamelist NAMELIST_ROMS_VERTICAL_STRATIFICATION_DIAGNOSTIC()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DESCRIPTION
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["GridFile"]="unset";
  ListIntValues1["N"] = 20;
  ListIntValues1["Vtransform"] = 20;
  ListIntValues1["Vstretching"] = 20;
  ListDoubleValues1["ThetaS"] = 7;
  ListDoubleValues1["ThetaB"] = 7;
  ListDoubleValues1["Tcline"] = 7;
  ListDoubleValues1["hc"] = 7;
  SingleBlock BlockDESC;
  BlockDESC.ListIntValues=ListIntValues1;
  BlockDESC.ListDoubleValues=ListDoubleValues1;
  BlockDESC.ListStringValues=ListStringValues1;
  ListBlock["DESCRIPTION"]=BlockDESC;
  //
  return {std::move(ListBlock), "undefined"};
}


double ComputeHydrostaticInconsistencyNumber(Eigen::Tensor<double,3> const& z_r, GridArray const& GrdArr)
{
  auto compute_rx1=[&](int const& iEta, int const& iXi, int const& jEta, int const& jXi, int k) -> double {
    double z_e1_kP = z_r(k  ,iEta,iXi);
    double z_e2_kP = z_r(k  ,jEta,jXi);
    double z_e1_kM = z_r(k-1,iEta,iXi);
    double z_e2_kM = z_r(k-1,jEta,jXi);
    double num = z_e1_kP - z_e2_kP + z_e1_kM - z_e2_kM;
    double den = z_e1_kP + z_e2_kP - z_e1_kM - z_e2_kM;
    double coef = fabs(num) / fabs(den);
    return coef;
  };
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho  = GrdArr.GrdArrRho.LON.cols();
  int N = GrdArr.ARVD.N;
  std::vector<std::vector<int>> LNeigh{{1,0},{0,1},{-1,0},{0,-1}};
  double rx1 = 0;
  for (int iEta=0; iEta<eta_rho; iEta++) {
    for (int iXi=0; iXi<xi_rho; iXi++) {
      if (GrdArr.GrdArrRho.MSK(iEta, iXi)) {
        for (auto & eNeigh : LNeigh) {
          int iEtaN = iEta + eNeigh[0];
          int iXiN = iXi + eNeigh[1];
          if (iEtaN >= 0 && iEtaN <eta_rho && iXiN >= 0 && iXiN < xi_rho) {
            if (GrdArr.GrdArrRho.MSK(iEtaN, iXiN)) {
              for (int k=1; k<N; k++) {
                double val = compute_rx1(iEta, iXi, iEtaN, iXiN, k);
                if (val > rx1)
                  rx1 = val;
              }
            }
          }
        }
      }
    }
  }
  return rx1;
}







void DiagnosticsVerticalStratificationDiagnostic(FullNamelist const& eFull)
{
  SingleBlock BlockDESC=eFull.ListBlock.at("DESCRIPTION");
  std::string GridFile = BlockDESC.ListStringValues.at("GridFile");
  int N = BlockDESC.ListIntValues.at("N");
  int Vtransform = BlockDESC.ListIntValues.at("Vtransform");
  int Vstretching = BlockDESC.ListIntValues.at("Vstretching");
  double ThetaS = BlockDESC.ListDoubleValues.at("ThetaS");
  double ThetaB = BlockDESC.ListDoubleValues.at("ThetaB");
  double Tcline = BlockDESC.ListDoubleValues.at("Tcline");
  double hc = BlockDESC.ListDoubleValues.at("hc");
  //
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  ARVDtyp ARVD = ROMSgetARrayVerticalDescription(N, Vtransform, Vstretching, Tcline, hc, ThetaS, ThetaB);
  GrdArr.ARVD = ARVD;
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho  = GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> zeta = ZeroMatrix<double>(eta_rho, xi_rho);
  Eigen::Tensor<double,3> z_r = ROMS_ComputeVerticalGlobalCoordinate(GrdArr, zeta);
  //
  double rx1 = ComputeHydrostaticInconsistencyNumber(z_r, GrdArr);
  std::cerr << "rx1=" << rx1 << "\n";
  int s_rho=z_r.dimension(0);
  MyMatrix<double> HighestLevel = DimensionExtraction(z_r, 0, s_rho-1);
  double minDep = std::numeric_limits<double>::max();
  double maxDep = std::numeric_limits<double>::min();
  double sumDep = 0;
  size_t nb = 0;
  for (int iEta=0; iEta<eta_rho; iEta++)
    for (int iXi=0; iXi<xi_rho; iXi++)
      if (GrdArr.GrdArrRho.MSK(iEta,iXi)) {
        double eDep = HighestLevel(iEta, iXi);
        minDep = std::min(minDep, eDep);
        maxDep = std::max(maxDep, eDep);
        sumDep += eDep;
        nb++;
      }
  double avgDep = sumDep / nb;
  std::cerr << "minDep=" << minDep << " maxDep=" << maxDep << " avgDep=" << avgDep << "\n";
}




MyMatrix<double> My_u2rho(MyMatrix<double> const& eVar_u, MyMatrix<uint8_t> const& MSK_u)
{
  int eta_u=eVar_u.rows();
  int xi_u=eVar_u.cols();
  int eta_rho = eta_u;
  int xi_rho = xi_u + 1;
  if (!IsEqualSizeMatrices(eVar_u, MSK_u)) {
    std::cerr << "|eVar_u|=" << eVar_u.rows() << " / " << eVar_u.cols() << "\n";
    std::cerr << "|MSK_u|=" << MSK_u.rows() << " / " << MSK_u.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the chosen history\n";
    throw TerminalException{1};
  }
  //
  MyMatrix<double> eVar_rho(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      double eSumVal=0;
      if (j<xi_u) {
	if (MSK_u(i,j) == 1) {
	  eSumMsk++;
	  eSumVal += eVar_u(i,j);
	}
      }
      if (j > 0) {
	if (MSK_u(i,j-1) == 1) {
	  eSumMsk++;
	  eSumVal += eVar_u(i,j-1);
	}
      }
      if (eSumMsk == 0) {
	eVar_rho(i,j)=0;
      } else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}


MyMatrix<double> My_v2rho(MyMatrix<double> const& eVar_v, MyMatrix<uint8_t> const& MSK_v)
{
  int eta_v=eVar_v.rows();
  int xi_v=eVar_v.cols();
  int xi_rho = xi_v;
  int eta_rho = eta_v + 1;
  if (!IsEqualSizeMatrices(eVar_v, MSK_v)) {
    std::cerr << "|eVar_v|=" << eVar_v.rows() << " / " << eVar_v.cols() << "\n";
    std::cerr << "|MSK_v|=" << MSK_v.rows() << " / " << MSK_v.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the chosen history\n";
    throw TerminalException{1};
  }
  MyMatrix<double> eVar_rho(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      double eSumVal=0;
      if (i < eta_v) {
	if (MSK_v(i,j) == 1) {
	  eSumMsk++;
	  eSumVal += eVar_v(i,j);
	}
      }
      if (i > 0) {
	if (MSK_v(i-1,j) == 1) {
	  eSumMsk++;
	  eSumVal += eVar_v(i-1,j);
	}
      }
      if (eSumMsk == 0) {
	eVar_rho(i,j)=0;
      } else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}



Eigen::Tensor<double,3> My_v2rho_3D(Eigen::Tensor<double,3> const& eVar_v, MyMatrix<uint8_t> const& MSK_v)
{
  auto LDim=eVar_v.dimensions();
  int s_vert=LDim[0];
  int eta_v=LDim[1];
  int xi_v=LDim[2];
  int xi_rho = xi_v;
  int eta_rho = eta_v + 1;
  if (eta_v != int(MSK_v.rows()) || xi_v != int(MSK_v.cols())) {
    std::cerr << "|eVar_v|=" << eta_v << " / " << xi_v << "\n";
    std::cerr << "|MSK_v|=" << MSK_v.rows() << " / " << MSK_v.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the chosen history\n";
    throw TerminalException{1};
  }
  std::vector<double> VertColumn(s_vert);
  Eigen::Tensor<double,3> eVar_rho(s_vert,eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      for (int k=0; k<s_vert; k++)
	VertColumn[k]=0;
      if (i < eta_v) {
	if (MSK_v(i,j) == 1) {
	  eSumMsk++;
	  for (int k=0; k<s_vert; k++)
	    VertColumn[k] += eVar_v(k,i,j);
	}
      }
      if (i > 0) {
	if (MSK_v(i-1,j) == 1) {
	  eSumMsk++;
	  for (int k=0; k<s_vert; k++)
	    VertColumn[k] += eVar_v(k,i-1,j);
	}
      }
      if (eSumMsk == 0) {
	for (int k=0; k<s_vert; k++)
	  eVar_rho(k,i,j)=0;
      } else {
	for (int k=0; k<s_vert; k++) {
	  double eVal=VertColumn[k]/double(eSumMsk);
	  eVar_rho(k,i,j)=eVal;
	}
      }
    }
  return eVar_rho;
}

MyMatrix<double> My_rho2u_2D(GridArray const& GrdArr, MyMatrix<double> const& M)
{
  int eta_rho=M.rows();
  int xi_rho =M.cols();
  int eta_u=eta_rho;
  int xi_u=xi_rho-1;
  MyMatrix<double> Mout(eta_u, xi_u);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u; j++) {
      int sumMSK=0;
      double sumVal=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	sumMSK++;
	sumVal += M(i,j);
      }
      if (GrdArr.GrdArrRho.MSK(i,j+1) == 1) {
	sumMSK++;
	sumVal += M(i,j+1);
      }
      if (sumMSK > 0)
	sumVal /= sumMSK;
      Mout(i,j)=sumVal;
    }
  return Mout;
}



MyMatrix<double> My_rho2v_2D(GridArray const& GrdArr, MyMatrix<double> const& M)
{
  int eta_rho=M.rows();
  int xi_rho =M.cols();
  int eta_v=eta_rho-1;
  int xi_v=xi_rho;
  MyMatrix<double> Mout(eta_v, xi_v);
  for (int i=0; i<eta_v; i++)
    for (int j=0; j<xi_v; j++) {
      int sumMSK=0;
      double sumVal=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	sumMSK++;
	sumVal += M(i,j);
      }
      if (GrdArr.GrdArrRho.MSK(i+1,j) == 1) {
	sumMSK++;
	sumVal += M(i+1,j);
      }
      if (sumMSK > 0)
	sumVal /= sumMSK;
      Mout(i,j)=sumVal;
    }
  return Mout;
}




Eigen::Tensor<double,3> My_rho2u_3D(GridArray const& GrdArr, Eigen::Tensor<double,3> const& M)
{
  auto LDim=M.dimensions();
  int s_vert=LDim[0];
  int eta_rho=LDim[1];
  int xi_rho=LDim[2];
  int eta_u=eta_rho;
  int xi_u=xi_rho-1;
  Eigen::Tensor<double,3> Mout(s_vert, eta_u, xi_u);
  MyVector<double> VertCol(s_vert);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u; j++) {
      int sumMSK=0;
      for (int k=0; k<s_vert; k++)
	VertCol(k)=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	sumMSK++;
	for (int k=0; k<s_vert; k++)
	  VertCol(k) += M(k,i,j);
      }
      if (GrdArr.GrdArrRho.MSK(i,j+1) == 1) {
	sumMSK++;
	for (int k=0; k<s_vert; k++)
	  VertCol(k) += M(k,i,j+1);
      }
      if (sumMSK > 0) {
	for (int k=0; k<s_vert; k++)
	  Mout(k,i,j)=VertCol(k) / sumMSK;
      } else {
	for (int k=0; k<s_vert; k++)
	  Mout(k,i,j) = 0;
      }
    }
  return Mout;
}



Eigen::Tensor<double,3> My_rho2v_3D(GridArray const& GrdArr, Eigen::Tensor<double,3> const& M)
{
  auto LDim=M.dimensions();
  int s_vert=LDim[0];
  int eta_rho=LDim[1];
  int xi_rho=LDim[2];
  int eta_v=eta_rho-1;
  int xi_v=xi_rho;
  Eigen::Tensor<double,3> Mout(s_vert, eta_v, xi_v);
  MyVector<double> VertCol(s_vert);
  for (int i=0; i<eta_v; i++)
    for (int j=0; j<xi_v; j++) {
      int sumMSK=0;
      for (int k=0; k<s_vert; k++)
	VertCol(k)=0;
      if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
	sumMSK++;
	for (int k=0; k<s_vert; k++)
	  VertCol(k) += M(k, i, j);
      }
      if (GrdArr.GrdArrRho.MSK(i+1,j) == 1) {
	sumMSK++;
	for (int k=0; k<s_vert; k++)
	  VertCol(k) += M(k, i+1, j);
      }
      if (sumMSK > 0) {
	for (int k=0; k<s_vert; k++)
	  Mout(k, i, j) = VertCol(k) / sumMSK;
      } else {
	for (int k=0; k<s_vert; k++)
	  Mout(k, i, j)=0;
      }
    }
  return Mout;
}


Eigen::Tensor<double,3> My_u2rho_3D(Eigen::Tensor<double,3> const& eVar_u, MyMatrix<uint8_t> const& MSK_u)
{
  auto LDim=eVar_u.dimensions();
  int s_vert=LDim[0];
  int eta_u=LDim[1];
  int xi_u=LDim[2];
  int eta_rho = eta_u;
  int xi_rho = xi_u + 1;
  if (eta_u != int(MSK_u.rows()) || xi_u != int(MSK_u.cols())) {
    std::cerr << "|eVar_u|=" << eta_u << " / " << xi_u << "\n";
    std::cerr << "|MSK_u|=" << MSK_u.rows() << " / " << MSK_u.cols() << "\n";
    std::cerr << "Most likely the grid array does not match the chosen history\n";
    throw TerminalException{1};
  }
  //
  std::vector<double> VertColumn(s_vert);
  Eigen::Tensor<double,3> eVar_rho(s_vert, eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      for (int k=0; k<s_vert; k++)
	VertColumn[k]=0;
      if (j<xi_u) {
	if (MSK_u(i,j) == 1) {
	  eSumMsk++;
	  for (int k=0; k<s_vert; k++)
	    VertColumn[k] += eVar_u(k,i,j);
	}
      }
      if (j > 0) {
	if (MSK_u(i,j-1) == 1) {
	  eSumMsk++;
	  for (int k=0; k<s_vert; k++)
	    VertColumn[k] += eVar_u(k,i,j-1);
	}
      }
      if (eSumMsk == 0) {
	for (int k=0; k<s_vert; k++)
	  eVar_rho(k,i,j)=0;
      } else {
	for (int k=0; k<s_vert; k++) {
	  double eVal=VertColumn[k]/double(eSumMsk);
	  eVar_rho(k,i,j)=eVal;
	}
      }
    }
  return eVar_rho;
}






void AngleRhoRot(MyMatrix<double> & U_rho, MyMatrix<double> & V_rho,
		 MyMatrix<double> const& ANG)
{
  if (U_rho.rows() != ANG.rows() || V_rho.rows() != ANG.rows() || U_rho.cols() != ANG.cols() || V_rho.cols() != ANG.cols()) {
    std::cerr << "Error in AngleRhoRot\n";
    std::cerr << "|U_rho|=" << U_rho.rows() << " / " << U_rho.cols() << "\n";
    std::cerr << "|V_rho|=" << V_rho.rows() << " / " << V_rho.cols() << "\n";
    std::cerr << "|ANG|=" << ANG.rows() << " / " << ANG.cols() << "\n";
    std::cerr << "Most likely the grid used does not match the history file used\n";
    throw TerminalException{1};
  }
  int eta_rho=U_rho.rows();
  int xi_rho=U_rho.cols();
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eU=U_rho(i,j);
      double eV=V_rho(i,j);
      double eAng=ANG(i,j);
      double NewU=eU*cos(eAng) - eV*sin(eAng);
      double NewV=eU*sin(eAng) + eV*cos(eAng);
      U_rho(i,j)=NewU;
      V_rho(i,j)=NewV;
    }
}



void AngleRhoRot_3D(Eigen::Tensor<double,3> & U_rho, Eigen::Tensor<double,3> & V_rho,
		    MyMatrix<double> const& ANG)
{
  auto LDim=U_rho.dimensions();
  int s_vert=LDim[0];
  int eta_rho=LDim[1];
  int xi_rho=LDim[2];
  auto LDimV=V_rho.dimensions();
  for (int i=0; i<3; i++)
    if (LDim[i] != LDimV[i]) {
      std::cerr << "Error in AngleRhoRot_3D\n";
      std::cerr << "Different dimension for U_rho and V_rho at posiiton i=" << i << "\n";
      throw TerminalException{1};
    }
  if (eta_rho != ANG.rows() || xi_rho != ANG.cols()) {
    std::cerr << "Error in AngleRhoRot_3D\n";
    std::cerr << "The dimensions of U_rho and ANG do not match\n";
    std::cerr << "|U_rho|=" << s_vert << " / " << eta_rho << " / " << xi_rho << "\n";
    std::cerr << "|ANG|=" << ANG.rows() << " / " << ANG.cols() << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eAng=ANG(i,j);
      double eCos=cos(eAng);
      double eSin=sin(eAng);
      for (int k=0; k<s_vert; k++) {
	double eU=U_rho(k, i,j);
	double eV=V_rho(k, i,j);
	double NewU=eU*eCos - eV*eSin;
	double NewV=eU*eSin + eV*eCos;
	U_rho(k, i,j)=NewU;
	V_rho(k, i,j)=NewV;
      }
    }
}



struct ROMSgridArray {
  ARVDtyp ARVD;
  MyMatrix<uint8_t> MSK_rho;
  MyMatrix<uint8_t> MSK_u;
  MyMatrix<uint8_t> MSK_v;
  MyMatrix<uint8_t> MSK_psi;
  MyMatrix<double> LON_u;
  MyMatrix<double> LAT_u;
  MyMatrix<double> LON_v;
  MyMatrix<double> LAT_v;
  MyMatrix<double> LON_rho;
  MyMatrix<double> LAT_rho;
  MyMatrix<double> LON_psi;
  MyMatrix<double> LAT_psi;
  //
  MyMatrix<double> pm;
  MyMatrix<double> pn;
  MyMatrix<double> h;
  MyMatrix<double> angle;
};


ROMSgridArray ReadFullROMSgridArray(std::string const& eFile)
{
  ROMSgridArray eA;
  eA.LON_rho=NC_Read2Dvariable(eFile, "lon_rho");
  eA.LAT_rho=NC_Read2Dvariable(eFile, "lat_rho");
  eA.LON_psi=NC_Read2Dvariable(eFile, "lon_psi");
  eA.LAT_psi=NC_Read2Dvariable(eFile, "lat_psi");
  eA.LON_u=NC_Read2Dvariable(eFile, "lon_u");
  eA.LAT_u=NC_Read2Dvariable(eFile, "lat_u");
  eA.LON_v=NC_Read2Dvariable(eFile, "lon_v");
  eA.LAT_v=NC_Read2Dvariable(eFile, "lat_v");
  //
  eA.h=NC_Read2Dvariable(eFile, "h");
  eA.pm=NC_Read2Dvariable(eFile, "pm");
  eA.pn=NC_Read2Dvariable(eFile, "pn");
  eA.angle=NC_Read2Dvariable(eFile, "angle");
  MyMatrix<double> MSK_rho_double=NC_Read2Dvariable(eFile, "mask_rho");
  eA.MSK_rho=UniversalMatrixConversion<uint8_t,double>(MSK_rho_double);
  //
  MyMatrix<double> MSK_u_double=NC_Read2Dvariable(eFile, "mask_u");
  eA.MSK_u=UniversalMatrixConversion<uint8_t,double>(MSK_u_double);
  //
  MyMatrix<double> MSK_v_double=NC_Read2Dvariable(eFile, "mask_v");
  eA.MSK_v=UniversalMatrixConversion<uint8_t,double>(MSK_v_double);
  //
  MyMatrix<double> MSK_psi_double=NC_Read2Dvariable(eFile, "mask_psi");
  eA.MSK_psi=UniversalMatrixConversion<uint8_t,double>(MSK_psi_double);
  //
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data_Vtrans=dataFile.getVar("Vtransform");
  netCDF::NcVar data_Vstret=dataFile.getVar("Vstretching");
  netCDF::NcVar data_theta_s=dataFile.getVar("theta_s");
  netCDF::NcVar data_theta_b=dataFile.getVar("theta_b");
  netCDF::NcVar data_Tcline=dataFile.getVar("Tcline");
  netCDF::NcVar data_hc=dataFile.getVar("hc");
  //
  int eValI;
  double eValD;
  //
  data_Vtrans.getVar(&eValI);
  eA.ARVD.Vtransform=eValI;
  data_Vstret.getVar(&eValI);
  eA.ARVD.Vstretching=eValI;
  data_theta_s.getVar(&eValD);
  eA.ARVD.theta_s=eValD;
  data_theta_b.getVar(&eValD);
  eA.ARVD.theta_b=eValD;
  data_Tcline.getVar(&eValD);
  eA.ARVD.Tcline=eValD;
  data_hc.getVar(&eValD);
  eA.ARVD.hc=eValD;
  //
  netCDF::NcVar data_Cs_r=dataFile.getVar("Cs_r");
  netCDF::NcVar data_Cs_w=dataFile.getVar("Cs_w");
  netCDF::NcVar data_s_r=dataFile.getVar("s_rho");
  netCDF::NcVar data_s_w=dataFile.getVar("s_w");
  netCDF::NcDim eDim=data_Cs_r.getDim(0);
  int dim_s_r=eDim.getSize();
  netCDF::NcDim fDim=data_Cs_w.getDim(0);
  int dim_s_w=fDim.getSize();
  std::vector<double> eVarR(dim_s_r), eVarW(dim_s_w);
  data_Cs_r.getVar(eVarR.data());
  for (int i=0; i<dim_s_r; i++)
    eA.ARVD.Cs_r(i)=eVarR[i];
  data_s_r.getVar(eVarR.data());
  for (int i=0; i<dim_s_r; i++)
    eA.ARVD.sc_r(i)=eVarR[i];
  data_Cs_w.getVar(eVarW.data());
  for (int i=0; i<dim_s_w; i++)
    eA.ARVD.Cs_w(i)=eVarW[i];
  data_s_w.getVar(eVarW.data());
  for (int i=0; i<dim_s_w; i++)
    eA.ARVD.sc_w(i)=eVarW[i];
  eA.ARVD.N=dim_s_r;
  //
  return eA;
}


FullNamelist Individual_Tracer_Variable_File()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DESCRIPTION
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["VariableName"]="unset";
  ListStringValues1["MethodSetting"]="unset"; // Possible values are "Constant", "Nearest", "VerticalProfile", "YearlyClimatology", etc.
  ListDoubleValues1["ConstantValue"] = -1.0;
  ListListDoubleValues1["ListVerticalLevels"] = {};
  ListListDoubleValues1["ListValues"] = {};
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
  return {std::move(ListBlock), "undefined"};
}

std::vector<std::string> GetListVariableBFM_full()
{
  std::vector<std::string> ListVar{"iOxyg", "iPO4_", "iNO3_", "iNH4_", "iO4n_", "iSiOH", "iN6r_", "iB1c_", "iB1n_", "iB1p_", "iP1c_", "iP1n_", "iP1p_", "iP1l_", "iP1s_", "iP2c_", "iP2n_", "iP2p_", "iP2l_", "iP3c_", "iP3n_", "iP3p_", "iP3l_", "iP4c_", "iP4n_", "iP4p_", "iP4l_", "iZ3c_", "iZ3n_", "iZ3p_", "iZ4c_", "iZ4n_", "iZ4p_", "iZ5c_", "iZ5n_", "iZ5p_", "iZ6c_", "iZ6n_", "iZ6p_", "iR1c_", "iR1n_", "iR1p_", "iR2c_", "iR3c_", "iR6c_", "iR6n_", "iR6p_", "iR6s_", "iO3c_", "iO3h_", "iEIR_", "iDIC_", "iChlo", "siP1_", "siP2_", "siP3_", "siP4_", "eiP1_", "eiP2_", "eiP3_", "eiP4_", "ruPTc", "ruZTc", "ixEPS"};
  return ListVar;
}

std::vector<std::string> GetListVariableBFM()
{
  std::vector<std::string> ListVar{"iOxyg", "iPO4_", "iNO3_", "iNH4_", "iO4n_", "iSiOH", "iN6r_", "iB1c_", "iB1n_", "iB1p_", "iP1c_", "iP1n_", "iP1p_", "iP1l_", "iP1s_", "iP2c_", "iP2n_", "iP2p_", "iP2l_", "iP3c_", "iP3n_", "iP3p_", "iP3l_", "iP4c_", "iP4n_", "iP4p_", "iP4l_", "iZ3c_", "iZ3n_", "iZ3p_", "iZ4c_", "iZ4n_", "iZ4p_", "iZ5c_", "iZ5n_", "iZ5p_", "iZ6c_", "iZ6n_", "iZ6p_", "iR1c_", "iR1n_", "iR1p_", "iR2c_", "iR3c_", "iR6c_", "iR6n_", "iR6p_", "iR6s_", "iO3c_", "iO3h_", "iDIC_", "iChlo"};
  return ListVar;
}

std::vector<std::string> GetListVariableDYE1()
{
  std::vector<std::string> ListVar{"dye_01"};
  return ListVar;
}



std::vector<std::string> GetListVariables(std::string const& eModelName)
{
  if (eModelName == "BFM_full")
    return GetListVariableBFM_full();
  if (eModelName == "BFM")
    return GetListVariableBFM();
  if (eModelName == "DYE1")
    return GetListVariableDYE1();
  std::cerr << "Error in GetListVariables\n";
  std::cerr << "eModelName=" << eModelName << "\n";
  std::cerr << "But only available model is BFM\n";
  throw TerminalException{1};
}




struct VarRomsDesc {
  std::string NetcdfName;
  std::string FullName;
  std::string Unit;
  double ScalarMult;
  std::string KeyStr;
  std::string shortStr;
};


std::vector<VarRomsDesc> GetFullVariablesNames(std::vector<std::string> const& ListVar, std::string const& VarInfoFile)
{
  int nbVar=ListVar.size();
  //
  std::cerr << "VarInfoFile=" << VarInfoFile << "\n";
  std::vector<std::string> ListLine = ReadFullFile(VarInfoFile);
  int nbLine=ListLine.size();
  //
  std::vector<int> LineHasString(nbLine);
  for (int iLine=0; iLine<nbLine; iLine++) {
    std::string eLine=ListLine[iLine];
    int HasStrOpCl=0;
    int len=eLine.size();
    std::string eOpCl="'";
    for (int i=0; i<len; i++) {
      std::string eSub = eLine.substr(i,1);
      if (eSub == eOpCl)
	HasStrOpCl=1;
    }
    LineHasString[iLine] = HasStrOpCl;
  }
  std::cerr << "We have ListHasString\n";
  //
  std::vector<int> ListLineFound(nbVar,-1);
  for (int iVar=0; iVar<nbVar; iVar++) {
    std::string eStrSearch = "idTvar(" + ListVar[iVar] + ")";
    std::cerr << "iVar=" << iVar << " / " << nbVar << " eStrSearch=" << eStrSearch << "\n";
    int iLineFound=-1;
    for (int iLine=0; iLine<nbLine; iLine++) {
      std::cerr << "Before split eLine=" << ListLine[iLine] << "\n";
      std::vector<std::string> LStr = STRING_Split(ListLine[iLine], eStrSearch);
      std::cerr << "After split\n";
      if (LStr.size() > 1)
	iLineFound=iLine;
    }
    if (iLineFound == -1) {
      std::cerr << "Failed to find a line containing the variable\n";
      std::cerr << "eStrSearch=" << eStrSearch << "\n";
      throw TerminalException{1};
    }
    ListLineFound[iVar] = iLineFound;
  }
  std::cerr << "We have ListLineFound\n";
  //
  auto ExtractSubString=[](std::string const& FullStr) -> std::string {
    bool IsInStr=false;
    std::string retStr;
    int eLen=FullStr.size();
    std::string DefStrChar="'";
    for (int iL=0; iL<eLen; iL++) {
      std::string eChar=FullStr.substr(iL,1);
      if (eChar != DefStrChar && IsInStr) {
	retStr += eChar;
      }
      if (eChar == DefStrChar) {
	IsInStr = !IsInStr;
      }
    }
    return retStr;
  };
  auto GetLineKeyString=[&](int const& iLine) -> std::string {
    if (LineHasString[iLine] == 0) {
      std::cerr << "The Line has no ' in it, therefore extraction is impossible\n";
      throw TerminalException{1};
    }
    return ExtractSubString(ListLine[iLine]);
  };
  //
  std::vector<VarRomsDesc> ListFullVarRomsDesc(nbVar);
  for (int iVar=0; iVar<nbVar; iVar++) {
    int iLineFound=ListLineFound[iVar];
    int iLineUnit=iLineFound - 3;
    int iLineFullName=iLineFound - 4;
    int iLineNetcdfName=iLineFound - 5;
    //
    std::string Unit = GetLineKeyString(iLineUnit);
    std::string FullName = GetLineKeyString(iLineFullName);
    std::string NetcdfName = GetLineKeyString(iLineNetcdfName);
    std::string KeyStr = "idTvar(" + ListVar[iVar] + ")";
    std::string shortStr = ListVar[iVar];
    //
    std::string scalStr = ListLine[iLineFound + 2];
    std::vector<std::string> LStr = STRING_Split(scalStr, "d");
    if (LStr.size() != 2) {
      std::cerr << "The line should be of something like 1.0d0\n";
      std::cerr << "Instead we have scalStr=" << scalStr << "\n";
      throw TerminalException{1};
    }
    std::string eStr=LStr[0];
    std::istringstream is(eStr);
    double ScalarMult;
    is >> ScalarMult;
    ListFullVarRomsDesc[iVar] = {NetcdfName, FullName, Unit, ScalarMult, KeyStr, shortStr};
  }
  return ListFullVarRomsDesc;
}


FullNamelist NAMELIST_SET_VARIABLE_INITIAL_ROMS()
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
  ListStringValues1["NetcdfInitialFile"]="unset";
  ListStringValues1["PrefixVariableDefinitions"]="unset";
  ListStringValues1["TracerModelName"] = "unset";
  ListStringValues1["VarInfoFile"] = "External/varinfo.dat";
  ListStringValues1["GridFile"] = "unset";
  ListStringValues1["FileDescARVD"] = "unset";
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
  return {std::move(ListBlock), "undefined"};
}

Eigen::Tensor<double,3> GetConditionsAccordingToDescription(GridArray const& GrdArr, ARVDtyp const& ARVD, FullNamelist const& eFullDesc, VarRomsDesc const& eVarRomsDesc)
{
  int N = ARVD.N;
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  Eigen::Tensor<double,3> eTens = ZeroTensor3<double>(N, eta_rho, xi_rho);
  SingleBlock BlDESC = eFullDesc.ListBlock.at("DESCRIPTION");
  std::string MethodSetting = BlDESC.ListStringValues.at("MethodSetting");
  if (MethodSetting == "Constant") {
    double eVal = BlDESC.ListDoubleValues.at("ConstantValue");
    for (int i=0; i<N; i++)
      for (int j=0; j<eta_rho; j++)
	for (int k=0; k<xi_rho; k++)
	  eTens(i,j,k) = eVal;
    return eTens;
  }
  if (MethodSetting == "VerticalProfile") {
    std::vector<double> ListDep = BlDESC.ListListDoubleValues.at("ListVerticalLevels");
    std::vector<double> ListVal = BlDESC.ListListDoubleValues.at("ListValues");
    MyVector<double> ListDep_V = VectorFromStdVector(ListDep);
    MyVector<double> ListVal_V = VectorFromStdVector(ListVal);
    std::vector<double> VectZ(eta_rho*xi_rho*N);
    int idx=0;
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	double eZeta=0;
	double eDep=GrdArr.GrdArrRho.DEP(i,j);
	MyVector<double> Zr_out=GetVertCoord_R(ARVD, eDep, eZeta);
	for (int k=0; k<N; k++) {
	  VectZ[idx] = Zr_out(k);
	  idx++;
	}
      }
    std::vector<double> retListVal = OneDimInterpolation_vector(ListVal_V, ListDep_V, VectZ);
    idx=0;
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	for (int k=0; k<N; k++) {
	  eTens(k,i,j) = retListVal[idx];
	  idx++;
	}
      }
    return eTens;
  }
  std::cerr << "Error in GetConditionsAccordingToDescription\n";
  std::cerr << "MethodSetting=" << MethodSetting << " has not been covered\n";
  throw TerminalException{1};
}



MyMatrix<double> My_Psi2Rho(MyMatrix<double> const& F_psi)
{
  int eta_psi = F_psi.rows();
  int xi_psi = F_psi.cols();
  int eta_rho = eta_psi + 1;
  int xi_rho = xi_psi + 1;
  MyMatrix<double> F_rho(eta_rho, xi_rho);
  for (int iEta=1; iEta<eta_psi; iEta++)
    for (int iXi=1; iXi<xi_psi; iXi++) {
      double val1 = F_psi(iEta-1,iXi-1);
      double val2 = F_psi(iEta-1,iXi);
      double val3 = F_psi(iEta,iXi-1);
      double val4 = F_psi(iEta,iXi);
      F_rho(iEta, iXi) = (val1 + val2 + val3 + val4) / 4;
    }
  for (int iEta=1; iEta<eta_psi; iEta++) {
    F_rho(iEta,0) = F_rho(iEta,1);
    F_rho(iEta,xi_rho-1) = F_rho(iEta,xi_rho-2);
  }
  for (int iXi=1; iXi<xi_psi; iXi++) {
    F_rho(0,iXi) = F_rho(1,iXi);
    F_rho(eta_rho-1,iXi) = F_rho(eta_rho-2, iXi);
  }
  F_rho(0,0) = (F_rho(0,1) + F_rho(1,0)) / 2;
  F_rho(eta_rho-1,0) = (F_rho(eta_rho-1,1) + F_rho(eta_rho-2,0)) / 2;
  F_rho(0,xi_rho-1) = (F_rho(0,xi_rho-2) + F_rho(1,xi_rho-1)) / 2;
  F_rho(eta_rho-1,xi_rho-1) = (F_rho(eta_rho-1,xi_rho-2) + F_rho(eta_rho-2,xi_rho-1)) / 2;
  return F_rho;
}



void print_matrix_info(std::string name, MyMatrix<double> const& M)
{
  double min_c = M.minCoeff();
  double max_c = M.maxCoeff();
  std::cerr << "name=" << name << " |M|=" << M.rows() << " / " << M.cols() << " min_c=" << min_c << " max_c=" << max_c << "\n";
}

MyMatrix<double> ZeroMask(MyMatrix<double> const& M, MyMatrix<uint8_t> const& msk)
{
  int eta_rho = M.rows();
  int xi_rho = M.cols();
  MyMatrix<double> Mret = M;
  for (int iEta=0; iEta<eta_rho; iEta++)
    for (int iXi=0; iXi<xi_rho; iXi++)
      if (msk(iEta,iXi) == 0)
        Mret(iEta, iXi) = 0;
  return Mret;
}



MyMatrix<double> GRID_VorticityRho(const GridArray& GrdArr, MyMatrix<double> const& U, MyMatrix<double> const& V)
{
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int eta_u = GrdArr.GrdArrU.LON.rows();
  int xi_u = GrdArr.GrdArrU.LON.cols();
  int eta_v = GrdArr.GrdArrV.LON.rows();
  int xi_v = GrdArr.GrdArrV.LON.cols();
  int eta_psi = eta_rho - 1;
  int xi_psi = xi_rho - 1;
  //
  const MyMatrix<double> & pm = GrdArr.GrdArrRho.pm;
  const MyMatrix<double> & pn = GrdArr.GrdArrRho.pn;
  print_matrix_info("pm", pm);
  print_matrix_info("pn", pn);
  print_matrix_info("U", U);
  print_matrix_info("V", V);
  MyMatrix<double> U_zero = ZeroMask(U, GrdArr.GrdArrU.MSK);
  MyMatrix<double> V_zero = ZeroMask(V, GrdArr.GrdArrV.MSK);
  print_matrix_info("U_zero", U_zero);
  print_matrix_info("V_zero", V_zero);
  MyMatrix<double> uom(eta_u, xi_u);
  for (int iEta=0; iEta<eta_u; iEta++)
    for (int iXi=0; iXi<xi_u; iXi++) {
      double s_pm = pm(iEta, iXi) + pm(iEta, iXi+1);
      uom(iEta, iXi) = 2 * U_zero(iEta,iXi) / s_pm;
    }
  print_matrix_info("uom", uom);
  MyMatrix<double> von(eta_v, xi_v);
  for (int iEta=0; iEta<eta_v; iEta++)
    for (int iXi=0; iXi<xi_v; iXi++) {
      double s_pn = pn(iEta, iXi) + pn(iEta+1, iXi);
      von(iEta, iXi) = 2 * V_zero(iEta,iXi) / s_pn;
    }
  print_matrix_info("von", von);
  //
  MyMatrix<double> mn_psi(eta_psi, xi_psi);
  for (int iEta=0; iEta<eta_psi; iEta++)
    for (int iXi=0; iXi<xi_psi; iXi++) {
      double val1 = pm(iEta,iXi) * pn(iEta,iXi);
      double val2 = pm(iEta,iXi+1) * pn(iEta,iXi+1);
      double val3 = pm(iEta+1,iXi) * pn(iEta+1,iXi);
      double val4 = pm(iEta+1,iXi+1) * pn(iEta+1,iXi+1);
      mn_psi(iEta, iXi) = (val1 + val2 + val3 + val4) / 4;
    }
  MyMatrix<double> TheVortPsi(eta_psi, xi_psi);
  for (int iEta=0; iEta<eta_psi; iEta++)
    for (int iXi=0; iXi<xi_psi; iXi++) {
      double val1 = von(iEta, iXi+1);
      double val2 = von(iEta, iXi);
      double val3 = uom(iEta+1, iXi);
      double val4 = uom(iEta, iXi);
      TheVortPsi(iEta,iXi) = mn_psi(iEta, iXi) * (val1 - val2 - val3 + val4);
    }
  MyMatrix<double> TheVort = My_Psi2Rho(TheVortPsi);
  return ZeroMask(TheVort, GrdArr.GrdArrRho.MSK);
}



void SetNetcdfInitial(FullNamelist const& eFull)
{
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  std::string NetcdfInitialFile = BlPROC.ListStringValues.at("NetcdfInitialFile");
  std::string PrefixVariableDefinitions = BlPROC.ListStringValues.at("PrefixVariableDefinitions");
  std::string TracerModelName = BlPROC.ListStringValues.at("TracerModelName");
  std::string VarInfoFile = BlPROC.ListStringValues.at("VarInfoFile");
  std::string GridFile = BlPROC.ListStringValues.at("GridFile");
  std::string FileDescARVD = BlPROC.ListStringValues.at("FileDescARVD");
  std::cerr << "The input file has been read\n";
  //
  std::vector<std::string> ListVar = GetListVariables(TracerModelName);
  std::cerr << "ListVar read\n";
  std::vector<VarRomsDesc> ListVarRomsDesc = GetFullVariablesNames(ListVar, VarInfoFile);
  std::cerr << "ListVarRomsDesc read\n";
  //
  ARVDtyp eARVD = ReadROMSverticalStratification(FileDescARVD);
  std::cerr << "eARVD read\n";
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  int N=eARVD.N;
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  int eDimTracer = N * eta_rho * xi_rho;
  std::cerr << "GrdArr and related read\n";
  //
  if (!IsExistingFile(NetcdfInitialFile)) {
    std::cerr << "Error the file NetcdfInitialFile=" << NetcdfInitialFile << " is missing\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(NetcdfInitialFile, netCDF::NcFile::write);
  std::vector<std::string> ListDimField{"ocean_time", "s_rho", "eta_rho", "xi_rho"};
  std::vector<double> A(eDimTracer);
  for (auto & eVarRomsDesc : ListVarRomsDesc) {
    std::string shortStr = eVarRomsDesc.shortStr;
    std::string eDescFile = PrefixVariableDefinitions + shortStr + ".nml";
    std::cerr << "Treating variable shortStr=" << shortStr << " eDescFile=" << eDescFile << " NetcdfName=" << eVarRomsDesc.NetcdfName << "\n";
    if (IsExistingFile(eDescFile)) {
      FullNamelist eFullDesc = Individual_Tracer_Variable_File();
      NAMELIST_ReadNamelistFile(eDescFile, eFullDesc);
      std::cerr << "eDescFile read\n";
      Eigen::Tensor<double,3> eTensTracer = GetConditionsAccordingToDescription(GrdArr, eARVD, eFullDesc, eVarRomsDesc);
      std::cerr << "eTensTracer read\n";
      std::vector<double> Vstat = GetMinMaxAvg(eTensTracer);
      std::cerr << "min=" << Vstat[0] << " max=" << Vstat[1] << " avg=" << Vstat[2] << "\n";
      if (!NC_IsVar(NetcdfInitialFile, eVarRomsDesc.NetcdfName)) {
        netCDF::NcVar eVarData = dataFile.addVar(eVarRomsDesc.NetcdfName, "float", ListDimField);
        eVarData.putAtt("long_name", eVarRomsDesc.FullName);
        eVarData.putAtt("units", eVarRomsDesc.Unit);
      }
      int idx=0;
      for (int k=0; k<N; k++)
        for (int i=0; i<eta_rho; i++)
          for (int j=0; j<xi_rho; j++) {
            A[idx] = eTensTracer(k,i,j);
            idx++;
          }
      netCDF::NcVar fVarData = dataFile.getVar(eVarRomsDesc.NetcdfName);
      fVarData.putVar(A.data());
    }
  }
}


FullNamelist NAMELIST_CREATE_DEFAULT_SETTING()
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
  ListStringValues1["PrefixVariableDefinitions"]="unset";
  ListStringValues1["TracerModelName"] = "unset";
  ListStringValues1["VarInfoFile"] = "External/varinfo.dat";
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
  return {std::move(ListBlock), "undefined"};
}



FullNamelist NAMELIST_ROMS_FIELD_COPY()
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
  ListStringValues1["PrefixVariableDefinitions"]="unset";
  ListStringValues1["TracerModelName"] = "unset";
  ListStringValues1["VarInfoFile"] = "External/varinfo.dat";
  ListStringValues1["NetcdfFileIn"] = "unset";
  ListStringValues1["NetcdfFileOut"] = "unset";
  ListStringValues1["RomsGridFile"] = "unset";
  ListIntValues1["NetcdfIdxIn"] = -1;
  ListIntValues1["NetcdfIdxOut"] = -1;
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
  return {std::move(ListBlock), "undefined"};
}



void CopyTracerFields(FullNamelist const& eFull)
{
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  std::string NetcdfFileIn = BlPROC.ListStringValues.at("NetcdfFileIn");
  int NetcdfIdxIn = BlPROC.ListIntValues.at("NetcdfIdxIn");
  std::string NetcdfFileOut = BlPROC.ListStringValues.at("NetcdfFileOut");
  int NetcdfIdxOut = BlPROC.ListIntValues.at("NetcdfIdxOut");
  std::string GridFile = BlPROC.ListStringValues.at("RomsGridFile");
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  //
  std::string PrefixVariableDefinitions = BlPROC.ListStringValues.at("PrefixVariableDefinitions");
  std::string TracerModelName = BlPROC.ListStringValues.at("TracerModelName");
  std::string VarInfoFile = BlPROC.ListStringValues.at("VarInfoFile");

  std::vector<std::string> ListVar = GetListVariables(TracerModelName);
  std::vector<VarRomsDesc> ListVarRomsDesc = GetFullVariablesNames(ListVar, VarInfoFile);
  for (auto & eVarRomsDesc : ListVarRomsDesc) {
    std::string FullFile = PrefixVariableDefinitions + eVarRomsDesc.shortStr + ".nml";
    //
    std::string eVarNetcdf = eVarRomsDesc.NetcdfName;
    std::cerr << "Treating eVarNetcdf = " << eVarNetcdf << "\n";
    //
    Eigen::Tensor<double,3> TheTens = NETCDF_Get3DvariableSpecEntry_FD(NetcdfFileIn, GrdArr, eVarNetcdf, NetcdfIdxIn);
    NETCDF_Write3DvariableSpecEntry(NetcdfFileOut, eVarNetcdf, NetcdfIdxOut, TheTens);
  }
}




void CreateDefaultInputFiles(FullNamelist const& eFull)
{
  SingleBlock BlPROC = eFull.ListBlock.at("PROC");
  std::string PrefixVariableDefinitions = BlPROC.ListStringValues.at("PrefixVariableDefinitions");
  std::string TracerModelName = BlPROC.ListStringValues.at("TracerModelName");
  std::string VarInfoFile = BlPROC.ListStringValues.at("VarInfoFile");

  std::vector<std::string> ListVar = GetListVariables(TracerModelName);
  std::vector<VarRomsDesc> ListVarRomsDesc = GetFullVariablesNames(ListVar, VarInfoFile);
  for (auto & eVarRomsDesc : ListVarRomsDesc) {
    std::string FullFile = PrefixVariableDefinitions + eVarRomsDesc.shortStr + ".nml";
    std::ofstream os(FullFile);
    os << "&DESCRIPTION\n";
    os << " VariableName = \"" << eVarRomsDesc.NetcdfName << "\",\n";
    os << " MethodSetting = \"Constant\",\n";
    os << " ConstantValue = 0.0,\n";
    os << "/\n";
  }
}



#endif
