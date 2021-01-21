#ifndef ROMS_FUNCTIONALITY_INCLUDE
#define ROMS_FUNCTIONALITY_INCLUDE


#include "Basic_netcdf.h"
#include "OneDimInterpolation.h"
#include "Namelist.h"
#include "Model_grids.h"







MyMatrix<double> My_u2rho(MyMatrix<double> const& eVar_u, MyMatrix<int> const& MSK_u)
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
      }
      else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}


MyMatrix<double> My_v2rho(MyMatrix<double> const& eVar_v, MyMatrix<int> const& MSK_v)
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
      }
      else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}



Eigen::Tensor<double,3> My_v2rho_3D(Eigen::Tensor<double,3> const& eVar_v, MyMatrix<int> const& MSK_v)
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
      }
      else {
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
      }
      else {
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
      }
      else {
	for (int k=0; k<s_vert; k++)
	  Mout(k, i, j)=0;
      }
    }
  return Mout;
}


Eigen::Tensor<double,3> My_u2rho_3D(Eigen::Tensor<double,3> const& eVar_u, MyMatrix<int> const& MSK_u)
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
      }
      else {
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
  MyMatrix<int> MSK_rho;
  MyMatrix<int> MSK_u;
  MyMatrix<int> MSK_v;
  MyMatrix<int> MSK_psi;
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
  eA.MSK_rho=ConvertMatrixUniversal<int,double>(MSK_rho_double);
  //
  MyMatrix<double> MSK_u_double=NC_Read2Dvariable(eFile, "mask_u");
  eA.MSK_u=ConvertMatrixUniversal<int,double>(MSK_u_double);
  //
  MyMatrix<double> MSK_v_double=NC_Read2Dvariable(eFile, "mask_v");
  eA.MSK_v=ConvertMatrixUniversal<int,double>(MSK_v_double);
  //
  MyMatrix<double> MSK_psi_double=NC_Read2Dvariable(eFile, "mask_psi");
  eA.MSK_psi=ConvertMatrixUniversal<int,double>(MSK_psi_double);
  //
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data_Vtrans=dataFile.getVar("Vtransform");
  netCDF::NcVar data_Vstret=dataFile.getVar("Vstretching");
  netCDF::NcVar data_theta_s=dataFile.getVar("theta_s");
  netCDF::NcVar data_theta_b=dataFile.getVar("theta_b");
  netCDF::NcVar data_Tcline=dataFile.getVar("Tcline");
  netCDF::NcVar data_hc=dataFile.getVar("hc");
  //
  int *eValI;
  eValI=new int[1];
  double *eValD;
  eValD=new double[1];
  //
  data_Vtrans.getVar(eValI);
  eA.ARVD.Vtransform=*eValI;
  data_Vstret.getVar(eValI);
  eA.ARVD.Vstretching=*eValI;
  data_theta_s.getVar(eValD);
  eA.ARVD.theta_s=*eValD;
  data_theta_b.getVar(eValD);
  eA.ARVD.theta_b=*eValD;
  data_Tcline.getVar(eValD);
  eA.ARVD.Tcline=*eValD;
  data_hc.getVar(eValD);
  eA.ARVD.hc=*eValD;
  //
  delete [] eValI;
  delete [] eValD;
  //
  double *eVarR, *eVarW;
  netCDF::NcVar data_Cs_r=dataFile.getVar("Cs_r");
  netCDF::NcVar data_Cs_w=dataFile.getVar("Cs_w");
  netCDF::NcVar data_s_r=dataFile.getVar("s_rho");
  netCDF::NcVar data_s_w=dataFile.getVar("s_w");
  netCDF::NcDim eDim=data_Cs_r.getDim(0);
  int dim_s_r=eDim.getSize();
  netCDF::NcDim fDim=data_Cs_w.getDim(0);
  int dim_s_w=fDim.getSize();
  eVarR=new double[dim_s_r];
  eVarW=new double[dim_s_w];
  data_Cs_r.getVar(eVarR);
  for (int i=0; i<dim_s_r; i++)
    eA.ARVD.Cs_r(i)=eVarR[i];
  data_s_r.getVar(eVarR);
  for (int i=0; i<dim_s_r; i++)
    eA.ARVD.sc_r(i)=eVarR[i];
  data_Cs_w.getVar(eVarW);
  for (int i=0; i<dim_s_w; i++)
    eA.ARVD.Cs_w(i)=eVarW[i];
  data_s_w.getVar(eVarW);
  for (int i=0; i<dim_s_w; i++)
    eA.ARVD.sc_w(i)=eVarW[i];
  delete [] eVarR;
  delete [] eVarW;
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

std::vector<std::string> GetListVariableBFM()
{
  std::vector<std::string> ListVar{"iOxyg", "iPO4_", "iNO3_", "iNH4_", "iO4n_", "iSiOH", "iN6r_", "iB1c_", "iB1n_", "iB1p_", "iP1c_", "iP1n_", "iP1p_", "iP1l_", "iP1s_", "iP2c_", "iP2n_", "iP2p_", "iP2l_", "iP3c_", "iP3n_", "iP3p_", "iP3l_", "iP4c_", "iP4n_", "iP4p_", "iP4l_", "iZ3c_", "iZ3n_", "iZ3p_", "iZ4c_", "iZ4n_", "iZ4p_", "iZ5c_", "iZ5n_", "iZ5p_", "iZ6c_", "iZ6n_", "iZ6p_", "iR1c_", "iR1n_", "iR1p_", "iR2c_", "iR3c_", "iR6c_", "iR6n_", "iR6p_", "iR6s_", "iO3c_", "iO3h_", "iEIR_", "iDIC_", "iChlo", "siP1_", "siP2_", "siP3_", "siP4_", "eiP1_", "eiP2_", "eiP3_", "eiP4_", "ruPTc", "ruZTc", "ixEPS"};
  return ListVar;
}

std::vector<std::string> GetListVariableDYE1()
{
  std::vector<std::string> ListVar{"dye_01"};
  return ListVar;
}



std::vector<std::string> GetListVariables(std::string const& eModelName)
{
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
  std::vector<std::string> ListLine = ReadFullFile(VarInfoFile);
  int nbLine=ListLine.size();
  std::vector<int> LineHasString(nbLine);
  //
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
  //
  std::vector<int> ListLineFound(nbVar,-1);
  for (int iVar=0; iVar<nbVar; iVar++) {
    std::string eStrSearch = "idTvar(" + ListVar[iVar] + ")";
    int iLineFound=-1;
    for (int iLine=0; iLine<nbLine; iLine++) {
      std::vector<std::string> LStr = STRING_Split(ListLine[iLine], eStrSearch);
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
  if (MethodSetting == "ConstantValue") {
    double eVal = BlDESC.ListDoubleValues.at("ConstantValue");
    for (int i=0; i<N; i++)
      for (int j=0; j<eta_rho; j++)
	for (int k=0; k<xi_rho; k++)
	  eTens(i,j,k) = eVal;
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
  }
  return eTens;
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
  //
  std::vector<std::string> ListVar = GetListVariables(TracerModelName);
  std::vector<VarRomsDesc> ListVarRomsDesc = GetFullVariablesNames(ListVar, VarInfoFile);
  //
  ARVDtyp eARVD = ReadROMSverticalStratification(FileDescARVD);
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  int N=eARVD.N;
  int eta_rho=GrdArr.GrdArrRho.LON.rows();
  int xi_rho =GrdArr.GrdArrRho.LON.cols();
  int eDimTracer = N * eta_rho * xi_rho;
  //
  if (!IsExistingFile(NetcdfInitialFile)) {
    std::cerr << "Error the file NetcdfInitialFile=" << NetcdfInitialFile << " is missing\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(NetcdfInitialFile, netCDF::NcFile::write);
  std::vector<std::string> ListDimField{"ocean_time", "s_rho", "eta_rho", "xi_rho"};
  double *A;
  A = new double[eDimTracer];
  for (auto & eVarRomsDesc : ListVarRomsDesc) {
    std::string shortStr = eVarRomsDesc.shortStr;
    std::cerr << "Treating variable shortStr=" << shortStr << "\n";
    std::string eDescFile = PrefixVariableDefinitions + shortStr + ".nml";
    FullNamelist eFullDesc = Individual_Tracer_Variable_File();
    NAMELIST_ReadNamelistFile(eDescFile, eFullDesc);
    
    Eigen::Tensor<double,3> eTensTracer = GetConditionsAccordingToDescription(GrdArr, eARVD, eFullDesc, eVarRomsDesc);
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
    fVarData.putVar(A);
  }
  delete [] A;
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
