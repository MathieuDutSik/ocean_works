#ifndef DEFINE_SEICHE_NUMERIC_EIGENVALUE
#include "MAT_Matrix.h"
#include "Model_grids.h"
#include "Namelist.h"
// Spectra code is available at https://github.com/yixuan/spectra
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>


FullNamelist NAMELIST_SEICHE_Eigen()
{
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["GridFile"]="unset";
  ListStringValues1["OutFile"]="unset_out";
  ListIntValues1["maxNbEig"] = 20;
  ListIntValues1["ncv"] = 20;
  ListBoolValues1["RescaleEigenvector"] = true;
  ListBoolValues1["UseSI_Bmatrix"] = false;
  ListDoubleValues1["h0"] = 0;
  SingleBlock BlockCOMP;
  BlockCOMP.ListDoubleValues=ListDoubleValues1;
  BlockCOMP.ListIntValues=ListIntValues1;
  BlockCOMP.ListBoolValues=ListBoolValues1;
  BlockCOMP.ListStringValues=ListStringValues1;
  ListBlock["COMP"]=BlockCOMP;
  //
  return {ListBlock, "undefined"};
}





struct PeriodicSolution {
  double lambda;
  double period;
  MyVector<double> Height;
};


std::string GetPairVectorInfo(MyVector<double> const& V1, MyVector<double> const& V2)
{
  double scal11 = V1.dot(V1);
  double scal12 = V2.dot(V1);
  double scal22 = V2.dot(V2);
  double scal_N = scal12 / sqrt(scal11 * scal22);
  double quot = scal22 / scal11;
  return "norm1=" + std::to_string(scal11) + " norm2=" + std::to_string(scal22) + " scal_N=" + std::to_string(scal_N) + " quot=" + std::to_string(quot);
}



/*
  The potential function considered is
  int_{Omega} (h0 - b) [(df/dx)^2 + (df/dy)^2]
  ---
  Here h0 is the surface level. b is the height and h0 - b is the total water height.
  Thus we set b = -h with h the depth.
  ---
  An eigenvalue lambda of the quadratic form will give us via
  lambda = omega^2 / g
 */
std::vector<PeriodicSolution> ComputeEigenvaluesSWE1(double const& h0, int const& maxNbEig, GridArray const& GrdArr, bool const& RescaleEigenvector, bool const& UseSI_Bmatrix, int const& ncv)
{
  int nb_elt = GrdArr.INE.rows();
  MyMatrix<double> X =  GrdArr.GrdArrRho.LON;
  MyMatrix<double> Y =  GrdArr.GrdArrRho.LAT;
  MyMatrix<double> B = -GrdArr.GrdArrRho.DEP;
  std::cerr << "X : max=" << X.maxCoeff() << " min=" << X.minCoeff() << "\n";
  std::cerr << "Y : max=" << Y.maxCoeff() << " min=" << Y.minCoeff() << "\n";
  int nb_point = X.size();
  //
  MyVector<double> ListArea(nb_elt);
  double Tot_Area = 0;
  for (int i_elt=0; i_elt<nb_elt; i_elt++) {
    int i0 = GrdArr.INE(i_elt,0);
    int i1 = GrdArr.INE(i_elt,1);
    int i2 = GrdArr.INE(i_elt,2);
    double x0 = X(i0,0);
    double x1 = X(i1,0);
    double x2 = X(i2,0);
    double y0 = Y(i0,0);
    double y1 = Y(i1,0);
    double y2 = Y(i2,0);
    double area_by_2 = (x2 - x0) * (y1 - y0) - (y2 - y0) * (x1 - x0);
    double area = std::abs(area_by_2 / 2);
    Tot_Area += area;
    ListArea(i_elt) = area;
  }
  std::cerr << "Tot_Area=" << Tot_Area << "\n";
  MyVector<double> SI = ZeroVector<double>(nb_point);
  MyVector<double> Coeff = ZeroVector<double>(nb_point);
  for (int i_elt=0; i_elt<nb_elt; i_elt++) {
    for (int idx=0; idx<3; idx++) {
      int ip = GrdArr.INE(i_elt, idx);
      SI(ip) += ListArea(i_elt) / 3;
    }
  }
  if (UseSI_Bmatrix) {
    for (int i_pt=0; i_pt<nb_point; i_pt++)
      Coeff(i_pt) = sqrt(SI(i_pt));
  }
  else {
    for (int i_pt=0; i_pt<nb_point; i_pt++)
      Coeff(i_pt) = 1;
  }
  GraphSparseImmutable eG = GetUnstructuredVertexAdjInfo(GrdArr.INE, nb_point);
  std::pair<std::vector<size_t>, std::vector<size_t>> ePair = eG.Get_ListStart_ListListAdj();
  int nb_adj = ePair.second.size();
  //
  // Now building the matrix
  //
  std::vector<double> ListDiagValue_A_SI(nb_point, 0);
  std::vector<double> ListDiagValue_A(nb_point, 0);
  std::vector<double> ListDiagValue_B(nb_point, 0);
  std::vector<double> ListDiagValue_B_SI(nb_point, 0);
  std::vector<double> ListOffDiagValue_A_SI(nb_adj, 0);
  std::vector<double> ListOffDiagValue_A(nb_adj, 0);
  std::vector<double> ListOffDiagValue_B(nb_adj, 0);
  double minWaterHeight = std::numeric_limits<double>::max();
  double maxWaterHeight = std::numeric_limits<double>::min();
  double sumDelta=0;
  std::vector<int> ListNbMatch(nb_adj, 0);
  for (int i_elt=0; i_elt<nb_elt; i_elt++) {
    int i0 = GrdArr.INE(i_elt,0);
    int i1 = GrdArr.INE(i_elt,1);
    int i2 = GrdArr.INE(i_elt,2);
    for (int idx=0; idx<3; idx++) {
      int i = GrdArr.INE(i_elt,idx);
      double e_val = h0 - B(i);
      minWaterHeight = std::min(minWaterHeight, e_val);
      maxWaterHeight = std::max(maxWaterHeight, e_val);
    }
    /*  The linear function is f(x,y) = a + bx + cy over the triangle.
        f0 = f(x0,y0) = a + b x0 + c y0
        f1 = f(x1,y1) = a + b x1 + c y1
        f2 = f(x2,y2) = a + b x2 + c y2
        This gives
        f1 - f0 = b (x1 - x0) + c (y1 - y0)
        f2 - f0 = b (x2 - x0) + c (y2 - y0)
        which gives the values
        b [(x1 - x0) (y2 - y0) - (x2 - x0) (y1 - y0)] = (y2 - y0) (f1 - f0) - (y1 - y0) (f2 - f0)
        c [(x1 - x0) (y2 - y0) - (x2 - x0) (y1 - y0)] = (x1 - x0) (f2 - f0) - (x2 - x0) (f1 - f0)
        Define Delta = (x1 - x0) (y2 - y0) - (x2 - x0) (y1 - y0)
        We thus have:
        b Delta = f0 (y1 - y2) + f1 (y2 - y0) + f2 (y0 - y1) = b0 f0 + b1 f1 + b2 f2
        c Delta = f0 (x2 - x1) + f1 (x0 - x2) + f2 (x1 - x0) = c0 f0 + c1 f1 + c2 f2
        The gradient term to put is thus:
        avg(h - B) * Area(Delta) * [b*b + c*c]
        The integral term \int_T f^2 is thus
        P = P0 + u (P1 - P0) + v (P2 - P0) with 0\leq u,v and u+v \leq 1
        Triangle expressed as 0\leq x,y and x+y\leq 1.
        I(f) = \int_{u=0}^{u=1}\int_{v=0}^{v=1-u} f(u,v)
        Thus I(1) = \int_{u=0}^{u=1} (1-u) = 1 - 1/2 = 1/2
        I(u) = \int_{u=0}^{u=1}\int_{v=0}^{v=1-u} u
             = \int_{u=0}^{u=1} u-u^2 = 1/2 - 1/3 = 1/6
        I(u^2) = 1/3 - 1/4 = 1/12
        I(v) = \int_{u=0}^{u=1}\int_{v=0}^{v=1-u} v
             = \int_{u=0}^{u=1} (1-u)^2 / 2 = \int_{u=0}^{u=1} u^2 / 2 = 1/6
        I(uv) = \int_{u=0}^{u=1}\int_{v=0}^{v=1-u} uv
              = \int_{u=0}^{u=1}u (1-u)^2/2 = \int_{u=0}^{u=1}u^2 (1-u)/2 = 
              = (1/3 - 1/4)/2 = 1/24
        I = int_T (a + bx + cy)^2
          = int_T a^2 + 2ab x + 2ac y + 2bc xy + b^2 x^2 + c^2 y^2
          = a^2 / 2 + ab / 3 + ac / 3 + bc / 12 + b^2 / 12 + c^2 / 12
          = [6a^2 + 4 ab + 4 ac + bc + b^2 + c^2] / 12
        I = [6f0^2 + 4f0 (f1 - f0) + 4f0 (f2 - f0) + (f1-f0)(f2- f0) + (f1-f0)^2 + (f2-f0)^2] / 12
          = [f0^2 + f1^2 + f2^2 + f0 f1 + f0 f2 + f1 f2] / 12
        The problem of this is that it requires solving a generalized eigenvalue problem.
        The dimensionality of the eigenvalues is [L^-1].
    */
    double val0 = h0 - B(i0);
    double val1 = h0 - B(i1);
    double val2 = h0 - B(i2);
    double avg_val = (val0 + val1 + val2) / 3.0;
    double b0 = Y(i1,0) - Y(i2,0);
    double b1 = Y(i2,0) - Y(i0,0);
    double b2 = Y(i0,0) - Y(i1,0);
    double c0 = X(i2,0) - X(i1,0);
    double c1 = X(i0,0) - X(i2,0);
    double c2 = X(i1,0) - X(i0,0);
    double delta = c2 * b1 - c1 * b2;
    sumDelta += delta;
    int idx01 = eG.GetIndex(i0, i1);
    int idx10 = eG.GetIndex(i1, i0);
    int idx02 = eG.GetIndex(i0, i2);
    int idx20 = eG.GetIndex(i2, i0);
    int idx12 = eG.GetIndex(i1, i2);
    int idx21 = eG.GetIndex(i2, i1);
    //    std::cerr << "i_elt=" << i_elt << " idx01=" << idx01 << " idx10=" << idx10 << "\n";
    //    std::cerr << "idx02=" << idx02 << " idx20=" << idx20 << " idx12=" << idx12 << " idx21=" << idx21 << "\n";
    // Construction of A matrix
    double tmp;
    double alpha = avg_val * ListArea(i_elt) / (delta * delta);
    // dimensionality of alpha is [L]
    ListDiagValue_A[i0] += (b0 * b0 + c0 * c0) * alpha;
    ListDiagValue_A[i1] += (b1 * b1 + c1 * c1) * alpha;
    ListDiagValue_A[i2] += (b2 * b2 + c2 * c2) * alpha;
    tmp = (b0 * b1 + c0 * c1) * alpha;
    ListOffDiagValue_A[idx01] += tmp;
    ListOffDiagValue_A[idx10] += tmp;
    tmp = (b0 * b2 + c0 * c2) * alpha;
    ListOffDiagValue_A[idx02] += tmp;
    ListOffDiagValue_A[idx20] += tmp;
    tmp = (b1 * b2 + c1 * c2) * alpha;
    ListOffDiagValue_A[idx12] += tmp;
    ListOffDiagValue_A[idx21] += tmp;
    // Construction of B matrix
    double beta = ListArea(i_elt) / 6;
    ListDiagValue_B[i0] += beta;
    ListDiagValue_B[i1] += beta;
    ListDiagValue_B[i2] += beta;
    ListOffDiagValue_B[idx01] += beta;
    ListOffDiagValue_B[idx10] += beta;
    ListOffDiagValue_B[idx02] += beta;
    ListOffDiagValue_B[idx20] += beta;
    ListOffDiagValue_B[idx12] += beta;
    ListOffDiagValue_B[idx21] += beta;
    ListNbMatch[idx01]++;
    ListNbMatch[idx10]++;
    ListNbMatch[idx02]++;
    ListNbMatch[idx20]++;
    ListNbMatch[idx12]++;
    ListNbMatch[idx21]++;
  }
  //  for (int i_adj=0; i_adj<nb_adj; i_adj++)
  //    std::cerr << "i_adj=" << i_adj << " ListNbMatch[i_adj]=" << ListNbMatch[i_adj] << "\n";
  std::cerr << "minWaterHeight = " << minWaterHeight << "  maxWaterHeight = " << maxWaterHeight << "\n";
  std::cerr << "sumDelta=" << sumDelta << "\n";
  using Ttrip = Eigen::Triplet<double>;
  std::vector<Ttrip> tripletList_A_SI, tripletList_A, tripletList_B, tripletList_B_SI;
  double e_diag, f_diag, e_si, f_si, e_val, f_val;
  for (int i_pt=0; i_pt<nb_point; i_pt++) {
    int eStart = ePair.first[i_pt];
    int eEnd = ePair.first[i_pt+1];
    // The A SI scaled matrix
    e_diag = ListDiagValue_A[i_pt];
    e_si = Coeff(i_pt);
    f_diag = e_diag / (e_si * e_si);
    tripletList_A_SI.push_back({i_pt, i_pt, f_diag});
    for (int i=eStart; i<eEnd; i++) {
      int j_pt = ePair.second[i];
      f_si = Coeff(j_pt);
      e_val = ListOffDiagValue_A[i];
      f_val = e_val / (e_si * f_si);
      tripletList_A_SI.push_back({i_pt, j_pt, f_val});
    }
    // The A matrix
    e_diag = ListDiagValue_A[i_pt];
    tripletList_A.push_back({i_pt, i_pt, e_diag});
    for (int i=eStart; i<eEnd; i++) {
      int j_pt = ePair.second[i];
      e_val = ListOffDiagValue_A[i];
      tripletList_A.push_back({i_pt, j_pt, e_val});
    }
    // The B matrix
    e_diag = ListDiagValue_B[i_pt];
    tripletList_B.push_back({i_pt, i_pt, e_diag});
    for (int i=eStart; i<eEnd; i++) {
      int j_pt = ePair.second[i];
      e_val = ListOffDiagValue_B[i];
      tripletList_B.push_back({i_pt, j_pt, e_val});
    }
    // The SI matrix
    tripletList_B_SI.push_back({i_pt, i_pt, SI(i_pt)});
  }
  MySparseMatrix<double> SpMat_A_SI(nb_point, nb_point);
  MySparseMatrix<double> SpMat_A(nb_point, nb_point);
  MySparseMatrix<double> SpMat_B(nb_point, nb_point);
  MySparseMatrix<double> SpMat_B_SI(nb_point, nb_point);
  SpMat_A_SI.setFromTriplets(tripletList_A_SI.begin(), tripletList_A_SI.end());
  SpMat_A.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
  SpMat_B.setFromTriplets(tripletList_B.begin(), tripletList_B.end());
  SpMat_B_SI.setFromTriplets(tripletList_B_SI.begin(), tripletList_B_SI.end());
  //
  // Now computing eigenvalues.
  //
  std::vector<double> ListEig;
  std::vector<MyVector<double>> ListVect;
  int nb_out = std::min(maxNbEig, nb_point);
  if (UseSI_Bmatrix) {
    std::cerr << "Using the diagonal SI matrix for the B matrix\n";
    Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(SpMat_A_SI);
    MyVector<double> ListEig_A=eig.eigenvalues();
    MyMatrix<double> ListVect_A=eig.eigenvectors();
    for (int i_eig=0; i_eig<nb_out; i_eig++) {
      MyVector<double> Height(nb_point);
      int pos = i_eig;
      for (int ipt=0; ipt<nb_point; ipt++)
        Height(ipt) = ListVect_A(ipt, pos) / Coeff(ipt);
      ListEig.push_back(ListEig_A[i_eig]);
      ListVect.push_back(Height);
    }
  } else {
    std::cerr << "Using the B matrix computed from triangle integration and Generalized eigenvalue\n";
    Spectra::SparseSymMatProd<double> op(SpMat_A);
    Spectra::SparseCholesky<double> Bop(SpMat_B);
    int nev = nb_out;
    Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> geigs(op, Bop, nev, ncv);

    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::SmallestAlge);
    std::cerr << "nconv=" << nconv << "\n";

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;
    auto result = geigs.info();
    if (result == Spectra::CompInfo::Successful) {
      evalues = geigs.eigenvalues();
      evecs = geigs.eigenvectors();
    } else {
      if (result == Spectra::CompInfo::NotComputed)
        std::cerr << "Error NOT_COMPUTED\n";
      if (result == Spectra::CompInfo::NotConverging)
        std::cerr << "Error NOT_CONVERGING\n";
      if (result == Spectra::CompInfo::NumericalIssue)
        std::cerr << "Error NUMERICAL_ISSUE\n";
      std::cerr << "The geigs algorithm failed\n";
      throw TerminalException{1};
    }
    for (int i_eig=0; i_eig<nb_out; i_eig++) {
      int pos = i_eig;
      MyVector<double> Height(nb_point);
      for (int ipt=0; ipt<nb_point; ipt++)
        Height(ipt) = evecs(ipt, pos);
      ListEig.push_back(evalues(i_eig));
      ListVect.push_back(Height);
    }
  }
  double gCst = 9.81;
  double piCst = 3.14159;
  std::vector<PeriodicSolution> ListSol;
  for (int i_eig=0; i_eig<nb_out; i_eig++) {
    // The lowest eigenvalue are the most interesting, so start from the lowest.
    int pos = i_eig;
    // The ListEig is the eigenvalues of the quadratic form.
    // We have the relation lambda_{quad} = |Omega| lambda_{operator}.
    double lambda_quad = ListEig[pos];
    double lambda_oper = lambda_quad;
    // Formula is lambda = omega^2 / g
    double omega = sqrt(lambda_oper * gCst);
    // Formula for period is omega = 2\pi / T
    double period = 2*piCst / omega;
    MyVector<double> Height = ListVect[i_eig];
    double e_max = Height.maxCoeff();
    double e_min = Height.minCoeff();
    double e_max_min = std::max(std::abs(e_max), std::abs(e_min));
    if (RescaleEigenvector)
      Height /= e_max_min;
    double VolSum=0, VolAbsSum=0;
    for (int ipt=0; ipt<nb_point; ipt++) {
      VolSum += SI(ipt) * Height(ipt);
      VolAbsSum += std::abs(SI(ipt) * Height(ipt));
    }
    std::cerr << "i_eig=" << i_eig << " VolSum=" << VolSum << " VolAbsSum=" << VolAbsSum << "\n";
    //
    double Norm_F_SI=0;
    for (int ipt=0; ipt<nb_point; ipt++) {
      Norm_F_SI += SI(ipt) * Height(ipt) * Height(ipt);
    }
    double Norm_F_Int = 0;
    for (int i_elt=0; i_elt<nb_elt; i_elt++) {
      int i0 = GrdArr.INE(i_elt,0);
      int i1 = GrdArr.INE(i_elt,1);
      int i2 = GrdArr.INE(i_elt,2);
      double f0 = Height(i0);
      double f1 = Height(i1);
      double f2 = Height(i2);
      double area = ListArea(i_elt) / 6;
      Norm_F_Int += area * (f0*f0 + f1*f1 + f2*f2 + f0*f1 + f0*f2 + f1*f2);
    }
    std::cerr << "  Norm_F_SI=" << Norm_F_SI << " Norm_F_Int=" << Norm_F_Int << "\n";
    //
    MyVector<double> V1 = SpMat_A * Height;
    MyVector<double> V2 = SpMat_B * Height;
    MyVector<double> V2_SI = SpMat_B_SI * Height;
    double scal_B = V2.dot(Height);
    double scal_B_SI = V2_SI.dot(Height);
    std::cerr << "  RES(V1,V2) : " << GetPairVectorInfo(V1, lambda_oper * V2) << "\n";
    std::cerr << "  RES(V2,V2_SI) : " << GetPairVectorInfo(V2_SI, V2) << "\n";
    std::cerr << "  scal_B=" << scal_B << " scal_B_SI=" << scal_B_SI << "\n";
    //
    PeriodicSolution eSol{lambda_oper, period, Height};
    ListSol.push_back(eSol);
  }
  return ListSol;
}


void WriteSeicheInfoAsNetcdfFile(std::vector<PeriodicSolution> const& ListSol, std::string const& FileName)
{
  int nb_sol = ListSol.size();
  int mnp = ListSol[0].Height.size();
  netCDF::NcFile dataFile(FileName, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  netCDF::NcDim eDimOne=dataFile.addDim("ocean_time", nb_sol);
  netCDF::NcDim eDimTwo=dataFile.addDim("mnp", mnp);
  //
  std::vector<std::string> ListDimTime{"ocean_time"};
  std::vector<std::string> ListDimH{"ocean_time", "mnp"};

  netCDF::NcVar eVarData_ocean=dataFile.addVar("ocean_time", "double", ListDimTime);
  eVarData_ocean.putAtt("units", "seconds since 1858-11-17 00:00:00");
  eVarData_ocean.putAtt("calendar", "gregorian");
  netCDF::NcVar eVarData_period=dataFile.addVar("periods", "double", ListDimTime);
  netCDF::NcVar eVarData_H    =dataFile.addVar("H", "double", ListDimH);

  int nb = nb_sol * mnp;
  std::vector<double> eFieldH(nb);
  std::vector<double> ListTime(nb_sol);
  std::vector<double> ListPeriod(nb_sol);
  for (int i_sol=0; i_sol<nb_sol; i_sol++) {
    ListTime[i_sol] = 86400 * i_sol;
    double e_period = ListSol[i_sol].period;
    ListPeriod[i_sol] = e_period;
    double lambda = ListSol[i_sol].lambda;
    double e_max = ListSol[i_sol].Height.maxCoeff();
    double e_min = ListSol[i_sol].Height.minCoeff();
    std::cerr << "i_sol=" << i_sol << " : lambda=" << lambda << " e_min=" << e_min << " e_max=" << e_max << " period = " << e_period << "\n";
    for (int ipt=0; ipt<mnp; ipt++) {
      int pos = mnp * i_sol + ipt;
      eFieldH[pos] = ListSol[i_sol].Height(ipt);
    }
  }
  eVarData_H.putVar(eFieldH.data());
  eVarData_ocean.putVar(ListTime.data());
  eVarData_period.putVar(ListPeriod.data());
}



#endif
