#ifndef DEFINE_SEICHE_NUMERIC_EIGENVALUE
#include "MAT_Matrix.h"
#include "Model_grids.h"
#include "Namelist.h"


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
  ListBoolValues1["RescaleEigenvector"] = true;
  ListBoolValues1["ApplySI"] = false;
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
std::vector<PeriodicSolution> ComputeEigenvaluesSWE1(double const& h0, int const& maxNbEig, GridArray const& GrdArr, bool const& RescaleEigenvector, bool const& ApplySI)
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
  if (ApplySI) {
    for (int i_elt=0; i_elt<nb_elt; i_elt++) {
      for (int idx=0; idx<3; idx++) {
        int ip = GrdArr.INE(i_elt, idx);
        SI(ip) += ListArea(i_elt);
      }
    }
    for (int i_pt=0; i_pt<nb_point; i_pt++)
      SI(i_pt) = sqrt(SI(i_pt));
  }
  else {
    for (int i_pt=0; i_pt<nb_point; i_pt++)
      SI(i_pt) = 1;
  }
  GraphSparseImmutable eG = GetUnstructuredVertexAdjInfo(GrdArr.INE, nb_point);
  std::pair<std::vector<int>, std::vector<int>> ePair = eG.Get_ListStart_ListListAdj();
  int nb_adj = ePair.second.size();
  //
  // Now building the matrix
  //
  std::vector<double> ListDiagValue(nb_point,0);
  std::vector<double> ListOffDiagValue(nb_adj, 0);
  double minWaterHeight = std::numeric_limits<double>::max();
  double maxWaterHeight = std::numeric_limits<double>::min();
  double sumDelta=0;
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
    //
    double tmp;
    double alpha = avg_val * ListArea(i_elt) / (delta * delta);
    // dimensionality of alpha is [L]
    ListDiagValue[i0] += (b0 * b0 + c0 * c0) * alpha;
    ListDiagValue[i1] += (b1 * b1 + c1 * c1) * alpha;
    ListDiagValue[i2] += (b2 * b2 + c2 * c2) * alpha;
    tmp = (b0 * b1 + c0 * c1) * alpha;
    ListOffDiagValue[idx01] += tmp;
    ListOffDiagValue[idx10] += tmp;
    tmp = (b0 * b2 + c0 * c2) * alpha;
    ListOffDiagValue[idx02] += tmp;
    ListOffDiagValue[idx20] += tmp;
    tmp = (b1 * b2 + c1 * c2) * alpha;
    ListOffDiagValue[idx12] += tmp;
    ListOffDiagValue[idx21] += tmp;
  }
  std::cerr << "minWaterHeight = " << minWaterHeight << "  maxWaterHeight = " << maxWaterHeight << "\n";
  std::cerr << "sumDelta=" << sumDelta << "\n";
  int tot_dim = nb_point;
  MySparseMatrix<double> eSP(tot_dim, tot_dim);
  using Ttrip = Eigen::Triplet<double>;
  std::vector<Ttrip> tripletList;
  for (int i_pt=0; i_pt<nb_point; i_pt++) {
    double e_diag = ListDiagValue[i_pt];
    double e_si = SI(i_pt);
    double f_diag = e_diag / (e_si * e_si);
    tripletList.push_back({i_pt, i_pt, f_diag});
    int eStart = ePair.first[i_pt];
    int eEnd = ePair.first[i_pt+1];
    for (int i=eStart; i<eEnd; i++) {
      int j_pt = ePair.second[i];
      double f_si = SI(j_pt);
      double e_val = ListOffDiagValue[i];
      double f_val = e_val / (e_si * f_si);
      tripletList.push_back({i_pt, j_pt, f_val});
    }
  }
  MySparseMatrix<double> SpMat=MySparseMatrix<double>(nb_point, nb_point);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  //
  // Now computing eigenvalues.
  //
  Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(SpMat);
  MyVector<double> ListEig=eig.eigenvalues();
  MyMatrix<double> ListVect=eig.eigenvectors();
  //  for (int i_eig=0; i_eig<tot_dim; i_eig++) {
  //    std::cerr << "i_eig=" << i_eig << " lambda=" << ListEig(i_eig) << "\n";
  //  }
  double gCst = 9.81;
  double piCst = 3.14159;
  std::vector<PeriodicSolution> ListSol;
  int nb_out = std::min(maxNbEig, nb_point);
  for (int i_eig=0; i_eig<nb_out; i_eig++) {
    // The lowest eigenvalue are the most interesting, so start from the lowest.
    int pos = i_eig;
    // The ListEig is the eigenvalues of the quadratic form.
    // We have the relation lambda_{quad} = |Omega| lambda_{operator}.
    double lambda_quad = ListEig(pos);
    double lambda_oper = lambda_quad;
    // Formula is lambda = omega^2 / g
    double omega = sqrt(lambda_oper * gCst);
    // Formula for period is omega = 2\pi / T
    double period = 2*piCst / omega;
    MyVector<double> Height(nb_point);
    for (int ipt=0; ipt<nb_point; ipt++)
      Height(ipt) = ListVect(ipt, pos);
    double e_max = Height.maxCoeff();
    double e_min = Height.minCoeff();
    double e_max_min = std::max(std::abs(e_max), std::abs(e_min));
    if (RescaleEigenvector)
      Height /= e_max_min;
    //
    double sumQuadF = 0;
    double sumdx2 = 0, sumdy2=0;
    double maxdx = 0, maxdy=0;
    double sumSqrF=0;
    for (int i_elt=0; i_elt<nb_elt; i_elt++) {
      int i0 = GrdArr.INE(i_elt,0);
      int i1 = GrdArr.INE(i_elt,1);
      int i2 = GrdArr.INE(i_elt,2);
      double f0 = Height(i0);
      double f1 = Height(i1);
      double f2 = Height(i2);
      double avg_f = (f0 + f1 + f2) / 3.0;
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
      double delta_x = (b0 * f0 + b1 * f1 + b2 * f2) / delta;
      double delta_y = (c0 * f0 + c1 * f1 + c2 * f2) / delta;
      sumdx2 += delta_x*delta_x;
      sumdy2 += delta_y*delta_y;
      maxdx = std::max(maxdx, std::abs(delta_x));
      maxdy = std::max(maxdy, std::abs(delta_y));
      sumQuadF += avg_val * ListArea(i_elt) * (delta_x*delta_x + delta_y*delta_y);
      sumSqrF += ListArea(i_elt) * avg_f * avg_f;
    }
    double sumF=0;
    for (int ipt=0; ipt<nb_point; ipt++) {
      sumF += Height(ipt) * Height(ipt);
    }
    double rayleigh_quot = sumQuadF / sumF;
    std::cerr << "i_eig=" << i_eig << " rayleigh_quot=" << rayleigh_quot << " lambda_quad=" << lambda_quad << "\n";
    std::cerr << "                      sumdx2=" << sumdx2 << " sumdy2=" << sumdy2 << "\n";
    std::cerr << "                       maxdx=" << maxdx  << "  maxdy=" << maxdy  << "\n";
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
    double e_sum = ListSol[i_sol].Height.sum();
    std::cerr << "i_sol=" << i_sol << " : lambda=" << lambda << " e_min=" << e_min << " e_max=" << e_max << " e_sum=" << e_sum << " period = " << e_period << "\n";
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
