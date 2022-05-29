// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_ROMSBOUNDARY_H_
#define SRC_OCEAN_ROMSBOUNDARY_H_

#include "NCL_Kernel.h"
#include "Plotting_fct.h"
#include "ROMSfunctionality.h"
#include <map>
#include <utility>
#include <vector>
#include <string>

FullNamelist NAMELIST_GetStandardPLOT_BOUNDARY() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["BEGTC"] = "20110915.000000";
  ListStringValues1["ENDTC"] = "20110925.000000";
  ListDoubleValues1["DELTC"] = 600;
  ListStringValues1["UNITC"] = "SEC";
  // KindSelect, possible values: direct, monthly, seasonal, yearly, specific
  ListStringValues1["KindSelect"] = "direct";
  ListStringValues1["GridFile"] = "UNK";
  ListStringValues1["BoundaryFile"] = "UNK";
  ListStringValues1["PicPrefix"] = "UNK";
  ListStringValues1["Extension"] = "png";
  ListStringValues1["__NaturePlot"] = "boundary";
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListIntValues1["NPROC"] = 1;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Lines_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues = ListIntValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  BlockPROC.ListDoubleValues = ListDoubleValues1;
  BlockPROC.ListListDoubleValues = ListListDoubleValues1;
  BlockPROC.ListListIntValues = ListListIntValues1;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListListStringValues = ListListStringValues1;
  ListBlock["PROC"] = BlockPROC;
  // PLOT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListListStringValues2["ListSides"] = {};
  ListBoolValues2["VariableRange"] = false;
  ListBoolValues2["PlotTemp"] = false;
  ListBoolValues2["PlotSalt"] = false;
  ListBoolValues2["PlotU"] = false;
  ListBoolValues2["PlotV"] = false;
  ListIntValues2["nbLevelSpa"] = 30;
  ListIntValues2["nbLabelStride"] = 10;
  ListBoolValues2["DrawAnnotation"] = false;
  ListDoubleValues2["AnnotationLon"] = 0;
  ListDoubleValues2["AnnotationLat"] = 0;
  ListStringValues2["AnnotationText"] = "something to write";
  ListStringValues2["FileDirectNCLins"] = "irrelevant";
  ListBoolValues2["DoTitle"] = true;
  ListDoubleValues2["vcRefLengthF"] = 0.02;
  ListBoolValues2["DoColorBar"] = true;
  ListStringValues2["cnFillMode"] = "RasterFill";
  ListBoolValues2["cnFillOn"] = true;
  ListBoolValues2["cnLinesOn"] = false;
  ListBoolValues2["cnLineLabelsOn"] = false;
  ListBoolValues2["cnSmoothingOn"] = true;
  ListStringValues2["ColorMap"] = "BlAqGrYeOrReVi200";
  ListBoolValues2["PrintMMA"] = false;
  ListBoolValues2["DoTitleString"] = true;
  ListStringValues2["LandPortr"] = "Landscape";
  ListStringValues2["optStatStr"] = "double";
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues = ListIntValues2;
  BlockPLOT.ListBoolValues = ListBoolValues2;
  BlockPLOT.ListDoubleValues = ListDoubleValues2;
  BlockPLOT.ListListDoubleValues = ListListDoubleValues2;
  BlockPLOT.ListStringValues = ListStringValues2;
  BlockPLOT.ListListStringValues = ListListStringValues2;
  ListBlock["PLOT"] = BlockPLOT;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

template <typename T>
MyVector<T> GetMatrixSide(MyMatrix<T> const &M, std::string const &eSide) {
  int eta = M.rows();
  int xi = M.cols();
  if (eSide == "South") {
    MyVector<T> V(xi);
    for (int i = 0; i < xi; i++)
      V(i) = M(0, i);
    return V;
  }
  if (eSide == "North") {
    MyVector<T> V(xi);
    for (int i = 0; i < xi; i++)
      V(i) = M(eta - 1, i);
    return V;
  }
  if (eSide == "East") {
    MyVector<T> V(eta);
    for (int i = 0; i < eta; i++)
      V(i) = M(i, xi - 1);
    return V;
  }
  if (eSide == "West") {
    MyVector<T> V(eta);
    for (int i = 0; i < eta; i++)
      V(i) = M(i, 0);
    return V;
  }
  std::cerr << "Failed to find Matching entry in GetMatrixSide\n";
  throw TerminalException{1};
}

void BOUND_Plotting_Function(FullNamelist const &eFull) {
  struct ArrSide {
    std::string InputName;
    std::string NcName;
    int eta_ts;
    int eta_u;
    int eta_v;
    MyVector<double> DEP_rho;
    MyVector<double> DEP_u;
    MyVector<double> DEP_v;
    MyVector<uint8_t> MSK_rho;
    MyVector<uint8_t> MSK_u;
    MyVector<uint8_t> MSK_v;
  };
  struct TypeVar {
    std::string VarName;
    std::string SystemName;
    std::string Nature;
  };
  auto GetVectorDEP_MSK = [&](std::string const &typeName, ArrSide const &eSide)
      -> std::pair<MyVector<double>, MyVector<uint8_t>> {
    if (typeName == "rho")
      return {eSide.DEP_rho, eSide.MSK_rho};
    if (typeName == "u")
      return {eSide.DEP_u, eSide.MSK_u};
    if (typeName == "v")
      return {eSide.DEP_v, eSide.MSK_v};
    std::cerr << "Failed to find Matching entry in GetVectorDEP_MSK\n";
    throw TerminalException{1};
  };
  std::string strSRho = "s_rho";
  //
  // PROC entries
  //
  SingleBlock eBlPROC = eFull.ListBlock.at("PROC");
  std::string strBEGTC = eBlPROC.ListStringValues.at("BEGTC");
  std::string strENDTC = eBlPROC.ListStringValues.at("ENDTC");
  std::string BoundaryFile = eBlPROC.ListStringValues.at("BoundaryFile");
  netCDF::NcFile dataFile(BoundaryFile, netCDF::NcFile::read);
  int s_rho = NC_ReadDimension(dataFile, strSRho);
  std::vector<double> ListTime = NC_ReadTimeFromFile(BoundaryFile, "zeta_time");
  std::pair<int, int> PairFirstLast =
      GetIdx_first_last(ListTime, strBEGTC, strENDTC);
  int idx_first = PairFirstLast.first;
  int idx_last = PairFirstLast.second;
  int idx_len = idx_last - idx_first;
  //  bool
  //  WriteITimeInFileName=eBlPROC.ListBoolValues.at("WriteITimeInFileName");
  std::string GridFile = eBlPROC.ListStringValues.at("GridFile");
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  ARVDtyp ARVD = ReadROMSverticalStratification(BoundaryFile);
  std::cerr << "PROC entries read\n";
  //
  // PLOT parameter
  //
  SingleBlock eBlPLOT = eFull.ListBlock.at("PLOT");
  bool PlotTemp = eBlPLOT.ListBoolValues.at("PlotTemp");
  bool PlotSalt = eBlPLOT.ListBoolValues.at("PlotSalt");
  bool PlotU = eBlPLOT.ListBoolValues.at("PlotU");
  bool PlotV = eBlPLOT.ListBoolValues.at("PlotV");
  bool VariableRange = eBlPLOT.ListBoolValues.at("VariableRange");
  std::vector<std::string> ListSides =
      eBlPLOT.ListListStringValues.at("ListSides");
  std::vector<std::string> ListSidesTot = {"South", "North", "West", "East"};
  std::vector<ArrSide> ListArrSide;
  for (auto &eStr : ListSides) {
    bool IsPresent = false;
    for (auto &eStrTot : ListSidesTot)
      if (eStr == eStrTot)
        IsPresent = true;
    if (!IsPresent) {
      std::cerr << "The variable eStr = " << eStr << " is not allowed\n";
      std::cerr << "Allowed variables are South, North, West, East\n";
      throw TerminalException{1};
    }
    ArrSide eArrSide;
    eArrSide.InputName = eStr;
    eArrSide.NcName = UpperCaseToLowerCase(eStr);
    eArrSide.DEP_rho = GetMatrixSide(GetDEP(GrdArr.GrdArrRho), eStr);
    eArrSide.DEP_u = GetMatrixSide(GetDEP(GrdArr.GrdArrU), eStr);
    eArrSide.DEP_v = GetMatrixSide(GetDEP(GrdArr.GrdArrV), eStr);
    eArrSide.MSK_rho = GetMatrixSide(GrdArr.GrdArrRho.MSK, eStr);
    eArrSide.MSK_u = GetMatrixSide(GrdArr.GrdArrU.MSK, eStr);
    eArrSide.MSK_v = GetMatrixSide(GrdArr.GrdArrV.MSK, eStr);
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

  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  std::cerr << "ePerm obtained\n";
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  std::cerr << "eCall obtained\n";

  for (int idx = 0; idx < idx_len; idx++) {
    int iTime = idx_first + idx;
    double eTimeDay = ListTime[iTime];
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    std::string strFile = DATE_ConvertMjd2mystringFile(eTimeDay);
    std::cerr << "idx=" << idx << "/" << idx_len << " iTime=" << iTime
              << " date=" << strPres << "\n";
    for (auto &eArrSide : ListArrSide) {
      for (auto &eTypeVar : ListTypeVar) {
        std::string varName = eTypeVar.VarName + "_" + eArrSide.NcName;
        netCDF::NcVar data = dataFile.getVar(varName);
        //
        // Reading the bathymetry
        //
        std::pair<MyVector<double>, MyVector<uint8_t>> DEP_MSK =
            GetVectorDEP_MSK(eTypeVar.Nature, eArrSide);
        MyVector<double> DEP = DEP_MSK.first;
        MyVector<uint8_t> MSK_grid = DEP_MSK.second;
        int siz = DEP.size();
        double maxDep = DEP.maxCoeff();
        int NbVert = 100;
        double DeltaZ = maxDep / static_cast<double>(NbVert);
        MyVector<double> ListVertPos(NbVert + 1);
        for (int i = 0; i <= NbVert; i++) {
          double eVertPos = -maxDep + i * DeltaZ;
          ListVertPos(i) = eVertPos;
        }
        //
        // Reading the data sets
        //
        std::vector<size_t> start{size_t(iTime), 0, 0};
        std::vector<size_t> count{1, size_t(s_rho), size_t(siz)};
        MyVector<double> eVal =
            NC_ReadVariable_data_start_count(data, start, count);
        MyMatrix<double> M(s_rho, siz);
        int idx = 0;
        for (int i = 0; i < s_rho; i++)
          for (int j = 0; j < siz; j++) {
            M(i, j) = eVal[idx];
            idx++;
          }
        //        std::cerr << "M assigned\n";
        //
        // The vertical coordinate
        //
        MyMatrix<double> F(NbVert + 1, siz);
        MyMatrix<uint8_t> MSK(NbVert + 1, siz);
        for (int i = 0; i < siz; i++) {
          double eZeta = 0;
          double eDep = DEP(i);
          MyVector<double> Zr_out = GetVertCoord_R(ARVD, eDep, eZeta);
          double eps = 0.00001;
          for (int iV = 0; iV <= NbVert; iV++) {
            double eVert = ListVertPos(iV);
            int eMSK = 0;
            double eF = 0;
            if (eVert >= Zr_out(0) && MSK_grid(i) == 1) {
              eMSK = 1;
              if (eVert > Zr_out(s_rho - 1) - eps) {
                eF = M(s_rho - 1, i);
                //                std::cerr << "i=" << i << " iV=" << iV << "
                //                eF=" << eF << "\n";
              } else {
                bool IsMatch = false;
                for (int iS = 0; iS < s_rho - 1; iS++) {
                  double dep1 = Zr_out(iS);
                  double dep2 = Zr_out(iS + 1);
                  if (dep1 - eps <= eVert && eVert <= dep2 + eps) {
                    IsMatch = true;
                    double alpha1 = (dep2 - eVert) / (dep2 - dep1);
                    double alpha2 = (eVert - dep1) / (dep2 - dep1);
                    eF = M(iS, i) * alpha1 + M(iS + 1, i) * alpha2;
                    //                    std::cerr << "i=" << i << " iV=" << iV
                    //                    << " iS=" << iS << " eF=" << eF <<
                    //                    "\n"; std::cerr << "   alpha1=" <<
                    //                    alpha1 << " alpha2=" << alpha2 << "
                    //                    M1=" << M(iS,i) << " M2=" << M(iS+1,i)
                    //                    << "\n";
                  }
                }
                if (!IsMatch) {
                  std::cerr << "Failed to find matching depth\n";
                  throw TerminalException{1};
                }
              }
            }
            MSK(iV, i) = eMSK;
            F(iV, i) = eF;
          }
        }
        std::cerr << " F(min/max)=" << F.minCoeff() << " / " << F.maxCoeff()
                  << "\n";
        //        std::cerr << "We have F and MSK\n";
        //
        // The eta/xi coordinate
        //
        MyMatrix<double> LON(NbVert + 1, siz);
        MyMatrix<double> LAT(NbVert + 1, siz);
        for (int i = 0; i < siz; i++)
          for (int iV = 0; iV <= NbVert; iV++) {
            LON(iV, i) = static_cast<double>(i);
            LAT(iV, i) = ListVertPos(iV);
          }
        //        std::cerr << "We have LON/LAT\n";
        //
        // The grid array
        //
        GridArray GrdArr;
        GrdArr.IsFE = 0;
        GrdArr.IsSpherical = false;
        GrdArr.GrdArrRho.LON = LON;
        GrdArr.GrdArrRho.LAT = LAT;
        GrdArr.GrdArrRho.MSK = MSK;
        //        std::cerr << "We have GrdArr\n";
        //
        // The plotting function
        //
        std::cerr << "  eName=" << eTypeVar.SystemName << "\n";
        RecVar eRecVarTriv = RetrieveTrivialRecVar(eTypeVar.SystemName);
        RecVar NewRecVar;
        NewRecVar.RecS = eRecVarTriv.RecS;
        NewRecVar.RecS.eTimeDay = eTimeDay;
        NewRecVar.RecS.strPres = strPres;
        NewRecVar.RecS.strFile = strFile;
        NewRecVar.RecS.iTime = iTime;
        if (VariableRange) {
          PairMinMax ePair = ComputeMinMaxMask(MSK, F);
          NewRecVar.RecS.mindiff = ePair.TheMin;
          NewRecVar.RecS.maxdiff = ePair.TheMax;
          NewRecVar.RecS.minval = ePair.TheMin;
          NewRecVar.RecS.maxval = ePair.TheMax;
        }
        std::cerr << "    mindiff=" << NewRecVar.RecS.mindiff
                  << " maxdiff=" << NewRecVar.RecS.maxdiff << "\n";
        std::cerr << "    minval=" << NewRecVar.RecS.minval
                  << " maxval=" << NewRecVar.RecS.maxval << "\n";
        NewRecVar.F = F;
        //        std::cerr << "We have NewRecVar\n";
        //
        DrawArr eDrawArr = ePerm.eDrawArr;
        eDrawArr.DrawContourBathy = false;
        eDrawArr.DoTitle = true;
        eDrawArr.TitleStr = varName + " " + NewRecVar.RecS.strPres;
        eDrawArr.nbLevelSpa = nbLevelSpa;
        //        std::cerr << "We have eDrawArr\n";
        //
        eDrawArr.VarNameUF = varName;
        std::string FileName = ePerm.eDir + varName + "_" +
                               NewRecVar.RecS.strAll + "_" +
                               StringNumber(iTime, 4);
        //
        PLOT_PCOLOR(FileName, GrdArr, eDrawArr, NewRecVar, eCall, ePerm);
      }
    }
  }
}

// clang-format off
#endif  // SRC_OCEAN_ROMSBOUNDARY_H_
// clang-format on
