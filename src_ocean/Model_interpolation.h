// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_MODEL_INTERPOLATION_H_
#define SRC_OCEAN_MODEL_INTERPOLATION_H_

#include "Basic_netcdf.h"
#include "CommonFuncModel.h"
#include "Interpolation.h"
#include "Model_grids.h"
#include "NamelistExampleOcean.h"
#include "Statistics.h"
#include "Triangulations.h"
#include "WW3_includes.h"
#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

FullNamelist NAMELIST_GetStandardMODEL_MERGING() {
  std::map<std::string, SingleBlock> ListBlock;
  // INPUT
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListListStringValues1["ListMODELNAME"] = {"UNK"};
  ListListStringValues1["ListGridFile"] = {"UNK"};
  ListListStringValues1["ListHisPrefix"] = {"UNK"};
  ListListIntValues1["ListSpongeSize"] = {1};
  ListListIntValues1["ListFatherGrid"] = {-1};
  ListBoolValues1["DoClimatology"] = false;
  ListBoolValues1["AllowExtrapolation"] = false;
  ListBoolValues1["PrintMMA"] = false;
  SingleBlock BlockINPUT;
  BlockINPUT.ListIntValues = ListIntValues1;
  BlockINPUT.ListBoolValues = ListBoolValues1;
  BlockINPUT.ListDoubleValues = ListDoubleValues1;
  BlockINPUT.ListListDoubleValues = ListListDoubleValues1;
  BlockINPUT.ListListIntValues = ListListIntValues1;
  BlockINPUT.ListStringValues = ListStringValues1;
  BlockINPUT.ListListStringValues = ListListStringValues1;
  ListBlock["INPUT"] = BlockINPUT;
  // OUTPUT
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListStringValues2["MODELNAME"] = "unset MODELNAME";
  ListDoubleValues2["MinLat"] = -1;
  ListDoubleValues2["MaxLat"] = -1;
  ListDoubleValues2["MinLon"] = -1;
  ListDoubleValues2["MaxLon"] = -1;
  ListDoubleValues2["deltaKM"] = -1;
  ListStringValues2["GridFile"] = "unset GridFile";
  ListStringValues2["BoundFile"] = "unset";
  ListStringValues2["HisPrefix"] = "unset HisPrefix";
  ListStringValues2["BEGTC"] = "20110915.000000";
  ListStringValues2["ENDTC"] = "20110925.000000";
  ListDoubleValues2["DELTC"] = 600;
  ListStringValues2["UNITC"] = "SEC";
  ListDoubleValues2["DEFINETC"] = 86400;
  ListStringValues2["KindSelect"] = "direct";
  ListDoubleValues2["TimeFrameDay"] = 1;
  ListStringValues2["WaveWatchFormat"]="UNFORMATTED or FORMATTED or PRNC";
  ListBoolValues2["DoNetcdfWrite"] = false;
  ListBoolValues2["DoGribWrite"] = false;
  ListBoolValues2["DoRomsWrite_Surface"] = false;
  ListBoolValues2["DoRomsWrite_InitialHistory"] = false;
  ListBoolValues2["DoRomsWrite_Boundary"] = false;
  ListBoolValues2["DoWaveWatchWrite"] = false;
  ListBoolValues2["DoSfluxWrite"] = false;
  SingleBlock BlockOUTPUT;
  BlockOUTPUT.ListIntValues = ListIntValues2;
  BlockOUTPUT.ListBoolValues = ListBoolValues2;
  BlockOUTPUT.ListDoubleValues = ListDoubleValues2;
  BlockOUTPUT.ListListDoubleValues = ListListDoubleValues2;
  BlockOUTPUT.ListStringValues = ListStringValues2;
  BlockOUTPUT.ListListStringValues = ListListStringValues2;
  ListBlock["OUTPUT"] = BlockOUTPUT;
  // ANALYTIC
  std::map<std::string, int> ListIntValues23;
  std::map<std::string, bool> ListBoolValues23;
  std::map<std::string, double> ListDoubleValues23;
  std::map<std::string, std::vector<double>> ListListDoubleValues23;
  std::map<std::string, std::string> ListStringValues23;
  std::map<std::string, std::vector<std::string>> ListListStringValues23;
  ListListStringValues23["AnalyticalListNameVariables"] = {};
  ListListDoubleValues23["AnalyticalListConstantValuesRho"] = {};
  ListListDoubleValues23["AnalyticalListConstantValuesU"] = {};
  ListListDoubleValues23["AnalyticalListConstantValuesV"] = {};
  SingleBlock BlockANALYTIC;
  BlockANALYTIC.ListIntValues = ListIntValues23;
  BlockANALYTIC.ListBoolValues = ListBoolValues23;
  BlockANALYTIC.ListDoubleValues = ListDoubleValues23;
  BlockANALYTIC.ListListDoubleValues = ListListDoubleValues23;
  BlockANALYTIC.ListStringValues = ListStringValues23;
  BlockANALYTIC.ListListStringValues = ListListStringValues23;
  ListBlock["ANALYTIC"] = BlockANALYTIC;
  // MEASUREMENT
  std::map<std::string, int> ListIntValues24;
  std::map<std::string, bool> ListBoolValues24;
  std::map<std::string, double> ListDoubleValues24;
  std::map<std::string, std::vector<double>> ListListDoubleValues24;
  std::map<std::string, std::string> ListStringValues24;
  std::map<std::string, std::vector<std::string>> ListListStringValues24;
  ListListStringValues24["ListFilesMeasurement"] = {};
  SingleBlock BlockMEAS;
  BlockMEAS.ListIntValues = ListIntValues24;
  BlockMEAS.ListBoolValues = ListBoolValues24;
  BlockMEAS.ListDoubleValues = ListDoubleValues24;
  BlockMEAS.ListListDoubleValues = ListListDoubleValues24;
  BlockMEAS.ListStringValues = ListStringValues24;
  BlockMEAS.ListListStringValues = ListListStringValues24;
  ListBlock["MEASUREMENT"] = BlockMEAS;
  // ROMS_SURFACE
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::vector<double>> ListListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListBoolValues3["IsRegrid"] = false;
  ListBoolValues3["SingleFile"] = true;
  ListStringValues3["RomsFile_surf"] = "unset";
  SingleBlock BlockROMS_SURFACE;
  BlockROMS_SURFACE.ListIntValues = ListIntValues3;
  BlockROMS_SURFACE.ListBoolValues = ListBoolValues3;
  BlockROMS_SURFACE.ListDoubleValues = ListDoubleValues3;
  BlockROMS_SURFACE.ListListDoubleValues = ListListDoubleValues3;
  BlockROMS_SURFACE.ListStringValues = ListStringValues3;
  BlockROMS_SURFACE.ListListStringValues = ListListStringValues3;
  ListBlock["ROMS_SURFACE"] = BlockROMS_SURFACE;
  // ROMS_INITIAL
  std::map<std::string, int> ListIntValues4;
  std::map<std::string, bool> ListBoolValues4;
  std::map<std::string, double> ListDoubleValues4;
  std::map<std::string, std::vector<double>> ListListDoubleValues4;
  std::map<std::string, std::string> ListStringValues4;
  std::map<std::string, std::vector<std::string>> ListListStringValues4;
  ListStringValues4["RomsFile_InitialHistory"] = "unset";
  ListIntValues4["ARVD_N"] = -1;
  ListIntValues4["ARVD_Vtransform"] = -1;
  ListIntValues4["ARVD_Vstretching"] = -1;
  ListDoubleValues4["ARVD_Tcline"] = -1;
  ListDoubleValues4["ARVD_hc"] = -1;
  ListDoubleValues4["ARVD_theta_s"] = -1;
  ListDoubleValues4["ARVD_theta_b"] = -1;
  SingleBlock BlockROMS_INIT_HIS;
  BlockROMS_INIT_HIS.ListIntValues = ListIntValues4;
  BlockROMS_INIT_HIS.ListBoolValues = ListBoolValues4;
  BlockROMS_INIT_HIS.ListDoubleValues = ListDoubleValues4;
  BlockROMS_INIT_HIS.ListListDoubleValues = ListListDoubleValues4;
  BlockROMS_INIT_HIS.ListStringValues = ListStringValues4;
  BlockROMS_INIT_HIS.ListListStringValues = ListListStringValues4;
  ListBlock["ROMS_INITIAL_HISTORY"] = BlockROMS_INIT_HIS;
  // ROMS_BOUND
  std::map<std::string, int> ListIntValues5;
  std::map<std::string, bool> ListBoolValues5;
  std::map<std::string, double> ListDoubleValues5;
  std::map<std::string, std::vector<double>> ListListDoubleValues5;
  std::map<std::string, std::string> ListStringValues5;
  std::map<std::string, std::vector<std::string>> ListListStringValues5;
  //  ListStringValues5["MODELNAME"]="unset MODELNAME";
  //  ListStringValues5["GridFile"]="unset GridFile";
  //  ListStringValues5["HisPrefix"]="unset HisPrefix";
  ListStringValues5["RomsFile_bound"] = "unset";
  ListListStringValues5["ListSides"] = {};
  ListIntValues5["ARVD_N"] = -1;
  ListIntValues5["ARVD_Vtransform"] = -1;
  ListIntValues5["ARVD_Vstretching"] = -1;
  ListDoubleValues5["ARVD_Tcline"] = -1;
  ListDoubleValues5["ARVD_hc"] = -1;
  ListDoubleValues5["ARVD_theta_s"] = -1;
  ListDoubleValues5["ARVD_theta_b"] = -1;
  SingleBlock BlockROMS_BOUND;
  BlockROMS_BOUND.ListIntValues = ListIntValues5;
  BlockROMS_BOUND.ListBoolValues = ListBoolValues5;
  BlockROMS_BOUND.ListDoubleValues = ListDoubleValues5;
  BlockROMS_BOUND.ListListDoubleValues = ListListDoubleValues5;
  BlockROMS_BOUND.ListStringValues = ListStringValues5;
  BlockROMS_BOUND.ListListStringValues = ListListStringValues5;
  ListBlock["ROMS_BOUND"] = BlockROMS_BOUND;
  // NETCDF STANDARD
  std::map<std::string, int> ListIntValues6;
  std::map<std::string, bool> ListBoolValues6;
  std::map<std::string, double> ListDoubleValues6;
  std::map<std::string, std::vector<double>> ListListDoubleValues6;
  std::map<std::string, std::string> ListStringValues6;
  std::map<std::string, std::vector<std::string>> ListListStringValues6;
  ListBoolValues6["WriteIFile"] = true;
  ListBoolValues6["WriteDate"] = false;
  ListStringValues6["HisPrefixOut"] = "FinalTarget_";
  SingleBlock BlockNETCDF_STANDARD;
  BlockNETCDF_STANDARD.ListIntValues = ListIntValues6;
  BlockNETCDF_STANDARD.ListBoolValues = ListBoolValues6;
  BlockNETCDF_STANDARD.ListDoubleValues = ListDoubleValues6;
  BlockNETCDF_STANDARD.ListListDoubleValues = ListListDoubleValues6;
  BlockNETCDF_STANDARD.ListStringValues = ListStringValues6;
  BlockNETCDF_STANDARD.ListListStringValues = ListListStringValues6;
  ListBlock["NETCDF_STANDARD"] = BlockNETCDF_STANDARD;
  // GRIB STANDARD
  std::map<std::string, int> ListIntValues7;
  std::map<std::string, bool> ListBoolValues7;
  std::map<std::string, double> ListDoubleValues7;
  std::map<std::string, std::vector<double>> ListListDoubleValues7;
  std::map<std::string, std::string> ListStringValues7;
  std::map<std::string, std::vector<std::string>> ListListStringValues7;
  ListBoolValues7["WriteFromStart"] = true;
  ListStringValues7["HisPrefixOut"] = "FinalTarget_";
  SingleBlock BlockGRIB_STANDARD;
  BlockGRIB_STANDARD.ListIntValues = ListIntValues7;
  BlockGRIB_STANDARD.ListBoolValues = ListBoolValues7;
  BlockGRIB_STANDARD.ListDoubleValues = ListDoubleValues7;
  BlockGRIB_STANDARD.ListListDoubleValues = ListListDoubleValues7;
  BlockGRIB_STANDARD.ListStringValues = ListStringValues7;
  BlockGRIB_STANDARD.ListListStringValues = ListListStringValues7;
  ListBlock["GRIB_STANDARD"] = BlockGRIB_STANDARD;
  // VARS
  std::map<std::string, bool> ListBoolValues100;
  for (auto &eVarName : GetAllPossibleVariables())
    ListBoolValues100[eVarName] = false;
  SingleBlock BlockVARS;
  BlockVARS.ListBoolValues = ListBoolValues100;
  ListBlock["VARS"] = BlockVARS;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

FullNamelist NAMELIST_GetStandard_CREATE_sflux() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["MODELNAME"] = "unset MODELNAME";
  ListStringValues1["GridFile"] = "unset GridFile";
  ListStringValues1["BoundFile"] = "unset";
  ListStringValues1["HisPrefix"] = "unset HisPrefix";
  ListStringValues1["BEGTC"] = "20110915.000000";
  ListDoubleValues1["DELTC"] = 600;
  ListDoubleValues1["TimeFrameDay"] = 1;
  ListStringValues1["UNITC"] = "SEC";
  ListStringValues1["ENDTC"] = "20110925.000000";
  ListStringValues1["OutPrefix"] = "Pictures/DIR_plot/";
  ListIntValues1["idxForc"] = 1;
  ListBoolValues1["DoAIRfiles"] = false;
  ListBoolValues1["DoPRCfiles"] = false;
  ListBoolValues1["DoRADfiles"] = false;
  ListBoolValues1["AnalyticWind"] = false;
  ListBoolValues1["AnalyticPRMSL"] = true;
  ListBoolValues1["AnalyticSPFH"] = true;
  ListBoolValues1["AnalyticSTMP"] = true;
  ListBoolValues1["AnalyticPRATE"] = true;
  ListBoolValues1["AnalyticDLWRF"] = true;
  ListBoolValues1["AnalyticDSWRF"] = true;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues = ListIntValues1;
  BlockPROC.ListBoolValues = ListBoolValues1;
  BlockPROC.ListDoubleValues = ListDoubleValues1;
  BlockPROC.ListStringValues = ListStringValues1;
  BlockPROC.ListListStringValues = ListListStringValues1;
  ListBlock["PROC"] = BlockPROC;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

void CREATE_sflux_files(FullNamelist const &eFull) {
  std::map<std::string, SingleBlock> ListBlock = eFull.ListBlock;
  SingleBlock eBlPROC = eFull.ListBlock.at("PROC");
  std::string eModelName = eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile = eBlPROC.ListStringValues.at("GridFile");
  std::string BoundFile = eBlPROC.ListStringValues.at("BoundFile");
  std::string HisPrefix = eBlPROC.ListStringValues.at("HisPrefix");
  TripleModelDesc eTriple{eModelName, GridFile, BoundFile, HisPrefix, {}};
  bool DoAIRfiles = eBlPROC.ListBoolValues.at("DoAIRfiles");
  bool DoPRCfiles = eBlPROC.ListBoolValues.at("DoPRCfiles");
  bool DoRADfiles = eBlPROC.ListBoolValues.at("DoRADfiles");
  bool AnalyticWind = eBlPROC.ListBoolValues.at("AnalyticWind");
  bool AnalyticPrmsl = eBlPROC.ListBoolValues.at("AnalyticPRMSL");
  bool AnalyticStmp = eBlPROC.ListBoolValues.at("AnalyticSTMP");
  bool AnalyticSpfh = eBlPROC.ListBoolValues.at("AnalyticSPFH");
  bool AnalyticDlwrf = eBlPROC.ListBoolValues.at("AnalyticDLWRF");
  bool AnalyticDswrf = eBlPROC.ListBoolValues.at("AnalyticDSWRF");
  bool AnalyticPrate = eBlPROC.ListBoolValues.at("AnalyticPRATE");
  int idxForc = eBlPROC.ListIntValues.at("idxForc");
  //
  // Retrieving the grid array
  //
  GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
  //
  // Setting up the timings.
  //
  ArrayHistory eArr = ReadArrayHistory(eTriple);
  std::vector<double> ListTime = GetIntervalGen(eBlPROC, {eArr});
  int nbTime = ListTime.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  //
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  double FirstTime = ListTime[0];
  double LastTime = ListTime[nbTime - 1];
  double eps = 0.0001;
  int iDayFirst = static_cast<int>(floor(FirstTime + eps));
  int iDayLast = static_cast<int>(floor(LastTime + eps));
  std::string OutPrefix = eBlPROC.ListStringValues.at("OutPrefix");
  std::string eDir = FILE_GetAbsoluteDirectory(OutPrefix);
  std::string eDirB = ExtractDirectoryFromFileString(eDir);
  CreateDirectory(eDirB);
  //
  // Now writing data
  //
  for (int iDay = iDayFirst; iDay <= iDayLast; iDay++) {
    std::cerr
        << "-------------------------------------------------------------\n";
    std::vector<int> ListITime;
    bool IsFirst = true;
    double MinTimeFirst = -1;
    for (int iTime = 0; iTime < nbTime; iTime++) {
      double eTime = ListTime[iTime];
      if (eTime > static_cast<double>(iDay) - eps &&
          eTime < static_cast<double>(iDay + 1) - eps) {
        ListITime.push_back(iTime);
        if (IsFirst) {
          MinTimeFirst = eTime;
          IsFirst = false;
        } else {
          if (eTime < MinTimeFirst)
            MinTimeFirst = eTime;
        }
      }
    }
    std::vector<int> eDateMinSix = DATE_ConvertMjd2six(MinTimeFirst);
    std::vector<int> eTimeFirstSix{
        eDateMinSix[0], eDateMinSix[1], eDateMinSix[2], 0, 0, 0};
    double eTimeFirst = DATE_ConvertSix2mjd(eTimeFirstSix);
    std::string strPresB = DATE_ConvertSix2mystringPres(eTimeFirstSix);
    int nbTimeRel = ListITime.size();
    std::cerr << "iDay=" << iDay << " date=" << strPresB
              << " nbTimeRel=" << nbTimeRel << "\n";
    //
    int idx = iDay - (iDayFirst - 1);
    //
    // Common variables
    //
    std::vector<std::string> ListDimTime{"time"};
    std::vector<std::string> ListDimTimeStr{"time", "dateString"};
    std::vector<std::string> ListDimABsize{"aSize", "bSize"};
    std::vector<std::string> ListDimField{"time", "aSize", "bSize"};
    int Data[4];
    std::vector<int> eDate = DATE_ConvertMjd2six(eTimeFirst);
    for (int i = 0; i < 4; i++)
      Data[i] = eDate[i];
    //
    struct InfoVarWrite {
      std::string VarName;
      std::string FullName;
      std::string Units;
      std::string PolyAccess;
      double DefaultValue;
      bool DoAnalytic;
    };
    auto WriteFullVariable =
        [&](netCDF::NcFile &recFile, netCDF::NcVar &eVarData_time,
            netCDF::NcVar &eVarData_timeStr, InfoVarWrite const &eV) -> void {
      netCDF::NcVar eVarData =
          recFile.addVar(eV.VarName, "float", ListDimField);
      eVarData.putAtt("long_name", eV.FullName);
      eVarData.putAtt("units", eV.Units);
      MyMatrix<double> VAR(eta_rho, xi_rho);
      std::vector<float> XfieldD(eta_rho * xi_rho);
      for (int iTimeRel = 0; iTimeRel < nbTimeRel; iTimeRel++) {
        //
        // Data loading
        //
        int iTime = ListITime[iTimeRel];
        double eTimeDay = ListTime[iTime];
        if (eV.DoAnalytic) {
          for (int i = 0; i < eta_rho; i++)
            for (int j = 0; j < xi_rho; j++)
              VAR(i, j) = eV.DefaultValue;
        } else {
          RecVar eRecVar =
              ModelSpecificVarSpecificTime(TotalArr, eV.PolyAccess, eTimeDay);
          VAR = eRecVar.F;
        }
        std::string strPresC = DATE_ConvertMjd2mystringPres(eTimeDay);
        std::vector<size_t> startpT{size_t(iTimeRel)};
        std::vector<size_t> countpT{1};
        double eTimePrint = eTimeDay - eTimeFirst;
        eVarData_time.putVar(startpT, countpT, &eTimePrint);
        std::vector<size_t> startpK{size_t(iTimeRel), 0};
        std::vector<size_t> countpK{1, 19};
        eVarData_timeStr.putVar(startpK, countpK, strPresC.c_str());
        //
        std::vector<size_t> startp{size_t(iTimeRel), 0, 0};
        std::vector<size_t> countp{1, size_t(eta_rho), size_t(xi_rho)};
        int posA = 0;
        for (int i = 0; i < eta_rho; i++)
          for (int j = 0; j < xi_rho; j++) {
            XfieldD[posA] = static_cast<float>(VAR(i, j));
            posA++;
          }
        eVarData.putVar(startp, countp, XfieldD.data());
      }
    };
    auto DoWriteDown = [&](bool const &DoAIR, bool const &DoPRC,
                           bool const &DoRAD) -> void {
      std::string eFileNC;
      std::string Postfix =
          IntToString(idxForc) + "." + StringNumber(idx, 3) + ".nc";
      if (DoAIR)
        eFileNC = OutPrefix + "sflux_air_" + Postfix;
      if (DoPRC)
        eFileNC = OutPrefix + "sflux_prc_" + Postfix;
      if (DoRAD)
        eFileNC = OutPrefix + "sflux_rad_" + Postfix;
      if (!FILE_IsFileMakeable(eFileNC)) {
        std::cerr << "Request to create file FileOut=" << eFileNC << "\n";
        std::cerr << "but the directory does not exist\n";
        throw TerminalException{1};
      }
      netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                              netCDF::NcFile::nc4);
      netCDF::NcDim eDimAsize = dataFile.addDim("aSize", eta_rho);
      netCDF::NcDim eDimBsize = dataFile.addDim("bSize", xi_rho);
      netCDF::NcDim eDimNb = dataFile.addDim("time");
      netCDF::NcDim eDimDate = dataFile.addDim("dateString", 19);
      //
      // Add variable of lon/lat and time
      //
      netCDF::NcVar eVarData_lon =
          dataFile.addVar("lon", "float", ListDimABsize);
      eVarData_lon.putAtt("long_name", std::string("longitude"));
      eVarData_lon.putAtt("units", std::string("degrees_east"));
      eVarData_lon.putAtt("standard_name", std::string("longitude"));
      netCDF::NcVar eVarData_lat =
          dataFile.addVar("lat", "float", ListDimABsize);
      eVarData_lat.putAtt("long_name", std::string("latitude"));
      eVarData_lat.putAtt("units", std::string("degrees_north"));
      eVarData_lat.putAtt("standard_name", std::string("latitude"));
      netCDF::NcVar eVarData_ang =
          dataFile.addVar("ang", "float", ListDimABsize);
      netCDF::NcVar eVarData_time =
          dataFile.addVar("time", "double", ListDimTime);
      eVarData_time.putAtt("long_name", std::string("Time"));
      eVarData_time.putAtt("standard_name", std::string("time"));
      std::string strPresFirst = DATE_ConvertMjd2mystringPres(eTimeFirst);
      std::string AttUnitTime = "days since " + strPresFirst;
      eVarData_time.putAtt("base_date", netCDF::NcType::nc_INT, 4, Data);
      eVarData_time.putAtt("units", AttUnitTime);
      netCDF::NcVar eVarData_timeStr =
          dataFile.addVar("time_str", "char", ListDimTimeStr);
      //
      // some var write
      //
      auto StdGridWrite = [&](netCDF::NcVar &eVarData,
                              MyMatrix<double> const &F) -> void {
        int posA = 0;
        std::vector<float> XfieldD(eta_rho * xi_rho);
        for (int i = 0; i < eta_rho; i++)
          for (int j = 0; j < xi_rho; j++) {
            double eLon = F(i, j);
            XfieldD[posA] = static_cast<float>(eLon);
            posA++;
          }
        eVarData.putVar(XfieldD.data());
      };
      StdGridWrite(eVarData_lon, GrdArr.GrdArrRho.LON);
      StdGridWrite(eVarData_lat, GrdArr.GrdArrRho.LAT);
      StdGridWrite(eVarData_ang, GrdArr.GrdArrRho.ANG);

      if (DoAIR) {
        InfoVarWrite eV_prmsl{"prmsl",  "Pressure reduced to MSL",
                              "Pa",     "SurfPres",
                              105384.9, AnalyticPrmsl};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_prmsl);
        InfoVarWrite eV_uwind{"uwind", "10 m U-wind", "m/s", "Uwind",
                              0,       AnalyticWind};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_uwind);
        InfoVarWrite eV_vwind{"vwind", "10 m V-wind", "m/s", "Vwind",
                              0,       AnalyticWind};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_vwind);
        InfoVarWrite eV_stmp{"stmp",      "Surface Air Temperature (2m AGL)",
                             "K",         "AIRT2K",
                             273.15 + 15, AnalyticStmp};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_stmp);
        InfoVarWrite eV_spfh{"spfh", "Surface Specific Humidity (2m AGL)",
                             "1",    "Rh2frac",
                             0,      AnalyticSpfh};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_spfh);
      }
      if (DoPRC) {
        InfoVarWrite eV_prate{"prate", "precipitation flux", "kg/m^2/s", "rain",
                              0,       AnalyticPrate};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_prate);
      }
      if (DoRAD) {
        InfoVarWrite eV_dlwrf{"dlwrf", "Downward Long Wave Radiation Flux",
                              "W/m^2", "lwrad",
                              0,       AnalyticDlwrf};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_dlwrf);
        InfoVarWrite eV_dswrf{"dswrf", "Downward Short Wave Radiation Flux",
                              "W/m^2", "swrad",
                              0,       AnalyticDswrf};
        WriteFullVariable(dataFile, eVarData_time, eVarData_timeStr, eV_dswrf);
      }
    };
    if (DoAIRfiles)
      DoWriteDown(true, false, false);
    if (DoPRCfiles)
      DoWriteDown(false, true, false);
    if (DoRADfiles)
      DoWriteDown(false, false, true);
  }
}

MyMatrix<double>
SingleInterpolationOfField_2D(SingleArrayInterpolation const &eInterp,
                              MyMatrix<double> const &Fin) {
  int eta_out = eInterp.eta_out;
  int xi_out = eInterp.xi_out;
  int eta_in = eInterp.eta_in;
  int xi_in = eInterp.xi_in;
  int dimOut = eta_out * xi_out;
  int dimIn = eta_in * xi_in;
  int nbRow = eInterp.SpMat.rows();
  int nbCol = eInterp.SpMat.cols();
  if (dimOut != nbRow) {
    std::cerr << "Output Dimension inconsistency\n";
    std::cerr << "dimOut=" << dimOut << " rows=" << nbRow << "\n";
    throw TerminalException{1};
  }
  if (dimIn != nbCol) {
    std::cerr << "Input Dimension inconsistency\n";
    std::cerr << "dimIn=" << dimIn << " cols=" << nbCol << "\n";
    throw TerminalException{1};
  }
  MyVector<double> xFin(eta_in * xi_in);
  int idxIn = 0;
  for (int iXi = 0; iXi < xi_in; iXi++)
    for (int iEta = 0; iEta < eta_in; iEta++) {
      xFin(idxIn) = Fin(iEta, iXi);
      idxIn++;
    }
  //  std::cerr << "eInterp.SpMat(rows/cols)=" << eInterp.SpMat.rows() << " / "
  //  << eInterp.SpMat.cols() << "\n"; std::cerr << "xFin assigned\n";
  MyVector<double> xFout = eInterp.SpMat * xFin;
  //  std::cerr << "xFout assigned\n";
  int idxOut = 0;
  MyMatrix<double> Fout(eta_out, xi_out);
  for (int iXi = 0; iXi < xi_out; iXi++)
    for (int iEta = 0; iEta < eta_out; iEta++) {
      Fout(iEta, iXi) = xFout(idxOut);
      idxOut++;
    }
  //  std::cerr << "Fout assigned\n";
  return Fout;
}

void LevelPrinting(std::string const &VarName,
                   Eigen::Tensor<double, 3> const &F) {
  auto LDim = F.dimensions();
  int s_rho = LDim[0];
  double eMin = minCoeff(F);
  double eMax = maxCoeff(F);
  double deltaMM = eMax - eMin;
  std::cerr << "LevelPrinting of " << VarName << " Global : min/max = " << eMin
            << " / " << eMax << " deltaMM=" << deltaMM << "\n";
  for (int iS = 0; iS < s_rho; iS++) {
    MyMatrix<double> Fmat = DimensionExtraction(F, 0, iS);
    std::cerr << "  " << VarName << " iS=" << iS
              << " min/max=" << Fmat.minCoeff() << " / " << Fmat.maxCoeff()
              << "\n";
  }
}

Eigen::Tensor<double, 3> SingleInterpolationOfField_3D_horizontal(
    SingleArrayInterpolationGen const &eInterp,
    Eigen::Tensor<double, 3> const &Fin) {
  LevelPrinting("Fin", Fin);
  auto LDim = Fin.dimensions();
  int Nvert = LDim[0];
  int eta_out = eInterp.e_arr.eta_out;
  int xi_out = eInterp.e_arr.xi_out;
  int eta_in = eInterp.e_arr.eta_in;
  int xi_in = eInterp.e_arr.xi_in;
  MyVector<double> xFin(eta_in * xi_in);
  MyVector<double> xFout(eta_out * xi_out);
  Eigen::Tensor<double, 3> Fout(Nvert, eta_out, xi_out);
  double threshold = 1000000000;
  for (int iVert = 0; iVert < Nvert; iVert++) {
    int idxIn = 0;
    int n_error = 0;
    for (int iXi = 0; iXi < xi_in; iXi++) {
      for (int iEta = 0; iEta < eta_in; iEta++) {
        xFin(idxIn) = Fin(iVert, iEta, iXi);
        if (eInterp.e_arr.ARVDin.TensMSKvert) {
          auto const &tens = *eInterp.e_arr.ARVDin.TensMSKvert;
          if (tens(iVert, iEta, iXi) == 1) {
            if (fabs(Fin(iVert, iEta, iXi)) > threshold)
              n_error += 1;
          }
        }
        idxIn++;
      }
    }
    std::cerr << "iVert=" << iVert << " n_error=" << n_error << "\n";
    if (!eInterp.e_arr.ARVDin.TensMSKvert) {
      xFout = eInterp.e_arr.SpMat * xFin;
    } else {
      xFout = eInterp.l_arr[iVert].SpMat * xFin;
    }
    int idxOut = 0;
    for (int iXi = 0; iXi < xi_out; iXi++)
      for (int iEta = 0; iEta < eta_out; iEta++) {
        Fout(iVert, iEta, iXi) = xFout(idxOut);
        idxOut++;
      }
  }
  LevelPrinting("Fout", Fout);
  return Fout;
}

Eigen::Tensor<double, 3>
SingleInterpolationOfField_3D(SingleArrayInterpolationGen const &eInterp,
                              Eigen::Tensor<double, 3> const &Fin) {
  std::cerr << "eInterp.GrdArrOut.ARVD.IsAssigned="
            << eInterp.e_arr.GrdArrOut.ARVD.IsAssigned
            << " eInterp.ARVDin.IsAssigned=" << eInterp.e_arr.ARVDin.IsAssigned
            << "\n";
  if (!eInterp.e_arr.GrdArrOut.ARVD.IsAssigned ||
      !eInterp.e_arr.ARVDin.IsAssigned)
    return SingleInterpolationOfField_3D_horizontal(eInterp, Fin);
  Eigen::Tensor<double, 3> Fhoriz =
      SingleInterpolationOfField_3D_horizontal(eInterp, Fin);
  //  std::cerr << "Before vertical interpolation\n";
  Eigen::Tensor<double, 3> Fret = VerticalInterpolationTensor_R(
      eInterp.e_arr.GrdArrOut, eInterp.e_arr.ARVDin, eInterp.e_arr.DEPinInterp,
      Fhoriz);
  LevelPrinting("Fret", Fret);
  //  std::cerr << "After vertical interpolation\n";
  /*
  LevelPrinting("Fin", Fin);
  LevelPrinting("Fhoriz", Fhoriz);
  LevelPrinting("Fret", Fret);
  */
  return Fret;
}

SingleArrayInterpolation
GetSingleArrayInterpolationTrivialCase(GridArray const &GrdArrOut,
                                       GridArray const &GrdArrIn) {
  std::cerr << "GetSingleArrayInterpolationTrivialCase, beginning\n";
  int eta_out = GrdArrOut.GrdArrRho.LON.rows();
  int xi_out = GrdArrOut.GrdArrRho.LON.cols();
  int eta_in = GrdArrIn.GrdArrRho.LON.rows();
  int xi_in = GrdArrIn.GrdArrRho.LON.cols();
  double deltaLON =
      GrdArrOut.GrdArrRho.LON(1, 0) - GrdArrOut.GrdArrRho.LON(0, 0);
  double deltaLAT =
      GrdArrOut.GrdArrRho.LAT(0, 1) - GrdArrOut.GrdArrRho.LAT(0, 0);
  double eQuad_MinLon = GrdArrOut.GrdArrRho.LON(0, 0);
  double eQuad_MinLat = GrdArrOut.GrdArrRho.LAT(0, 0);
  MyMatrix<uint8_t> const &MSK_in = GrdArrIn.GrdArrRho.MSK;
  //  MyMatrix<uint8_t> const& MSK_out=GrdArrOut.GrdArrRho.MSK;
  int nbNodeIn = eta_in * xi_in;
  int nbNodeFD = eta_out * xi_out;
  //  std::cerr << "eta_out=" << eta_out << " xi_out=" << xi_out << "\n";
  //  std::cerr << "eta_in =" << eta_in  << " xi_in =" << xi_in  << "\n";
  //
  // Now computing the sparse matrix for interpolation and the mask.
  //
  double THR = 1e-10;
  typedef Eigen::Triplet<double> T2;
  auto ComputeInterpolationArray =
      [&](MyMatrix<int> const &INE, MyMatrix<double> const &LON,
          MyMatrix<double> const &LAT,
          int const &nbNodeInput) -> SingleArrayInterpolation {
    int nbTrig = INE.rows();
    MyMatrix<int> MatITrig(eta_out, xi_out);
    for (int i = 0; i < eta_out; i++)
      for (int j = 0; j < xi_out; j++)
        MatITrig(i, j) = -1;
    for (int iTrig = 0; iTrig < nbTrig; iTrig++) {
      double MinLon = 0, MaxLon = 0, MinLat = 0, MaxLat = 0;
      for (int i = 0; i < 3; i++) {
        int ip = INE(iTrig, i);
        double eLon = LON(ip, 0);
        double eLat = LAT(ip, 0);
        if (i == 0) {
          MinLon = eLon;
          MaxLon = eLon;
          MinLat = eLat;
          MaxLat = eLat;
        } else {
          if (eLon > MaxLon)
            MaxLon = eLon;
          if (eLon < MinLon)
            MinLon = eLon;
          if (eLat > MaxLat)
            MaxLat = eLat;
          if (eLat < MinLat)
            MinLat = eLat;
        }
      }
      int iLonMin = static_cast<int>(floor((MinLon - eQuad_MinLon) / deltaLON));
      int iLatMin = static_cast<int>(floor((MinLat - eQuad_MinLat) / deltaLAT));
      int iLonMax = static_cast<int>(ceil((MaxLon - eQuad_MinLon) / deltaLON));
      int iLatMax = static_cast<int>(ceil((MaxLat - eQuad_MinLat) / deltaLAT));
      iLonMin = std::max(iLonMin, 0);
      iLonMax = std::min(iLonMax, eta_out - 1);
      iLatMin = std::max(iLatMin, 0);
      iLatMax = std::min(iLatMax, xi_out - 1);
      auto IsCorrect = [&](double const &Xp, double const &Yp) -> bool {
        int ki = INE(iTrig, 0);
        int kj = INE(iTrig, 1);
        int kk = INE(iTrig, 2);
        double xi = LON(ki);
        double yi = LAT(ki);
        double xj = LON(kj);
        double yj = LAT(kj);
        double xk = LON(kk);
        double yk = LAT(kk);
        double f1, f2, f3;
        f1 = xi * (yj - Yp) + xj * (Yp - yi) + Xp * (yi - yj);
        f2 = xj * (yk - Yp) + xk * (Yp - yj) + Xp * (yj - yk);
        f3 = xk * (yi - Yp) + xi * (Yp - yk) + Xp * (yk - yi);
        if (f1 > -THR && f2 > -THR && f3 > -THR)
          return true;
        return false;
      };
      for (int iLon = iLonMin; iLon <= iLonMax; iLon++)
        for (int iLat = iLatMin; iLat <= iLatMax; iLat++) {
          double eLon = GrdArrOut.GrdArrRho.LON(iLon, iLat);
          double eLat = GrdArrOut.GrdArrRho.LAT(iLon, iLat);
          bool test = IsCorrect(eLon, eLat);
          if (test) {
            MatITrig(iLon, iLat) = iTrig;
          }
        }
    }
    int nbWet = 0;
    for (int iLon = 0; iLon < eta_out; iLon++)
      for (int iLat = 0; iLat < xi_out; iLat++) {
        int iTrig = MatITrig(iLon, iLat);
        if (iTrig != -1)
          nbWet++;
      }
    int nnz = 3 * nbWet;
    std::vector<T2> ListTr(nnz);
    int iNnz = 0;
    for (int iLon = 0; iLon < eta_out; iLon++)
      for (int iLat = 0; iLat < xi_out; iLat++) {
        int iTrig = MatITrig(iLon, iLat);
        if (iTrig != -1) {
          double Xp = GrdArrOut.GrdArrRho.LON(iLon, iLat);
          double Yp = GrdArrOut.GrdArrRho.LAT(iLon, iLat);
          std::vector<double> Xcall(3), Ycall(3);
          std::vector<int> LEta(3);
          for (int i = 0; i < 3; i++) {
            int IP = INE(iTrig, i);
            double eLon = LON(IP);
            double eLat = LAT(IP);
            LEta[i] = IP;
            Xcall[i] = eLon;
            Ycall[i] = eLat;
          }
          std::vector<double> LCoeff =
              DetermineCoefficient(Xcall, Ycall, Xp, Yp);
          int iNodeFD = iLon + eta_out * iLat;
          for (int i = 0; i < 3; i++) {
            T2 eTr{iNodeFD, LEta[i], LCoeff[i]};
            ListTr[iNnz] = eTr;
            iNnz++;
          }
        }
      }
    MySparseMatrix<double> SpMat(nbNodeFD, nbNodeInput);
    SpMat.setFromTriplets(ListTr.begin(), ListTr.end());
    return {eta_out,
            xi_out,
            eta_in,
            xi_in,
            std::move(GrdArrOut),
            std::move(GrdArrIn.ARVD),
            std::move(SpMat),
            {}};
  };
  std::cerr << "GrdArrIn.IsFE=" << GrdArrIn.IsFE << "\n";
  if (GrdArrIn.IsFE == 1) {
    return ComputeInterpolationArray(GrdArrIn.INE, GrdArrIn.GrdArrRho.LON,
                                     GrdArrIn.GrdArrRho.LAT, nbNodeIn);
  } else {
    int nbWet = GrdArrIn.GrdArrRho.nbWet;
    std::cerr << "nbWet=" << nbWet << "\n";
    std::cerr << "eta_in=" << eta_in << " xi_in=" << xi_in << "\n";
    MyMatrix<int> MatDirect(eta_in, xi_in);
    MyVector<int> EtaRev(nbWet), XiRev(nbWet);
    int idx = 0;
    MyMatrix<double> LON(nbWet, 1);
    MyMatrix<double> LAT(nbWet, 1);
    for (int i = 0; i < eta_in; i++)
      for (int j = 0; j < xi_in; j++) {
        int eVal = -1;
        if (MSK_in(i, j) == 1) {
          eVal = idx;
          EtaRev(idx) = i;
          XiRev(idx) = j;
          LON(idx, 0) = GrdArrIn.GrdArrRho.LON(i, j);
          LAT(idx, 0) = GrdArrIn.GrdArrRho.LAT(i, j);
          idx++;
        }
        MatDirect(i, j) = eVal;
      }
    std::vector<MyVector<int>> ListVectINE;
    for (int i = 0; i < eta_in - 1; i++)
      for (int j = 0; j < xi_in - 1; j++) {
        int sumMSK = MSK_in(i, j) + MSK_in(i + 1, j) + MSK_in(i, j + 1) +
                     MSK_in(i + 1, j + 1);
        if (sumMSK == 4) {
          MyVector<int> eVect1(3);
          eVect1(0) = MatDirect(i + 1, j);
          eVect1(1) = MatDirect(i, j);
          eVect1(2) = MatDirect(i, j + 1);
          ListVectINE.push_back(eVect1);
          //
          MyVector<int> eVect2(3);
          eVect2(0) = MatDirect(i, j + 1);
          eVect2(1) = MatDirect(i + 1, j + 1);
          eVect2(2) = MatDirect(i + 1, j);
          ListVectINE.push_back(eVect2);
        }
        if (sumMSK == 3 && MSK_in(i, j) == 0) {
          MyVector<int> eVect1(3);
          eVect1(0) = MatDirect(i, j + 1);
          eVect1(1) = MatDirect(i + 1, j + 1);
          eVect1(2) = MatDirect(i + 1, j);
          ListVectINE.push_back(eVect1);
        }
        if (sumMSK == 3 && MSK_in(i + 1, j + 1) == 0) {
          MyVector<int> eVect1(3);
          eVect1(0) = MatDirect(i + 1, j);
          eVect1(1) = MatDirect(i, j);
          eVect1(2) = MatDirect(i, j + 1);
          ListVectINE.push_back(eVect1);
        }
        if (sumMSK == 3 && MSK_in(i + 1, j) == 0) {
          MyVector<int> eVect1(3);
          eVect1(0) = MatDirect(i + 1, j + 1);
          eVect1(1) = MatDirect(i, j);
          eVect1(2) = MatDirect(i, j + 1);
          ListVectINE.push_back(eVect1);
        }
        if (sumMSK == 3 && MSK_in(i, j + 1) == 0) {
          MyVector<int> eVect1(3);
          eVect1(0) = MatDirect(i + 1, j + 1);
          eVect1(1) = MatDirect(i + 1, j);
          eVect1(2) = MatDirect(i, j);
          ListVectINE.push_back(eVect1);
        }
      }
    int nbTrig = ListVectINE.size();
    MyMatrix<int> INE(nbTrig, 3);
    for (int iTrig = 0; iTrig < nbTrig; iTrig++)
      for (int i = 0; i < 3; i++)
        INE(iTrig, i) = ListVectINE[iTrig](i);
    MySparseMatrix<double> SpMat =
        ComputeInterpolationArray(INE, LON, LAT, nbWet).SpMat;
    int nnz = SpMat.nonZeros();
    int iNNZ = 0;
    std::vector<T2> ListTr(nnz);
    for (int k = 0; k < SpMat.outerSize(); k++)
      for (typename MySparseMatrix<double>::InnerIterator it(SpMat, k); it;
           ++it) {
        double eVal = it.value();
        int iRow = it.row();
        int iCol = it.col();
        int iEta = EtaRev(iCol);
        int iXi = XiRev(iCol);
        int iColNew = iEta + eta_in * iXi;
        ListTr[iNNZ] = T2(iRow, iColNew, eVal);
        iNNZ++;
      }
    MySparseMatrix<double> NewSpMat(nbNodeFD, nbNodeIn);
    NewSpMat.setFromTriplets(ListTr.begin(), ListTr.end());
    return {eta_out,
            xi_out,
            eta_in,
            xi_in,
            std::move(GrdArrOut),
            std::move(GrdArrIn.ARVD),
            std::move(NewSpMat),
            {}};
  }
}

SingleArrayInterpolation
ConvertToArrayInt(int const &eta_out, int const &xi_out, int const &eta_in,
                  int const &xi_in, std::vector<int> const &LEta,
                  std::vector<int> const &LXi,
                  std::vector<SingleRecInterp> const &LRec,
                  GridArray const &GrdArrOut, ARVDtyp const &ARVDin) {
  typedef Eigen::Triplet<double> T2;
  std::vector<T2> tripletList;
  int nbEnt = LEta.size();
  std::vector<std::pair<int, int>> ListMiss;
  for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
    int iEta_out = LEta[iEnt];
    int iXi_out = LXi[iEnt];
    int iNodeOut = iEta_out + eta_out * iXi_out;
    SingleRecInterp eRec = LRec[iEnt];
    if (eRec.status) {
      int len = eRec.LPart.size();
      for (int i = 0; i < len; i++) {
        SinglePartInterp ePart = eRec.LPart[i];
        int iEta_in = ePart.eEta;
        int iXi_in = ePart.eXi;
        double eCoeff = ePart.eCoeff;
        int iNodeIn = iEta_in + eta_in * iXi_in;
        int iRow = iNodeOut;
        int iCol = iNodeIn;
        T2 eTr(iRow, iCol, eCoeff);
        tripletList.push_back(eTr);
      }
    } else {
      ListMiss.push_back({iEta_out, iXi_out});
    }
  }
  std::cerr << "|ListMiss|=" << ListMiss.size() << "\n";
  for (auto &eMiss : ListMiss) {
    int i = eMiss.first;
    int j = eMiss.second;
    std::cerr << " eMiss=[" << i << "," << j
              << "] lon=" << GrdArrOut.GrdArrRho.LON(i, j)
              << " lat=" << GrdArrOut.GrdArrRho.LAT(i, j) << "\n";
  }
  int nbRow = eta_out * xi_out;
  int nbCol = eta_in * xi_in;
  MySparseMatrix<double> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return {eta_out, xi_out, eta_in, xi_in, std::move(GrdArrOut),
    std::move(ARVDin), std::move(SpMat), {}};
}

SingleArrayInterpolation
INTERPOL_CreateSingleRecVarInterpol(GridArray const &GrdArrOut,
                                    GridArray const &GrdArrIn,
                                    bool const &AllowExtrapolation) {
  //  std::cerr << "Begining of INTERPOL_CreateSingleRecVarInterpol\n";
  int eta_in = GrdArrIn.GrdArrRho.LON.rows();
  int xi_in = GrdArrIn.GrdArrRho.LON.cols();
  int eta_out = GrdArrOut.GrdArrRho.LON.rows();
  int xi_out = GrdArrOut.GrdArrRho.LON.cols();
  SingleArrayInterpolation RecArr;
  if (GrdArrOut.IsFE == 1) {
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 1\n";
    MyMatrix<double> ListXY(2, eta_out);
    std::vector<int> LEta(eta_out), LXi(eta_out);
    for (int i = 0; i < eta_out; i++) {
      LEta[i] = i;
      LXi[i] = 0;
      ListXY(0, i) = GrdArrOut.GrdArrRho.LON(i, 0);
      ListXY(1, i) = GrdArrOut.GrdArrRho.LAT(i, 0);
    }
    std::vector<SingleRecInterp> LSingle =
        General_FindInterpolationWeight(GrdArrIn, ListXY, AllowExtrapolation);
    RecArr = ConvertToArrayInt(eta_out, xi_out, eta_in, xi_in, LEta, LXi,
                               LSingle, GrdArrOut, GrdArrIn.ARVD);
  } else {
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2\n";
    if (GrdArrOut.ModelName == "RECTANGULAR" && GrdArrIn.IsFE == 1) {
      return GetSingleArrayInterpolationTrivialCase(GrdArrOut, GrdArrIn);
    }
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.1\n";
    int nbWet = GrdArrOut.GrdArrRho.nbWet;
    MyMatrix<double> ListXY(2, nbWet);
    std::vector<int> LEta(nbWet), LXi(nbWet);
    int idx = 0;
    for (int i = 0; i < eta_out; i++)
      for (int j = 0; j < xi_out; j++)
        if (GrdArrOut.GrdArrRho.MSK(i, j) == 1) {
          LEta[idx] = i;
          LXi[idx] = j;
          ListXY(0, idx) = GrdArrOut.GrdArrRho.LON(i, j);
          ListXY(1, idx) = GrdArrOut.GrdArrRho.LAT(i, j);
          idx++;
        }
    if (nbWet != idx) {
      std::cerr << "Inconsistency in number of wet points\n";
      std::cerr << "nbWet=" << nbWet << " idx=" << idx << "\n";
      throw TerminalException{1};
    }
    std::cerr << "|ListXY|=" << ListXY.cols() << " / " << ListXY.rows() << "\n";
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.2\n";
    std::vector<SingleRecInterp> LSingle =
        General_FindInterpolationWeight(GrdArrIn, ListXY, AllowExtrapolation);
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.3\n";
    RecArr = ConvertToArrayInt(eta_out, xi_out, eta_in, xi_in, LEta, LXi,
                               LSingle, GrdArrOut, GrdArrIn.ARVD);
    //    std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.4\n";
  }
  //  std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.5\n";
  if (GrdArrIn.GrdArrRho.DEP)
    RecArr.DEPinInterp =
        SingleInterpolationOfField_2D(RecArr, GetDEP(GrdArrIn.GrdArrRho));
  //  std::cerr << "INTERPOL_CreateSingleRecVarInterpol, case 2.6\n";
  return RecArr;
}

SingleArrayInterpolationGen
INTERPOL_CreateSingleRecVarInterpolGen(GridArray const &GrdArrOut,
                                       GridArray const &GrdArrIn,
                                       bool const &AllowExtrapolation) {
  if (!GrdArrIn.ARVD.TensMSKvert)
    return {INTERPOL_CreateSingleRecVarInterpol(GrdArrOut, GrdArrIn,
                                                AllowExtrapolation),
            {}};
  if (!AllowExtrapolation) {
    std::cerr
        << "In the case of Z-coordinates, we need to allow for Extrapolation\n";
    throw TerminalException{1};
  }
  std::vector<SingleArrayInterpolation> l_arr;
  Eigen::Tensor<uint8_t, 3> const &TensMSKvert = *GrdArrIn.ARVD.TensMSKvert;
  auto LDim = TensMSKvert.dimensions();
  int Nvert = LDim[0];
  int eta = LDim[1];
  int xi = LDim[2];
  GridArray GrdArrIn_modif = GrdArrIn;
  for (int iVert = 0; iVert < Nvert; iVert++) {
    MyMatrix<uint8_t> MSK(eta, xi);
    int nWet = 0;
    for (int i = 0; i < eta; i++)
      for (int j = 0; j < xi; j++) {
        MSK(i, j) = TensMSKvert(iVert, i, j);
        nWet += static_cast<int>(TensMSKvert(iVert, i, j));
      }
    std::cerr << "iVert=" << iVert << " / " << Nvert << " nWet=" << nWet
              << "\n";
    GrdArrIn_modif.GrdArrRho.MSK = MSK;
    l_arr.push_back(INTERPOL_CreateSingleRecVarInterpol(
        GrdArrOut, GrdArrIn_modif, AllowExtrapolation));
  }
  SingleArrayInterpolation e_arr = INTERPOL_CreateSingleRecVarInterpol(
      GrdArrOut, GrdArrIn, AllowExtrapolation);
  return {e_arr, l_arr};
}

MyMatrix<double> RegionAveraging_2D(SingleArrayRegionAveraging const &eInterp,
                                    MyMatrix<double> const &F) {
  int len = eInterp.ListListEtaXi.size();
  MyMatrix<double> Fret(len, 1);
  for (int i = 0; i < len; i++) {
    double sum = 0;
    for (auto &ePair : eInterp.ListListEtaXi[i]) {
      sum += F(ePair.first, ePair.second);
    }
    Fret(i, 0) = sum / eInterp.ListListEtaXi[i].size();
  }
  return Fret;
}

Eigen::Tensor<double, 3>
RegionAveraging_3D(SingleArrayRegionAveraging const &eInterp,
                   Eigen::Tensor<double, 3> const &F) {
  int len = eInterp.ListListEtaXi.size();
  auto LDim = F.dimensions();
  int Nvert = LDim[0];
  Eigen::Tensor<double, 3> Fret(Nvert, len, 1);
  for (int i = 0; i < len; i++) {
    MyVector<double> Vsum = ZeroVector<double>(Nvert);
    for (auto &ePair : eInterp.ListListEtaXi[i]) {
      for (int iV = 0; iV < Nvert; iV++)
        Vsum(iV) += F(iV, ePair.first, ePair.second);
    }
    for (int iV = 0; iV < Nvert; iV++)
      Fret(iV, i, 0) = Vsum(iV) / eInterp.ListListEtaXi[i].size();
  }
  return Fret;
}

RecVar REGAVE_SingleRecVarAveraging(SingleArrayRegionAveraging const &eInterp,
                                    RecVar const &fRecVar) {
  //  std::cerr << "REGAVE_SingleRecVarAveraging, begin\n";
  RecVar eRecVar;
  eRecVar.RecS = fRecVar.RecS;
  //  std::cerr << "VarNature=" << fRecVar.RecS.VarNature << "\n";
  if (fRecVar.RecS.VarNature == "rho") {
    //    std::cerr << "fRecVar.F min/max=" << fRecVar.F.minCoeff() << " / " <<
    //    fRecVar.F.maxCoeff() << "\n";
    eRecVar.F = RegionAveraging_2D(eInterp, fRecVar.F);
    //    std::cerr << "eRecVar.F min/max=" << eRecVar.F.minCoeff() << " / " <<
    //    eRecVar.F.maxCoeff() << "\n";
  }
  if (fRecVar.RecS.VarNature == "uv") {
    eRecVar.U = RegionAveraging_2D(eInterp, fRecVar.U);
    eRecVar.V = RegionAveraging_2D(eInterp, fRecVar.V);
    eRecVar.F = RegionAveraging_2D(eInterp, fRecVar.F);
  }
  if (fRecVar.RecS.VarNature == "3Drho") {
    eRecVar.Tens3 = RegionAveraging_3D(eInterp, fRecVar.Tens3);
  }
  if (fRecVar.RecS.VarNature == "3Duv") {
    eRecVar.Uthree = RegionAveraging_3D(eInterp, fRecVar.Uthree);
    eRecVar.Vthree = RegionAveraging_3D(eInterp, fRecVar.Vthree);
    eRecVar.Tens3 = RegionAveraging_3D(eInterp, fRecVar.Tens3);
  }
  //  std::cerr << "REGAVE_SingleRecVarAveraging, end\n";
  return eRecVar;
}

RecVar
INTERPOL_SingleRecVarInterpolation(SingleArrayInterpolationGen const &eInterp,
                                   RecVar const &fRecVar) {
  RecVar eRecVar;
  eRecVar.RecS = fRecVar.RecS;
  if (fRecVar.RecS.VarNature == "rho") {
    //    std::cerr << "fRecVar.F min/max=" << fRecVar.F.minCoeff() << " / " <<
    //    fRecVar.F.maxCoeff() << "\n";
    eRecVar.F = SingleInterpolationOfField_2D(eInterp.e_arr, fRecVar.F);
    //    std::cerr << "eRecVar.F min/max=" << eRecVar.F.minCoeff() << " / " <<
    //    eRecVar.F.maxCoeff() << "\n";
  }
  if (fRecVar.RecS.VarNature == "uv") {
    eRecVar.U = SingleInterpolationOfField_2D(eInterp.e_arr, fRecVar.U);
    eRecVar.V = SingleInterpolationOfField_2D(eInterp.e_arr, fRecVar.V);
    eRecVar.F = SingleInterpolationOfField_2D(eInterp.e_arr, fRecVar.F);
  }
  if (fRecVar.RecS.VarNature == "3Drho") {
    eRecVar.Tens3 = SingleInterpolationOfField_3D(eInterp, fRecVar.Tens3);
  }
  if (fRecVar.RecS.VarNature == "3Duv") {
    eRecVar.Uthree = SingleInterpolationOfField_3D(eInterp, fRecVar.Uthree);
    eRecVar.Vthree = SingleInterpolationOfField_3D(eInterp, fRecVar.Vthree);
    eRecVar.Tens3 = SingleInterpolationOfField_3D(eInterp, fRecVar.Tens3);
  }
  return eRecVar;
}

MyMatrix<uint8_t> ComputeInsideMask(SingleArrayInterpolation const &eSingArr) {
  int eta_out = eSingArr.eta_out;
  int xi_out = eSingArr.xi_out;
  MyMatrix<uint8_t> F;
  F.setZero(eta_out, xi_out);
  int nnz = eSingArr.SpMat.nonZeros();
  int nb = 0;
  for (int k = 0; k < eSingArr.SpMat.outerSize(); ++k)
    for (typename MySparseMatrix<double>::InnerIterator it(eSingArr.SpMat, k);
         it; ++it) {
      int iRow = it.row();
      int iEta_out = iRow % eta_out;
      int iXi_out = (iRow - iEta_out) / eta_out;
      F(iEta_out, iXi_out) = 1;
      nb++;
    }
  if (nnz != nb) {
    std::cerr << "Counting error in ComputeInsideMask\n";
    throw TerminalException{1};
  }
  return F;
}

struct TotalArrayInterpolation {
  bool NeedInterp;
  int eta_rho, xi_rho;
  MyMatrix<int> MSK;
  int nbGrid;
  ARVDtyp ARVD;
  std::vector<MyMatrix<double>> ListHatFunction;
  std::vector<SingleArrayInterpolationGen> ListSingleArrayInterpolationGen;
  std::vector<TotalArrGetData> ListTotalArr;
  double StartTime;
  double EndTime;
};

RecVar INTERPOL_GetHatFunction(GridArray const &GrdArrOut,
                               double const &TheSize) {
  int eta_rho = GrdArrOut.GrdArrRho.LON.rows();
  int xi_rho = GrdArrOut.GrdArrRho.LON.cols();
  MyMatrix<double> F(eta_rho, xi_rho);
  RecVar eRecVar;
  eRecVar.RecS.VarName1 = "hatfct";
  eRecVar.RecS.VarName2 = "Hat function";
  eRecVar.RecS.minval = 0;
  eRecVar.RecS.maxval = 1;
  eRecVar.RecS.mindiff = -1;
  eRecVar.RecS.maxdiff = 1;
  eRecVar.RecS.Unit = "nondimensional";
  eRecVar.RecS.VarNature = "rho";
  if (GrdArrOut.IsFE == 1) {
    std::cerr << "Need to write the code for the hat function in\n";
    std::cerr << "INTERPOL_GetHatFunction\n";
    throw TerminalException{1};
  } else {
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        std::vector<double> LValX{
            static_cast<double>(1), static_cast<double>(i) / TheSize,
            static_cast<double>(eta_rho - 1 - i) / TheSize};
        std::vector<double> LValY{
            static_cast<double>(1), static_cast<double>(j) / TheSize,
            static_cast<double>(xi_rho - 1 - j) / TheSize};
        double eValX = *std::min(LValX.begin(), LValX.end());
        double eValY = *std::min(LValY.begin(), LValY.end());
        F(i, j) = eValX * eValY;
      }
  }
  eRecVar.F = F;
  return eRecVar;
}

MyMatrix<double> HatFunctionFromMask(MyMatrix<uint8_t> const &MSKinput,
                                     GridArray const &GrdArr,
                                     GraphSparseImmutable const &eGR,
                                     int const &SpongeSize) {
  int eta_rho = MSKinput.rows();
  int xi_rho = MSKinput.cols();
  std::cerr << "HatFunctionFromMask computation\n";
  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
  std::vector<std::vector<int>> ListNeigh{{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  MyMatrix<int> TheMSKwork = ZeroMatrix<int>(eta_rho, xi_rho);
  struct Pair {
    size_t i;
    size_t j;
  };
  auto GetListAdjacent = [&](int const &i, int const &j) -> std::vector<Pair> {
    std::vector<Pair> TheRet;
    if (GrdArr.IsFE == 0) {
      for (int inei = 0; inei < 4; inei++) {
        int iN = i + ListNeigh[inei][0];
        int jN = j + ListNeigh[inei][1];
        if (iN >= 0 && iN < eta_rho && jN >= 0 && jN < xi_rho &&
            MSKinput(iN, jN) == 0)
          TheRet.push_back({size_t(iN), size_t(jN)});
      }
    } else {
      for (auto &eAdj : eGR.Adjacency(i))
        TheRet.push_back({size_t(eAdj), 0});
    }
    return TheRet;
  };
  std::cerr << "HatFunctionFromMask, step 1\n";
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      if (MSKinput(i, j) == 1) {
        int IsNeighZero = 0;
        std::vector<Pair> LPair = GetListAdjacent(i, j);
        for (auto &ePair : LPair)
          if (MSKinput(ePair.i, ePair.j) == 0)
            IsNeighZero = 1;
        if (IsNeighZero == 1)
          TheMSKwork(i, j) = 1;
      }
  std::cerr << "HatFunctionFromMask, step 2\n";
  for (int iVal = 2; iVal < SpongeSize; iVal++) {
    std::cerr << "  iVal=" << iVal << " / " << SpongeSize << "\n";
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        if (MSKinput(i, j) == 1 && TheMSKwork(i, j) == 0) {
          int IsNeighLower = 0;
          std::vector<Pair> LPair = GetListAdjacent(i, j);
          for (auto &ePair : LPair)
            if (TheMSKwork(ePair.i, ePair.j) == iVal - 1)
              IsNeighLower = 1;
          if (IsNeighLower == 1)
            TheMSKwork(i, j) = iVal;
        }
  }
  std::cerr << "HatFunctionFromMask, step 3\n";
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      if (MSKinput(i, j) == 1 && TheMSKwork(i, j) == 0)
        TheMSKwork(i, j) = SpongeSize;
  std::cerr << "HatFunctionFromMask, step 4\n";
  MyMatrix<double> TheFCT(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      TheFCT(i, j) = static_cast<double>(TheMSKwork(i, j)) /
                     static_cast<double>(SpongeSize);
  std::cerr << "HatFunctionFromMask, step 5\n";
  return TheFCT;
}

std::pair<GraphSparseImmutable, std::vector<std::pair<int, int>>>
GetGraphSparseVertexAdjacency(GridArray const &GrdArr) {
  if (GrdArr.IsFE == 1) {
    std::cerr << "GetGraphSparseVertexAdjacency : Unstructured scheme\n";
    int nbNode = GrdArr.GrdArrRho.LON.rows();
    std::vector<std::pair<int, int>> ListPoint(nbNode);
    for (int iNode = 0; iNode < nbNode; iNode++)
      ListPoint[iNode] = {iNode, 0};
    return {GetUnstructuredVertexAdjInfo(GrdArr.INE, nbNode),
            std::move(ListPoint)};
  } else {
    std::cerr << "GetGraphSparseVertexAdjacency : Structured scheme\n";
    // Determining the list of wet points.
    size_t eta_rho = GrdArr.GrdArrRho.LON.rows();
    size_t xi_rho = GrdArr.GrdArrRho.LON.cols();
    size_t miss_val = std::numeric_limits<size_t>::max();
    std::vector<std::pair<int, int>> ListPoint;
    MyMatrix<size_t> MappingIndex(eta_rho, xi_rho);
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++)
        MappingIndex(iEta, iXi) = miss_val;
    size_t index = 0;
    for (size_t iEta = 0; iEta < eta_rho; iEta++)
      for (size_t iXi = 0; iXi < xi_rho; iXi++)
        if (GrdArr.GrdArrRho.MSK(iEta, iXi) == 1) {
          std::pair<int, int> ePair{static_cast<int>(iEta),
                                    static_cast<int>(iXi)};
          ListPoint.push_back(ePair);
          MappingIndex(iEta, iXi) = index;
          index++;
        }
    size_t nb_point = index;
    std::cerr << "  nb_point=" << nb_point << " eta_rho=" << eta_rho
              << " xi_rho=" << xi_rho << "\n";
    // Building the adjacencies
    std::vector<size_t> NbAdj(nb_point, 0);
    std::vector<size_t> VectAdj(4 * nb_point, 0);
    auto GetADJ_index = [&](size_t const &iEta, size_t const &iXi) -> int {
      if (iEta >= eta_rho || iXi >= xi_rho)
        return miss_val;
      if (GrdArr.GrdArrRho.MSK(iEta, iXi) == 0)
        return miss_val;
      return MappingIndex(iEta, iXi);
    };
    auto GetADJ = [&](size_t const &idx, size_t const &iPoint) -> size_t {
      size_t iEta = ListPoint[iPoint].first;
      size_t iXi = ListPoint[iPoint].second;
      if (idx == 0)
        return GetADJ_index(iEta - 1, iXi);
      if (idx == 1)
        return GetADJ_index(iEta + 1, iXi);
      if (idx == 2)
        return GetADJ_index(iEta, iXi - 1);
      if (idx == 3)
        return GetADJ_index(iEta, iXi + 1);
      return miss_val;
    };
    for (size_t iPoint = 0; iPoint < nb_point; iPoint++) {
      for (size_t idx = 0; idx < 4; idx++) {
        size_t jPoint = GetADJ(idx, iPoint);
        if (jPoint != miss_val) {
          size_t eNB = NbAdj[iPoint];
          VectAdj[4 * iPoint + eNB] = jPoint;
          NbAdj[iPoint] = eNB + 1;
        }
      }
    }
    size_t nb_adj = 0;
    for (size_t iPoint = 0; iPoint < nb_point; iPoint++)
      nb_adj += NbAdj[iPoint];
    std::cerr << "  nb_adj=" << nb_adj << "\n";
    std::vector<size_t> ListStart(nb_point + 1, 0);
    for (size_t iPoint = 0; iPoint < nb_point; iPoint++)
      ListStart[iPoint + 1] = ListStart[iPoint] + NbAdj[iPoint];
    std::vector<size_t> ListListAdj(nb_adj);
    size_t pos = 0;
    for (size_t iPoint = 0; iPoint < nb_point; iPoint++) {
      size_t eNB = NbAdj[iPoint];
      for (size_t i = 0; i < eNB; i++) {
        ListListAdj[pos] = VectAdj[4 * iPoint + i];
        pos++;
      }
    }
    return {GraphSparseImmutable(nb_point, ListStart, ListListAdj),
            std::move(ListPoint)};
  }
}

TotalArrayInterpolation
INTERPOL_ConstructTotalArray(std::vector<TotalArrGetData> const &ListTotalArr,
                             std::vector<int> const &ListSpongeSize,
                             std::vector<int> const &ListFatherGrid,
                             GridArray const &GrdArrOut,
                             bool const &AllowExtrapolation) {
  TotalArrayInterpolation TotalArrInt;
  TotalArrInt.ARVD = GrdArrOut.ARVD;
  int nbGrid = ListTotalArr.size();
  if (nbGrid == 1) {
    bool testEq = TestEqualityGridArray(GrdArrOut, ListTotalArr[0].GrdArr);
    if (testEq) {
      TotalArrInt.NeedInterp = false;
      TotalArrInt.ListTotalArr = ListTotalArr;
      return TotalArrInt;
    }
  }
  int eta_rho = GrdArrOut.GrdArrRho.LON.rows();
  int xi_rho = GrdArrOut.GrdArrRho.LON.cols();
  MyMatrix<uint8_t> MSKatt;
  MSKatt.setZero(eta_rho, xi_rho);
  std::vector<SingleArrayInterpolationGen> ListSingleArrayInterpolationGen(
      nbGrid);
  std::vector<MyMatrix<int>> ListInsideMask(nbGrid);
  std::vector<MyMatrix<double>> ListHatFunction1(nbGrid);
  GraphSparseImmutable eGR = GetGraphSparseVertexAdjacency(GrdArrOut).first;
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    ListSingleArrayInterpolationGen[iGrid] =
        INTERPOL_CreateSingleRecVarInterpolGen(
            GrdArrOut, ListTotalArr[iGrid].GrdArr, AllowExtrapolation);
    MyMatrix<uint8_t> MSKinside =
        ComputeInsideMask(ListSingleArrayInterpolationGen[iGrid].e_arr);
    MSKatt += MSKinside;
    if (ListSpongeSize[iGrid] <= 0) {
      std::cerr << "The value of ListSpongeSize should be non-negative\n";
      throw TerminalException{1};
    }
    ListHatFunction1[iGrid] =
        HatFunctionFromMask(MSKinside, GrdArrOut, eGR, ListSpongeSize[iGrid]);
  }
  //  std::cerr << "HatFunctions have been computed\n";
  MyMatrix<int> MSK;
  MSK.setZero(eta_rho, xi_rho);
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      if (MSKatt(i, j) > 0)
        MSK(i, j) = 1;
  std::vector<std::vector<int>> ListChildren(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    int eFather = ListFatherGrid[iGrid];
    if (eFather != -1)
      ListChildren[eFather].push_back(iGrid);
  }
  //  std::cerr << "ListChildren has been computed\n";
  MyMatrix<double> TotalSumHat = ZeroMatrix<double>(eta_rho, xi_rho);
  std::vector<MyMatrix<double>> ListHatFunction2(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    MyMatrix<double> TheHatSma = ListHatFunction1[iGrid];
    for (auto &eGrid : ListChildren[iGrid])
      for (int i = 0; i < eta_rho; i++)
        for (int j = 0; j < xi_rho; j++)
          TheHatSma(i, j) = TheHatSma(i, j) * (static_cast<double>(1) -
                                               ListHatFunction1[eGrid](i, j));
    ListHatFunction2[iGrid] = TheHatSma;
    TotalSumHat += TheHatSma;
  }
  //  std::cerr << "TotalSumHat and ListHatFunction2 have been computed\n";
  std::vector<MyMatrix<double>> ListHatFunction3(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    MyMatrix<double> TheHatSma = ListHatFunction2[iGrid];
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        if (TotalSumHat(i, j) > 0) {
          TheHatSma(i, j) = TheHatSma(i, j) / TotalSumHat(i, j);
        } else {
          TheHatSma(i, j) = 0;
        }
      }
    ListHatFunction3[iGrid] = TheHatSma;
  }
  //  std::cerr << "ListHatFunction3 have been computed\n";

  //
  // Computing information for the climatological reading of data.
  //
  std::vector<double> ListMinTime(nbGrid), ListMaxTime(nbGrid);
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    double BeginTime = MinimumTimeHistoryArray(ListTotalArr[iGrid].eArr);
    double EndTime = MaximumTimeHistoryArray(ListTotalArr[iGrid].eArr);
    ListMinTime[iGrid] = BeginTime;
    ListMaxTime[iGrid] = EndTime;
  }
  double StartTime = VectorMax(ListMinTime);
  double EndTime = VectorMin(ListMaxTime);

  TotalArrInt.NeedInterp = true;
  TotalArrInt.eta_rho = eta_rho;
  TotalArrInt.xi_rho = xi_rho;
  TotalArrInt.MSK = MSK;
  TotalArrInt.nbGrid = nbGrid;
  TotalArrInt.ListHatFunction = ListHatFunction3;
  TotalArrInt.ListSingleArrayInterpolationGen = ListSingleArrayInterpolationGen;
  TotalArrInt.ListTotalArr = ListTotalArr;
  TotalArrInt.StartTime = StartTime;
  TotalArrInt.EndTime = EndTime;
  return TotalArrInt;
}

RecVar INTERPOL_MultipleRecVarInterpolation(
    TotalArrayInterpolation const &TotalArrInt, GridArray const &GrdArrOut,
    std::string const &eVarName, double const &eTimeDay) {
  if (!TotalArrInt.NeedInterp)
    return ModelSpecificVarSpecificTime(TotalArrInt.ListTotalArr[0], eVarName,
                                        eTimeDay);
  RecVar eRecVar;
  int nbGrid = TotalArrInt.nbGrid;
  int eta_rho = TotalArrInt.eta_rho;
  int xi_rho = TotalArrInt.xi_rho;
  TotalArrGetData TotalArrTrivial;
  TotalArrTrivial.GrdArr.ModelName = "TRIVIAL";
  RecSymbolic RecS =
      ModelSpecificVarSpecificTime(TotalArrTrivial, eVarName, eTimeDay).RecS;
  std::string const &VarNature = RecS.VarNature;
  //
  // allocating needed variables.
  //
  int Nvert;
  if (TotalArrInt.ARVD.IsAssigned)
    Nvert = TotalArrInt.ARVD.N;
  MyMatrix<double> F, U, V;
  Eigen::Tensor<double, 3> Uthree, Vthree, Tens3;
  if (VarNature == "rho" || VarNature == "uv") {
    F = ZeroMatrix<double>(eta_rho, xi_rho);
  }
  if (VarNature == "uv") {
    U = ZeroMatrix<double>(eta_rho, xi_rho);
    V = ZeroMatrix<double>(eta_rho, xi_rho);
  }
  if (VarNature == "3Drho" || VarNature == "3Duv") {
    Tens3 = ZeroTensor3<double>(Nvert, eta_rho, xi_rho);
  }
  if (VarNature == "3Duv") {
    Uthree = ZeroTensor3<double>(Nvert, eta_rho, xi_rho);
    Vthree = ZeroTensor3<double>(Nvert, eta_rho, xi_rho);
  }
  MyMatrix<double> Unity = ZeroMatrix<double>(eta_rho, xi_rho);
  bool PrintDebug = true;
  auto PrintRecVarInfo = [&](RecVar const &uRecVar, GridArray const &GrdArr,
                             std::string const &strO) -> void {
    if (PrintDebug) {
      std::cerr << "  Begin debugging --------- " << strO
                << " ------------------------------\n";
      if (VarNature == "rho") {
        PrintMMA_FCT(uRecVar.F, GrdArr.GrdArrRho.MSK, "F", "unitblk");
      }
      if (VarNature == "uv") {
        PrintMMA_FCT(uRecVar.F, GrdArr.GrdArrRho.MSK, "F", "unitblk");
        PrintMMA_FCT(uRecVar.U, GrdArr.GrdArrRho.MSK, "U", "uniblk");
        PrintMMA_FCT(uRecVar.V, GrdArr.GrdArrRho.MSK, "V", "unitblk");
      }
      if (VarNature == "3Drho") {
        for (int k = 0; k < Nvert; k++) {
          MyMatrix<double> Fred = DimensionExtraction(uRecVar.Tens3, 0, k);
          std::string strF = "Fred k" + std::to_string(k);
          PrintMMA_FCT(Fred, GrdArr.GrdArrRho.MSK, strF, "unitblk");
        }
      }
      std::cerr << "  End debugging ---------- " << strO
                << " -----------------------------\n";
    }
  };

  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    MyMatrix<double> eHatFunction = TotalArrInt.ListHatFunction[iGrid];
    std::cerr << "   eVarName=" << eVarName << "\n";
    RecVar fRecVar = ModelSpecificVarSpecificTime(
        TotalArrInt.ListTotalArr[iGrid], eVarName, eTimeDay);
    PrintRecVarInfo(fRecVar, TotalArrInt.ListTotalArr[iGrid].GrdArr, "fRecVar");
    RecVar gRecVar = INTERPOL_SingleRecVarInterpolation(
        TotalArrInt.ListSingleArrayInterpolationGen[iGrid], fRecVar);
    PrintRecVarInfo(gRecVar, GrdArrOut, "gRecVar");

    Unity += eHatFunction;
    if (VarNature == "rho") {
      F += eHatFunction.cwiseProduct(gRecVar.F);
    }
    if (VarNature == "uv") {
      U += eHatFunction.cwiseProduct(gRecVar.U);
      V += eHatFunction.cwiseProduct(gRecVar.V);
      F += eHatFunction.cwiseProduct(gRecVar.F);
    }
    if (VarNature == "3Drho") {
      for (int i = 0; i < eta_rho; i++)
        for (int j = 0; j < xi_rho; j++)
          for (int k = 0; k < Nvert; k++)
            Tens3(k, i, j) += gRecVar.Tens3(k, i, j) * eHatFunction(i, j);
    }
    if (VarNature == "3Duv") {
      for (int i = 0; i < eta_rho; i++)
        for (int j = 0; j < xi_rho; j++)
          for (int k = 0; k < Nvert; k++) {
            Uthree(k, i, j) += gRecVar.Uthree(k, i, j) * eHatFunction(i, j);
            Vthree(k, i, j) += gRecVar.Vthree(k, i, j) * eHatFunction(i, j);
            Tens3(k, i, j) += gRecVar.Tens3(k, i, j) * eHatFunction(i, j);
          }
    }
  }
  double TotalErr = 0;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      if (TotalArrInt.MSK(i, j) == 1)
        TotalErr += fabs(Unity(i, j) - static_cast<double>(1));
    }
  if (TotalErr > 1)
    std::cerr << "TotalErr=" << TotalErr << "\n";
  eRecVar.RecS = RecS;
  //
  // assigning arrays
  //
  if (VarNature == "rho") {
    eRecVar.F = F;
  }
  if (VarNature == "uv") {
    eRecVar.U = U;
    eRecVar.V = V;
    eRecVar.F = F;
  }
  if (VarNature == "3Drho") {
    eRecVar.Tens3 = Tens3;
  }
  if (VarNature == "3Duv") {
    eRecVar.Uthree = Uthree;
    eRecVar.Vthree = Vthree;
    eRecVar.Tens3 = Tens3;
  }
  PrintRecVarInfo(eRecVar, GrdArrOut, "eRecVar");
  return eRecVar;
}

std::string ConvertTypename(std::string const &name1,
                            std::string const &name2) {
  std::vector<std::string> LStr = STRING_Split(name2, ":");
  int siz = LStr.size();
  if (siz == 1)
    return name1;
  if (siz > 2) {
    std::cerr << "|LStr|=" << siz << "\n";
    std::cerr << "Allowed sizes are 1 or 2\n";
    std::cerr << "name1 = " << name1 << "\n";
    std::cerr << "name2 = " << name2 << "\n";
    throw TerminalException{1};
  }
  return name1 + "_" + LStr[1];
}

RecTime INTERPOL_NetcdfInitialize(std::string const &eFileNC,
                                  GridArray const &GrdArr,
                                  std::vector<std::string> const &ListVarName) {
  if (!FILE_IsFileMakeable(eFileNC)) {
    std::cerr << "We cannot create the file eFileNC=" << eFileNC << "\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  netCDF::NcDim dateStrDim = dataFile.addDim("dateString", 19);
  dataFile.putAtt("Conventions", "CF-1.4");
  RecTime eRec = AddTimeArray(dataFile, "ocean_time", static_cast<double>(0));

  std::vector<std::string> LDim;
  if (GrdArr.IsFE) {
    int nbWet = GrdArr.GrdArrRho.LON.rows();
    int nbEle = GrdArr.INE.rows();
    netCDF::NcDim eDim = dataFile.addDim("nbNode", nbWet);
    netCDF::NcDim eDim3 = dataFile.addDim("three", 3);
    netCDF::NcDim eDim1 = dataFile.addDim("one", 1);
    netCDF::NcDim eDimEle = dataFile.addDim("nbEle", nbEle);
    std::vector<std::string> LDimNode{"nbNode"};
    netCDF::NcVar eVAR_lon = dataFile.addVar("lon", "double", LDimNode);
    eVAR_lon.putAtt("long_name", std::string("longitude"));
    eVAR_lon.putAtt("units", std::string("degrees_east"));
    eVAR_lon.putAtt("standard_name", std::string("longitude"));
    netCDF::NcVar eVAR_lat = dataFile.addVar("lat", "double", LDimNode);
    eVAR_lat.putAtt("long_name", std::string("latitude"));
    eVAR_lat.putAtt("units", std::string("degrees_north"));
    eVAR_lat.putAtt("standard_name", std::string("latitude"));
    netCDF::NcVar eVAR_dep = dataFile.addVar("depth", "double", LDimNode);
    std::vector<std::string> LDimOne{"one"};
    netCDF::NcVar eVAR_lsph = dataFile.addVar("LSPHE", "int", LDimOne);
    std::vector<std::string> LDimINE{"nbEle", "three"};
    netCDF::NcVar eVAR_ine = dataFile.addVar("ele", "int", LDimINE);
    //
    std::vector<double> A(nbWet);
    for (int i = 0; i < nbWet; i++)
      A[i] = GrdArr.GrdArrRho.LON(i, 0);
    eVAR_lon.putVar(A.data());
    for (int i = 0; i < nbWet; i++)
      A[i] = GrdArr.GrdArrRho.LAT(i, 0);
    eVAR_lat.putVar(A.data());
    const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
    for (int i = 0; i < nbWet; i++)
      A[i] = DEP(i, 0);
    eVAR_dep.putVar(A.data());
    //
    int eSphe = 1;
    eVAR_lsph.putVar(&eSphe);
    //
    std::vector<int> ine(3 * nbEle);
    int idx = 0;
    for (int ie = 0; ie < nbEle; ie++)
      for (int i = 0; i < 3; i++) {
        ine[idx] = GrdArr.INE(ie, i) + 1;
        idx++;
      }
    eVAR_ine.putVar(ine.data());
    //
    LDim = {"ocean_time", "nbNode"};
  } else {
    int eta_rho = GrdArr.GrdArrRho.LON.rows();
    int xi_rho = GrdArr.GrdArrRho.LON.cols();
    netCDF::NcDim eDim1 = dataFile.addDim("eta_rho", eta_rho);
    netCDF::NcDim eDim2 = dataFile.addDim("xi_rho", xi_rho);
    auto InsertVar = [&](std::string const &name,
                         MyMatrix<double> const &eVar) -> void {
      std::vector<double> A(eta_rho * xi_rho);
      int idx = 0;
      for (int i = 0; i < eta_rho; i++)
        for (int j = 0; j < xi_rho; j++) {
          A[idx] = eVar(i, j);
          idx++;
        }
      std::vector<std::string> LDim{"eta_rho", "xi_rho"};
      CFinformation eRecCF = GetCFnames(name);
      netCDF::NcVar VAR = dataFile.addVar(name, "double", LDim);
      VAR.putAtt("long_name", eRecCF.LongName);
      VAR.putAtt("units", eRecCF.Units);
      VAR.putAtt("standard_name", eRecCF.StdName);
      VAR.putVar(A.data());
    };
    InsertVar("lon", GrdArr.GrdArrRho.LON);
    InsertVar("lat", GrdArr.GrdArrRho.LAT);
    InsertVar("ang", GrdArr.GrdArrRho.ANG);
    MyMatrix<double> MSK_double(eta_rho, xi_rho);
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        MSK_double(i, j) = static_cast<double>(GrdArr.GrdArrRho.MSK(i, j));
    InsertVar("mask", MSK_double);
    LDim = {"ocean_time", "eta_rho", "xi_rho"};
  }
  //
  TotalArrGetData TotalArr;
  TotalArr.GrdArr.ModelName = "TRIVIAL";
  double eTimeDay = 0;
  std::cerr << "|ListVarName|=" << ListVarName.size() << "\n";
  for (auto &eVarName : ListVarName) {
    std::cerr << "Init for variable eVarName=" << eVarName << "\n";
    RecVar eRecVar = ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
    if (eRecVar.RecS.VarNature == "rho") {
      netCDF::NcVar eVAR_rho =
          dataFile.addVar(eRecVar.RecS.VarName1, "float", LDim);
      eVAR_rho.putAtt("long_name", eRecVar.RecS.VarName2);
      eVAR_rho.putAtt("units", eRecVar.RecS.Unit);
      eVAR_rho.putAtt("coordinates", "lon lat");
    } else {
      std::string nameU = ConvertTypename(eRecVar.RecS.nameU, eVarName);
      std::string nameV = ConvertTypename(eRecVar.RecS.nameV, eVarName);
      netCDF::NcVar eVAR_u = dataFile.addVar(nameU, "float", LDim);
      netCDF::NcVar eVAR_v = dataFile.addVar(nameV, "float", LDim);
      eVAR_u.putAtt("long_name", eRecVar.RecS.VarName2);
      eVAR_v.putAtt("long_name", eRecVar.RecS.VarName2);
      eVAR_u.putAtt("units", eRecVar.RecS.Unit);
      eVAR_v.putAtt("units", eRecVar.RecS.Unit);
      eVAR_u.putAtt("coordinates", "lon lat");
      eVAR_v.putAtt("coordinates", "lon lat");
    }
  }
  return eRec;
}

std::string ROMS_Surface_NetcdfInitialize_SingleVar(
    netCDF::NcFile &dataFile, bool const &ConstantDefinition,
    bool const &PutLonLatAngArray, bool const &IsRegrid,
    GridArray const &GrdArr, std::string const &eVarName) {
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  TotalArrGetData TotalArr;
  TotalArr.GrdArr.ModelName = "TRIVIAL";
  double eTimeDay = 0;
  RecVar eRecVar = ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
  //
  double RefTimeROMS = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  std::string strTime_ROMS = eRecVar.RecS.strTime_ROMS;
  if (strTime_ROMS == "unset") {
    std::cerr << "We require the variable to be one of the standard ROMS "
                 "supported ones.\n";
    std::cerr << "Those are WIND10, cloud, rain, swrad, lwrad, AIRT2, Rh2\n";
    std::cerr << "eVarName=" << eVarName << "\n";
    throw TerminalException{1};
  }
  if (ConstantDefinition) {
    netCDF::NcDim dateStrDim = dataFile.addDim("dateString", 19);
    netCDF::NcDim eDim1 = dataFile.addDim("eta_rho", eta_rho);
    netCDF::NcDim eDim2 = dataFile.addDim("xi_rho", xi_rho);
  }
  AddTimeArrayROMS(dataFile, strTime_ROMS, RefTimeROMS);
  auto InsertVar = [&](std::string const &name,
                       MyMatrix<double> const &eVar) -> void {
    std::vector<std::string> LDim{"eta_rho", "xi_rho"};
    CFinformation eRecCF = GetCFnames(name);
    netCDF::NcVar VAR = dataFile.addVar(name, "double", LDim);
    if (VAR.isNull()) {
      std::cerr << "We have VAR that is NULL\n";
      throw TerminalException{1};
    }
    VAR.putAtt("long_name", eRecCF.LongName);
    VAR.putAtt("units", eRecCF.Units);
    VAR.putAtt("standard_name", eRecCF.StdName);
    std::vector<double> A(eta_rho * xi_rho);
    int idx = 0;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A[idx] = eVar(i, j);
        idx++;
      }
    VAR.putVar(A.data());
  };
  if (PutLonLatAngArray && ConstantDefinition) {
    InsertVar("lon", GrdArr.GrdArrRho.LON);
    InsertVar("lat", GrdArr.GrdArrRho.LAT);
    InsertVar("ang", GrdArr.GrdArrRho.ANG);
    MyMatrix<double> MSK_double(eta_rho, xi_rho);
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        MSK_double(i, j) = static_cast<double>(GrdArr.GrdArrRho.MSK(i, j));
    InsertVar("mask", MSK_double);
  }
  std::vector<std::string> LDimVar{strTime_ROMS, "eta_rho", "xi_rho"};
  //
  if (eRecVar.RecS.VarNature == "rho") {
    if (!eRecVar.RecS.varName_ROMS) {
      std::cerr << "varName_ROMS has not been assigned 1\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
      std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &varname = *eRecVar.RecS.varName_ROMS;
    netCDF::NcVar eVAR_rho = dataFile.addVar(varname, "float", LDimVar);
    eVAR_rho.putAtt("long_name", eRecVar.RecS.VarName2);
    eVAR_rho.putAtt("units", eRecVar.RecS.Unit);
    if (IsRegrid)
      eVAR_rho.putAtt("coordinates", "lon lat");
  } else {
    if (!eRecVar.RecS.varName_ROMS_U || !eRecVar.RecS.varName_ROMS_V) {
      std::cerr << "varName_ROMS_U of varName_ROMS_V has not been assigned\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1
                << "  VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &varname_U = *eRecVar.RecS.varName_ROMS_U;
    std::string const &varname_V = *eRecVar.RecS.varName_ROMS_V;
    std::string nameU = ConvertTypename(varname_U, eVarName);
    std::string nameV = ConvertTypename(varname_V, eVarName);
    netCDF::NcVar eVAR_u = dataFile.addVar(nameU, "float", LDimVar);
    netCDF::NcVar eVAR_v = dataFile.addVar(nameV, "float", LDimVar);
    eVAR_u.putAtt("long_name", eRecVar.RecS.VarName2);
    eVAR_v.putAtt("long_name", eRecVar.RecS.VarName2);
    eVAR_u.putAtt("units", eRecVar.RecS.Unit);
    eVAR_v.putAtt("units", eRecVar.RecS.Unit);
    if (IsRegrid) {
      eVAR_u.putAtt("coordinates", "lon lat");
      eVAR_v.putAtt("coordinates", "lon lat");
    }
  }
  return strTime_ROMS;
}

struct ROMS_NC_VarInfo {
  std::string strTime;
  std::string eFileNC;
  std::string eVarName;
};

std::vector<ROMS_NC_VarInfo>
ROMS_Surface_NetcdfInitialize(std::string const &eFileNC, bool const &IsRegrid,
                              bool const &SingleFile, GridArray const &GrdArr,
                              std::vector<std::string> const &ListVarName) {
  bool IsFirst = true;
  std::vector<std::string> ListFileNC;
  for (auto &eVarName : ListVarName) {
    std::string eFileNC_real;
    if (SingleFile) {
      eFileNC_real = eFileNC;
    } else {
      int len = eFileNC.size();
      std::string LastChar = eFileNC.substr(len - 3, 3);
      if (LastChar != ".nc") {
        std::cerr << "Last 3 characters of eFileNC must be .nc\n";
        std::cerr << "LastChar=" << LastChar << "\n";
        throw TerminalException{1};
      }
      eFileNC_real = eFileNC.substr(0, len - 3) + "_" + eVarName + ".nc";
    }
    ListFileNC.push_back(eFileNC_real);
    if (IsFirst || !SingleFile) {
      if (!FILE_IsFileMakeable(eFileNC_real)) {
        std::cerr << "Request to create file eFileNC_real=" << eFileNC_real
                  << "\n";
        std::cerr << "but the directory does not exist\n";
        throw TerminalException{1};
      }
      netCDF::NcFile dataFile(eFileNC_real, netCDF::NcFile::replace,
                              netCDF::NcFile::nc4);
    }
    IsFirst = false;
  }
  bool PutLonLatAngArray = IsRegrid;
  int nbVar = ListVarName.size();
  //
  std::vector<ROMS_NC_VarInfo> ListArrROMS(nbVar);
  bool ConstantDefinition = true;
  for (int iVar = 0; iVar < nbVar; iVar++) {
    std::string eFileNC_real = ListFileNC[iVar];
    std::string eVarName = ListVarName[iVar];
    netCDF::NcFile dataFile(eFileNC_real, netCDF::NcFile::write);
    std::string strTime = ROMS_Surface_NetcdfInitialize_SingleVar(
        dataFile, ConstantDefinition, PutLonLatAngArray, IsRegrid, GrdArr,
        eVarName);
    ListArrROMS[iVar] = {strTime, eFileNC_real, eVarName};
    if (SingleFile)
      ConstantDefinition = false;
  }
  return ListArrROMS;
}

std::vector<RecVar>
GetListArrayTracerTrivial(std::vector<std::string> const &ListVarName) {
  std::vector<RecVar> ListArrayTracer;
  for (auto &eVarName : ListVarName) {
    int posVect = PositionVect(
        {"ZetaOcean", "Temp", "Salt", "Curr", "CurrBaro"}, eVarName);
    if (posVect == -1) {
      RecVar eRecVar = RetrieveTrivialRecVar(eVarName);
      ListArrayTracer.push_back(eRecVar);
    }
  }
  return ListArrayTracer;
}

struct ROMSstate {
  double eTimeDay;
  MyMatrix<double> ZETA;
  Eigen::Tensor<double, 3> Temp;
  Eigen::Tensor<double, 3> Salt;
  MyMatrix<double> Ubar;
  MyMatrix<double> Vbar;
  Eigen::Tensor<double, 3> U;
  Eigen::Tensor<double, 3> V;
  std::vector<RecVar> ListAddiTracer;
};

void ROMS_InitialHistory_NetcdfInitialize(
    std::string const &FileOut, GridArray const &GrdArr,
    std::vector<std::string> const &ListAddiVarnameROMS) {
  if (!FILE_IsFileMakeable(FileOut)) {
    std::cerr << "Request to create file FileOut=" << FileOut << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 1\n";
  netCDF::NcFile dataFile(FileOut, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  //  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
  //  netCDF::NcFile::nc4);
  double RefTimeROMS = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int s_rho = GrdArr.ARVD.N;
  int s_w = GrdArr.ARVD.N + 1;
  int eta_u = eta_rho;
  int eta_v = eta_rho - 1;
  int xi_u = xi_rho - 1;
  int xi_v = xi_rho;
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 2\n";
  netCDF::NcDim dateStrDim = dataFile.addDim("dateString", 19);
  netCDF::NcDim eDim_eta_rho = dataFile.addDim("eta_rho", eta_rho);
  netCDF::NcDim eDim_xi_rho = dataFile.addDim("xi_rho", xi_rho);
  netCDF::NcDim eDim_eta_u = dataFile.addDim("eta_u", eta_u);
  netCDF::NcDim eDim_xi_u = dataFile.addDim("xi_u", xi_u);
  netCDF::NcDim eDim_eta_v = dataFile.addDim("eta_v", eta_v);
  netCDF::NcDim eDim_xi_v = dataFile.addDim("xi_v", xi_v);
  netCDF::NcDim eDim_s_rho = dataFile.addDim("s_rho", s_rho);
  netCDF::NcDim eDim_s_w = dataFile.addDim("s_w", s_w);
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 3\n";
  //
  (void)AddTimeArray(dataFile, "ocean_time", RefTimeROMS);
  std::string strOceanTime = "ocean_time";
  std::string strEtaRho = "eta_rho";
  std::string strEtaU = "eta_u";
  std::string strEtaV = "eta_v";
  std::string strXiRho = "xi_rho";
  std::string strXiU = "xi_u";
  std::string strXiV = "xi_v";
  std::string strSRho = "s_rho";
  std::string strSW = "s_w";
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 4\n";
  netCDF::NcVar eVAR_zeta =
      dataFile.addVar("zeta", "float", {strOceanTime, strEtaRho, strXiRho});
  eVAR_zeta.putAtt("long_name", "free-surface");
  eVAR_zeta.putAtt("units", "meter");
  eVAR_zeta.putAtt("time", "ocean_time");
  eVAR_zeta.putAtt("grid", "grid");
  eVAR_zeta.putAtt("location", "face");
  eVAR_zeta.putAtt("coordinates", "lon_rho lat_rho ocean_time");
  eVAR_zeta.putAtt("field", "free-surface, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 5\n";
  netCDF::NcVar eVAR_temp = dataFile.addVar(
      "temp", "float", {strOceanTime, strSRho, strEtaRho, strXiRho});
  eVAR_temp.putAtt("long_name", "potential temperature");
  eVAR_temp.putAtt("units", "Celsius");
  eVAR_temp.putAtt("time", "ocean_time");
  eVAR_temp.putAtt("grid", "grid");
  eVAR_temp.putAtt("location", "face");
  eVAR_temp.putAtt("coordinates", "lon_rho lat_rho s_rho ocean_time");
  eVAR_temp.putAtt("field", "temperature, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 6\n";
  netCDF::NcVar eVAR_salt = dataFile.addVar(
      "salt", "float", {strOceanTime, strSRho, strEtaRho, strXiRho});
  eVAR_salt.putAtt("long_name", "salinity");
  eVAR_salt.putAtt("time", "ocean_time");
  eVAR_salt.putAtt("grid", "grid");
  eVAR_salt.putAtt("location", "face");
  eVAR_salt.putAtt("coordinates", "lon_rho lat_rho s_rho ocean_time");
  eVAR_salt.putAtt("field", "salinity, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 7\n";
  netCDF::NcVar eVAR_ubar =
      dataFile.addVar("ubar", "float", {strOceanTime, strEtaU, strXiU});
  eVAR_ubar.putAtt("long_name", "vertically integrated u-momentum component");
  eVAR_ubar.putAtt("units", "meter second-1");
  eVAR_ubar.putAtt("time", "ocean_time");
  eVAR_ubar.putAtt("grid", "grid");
  eVAR_ubar.putAtt("location", "edge1");
  eVAR_ubar.putAtt("coordinates", "lon_u lat_u ocean_time");
  eVAR_ubar.putAtt("field", "ubar-velocity, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 8\n";
  netCDF::NcVar eVAR_vbar =
      dataFile.addVar("vbar", "float", {strOceanTime, strEtaV, strXiV});
  eVAR_vbar.putAtt("long_name", "vertically integrated v-momentum component");
  eVAR_vbar.putAtt("units", "meter second-1");
  eVAR_vbar.putAtt("time", "ocean_time");
  eVAR_vbar.putAtt("grid", "grid");
  eVAR_vbar.putAtt("location", "edge2");
  eVAR_vbar.putAtt("coordinates", "lon_v lat_v ocean_time");
  eVAR_vbar.putAtt("field", "vbar-velocity, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 9\n";
  netCDF::NcVar eVAR_u =
      dataFile.addVar("u", "float", {strOceanTime, strSRho, strEtaU, strXiU});
  eVAR_u.putAtt("long_name", "u-momentum component");
  eVAR_u.putAtt("units", "meter second-1");
  eVAR_u.putAtt("time", "ocean_time");
  eVAR_u.putAtt("grid", "grid");
  eVAR_u.putAtt("location", "edge1");
  eVAR_u.putAtt("coordinates", "lon_u lat_u s_rho ocean_time");
  eVAR_u.putAtt("field", "u-velocity, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 10\n";
  netCDF::NcVar eVAR_v =
      dataFile.addVar("v", "float", {strOceanTime, strSRho, strEtaV, strXiV});
  eVAR_v.putAtt("long_name", "v-momentum component");
  eVAR_v.putAtt("units", "meter second-1");
  eVAR_v.putAtt("time", "ocean_time");
  eVAR_v.putAtt("grid", "grid");
  eVAR_v.putAtt("location", "edge2");
  eVAR_v.putAtt("coordinates", "lon_v lat_v s_rho ocean_time");
  eVAR_v.putAtt("field", "v-velocity, scalar, series");
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 11\n";
  for (auto &VarNameRoms : ListAddiVarnameROMS) {
    std::cerr << "VarNameRoms=" << VarNameRoms << "\n";
    netCDF::NcVar eVAR_tracer = dataFile.addVar(
        VarNameRoms, "float", {strOceanTime, strSRho, strEtaRho, strXiRho});
  }
  //
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 12\n";
  WriteROMSverticalStratification(dataFile, GrdArr.ARVD);
  std::cerr << "ROMS_InitialHistory_NetcdfInitialize, step 13\n";
}

// This is for
// * ROMS surface forcing
// * ROMS boundary forcing
void ROMS_WRITE_TIME(netCDF::NcFile &dataFile, std::string const &strTimeName,
                     int const &pos, double const &eTimeDay) {
  std::string strTimeNameDay = strTimeName;
  std::string strTimeNameStr = strTimeName + "_str";
  netCDF::NcVar timeVarDay = dataFile.getVar(strTimeNameDay);
  if (timeVarDay.isNull()) {
    std::cerr << "strTimeNameDay = " << strTimeNameDay << "\n";
    std::cerr << "timeVarDay is null\n";
    throw TerminalException{1};
  }
  netCDF::NcVar timeVarStr = dataFile.getVar(strTimeNameStr);
  if (timeVarStr.isNull()) {
    std::cerr << "strTimeNameStr = " << strTimeNameStr << "\n";
    std::cerr << "timeVarStr is null\n";
    throw TerminalException{1};
  }
  double RefTimeROMS = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
  std::vector<size_t> start2{size_t(pos)};
  std::vector<size_t> count2{1};
  double eTimeWrite = eTimeDay - RefTimeROMS;
  timeVarDay.putVar(start2, count2, &eTimeWrite);
  std::vector<size_t> start3{size_t(pos), 0};
  std::vector<size_t> count3{1, 19};
  timeVarStr.putVar(start3, count3, strPres.c_str());
}

// This is for
// * ROMS history file
// * ROMS initial file
void ROMS_WRITE_TIME_HISTORY_INITIAL(netCDF::NcFile &dataFile,
                                     std::string const &strTimeName,
                                     int const &pos, double const &eTimeDay) {
  std::string strTimeNameSec = strTimeName;
  std::string strTimeNameDay = strTimeName + "_day";
  std::string strTimeNameStr = strTimeName + "_str";
  netCDF::NcVar timeVarSec = dataFile.getVar(strTimeNameSec);
  if (timeVarSec.isNull()) {
    std::cerr << "strTimeNameSec = " << strTimeNameSec << "\n";
    std::cerr << "timeVarDay is null\n";
    throw TerminalException{1};
  }
  netCDF::NcVar timeVarDay = dataFile.getVar(strTimeNameDay);
  if (timeVarDay.isNull()) {
    std::cerr << "strTimeNameDay = " << strTimeNameDay << "\n";
    std::cerr << "timeVarDay is null\n";
    throw TerminalException{1};
  }
  netCDF::NcVar timeVarStr = dataFile.getVar(strTimeNameStr);
  if (timeVarStr.isNull()) {
    std::cerr << "strTimeNameStr = " << strTimeNameStr << "\n";
    std::cerr << "timeVarStr is null\n";
    throw TerminalException{1};
  }
  double RefTimeROMS = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
  std::vector<size_t> start2{size_t(pos)};
  std::vector<size_t> count2{1};
  double eTimeWriteDay = eTimeDay - RefTimeROMS;
  double eTimeWriteSec = eTimeWriteDay * 86400;
  timeVarSec.putVar(start2, count2, &eTimeWriteSec);
  timeVarDay.putVar(start2, count2, &eTimeWriteDay);
  std::vector<size_t> start3{size_t(pos), 0};
  std::vector<size_t> count3{1, 19};
  timeVarStr.putVar(start3, count3, strPres.c_str());
}

void ROMS_InitialHistory_NetcdfAppend(std::string const &FileOut,
                                      ROMSstate const &eState,
                                      GridArray const &GrdArr,
                                      size_t const &pos) {
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int s_rho = GrdArr.ARVD.N;
  int eta_u = eta_rho;
  int eta_v = eta_rho - 1;
  int xi_u = xi_rho - 1;
  int xi_v = xi_rho;
  netCDF::NcFile dataFile(FileOut, netCDF::NcFile::write, netCDF::NcFile::nc4);
  netCDF::NcVar eVAR_zeta = dataFile.getVar("zeta");
  netCDF::NcVar eVAR_temp = dataFile.getVar("temp");
  netCDF::NcVar eVAR_salt = dataFile.getVar("salt");
  netCDF::NcVar eVAR_ubar = dataFile.getVar("ubar");
  netCDF::NcVar eVAR_vbar = dataFile.getVar("vbar");
  netCDF::NcVar eVAR_u = dataFile.getVar("u");
  netCDF::NcVar eVAR_v = dataFile.getVar("v");
  //
  ROMS_WRITE_TIME_HISTORY_INITIAL(dataFile, "ocean_time", pos, eState.eTimeDay);
  std::vector<size_t> start, count;
  int idx;
  //
  std::vector<float> A(eta_rho * xi_rho);
  start = {pos, 0, 0};
  count = {1, size_t(eta_rho), size_t(xi_rho)};
  idx = 0;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++) {
      A[idx] = static_cast<float>(eState.ZETA(i, j));
      idx++;
    }
  netCDF::NcVar eVar1 = dataFile.getVar("zeta");
  eVAR_zeta.putVar(start, count, A.data());
  //
  // The tracers
  //
  std::vector<float> Atr(s_rho * eta_rho * xi_rho);
  start = {pos, 0, 0, 0};
  count = {1, size_t(s_rho), size_t(eta_rho), size_t(xi_rho)};
  idx = 0;
  for (int i = 0; i < s_rho; i++)
    for (int j = 0; j < eta_rho; j++)
      for (int k = 0; k < xi_rho; k++) {
        Atr[idx] = static_cast<float>(eState.Temp(i, j, k));
        idx++;
      }
  netCDF::NcVar eVar2 = dataFile.getVar("temp");
  eVAR_temp.putVar(start, count, Atr.data());
  //
  idx = 0;
  for (int i = 0; i < s_rho; i++)
    for (int j = 0; j < eta_rho; j++)
      for (int k = 0; k < xi_rho; k++) {
        Atr[idx] = static_cast<float>(eState.Salt(i, j, k));
        idx++;
      }
  netCDF::NcVar eVar3 = dataFile.getVar("salt");
  eVAR_salt.putVar(start, count, Atr.data());
  //
  for (auto &eRecVar : eState.ListAddiTracer) {
    if (!eRecVar.RecS.varName_ROMS) {
      std::cerr << "varName_ROMS has not been assigned 2\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
      std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &strNameROMS = *eRecVar.RecS.varName_ROMS;
    netCDF::NcVar eVAR_tracer = dataFile.getVar(strNameROMS);
    //
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_rho; j++)
        for (int k = 0; k < xi_rho; k++) {
          Atr[idx] = static_cast<float>(eRecVar.Tens3(i, j, k));
          idx++;
        }
    eVAR_tracer.putVar(start, count, Atr.data());
  }
  //
  std::vector<float> Au(s_rho * eta_u * xi_u);
  start = {pos, 0, 0, 0};
  count = {1, size_t(s_rho), size_t(eta_u), size_t(xi_u)};
  idx = 0;
  for (int i = 0; i < s_rho; i++)
    for (int j = 0; j < eta_u; j++)
      for (int k = 0; k < xi_u; k++) {
        Au[idx] = static_cast<float>(eState.U(i, j, k));
        idx++;
      }
  netCDF::NcVar eVar4 = dataFile.getVar("u");
  eVAR_u.putVar(start, count, Au.data());
  //
  std::vector<float> Av(s_rho * eta_v * xi_v);
  start = {pos, 0, 0, 0};
  count = {1, size_t(s_rho), size_t(eta_v), size_t(xi_v)};
  idx = 0;
  for (int i = 0; i < s_rho; i++)
    for (int j = 0; j < eta_v; j++)
      for (int k = 0; k < xi_v; k++) {
        Av[idx] = static_cast<float>(eState.V(i, j, k));
        idx++;
      }
  netCDF::NcVar eVar5 = dataFile.getVar("v");
  eVAR_v.putVar(start, count, Av.data());
  //
  std::vector<float> Aubar(eta_u * xi_u);
  start = {pos, 0, 0};
  count = {1, size_t(eta_u), size_t(xi_u)};
  idx = 0;
  for (int i = 0; i < eta_u; i++)
    for (int j = 0; j < xi_u; j++) {
      Aubar[idx] = static_cast<float>(eState.Ubar(i, j));
      idx++;
    }
  netCDF::NcVar eVar6 = dataFile.getVar("ubar");
  eVAR_ubar.putVar(start, count, Aubar.data());
  //
  std::vector<float> Avbar(eta_v * xi_v);
  start = {pos, 0, 0};
  count = {1, size_t(eta_v), size_t(xi_v)};
  idx = 0;
  for (int i = 0; i < eta_v; i++)
    for (int j = 0; j < xi_v; j++) {
      Avbar[idx] = static_cast<float>(eState.Vbar(i, j));
      idx++;
    }
  netCDF::NcVar eVar7 = dataFile.getVar("vbar");
  eVAR_vbar.putVar(start, count, Avbar.data());
}

void ROMS_BOUND_NetcdfInitialize(std::string const &eFileNC,
                                 GridArray const &GrdArr,
                                 std::vector<std::string> const &ListSides,
                                 std::vector<RecVar> const &ListArrayTracer) {
  int posSouth = PositionVect(ListSides, std::string("South"));
  int posNorth = PositionVect(ListSides, std::string("North"));
  int posWest = PositionVect(ListSides, std::string("West"));
  int posEast = PositionVect(ListSides, std::string("East"));
  for (auto &eStr : ListSides) {
    if (eStr != "South" && eStr != "North" && eStr != "West" &&
        eStr != "East") {
      std::cerr << "Possible allowed values are North, South, West and East\n";
      std::cerr << "eStr=" << eStr << "\n";
      throw TerminalException{1};
    }
  }
  if (ListSides.size() == 0) {
    std::cerr << "|ListSides| = 0\n";
    std::cerr << "If calling for netcdf boundary conditions\n";
    std::cerr << "you need to select in ListSides what you want\n";
    std::cerr
        << "e.g. ListSides = \"South\", \"North\", \"West\", \"East\",    \n";
    throw TerminalException{1};
  }
  double RefTimeROMS = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int s_rho = GrdArr.ARVD.N;
  int s_w = GrdArr.ARVD.N + 1;
  int eta_u = eta_rho;
  int eta_v = eta_rho - 1;
  int xi_u = xi_rho - 1;
  int xi_v = xi_rho;
  if (!FILE_IsFileMakeable(eFileNC)) {
    std::cerr << "Request to create file eFileNC=" << eFileNC << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  netCDF::NcDim dateStrDim = dataFile.addDim("dateString", 19);
  netCDF::NcDim eDim_eta_rho = dataFile.addDim("eta_rho", eta_rho);
  netCDF::NcDim eDim_xi_rho = dataFile.addDim("xi_rho", xi_rho);
  netCDF::NcDim eDim_eta_u = dataFile.addDim("eta_u", eta_u);
  netCDF::NcDim eDim_xi_u = dataFile.addDim("xi_u", xi_u);
  netCDF::NcDim eDim_eta_v = dataFile.addDim("eta_v", eta_v);
  netCDF::NcDim eDim_xi_v = dataFile.addDim("xi_v", xi_v);
  netCDF::NcDim eDim_s_rho = dataFile.addDim("s_rho", s_rho);
  netCDF::NcDim eDim_s_w = dataFile.addDim("s_w", s_w);
  //
  AddTimeArrayRomsBound(dataFile, "zeta_time", RefTimeROMS);
  AddTimeArrayRomsBound(dataFile, "temp_time", RefTimeROMS);
  AddTimeArrayRomsBound(dataFile, "salt_time", RefTimeROMS);
  AddTimeArrayRomsBound(dataFile, "v2d_time", RefTimeROMS);
  AddTimeArrayRomsBound(dataFile, "v3d_time", RefTimeROMS);
  std::string strZetaTime = "zeta_time";
  std::string strTempTime = "temp_time";
  std::string strSaltTime = "salt_time";
  std::string strV2dTime = "v2d_time";
  std::string strV3dTime = "v3d_time";
  std::string strEtaRho = "eta_rho";
  std::string strEtaU = "eta_u";
  std::string strEtaV = "eta_v";
  std::string strXiRho = "xi_rho";
  std::string strXiU = "xi_u";
  std::string strXiV = "xi_v";
  std::string strSRho = "s_rho";
  std::string strSW = "s_w";
  if (posEast != -1) {
    std::cerr << "Writing EAST from Lon/LAT = ("
              << GrdArr.GrdArrRho.LON(0, xi_rho - 1) << "/"
              << GrdArr.GrdArrRho.LAT(0, xi_rho - 1) << ") to ("
              << GrdArr.GrdArrRho.LON(eta_rho - 1, xi_rho - 1) << "/"
              << GrdArr.GrdArrRho.LAT(eta_rho - 1, xi_rho - 1) << ")\n";
    netCDF::NcVar eVAR1 =
        dataFile.addVar("zeta_east", "float", {strZetaTime, strEtaRho});
    netCDF::NcVar eVAR2 = dataFile.addVar("temp_east", "float",
                                          {strTempTime, strSRho, strEtaRho});
    netCDF::NcVar eVAR3 = dataFile.addVar("salt_east", "float",
                                          {strSaltTime, strSRho, strEtaRho});
    netCDF::NcVar eVAR4 =
        dataFile.addVar("ubar_east", "float", {strV2dTime, strEtaU});
    netCDF::NcVar eVAR5 =
        dataFile.addVar("vbar_east", "float", {strV2dTime, strEtaV});
    netCDF::NcVar eVAR6 =
        dataFile.addVar("u_east", "float", {strV3dTime, strSRho, strEtaU});
    netCDF::NcVar eVAR7 =
        dataFile.addVar("v_east", "float", {strV3dTime, strSRho, strEtaV});
  }
  if (posWest != -1) {
    std::cerr << "Writing EAST from Lon/LAT = (" << GrdArr.GrdArrRho.LON(0, 0)
              << "/" << GrdArr.GrdArrRho.LAT(0, 0) << ") to ("
              << GrdArr.GrdArrRho.LON(eta_rho - 1, 0) << "/"
              << GrdArr.GrdArrRho.LAT(eta_rho - 1, 0) << ")\n";
    netCDF::NcVar eVAR1 =
        dataFile.addVar("zeta_west", "float", {strZetaTime, strEtaRho});
    netCDF::NcVar eVAR2 = dataFile.addVar("temp_west", "float",
                                          {strTempTime, strSRho, strEtaRho});
    netCDF::NcVar eVAR3 = dataFile.addVar("salt_west", "float",
                                          {strSaltTime, strSRho, strEtaRho});
    netCDF::NcVar eVAR4 =
        dataFile.addVar("ubar_west", "float", {strV2dTime, strEtaU});
    netCDF::NcVar eVAR5 =
        dataFile.addVar("vbar_west", "float", {strV2dTime, strEtaV});
    netCDF::NcVar eVAR6 =
        dataFile.addVar("u_west", "float", {strV3dTime, strSRho, strEtaU});
    netCDF::NcVar eVAR7 =
        dataFile.addVar("v_west", "float", {strV3dTime, strSRho, strEtaV});
  }
  if (posNorth != -1) {
    std::cerr << "Writing EAST from Lon/LAT = ("
              << GrdArr.GrdArrRho.LON(eta_rho - 1, 0) << "/"
              << GrdArr.GrdArrRho.LAT(eta_rho - 1, 0) << ") to ("
              << GrdArr.GrdArrRho.LON(eta_rho - 1, xi_rho - 1) << "/"
              << GrdArr.GrdArrRho.LAT(eta_rho - 1, xi_rho - 1) << ")\n";
    netCDF::NcVar eVAR1 =
        dataFile.addVar("zeta_north", "float", {strZetaTime, strXiRho});
    netCDF::NcVar eVAR2 = dataFile.addVar("temp_north", "float",
                                          {strTempTime, strSRho, strXiRho});
    netCDF::NcVar eVAR3 = dataFile.addVar("salt_north", "float",
                                          {strSaltTime, strSRho, strXiRho});
    netCDF::NcVar eVAR4 =
        dataFile.addVar("ubar_north", "float", {strV2dTime, strXiU});
    netCDF::NcVar eVAR5 =
        dataFile.addVar("vbar_north", "float", {strV2dTime, strXiV});
    netCDF::NcVar eVAR6 =
        dataFile.addVar("u_north", "float", {strV3dTime, strSRho, strXiU});
    netCDF::NcVar eVAR7 =
        dataFile.addVar("v_north", "float", {strV3dTime, strSRho, strXiV});
  }
  if (posSouth != -1) {
    std::cerr << "Writing EAST from Lon/LAT = (" << GrdArr.GrdArrRho.LON(0, 0)
              << "/" << GrdArr.GrdArrRho.LAT(0, 0) << ") to ("
              << GrdArr.GrdArrRho.LON(0, xi_rho - 1) << "/"
              << GrdArr.GrdArrRho.LAT(0, xi_rho - 1) << ")\n";
    netCDF::NcVar eVAR1 =
        dataFile.addVar("zeta_south", "float", {strZetaTime, strXiRho});
    netCDF::NcVar eVAR2 = dataFile.addVar("temp_south", "float",
                                          {strTempTime, strSRho, strXiRho});
    netCDF::NcVar eVAR3 = dataFile.addVar("salt_south", "float",
                                          {strSaltTime, strSRho, strXiRho});
    netCDF::NcVar eVAR4 =
        dataFile.addVar("ubar_south", "float", {strV2dTime, strXiU});
    netCDF::NcVar eVAR5 =
        dataFile.addVar("vbar_south", "float", {strV2dTime, strXiV});
    netCDF::NcVar eVAR6 =
        dataFile.addVar("u_south", "float", {strV3dTime, strSRho, strXiU});
    netCDF::NcVar eVAR7 =
        dataFile.addVar("v_south", "float", {strV3dTime, strSRho, strXiV});
  }
  //
  // Writing the ARVD vertical discretization
  //
  std::vector<std::string> ListDimEmpty;
  netCDF::NcVar eVAR_1 = dataFile.addVar("Vtransform", "int", ListDimEmpty);
  eVAR_1.putVar(&GrdArr.ARVD.Vtransform);
  netCDF::NcVar eVAR_2 = dataFile.addVar("Vstretching", "int", ListDimEmpty);
  eVAR_2.putVar(&GrdArr.ARVD.Vstretching);
  netCDF::NcVar eVAR_3 = dataFile.addVar("theta_s", "double", ListDimEmpty);
  eVAR_3.putVar(&GrdArr.ARVD.theta_s);
  netCDF::NcVar eVAR_4 = dataFile.addVar("theta_b", "double", ListDimEmpty);
  eVAR_4.putVar(&GrdArr.ARVD.theta_b);
  netCDF::NcVar eVAR_5 = dataFile.addVar("Tcline", "double", ListDimEmpty);
  eVAR_5.putVar(&GrdArr.ARVD.Tcline);
  netCDF::NcVar eVAR_6 = dataFile.addVar("hc", "double", ListDimEmpty);
  eVAR_6.putVar(&GrdArr.ARVD.hc);
  //
  // Now the variable
  //
  netCDF::NcVar eVAR_7 = dataFile.addVar("Cs_r", "double", {strSRho});
  eVAR_7.putVar(GrdArr.ARVD.Cs_r.data());
  netCDF::NcVar eVAR_8 = dataFile.addVar("Cs_w", "double", {strSW});
  eVAR_8.putVar(GrdArr.ARVD.Cs_w.data());
  netCDF::NcVar eVAR_9 = dataFile.addVar("s_rho", "double", {strSRho});
  eVAR_9.putVar(GrdArr.ARVD.sc_r.data());
  netCDF::NcVar eVAR_10 = dataFile.addVar("s_w", "double", {strSW});
  eVAR_10.putVar(GrdArr.ARVD.sc_w.data());
  //
  // Now the additional tracers on output
  //
  std::cerr << "|ListArrayTracer|=" << ListArrayTracer.size() << "\n";
  for (auto &eRecVar : ListArrayTracer) {
    if (!eRecVar.RecS.varName_ROMS) {
      std::cerr << "varName_ROMS has not been assigned 3\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
      std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &strNameROMS = *eRecVar.RecS.varName_ROMS;
    if (posEast != -1) {
      std::string str2 = strNameROMS + "_east";
      netCDF::NcVar eVAR2 =
          dataFile.addVar(str2, "float", {strTempTime, strSRho, strEtaRho});
    }
    if (posWest != -1) {
      std::string str2 = strNameROMS + "_west";
      netCDF::NcVar eVAR2 =
          dataFile.addVar(str2, "float", {strTempTime, strSRho, strEtaRho});
    }
    if (posNorth != -1) {
      std::string str2 = strNameROMS + "_north";
      netCDF::NcVar eVAR2 =
          dataFile.addVar(str2, "float", {strTempTime, strSRho, strXiRho});
    }
    if (posSouth != -1) {
      std::string str2 = strNameROMS + "_south";
      netCDF::NcVar eVAR2 =
          dataFile.addVar(str2, "float", {strTempTime, strSRho, strXiRho});
    }
  }
}

void ROMS_BOUND_NetcdfAppend(std::string const &eFileNC,
                             ROMSstate const &eState,
                             std::vector<std::string> const &ListSides,
                             int const &pos) {
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 1\n";
  int posSouth = PositionVect(ListSides, std::string("South"));
  int posNorth = PositionVect(ListSides, std::string("North"));
  int posWest = PositionVect(ListSides, std::string("West"));
  int posEast = PositionVect(ListSides, std::string("East"));
  int eta_rho = eState.ZETA.rows();
  int xi_rho = eState.ZETA.cols();
  int eta_u = eState.Ubar.rows();
  int xi_u = eState.Ubar.cols();
  int eta_v = eState.Vbar.rows();
  int xi_v = eState.Vbar.cols();
  auto LDim = eState.Temp.dimensions();
  int s_rho = LDim[0];
  std::vector<size_t> start, count;
  int idx;
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::write);
  ROMS_WRITE_TIME(dataFile, "zeta_time", pos, eState.eTimeDay);
  ROMS_WRITE_TIME(dataFile, "temp_time", pos, eState.eTimeDay);
  ROMS_WRITE_TIME(dataFile, "salt_time", pos, eState.eTimeDay);
  ROMS_WRITE_TIME(dataFile, "v3d_time", pos, eState.eTimeDay);
  ROMS_WRITE_TIME(dataFile, "v2d_time", pos, eState.eTimeDay);
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 2\n";
  if (posSouth != -1) {
    std::vector<float> A1(xi_rho);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_rho)};
    for (int i = 0; i < xi_rho; i++)
      A1[i] = static_cast<float>(eState.ZETA(0, i));
    netCDF::NcVar eVar1 = dataFile.getVar("zeta_south");
    eVar1.putVar(start, count, A1.data());
    //
    std::vector<float> A2(s_rho * xi_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A2[idx] = static_cast<float>(eState.Temp(i, 0, j));
        idx++;
      }
    netCDF::NcVar eVar2 = dataFile.getVar("temp_south");
    eVar2.putVar(start, count, A2.data());
    //
    std::vector<float> A3(s_rho * xi_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A3[idx] = static_cast<float>(eState.Salt(i, 0, j));
        idx++;
      }
    netCDF::NcVar eVar3 = dataFile.getVar("salt_south");
    eVar3.putVar(start, count, A3.data());
    //
    std::vector<float> A4(s_rho * xi_u);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_u)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_u; j++) {
        A4[idx] = static_cast<float>(eState.U(i, 0, j));
        idx++;
      }
    netCDF::NcVar eVar4 = dataFile.getVar("u_south");
    eVar4.putVar(start, count, A4.data());
    //
    std::vector<float> A5(s_rho * xi_v);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_v)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_v; j++) {
        A5[idx] = static_cast<float>(eState.V(i, 0, j));
        idx++;
      }
    netCDF::NcVar eVar5 = dataFile.getVar("v_south");
    eVar5.putVar(start, count, A5.data());
    //
    std::vector<float> A6(xi_u);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_u)};
    idx = 0;
    for (int j = 0; j < xi_u; j++) {
      A6[idx] = static_cast<float>(eState.Ubar(0, j));
      idx++;
    }
    netCDF::NcVar eVar6 = dataFile.getVar("ubar_south");
    eVar6.putVar(start, count, A6.data());
    //
    std::vector<float> A7(xi_v);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_v)};
    idx = 0;
    for (int j = 0; j < xi_v; j++) {
      A7[idx] = static_cast<float>(eState.Vbar(0, j));
      idx++;
    }
    netCDF::NcVar eVar7 = dataFile.getVar("vbar_south");
    eVar7.putVar(start, count, A7.data());
    //
  }
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 3\n";
  if (posNorth != -1) {
    std::vector<float> A1(xi_rho);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_rho)};
    for (int i = 0; i < xi_rho; i++)
      A1[i] = static_cast<float>(eState.ZETA(eta_rho - 1, i));
    netCDF::NcVar eVar1 = dataFile.getVar("zeta_north");
    eVar1.putVar(start, count, A1.data());
    //
    std::vector<float> A2(s_rho * xi_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A2[idx] = static_cast<float>(eState.Temp(i, eta_rho - 1, j));
        idx++;
      }
    netCDF::NcVar eVar2 = dataFile.getVar("temp_north");
    eVar2.putVar(start, count, A2.data());
    //
    std::vector<float> A3(s_rho * xi_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A3[idx] = static_cast<float>(eState.Salt(i, eta_rho - 1, j));
        idx++;
      }
    netCDF::NcVar eVar3 = dataFile.getVar("salt_north");
    eVar3.putVar(start, count, A3.data());
    //
    std::vector<float> A4(s_rho * xi_u);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_u)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_u; j++) {
        A4[idx] = static_cast<float>(eState.U(i, eta_u - 1, j));
        idx++;
      }
    netCDF::NcVar eVar4 = dataFile.getVar("u_north");
    eVar4.putVar(start, count, A4.data());
    //
    std::vector<float> A5(s_rho * xi_v);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(xi_v)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < xi_v; j++) {
        A5[idx] = static_cast<float>(eState.V(i, eta_v - 1, j));
        idx++;
      }
    netCDF::NcVar eVar5 = dataFile.getVar("v_north");
    eVar5.putVar(start, count, A5.data());
    //
    std::vector<float> A6(xi_u);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_u)};
    idx = 0;
    for (int j = 0; j < xi_u; j++) {
      A6[idx] = static_cast<float>(eState.Ubar(eta_u - 1, j));
      idx++;
    }
    netCDF::NcVar eVar6 = dataFile.getVar("ubar_north");
    eVar6.putVar(start, count, A6.data());
    //
    std::vector<float> A7(xi_v);
    start = {size_t(pos), 0};
    count = {1, size_t(xi_v)};
    idx = 0;
    for (int j = 0; j < xi_v; j++) {
      A7[idx] = static_cast<float>(eState.Vbar(eta_v - 1, j));
      idx++;
    }
    netCDF::NcVar eVar7 = dataFile.getVar("vbar_north");
    eVar7.putVar(start, count, A7.data());
    //
  }
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 4\n";
  if (posEast != -1) {
    //    std::cerr << "Doing posEast\n";
    std::vector<float> A1(eta_rho);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_rho)};
    for (int i = 0; i < eta_rho; i++)
      A1[i] = static_cast<float>(eState.ZETA(i, xi_rho - 1));
    /*
    std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
    std::cerr << "ZETA : A1(min/max)=" << VectorMin(A1) << " / " <<
    VectorMax(A1) << "\n"; std::cerr << "eState.ZETA : (min/max)=" <<
    eState.ZETA.minCoeff() << " / " << eState.ZETA.maxCoeff() << "\n"; for (int
    j=0; j<xi_rho; j++) { double sumAbsZeta = 0; for (int i=0; i<eta_rho; i++)
        sumAbsZeta += T_abs(eState.ZETA(i,j));
      std::cerr << " j=" << j << " |zeta|=" << sumAbsZeta << "\n";
      }*/
    netCDF::NcVar eVar1 = dataFile.getVar("zeta_east");
    eVar1.putVar(start, count, A1.data());
    //
    std::vector<float> A2(s_rho * eta_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_rho; j++) {
        A2[idx] = static_cast<float>(eState.Temp(i, j, xi_rho - 1));
        idx++;
      }
    netCDF::NcVar eVar2 = dataFile.getVar("temp_east");
    eVar2.putVar(start, count, A2.data());
    //
    std::vector<float> A3(s_rho * eta_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_rho; j++) {
        A3[idx] = static_cast<float>(eState.Salt(i, j, xi_rho - 1));
        idx++;
      }
    netCDF::NcVar eVar3 = dataFile.getVar("salt_east");
    eVar3.putVar(start, count, A3.data());
    //
    std::vector<float> A4(s_rho * eta_u);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_u)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_u; j++) {
        A4[idx] = static_cast<float>(eState.U(i, j, xi_u - 1));
        idx++;
      }
    netCDF::NcVar eVar4 = dataFile.getVar("u_east");
    eVar4.putVar(start, count, A4.data());
    //
    std::vector<float> A5(s_rho * eta_v);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_v)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_v; j++) {
        A5[idx] = static_cast<float>(eState.V(i, j, xi_v - 1));
        idx++;
      }
    netCDF::NcVar eVar5 = dataFile.getVar("v_east");
    eVar5.putVar(start, count, A5.data());
    //
    std::vector<float> A6(eta_u);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_u)};
    idx = 0;
    for (int j = 0; j < eta_u; j++) {
      A6[idx] = static_cast<float>(eState.Ubar(j, xi_u - 1));
      idx++;
    }
    netCDF::NcVar eVar6 = dataFile.getVar("ubar_east");
    eVar6.putVar(start, count, A6.data());
    //
    std::vector<float> A7(eta_v);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_v)};
    idx = 0;
    for (int j = 0; j < eta_v; j++) {
      A7[idx] = static_cast<float>(eState.Vbar(j, xi_v - 1));
      idx++;
    }
    netCDF::NcVar eVar7 = dataFile.getVar("vbar_east");
    eVar7.putVar(start, count, A7.data());
    //
  }
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 5\n";
  if (posWest != -1) {
    std::vector<float> A1(eta_rho);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_rho)};
    for (int i = 0; i < eta_rho; i++)
      A1[i] = static_cast<float>(eState.ZETA(i, 0));
    netCDF::NcVar eVar1 = dataFile.getVar("zeta_west");
    eVar1.putVar(start, count, A1.data());
    //
    std::vector<float> A2(s_rho * eta_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_rho; j++) {
        A2[idx] = static_cast<float>(eState.Temp(i, j, 0));
        idx++;
      }
    netCDF::NcVar eVar2 = dataFile.getVar("temp_west");
    eVar2.putVar(start, count, A2.data());
    //
    std::vector<float> A3(s_rho * eta_rho);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_rho)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_rho; j++) {
        A3[idx] = static_cast<float>(eState.Salt(i, j, 0));
        idx++;
      }
    netCDF::NcVar eVar3 = dataFile.getVar("salt_west");
    eVar3.putVar(start, count, A3.data());
    //
    std::vector<float> A4(s_rho * eta_u);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_u)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_u; j++) {
        A4[idx] = static_cast<float>(eState.U(i, j, 0));
        idx++;
      }
    netCDF::NcVar eVar4 = dataFile.getVar("u_west");
    eVar4.putVar(start, count, A4.data());
    //
    std::vector<float> A5(s_rho * eta_v);
    start = {size_t(pos), 0, 0};
    count = {1, size_t(s_rho), size_t(eta_v)};
    idx = 0;
    for (int i = 0; i < s_rho; i++)
      for (int j = 0; j < eta_v; j++) {
        A5[idx] = static_cast<float>(eState.V(i, j, 0));
        idx++;
      }
    netCDF::NcVar eVar5 = dataFile.getVar("v_west");
    eVar5.putVar(start, count, A5.data());
    //
    std::vector<float> A6(eta_u);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_u)};
    idx = 0;
    for (int j = 0; j < eta_u; j++) {
      A6[idx] = static_cast<float>(eState.Ubar(j, 0));
      idx++;
    }
    netCDF::NcVar eVar6 = dataFile.getVar("ubar_west");
    eVar6.putVar(start, count, A6.data());
    //
    std::vector<float> A7(eta_v);
    start = {size_t(pos), 0};
    count = {1, size_t(eta_v)};
    idx = 0;
    for (int j = 0; j < eta_v; j++) {
      A7[idx] = static_cast<float>(eState.Vbar(j, 0));
      idx++;
    }
    netCDF::NcVar eVar7 = dataFile.getVar("vbar_west");
    eVar7.putVar(start, count, A7.data());
    //
  }
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 6\n";
  //
  // Additional tracers
  //
  for (auto &eRecVar : eState.ListAddiTracer) {
    if (!eRecVar.RecS.varName_ROMS) {
      std::cerr << "varName_ROMS has not been assigned 4\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
      std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &strNameROMS = *eRecVar.RecS.varName_ROMS;
    if (posEast != -1) {
      std::string str2 = strNameROMS + "_east";
      std::vector<float> A(s_rho * eta_rho);
      start = {size_t(pos), 0, 0};
      count = {1, size_t(s_rho), size_t(eta_rho)};
      idx = 0;
      for (int i = 0; i < s_rho; i++)
        for (int j = 0; j < eta_rho; j++) {
          A[idx] = static_cast<float>(eRecVar.Tens3(i, j, xi_rho - 1));
          idx++;
        }
      netCDF::NcVar eVar2 = dataFile.getVar(str2);
      eVar2.putVar(start, count, A.data());
    }
    if (posWest != -1) {
      std::string str2 = strNameROMS + "_west";
      std::vector<float> A(s_rho * eta_rho);
      start = {size_t(pos), 0, 0};
      count = {1, size_t(s_rho), size_t(eta_rho)};
      idx = 0;
      for (int i = 0; i < s_rho; i++)
        for (int j = 0; j < eta_rho; j++) {
          A[idx] = static_cast<float>(eRecVar.Tens3(i, j, 0));
          idx++;
        }
      netCDF::NcVar eVar2 = dataFile.getVar(str2);
      eVar2.putVar(start, count, A.data());
    }
    if (posNorth != -1) {
      std::string str2 = strNameROMS + "_north";
      std::vector<float> A(s_rho * xi_rho);
      start = {size_t(pos), 0, 0};
      count = {1, size_t(s_rho), size_t(xi_rho)};
      idx = 0;
      for (int i = 0; i < s_rho; i++)
        for (int j = 0; j < xi_rho; j++) {
          A[idx] = static_cast<float>(eRecVar.Tens3(i, eta_rho - 1, j));
          idx++;
        }
      netCDF::NcVar eVar2 = dataFile.getVar(str2);
      eVar2.putVar(start, count, A.data());
    }
    if (posSouth != -1) {
      std::string str2 = strNameROMS + "_south";
      std::vector<float> A(s_rho * xi_rho);
      start = {size_t(pos), 0, 0};
      count = {1, size_t(s_rho), size_t(xi_rho)};
      idx = 0;
      for (int i = 0; i < s_rho; i++)
        for (int j = 0; j < xi_rho; j++) {
          A[idx] = static_cast<float>(eRecVar.Tens3(i, 0, j));
          idx++;
        }
      netCDF::NcVar eVar2 = dataFile.getVar(str2);
      eVar2.putVar(start, count, A.data());
    }
  }
  std::cerr << "ROMS_BOUND_NetcdfAppend, step 7\n";
}

Eigen::Tensor<double, 3> ZeroThreeTensor(int const &dim0, int const &dim1,
                                         int const &dim2) {
  Eigen::Tensor<double, 3> Tens3(dim0, dim1, dim2);
  for (int i0 = 0; i0 < dim0; i0++)
    for (int i1 = 0; i1 < dim1; i1++)
      for (int i2 = 0; i2 < dim2; i2++)
        Tens3(i0, i1, i2) = 0;
  return Tens3;
}

Eigen::Tensor<double, 3> ConstantThreeTensor(int const &dim0, int const &dim1,
                                             int const &dim2,
                                             double const &val) {
  Eigen::Tensor<double, 3> Tens3(dim0, dim1, dim2);
  for (int i0 = 0; i0 < dim0; i0++)
    for (int i1 = 0; i1 < dim1; i1++)
      for (int i2 = 0; i2 < dim2; i2++)
        Tens3(i0, i1, i2) = val;
  return Tens3;
}

ROMSstate GetRomsStateFromVariables(GridArray const &GrdArr,
                                    std::vector<RecVar> const &ListRecVar) {
  std::cerr << "GetRomsStateFromVariables, step 1\n";
  bool HasZeta = false, HasTemp = false, HasSalt = false, HasCurr = false,
       HasCurrBaro = false;
  ROMSstate eState;
  Eigen::Tensor<double, 3> Ufield, Vfield;
  std::vector<RecVar> ListAddiTracer;
  std::cerr << "GetRomsStateFromVariables, step 2\n";
  MyMatrix<double> UbarMod, VbarMod;
  for (auto &eRecVar : ListRecVar) {
    bool IsMatch = false;
    std::string VarName1 = eRecVar.RecS.VarName1;
    //    std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
    if (VarName1 == "ZetaOcean") {
      eState.eTimeDay = eRecVar.RecS.eTimeDay;
      eState.ZETA = eRecVar.F;
      HasZeta = true;
      IsMatch = true;
    }
    if (VarName1 == "Temp") {
      eState.eTimeDay = eRecVar.RecS.eTimeDay;
      eState.Temp = eRecVar.Tens3;
      HasTemp = true;
      IsMatch = true;
    }
    if (VarName1 == "Salt") {
      eState.eTimeDay = eRecVar.RecS.eTimeDay;
      eState.Salt = eRecVar.Tens3;
      HasSalt = true;
      IsMatch = true;
    }
    if (VarName1 == "Curr") {
      eState.eTimeDay = eRecVar.RecS.eTimeDay;
      Ufield = eRecVar.Uthree;
      Vfield = eRecVar.Vthree;
      HasCurr = true;
      IsMatch = true;
    }
    if (VarName1 == "CurrBaro") {
      eState.eTimeDay = eRecVar.RecS.eTimeDay;
      UbarMod = eRecVar.U;
      VbarMod = eRecVar.V;
      HasCurrBaro = true;
      IsMatch = true;
    }
    if (!IsMatch) {
      std::cerr << "Inserting tracer named " << VarName1 << "\n";
      ListAddiTracer.push_back(eRecVar);
    }
  }
  if (!HasZeta || !HasTemp || !HasSalt || !HasCurr || !HasCurrBaro) {
    std::cerr
        << "For the ROMS boundary forcing, we need Zeta, Temp, Salt and Curr\n";
    std::cerr << "    HasZeta=" << HasZeta << "\n";
    std::cerr << "    HasTemp=" << HasTemp << "\n";
    std::cerr << "    HasSalt=" << HasSalt << "\n";
    std::cerr << "    HasCurr=" << HasCurr << "\n";
    std::cerr << "HasCurrBaro=" << HasCurrBaro << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "eState.eTimeDay=" << eState.eTimeDay << "\n";
  std::cerr << "GetRomsStateFromVariables, step 3\n";
  MyMatrix<double> ANGrotation = -GrdArr.GrdArrRho.ANG;
  AngleRhoRot_3D(Ufield, Vfield, ANGrotation);
  MyMatrix<double> UbarInt =
      ConvertBaroclinic_to_Barotropic(Ufield, eState.ZETA, GrdArr);
  MyMatrix<double> VbarInt =
      ConvertBaroclinic_to_Barotropic(Vfield, eState.ZETA, GrdArr);
  //
  bool DoDebug = true;
  if (DoDebug) {
    MyMatrix<double> Ubar_diff = (UbarInt - UbarMod).cwiseAbs();
    MyMatrix<double> Vbar_diff = (VbarInt - VbarMod).cwiseAbs();
    double Li_UbarInt = UbarInt.cwiseAbs().maxCoeff();
    double Li_VbarInt = VbarInt.cwiseAbs().maxCoeff();
    double Li_UbarMod = UbarMod.cwiseAbs().maxCoeff();
    double Li_VbarMod = VbarMod.cwiseAbs().maxCoeff();
    std::cerr << "Li(UbarInt - UbarMod)=" << Ubar_diff.maxCoeff()
              << " Li(VbarInt - VbarMod)=" << Vbar_diff.maxCoeff() << "\n";
    std::cerr << "L1(UbarInt)=" << Li_UbarInt << " Li(VbarInt)=" << Li_VbarInt
              << "\n";
    std::cerr << "L1(UbarMod)=" << Li_UbarMod << " Li(VbarMod)=" << Li_VbarMod
              << "\n";
  }
  //
  std::string CurrentMethod = "Rescaling";
  if (CurrentMethod == "DirectInterp") {
    eState.U = My_rho2u_3D(GrdArr, Ufield);
    eState.V = My_rho2v_3D(GrdArr, Vfield);
    eState.Ubar = My_rho2u_2D(GrdArr, UbarInt);
    eState.Vbar = My_rho2v_2D(GrdArr, VbarInt);
  }
  if (CurrentMethod == "Rescaling") {
    auto get_norm = [&](double const &dx, double const &dy) -> double {
      return sqrt(dx * dx + dy * dy);
    };
    int N = GrdArr.ARVD.N;
    int eta_rho = GrdArr.GrdArrRho.LON.rows();
    int xi_rho = GrdArr.GrdArrRho.LON.cols();
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
          double norm_uvbar_mod = get_norm(UbarMod(i, j), VbarMod(i, j));
          double norm_uvbar_int = get_norm(UbarInt(i, j), VbarInt(i, j));
          if (norm_uvbar_int > norm_uvbar_mod) {
            double coeff = norm_uvbar_mod / norm_uvbar_int;
            for (int iS = 0; iS < N; iS++) {
              Ufield(iS, i, j) *= coeff;
              Vfield(iS, i, j) *= coeff;
            }
          }
        }
    UbarInt = ConvertBaroclinic_to_Barotropic(Ufield, eState.ZETA, GrdArr);
    VbarInt = ConvertBaroclinic_to_Barotropic(Vfield, eState.ZETA, GrdArr);
    //
    eState.U = My_rho2u_3D(GrdArr, Ufield);
    eState.V = My_rho2v_3D(GrdArr, Vfield);
    eState.Ubar = My_rho2u_2D(GrdArr, UbarInt);
    eState.Vbar = My_rho2v_2D(GrdArr, VbarInt);
  }
  eState.ListAddiTracer = ListAddiTracer;
  std::cerr << "GetRomsStateFromVariables, step 4\n";
  return eState;
}

void ROMS_Surface_NetcdfAppendVarName_SingleVar(netCDF::NcFile &dataFile,
                                                GridArray const &GrdArr,
                                                RecVar const &eRecVar,
                                                ROMS_NC_VarInfo const &eArr,
                                                bool const &PrintMMA) {
  std::string eTime = eArr.strTime;
  //  std::cerr << "eTime=" << eTime << "\n";
  std::multimap<std::string, netCDF::NcDim> MapDims = dataFile.getDims();
  std::multimap<std::string, netCDF::NcDim>::iterator iter = MapDims.begin();
  netCDF::NcDim timeDim;
  bool WeFound = false;
  while (iter != MapDims.end()) {
    //    std::cerr << "iter->first=" << iter->first << "\n";
    if (iter->first == eTime) {
      WeFound = true;
      timeDim = iter->second;
    }
    iter++;
  }
  if (!WeFound) {
    std::cerr << "We fail to find the variable eTime = " << eTime << "\n";
    throw TerminalException{1};
  }
  if (!timeDim.isUnlimited()) {
    std::cerr << "Error the time dimension should be unlimited\n";
    throw TerminalException{1};
  }
  bool IsPressure = false;
  if (eRecVar.RecS.FullVarName == "SurfPres")
    IsPressure = true;
  MyMatrix<double> F_raw;
  if (IsPressure)
    F_raw = eRecVar.F / 100;
  else
    F_raw = eRecVar.F;
  int siz = timeDim.getSize();
  double eTimeDay = eRecVar.RecS.eTimeDay;
  ROMS_WRITE_TIME(dataFile, eArr.strTime, siz, eTimeDay);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  std::vector<float> A(eta_rho * xi_rho);
  std::vector<size_t> start{size_t(siz), 0, 0};
  std::vector<size_t> count{1, size_t(eta_rho), size_t(xi_rho)};
  auto print_avg_max_min = [&](MyMatrix<double> const &F,
                               std::string const &name) -> void {
    double sumVal = 0;
    double maxVal = std::numeric_limits<double>::min();
    double minVal = std::numeric_limits<double>::max();
    size_t nb = 0;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
          double val = F(i, j);
          maxVal = std::max(maxVal, val);
          minVal = std::min(minVal, val);
          sumVal += val;
          nb++;
        }
    double avgVal = sumVal / nb;
    std::cerr << "MMA : " << name << " avg=" << avgVal << " min=" << minVal
              << " max=" << maxVal << "\n";
  };
  if (eRecVar.RecS.VarNature == "rho") {
    if (!eRecVar.RecS.varName_ROMS) {
      std::cerr << "varName_ROMS has not been assigned 5\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
      std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &varName_ROMS = *eRecVar.RecS.varName_ROMS;
    if (PrintMMA) {
      print_avg_max_min(F_raw, varName_ROMS);
    }
    //    std::cerr << "varName_ROMS=" << eRecVar.RecS.varName_ROMS << "\n";
    netCDF::NcVar eVar_F = dataFile.getVar(varName_ROMS);
    //    std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
    //    std::cerr << "|eRecVar.F|=" << eRecVar.F.rows() << " / " <<
    //    eRecVar.F.cols() << "\n";
    int sumWet = 0;
    int sumLand = 0;
    double sumWet_d = 0;
    double sumLand_d = 0;
    bool PrintDebug = false;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
          sumWet++;
          sumWet_d += F_raw(i, j);
        } else {
          if (PrintDebug) {
            sumLand++;
            sumLand_d += F_raw(i, j);
          }
        }
      }
    double avgWet = sumWet_d / static_cast<double>(sumWet);
    if (PrintDebug) {
      double avgLand = sumLand_d / double(sumLand);
      std::cerr << "avgWet = " << avgWet << " avgLand = " << avgLand << "\n";
    }
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++)
        if (GrdArr.GrdArrRho.MSK(i, j) == 0)
          F_raw(i, j) = avgWet;
    int idx = 0;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A[idx] = static_cast<float>(F_raw(i, j));
        idx++;
      }
    eVar_F.putVar(start, count, A.data());
  } else {
    if (!eRecVar.RecS.varName_ROMS_U || !eRecVar.RecS.varName_ROMS_V) {
      std::cerr << "varName_ROMS_U or arName_ROMS_V has not been assigned\n";
      std::cerr << "VarName1=" << eRecVar.RecS.VarName1
                << "  VarName2=" << eRecVar.RecS.VarName2 << "\n";
      throw TerminalException{1};
    }
    std::string const &varName_ROMS_U = *eRecVar.RecS.varName_ROMS_U;
    std::string const &varName_ROMS_V = *eRecVar.RecS.varName_ROMS_V;
    if (PrintMMA) {
      print_avg_max_min(eRecVar.U, varName_ROMS_U);
      print_avg_max_min(eRecVar.V, varName_ROMS_V);
    }
    netCDF::NcVar eVar_U = dataFile.getVar(varName_ROMS_U);
    netCDF::NcVar eVar_V = dataFile.getVar(varName_ROMS_V);
    int idx = 0;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A[idx] = static_cast<float>(eRecVar.U(i, j));
        idx++;
      }
    eVar_U.putVar(start, count, A.data());
    idx = 0;
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        A[idx] = static_cast<float>(eRecVar.V(i, j));
        idx++;
      }
    eVar_V.putVar(start, count, A.data());
  }
}

void ROMS_Surface_NetcdfAppendVarName(GridArray const &GrdArr,
                                      std::vector<RecVar> const &ListRecVar,
                                      std::vector<ROMS_NC_VarInfo> &ListArrROMS,
                                      bool const &PrintMMA) {
  int nbVar = ListRecVar.size();
  for (int iVar = 0; iVar < nbVar; iVar++) {
    netCDF::NcFile dataFile(ListArrROMS[iVar].eFileNC, netCDF::NcFile::write);
    ROMS_Surface_NetcdfAppendVarName_SingleVar(
        dataFile, GrdArr, ListRecVar[iVar], ListArrROMS[iVar], PrintMMA);
  }
}

void INTERPOL_NetcdfAppendVarName(std::string const &eFileNC,
                                  GridArray const &GrdArr,
                                  std::vector<RecVar> const &ListRecVar,
                                  RecTime &eRec) {
  std::cerr << "INTERPOL_NetcdfAppendVarName : eFileNC = " << eFileNC << "\n";
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::write);
  std::cerr << "   we have dataFile\n";
  std::multimap<std::string, netCDF::NcDim> MapDims = dataFile.getDims();
  std::multimap<std::string, netCDF::NcDim>::iterator iter = MapDims.begin();
  netCDF::NcDim timeDim;
  bool timeDimAssigned = false;
  while (iter != MapDims.end()) {
    if (iter->first == "ocean_time") {
      timeDim = iter->second;
      timeDimAssigned = true;
    }
    iter++;
  }
  if (!timeDimAssigned) {
    std::cerr << "Error, timeDim has not been assigned\n";
    throw TerminalException{1};
  }
  if (!timeDim.isUnlimited()) {
    std::cerr << "Error the dimension should be unlimited\n";
    throw TerminalException{1};
  }
  size_t siz = eRec.timeDim.getSize();
  // Putting the time
  if (ListRecVar.size() == 0) {
    std::cerr << "What is the point of writing\n";
    std::cerr << "ZERO variables\n";
    throw TerminalException{1};
  }
  double eTimeDay = ListRecVar[0].RecS.eTimeDay;
  PutTimeDay(eRec, siz, eTimeDay);
  for (auto &eRecVar : ListRecVar) {
    if (GrdArr.IsFE == 1) {
      int nbWet = GrdArr.GrdArrRho.LON.rows();
      std::vector<float> A(nbWet);
      std::vector<size_t> start{siz, 0};
      std::vector<size_t> count{1, size_t(nbWet)};
      if (eRecVar.RecS.VarNature == "rho") {
        netCDF::NcVar eVar_F = dataFile.getVar(eRecVar.RecS.VarName1);
        for (int i = 0; i < nbWet; i++)
          A[i] = static_cast<float>(eRecVar.F(i, 0));
        eVar_F.putVar(start, count, A.data());
      } else {
        std::string nameU =
            ConvertTypename(eRecVar.RecS.nameU, eRecVar.RecS.FullVarName);
        std::string nameV =
            ConvertTypename(eRecVar.RecS.nameV, eRecVar.RecS.FullVarName);
        netCDF::NcVar eVar_U = dataFile.getVar(nameU);
        netCDF::NcVar eVar_V = dataFile.getVar(nameV);
        for (int i = 0; i < nbWet; i++)
          A[i] = static_cast<float>(eRecVar.U(i, 0));
        eVar_U.putVar(start, count, A.data());
        for (int i = 0; i < nbWet; i++)
          A[i] = static_cast<float>(eRecVar.V(i, 0));
        eVar_V.putVar(start, count, A.data());
      }
    } else {
      int eta_rho = GrdArr.GrdArrRho.LON.rows();
      int xi_rho = GrdArr.GrdArrRho.LON.rows();
      std::vector<float> A(eta_rho * xi_rho);
      std::vector<size_t> start{siz, 0, 0};
      std::vector<size_t> count{1, size_t(eta_rho), size_t(xi_rho)};
      if (eRecVar.RecS.VarNature == "rho") {
        netCDF::NcVar eVar_F = dataFile.getVar(eRecVar.RecS.VarName1);
        int idx = 0;
        for (int i = 0; i < eta_rho; i++)
          for (int j = 0; j < xi_rho; j++) {
            A[idx] = static_cast<float>(eRecVar.F(i, j));
            idx++;
          }
        eVar_F.putVar(start, count, A.data());
      } else {
        std::string nameU =
            ConvertTypename(eRecVar.RecS.nameU, eRecVar.RecS.FullVarName);
        std::string nameV =
            ConvertTypename(eRecVar.RecS.nameV, eRecVar.RecS.FullVarName);
        netCDF::NcVar eVar_U = dataFile.getVar(nameU);
        netCDF::NcVar eVar_V = dataFile.getVar(nameV);
        int idx = 0;
        for (int i = 0; i < eta_rho; i++)
          for (int j = 0; j < xi_rho; j++) {
            A[idx] = static_cast<float>(eRecVar.U(i, j));
            idx++;
          }
        eVar_U.putVar(start, count, A.data());
        idx = 0;
        for (int i = 0; i < eta_rho; i++)
          for (int j = 0; j < xi_rho; j++) {
            A[idx] = static_cast<float>(eRecVar.V(i, j));
            idx++;
          }
        eVar_V.putVar(start, count, A.data());
      }
    }
  }
}

void WaveWatch_WriteData_direct(GridArray const &GrdArrOut,
                                std::vector<RecVar> const &ListRecVar,
                                int &WWIII_nbWritten,
                                int const& IsFormatted_inp) {
  if (ListRecVar.size() == 0) {
    std::cerr << "|ListRecVar| = 0 but it should not\n";
    std::cerr << "We should select at least one variable\n";
    throw TerminalException{1};
  }
  int IsFormatted = IsFormatted_inp;
  std::vector<std::string> ListVarName;
  double eTimeDay = 0;
  for (auto &eRecVar : ListRecVar) {
    ListVarName.push_back(eRecVar.RecS.VarName1);
    eTimeDay = eRecVar.RecS.eTimeDay;
  }
  int WWIII_posWind10 = PositionVect(ListVarName, std::string("WIND10"));
  int WWIII_posSurfCurr = PositionVect(ListVarName, std::string("SurfCurr"));
  int WWIII_posZetaOcean = PositionVect(ListVarName, std::string("ZetaOcean"));
  if (WWIII_posWind10 == -1 && WWIII_posSurfCurr == -1 &&
      WWIII_posZetaOcean == -1) {
    std::cerr << "It is erroneous to call the WaveWatch III output\n";
    std::cerr << "if neither WIND10 nor surface current nor free surface are "
                 "not selected\n";
    throw TerminalException{1};
  }
  //
  //
  //
  std::vector<int> tfnvect = DATE_ConvertMjd2tfn(eTimeDay);
  std::vector<int> TFN(2);
  TFN[0] = tfnvect[0];
  TFN[1] = tfnvect[1];
  int nx = GrdArrOut.GrdArrRho.LON.rows();
  int ny = GrdArrOut.GrdArrRho.LON.cols();
  if (WWIII_nbWritten == 0) {
    int UNGTYPE = 3;
    int GTYPE = UNGTYPE;
    std::cout << "Creation of files, step 1\n";
    std::cout << "Creation of files, nx=" << nx << " ny=" << ny << "\n";
    if (WWIII_posWind10 != -1) {
      int ChoiceFile = 1;
      write_wavewatch_header_(&ChoiceFile, &nx, &ny, &GTYPE, &IsFormatted);
    }
    std::cout << "Creation of files, step 2\n";
    if (WWIII_posSurfCurr != -1) {
      int ChoiceFile = 2;
      write_wavewatch_header_(&ChoiceFile, &nx, &ny, &GTYPE, &IsFormatted);
    }
    std::cout << "Creation of files, step 3\n";
    if (WWIII_posZetaOcean != -1) {
      int ChoiceFile = 3;
      write_wavewatch_header_(&ChoiceFile, &nx, &ny, &GTYPE, &IsFormatted);
    }
    std::cout << "Creation of files, step 4\n";
  }
  if (WWIII_posWind10 != -1) {
    std::vector<float> Uvect(nx * ny);
    std::vector<float> Vvect(nx * ny);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++) {
        Uvect[i + nx * j] =
            static_cast<float>(ListRecVar[WWIII_posWind10].U(i, j));
        Vvect[i + nx * j] =
            static_cast<float>(ListRecVar[WWIII_posWind10].V(i, j));
      }
    std::cerr << "Before call to two entry files, wind\n";
    write_wavewatch_entry_two_field_("wind.ww3", TFN.data(), &nx, &ny,
                                     Uvect.data(), Vvect.data(), &IsFormatted);
    std::cerr << " After call to two entry files, wind\n";
  }
  if (WWIII_posSurfCurr != -1) {
    std::vector<float> Uvect(nx * ny);
    std::vector<float> Vvect(nx * ny);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++) {
        Uvect[i + nx * j] =
            static_cast<float>(ListRecVar[WWIII_posSurfCurr].U(i, j));
        Vvect[i + nx * j] =
            static_cast<float>(ListRecVar[WWIII_posSurfCurr].V(i, j));
      }
    std::cerr << "Before call to two entry files, current\n";
    write_wavewatch_entry_two_field_("current.ww3", TFN.data(), &nx, &ny,
                                     Uvect.data(), Vvect.data(), &IsFormatted);
    std::cerr << " After call to two entry files, current\n";
  }
  if (WWIII_posZetaOcean != -1) {
    std::vector<float> Fvect(nx * ny);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        Fvect[i + nx * j] =
            static_cast<float>(ListRecVar[WWIII_posZetaOcean].F(i, j));
    std::cerr << "Before call to two entry files, level\n";
    write_wavewatch_entry_one_field_("level.ww3", TFN.data(), &nx, &ny,
                                     Fvect.data(), &IsFormatted);
    std::cerr << " After call to two entry files, level\n";
  }
  WWIII_nbWritten++;
}

void WaveWatch_WriteData_nc(GridArray const &GrdArrOut,
                            std::vector<RecVar> const &ListRecVar,
                            int &WWIII_nbWritten) {
  if (ListRecVar.size() != 1) {
    std::cerr << "Writing only one entry at a time\n";
    throw TerminalException{1};
  }
  RecVar const& eRecVar = ListRecVar[0];
  double eTimeDay = eRecVar.RecS.eTimeDay;
  std::string const& eVarName = eRecVar.RecS.VarName1;
  std::string eFile = "/irrelevant/file.nc";
  std::vector<std::string> LVar;
  int nx = GrdArrOut.GrdArrRho.LON.rows();
  int ny = GrdArrOut.GrdArrRho.LON.cols();
  if (eVarName == "WIND10") {
    eFile = "wind.nc";
    LVar = {"uwnd", "vwnd"};
  }
  if (eVarName == "SurfCurr") {
    eFile = "curr.nc";
    LVar = {"ucur", "vcur"};
  }
  if (eVarName == "ZetaOcean") {
    eFile = "zeta.nc";
    LVar = {"wlv"};
  }
  std::string strTime = "time";
  double RefTime = DATE_ConvertSix2mjd({1968, 05, 23, 0, 0, 0});
  if (WWIII_nbWritten == 0) {
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::write);
    AddTimeArray(dataFile, strTime, RefTime);
    dataFile.addDim("nx", nx);
    dataFile.addDim("ny", ny);
    for (auto & eVar : LVar) {
      netCDF::NcVar ncvar = dataFile.addVar(eVar, "float", {"time", "nx", "ny"});
    }
  }
  //
  // Now writing
  //
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::write);
  std::vector<size_t> start_var(3), count_var(3);
  start_var[0] = WWIII_nbWritten;
  start_var[1] = 0;
  start_var[2] = 0;
  count_var[0] = 1;
  count_var[1] = nx;
  count_var[2] = ny;
  std::vector<float> FillVector(nx * ny);
  auto write_array=[&](MyMatrix<double> const& M, std::string const& the_var) -> void {
    size_t pos = 0;
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        FillVector[pos] = M(i,j);
        pos++;
      }
    }
    netCDF::NcVar ncvar = dataFile.getVar(the_var);
    ncvar.putVar(start_var, count_var, FillVector.data());
  };
  // Writing of the time
  ROMS_WRITE_TIME_HISTORY_INITIAL(dataFile, strTime, WWIII_nbWritten, eTimeDay);
  // The wind.
  if (eVarName == "WIND10") {
    write_array(eRecVar.U, "uwnd");
    write_array(eRecVar.V, "vwnd");
  }
  // The current
  if (eVarName == "SurfCurr") {
    write_array(eRecVar.U, "ucur");
    write_array(eRecVar.V, "vcur");
  }
  // The water level
  if (eVarName == "ZetaOcean") {
    write_array(eRecVar.F, "wlv");
  }
  WWIII_nbWritten++;
}



void WaveWatch_WriteData(GridArray const &GrdArrOut,
                         std::vector<RecVar> const &ListRecVar,
                         int &WWIII_nbWritten, std::string const& OutFormat) {

  if (OutFormat == "UNFORMATTED") {
    int IsFormatted = 0;
    return WaveWatch_WriteData_direct(GrdArrOut, ListRecVar, WWIII_nbWritten, IsFormatted);
  }
  if (OutFormat == "FORMATTED") {
    int IsFormatted = 1;
    return WaveWatch_WriteData_direct(GrdArrOut, ListRecVar, WWIII_nbWritten, IsFormatted);
  }
  if (OutFormat == "PRNC") {
  }
  std::cerr << "Failed to find a matching entry in WaveWatch_WriteData\n";
  throw TerminalException{1};
}



struct recNetcdfOutput {
  int iFile;
  int nbWritten;
  double eTimeDay;
  RecTime eRec;
};

void INTERPOL_NetcdfOutput(GridArray const &GrdArrOut,
                           std::vector<RecVar> const &ListRecVar,
                           FullNamelist const &eFull, int const &eMult,
                           recNetcdfOutput &recNO) {
  SingleBlock eBlNETCDF_STANDARD = eFull.ListBlock.at("NETCDF_STANDARD");
  bool WriteIFile = eBlNETCDF_STANDARD.ListBoolValues.at("WriteIFile");
  bool WriteDate = eBlNETCDF_STANDARD.ListBoolValues.at("WriteDate");
  std::string HisPrefixOut =
      eBlNETCDF_STANDARD.ListStringValues.at("HisPrefixOut");
  std::vector<std::string> ListVarName;
  double eTimeDay = 0;
  for (auto &eRecVar : ListRecVar) {
    ListVarName.push_back(eRecVar.RecS.VarName1);
    eTimeDay = eRecVar.RecS.eTimeDay;
  }
  if (recNO.nbWritten == 0)
    recNO.eTimeDay = eTimeDay;
  //
  auto GetFileNC = [&]() -> std::string {
    std::string eFileNC = HisPrefixOut;
    if (WriteIFile)
      eFileNC += "_" + StringNumber(recNO.iFile, 4);
    if (WriteDate) {
      std::string strFile = DATE_ConvertMjd2mystringFile(recNO.eTimeDay);
      eFileNC += "_" + strFile;
    }
    eFileNC += ".nc";
    return eFileNC;
  };
  std::string eFileNC = GetFileNC();
  if (recNO.nbWritten == 0)
    recNO.eRec = INTERPOL_NetcdfInitialize(eFileNC, GrdArrOut, ListVarName);
  INTERPOL_NetcdfAppendVarName(eFileNC, GrdArrOut, ListRecVar, recNO.eRec);
  recNO.nbWritten++;
  if (recNO.nbWritten == eMult) {
    recNO.iFile++;
    recNO.nbWritten = 0;
  }
}

struct recGribOutput {
  int nbWritten;
  double eTimeDay;
  double StartDate_mjd;
};

void INTERPOL_GribAppendVarName(std::string const &eFileGrib,
                                GridArray const &GrdArrOut,
                                std::vector<RecVar> const &ListRecVar) {
  grib_multi_handle *mh = NULL;
  mh = grib_multi_handle_new(0);
  if (!mh) {
    std::cerr << "ERROR: Unable to create multi field handle\n";
    throw TerminalException{1};
  }
  grib_handle *h = NULL;
  h = grib_handle_new_from_samples(NULL, "GRIB2");
  char gridType_str[] = "regular_ll";
  size_t gridType_len = strlen(gridType_str);
  GRIB_CHECK(grib_set_string(h, "gridType", gridType_str, &gridType_len), 0);
  //
  RecVar eRecVar = ListRecVar[0];
  //
  double eTimeDay = eRecVar.RecS.eTimeDay;
  std::vector<int> Date = DATE_ConvertMjd2six(eTimeDay);
  long eYear = Date[0];
  long eMonth = Date[1];
  long eDay = Date[2];
  long eHour = Date[3];
  long eMin = Date[4];
  long eSec = Date[5];
  GRIB_CHECK(grib_set_long(h, "year", eYear), 0);
  GRIB_CHECK(grib_set_long(h, "month", eMonth), 0);
  GRIB_CHECK(grib_set_long(h, "day", eDay), 0);
  GRIB_CHECK(grib_set_long(h, "hour", eHour), 0);
  GRIB_CHECK(grib_set_long(h, "minute", eMin), 0);
  GRIB_CHECK(grib_set_long(h, "second", eSec), 0);
  long dataDate = eYear * 10000 + eMonth * 100 + eDay;
  GRIB_CHECK(grib_set_long(h, "dataDate", dataDate), 0);
  long dataTime = eHour * 100 + eMin;
  GRIB_CHECK(grib_set_long(h, "dataTime", dataTime), 0);
  //
  if (!eRecVar.RecS.varName_GRIB) {
    std::cerr << "varName_GRIB has not been assigned\n";
    std::cerr << "VarName1=" << eRecVar.RecS.VarName1
              << "  VarName2=" << eRecVar.RecS.VarName2 << "\n";
    throw TerminalException{1};
  }
  std::string const &vn_GRIB = *eRecVar.RecS.varName_GRIB;
  size_t shortName_len = strlen(vn_GRIB.c_str());
  GRIB_CHECK(grib_set_string(h, "shortName", vn_GRIB.c_str(), &shortName_len),
             0);
  //
  MyMatrix<double> LON = GrdArrOut.GrdArrRho.LON;
  MyMatrix<double> LAT = GrdArrOut.GrdArrRho.LAT;
  int nbRow = LON.rows();
  int nbCol = LON.cols();
  int Ni = nbRow;
  int Nj = nbCol;
  double lonFirst = LON.minCoeff();
  double lonLast = LON.maxCoeff();
  double latFirst = LAT.minCoeff();
  double latLast = LAT.maxCoeff();
  GRIB_CHECK(grib_set_long(h, "Ni", Ni), 0);
  GRIB_CHECK(grib_set_long(h, "Nj", Nj), 0);
  GRIB_CHECK(grib_set_double(h, "longitudeOfFirstGridPoint", lonFirst), 0);
  GRIB_CHECK(grib_set_double(h, "longitudeOfFirstGridPointInDegrees", lonFirst),
             0);
  GRIB_CHECK(grib_set_double(h, "latitudeOfFirstGridPoint", latFirst), 0);
  GRIB_CHECK(grib_set_double(h, "latitudeOfFirstGridPointInDegrees", latFirst),
             0);
  GRIB_CHECK(grib_set_double(h, "longitudeOfLastGridPoint", lonLast), 0);
  GRIB_CHECK(grib_set_double(h, "longitudeOfLastGridPointInDegrees", lonLast),
             0);
  GRIB_CHECK(grib_set_double(h, "latitudeOfLastGridPoint", latLast), 0);
  GRIB_CHECK(grib_set_double(h, "latitudeOfLastGridPointInDegrees", latLast),
             0);
  double *A;
  int siz = Ni * Nj;
  A = (double *)malloc(siz * sizeof(double));
  int idx = 0;
  for (int i = 0; i < Ni; i++)
    for (int j = 0; j < Nj; j++) {
      A[idx] = eRecVar.F(i, j);
      idx++;
    }
  GRIB_CHECK(grib_set_double_array(h, "values", A, siz), 0);
  //
  //
  //
  FILE *of = NULL;
  of = fopen(eFileGrib.c_str(), "w");
  if (!of) {
    std::cerr << "ERROR: unable to open output file " << eFileGrib << "\n";
    throw TerminalException{1};
  }
  grib_multi_handle_write(mh, of);
  const int start_section = 4; /* Grib2 Product Definition Section */
  grib_multi_handle_append(h, start_section, mh);
  grib_multi_handle_write(mh, of);
  fclose(of);
  grib_handle_delete(h);
  grib_multi_handle_delete(mh);
}

void INTERPOL_GribOutput(GridArray const &GrdArrOut,
                         std::vector<RecVar> const &ListRecVar,
                         FullNamelist const &eFull, int const &eMult,
                         recGribOutput &recGO) {
  SingleBlock eBlGRIB_STANDARD = eFull.ListBlock.at("GRIB_STANDARD");
  bool WriteFromStart = eBlGRIB_STANDARD.ListBoolValues.at("WriteFromStart");
  std::string HisPrefixOut =
      eBlGRIB_STANDARD.ListStringValues.at("HisPrefixOut");
  std::vector<std::string> ListVarName;
  std::cerr << "recGO.nbWritten=" << recGO.nbWritten << "\n";
  double eTimeDay = 0;
  for (auto &eRecVar : ListRecVar) {
    ListVarName.push_back(eRecVar.RecS.VarName1);
    eTimeDay = eRecVar.RecS.eTimeDay;
  }
  if (recGO.nbWritten == 0) {
    recGO.StartDate_mjd = eTimeDay;
  }
  //
  auto GetFileGrib = [&]() -> std::string {
    std::string eFileGrib = HisPrefixOut;
    if (WriteFromStart) {
      double deltaTime =
          (eTimeDay - recGO.StartDate_mjd) * static_cast<double>(24);
      int deltaTime_i = static_cast<int>(round(deltaTime));
      int nbDigit = 2;
      if (deltaTime_i >= 100)
        nbDigit = GetNumberDigit(deltaTime_i);
      eFileGrib += DATE_ConvertMjd2dhmz(recGO.StartDate_mjd) + "+" +
                   StringNumber(deltaTime_i, nbDigit);
    } else {
      eFileGrib += DATE_ConvertMjd2dhmz(eTimeDay);
    }
    eFileGrib += ".grb";
    return eFileGrib;
  };
  std::string eFileGrib = GetFileGrib();
  INTERPOL_GribAppendVarName(eFileGrib, GrdArrOut, ListRecVar);
  recGO.nbWritten++;
}

bool HasField(AnalyticalAlgorithm const &AnalField,
              std::string const &eVarName) {
  if (PositionVect(AnalField.ListNameVariables, eVarName) != -1)
    return true;
  return false;
}

RecVar GetRecVarAnalytical(GridArray const &GrdArr, std::string const &eVarName,
                           double const &eTimeDay,
                           AnalyticalAlgorithm const &AnalField) {
  RecVar eRecVar = RetrieveTrivialRecVar(eVarName);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int N = GrdArr.ARVD.N;
  int pos = PositionVect(AnalField.ListNameVariables, eVarName);
  if (pos == -1) {
    std::cerr << "For the model ANALYTICAL the variable eVarName = " << eVarName
              << " was requested\n";
    std::cerr << "But it is absent from the \"ListNameVariables\" field\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  std::string const &VarNature = eRecVar.RecS.VarNature;
  if (VarNature == "rho") {
    double eValRho = AnalField.ListConstantValuesRho[pos];
    eRecVar.F.setConstant(eta_rho, xi_rho, eValRho);
    return eRecVar;
  }
  if (VarNature == "uv") {
    double eValU = AnalField.ListConstantValuesU[pos];
    double eValV = AnalField.ListConstantValuesV[pos];
    eRecVar.U.setConstant(eta_rho, xi_rho, eValU);
    eRecVar.V.setConstant(eta_rho, xi_rho, eValV);
    eRecVar.F = GetNormMatrix(eRecVar.U, eRecVar.V);
    return eRecVar;
  }
  if (VarNature == "3Drho") {
    double eValRho = AnalField.ListConstantValuesRho[pos];
    eRecVar.Tens3 = ConstantThreeTensor(N, eta_rho, xi_rho, eValRho);
    return eRecVar;
  }
  if (VarNature == "3Duv") {
    double eValU = AnalField.ListConstantValuesU[pos];
    double eValV = AnalField.ListConstantValuesV[pos];
    eRecVar.Uthree = ConstantThreeTensor(N, eta_rho, xi_rho, eValU);
    eRecVar.Vthree = ConstantThreeTensor(N, eta_rho, xi_rho, eValV);
    eRecVar.Tens3 = ComputeNormPairOfTensor(eRecVar.Uthree, eRecVar.Vthree);
    return eRecVar;
  }
  std::cerr << "We have VarNature = " << VarNature << "\n";
  std::cerr << "Which is not matched\n";
  throw TerminalException{1};
}

std::pair<double, double> get_year_day(double const &eDate) {
  std::vector<int> LDate = DATE_ConvertMjd2six(eDate);
  double year = LDate[0];
  std::vector<int> LDateB{LDate[0], 1, 1, 0, 0, 0};
  double dateStartYear = DATE_ConvertSix2mjd(LDateB);
  double delta = eDate - dateStartYear;
  return {year, delta};
}

RecVar GetRecVarMeasurement(GridArray const &GrdArr,
                            std::string const &eVarName, double const &eTimeMjd,
                            CompleteMeasurementData const &cmd) {
  RecVar eRecVar = RetrieveTrivialRecVar(eVarName);
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  int N = GrdArr.ARVD.N;
  if (cmd.arr_meas.count(eVarName) == 0) {
    std::cerr << "For the model MEASUREMENT the variable eVarName = "
              << eVarName << " was requested\n";
    std::cerr << "But it is absent from the measurement data array\n";
    std::cerr << "Please correct\n";
    throw TerminalException{1};
  }
  const MeasurementData &md = cmd.arr_meas.at(eVarName);
  if (md.l_pos.size() == 0) {
    std::cerr << "We should hqve l_lon of non-trivial length\n";
    throw TerminalException{1};
  }
  std::string const &VarNature = eRecVar.RecS.VarNature;
  if (VarNature != "3Drho") {
    std::cerr << "Right now, only 3Drho is supported\n";
    throw TerminalException{1};
  }
  if (GrdArr.ModelName != "ROMS") {
    std::cerr << "Right now, only ROMS is supported\n";
    throw TerminalException{1};
  }
  MyMatrix<double> zeta = ZeroMatrix<double>(eta_rho, xi_rho);
  Eigen::Tensor<double, 3> Zr =
      ROMS_ComputeVerticalGlobalCoordinate_r(GrdArr, zeta);
  auto get_dist = [&](MeasurementSingPos const &e_msp,
                      MeasurementSingPos const &f_msp) -> double {
    double dx = e_msp.x - f_msp.x;
    double dy = e_msp.y - f_msp.y;
    double dz = e_msp.z - f_msp.z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double EarthRadius = 6371;
    double dist1 = 800;
    double coef1 = EarthRadius * dist / dist1;
    //
    double rel_dep = 2;
    double coef2 = T_abs(e_msp.dep - f_msp.dep) / rel_dep;
    //
    double rel_year = 2;
    double coef3 = T_abs(e_msp.year - f_msp.year) / rel_year;
    //
    double rel_day = 180;
    double coef4 = T_abs(e_msp.day - e_msp.day);
    if (coef4 > 183)
      coef4 = 365 - coef4;
    coef4 /= rel_day;
    //
    return coef1 + coef2 + coef3 + coef4;
  };
  auto get_value = [&](MeasurementSingPos const &f_msp) -> double {
    double sum = 0;
    double sum_val = 0;
    for (size_t i = 0; i < md.l_pos.size(); i++) {
      double eD = get_dist(md.l_pos[i], f_msp);
      double eW = 1 / eD;
      sum += eW;
      sum_val += eW * md.l_data[i];
    }
    return sum_val / sum;
  };
  std::pair<double, double> epair = get_year_day(eTimeMjd);
  Eigen::Tensor<double, 3> Tens3(N, eta_rho, xi_rho);
  TripleXYZ trip = ComputeTripleXYZ(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
  for (int iEta = 0; iEta < eta_rho; iEta++)
    for (int iXi = 0; iXi < xi_rho; iXi++)
      if (GrdArr.GrdArrRho.MSK(iEta, iXi) == 1) {
        double x = trip.X(iEta, iXi);
        double y = trip.Y(iEta, iXi);
        double z = trip.Z(iEta, iXi);
        for (int i = 0; i < N; i++) {
          double dep = Zr(i, iEta, iXi);
          MeasurementSingPos f_msp{x, y, z, dep, epair.first, epair.second};
          double interVal = get_value(f_msp);
          Tens3(i, iEta, iXi) = interVal;
        }
      }
  eRecVar.Tens3 = Tens3;
  return eRecVar;
}

bool HasField(CompleteMeasurementData const &cmd, std::string const &eVarName) {
  if (cmd.arr_meas.count(eVarName) == 1)
    return true;
  return false;
}

CompleteMeasurementData
ReadCompleteMeasurementData(std::vector<std::string> const &List_files) {
  std::unordered_map<std::string, MeasurementData> arr_meas;
  for (auto &eFile : List_files) {
    std::vector<std::string> ListLines = ReadFullFile(eFile);
    for (auto &eLine : ListLines) {
      std::vector<std::string> LStr = STRING_Split(eLine, " ");
      // Format is: VarName date lon lat dep meas
      std::string VarName = LStr[0];
      double date = DATE_ConvertString2mjd(LStr[1]);
      std::pair<double, double> epair = get_year_day(date);

      double lon = ParseScalar<double>(LStr[2]);
      double lat = ParseScalar<double>(LStr[3]);
      double dep = ParseScalar<double>(LStr[4]);
      double meas = ParseScalar<double>(LStr[5]);
      MeasurementData &md = arr_meas[VarName];
      std::vector<double> V = GetXYZcoordinateLL(lon, lat);
      double x = V[0];
      double y = V[1];
      double z = V[2];
      md.l_pos.push_back({x, y, z, dep, epair.first, epair.second});
      md.l_data.push_back(meas);
    }
  }
  return {arr_meas};
}

void Average_field_Function(FullNamelist const &eFull) {
  const std::map<std::string, SingleBlock> &ListBlock = eFull.ListBlock;
  SingleBlock eBlPROC = ListBlock.at("PROC");
  std::string ModelName = eBlPROC.ListStringValues.at("ModelName");
  std::string GridFile = eBlPROC.ListStringValues.at("GridFile");
  std::string HisPrefix = eBlPROC.ListStringValues.at("HisPrefix");
  //
  // Reading grid information.
  //
  TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
  GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
  ArrayHistory eArr = ReadArrayHistory(eTriple);
  //
  // Now the periods
  //
  SingleBlock eBlSELECT = ListBlock.at("SELECT");
  std::string Prefix = eBlSELECT.ListStringValues.at("Prefix");
  std::vector<std::string> ListNamesFile =
      eBlSELECT.ListListStringValues.at("ListNamesFile");
  std::vector<std::string> ListStartTime_str =
      eBlSELECT.ListListStringValues.at("ListStartTime");
  std::vector<std::string> ListEndTime_str =
      eBlSELECT.ListListStringValues.at("ListEndTime");
  std::vector<double> ListStartTime, ListEndTime;
  size_t n_ent = ListNamesFile.size();
  if (ListStartTime_str.size() == n_ent && ListEndTime_str.size() == n_ent) {
    for (auto &eTime_str : ListStartTime_str) {
      double eTime = CT2MJD(eTime_str);
      ListStartTime.push_back(eTime);
    }
    for (auto &eTime_str : ListEndTime_str) {
      double eTime = CT2MJD(eTime_str);
      ListEndTime.push_back(eTime);
    }
  } else {
    if (ListStartTime_str.size() != 0 || ListEndTime_str.size() != 0) {
      std::cerr << "If the ListStartTime and ListEndTime are not fully set\n";
      std::cerr << "then they have to be empty\n";
      throw TerminalException{1};
    }
    for (auto &eStr : ListNamesFile) {
      if (eStr.size() != 7) {
        std::cerr << "The entries of ListNamesFile need to be of the form "
                     "May2019 or such\n";
        throw TerminalException{1};
      }
      std::string eStrMonth = eStr.substr(0, 3);
      std::string eStrYear = eStr.substr(3, 4);
      int iMonth = GetIMonth(eStrMonth);
      int iYear = ParseScalar<int>(eStrYear);
      int month_len = MONTH_LEN(iYear, iMonth);
      double date_start = DATE_ConvertSix2mjd({iYear, iMonth, 1, 0, 0, 0});
      double date_end =
          DATE_ConvertSix2mjd({iYear, iMonth, month_len, 0, 0, 0});
      ListStartTime.push_back(date_start);
      ListEndTime.push_back(date_end);
    }
  }
  double DeltaT = 1 / static_cast<double>(24);
  for (size_t i_ent = 0; i_ent < n_ent; i_ent++) {
    std::string FullOutFile = Prefix + ListNamesFile[i_ent] + ".nc";
    RemoveFileIfExist(FullOutFile);
    std::cerr << "i_ent=" << i_ent << " FullOutFile=" << FullOutFile << "\n";
    std::vector<std::pair<std::string, size_t>> ListEnt;
    double eStartTime = ListStartTime[i_ent];
    double eEndTime = ListEndTime[i_ent];
    double eTime = eStartTime;
    while (true) {
      InterpInfo eInterp = GetTimeInterpolationInfoGeneralized(eArr, eTime);
      int iTimeLow = eInterp.iTimeLow;
      //
      std::vector<int> eRecLow = GetIFileIRec(eArr, iTimeLow);
      int iFile = eRecLow[0];
      int iRec = eRecLow[1];
      std::string eVar = "Uwind";
      std::string HisFile = ARR_GetHisFileName(eArr, eVar, iFile);
      std::cerr << "  iTimeLow=" << iTimeLow << " HisFile=" << HisFile
                << " iRec=" << iRec << "\n";
      ListEnt.push_back({HisFile, iRec});
      eTime += DeltaT;
      if (eTime >= eEndTime)
        break;
    }
    size_t n_part = ListEnt.size();
    std::cerr << "i_ent=" << i_ent << " |ListEnt|=" << n_part << "\n";
    //
    // Selecting entries in the database
    //
    std::vector<std::string> ListBlock;
    for (size_t i_part = 0; i_part < n_part; i_part++) {
      std::string command = "ncks";
      std::string blk_name = "/tmp/block_" + std::to_string(i_part) + ".nc";
      RemoveFileIfExist(blk_name);
      std::string eFile = ListEnt[i_part].first;
      int pos = ListEnt[i_part].second;
      std::string pos_s = std::to_string(pos);
      std::string order = command + " -d ocean_time," + pos_s + "," + pos_s +
                          " " + eFile + " " + blk_name;
      ListBlock.push_back(blk_name);
      //
      std::cerr << "i_part=" << i_part << " order=" << order << "\n";
      int iret1 = system(order.c_str());
      if (iret1 != 0) {
        std::cerr << "Error at ncks operation\n";
        throw TerminalException{1};
      }
    }
    //
    // Merging the separate netcdf files
    //
    std::string FileOut = "/tmp/Merge_" + std::to_string(i_ent) + ".nc";
    RemoveFileIfExist(FileOut);
    std::string order_concat = "ncrcat";
    for (size_t i_part = 0; i_part < n_part; i_part++) {
      order_concat += " " + ListBlock[i_part];
    }
    order_concat += " -o " + FileOut;
    std::cerr << "order_concat=" << order_concat << "\n";
    int iret2 = system(order_concat.c_str());
    if (iret2 != 0) {
      std::cerr << "Error at ncrcat operation\n";
      throw TerminalException{1};
    }
    //
    // Removing the files
    //
    for (size_t i_part = 0; i_part < n_part; i_part++)
      RemoveFileIfExist(ListBlock[i_part]);
    //
    // Computing the average
    //
    std::string order_avg =
        "ncwa -a ocean_time -b " + FileOut + " " + FullOutFile;
    std::cerr << "order_avg=" << order_avg << "\n";
    int iret3 = system(order_avg.c_str());
    if (iret3 != 0) {
      std::cerr << "Error at ncwa operation\n";
      throw TerminalException{1};
    }
    //
    // Removing the input
    //
    RemoveFileIfExist(FileOut);
  }
}

void INTERPOL_field_Function(FullNamelist const &eFull) {
  //
  // Construction of global arrays
  //
  std::map<std::string, SingleBlock> ListBlock = eFull.ListBlock;
  SingleBlock eBlINPUT = ListBlock.at("INPUT");
  std::vector<std::string> ListModelName =
      eBlINPUT.ListListStringValues.at("ListMODELNAME");
  std::vector<std::string> ListGridFile =
      eBlINPUT.ListListStringValues.at("ListGridFile");
  std::vector<std::string> ListHisPrefix =
      eBlINPUT.ListListStringValues.at("ListHisPrefix");
  std::vector<int> ListSpongeSize =
      eBlINPUT.ListListIntValues.at("ListSpongeSize");
  std::vector<int> ListFatherGrid =
      eBlINPUT.ListListIntValues.at("ListFatherGrid");
  bool PrintMMA = eBlINPUT.ListBoolValues.at("PrintMMA");
  int nbGrid = ListGridFile.size();
  std::cerr << "nbGrid=" << nbGrid << "\n";
  size_t nbGrid_s = nbGrid;
  if (ListModelName.size() != nbGrid_s || ListHisPrefix.size() != nbGrid_s ||
      ListFatherGrid.size() != nbGrid_s || ListSpongeSize.size() != nbGrid_s) {
    std::cerr << "Incoherent lengths of arrays\n";
    std::cerr << "|ListGridFile|   = " << ListGridFile.size() << "\n";
    std::cerr << "|ListModelName|  = " << ListModelName.size() << "\n";
    std::cerr << "|ListHisPrefix|  = " << ListHisPrefix.size() << "\n";
    std::cerr << "|ListFatherGrid| = " << ListFatherGrid.size() << "\n";
    std::cerr << "|ListSpongeSize| = " << ListSpongeSize.size() << "\n";
    throw TerminalException{1};
  }
  std::vector<GridArray> ListGrdArr;
  std::vector<ArrayHistory> ListArrayHistory;
  std::vector<TotalArrGetData> ListTotalArr;
  for (int iGrid = 0; iGrid < nbGrid; iGrid++) {
    std::string eModelName = ListModelName[iGrid];
    std::string GridFile = ListGridFile[iGrid];
    std::string HisPrefix = ListHisPrefix[iGrid];
    TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
    GridArray GrdArr = RETRIEVE_GRID_ARRAY(eTriple);
    ListGrdArr.push_back(GrdArr);
    ArrayHistory eArr = ReadArrayHistory(eTriple);
    ListArrayHistory.push_back(eArr);
    TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
    ListTotalArr.push_back(TotalArr);
  }
  bool DoClimatology = eBlINPUT.ListBoolValues.at("DoClimatology");
  bool AllowExtrapolation = eBlINPUT.ListBoolValues.at("AllowExtrapolation");
  std::cerr << "DoClimatology=" << DoClimatology
            << " AllowExtrapolation=" << AllowExtrapolation << "\n";
  std::cerr << "Arrays ListTotalArr, ListGrdArr and ListArrayHistory have been "
               "read\n";
  //
  // The target grid for the interpolation and the total array for interpolation
  //
  SingleBlock eBlOUTPUT = ListBlock.at("OUTPUT");
  std::string eModelName = eBlOUTPUT.ListStringValues.at("MODELNAME");
  std::string GridFile = eBlOUTPUT.ListStringValues.at("GridFile");
  std::string BoundFile = eBlOUTPUT.ListStringValues.at("BoundFile");
  std::string HisPrefix = eBlOUTPUT.ListStringValues.at("HisPrefix");
  //
  // The analytical arrays if needed
  //
  SingleBlock eBlANALYTIC = ListBlock.at("ANALYTIC");
  std::vector<std::string> AnalyticalListNameVariables =
      eBlANALYTIC.ListListStringValues.at("AnalyticalListNameVariables");
  std::vector<double> AnalyticalListConstantValuesRho =
      eBlANALYTIC.ListListDoubleValues.at("AnalyticalListConstantValuesRho");
  std::vector<double> AnalyticalListConstantValuesU =
      eBlANALYTIC.ListListDoubleValues.at("AnalyticalListConstantValuesU");
  std::vector<double> AnalyticalListConstantValuesV =
      eBlANALYTIC.ListListDoubleValues.at("AnalyticalListConstantValuesV");
  AnalyticalAlgorithm AnalField{
      AnalyticalListNameVariables, AnalyticalListConstantValuesRho,
      AnalyticalListConstantValuesU, AnalyticalListConstantValuesV};
  std::cerr << "Analytical fields have been read\n";
  //
  // The measurement data arrays if needed
  //
  SingleBlock eBlMEAS = ListBlock.at("MEASUREMENT");
  std::vector<std::string> List_files =
      eBlMEAS.ListListStringValues.at("ListFilesMeasurement");
  CompleteMeasurementData cmd = ReadCompleteMeasurementData(List_files);
  //
  // ROMS boundary related stuff
  //
  bool DoSfluxWrite = eBlOUTPUT.ListBoolValues.at("DoSfluxWrite");
  bool DoNetcdfWrite = eBlOUTPUT.ListBoolValues.at("DoNetcdfWrite");
  bool DoGribWrite = eBlOUTPUT.ListBoolValues.at("DoGribWrite");
  bool DoRomsWrite_Surface = eBlOUTPUT.ListBoolValues.at("DoRomsWrite_Surface");
  bool DoRomsWrite_InitialHistory =
      eBlOUTPUT.ListBoolValues.at("DoRomsWrite_InitialHistory");
  bool DoRomsWrite_Boundary =
      eBlOUTPUT.ListBoolValues.at("DoRomsWrite_Boundary");
  bool DoWaveWatchWrite = eBlOUTPUT.ListBoolValues.at("DoWaveWatchWrite");
  std::string WaveWatchFormat = eBlOUTPUT.ListStringValues.at("WaveWatchFormat");
  int nbTypeOutput = 0;
  if (DoSfluxWrite)
    nbTypeOutput++;
  if (DoNetcdfWrite)
    nbTypeOutput++;
  if (DoGribWrite)
    nbTypeOutput++;
  if (DoRomsWrite_Surface)
    nbTypeOutput++;
  if (DoRomsWrite_InitialHistory)
    nbTypeOutput++;
  if (DoRomsWrite_Boundary)
    nbTypeOutput++;
  if (DoWaveWatchWrite)
    nbTypeOutput++;
  if (nbTypeOutput != 1) {
    std::cerr << "We have DoSfluxWrite = " << DoSfluxWrite << "\n";
    std::cerr << "We have DoNetcdfWrite = " << DoNetcdfWrite << "\n";
    std::cerr << "We have DoGribWrite = " << DoGribWrite << "\n";
    std::cerr << "We have DoRomsWrite_Surface = " << DoRomsWrite_Surface
              << "\n";
    std::cerr << "We have DoRomsWrite_InitialHistory = "
              << DoRomsWrite_InitialHistory << "\n";
    std::cerr << "We have DoRomsWrite_Boundary = " << DoRomsWrite_Boundary
              << "\n";
    std::cerr << "We have DoWaveWatchWrite = " << DoWaveWatchWrite << "\n";
    std::cerr << "We have nbTypeOutput = " << nbTypeOutput << "\n";
    std::cerr << "We should select exactly 1 output\n";
    throw TerminalException{1};
  }
  std::cerr << "Selection of output type done\n";
  //
  // The output grid
  //
  double MinLat = eBlOUTPUT.ListDoubleValues.at("MinLat");
  double MaxLat = eBlOUTPUT.ListDoubleValues.at("MaxLat");
  double MinLon = eBlOUTPUT.ListDoubleValues.at("MinLon");
  double MaxLon = eBlOUTPUT.ListDoubleValues.at("MaxLon");
  double deltaKM = eBlOUTPUT.ListDoubleValues.at("deltaKM");
  std::string Sphericity = "Spherical";
  GridSymbolic RecGridSymb(Sphericity, false, false, 0, 0, MinLat, MaxLat,
                           MinLon, MaxLon, deltaKM);
  TripleModelDesc eTripleOut{eModelName, GridFile, BoundFile, HisPrefix,
                             RecGridSymb};
  std::cerr
      << "INTERPOL_field_Function : Before RETRIEVE_GRID_ARRAY eModelName="
      << eModelName << "\n";
  GridArray GrdArrOut = RETRIEVE_GRID_ARRAY(eTripleOut);
  std::cerr << "INTERPOL_field_Function : After RETRIEVE_GRID_ARRAY\n";
  //
  // Reading the relevant variables of the output
  //
  std::vector<std::string> ListVarName =
      NAMELIST_ListTrueEntryBool(eFull, "VARS");
  //
  // The ROMS functionality for surface forcing
  //
  std::vector<ROMS_NC_VarInfo> ListArrROMS;
  if (DoRomsWrite_Surface) {
    SingleBlock eBlSURF = ListBlock.at("ROMS_SURFACE");
    std::string RomsFileNC_surf = eBlSURF.ListStringValues.at("RomsFile_surf");
    bool IsRegrid = eBlSURF.ListBoolValues.at("IsRegrid");
    bool SingleFile = eBlSURF.ListBoolValues.at("SingleFile");
    if (IsRegrid) {
      if (nbGrid != 1) {
        std::cerr << "nbGrid = " << nbGrid << " and it should be equal to 1\n";
        throw TerminalException{1};
      }
      GrdArrOut = ListGrdArr[0];
    }
    ListArrROMS = ROMS_Surface_NetcdfInitialize(
        RomsFileNC_surf, IsRegrid, SingleFile, GrdArrOut, ListVarName);
    std::cerr << "After DoRomsWrite_Surface initialization\n";
  }
  //
  // The ROMS initial file
  //
  std::string RomsFileNC_InitialHistory;
  if (DoRomsWrite_InitialHistory) {
    SingleBlock eBlROMS_INIT_HIS = ListBlock.at("ROMS_INITIAL_HISTORY");
    RomsFileNC_InitialHistory =
        eBlROMS_INIT_HIS.ListStringValues.at("RomsFile_InitialHistory");
    int N = eBlROMS_INIT_HIS.ListIntValues.at("ARVD_N");
    int Vtransform = eBlROMS_INIT_HIS.ListIntValues.at("ARVD_Vtransform");
    int Vstretching = eBlROMS_INIT_HIS.ListIntValues.at("ARVD_Vstretching");
    double Tcline = eBlROMS_INIT_HIS.ListDoubleValues.at("ARVD_Tcline");
    double hc = eBlROMS_INIT_HIS.ListDoubleValues.at("ARVD_hc");
    double theta_s = eBlROMS_INIT_HIS.ListDoubleValues.at("ARVD_theta_s");
    double theta_b = eBlROMS_INIT_HIS.ListDoubleValues.at("ARVD_theta_b");
    GrdArrOut.ARVD = ROMSgetARrayVerticalDescription(
        N, Vtransform, Vstretching, Tcline, hc, theta_s, theta_b);
    std::vector<std::string> ListAddiVarnameROMS;
    std::vector<std::string> ListClassic = {"Temp", "Salt", "ZetaOcean", "Curr",
                                            "CurrBaro"};
    for (auto &eVarName : ListVarName) {
      if (PositionVect(ListClassic, eVarName) == -1) {
        RecVar eRecVar = RetrieveTrivialRecVar(eVarName);
        if (!eRecVar.RecS.varName_ROMS) {
          std::cerr << "varName_ROMS has not been assigned 6\n";
          std::cerr << "eVarName=" << eVarName << "\n";
          std::cerr << "VarName1=" << eRecVar.RecS.VarName1 << "\n";
          std::cerr << "VarName2=" << eRecVar.RecS.VarName2 << "\n";
          throw TerminalException{1};
        }
        std::string const &VarNameRoms = *eRecVar.RecS.varName_ROMS;
        ListAddiVarnameROMS.push_back(VarNameRoms);
        std::cerr << "eVarName=" << eVarName << " VarNameRoms=" << VarNameRoms
                  << "\n";
      }
    }
    std::cerr << "Before ROMS_InitialHistory_NetcdfInitialize\n";
    ROMS_InitialHistory_NetcdfInitialize(RomsFileNC_InitialHistory, GrdArrOut,
                                         ListAddiVarnameROMS);
    std::cerr << "After  ROMS_InitialHistory_NetcdfInitialize\n";
  }
  //
  // The ROMS functionality for boundary forcing
  //
  SingleBlock eBLROMS_BOUND = ListBlock.at("ROMS_BOUND");
  std::vector<std::string> ListSides =
      eBLROMS_BOUND.ListListStringValues.at("ListSides");
  std::string RomsFileNC_bound =
      eBLROMS_BOUND.ListStringValues.at("RomsFile_bound");
  std::cerr << "DoRomsWrite_Boundary=" << DoRomsWrite_Boundary << "\n";
  if (DoRomsWrite_Boundary) {
    int N = eBLROMS_BOUND.ListIntValues.at("ARVD_N");
    int Vtransform = eBLROMS_BOUND.ListIntValues.at("ARVD_Vtransform");
    int Vstretching = eBLROMS_BOUND.ListIntValues.at("ARVD_Vstretching");
    double Tcline = eBLROMS_BOUND.ListDoubleValues.at("ARVD_Tcline");
    double hc = eBLROMS_BOUND.ListDoubleValues.at("ARVD_hc");
    double theta_s = eBLROMS_BOUND.ListDoubleValues.at("ARVD_theta_s");
    double theta_b = eBLROMS_BOUND.ListDoubleValues.at("ARVD_theta_b");
    GrdArrOut.ARVD = ROMSgetARrayVerticalDescription(
        N, Vtransform, Vstretching, Tcline, hc, theta_s, theta_b);
    int eta_rho = GrdArrOut.GrdArrRho.MSK.rows();
    int xi_rho = GrdArrOut.GrdArrRho.MSK.cols();
    auto ComputeBndStat =
        [&](std::string const &eSide,
            std::vector<std::pair<int, int>> const &ListPairIdx) -> void {
      int sumMSK = 0;
      double sumLon = 0, sumLat = 0;
      std::vector<double> LLon, LLat;
      for (auto epair : ListPairIdx) {
        sumMSK += GrdArrOut.GrdArrRho.MSK(epair.first, epair.second);
        double eLon = GrdArrOut.GrdArrRho.LON(epair.first, epair.second);
        double eLat = GrdArrOut.GrdArrRho.LAT(epair.first, epair.second);
        LLon.push_back(eLon);
        LLat.push_back(eLat);
        sumLon += eLon;
        sumLat += eLat;
      }
      int pos = PositionVect(ListSides, eSide);
      if (pos == -1 && sumMSK > 0) {
        std::cerr << "The side " << eSide
                  << " is not included in the ListSides.\n";
        std::cerr
            << "However, we have some boundary point. This could be an error\n";
      }
      if (pos != -1 && sumMSK == 0) {
        std::cerr << "The side " << eSide << " is selected in ListSides.\n";
        std::cerr << "However, the boundary has no masking point.\n";
        std::cerr << "There is no scenario n in which this makes sense\n";
        throw TerminalException{1};
      }
      double avgLon = sumLon / ListPairIdx.size();
      double avgLat = sumLat / ListPairIdx.size();
      std::cerr << "  |Side|=" << ListPairIdx.size() << " name=" << eSide
                << " sumMSK=" << sumMSK << " avgLon=" << avgLon
                << " avgLat=" << avgLat << "\n";
      std::cerr << "     LON(min/max)=" << VectorMin(LLon) << " / "
                << VectorMax(LLon) << " LAT(min/max)=" << VectorMin(LLat)
                << " / " << VectorMax(LLat) << "\n";
    };
    std::vector<std::pair<int, int>> ListPairEast, ListPairWest, ListPairSouth,
        ListPairNorth;
    for (int i = 0; i < xi_rho; i++) {
      ListPairSouth.push_back({0, i});
      ListPairNorth.push_back({eta_rho - 1, i});
    }
    for (int i = 0; i < eta_rho; i++) {
      ListPairEast.push_back({i, xi_rho - 1});
      ListPairWest.push_back({i, 0});
    }
    ComputeBndStat("South", ListPairSouth);
    ComputeBndStat("North", ListPairNorth);
    ComputeBndStat("West", ListPairWest);
    ComputeBndStat("East", ListPairEast);
    //
    std::cerr << "Before call to ROMS_BOUND_NetcdfInitialize\n";
    std::vector<RecVar> ListArrayTracer =
        GetListArrayTracerTrivial(ListVarName);
    ROMS_BOUND_NetcdfInitialize(RomsFileNC_bound, GrdArrOut, ListSides,
                                ListArrayTracer);
    std::cerr << " After call to ROMS_BOUND_NetcdfInitialize\n";
  }
  std::cerr << "After DoRomsWrite_Boundary initialization\n";
  //
  // The interpolation arrays
  //
  TotalArrayInterpolation TotalArrInt =
      INTERPOL_ConstructTotalArray(ListTotalArr, ListSpongeSize, ListFatherGrid,
                                   GrdArrOut, AllowExtrapolation);
  std::cerr << "NeedInterp=" << TotalArrInt.NeedInterp << "\n";
  std::cerr << "We have the interpolation array TotalArrInt\n";
  auto GetRecVarInterpolate = [&](std::string const &eVarName,
                                  double const &eTimeDay) -> RecVar {
    if (!DoClimatology)
      return INTERPOL_MultipleRecVarInterpolation(TotalArrInt, GrdArrOut,
                                                  eVarName, eTimeDay);
    // Now doing the climatological computation
    int eYearStart = DATE_ConvertMjd2six(TotalArrInt.StartTime)[0];
    int eYearEnd = DATE_ConvertMjd2six(TotalArrInt.EndTime)[0];
    int YearBegin = eYearStart - 1;
    int YearEnd = eYearEnd + 1;
    std::vector<RecVar> ListRecVar;
    auto FuncInsert = [&](double const &eTimeDayIns) -> void {
      std::vector<int> eDateIns = DATE_ConvertMjd2six(eTimeDayIns);
      for (int eYearW = YearBegin; eYearW <= YearEnd; eYearW++) {
        std::vector<int> eDateW = eDateIns;
        eDateW[0] = eYearW;
        if (TestCorrectnessVectorTime(eDateW).first) {
          double eTimeW = DATE_ConvertSix2mjd(eDateW);
          if (TotalArrInt.StartTime <= eTimeW && eTimeW <= TotalArrInt.EndTime)
            ListRecVar.push_back(INTERPOL_MultipleRecVarInterpolation(
                TotalArrInt, GrdArrOut, eVarName, eTimeW));
        }
      }
    };
    //
    FuncInsert(eTimeDay);
    if (ListRecVar.size() == 0) {
      FuncInsert(eTimeDay - 1);
      FuncInsert(eTimeDay + 1);
    }
    //
    if (ListRecVar.size() == 0) {
      std::cerr << "We found zero relevant record, therefore we cannot compute "
                   "the average of them\n";
      throw TerminalException{1};
    }
    return Average_RecVar(ListRecVar);
  };
  //
  // timings for all the model output
  //
  std::vector<double> ListTime = GetIntervalGen(eBlOUTPUT, ListArrayHistory);
  int nbTime = ListTime.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  double DEFINETC = eBlOUTPUT.ListDoubleValues.at("DEFINETC");
  double DELTC = eBlOUTPUT.ListDoubleValues.at("DELTC");
  int eMult = static_cast<int>(round(DEFINETC / DELTC));
  double eDiff = DEFINETC - static_cast<double>(eMult) * DELTC;
  if (fabs(eDiff) > 1) {
    std::cerr << "nbTime=" << nbTime << "\n";
    std::cerr << "eDiff=" << eDiff << "\n";
    std::cerr << "DEFINETC=" << DEFINETC << " DELTC=" << DELTC << "\n";
    std::cerr << "DEFINETC should be an integer multiple of\n";
    std::cerr << "of DELTC\n";
    throw TerminalException{1};
  }
  std::cerr << "Creation of defining times\n";
  //
  // WaveWatch III related variable
  //
  int WWIII_nbWritten = 0;
  //
  // The relevant variable for the netcdf output
  //
  recNetcdfOutput recNO;
  recNO.iFile = 1;
  recNO.nbWritten = 0;
  recGribOutput recGO;
  recGO.nbWritten = 0;
  auto get_rec_var = [&](std::string const &eVarName,
                         double const &eTimeDay) -> RecVar {
    if (HasField(AnalField, eVarName)) {
      std::cerr << "Retrieving eVarName=" << eVarName << " analytically\n";
      return GetRecVarAnalytical(GrdArrOut, eVarName, eTimeDay, AnalField);
    }
    if (HasField(cmd, eVarName)) {
      std::cerr << "Retrieving eVarName=" << eVarName << " from measurements\n";
      return GetRecVarMeasurement(GrdArrOut, eVarName, eTimeDay, cmd);
    }
    std::cerr << "Retrieving eVarName=" << eVarName
              << " from interpolation of model fields\n";
    return GetRecVarInterpolate(eVarName, eTimeDay);
  };
  //
  // The big bad loop
  //
  for (int iTime = 0; iTime < nbTime; iTime++) {
    double eTimeDay = ListTime[iTime];
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << "  date=" << strPres
              << "\n";
    //
    // Data retrieval
    //
    std::vector<RecVar> ListRecVar;
    for (auto &eVarName : ListVarName) {
      RecVar eRecVar = get_rec_var(eVarName, eTimeDay);
      Set_iTime_eTimeDay(eRecVar, iTime, eTimeDay);
      ListRecVar.push_back(eRecVar);
    }
    //
    // Write standard netcdf output
    //
    if (DoNetcdfWrite)
      INTERPOL_NetcdfOutput(GrdArrOut, ListRecVar, eFull, eMult, recNO);
    //
    // Write standard grib output
    //
    if (DoGribWrite)
      INTERPOL_GribOutput(GrdArrOut, ListRecVar, eFull, eMult, recGO);
    //
    // Write .ww3 wavewatch III forcing
    //
    if (DoWaveWatchWrite)
      WaveWatch_WriteData(GrdArrOut, ListRecVar, WWIII_nbWritten, WaveWatchFormat);
    //
    // Write SFLUX files
    //
    if (DoSfluxWrite) {
      std::cerr << "The code for SFLUX has to be written\n";
      throw TerminalException{1};
    }
    //
    // Write ROMS surface forcing file
    //
    if (DoRomsWrite_Surface)
      ROMS_Surface_NetcdfAppendVarName(GrdArrOut, ListRecVar, ListArrROMS,
                                       PrintMMA);
    //
    // Write ROMS initial file (just one entry) or depending on the viewpoint
    // the history file (several entries)
    //
    if (DoRomsWrite_InitialHistory) {
      ROMSstate eState = GetRomsStateFromVariables(GrdArrOut, ListRecVar);
      std::cerr << "Before ROMS_InitialHistory_NetcdfAppend\n";
      ROMS_InitialHistory_NetcdfAppend(RomsFileNC_InitialHistory, eState,
                                       GrdArrOut, iTime);
      std::cerr << " After ROMS_InitialHistory_NetcdfAppend\n";
    }
    //
    // Write ROMS boundary forcing
    //
    if (DoRomsWrite_Boundary) {
      ROMSstate eState = GetRomsStateFromVariables(GrdArrOut, ListRecVar);
      ROMS_BOUND_NetcdfAppend(RomsFileNC_bound, eState, ListSides, iTime);
    }
  }
}

// clang-format off
#endif  // SRC_OCEAN_MODEL_INTERPOLATION_H_
// clang-format on
