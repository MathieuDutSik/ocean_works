#ifndef INCLUDE_SST_H
#define INCLUDE_SST_H

#include "NamelistExampleOcean.h"


int GetSmallestIndex(MyVector<double> const& V, double eVal)
{
  int pos = -1;
  double MinDist=100000000;
  for (int idx=0; idx<int(V.size()); idx++) {
    double eDist = fabs(V(idx) - eVal);
    if (eDist < MinDist) {
      MinDist = eDist;
      pos=idx;
    }
  }
  return pos;
}


void Process_sst_Comparison_Request(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  SingleBlock eBlSTAT=eFull.ListBlock.at("STAT");
  //
  // Reading the model
  //
  std::string ModelName = eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile  = eBlPROC.ListStringValues.at("GridFile");
  std::string HisPrefix = eBlPROC.ListStringValues.at("HisPrefix");
  TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, {}};
  TotalArrGetData TotalArr = RetrieveTotalArr(eTriple);
  //
  // Reading the list of files and times.
  //
  std::string SST_files_prefix = eBlPROC.ListStringValues.at("SST_files_prefix");
  std::vector<std::string> ListFile = FILE_DirectoryMatchingPrefixExtension(SST_files_prefix, "nc");
  struct SingEnt {
    size_t iFile;
    size_t iTime;
    double eTime;
  };
  std::vector<SingEnt> ListSingEnt;
  for (size_t iFile=0; iFile<ListFile.size(); iFile++) {
    std::string eFile = ListFile[iFile];
    std::vector<double> LTime = NC_ReadTimeFromFile(eFile, "time");
    for (size_t iTime=0; iTime<LTime.size(); iTime++)
      ListEntries.push_back({iFile, iTime, LTime[iTime]});
  }
  //
  // Determining the beginning and ending of time for comparison
  //
  std::string strBEGTC = eBlPROC.ListStringValues.at("BEGTC");
  std::string strENDTC = eBlPROC.ListStringValues.at("ENDTC");
  double BeginTime=0, EndTime=0;
  if (strBEGTC == "earliest") {
    BeginTime = MinimumTimeHistoryArray(TotalArr.eArr);
  } else {
    BeginTime = CT2MJD(strBEGTC);
  }
  if (strENDTC == "latest") {
    EndTime = MaximumTimeHistoryArray(TotalArr.eArr);
  } else {
    EndTime = CT2MJD(strBEGTC);
  }
  //
  // Reading the SST grid and computing interpolation arrays
  //
  std::string eFileSST = ListFile[0];
  MyVector<double> VectLat = NC_Read1Dvariable(eFileSST, "lat");
  MyVector<double> VectLon = NC_Read1Dvariable(eFileSST, "lon");
  size_t n_lat = VectLat.size();
  size_t n_lon = VectLon.size();
  MyMatrix<double> const& LON = TotalArr.GrdArr.GrdArrRho.LON;
  MyMatrix<double> const& LAT = TotalArr.GrdArr.GrdArrRho.LAT;
  MyMatrix<uint8_t> const& MSK = TotalArr.GrdArr.GrdArrRho.LAT;
  size_t eta_rho = LON.rows();
  size_t xi_rho = LON.cols();
  MyMatrix<int> MatIdxLat(eta_rho, xi_rho);
  MyMatrix<int> MatIdxLon(eta_rho, xi_rho);
  for (size_t iEta=0; iEta<eta_rho; iEta++)
    for (size_t iXi=0; iXi<xi_rho; iXi++) {
      double eLon = LON(iEta, iXi);
      double eLat = LAT(iEta, iXi);
      //
      MatIdxLat(iEta, iXi) = GetSmallestIndex(VectLat, eLat);
      MatIdxLon(iEta, iXi) = GetSmallestIndex(VectLon, eLon);
    }
  //
  // Now processing the comparison
  //
  auto GetEntry=[&](double const& eTime) -> std::pair<size_t, size_t> {
    for (auto & eSingEnt : ListSingEnt) {
      if (fabs(eSingEnt.eTime - eTime) < 0.0001) {
        return {eSingEnt.iFile, eSingEnt.iTime};
      }
    }
    return {-1, -1};
  };
  auto ReadSST_entry=[&](std::pair<size_t, size_t> const& ePair, std::string const& eVar) -> MyMatrix<double> {
    std::string eFile = ListFile[ePair.first];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar data=dataFile.getVar(eVar);
    //
    std::vector<size_t> start{ePair.second, 0, 0};
    std::vector<size_t> count{1,n_lat, n_lon};
    MyVector<double> V = NC_ReadVariable_data_start_count(data, start, count);
    MyMatrix<double> M(n_lat, n_lon);
    size_t idx=0;
    for (size_t i_lat=0; i_lat<n_lat; i_lat++)
      for (size_t i_lon=0; i_lon<n_lon; i_lon++) {
        M(i_lat, i_lon) = V(idx);
        idx++;
      }
    return M;
  };
  double MaxErr_L4 = eBlSTAT.ListDoubleValues.at("MaxErr_L4");
  for (double eTime = BeginTime; eTime <= EndTime; eTime += 1.0) {
    RecVar eRecVar = ModelSpecificVarSpecificTime_Kernel(TotalArr, "TempSurf", eTime);
    //
    std::pair<size_t, size_t> ePair = GetEntry(eTime);
    MyMatrix<double> Mat_ERR = ReadSST_entry(ePair, "analysis_error");
    MyMatrix<double> Mat_SST = ReadSST_entry(ePair, "analysed_sst");
    //
    std::vector<double> V_meas;
    std::vector<double> V_model;
    for (size_t iEta=0; iEta<eta_rho; iEta++)
      for (size_t iXi=0; iXi<xi_rho; iXi++) {
        double eValModel = eRecVar.F(iEta, iXi);
        if (MSK(iEta, iXi) == 1) {
          int i_lat = MatIdxLat(iEta, iXi);
          int i_lon = MatIdxLon(iEta, iXi);
          if (Mat_ERR(i_lat, i_lon) < MaxErr_L4) {
            double eValMeas = Mat_SST(i_lat, i_lon);
            V_meas.push_back(eValMeas);
            V_model.push_back(eValModel);
          }
        }
      }
    // Now computing the stats
    T_stat estat = ComputeStatistics_vector(V_meas, V_model);
    T_statString estatstr = ComputeStatisticString_from_Statistics(estat, "4dot2f");
    std::cerr << "stat=" << estatstr.str << "\n";

  }
}







#endif


