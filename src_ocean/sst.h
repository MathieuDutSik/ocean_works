#ifndef INCLUDE_SST_H
#define INCLUDE_SST_H

#include "NamelistExampleOcean.h"
#include "Statistics.h"
#include "Plotting_fct.h"


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

void RAW_SCATTER_SST(std::vector<double> const& V_meas, std::vector<double> const& V_model,
                     NCLcaller<GeneralType> & eCall,
                     PermanentInfoDrawing const& ePerm)
{
  std::cerr << "Running RAW_SCATTER_SST\n";
  SingleBlock eBlPLOT=ePerm.eFull.ListBlock.at("PLOT");
  bool UseDynamicRangeInScatter=eBlPLOT.ListBoolValues.at("UseDynamicRangeInScatter");
  //
  DrawScatterArr eDrw;
  int siz=V_meas.size();
  MyVector<double> eVectA(siz);
  MyVector<double> eVectB(siz);
  for (int i=0; i<siz; i++) {
    eVectA[i]=V_meas[i];
    eVectB[i]=V_model[i];
  }
  std::vector<double> data_rangeA(2);
  std::vector<double> data_rangeB(2);
  //
  std::string VarType="SST";
  std::string eUnit="deg C";
  double TheMin = 10;
  double TheMax = 25;
  if (UseDynamicRangeInScatter && siz > 0) {
    double maxA = eVectA.maxCoeff();
    double maxB = eVectB.maxCoeff();
    TheMax = std::max(maxA, maxB);
    double minA = eVectA.maxCoeff();
    double minB = eVectB.maxCoeff();
    TheMin = std::max(minA, minB);
  }
  data_rangeA[0]=TheMin;
  data_rangeA[1]=TheMax;
  data_rangeB[0]=TheMin;
  data_rangeB[1]=TheMax;
  eDrw.VarNameAB_file="Scatter_iGrid_" + VarType;
  //
  eDrw.DoTitle=false;
  eDrw.AddStatMeasModel=true;
  eDrw.NameA_plot="Data (" + eUnit + ")";
  eDrw.NameB_plot="Model (" + eUnit + ")";
  eDrw.data_rangeA=data_rangeA;
  eDrw.data_rangeB=data_rangeB;
  eDrw.eVectA=eVectA;
  eDrw.eVectB=eVectB;
  eDrw.aSize=100;
  eDrw.bSize=100;
  PLOT_SCATTER(eDrw, eCall, ePerm);
}







void Process_sst_Comparison_Request(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  SingleBlock eBlSTAT=eFull.ListBlock.at("STAT");
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  //
  // Now basic definitions
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  ePerm.eDrawArr = CommonAssignation_DrawArr(ePerm.eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC); // has to be after ePerm
  //
  // Reading the model
  //
  std::string ModelName = eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile  = eBlPROC.ListStringValues.at("GridFile");
  std::string HisPrefix = eBlPROC.ListStringValues.at("HisPrefix");
  TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
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
    std::cerr << "iFile=" << iFile << " eFile=" << eFile << "\n";
    std::vector<double> LTime = NC_ReadTimeFromFile(eFile, "time");
    for (size_t iTime=0; iTime<LTime.size(); iTime++)
      ListSingEnt.push_back({iFile, iTime, LTime[iTime]});
  }
  std::cerr << "|ListSingEnt|=" << ListSingEnt.size() << "\n";
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
  double PreDawnHour = eBlPROC.ListDoubleValues.at("PreDawnHour");
  if (!IsZeroHour(BeginTime)) {
    std::string strPresBegin = DATE_ConvertMjd2mystringPres(BeginTime);
    std::cerr << "The initial date should be a zero hour\n";
    std::cerr << "That is hour=0 , min=0 , sec=0\n";
    std::cerr << "On the other hand we have BeginTime=" << strPresBegin << "\n";
    throw TerminalException{1};
  }
  std::cerr << "We have BeginTime, EndTime, PreDawnTime\n";
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
  MyMatrix<uint8_t> const& MSK = TotalArr.GrdArr.GrdArrRho.MSK;
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
  std::cerr << "We have MatIdxLat, MatIdxLon\n";
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
  auto ReadSST_Pair=[&](std::pair<size_t, size_t> const& ePair) -> std::pair<MyMatrix<double>, MyMatrix<double>> {
    MyMatrix<double> Mat_SST = ReadSST_entry(ePair, "analysed_sst");
    MyMatrix<double> Mat_ERR = ReadSST_entry(ePair, "analysis_error");
    int eta_rho = Mat_ERR.rows();
    int xi_rho = Mat_ERR.rows();
    for (int iEta=0; iEta<eta_rho; iEta++)
      for (int iXi=0; iXi<xi_rho; iXi++) {
        if (Mat_SST(iEta, iXi) < -50)
          Mat_ERR(iEta, iXi) = 300;
      }
    return {std::move(Mat_SST), std::move(Mat_ERR)};
  };
  double MaxErr_L4 = eBlSTAT.ListDoubleValues.at("MaxErr_L4");
  std::cerr << "MaxErr_L4=" << MaxErr_L4 << "\n";
  std::string OutPrefix = eBlSTAT.ListStringValues.at("OutPrefix");
  std::string FileStatDaily = OutPrefix + "Statistics_Daily.txt";
  std::ofstream os(FileStatDaily);
  std::vector<double> V_meas_total;
  std::vector<double> V_model_total;
  for (double eTime = BeginTime; eTime <= EndTime; eTime += 1.0) {
    std::string strPres = DATE_ConvertMjd2mystringPres(eTime);
    std::cerr << "eTime=" << eTime << " date=" << strPres << "\n";
    double eTimeCall = eTime + PreDawnHour;
    RecVar eRecVar = ModelSpecificVarSpecificTime_Kernel(TotalArr, "TempSurf", eTimeCall);
    //
    std::pair<size_t, size_t> ePair = GetEntry(eTime);
    auto ePairSST_ERR = ReadSST_Pair(ePair);
    MyMatrix<double> const& Mat_SST = ePairSST_ERR.first;
    MyMatrix<double> const& Mat_ERR = ePairSST_ERR.second;
    //
    std::vector<double> V_meas;
    std::vector<double> V_model;
    double CorrKelvin_Celsius = 273.15;
    for (size_t iEta=0; iEta<eta_rho; iEta++)
      for (size_t iXi=0; iXi<xi_rho; iXi++) {
        double eValModel = eRecVar.F(iEta, iXi);
        if (MSK(iEta, iXi) == 1) {
          int i_lat = MatIdxLat(iEta, iXi);
          int i_lon = MatIdxLon(iEta, iXi);
          if (Mat_ERR(i_lat, i_lon) < MaxErr_L4) {
            double eValMeas = Mat_SST(i_lat, i_lon) - CorrKelvin_Celsius;
            V_meas.push_back(eValMeas);
            V_model.push_back(eValModel);
            V_meas_total.push_back(eValMeas);
            V_model_total.push_back(eValModel);
          }
        }
      }
    // Now computing the stats
    T_stat estat = ComputeStatistics_vector(V_meas, V_model);
    T_statString estatstr = ComputeStatisticString_from_Statistics(estat, "4dot2f");
    std::cerr << "    date=" << strPres << " stat=" << estatstr.str << "\n";
    std::cerr << "    nbMeas=" << estat.nbMeas << " MeanMeas=" << estat.MeanMeas << " MeanModel=" << estat.MeanModel << "\n";
    std::cerr << "    MinMeas=" << estat.MinMeas << " MaxMeas=" << estat.MaxMeas << "\n";
  }
  RAW_SCATTER_SST(V_meas_total, V_model_total,
                  eCall, ePerm);
}







#endif


