// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_RIVER_H_
#define SRC_OCEAN_RIVER_H_

#include "Basic_file.h"
#include "Basic_netcdf.h"
#include "Basic_plot.h"
#include "Model_grids.h"
#include "Namelist.h"
#include "ROMSfunctionality.h"
#include "SVGfunctions.h"
#include "SphericalGeom.h"
#include "mjdv2.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

struct PairTimeMeas {
  double time;
  double meas;
};

std::vector<PairTimeMeas>
ReadFileInterpolationInformation(std::string const &eFile) {
  std::vector<PairTimeMeas> l_pairs;
  std::vector<std::string> ListLines = ReadFullFile(eFile);
  for (auto &eLine : ListLines) {
    std::vector<std::string> LStrA = STRING_Split(eLine, " ");
    if (LStrA.size() != 2) {
      std::cerr << "eLine=" << eLine << "\n";
      std::cerr << "Required format is 20160120.000000 45.0\n";
      throw TerminalException{1};
    }
    double eTime = DATE_ConvertString2mjd(LStrA[0]);
    double value = ParseScalar<double>(LStrA[1]);
    PairTimeMeas eP{eTime, value};
    l_pairs.push_back(eP);
  }
  size_t len = l_pairs.size();
  for (size_t i = 1; i < len; i++) {
    double time0 = l_pairs[i - 1].time;
    double time1 = l_pairs[i].time;
    if (time0 > time1) {
      std::cerr << "The entries in the file are not correctly ordered\n";
      std::cerr << "i=" << i << " time0=" << time0 << " time1=" << time1
                << "\n";
      std::cerr << "line0=" << ListLines[i - 1] << "\n";
      std::cerr << "line1=" << ListLines[i] << "\n";
      throw TerminalException{1};
    }
  }
  return l_pairs;
}

double InterpolateMeasurement(std::vector<PairTimeMeas> const &ListPairTimeMeas,
                              double const &eTime,
                              double const &maxAllowedTimeInterval) {
  int siz = ListPairTimeMeas.size();
  if (siz == 0) {
    std::cerr << "We have |ListPairTimeMeas| = 0\n";
    std::cerr << "InterpolateMeasurement cannot be run correctly\n";
    throw TerminalException{1};
  }
  double epsilon = 0.00001;
  for (int i = 1; i < siz; i++) {
    double time0 = ListPairTimeMeas[i - 1].time;
    double time1 = ListPairTimeMeas[i].time;
    if (time0 - epsilon <= eTime && eTime <= time1 + epsilon) {
      double delta_time = time1 - time0;
      if (delta_time > maxAllowedTimeInterval) {
        std::cerr
            << "The time interval is too large compared to what we allow\n";
        std::cerr << "delta_time=" << delta_time
                  << " maxAllowedTimeInterval=" << maxAllowedTimeInterval
                  << "\n";
        std::cerr << "The complete list of missing intervals is:\n";
        for (int j = 1; j < siz; j++) {
          double ftime0 = ListPairTimeMeas[j - 1].time;
          double ftime1 = ListPairTimeMeas[j].time;
          double fdelta_time = ftime1 - ftime0;
          if (fdelta_time > maxAllowedTimeInterval) {
            std::string strTime0 = DATE_ConvertMjd2mystringPres(ftime0);
            std::string strTime1 = DATE_ConvertMjd2mystringPres(ftime1);
            std::cerr << "j=" << j << " delta=" << fdelta_time
                      << " time0=" << strTime0 << " time1=" << strTime1 << "\n";
          }
        }
        throw TerminalException{1};
      }
      double alpha0 = (time1 - eTime) / (time1 - time0);
      double alpha1 = (eTime - time0) / (time1 - time0);
      double meas0 = ListPairTimeMeas[i - 1].meas;
      double meas1 = ListPairTimeMeas[i].meas;
      double eValInterp = alpha0 * meas0 + alpha1 * meas1;
      return eValInterp;
    }
  }
  std::cerr << "siz=" << siz << "\n";
  double minTime = ListPairTimeMeas[0].time;
  double maxTime = ListPairTimeMeas[siz - 1].time;
  std::cerr << "Failed to find the right entry in the list of values\n";
  std::cerr << "eTime=" << eTime
            << " strPres=" << DATE_ConvertMjd2mystringPres(eTime) << "\n";
  std::cerr << "minTime=" << minTime
            << " strPres=" << DATE_ConvertMjd2mystringPres(minTime) << "\n";
  std::cerr << "maxTime=" << maxTime
            << " strPres=" << DATE_ConvertMjd2mystringPres(maxTime) << "\n";
  throw TerminalException{1};
}

// The other code

FullNamelist Individual_Tracer() {
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  // possibilities: PO4, dye_1
  ListStringValues1["TracerName"] = "unset";
  // possibilities: Constant, Monthly, Seasonal, Interpolation
  ListStringValues1["TypeVariation"] = "unset";
  ListListDoubleValues1["ListMonthlyValue"] = {};
  ListListDoubleValues1["ListSeasonalValue"] = {};
  ListStringValues1["FileInterpolation"] = "unset";
  ListDoubleValues1["ConstantValue"] = -1;
  SingleBlock BlockDESC;
  BlockDESC.setListDoubleValues(ListDoubleValues1);
  BlockDESC.setListListDoubleValues(ListListDoubleValues1);
  BlockDESC.setListStringValues(ListStringValues1);
  ListBlock["DESCRIPTION"] = BlockDESC;
  //
  return FullNamelist(ListBlock);
}

struct TracerTimeVariability {
  std::string TracerName;
  std::string TypeVariation;
  double ConstantValue;
  std::vector<double> ListMonthlyValue;
  std::vector<double> ListSeasonalValue;
  // For interpolation
  std::vector<PairTimeMeas> ListPairTimeValue;
};

TracerTimeVariability ReadIndividualTracer(FullNamelist const &eFull) {
  TracerTimeVariability ttv;
  //
  SingleBlock const &eBlDESC = eFull.get_block("DESCRIPTION");
  ttv.TracerName = eBlDESC.get_string("TracerName");
  ttv.TypeVariation = eBlDESC.get_string("TypeVariation");
  ttv.ConstantValue = eBlDESC.get_double("ConstantValue");
  ttv.ListMonthlyValue = eBlDESC.get_list_double("ListMonthlyValue");
  ttv.ListSeasonalValue = eBlDESC.get_list_double("ListSeasonalValue");
  if (ttv.TypeVariation == "Interpolation") {
    std::string FileInterpolation =
      eBlDESC.get_string("FileInterpolation");
    ttv.ListPairTimeValue = ReadFileInterpolationInformation(FileInterpolation);
  }
  if (ttv.TypeVariation == "Seasonal") {
    if (ttv.ListSeasonalValue.size() != 4) {
      eFull.NAMELIST_WriteNamelistFile(std::cerr, false);
      std::cerr << "|ttv.ListSeasonalValue|=" << ttv.ListSeasonalValue.size()
                << "\n";
      std::cerr << "ListSeasonal should have length 4 if option Seasonal is "
                   "selected\n";
      throw TerminalException{1};
    }
  }
  if (ttv.TypeVariation == "Monthly") {
    if (ttv.ListMonthlyValue.size() != 12) {
      eFull.NAMELIST_WriteNamelistFile(std::cerr, false);
      std::cerr << "|ttv.ListMonthlyValue|=" << ttv.ListMonthlyValue.size()
                << "\n";
      std::cerr << "ListSeasonal should have length 12 if option Monthly is "
                   "selected\n";
      throw TerminalException{1};
    }
  }
  return ttv;
}

double RetrieveTracerValue(TracerTimeVariability const &ttv,
                           double const &CurrentTime,
                           double const &maxAllowedTimeInterval) {
  std::vector<int> eDate = DATE_ConvertMjd2six(CurrentTime);
  if (ttv.TypeVariation == "Constant") {
    return ttv.ConstantValue;
  }
  if (ttv.TypeVariation == "Monthly") {
    int iMonth = eDate[1];
    if (iMonth < 1 || iMonth > 12) {
      std::cerr << "iMonth should be between 1 and 12\n";
      std::cerr << "The month is incorrect iMonth=" << iMonth << "\n";
      throw TerminalException{1};
    }
    return ttv.ListMonthlyValue[iMonth - 1];
  }
  if (ttv.TypeVariation == "Seasonal") {
    int iMonth = eDate[1];
    int iSeason = (iMonth - 1) / 3;
    if (iSeason < 0 || iSeason > 3) {
      std::cerr << "iSeason =" << iSeason << " and should be between 0 and 3\n";
      std::cerr << "The season is incorrect iMonth=" << iMonth << "\n";
      throw TerminalException{1};
    }
    return ttv.ListSeasonalValue[iSeason];
  }
  if (ttv.TypeVariation == "Interpolation") {
    return InterpolateMeasurement(ttv.ListPairTimeValue, CurrentTime,
                                  maxAllowedTimeInterval);
  }
  std::cerr << "Missing code for the method you choose TypeVariation = "
            << ttv.TypeVariation << "\n";
  std::cerr << "Allowed methods are Constant, Monthly, Seasonal, "
               "Interpolation\n";
  throw TerminalException{1};
}

struct TracerDescription {
  std::string namefull;
  std::string name;
  std::string long_name;
  std::string units;
  std::string field;
};

TracerDescription RetrieveTracerDescription(std::string const &eFileExternal,
                                            std::string const &TracerName) {
  std::string TracerNameFull = "river_" + TracerName;
  std::vector<std::string> ListLines = ReadFullFile(eFileExternal);
  int nbLine = ListLines.size();
  auto RemovalParenthesis = [](std::string const &estr) -> std::string {
    std::string eSep = "'";
    std::string eRET;
    int len = estr.size();
    bool IsIN = false;
    for (int i = 0; i < len; i++) {
      std::string eChar = estr.substr(i, 1);
      if (eChar == eSep) {
        IsIN = !IsIN;
      } else {
        if (IsIN)
          eRET += eChar;
      }
    }
    if (IsIN || eRET.size() == 0) {
      std::cerr << "Failed to find matching entry\n";
      std::cerr << "estr=" << estr << "\n";
      std::cerr << "IsIN=" << IsIN << "\n";
      std::cerr << "eRET=" << eRET << "\n";
      throw TerminalException{1};
    }
    return eRET;
  };
  std::vector<TracerDescription> ListMatch;
  for (int iLine = 0; iLine < nbLine; iLine++) {
    std::string eLine = ListLines[iLine];
    std::vector<std::string> LStr = STRING_Split(eLine, TracerNameFull);
    if (LStr.size() > 1) {
      if (RemovalParenthesis(ListLines[iLine]) == TracerNameFull &&
          iLine < nbLine - 4) {
        std::string eLine_longname = ListLines[iLine + 1];
        std::string eLine_units = ListLines[iLine + 2];
        std::string eLine_field = ListLines[iLine + 3];
        TracerDescription eTracer{
            TracerNameFull, TracerName, RemovalParenthesis(eLine_longname),
            RemovalParenthesis(eLine_units), RemovalParenthesis(eLine_field)};
        ListMatch.push_back(eTracer);
      }
    }
  }
  if (ListMatch.size() > 1) {
    std::cerr << "|ListMatch|=" << ListMatch.size() << "\n";
    std::cerr << "We found several matching entries for TracerName = "
              << TracerName << "\n";
    std::cerr << "eFileExternal = " << eFileExternal << "\n";
    for (auto &etd : ListMatch) {
      std::cerr << "name=" << etd.name << " long_name=" << etd.long_name
                << " units=" << etd.units << " field=" << etd.field << "\n";
    }
    throw TerminalException{1};
  }
  if (ListMatch.size() == 0) {
    std::cerr << "Failed to find the Tracer information for TracerName = "
              << TracerName << "\n";
    std::cerr << "eFileExternal = " << eFileExternal << "\n";
    throw TerminalException{1};
  }
  return ListMatch[0];
}

FullNamelist Individual_River_File() {
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["name"] = "";
  ListStringValues1["TypeVaryingTransport"] =
      "Please select ConstantFlux, MonthlyFlux, InterpolationFlux";
  ListStringValues1["TypeVaryingTemperature"] =
      "Please select ConstantTemp, MonthlyTemp, InterpolationTemp";
  ListStringValues1["TypeVaryingSalinity"] = "Only possibilities, ConstantSalt";
  ListStringValues1["verticalShapeOption"] = "UpperLayer";
  ListListDoubleValues1["ListMonthlyFlux"] = {};
  ListListDoubleValues1["ListMonthlyTemp"] = {};
  ListListDoubleValues1["ListMonthlySalt"] = {};
  ListDoubleValues1["ConstantFlux"] = -1;
  ListDoubleValues1["ConstantFactorFlux"] = 1.0;
  ListStringValues1["FileRiverFlux"] = "unset";
  ListStringValues1["FileRiverTemp"] = "unset";
  ListStringValues1["FileRiverSalt"] = "unset";
  ListStringValues1["WScase"] = "River";
  ListBoolValues1["SetRiverTemperature"] = true;
  ListBoolValues1["SetRiverSalinity"] = true;
  ListDoubleValues1["ConstantRiverTemperature"] = 14;
  ListDoubleValues1["ConstantRiverSalinity"] = 0;
  ListDoubleValues1["lon"] = -400;
  ListDoubleValues1["lat"] = -400;
  ListDoubleValues1["direction"] = -400;
  ListDoubleValues1["MaxDepth"] = 2;
  ListDoubleValues1["targetDepth"] = 2;
  ListDoubleValues1["FrequencyDay"] = 7;
  ListDoubleValues1["DurationHour"] = 2;
  ListDoubleValues1["TotalFlux"] = 15300;
  ListListStringValues1["ListTracerFile"] = {};
  ListStringValues1["PrefixTracer"] = {};
  ListIntValues1["iSelect"] = -1;
  ListIntValues1["jSelect"] = -1;
  ListIntValues1["SignSelect"] = -400;
  ListIntValues1["DirSelect"] = -1;
  ListIntValues1["ChoiceSelect"] = -1;
  SingleBlock BlockDESC;
  BlockDESC.setListIntValues(ListIntValues1);
  BlockDESC.setListBoolValues(ListBoolValues1);
  BlockDESC.setListDoubleValues(ListDoubleValues1);
  BlockDESC.setListListDoubleValues(ListListDoubleValues1);
  BlockDESC.setListListIntValues(ListListIntValues1);
  BlockDESC.setListStringValues(ListStringValues1);
  BlockDESC.setListListStringValues(ListListStringValues1);
  ListBlock["DESCRIPTION"] = BlockDESC;
  //
  return FullNamelist(ListBlock);
}

struct DescriptionRiver {
  double lon;
  double lat;
  double direction;
  std::string verticalShapeOption;
  double MaxDepth;
  double targetDepth;
  double ConstantRiverTemperature;
  double ConstantRiverSalinity;
  bool SetRiverTemperature;
  bool SetRiverSalinity;
  std::vector<double> ListMonthlyFlux;
  std::vector<double> ListMonthlyTemp;
  std::vector<double> ListMonthlySalt;
  double ConstantFlux;
  double ConstantFactorFlux;
  std::vector<PairTimeMeas> ListPairTimeFlux;
  std::vector<PairTimeMeas> ListPairTimeTemp;
  std::vector<PairTimeMeas> ListPairTimeSalt;
  std::string TypeVaryingTransport;
  std::string TypeVaryingTemperature;
  std::string TypeVaryingSalinity;
  std::string name;
  double FrequencyDay;
  double DurationHour;
  double TotalFlux;
  std::string WScase;
  int iSelect;
  int jSelect;
  int SignSelect;
  int DirSelect;
  int ChoiceSelect;
  std::map<std::string, TracerTimeVariability> MapTracerDesc;
};

DescriptionRiver ReadRiverDescription(std::string const &RiverDescriptionFile) {
  //  std::cerr << "ReadRiverDescription, step 1\n";
  FullNamelist eFull = Individual_River_File();
  //  std::cerr << "ReadRiverDescription, step 2\n";
  NAMELIST_ReadNamelistFile(RiverDescriptionFile, eFull);
  //  std::cerr << "ReadRiverDescription, step 3\n";
  SingleBlock eBlDESC = eFull.get_block("DESCRIPTION");
  //  std::cerr << "ReadRiverDescription, step 4\n";
  DescriptionRiver eDesc;
  //  std::cerr << "ReadRiverDescription, step 5\n";
  eDesc.lon = eBlDESC.get_double("lon");
  eDesc.lat = eBlDESC.get_double("lat");
  eDesc.direction = eBlDESC.get_double("direction");
  eDesc.MaxDepth = eBlDESC.get_double("MaxDepth");
  eDesc.targetDepth = eBlDESC.get_double("targetDepth");
  eDesc.SetRiverTemperature = eBlDESC.get_bool("SetRiverTemperature");
  eDesc.SetRiverSalinity = eBlDESC.get_bool("SetRiverSalinity");
  eDesc.ConstantRiverTemperature =
    eBlDESC.get_double("ConstantRiverTemperature");
  eDesc.ConstantRiverSalinity =
    eBlDESC.get_double("ConstantRiverSalinity");
  eDesc.TypeVaryingTransport =
    eBlDESC.get_string("TypeVaryingTransport");
  eDesc.TypeVaryingTemperature =
    eBlDESC.get_string("TypeVaryingTemperature");
  eDesc.TypeVaryingSalinity =
    eBlDESC.get_string("TypeVaryingSalinity");
  eDesc.verticalShapeOption =
    eBlDESC.get_string("verticalShapeOption");
  eDesc.FrequencyDay = eBlDESC.get_double("FrequencyDay");
  eDesc.DurationHour = eBlDESC.get_double("DurationHour");
  eDesc.TotalFlux = eBlDESC.get_double("TotalFlux");
  eDesc.WScase = eBlDESC.get_string("WScase");
  eDesc.name = eBlDESC.get_string("name");
  eDesc.ListMonthlyFlux = eBlDESC.get_list_double("ListMonthlyFlux");
  eDesc.ListMonthlyTemp = eBlDESC.get_list_double("ListMonthlyTemp");
  eDesc.ListMonthlySalt = eBlDESC.get_list_double("ListMonthlySalt");
  eDesc.ConstantFlux = eBlDESC.get_double("ConstantFlux");
  eDesc.ConstantFactorFlux = eBlDESC.get_double("ConstantFactorFlux");
  //  std::cerr << "ReadRiverDescription, step 6\n";
  auto CheckListMonth = [&](std::vector<double> const &ListMon) -> void {
    int nbMonth = ListMon.size();
    if (nbMonth != 12 && nbMonth != 0) {
      std::cerr << "RiverDescriptionFile = " << RiverDescriptionFile << "\n";
      std::cerr << "We have nbMonth=" << nbMonth << "\n";
      std::cerr << "It should be 12 or not set\n";
      throw TerminalException{1};
    }
  };
  CheckListMonth(eDesc.ListMonthlyFlux);
  CheckListMonth(eDesc.ListMonthlyTemp);
  CheckListMonth(eDesc.ListMonthlySalt);
  if (eDesc.TypeVaryingTransport == "InterpolationFlux") {
    std::string FileRiverFlux = eBlDESC.get_string("FileRiverFlux");
    eDesc.ListPairTimeFlux = ReadFileInterpolationInformation(FileRiverFlux);
  }
  if (eDesc.TypeVaryingTemperature == "InterpolationTemp") {
    std::string FileRiverTemp = eBlDESC.get_string("FileRiverTemp");
    eDesc.ListPairTimeTemp = ReadFileInterpolationInformation(FileRiverTemp);
  }
  if (eDesc.TypeVaryingTemperature == "InterpolationSalt") {
    std::string FileRiverSalt = eBlDESC.get_string("FileRiverSalt");
    eDesc.ListPairTimeSalt = ReadFileInterpolationInformation(FileRiverSalt);
  }
  eDesc.iSelect = eBlDESC.get_int("iSelect");
  eDesc.jSelect = eBlDESC.get_int("jSelect");
  eDesc.SignSelect = eBlDESC.get_int("SignSelect");
  eDesc.DirSelect = eBlDESC.get_int("DirSelect");
  eDesc.ChoiceSelect = eBlDESC.get_int("ChoiceSelect");
  //
  // The additional tracers
  //
  std::vector<std::string> ListTracerFile =
    eBlDESC.get_list_string("ListTracerFile");
  std::string PrefixTracer = eBlDESC.get_string("PrefixTracer");
  std::map<std::string, TracerTimeVariability> MapTracerDesc;
  for (auto &eTracerFile : ListTracerFile) {
    FullNamelist eFullTracer = Individual_Tracer();
    std::string FullFile = PrefixTracer + eTracerFile;
    NAMELIST_ReadNamelistFile(FullFile, eFullTracer);
    TracerTimeVariability ttv = ReadIndividualTracer(eFullTracer);
    MapTracerDesc[ttv.TracerName] = ttv;
  }
  eDesc.MapTracerDesc = MapTracerDesc;
  return eDesc;
}

FullNamelist NAMELIST_PLOT_River() {
  std::map<std::string, SingleBlock> ListBlock;
  //
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["GridFile"] = "unset";
  ListStringValues1["RiverFile"] = "unset";
  ListStringValues1["SVGfile"] = "unset";
  ListStringValues1["PicPrefix"] = "unset";
  ListStringValues1["Extension"] = "png";
  ListStringValues1["BEGTC"] = "20110915.000000";
  ListStringValues1["ENDTC"] = "20110915.000000";
  ListStringValues1["__NaturePlot"] = "RIVER";
  ListBoolValues1["FirstCleanDirectory"] = true;
  ListBoolValues1["KeepNC_NCL"] = false;
  ListBoolValues1["InPlaceRun"] = false;
  ListBoolValues1["PrintDebugInfo"] = false;
  ListBoolValues1["OnlyCreateFiles"] = false;
  ListIntValues1["NPROC"] = 1;
  ListStringValues1["Pcolor_method"] = "ncl";
  ListStringValues1["Quiver_method"] = "ncl";
  ListStringValues1["Lines_method"] = "ncl";
  ListStringValues1["Scatter_method"] = "ncl";
  SingleBlock BlockPROC;
  BlockPROC.setListIntValues(ListIntValues1);
  BlockPROC.setListBoolValues(ListBoolValues1);
  BlockPROC.setListDoubleValues(ListDoubleValues1);
  BlockPROC.setListListDoubleValues(ListListDoubleValues1);
  BlockPROC.setListListIntValues(ListListIntValues1);
  BlockPROC.setListStringValues(ListStringValues1);
  BlockPROC.setListListStringValues(ListListStringValues1);
  ListBlock["PROC"] = BlockPROC;
  //
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::vector<double>> ListListDoubleValues2;
  std::map<std::string, std::vector<int>> ListListIntValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["PlotFlux"] = false;
  ListStringValues2["BEGTC"] = "20090915.000000";
  ListStringValues2["ENDTC"] = "20110915.000000";
  ListDoubleValues2["WidthLine"] = 0.4;
  ListDoubleValues2["WidthLineRiver"] = 0.4;
  ListIntValues2["SizeText"] = 2;
  ListIntValues2["nbLabel"] = 5;
  ListStringValues2["StyleDate"] = "unset";
  ListStringValues2["StyleTitle"] = "unset";
  ListStringValues2["TitleString"] = "unset";
  ListStringValues2["StyleTitle"] = "unset";
  ListDoubleValues2["ShiftLonText"] = 0.05;
  SingleBlock BlockPLOT;
  BlockPLOT.setListIntValues(ListIntValues2);
  BlockPLOT.setListBoolValues(ListBoolValues2);
  BlockPLOT.setListDoubleValues(ListDoubleValues2);
  BlockPLOT.setListListDoubleValues(ListListDoubleValues2);
  BlockPLOT.setListListIntValues(ListListIntValues2);
  BlockPLOT.setListStringValues(ListStringValues2);
  BlockPLOT.setListListStringValues(ListListStringValues2);
  ListBlock["PLOT"] = BlockPLOT;
  //
  return FullNamelist(ListBlock);
}

struct ijSeaLand {
  int iSea;
  int jSea;
  int iLand;
  int jLand;
};

struct ijdsInfo {
  int eStatus;
  int iSelect;
  int jSelect;
  int DirSelect;
  int SignSelect;
};

ijSeaLand GetArrayIjSeaLand(ijdsInfo const &RecInfo) {
  int iSelect = RecInfo.iSelect;
  int jSelect = RecInfo.jSelect;
  int DirSelect = RecInfo.DirSelect;
  int SignSelect = RecInfo.SignSelect;
  int iSea = -1, jSea = -1, iLand = -1, jLand = -1;
  if (DirSelect == 0 && SignSelect == 1) {
    // Right case in map_rivers.m
    iSea = iSelect + 1;
    jSea = jSelect + 1;
    iLand = iSea;
    jLand = jSea - 1;
  }
  if (DirSelect == 1 && SignSelect == 1) {
    // Up case in map_rivers.m
    iSea = iSelect + 1;
    jSea = jSelect + 1;
    iLand = iSea - 1;
    jLand = jSea;
  }
  if (DirSelect == 0 && SignSelect == -1) {
    // Left case in map_rivers.m
    iSea = iSelect + 1;
    jSea = jSelect;
    iLand = iSea;
    jLand = jSea + 1;
  }
  if (DirSelect == 1 && SignSelect == -1) {
    // Down case in map_rivers.m
    iSea = iSelect;
    jSea = jSelect + 1;
    iLand = iSea + 1;
    jLand = jSea;
  }
  return {iSea - 1, jSea - 1, iLand - 1, jLand - 1};
}

/*
  Choice=0 means right
  Choice=1 means up
  Choice=2 means left
  Choice=3 means down
 */
ijdsInfo RetrieveIJDSarray(int const &eEtaSea, int const &eXiSea,
                           int const &iChoice) {
  int iSelect = -1, jSelect = -1, DirSelect = -1, SignSelect = -1;
  if (iChoice == 0) {
    iSelect = eEtaSea - 1;
    jSelect = eXiSea - 1;
    DirSelect = 0;
    SignSelect = 1;
  }
  if (iChoice == 1) {
    iSelect = eEtaSea - 1;
    jSelect = eXiSea - 1;
    DirSelect = 1;
    SignSelect = 1;
  }
  if (iChoice == 2) {
    iSelect = eEtaSea - 1;
    jSelect = eXiSea;
    DirSelect = 0;
    SignSelect = -1;
  }
  if (iChoice == 3) {
    iSelect = eEtaSea;
    jSelect = eXiSea - 1;
    DirSelect = 1;
    SignSelect = -1;
  }
  if (iChoice == -1) {
    std::cerr << "RetrieveIJDSarray INVALID ICHOICE iChoice=-1\n";
    return {0, -1, -1, -1, -1};
  }
  // increment conversion between matlab and C++ indexing:
  iSelect++;
  jSelect++;
  return {1, iSelect, jSelect, DirSelect, SignSelect};
}

template <typename T> T AverageValue(MyMatrix<T> const &M) {
  return M.sum() / M.size();
}

template <typename T> T AverageValue(MyVector<T> const &M) {
  return M.sum() / M.size();
}

void PlotRiverInformation(FullNamelist const &eFull) {
  SingleBlock const& eBlPLOT = eFull.get_block("PLOT");
  SingleBlock const& eBlPROC = eFull.get_block("PROC");
  //
  std::string GridFile = eBlPROC.get_string("GridFile");
  std::string RiverFile = eBlPROC.get_string("RiverFile");
  std::string SVGfile = eBlPROC.get_string("SVGfile");
  //
  double SizeLine = eBlPLOT.get_double("WidthLine");
  double SizeLineRiver = eBlPLOT.get_double("WidthLineRiver");
  int SizeText = eBlPLOT.get_int("SizeText");
  double ShiftLonText = eBlPLOT.get_double("ShiftLonText");
  //
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  MyVector<double> ListETA_v = NC_Read1Dvariable(RiverFile, "river_Eposition");
  std::cerr << "We have ListETA_v\n";
  MyVector<double> ListXI_v = NC_Read1Dvariable(RiverFile, "river_Xposition");
  std::cerr << "We have ListXI_v\n";
  MyVector<double> ListDir_v = NC_Read1Dvariable(RiverFile, "river_direction");
  std::cerr << "We have ListDir_v\n";
  Eigen::Tensor<double, 3> DATA_Salt =
      NC_Read3Dvariable(RiverFile, "river_salt");
  std::cerr << "We have DATA_Salt\n";
  Eigen::Tensor<double, 3> DATA_Temp =
      NC_Read3Dvariable(RiverFile, "river_temp");
  std::cerr << "We have DATA_Temp\n";
  std::vector<Eigen::Tensor<double, 3>> ListDATA_dye;
  int iDye = 0;
  while (true) {
    std::string eVar = "river_dye_" + StringNumber(iDye + 1, 2);
    if (!NC_IsVar(RiverFile, eVar))
      break;
    //
    Eigen::Tensor<double, 3> DATA_dye = NC_Read3Dvariable(RiverFile, eVar);
    std::cerr << "We have DATA_dye\n";
    ListDATA_dye.emplace_back(DATA_dye);
    iDye++;
  }
  MyMatrix<double> MatTransport =
      NC_Read2Dvariable(RiverFile, "river_transport");
  std::vector<double> ListRiverTime =
      NC_ReadTimeFromFile(RiverFile, "river_time");
  int nbRiver = ListETA_v.size();
  //
  PermanentInfoDrawing ePerm = GET_PERMANENT_INFO(eFull);
  NCLcaller<GeneralType> eCall(ePerm.NPROC);
  //
  std::vector<std::string> ListRiverName;
  try {
    netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::read);
    netCDF::NcGroupAtt RiversAtt = dataFile.getAtt("rivers");
    std::string strRiverName;
    RiversAtt.getValues(strRiverName);
    ListRiverName = STRING_Split(strRiverName, ",");
  } catch (...) {
    for (int i = 0; i < nbRiver; i++)
      ListRiverName.push_back(IntToString(i));
  }
  //  std::cerr << "strRiverName = " << strRiverName << "\n";
  //
  int nbTime = MatTransport.rows();
  std::cerr << "MatTransport Rows=" << MatTransport.rows()
            << " Cols=" << MatTransport.cols() << "\n";
  std::vector<int> ListETA(nbRiver), ListXI(nbRiver), ListDir(nbRiver),
      ListSign(nbRiver);
  std::vector<double> ListAvgFlux(nbRiver);
  for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
    ListETA[iRiver] = ListETA_v(iRiver);
    ListXI[iRiver] = ListXI_v(iRiver);
    ListDir[iRiver] = ListDir_v(iRiver);
    double TheSum = 0;
    for (int iTime = 0; iTime < nbTime; iTime++)
      TheSum += MatTransport(iTime, iRiver);
    double TheAvgFlux = TheSum / static_cast<double>(nbTime);
    int eSign;
    if (TheSum > 0)
      eSign = 1;
    else
      eSign = -1;
    ListSign[iRiver] = eSign;
    ListAvgFlux[iRiver] = TheAvgFlux;
  }
  //
  // Some standard definitions
  //
  MyMatrix<int> MatDir(4, 2);
  MatDir(0, 0) = 0;
  MatDir(0, 1) = 0;
  MatDir(1, 0) = 1;
  MatDir(1, 1) = 0;
  MatDir(2, 0) = 1;
  MatDir(2, 1) = 1;
  MatDir(3, 0) = 0;
  MatDir(3, 1) = 1;
  CoordGridArrayFD RecPsi2 = GRID_ExtendedPsiThi(
      GrdArr.GrdArrRho, GrdArr.GrdArrU, GrdArr.GrdArrV, GrdArr.GrdArrPsi);
  int eta_psi2 = RecPsi2.LON.rows();
  int xi_psi2 = RecPsi2.LON.cols();
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  MyMatrix<double> MatRadius = GetMatrixRadiusROMS(GrdArr);
  //
  // The list of entries to be printed.
  //
  std::vector<SVGgeneral> ListGeneral;
  //
  // Computing the coloring of the points
  //
  bool PolylineColoring = true;
  if (PolylineColoring) {
    for (int i = 0; i < eta_rho; i++)
      for (int j = 0; j < xi_rho; j++) {
        std::vector<coor> ListCoor(4);
        for (int u = 0; u < 4; u++) {
          int i2 = i + MatDir(u, 0);
          int j2 = j + MatDir(u, 1);
          coor c{RecPsi2.LON(i2, j2), RecPsi2.LAT(i2, j2)};
          ListCoor[u] = c;
        }
        std::string MarkerEnd = "";
        std::string clip = "full";
        double Size = 0;
        std::vector<int> colorstroke{0, 0, 0};
        std::vector<int> colorfill;
        if (GrdArr.GrdArrRho.MSK(i, j) == 0)
          colorfill = {255, 0, 0};
        else
          colorfill = {0, 0, 255};
        SVGqualInfoPolyline eQual{colorfill, colorstroke, Size, MarkerEnd,
                                  clip};
        SVGpolyline ePolyline{ListCoor, eQual};
        ListGeneral.push_back(SVGgeneral(ePolyline));
      }
    std::cerr << "ListPolyline inserted\n";
  }
  //
  // Inserting the lines of the grid
  //
  auto InsertGridLine = [&](int const &i1, int const &j1, int const &i2,
                            int const &j2) -> void {
    //    std::cerr << "InsertGridLine,  step 1, i1=" << i1 << " j1=" << j1 << "
    //    i2=" << i2 << " j2=" << j2 << "\n";
    coor ePt{RecPsi2.LON(i1, j1), RecPsi2.LAT(i1, j1)};
    //    std::cerr << "InsertGridLine,  step 2\n";
    coor fPt{RecPsi2.LON(i2, j2), RecPsi2.LAT(i2, j2)};
    //    std::cerr << "InsertGridLine,  step 3\n";
    std::string MarkerEnd = "";
    std::string clip = "";
    double Size = SizeLine;
    std::vector<int> color{0, 0, 0};
    SVGqualInfo eQual{color, Size, MarkerEnd, clip};
    SVGline eLine{ePt, fPt, eQual};
    ListGeneral.push_back(SVGgeneral(eLine));
  };
  for (int i = 0; i < eta_psi2 - 1; i++)
    for (int j = 0; j < xi_psi2; j++)
      InsertGridLine(i, j, i + 1, j);
  std::cerr << "Block 1 inserted\n";
  for (int i = 0; i < eta_psi2; i++)
    for (int j = 0; j < xi_psi2 - 1; j++)
      InsertGridLine(i, j, i, j + 1);
  std::cerr << "Block 2 inserted\n";
  //
  // Computing the river information
  //
  bool PlotRiver = true;
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  if (PlotRiver) {
    for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
      int iSelect = ListETA[iRiver];
      int jSelect = ListXI[iRiver];
      int DirSelect = ListDir[iRiver];
      int SignSelect = ListSign[iRiver];
      ijSeaLand recIJSL =
          GetArrayIjSeaLand({1, iSelect, jSelect, DirSelect, SignSelect});
      for (int u = 0; u < 2; u++) {
        std::vector<int> color;
        int i, j;
        if (u == 0) {
          color = {255, 0, 0};
          i = recIJSL.iSea;
          j = recIJSL.jSea;
        } else {
          color = {0, 0, 255};
          i = recIJSL.iLand;
          j = recIJSL.jLand;
        }
        double coef = 0.3;
        double eRad = coef * MatRadius(i, j);
        coor r{eRad, eRad};
        double lon = GrdArr.GrdArrRho.LON(i, j);
        double lat = GrdArr.GrdArrRho.LAT(i, j);
        coor c{lon, lat};
        std::string MarkerEnd = "";
        std::string clip = "";
        double Size = SizeLineRiver;
        SVGqualInfo eQual{color, Size, MarkerEnd, clip};
        SVGellipse eEll{c, r, eQual};
        ListGeneral.push_back(SVGgeneral(eEll));
      }
      int iSea = recIJSL.iSea;
      int jSea = recIJSL.jSea;
      int iLand = recIJSL.iLand;
      int jLand = recIJSL.jLand;
      double lonSea = GrdArr.GrdArrRho.LON(iSea, jSea);
      double latSea = GrdArr.GrdArrRho.LAT(iSea, jSea);
      double lonLand = GrdArr.GrdArrRho.LON(iLand, jLand);
      double latLand = GrdArr.GrdArrRho.LAT(iLand, jLand);
      double lonMid = (lonSea + lonLand) / 2;
      double latMid = (latSea + latLand) / 2;
      coor ePt{lonSea, latSea};
      coor fPt{lonLand, latLand};
      int eMSKsea = GrdArr.GrdArrRho.MSK(iSea, jSea);
      int eMSKland = GrdArr.GrdArrRho.MSK(iLand, jLand);
      if (eMSKsea != 1 || eMSKland != 0) {
        std::cerr << "INCONSISTENCY iRiver=" << iRiver
                  << " has eMSKsea=" << eMSKsea << " eMSKland=" << eMSKland
                  << "\n";
      }
      int eDEPsea = DEP(iSea, jSea);
      std::cerr << "iRiver=" << iRiver << " Land(lon/lat)=" << lonLand << " / "
                << latLand << " SeaLand(Lon/Lat)=" << lonMid << " / " << latMid
                << "\n";
      std::cerr << "        Sea(lon/lat/dep)=" << lonSea << " / " << latSea
                << " / " << eDEPsea << "\n";
      std::cerr << "avgFlux=" << ListAvgFlux[iRiver] << "\n";
      MyVector<double> VectTrans = GetMatrixCol(MatTransport, iRiver);
      MyMatrix<double> ArrSalt = DimensionExtraction(DATA_Salt, 2, iRiver);
      MyMatrix<double> ArrTemp = DimensionExtraction(DATA_Temp, 2, iRiver);
      std::string str1 = std::string("  transport: ") +
                         std::to_string(AverageValue(VectTrans));
      std::string str2 =
          std::string("salt: ") + std::to_string(AverageValue(ArrSalt));
      std::string str3 =
          std::string("temp: ") + std::to_string(AverageValue(ArrTemp));
      std::string strO = str1 + " " + str2 + " " + str3;
      int iDye = 0;
      for (auto &eDATA : ListDATA_dye) {
        MyMatrix<double> eArr = DimensionExtraction(eDATA, 2, iRiver);
        strO += std::string(" dye") + StringNumber(iDye + 1, 2) +
                std::string(": ") + std::to_string(AverageValue(eArr));
      }
      std::cerr << strO << "\n";
      //
      {
        std::string MarkerEnd = "";
        std::string clip = "";
        double Size = SizeLineRiver;
        std::vector<int> color{255, 255, 0};
        SVGqualInfo eQual{color, Size, MarkerEnd, clip};
        SVGline eLine{ePt, fPt, eQual};
        ListGeneral.push_back(SVGgeneral(eLine));
      }
      //
      {
        double lonLnd = GrdArr.GrdArrRho.LON(iLand, jLand);
        double latLnd = GrdArr.GrdArrRho.LAT(iLand, jLand);
        double eLon = (lonSea + lonLnd) / 2;
        double eLat = (latSea + latLnd) / 2;
        coor point{eLon + ShiftLonText, eLat};
        std::vector<int> color{0, 255, 255};
        std::string str =
            ListRiverName[iRiver] + " : " + DoubleToString(ListAvgFlux[iRiver]);
        SVGtext eText{point, str, color, SizeText};
        ListGeneral.push_back(SVGgeneral(eText));
      }
    }
    std::cerr << "River lines inserted\n";
  }
  //
  // Plot the flux
  //
  bool PlotFlux = eBlPLOT.get_bool("PlotFlux");
  std::cerr << "PlotFlux=" << PlotFlux << "\n";
  if (PlotFlux) {
    std::string strBEGTC = eBlPROC.get_string("BEGTC");
    std::string strENDTC = eBlPROC.get_string("ENDTC");
    double BeginTime = CT2MJD(strBEGTC);
    double EndTime = CT2MJD(strENDTC);
    std::cerr << "BeginTime=" << BeginTime << " EndTime=" << EndTime << "\n";
    std::vector<int> ListIdx;
    int nbTime = ListRiverTime.size();
    for (int iTime = 0; iTime < nbTime; iTime++) {
      double eTime = ListRiverTime[iTime];
      if (BeginTime <= eTime && eTime <= EndTime)
        ListIdx.push_back(iTime);
    }
    int nbCorr = ListIdx.size();
    std::cerr << "nbTime=" << nbTime << " nbCorr=" << nbCorr << "\n";
    MyVector<double> ListRiverTime_sel(nbCorr);
    for (int iCorr = 0; iCorr < nbCorr; iCorr++) {
      int idx = ListIdx[iCorr];
      ListRiverTime_sel(iCorr) = ListRiverTime[idx];
    }
    std::cerr << "nbRiver=" << nbRiver << "\n";
    int nbLabel = eBlPLOT.get_int("nbLabel");
    std::string StyleDate = eBlPLOT.get_string("StyleDate");
    std::string StyleTitle = eBlPLOT.get_string("StyleTitle");
    for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
      MyVector<double> ListTransport_sel(nbCorr);
      for (int iCorr = 0; iCorr < nbCorr; iCorr++) {
        int idx = ListIdx[iCorr];
        ListTransport_sel(iCorr) = MatTransport(idx, iRiver);
      }
      double TheMax = ListTransport_sel.maxCoeff();
      DrawLinesArr eDrawArr;
      eDrawArr.DoTitle = true;
      if (StyleTitle == "style0")
        eDrawArr.TitleStr =
            "River " + IntToString(iRiver) + " name=" + ListRiverName[iRiver];
      if (StyleTitle == "style1")
        eDrawArr.TitleStr = ListRiverName[iRiver];
      if (StyleTitle == "style2")
        eDrawArr.TitleStr = eBlPLOT.get_string("TitleString");
      eDrawArr.IsTimeSeries = true;
      eDrawArr.PairComparison = false;
      eDrawArr.DoExplicitLabel = true;
      eDrawArr.DrawHorizVertLines = false;
      eDrawArr.nbLabel = nbLabel;
      eDrawArr.StyleDate = StyleDate;
      eDrawArr.VarName = "River_flux";
      eDrawArr.TheMax = TheMax;
      eDrawArr.TheMin = 0;
      eDrawArr.ListX = ListRiverTime_sel;
      eDrawArr.ListListVect = {ListTransport_sel};
      eDrawArr.YAxisString = "River flux [m3/s]";
      //      eDrawArr.ListName_plot={"flux"};
      std::string FileName = ePerm.eDir + "FLUX_" + IntToString(iRiver);
      LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
    }
  }
  //
  // Computing the river information
  //
  SVGplotDescription SVGplot;
  SVGplot.ListGeneral = ListGeneral;
  SVGplot.FrameOption = 1;
  SVGplot.height = 400;
  SVGplot.width = 400;
  SVGplot.scale_factor = 1;
  SVGplot.add_offsetX = 0;
  SVGplot.add_offsetY = 0;
  SVGplot.RoundMethod = 2;
  GeneralWriteSVGfile(SVGfile, SVGplot);
}

FullNamelist NAMELIST_GetStandard_ComputeRiverForcing_ROMS() {
  std::map<std::string, SingleBlock> ListBlock;
  // INPUT
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListListStringValues1["ListRiverName"] = {};
  ListStringValues1["RiverPrefix"] = "unset";
  ListStringValues1["RiverSuffix"] = "unset";
  ListStringValues1["GridFile"] = "UNK";
  ListStringValues1["BEGTC"] = "20110915.000000";
  ListStringValues1["ENDTC"] = "20110925.000000";
  ListStringValues1["RefTime"] = "19680523.000000";
  ListDoubleValues1["DELTC"] = 600;
  ListDoubleValues1["maxAllowedTimeInterval"] = 10000;
  ListStringValues1["UNITC"] = "SEC";
  ListIntValues1["ARVD_N"] = -1;
  ListIntValues1["ARVD_Vtransform"] = -1;
  ListIntValues1["ARVD_Vstretching"] = -1;
  ListDoubleValues1["ARVD_Tcline"] = -1;
  ListDoubleValues1["ARVD_hc"] = -1;
  ListDoubleValues1["ARVD_theta_s"] = -1;
  ListDoubleValues1["ARVD_theta_b"] = -1;
  ListStringValues1["RiverFile"] = "unset.nc";
  ListStringValues1["ExternalInfoFile"] = "unset";
  ListListStringValues1["ListAdditionalTracers"] = {};
  SingleBlock BlockINPUT;
  BlockINPUT.setListIntValues(ListIntValues1);
  BlockINPUT.setListBoolValues(ListBoolValues1);
  BlockINPUT.setListDoubleValues(ListDoubleValues1);
  BlockINPUT.setListListDoubleValues(ListListDoubleValues1);
  BlockINPUT.setListListIntValues(ListListIntValues1);
  BlockINPUT.setListStringValues(ListStringValues1);
  BlockINPUT.setListListStringValues(ListListStringValues1);
  ListBlock["INPUT"] = BlockINPUT;
  //
  return FullNamelist(ListBlock);
}

FullNamelist NAMELIST_RetrieveData() {
  std::map<std::string, SingleBlock> ListBlock;
  // INPUT
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::vector<double>> ListListDoubleValues1;
  std::map<std::string, std::vector<int>> ListListIntValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["RiverDescriptionFile"] = "unset";
  ListListStringValues1["ListTimes"] = {};
  ListIntValues1["StylePrint"] = 1;
  SingleBlock BlockINPUT;
  BlockINPUT.setListIntValues(ListIntValues1);
  BlockINPUT.setListBoolValues(ListBoolValues1);
  BlockINPUT.setListDoubleValues(ListDoubleValues1);
  BlockINPUT.setListListDoubleValues(ListListDoubleValues1);
  BlockINPUT.setListListIntValues(ListListIntValues1);
  BlockINPUT.setListStringValues(ListStringValues1);
  BlockINPUT.setListListStringValues(ListListStringValues1);
  ListBlock["INPUT"] = BlockINPUT;
  //
  return FullNamelist(ListBlock);
}

double MonthlyInterpolation(std::vector<double> const &ListMonthly,
                            double const &eTimeDay) {
  std::vector<int> LDate = DATE_ConvertMjd2six(eTimeDay);
  int eYear = LDate[0];
  int eMonth = LDate[1];
  double eTimeDayMid = DATE_ConvertSix2mjd({eYear, eMonth, 15, 0, 0, 0});
  double eTimeDayMidSEQapprox;
  if (eTimeDay < eTimeDayMid) {
    eTimeDayMidSEQapprox = eTimeDayMid - 30;
  } else {
    eTimeDayMidSEQapprox = eTimeDayMid + 30;
  }
  std::vector<int> LDate2 = DATE_ConvertMjd2six(eTimeDayMidSEQapprox);
  int eYear2 = LDate2[0];
  int eMonth2 = LDate2[1];
  double eTimeDayMidSEQ = DATE_ConvertSix2mjd({eYear2, eMonth2, 15, 0, 0, 0});
  //
  double alpha2 = (eTimeDay - eTimeDayMidSEQ) / (eTimeDayMid - eTimeDayMidSEQ);
  double alpha1 = (eTimeDayMid - eTimeDay) / (eTimeDayMid - eTimeDayMidSEQ);
  double epsilon = 0.0001;
  if (alpha1 < -epsilon || alpha2 < -epsilon) {
    std::cerr << "--------------------------------\n";
    std::string strPres = DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "eTimeDay=" << eTimeDay << " date=" << strPres << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMid = DATE_ConvertMjd2mystringPres(eTimeDayMid);
    std::cerr << "eTimeDayMid=" << eTimeDayMid << " date=" << strPresMid
              << "\n";
    std::cerr << "eYear=" << eYear << " eMonth=" << eMonth << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMidSEQ = DATE_ConvertMjd2mystringPres(eTimeDayMidSEQ);
    std::cerr << "eTimeDayMidSEQ=" << eTimeDayMidSEQ
              << " date=" << strPresMidSEQ << "\n";
    std::cerr << "eYear2=" << eYear2 << " eMonth2=" << eMonth2 << "\n";
    std::cerr << "--------------------------------\n";
    std::string strPresMidSEQapprox =
        DATE_ConvertMjd2mystringPres(eTimeDayMidSEQapprox);
    std::cerr << "eTimeDayMidSEQapprox=" << eTimeDayMidSEQapprox
              << " date=" << strPresMidSEQapprox << "\n";
    std::cerr << "--------------------------------\n";
    std::cerr << "eTimeDayMid=" << eTimeDayMid
              << " eTimeDayMidSEQ=" << eTimeDayMidSEQ << "\n";
    std::cerr << "alpha1=" << alpha1 << " alpha2=" << alpha2 << "\n";
    std::cerr << "Error on alpha1 / alpha2\n";
    throw TerminalException{1};
  }
  double Flux1 = ListMonthly[eMonth2 - 1];
  double Flux2 = ListMonthly[eMonth - 1];
  return alpha1 * Flux1 + alpha2 * Flux2;
}

struct TransTempSalt {
  double eTransport;
  double eTemp;
  double eSalt;
};

TransTempSalt RetrieveTTS(DescriptionRiver const &eDescRiv,
                          double const &eTimeDay,
                          double const &maxAllowedTimeInterval) {
  double eTransport = 0, eTemp = 0, eSalt = 0;
  bool HasTransport = false, HasTemp = false, HasSalt = false;
  if (eDescRiv.TypeVaryingTransport == "InterpolationFlux") {
    eTransport = InterpolateMeasurement(eDescRiv.ListPairTimeFlux, eTimeDay,
                                        maxAllowedTimeInterval);
    HasTransport = true;
  }
  if (eDescRiv.TypeVaryingTransport == "ConstantFlux") {
    eTransport = eDescRiv.ConstantFlux;
    HasTransport = true;
  }
  if (eDescRiv.TypeVaryingTransport == "MonthlyFlux") {
    eTransport = MonthlyInterpolation(eDescRiv.ListMonthlyFlux, eTimeDay);
    HasTransport = true;
  }
  if (eDescRiv.TypeVaryingTransport == "RegularBurstOutflow") {
    double FrequencyDay = eDescRiv.FrequencyDay;
    double DurationHour = eDescRiv.DurationHour;
    // Flux is in cubic meter
    double TotalFlux = eDescRiv.TotalFlux;
    //
    double DurationSec = DurationHour * 3600;
    double DurationDay = DurationHour / 24;
    double InstantFlux = TotalFlux / DurationSec;
    //
    double eQuot = round(eTimeDay / FrequencyDay);
    double eDiff = eTimeDay - eQuot * FrequencyDay;
    if (fabs(eDiff) <= DurationDay / 2) {
      eTransport = InstantFlux;
    } else {
      eTransport = 0;
    }
    HasTransport = true;
  }
  //
  // Temperature
  //
  if (eDescRiv.TypeVaryingTemperature == "InterpolationTemp") {
    eTemp = InterpolateMeasurement(eDescRiv.ListPairTimeTemp, eTimeDay,
                                   maxAllowedTimeInterval);
    HasTemp = true;
  }
  if (eDescRiv.TypeVaryingTemperature == "MonthlyTemp") {
    eTemp = MonthlyInterpolation(eDescRiv.ListMonthlyTemp, eTimeDay);
    HasTemp = true;
  }
  if (eDescRiv.TypeVaryingTemperature == "ConstantTemp") {
    eTemp = eDescRiv.ConstantRiverTemperature;
    HasTemp = true;
  }
  //
  // Salt
  //
  if (eDescRiv.TypeVaryingTemperature == "InterpolationSalt") {
    eSalt = InterpolateMeasurement(eDescRiv.ListPairTimeSalt, eTimeDay,
                                   maxAllowedTimeInterval);
    HasSalt = true;
  }
  if (eDescRiv.TypeVaryingTemperature == "MonthlySalt") {
    eSalt = MonthlyInterpolation(eDescRiv.ListMonthlySalt, eTimeDay);
    HasSalt = true;
  }
  if (eDescRiv.TypeVaryingSalinity == "ConstantSalt") {
    eSalt = eDescRiv.ConstantRiverSalinity;
    HasSalt = true;
  }
  //
  if (!HasTransport) {
    std::cerr << "eDescRiv.name=" << eDescRiv.name << "\n";
    std::cerr << "We have HasTransport = " << HasTransport << "\n";
    std::cerr << "TypeVaryingTransport = " << eDescRiv.TypeVaryingTransport
              << "\n";
    std::cerr << "AllowedValues are : YearlyClimatology, RegularBurstOutflow\n";
    throw TerminalException{1};
  }
  if (!HasTemp) {
    std::cerr << "eDescRiv.name=" << eDescRiv.name << "\n";
    std::cerr << "eDescRiv.TypeVaryingTemperature = "
              << eDescRiv.TypeVaryingTemperature << "\n";
    std::cerr << "We have HasTemp = " << HasTemp << "\n";
    std::cerr << "Available parametrization: ConstantTemp, MonthlyTemp, InterpolationTemp\n";
    throw TerminalException{1};
  }
  if (!HasSalt) {
    std::cerr << "eDescRiv.name=" << eDescRiv.name << "\n";
    std::cerr << "eDescRiv.TypeVaryingSalinity = "
              << eDescRiv.TypeVaryingSalinity << "\n";
    std::cerr << "We have HasSalt = " << HasSalt << "\n";
    std::cerr << "Available parametrization: ConstantSalt, MonthlySalt, InterpolationSalt\n";
    throw TerminalException{1};
  }
  eTransport *= eDescRiv.ConstantFactorFlux;
  return {eTransport, eTemp, eSalt};
}

int RIVER_ExtendedMask(MyMatrix<uint8_t> const &MSK_rho, int const &iEta,
                       int const &iXi) {
  int eta_rho = MSK_rho.rows();
  int xi_rho = MSK_rho.cols();
  if (iEta >= 0 && iEta < eta_rho && iXi >= 0 && iXi < xi_rho)
    return MSK_rho(iEta, iXi);
  else
    return 0;
}

double RIVER_AngleReduction(double const &TheAng) {
  double pi = GetPI();
  double RetAng = TheAng;
  while (1) {
    //    std::cerr << "RetAng=" << RetAng << "\n";
    if (RetAng > pi) {
      RetAng -= 2 * pi;
    } else {
      if (RetAng < -pi)
        RetAng += 2 * pi;
      else
        break;
    }
  }
  return RetAng;
}

struct RecordAngleStatusRiver {
  Eigen::Tensor<int, 3> TotalListStatus;
  Eigen::Tensor<double, 3> TotalListAngle;
  MyMatrix<int> StatusRULD;
};

/*
  StatusRULD is a matrix of size (eta_rho, xi_rho)
  It is 1 if a point is wet but has one of Right, Up, Left or Down point land.

 */
RecordAngleStatusRiver
DetermineRiverPossibleCandidates(MyMatrix<uint8_t> const &MSK_rho,
                                 MyMatrix<double> const &ANG_rho) {
  double pi = 3.141592653589793;
  int eta_rho = ANG_rho.rows();
  int xi_rho = ANG_rho.cols();
  Eigen::Tensor<int, 3> TotalListStatus(4, eta_rho, xi_rho);
  Eigen::Tensor<double, 3> TotalListAngle(4, eta_rho, xi_rho);
  MyMatrix<int> StatusRULD(eta_rho, xi_rho);
  int sumMSK_rho = 0;
  for (int iEta = 0; iEta < eta_rho; iEta++)
    for (int iXi = 0; iXi < xi_rho; iXi++) {
      int iRight = 1, iLeft = 1, iUp = 1, iDown = 1;
      double dRight = 0, dLeft = 0, dUp = 0, dDown = 0;
      if (MSK_rho(iEta, iXi) == 1) {
        sumMSK_rho++;
        int iEtaRight = iEta;
        int iXiRight = iXi - 1;
        iRight = RIVER_ExtendedMask(MSK_rho, iEtaRight, iXiRight);
        dRight = ANG_rho(iEta, iXi);
        //
        int iEtaUp = iEta - 1;
        int iXiUp = iXi;
        iUp = RIVER_ExtendedMask(MSK_rho, iEtaUp, iXiUp);
        dUp = ANG_rho(iEta, iXi) + pi / 2;
        //
        int iEtaLeft = iEta;
        int iXiLeft = iXi + 1;
        iLeft = RIVER_ExtendedMask(MSK_rho, iEtaLeft, iXiLeft);
        dLeft = ANG_rho(iEta, iXi) + pi;
        //
        int iEtaDown = iEta + 1;
        int iXiDown = iXi;
        iDown = RIVER_ExtendedMask(MSK_rho, iEtaDown, iXiDown);
        dDown = ANG_rho(iEta, iXi) + 3 * (pi / 2);
      }
      TotalListStatus(0, iEta, iXi) = iRight;
      TotalListStatus(1, iEta, iXi) = iUp;
      TotalListStatus(2, iEta, iXi) = iLeft;
      TotalListStatus(3, iEta, iXi) = iDown;
      TotalListAngle(0, iEta, iXi) = dRight;
      TotalListAngle(1, iEta, iXi) = dUp;
      TotalListAngle(2, iEta, iXi) = dLeft;
      TotalListAngle(3, iEta, iXi) = dDown;
      StatusRULD(iEta, iXi) = 1 - iRight * iUp * iLeft * iDown;
    }
  int eProd = eta_rho * xi_rho;
  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho
            << " eProd=" << eProd << "\n";
  std::cerr << "sum(MSK_rho)=" << sumMSK_rho << "\n";
  std::cerr << "sum(StatusRULD)=" << StatusRULD.sum()
            << " min/max=" << StatusRULD.minCoeff() << " / "
            << StatusRULD.maxCoeff() << "\n";
  for (int i = 0; i < 4; i++) {
    MyMatrix<int> ArrStatus = DimensionExtraction(TotalListStatus, 0, i);
    MyMatrix<double> ArrAngle = DimensionExtraction(TotalListAngle, 0, i);
    std::cerr << "i=" << i << " sum(ArrStatus)=" << ArrStatus.sum()
              << " sum(ArrAngle)=" << ArrAngle.sum() << "\n";
  }
  return {TotalListStatus, TotalListAngle, StatusRULD};
}

MyVector<double> RetrieveListOfWeight(MyVector<double> const &Zr,
                                      MyVector<double> const &Zw,
                                      DescriptionRiver const &eDescRiv) {
  int N = Zr.size();
  bool DoPrint = false;
  if (DoPrint) {
    std::cerr << "N = " << N << "\n";
    for (int i = 0; i < N; i++) {
      std::cerr << "i=" << i << " Zr=" << Zr(i) << "\n";
    }
    int Np1 = Zw.size();
    std::cerr << "Np1 = " << Np1 << "\n";
    for (int i = 0; i < Np1; i++) {
      std::cerr << "i=" << i << " Zw=" << Zw(i) << "\n";
    }
  }
  MyVector<double> PreListWeight(N);
  if (eDescRiv.verticalShapeOption == "UniformUpperLayer") {
    for (int iS = 0; iS < N; iS++) {
      double eDepLevel = Zr(iS);
      double eHeight = Zw(iS + 1) - Zw(iS);
      double eWeight;
      if (eDepLevel > -eDescRiv.MaxDepth) {
        eWeight = eHeight;
      } else {
        eWeight = 0;
      }
      PreListWeight(iS) = eWeight;
    }
    double eSum = PreListWeight.sum();
    return PreListWeight / eSum;
  }
  if (eDescRiv.verticalShapeOption == "FixedDepth") {
    for (int iS = 0; iS < N; iS++)
      PreListWeight(iS) = 0;
    int iSfound = -1;
    double eDep = -eDescRiv.targetDepth;
    double thr = 0.00001;
    for (int iS = 0; iS < N; iS++) {
      if (Zw(iS + 1) > eDep - thr && thr + eDep > Zw(iS))
        iSfound = iS;
    }
    if (iSfound == -1) {
      std::cerr << "We did not find wanted depth\n";
      std::cerr << "Requested eDep=" << eDep << " N=" << N << " list of available depths:\n";
      for (int iS=0; iS<=N; iS++)
        std::cerr << "iS=" << iS << " Zw=" << Zw(iS) << "\n";
      throw TerminalException{1};
    }
    PreListWeight(iSfound) = 1;
    return PreListWeight;
  }
  std::cerr << "No option chosen for the vertical shape\n";
  throw TerminalException{1};
}

void CreateRiverFile(FullNamelist const &eFull) {
  SingleBlock eBlINPUT = eFull.get_block("INPUT");
  //
  // List of river names
  //
  std::vector<std::string> ListRiverName =
    eBlINPUT.get_list_string("ListRiverName");
  int nbRiver = ListRiverName.size();
  std::string RiverPrefix = eBlINPUT.get_string("RiverPrefix");
  std::string RiverSuffix = eBlINPUT.get_string("RiverSuffix");
  std::string RiverFile = eBlINPUT.get_string("RiverFile");
  int N = eBlINPUT.get_int("ARVD_N");
  //
  // Timings
  //
  std::string RefTimeStr = eBlINPUT.get_string("RefTime");
  double RefTime = CT2MJD(RefTimeStr);
  std::string strBeginTime = eBlINPUT.get_string("BEGTC");
  double BeginTime = CT2MJD(strBeginTime);
  std::string strEndTime = eBlINPUT.get_string("ENDTC");
  double EndTime = CT2MJD(strEndTime);
  std::string UNITC = eBlINPUT.get_string("UNITC");
  double DELTC = eBlINPUT.get_double("DELTC");
  double DeltaTime = GetIntervalSize(DELTC, UNITC);
  double maxAllowedTimeInterval =
    eBlINPUT.get_double("maxAllowedTimeInterval");
  //
  // Now reading the vertical stratification
  //
  int Vtransform = eBlINPUT.get_int("ARVD_Vtransform");
  int Vstretching = eBlINPUT.get_int("ARVD_Vstretching");
  double Tcline = eBlINPUT.get_double("ARVD_Tcline");
  double hc = eBlINPUT.get_double("ARVD_hc");
  double theta_s = eBlINPUT.get_double("ARVD_theta_s");
  double theta_b = eBlINPUT.get_double("ARVD_theta_b");
  ARVDtyp ARVD = ROMSgetARrayVerticalDescription(N, Vtransform, Vstretching,
                                                 Tcline, hc, theta_s, theta_b);
  //
  // Reading Tracers different from the temp/salt
  // We have to treat the temp/salt separately since they use the
  // SetRiverSalinity / SetRiverTemp
  //
  std::string eFileExternal = eBlINPUT.get_string("ExternalInfoFile");
  std::vector<std::string> ListAdditionalTracers =
    eBlINPUT.get_list_string("ListAdditionalTracers");
  std::vector<TracerDescription> ListTracerStringDescription;
  for (auto eTracerName : ListAdditionalTracers) {
    TracerDescription eTracer =
        RetrieveTracerDescription(eFileExternal, eTracerName);
    ListTracerStringDescription.push_back(eTracer);
  }
  int nbAdditionalTracer = ListTracerStringDescription.size();
  //
  // Now reading the input file for the rivers
  //
  std::vector<DescriptionRiver> ListRiverDescription(nbRiver);
  std::cerr << "nbRiver=" << nbRiver << "\n";
  for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
    std::string eRiverName = ListRiverName[iRiver];
    std::string eRiverNameFull = RiverPrefix + eRiverName + RiverSuffix;
    std::cerr << "iRiver=" << iRiver << "/" << nbRiver
              << "   RiverDescriptionFile=" << eRiverNameFull << "\n";
    ListRiverDescription[iRiver] = ReadRiverDescription(eRiverNameFull);
  }
  //
  // Now reading the grid arrays and related stuff
  //
  std::string GridFile = eBlINPUT.get_string("GridFile");
  GridArray GrdArr = NC_ReadRomsGridFile(GridFile);
  std::cerr << "We have GrdArr\n";
  int eta_rho = GrdArr.GrdArrRho.LON.rows();
  int xi_rho = GrdArr.GrdArrRho.LON.cols();
  struct ijLL {
    int i;
    int j;
    double lon;
    double lat;
  };
  RecordAngleStatusRiver RecAngStatRiv = DetermineRiverPossibleCandidates(
      GrdArr.GrdArrRho.MSK, GrdArr.GrdArrRho.ANG);
  std::cerr << "We have RecAngStatRiv\n";
  std::vector<ijLL> ListIJLLruld;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      if (RecAngStatRiv.StatusRULD(i, j) == 1) {
        double lon = GrdArr.GrdArrRho.LON(i, j);
        double lat = GrdArr.GrdArrRho.LAT(i, j);
        ListIJLLruld.push_back({i, j, lon, lat});
      }
  std::vector<ijLL> ListIJLLwet;
  for (int i = 0; i < eta_rho; i++)
    for (int j = 0; j < xi_rho; j++)
      if (GrdArr.GrdArrRho.MSK(i, j) == 1) {
        double lon = GrdArr.GrdArrRho.LON(i, j);
        double lat = GrdArr.GrdArrRho.LAT(i, j);
        ListIJLLwet.push_back({i, j, lon, lat});
      }
  std::cerr << "sum(RecAngStatRiv.StatusRULD)="
            << RecAngStatRiv.StatusRULD.sum() << "\n";
  auto GetNearest = [&](double const &lon, double const &lat,
                        std::vector<ijLL> const &ListIJLL) -> ijLL {
    bool IsFirst = true;
    ijLL TheSel;
    double MinNorm = -1;
    for (auto &eEnt : ListIJLL) {
      double dist = GeodesicDistanceKM(eEnt.lon, eEnt.lat, lon, lat);
      if (IsFirst) {
        TheSel = eEnt;
        IsFirst = false;
        MinNorm = dist;
      } else {
        if (dist < MinNorm) {
          MinNorm = dist;
          TheSel = eEnt;
        }
      }
    }
    return TheSel;
  };
  MyMatrix<double> InflMatrixSph_rho = GRID_GetMaxRadiusInfluence_kernel(
      GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT, "spherical");
  //
  // Now creating the list of signs, ETA, XI and direction
  //
  std::vector<int> ListSign, ListETA, ListXI, ListDirection, ListIRiver;
  std::vector<int> ListETAland, ListXIland, ListETAsea, ListXIsea;
  std::vector<double> ListDepArrival;
  std::vector<std::string> ListStringName;
  double pi = GetPI();
  auto GetIJDSriverCase = [&](double const &lon, double const &lat,
                              double const &direction) -> ijdsInfo {
    double TheAngOrig = direction * (pi / 180);
    ijLL ijllNear = GetNearest(lon, lat, ListIJLLruld);
    double distance = GeodesicDistanceKM(ijllNear.lon, ijllNear.lat, lon, lat);
    int eEtaSea = ijllNear.i;
    int eXiSea = ijllNear.j;
    double dist1 = InflMatrixSph_rho(eEtaSea, eXiSea);
    if (distance > dist1) {
      std::cerr << "GetIJDSriverCase distance=" << distance
                << " dist1=" << dist1 << " INVALID distance.\n";
      return {0, -1, -1, -1, -1};
    }
    int iChoice = -1;
    double minDeltaAng = 4545;
    for (int i = 0; i < 4; i++) {
      if (RecAngStatRiv.TotalListStatus(i, eEtaSea, eXiSea) == 0) {
        double TheAngNew = RecAngStatRiv.TotalListAngle(i, eEtaSea, eXiSea);
        double deltaAng = TheAngNew - TheAngOrig;
        double deltaAngRed = abs(RIVER_AngleReduction(deltaAng));
        if (deltaAngRed < minDeltaAng) {
          iChoice = i;
          minDeltaAng = deltaAngRed;
        }
      }
    }
    return RetrieveIJDSarray(eEtaSea, eXiSea, iChoice);
  };
  auto GetIJDSwetCase = [&](double const &lon, double const &lat,
                            int const &ChoiceSelect) -> ijdsInfo {
    ijLL ijllNear = GetNearest(lon, lat, ListIJLLwet);
    double MinDist = GeodesicDistanceKM(ijllNear.lon, ijllNear.lat, lon, lat);
    int eEtaSea = ijllNear.i;
    int eXiSea = ijllNear.j;
    double distInfl = InflMatrixSph_rho(eEtaSea, eXiSea);
    if (MinDist > distInfl) {
      std::cerr << "GetIJDSwetCase : distInfl=" << distInfl
                << " MinDist=" << MinDist << " INVALID distance\n";
      return {0, -1, -1, -1, -1};
    }
    return RetrieveIJDSarray(eEtaSea, eXiSea, ChoiceSelect);
  };
  const MyMatrix<double> &DEP = GetDEP(GrdArr.GrdArrRho);
  for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
    DescriptionRiver eDescRiv = ListRiverDescription[iRiver];
    //
    ijdsInfo RecordInfo;
    bool IsDone = false;
    if (eDescRiv.WScase == "River") {
      RecordInfo =
          GetIJDSriverCase(eDescRiv.lon, eDescRiv.lat, eDescRiv.direction);
      IsDone = true;
    }
    if (eDescRiv.WScase == "Wet") {
      RecordInfo =
          GetIJDSwetCase(eDescRiv.lon, eDescRiv.lat, eDescRiv.ChoiceSelect);
      IsDone = true;
    }
    if (eDescRiv.WScase == "Direct") {
      RecordInfo = {1, eDescRiv.iSelect, eDescRiv.jSelect, eDescRiv.DirSelect,
                    eDescRiv.SignSelect};
      IsDone = true;
    }
    if (!IsDone) {
      std::cerr << "We failed to find matching algorithm\n";
      throw TerminalException{1};
    }
    std::cerr << "iRiver=" << iRiver << " name=" << eDescRiv.name
              << " status=" << RecordInfo.eStatus << "\n";
    if (RecordInfo.eStatus == 0) {
      std::cerr << "  No valid choice found WScase = " << eDescRiv.WScase
                << "\n";
    } else {
      if (RecordInfo.DirSelect != 0 && RecordInfo.DirSelect != 1) {
        std::cerr << "DirSelect=" << RecordInfo.DirSelect << "\n";
        std::cerr << "DirSelect should be 0 or 1\n";
        throw TerminalException{1};
      }
      if (RecordInfo.SignSelect != -1 && RecordInfo.SignSelect != 1) {
        std::cerr << "SignSelect=" << RecordInfo.SignSelect << "\n";
        std::cerr << "SignSelect should be -1 or 1\n";
        throw TerminalException{1};
      }
      ijSeaLand recIJSL = GetArrayIjSeaLand(RecordInfo);
      std::cerr << "recIJSL Sea(i/j)=" << recIJSL.iSea << " / " << recIJSL.jSea << " Land(i/j)=" << recIJSL.iLand << " / " << recIJSL.jLand << "\n";
      ListETAsea.push_back(recIJSL.iSea);
      ListXIsea.push_back(recIJSL.jSea);
      ListETAland.push_back(recIJSL.iLand);
      ListXIland.push_back(recIJSL.jLand);
      double eDEP = DEP(recIJSL.iSea, recIJSL.jSea);
      double eLONsea = GrdArr.GrdArrRho.LON(recIJSL.iSea, recIJSL.jSea);
      double eLATsea = GrdArr.GrdArrRho.LAT(recIJSL.iSea, recIJSL.jSea);
      std::string strLONland = "unset", strLATland = "unset";
      if (recIJSL.iLand >= 0 && recIJSL.jLand >= 0 && recIJSL.iLand < eta_rho && recIJSL.jLand < xi_rho) {
        strLONland = std::to_string(GrdArr.GrdArrRho.LON(recIJSL.iLand, recIJSL.jLand));
        strLATland = std::to_string(GrdArr.GrdArrRho.LAT(recIJSL.iLand, recIJSL.jLand));
      }
      ListDepArrival.push_back(eDEP);
      //
      ListETA.push_back(RecordInfo.iSelect);
      ListXI.push_back(RecordInfo.jSelect);
      ListDirection.push_back(RecordInfo.DirSelect);
      ListSign.push_back(RecordInfo.SignSelect);
      ListIRiver.push_back(iRiver);
      ListStringName.push_back(eDescRiv.name);
      std::cerr << "    Found river iS=" << RecordInfo.iSelect
                << " jS=" << RecordInfo.jSelect
                << " iD=" << RecordInfo.DirSelect
                << " iN=" << RecordInfo.SignSelect << "\n";
      std::cerr << "    Land(i,j)=" << recIJSL.iLand << "/" << recIJSL.jLand
                << " " << strLONland << "/" << strLATland
                << "   Sea(i,j)=" << recIJSL.iSea << "/" << recIJSL.jSea << " "
                << eLONsea << "/" << eLATsea << " dep=" << eDEP << "\n";
    }
  }
  int nbRiverReal = ListIRiver.size();
  std::cerr << "nbRiver=" << nbRiver << " nbRiverReal=" << nbRiverReal << "\n";
  MyMatrix<int> MatrixIdx(N, nbRiverReal);
  int idx = 0;
  for (int iN = 0; iN < N; iN++)
    for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++) {
      MatrixIdx(iN, iRiverReal) = idx;
      idx++;
    }
  //
  // Basic definition of the river file
  //
  if (!FILE_IsFileMakeable(RiverFile)) {
    std::cerr << "Request to create file RiverFile=" << RiverFile << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  // Now dimensions
  netCDF::NcDim eDimRiver = dataFile.addDim("river", nbRiverReal);
  netCDF::NcDim eDimSvert = dataFile.addDim("s_rho", N);
  netCDF::NcDim eDimRiverTime = dataFile.addDim("river_time");
  netCDF::NcDim eDim19 = dataFile.addDim("dateString", 19);
  // Now variables
  std::vector<std::string> LDim1{"river"};
  std::vector<std::string> LDim2{"s_rho", "river"};
  std::vector<std::string> LDim3{"river_time", "s_rho", "river"};
  std::vector<std::string> LDim4{"river_time"};
  std::vector<std::string> LDim5{"river_time", "dateString"};
  std::vector<std::string> LDim6{"river_time", "river"};
  netCDF::NcVar varRiver = dataFile.addVar("river", "double", LDim1);
  varRiver.putAtt("long_name", "river runoff identification number");
  varRiver.putAtt("units", "nondimensional");
  varRiver.putAtt("field", "river, scalar");
  netCDF::NcVar varEposition =
      dataFile.addVar("river_Eposition", "double", LDim1);
  varEposition.putAtt("long_name", "river ETA-position at RHO-points");
  varEposition.putAtt("units", "nondimensional");
  varEposition.putAtt("field", "river_Eposition, scalar");
  netCDF::NcVar varXposition =
      dataFile.addVar("river_Xposition", "double", LDim1);
  varXposition.putAtt("long_name", "river XI-position at RHO-points");
  varXposition.putAtt("units", "nondimensional");
  varXposition.putAtt("field", "river_Xposition, scalar");
  netCDF::NcVar varVshape = dataFile.addVar("river_Vshape", "double", LDim2);
  varVshape.putAtt("long_name", "river runoff mass transport vertical profile");
  varVshape.putAtt("units", "nondimensional");
  varVshape.putAtt("field", "river_Vshape, scalar");
  netCDF::NcVar varDirection =
      dataFile.addVar("river_direction", "double", LDim1);
  varDirection.putAtt("long_name", "river runoff direction");
  varDirection.putAtt("units", "nondimensional");
  varDirection.putAtt("field", "river_direction, scalar");
  netCDF::NcVar varFlag = dataFile.addVar("river_flag", "double", LDim1);
  varFlag.putAtt("long_name", "river runoff tracer flag");
  varFlag.putAtt("units", "nondimensional");
  varFlag.putAtt("option_0", "all tracers are off");
  varFlag.putAtt("option_1", "only temperature is on");
  varFlag.putAtt("option_2", "only salinity is on");
  varFlag.putAtt("option_3", "only both are on");
  varFlag.putAtt("field", "river_flag, scalar");
  //
  // The tracer variables. First salt/temp and then other tracers.
  //
  netCDF::NcVar varSalt = dataFile.addVar("river_salt", "double", LDim3);
  varSalt.putAtt("long_name", "river runoff salinity");
  varSalt.putAtt("units", "PSU");
  varSalt.putAtt("field", "river_salt, scalar, series");
  netCDF::NcVar varTemp = dataFile.addVar("river_temp", "double", LDim3);
  varTemp.putAtt("long_name", "river runoff potential temperature");
  varTemp.putAtt("units", "Celsius");
  varTemp.putAtt("field", "river_temp, scalar, series");
  std::vector<netCDF::NcVar> ListVarTracer;
  for (auto eTracerDesc : ListTracerStringDescription) {
    netCDF::NcVar varTracer =
        dataFile.addVar(eTracerDesc.namefull, "double", LDim3);
    varTracer.putAtt("long_name", eTracerDesc.long_name);
    varTracer.putAtt("units", eTracerDesc.units);
    varTracer.putAtt("field", eTracerDesc.field);
    ListVarTracer.push_back(varTracer);
  }
  netCDF::NcVar varRiverTime = dataFile.addVar("river_time", "double", LDim4);
  varRiverTime.putAtt("long_name", "river runoff time");
  std::string dateStr = DATE_ConvertMjd2mystringPres(RefTime);
  std::string attTime = "days since " + dateStr;
  varRiverTime.putAtt("units", attTime);
  varRiverTime.putAtt("calendar", "gregorian");
  varRiverTime.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTimeStr =
      dataFile.addVar("river_time_str", "char", LDim5);
  netCDF::NcVar varTransport =
      dataFile.addVar("river_transport", "double", LDim6);
  varTransport.putAtt("long_name",
                      "river runoff vertically integrated mass transport");
  varTransport.putAtt("units", "meter3 second-1");
  varTransport.putAtt("field", "river_transport, scalar, series");
  //
  // List of all names
  //
  std::string stringListNameRiver;
  for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++) {
    if (iRiverReal > 0)
      stringListNameRiver += ",";
    stringListNameRiver += ListStringName[iRiverReal];
  }
  dataFile.putAtt("rivers", stringListNameRiver);
  //
  // Function for write downs
  //
  auto WriteDownNbRiver = [&](netCDF::NcVar &eVAR,
                              std::vector<int> const &ListVal) -> void {
    std::vector<double> A(nbRiverReal);
    for (int i = 0; i < nbRiverReal; i++)
      A[i] = static_cast<double>(ListVal[i]);
    eVAR.putVar(A.data());
  };
  //
  // Now easy definitions
  //
  std::vector<int> ListIdxRiver(nbRiverReal);
  for (int i = 0; i < nbRiverReal; i++)
    ListIdxRiver[i] = i + 1;
  WriteDownNbRiver(varRiver, ListIdxRiver);
  WriteDownNbRiver(varEposition, ListETA);
  WriteDownNbRiver(varXposition, ListXI);
  WriteDownNbRiver(varDirection, ListDirection);
  //
  // Now creating whether we assign salinity and/or temperature
  //
  std::vector<int> ListFlag(nbRiver);
  for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++) {
    int iRiver = ListIRiver[iRiverReal];
    int eFlag = 0;
    bool SetTemp = ListRiverDescription[iRiver].SetRiverTemperature;
    bool SetSalt = ListRiverDescription[iRiver].SetRiverSalinity;
    if (SetTemp && SetSalt)
      eFlag = 3;
    if (!SetTemp && SetSalt)
      eFlag = 2;
    if (SetTemp && !SetSalt)
      eFlag = 1;
    if (!SetTemp && !SetSalt)
      eFlag = 0;
    std::cerr << "iRiverReal=" << iRiverReal << " iRiver=" << iRiver
              << " name=" << ListRiverDescription[iRiver].name
              << " SetTemp=" << SetTemp << " SetSalt=" << SetSalt
              << " eFlag=" << eFlag << "\n";
    ListFlag[iRiverReal] = eFlag;
  }
  WriteDownNbRiver(varFlag, ListFlag);
  //
  // Now we write down the vertical shape
  //
  std::vector<double> Ashape(N * nbRiverReal);
  for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++) {
    int iRiver = ListIRiver[iRiverReal];
    DescriptionRiver eDescRiv = ListRiverDescription[iRiver];
    double eDep = ListDepArrival[iRiverReal];
    double eZeta = 0;
    MyVector<double> Zr_out = GetVertCoord_R(ARVD, eDep, eZeta);
    MyVector<double> Zw_out = GetVertCoord_W(ARVD, eDep, eZeta);
    MyVector<double> ListWeight =
        RetrieveListOfWeight(Zr_out, Zw_out, eDescRiv);
    for (int iS = 0; iS < N; iS++) {
      int idx = MatrixIdx(iS, iRiverReal);
      Ashape[idx] = ListWeight(iS);
    }
  }
  varVshape.putVar(Ashape.data());
  //
  // The time loop of creating the data
  //
  double CurrentTime = BeginTime;
  double epsilon = 0.00001;
  int pos = 0;
  std::vector<double> Atransport(nbRiverReal);
  std::cerr << "Now writing flux and tracers\n";
  while (1) {
    std::string strPres = DATE_ConvertMjd2mystringPres(CurrentTime);
    std::cerr << "pos=" << pos << " time=" << strPres << "\n";
    // First writing the time
    std::vector<size_t> start1{size_t(pos)};
    std::vector<size_t> count1{1};
    double DiffTime = CurrentTime - RefTime;
    varRiverTime.putVar(start1, count1, &DiffTime);
    std::vector<size_t> start2{size_t(pos), 0};
    std::vector<size_t> count2{1, 19};
    varRiverTimeStr.putVar(start2, count2, strPres.c_str());
    //
    std::vector<std::vector<double>> ListListTracerValue(
        nbAdditionalTracer, std::vector<double>(nbRiverReal));
    std::vector<double> ListTemp(nbRiverReal), ListSalt(nbRiverReal);
    for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++) {
      int iRiver = ListIRiver[iRiverReal];
      TransTempSalt eTTS = RetrieveTTS(ListRiverDescription[iRiver],
                                       CurrentTime, maxAllowedTimeInterval);
      ListTemp[iRiverReal] = eTTS.eTemp;
      ListSalt[iRiverReal] = eTTS.eSalt;
      std::vector<double> ListTracerVal(nbAdditionalTracer);
      for (int iAdditionalTracer = 0; iAdditionalTracer < nbAdditionalTracer;
           iAdditionalTracer++) {
        std::string TracerName =
            ListTracerStringDescription[iAdditionalTracer].name;
        try {
          TracerTimeVariability const &ttv =
              ListRiverDescription[iRiver].MapTracerDesc.at(TracerName);
          double eValue =
              RetrieveTracerValue(ttv, CurrentTime, maxAllowedTimeInterval);
          ListTracerVal[iAdditionalTracer] = eValue;
        } catch (...) {
          std::cerr << "Failed to find a matching entry for this tracer and "
                       "iRiverReal="
                    << iRiverReal << "\n";
          std::cerr << "For TracerName=" << TracerName << "\n";
          throw TerminalException{1};
        }
      }
      for (int iTracer = 0; iTracer < nbAdditionalTracer; iTracer++) {
        std::cerr << "  iRiverReal=" << iRiverReal << " iTracer=" << iTracer
                  << "   TracerValue=" << ListTracerVal[iTracer] << "\n";
        ListListTracerValue[iTracer][iRiverReal] = ListTracerVal[iTracer];
      }
      Atransport[iRiverReal] = ListSign[iRiverReal] * eTTS.eTransport;
    }
    // Now writing the temp/salt
    std::vector<size_t> startTracer{size_t(pos), 0, 0};
    std::vector<size_t> countTracer{1, size_t(N), size_t(nbRiverReal)};
    auto PutTracer = [&](netCDF::NcVar &varTracer,
                         std::vector<double> const &ListValue) -> void {
      std::vector<double> A(N * nbRiverReal);
      for (int iRiverReal = 0; iRiverReal < nbRiverReal; iRiverReal++)
        for (int iS = 0; iS < N; iS++) {
          int idx = MatrixIdx(iS, iRiverReal);
          A[idx] = ListValue[iRiverReal];
        }
      varTracer.putVar(startTracer, countTracer, A.data());
    };
    PutTracer(varSalt, ListSalt);
    PutTracer(varTemp, ListTemp);
    // Now writing the other tracers
    for (int iTracer = 0; iTracer < nbAdditionalTracer; iTracer++)
      PutTracer(ListVarTracer[iTracer], ListListTracerValue[iTracer]);
    // Now writing the transport
    std::vector<size_t> startTrans{size_t(pos), 0};
    std::vector<size_t> countTrans{1, size_t(nbRiverReal)};
    varTransport.putVar(startTrans, countTrans, Atransport.data());
    // Now maybe leaving
    CurrentTime += DeltaTime;
    pos++;
    if (CurrentTime > EndTime + epsilon)
      break;
  }
  std::cerr << "Finished creation of   RiverFile = " << RiverFile << "\n";
}

void MergeRiverFile(std::string const &RiverFile,
                    std::vector<std::string> const &ListRiverFile) {
  int nbFile = ListRiverFile.size();
  if (nbFile == 0) {
    std::cerr << "We have |ListRiverFile| = 0\n";
    throw TerminalException{1};
  }
  std::vector<int> ListNbRiver(nbFile);
  std::vector<int> ListNbTime(nbFile);
  std::vector<int> ListSrho(nbFile);
  for (int iFile = 0; iFile < nbFile; iFile++) {
    std::string eFile = ListRiverFile[iFile];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    std::string strRiver = "river";
    std::string strRiverTime = "river_time";
    std::string strSrho = "s_rho";
    ListNbRiver[iFile] = NC_ReadDimension(dataFile, strRiver);
    ListNbTime[iFile] = NC_ReadDimension(dataFile, strRiverTime);
    ListSrho[iFile] = NC_ReadDimension(dataFile, strSrho);
  }
  auto CheckIdentityDimension = [&](std::vector<int> const &ListVal) -> void {
    int eMin = VectorMin(ListVal);
    int eMax = VectorMax(ListVal);
    if (eMin != eMax) {
      std::cerr << "Wrong identity of the dimensions\n";
      std::cerr << "eMin=" << eMin << " eMax=" << eMax << "\n";
      throw TerminalException{1};
    }
  };
  CheckIdentityDimension(ListNbTime);
  CheckIdentityDimension(ListSrho);
  int nbRiver = VectorSum(ListNbRiver);
  int nbTime = ListNbTime[0];
  int s_rho = ListSrho[0];
  //
  // Reading time from Single river file
  //
  MyVector<double> ListTime;
  auto GetRefTime = [&]() -> double {
    std::string SingleRiverFile = ListRiverFile[0];
    netCDF::NcFile dataFile(SingleRiverFile, netCDF::NcFile::read);
    netCDF::NcVar data = dataFile.getVar("river_time");
    ListTime = NC_ReadVariable_data(data);
    return CT2MJD("19680523.000000");
  };
  double RefTime = GetRefTime();
  //
  // Basic definition of the river file
  //
  MyMatrix<int> MatrixIdx(s_rho, nbRiver);
  int idxM = 0;
  for (int iS = 0; iS < s_rho; iS++) {
    for (int iRiver = 0; iRiver < nbRiver; iRiver++) {
      MatrixIdx(iS, iRiver) = idxM;
      idxM++;
    }
  }
  //
  // Basic definition of the river file
  //
  if (!FILE_IsFileMakeable(RiverFile)) {
    std::cerr << "Request to create file RiverFile=" << RiverFile << "\n";
    std::cerr << "but the directory does not exist\n";
    throw TerminalException{1};
  }
  netCDF::NcFile dataFile(RiverFile, netCDF::NcFile::replace,
                          netCDF::NcFile::nc4);
  // Now dimensions
  netCDF::NcDim eDimRiver = dataFile.addDim("river", nbRiver);
  netCDF::NcDim eDimSvert = dataFile.addDim("s_rho", s_rho);
  netCDF::NcDim eDimRiverTime = dataFile.addDim("river_time");
  netCDF::NcDim eDim19 = dataFile.addDim("dateString", 19);
  // Now variables
  std::vector<std::string> LDim1{"river"};
  std::vector<std::string> LDim2{"s_rho", "river"};
  std::vector<std::string> LDim3{"river_time", "s_rho", "river"};
  std::vector<std::string> LDim4{"river_time"};
  std::vector<std::string> LDim5{"river_time", "dateString"};
  std::vector<std::string> LDim6{"river_time", "river"};
  netCDF::NcVar varRiver = dataFile.addVar("river", "double", LDim1);
  varRiver.putAtt("long_name", "river runoff identification number");
  varRiver.putAtt("units", "nondimensional");
  varRiver.putAtt("field", "river, scalar");
  netCDF::NcVar varEposition =
      dataFile.addVar("river_Eposition", "double", LDim1);
  varEposition.putAtt("long_name", "river ETA-position at RHO-points");
  varEposition.putAtt("units", "nondimensional");
  varEposition.putAtt("field", "river_Eposition, scalar");
  netCDF::NcVar varXposition =
      dataFile.addVar("river_Xposition", "double", LDim1);
  varXposition.putAtt("long_name", "river XI-position at RHO-points");
  varXposition.putAtt("units", "nondimensional");
  varXposition.putAtt("field", "river_Xposition, scalar");
  netCDF::NcVar varVshape = dataFile.addVar("river_Vshape", "double", LDim2);
  varVshape.putAtt("long_name", "river runoff mass transport vertical profile");
  varVshape.putAtt("units", "nondimensional");
  varVshape.putAtt("field", "river_Vshape, scalar");
  netCDF::NcVar varDirection =
      dataFile.addVar("river_direction", "double", LDim1);
  varDirection.putAtt("long_name", "river runoff direction");
  varDirection.putAtt("units", "nondimensional");
  varDirection.putAtt("field", "river_direction, scalar");
  netCDF::NcVar varFlag = dataFile.addVar("river_flag", "double", LDim1);
  varFlag.putAtt("long_name", "river runoff tracer flag");
  varFlag.putAtt("units", "nondimensional");
  varFlag.putAtt("option_0", "all tracers are off");
  varFlag.putAtt("option_1", "only temperature is on");
  varFlag.putAtt("option_2", "only salinity is on");
  varFlag.putAtt("option_3", "only both are on");
  varFlag.putAtt("field", "river_flag, scalar");
  netCDF::NcVar varSalt = dataFile.addVar("river_salt", "double", LDim3);
  varSalt.putAtt("long_name", "river runoff salinity");
  varSalt.putAtt("units", "PSU");
  varSalt.putAtt("field", "river_salt, scalar, series");
  netCDF::NcVar varTemp = dataFile.addVar("river_temp", "double", LDim3);
  varTemp.putAtt("long_name", "river runoff potential temperature");
  varTemp.putAtt("units", "Celsius");
  varTemp.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTime = dataFile.addVar("river_time", "double", LDim4);
  varRiverTime.putAtt("long_name", "river runoff time");
  std::string dateStr = DATE_ConvertMjd2mystringPres(RefTime);
  std::string attTime = "days since " + dateStr;
  varRiverTime.putAtt("units", attTime);
  varRiverTime.putAtt("calendar", "gregorian");
  varRiverTime.putAtt("field", "river_temp, scalar, series");
  netCDF::NcVar varRiverTimeStr =
      dataFile.addVar("river_time_str", "char", LDim5);
  netCDF::NcVar varTransport =
      dataFile.addVar("river_transport", "double", LDim6);
  varTransport.putAtt("long_name",
                      "river runoff vertically integrated mass transport");
  varTransport.putAtt("units", "meter3 second-1");
  varTransport.putAtt("field", "river_transport, scalar, series");
  //
  // List of all names
  //
  std::string stringListNameRiver;
  for (int iFile = 0; iFile < nbFile; iFile++) {
    std::string SingleRiverFile = ListRiverFile[0];
    netCDF::NcFile dataFile(SingleRiverFile, netCDF::NcFile::read);
    netCDF::NcGroupAtt RiversAtt = dataFile.getAtt("rivers");
    std::string strRiverName;
    RiversAtt.getValues(strRiverName);
    if (iFile > 0)
      stringListNameRiver += ",";
    stringListNameRiver += strRiverName;
  }
  dataFile.putAtt("rivers", stringListNameRiver);
  //
  // Function for write downs
  //
  auto ExtractGlobalField =
      [&](std::string const &VarName) -> std::vector<double> {
    std::vector<double> retVal;
    for (int iFile = 0; iFile < nbFile; iFile++) {
      std::string eFile = ListRiverFile[iFile];
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar data = dataFile.getVar(VarName);
      MyVector<double> ListVal = NC_ReadVariable_data(data);
      for (int iRiver = 0; iRiver < ListNbRiver[iFile]; iRiver++) {
        retVal.push_back(ListVal(iRiver));
      }
    }
    return retVal;
  };
  auto WriteDownEntry = [&](netCDF::NcVar &eVAR,
                            std::vector<double> ListVal) -> void {
    std::vector<double> A(nbRiver);
    for (int iRiver = 0; iRiver < nbRiver; iRiver++)
      A[iRiver] = ListVal[iRiver];
    eVAR.putVar(A.data());
  };
  auto MergeAndWriteDownEntry = [&](netCDF::NcVar &eVAR,
                                    std::string const &VarName) -> void {
    std::vector<double> ListVal = ExtractGlobalField(VarName);
    WriteDownEntry(eVAR, ListVal);
  };
  //
  // Now easy definitions
  //
  std::vector<double> ListIdxRiver(nbRiver);
  for (int i = 0; i < nbRiver; i++)
    ListIdxRiver[i] = i + 1;
  WriteDownEntry(varRiver, ListIdxRiver);
  MergeAndWriteDownEntry(varEposition, "river_Eposition");
  MergeAndWriteDownEntry(varXposition, "river_Xposition");
  MergeAndWriteDownEntry(varDirection, "river_direction");
  MergeAndWriteDownEntry(varFlag, "river_flag");
  //
  // Now we write down the vertical shape
  //
  std::vector<double> Ashape(s_rho * nbRiver);
  int idx2 = 0;
  for (int iFile = 0; iFile < nbFile; iFile++) {
    std::string eFile = ListRiverFile[iFile];
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar data = dataFile.getVar("river_Vshape");
    MyMatrix<double> eMat = NC_Read2Dvariable_data(data);
    int nbRiverLoc = ListNbRiver[iFile];
    for (int iRiverLoc = 0; iRiverLoc < nbRiverLoc; iRiverLoc++) {
      for (int iS = 0; iS < s_rho; iS++) {
        int idxB = MatrixIdx(iS, idx2);
        Ashape[idxB] = eMat(iS, iRiverLoc);
      }
      idx2++;
    }
  }
  varVshape.putVar(Ashape.data());
  //
  // The time loop of creating the data
  //
  std::vector<double> Asalt(s_rho * nbRiver);
  std::vector<double> Atemp(s_rho * nbRiver);
  std::vector<double> Atransport(nbRiver);
  for (int iTime = 0; iTime < nbTime; iTime++) {
    // First writing the time
    double CurrentTime = ListTime(iTime);
    std::string strPres = DATE_ConvertMjd2mystringPres(RefTime + CurrentTime);
    std::vector<size_t> start1{size_t(iTime)};
    std::vector<size_t> count1{1};
    varRiverTime.putVar(start1, count1, &CurrentTime);
    std::vector<size_t> start2{size_t(iTime), 0};
    std::vector<size_t> count2{1, 19};
    varRiverTimeStr.putVar(start2, count2, strPres.c_str());
    //
    int idx3 = 0;
    for (int iFile = 0; iFile < nbFile; iFile++) {
      int nbRiverLoc = ListNbRiver[iFile];
      std::vector<size_t> startTS{size_t(iTime), 0, 0};
      std::vector<size_t> countTS{1, size_t(s_rho), size_t(nbRiverLoc)};
      std::vector<size_t> startTrans{size_t(iTime), 0};
      std::vector<size_t> countTrans{1, size_t(nbRiverLoc)};
      //
      std::string eFile = ListRiverFile[iFile];
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar dataSalt = dataFile.getVar("river_salt");
      netCDF::NcVar dataTemp = dataFile.getVar("river_temp");
      netCDF::NcVar dataTrans = dataFile.getVar("river_transport");
      MyVector<double> ArrSalt =
          NC_ReadVariable_data_start_count(dataSalt, startTS, countTS);
      MyVector<double> ArrTemp =
          NC_ReadVariable_data_start_count(dataTemp, startTS, countTS);
      MyVector<double> ArrTrans =
          NC_ReadVariable_data_start_count(dataTrans, startTrans, countTrans);
      MyMatrix<int> MatrixIdxLoc(s_rho, nbRiverLoc);
      int idxM2 = 0;
      for (int iS = 0; iS < s_rho; iS++)
        for (int iRiver = 0; iRiver < nbRiverLoc; iRiver++) {
          MatrixIdxLoc(iS, iRiver) = idxM2;
          idxM2++;
        }
      for (int iRiverLoc = 0; iRiverLoc < nbRiverLoc; iRiverLoc++) {
        for (int iS = 0; iS < s_rho; iS++) {
          int idx4 = MatrixIdxLoc(iS, iRiverLoc);
          Asalt[s_rho * idx3 + iS] = ArrSalt(idx4);
          Atemp[s_rho * idx3 + iS] = ArrTemp(idx4);
        }
        Atransport[idx3] = ArrTrans(iRiverLoc);
      }
    }
    // Now writing your data
    std::vector<size_t> startTS{size_t(iTime), 0, 0};
    std::vector<size_t> countTS{1, size_t(s_rho), size_t(nbRiver)};
    varSalt.putVar(startTS, countTS, Asalt.data());
    varTemp.putVar(startTS, countTS, Atemp.data());
    std::vector<size_t> startTrans{size_t(iTime), 0};
    std::vector<size_t> countTrans{1, size_t(nbRiver)};
    varTransport.putVar(startTrans, countTrans, Atransport.data());
  }
}

void PrintRiverInformation(FullNamelist const &eFull) {
  SingleBlock eBlINPUT = eFull.get_block("INPUT");
  std::string RiverDescriptionFile =
    eBlINPUT.get_string("RiverDescriptionFile");
  std::vector<std::string> ListTimes =
    eBlINPUT.get_list_string("ListTimes");
  int StylePrint = eBlINPUT.get_int("StylePrint");
  double maxAllowedTimeInterval =
    eBlINPUT.get_double("maxAllowedTimeInterval");

  DescriptionRiver eDescRiv = ReadRiverDescription(RiverDescriptionFile);
  int nbTime = ListTimes.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  for (int iTime = 0; iTime < nbTime; iTime++) {
    std::string eTimeStr = ListTimes[iTime];
    double eTime = CT2MJD(eTimeStr);
    std::string strPres = DATE_ConvertMjd2mystringPres(eTime);
    TransTempSalt eTTS = RetrieveTTS(eDescRiv, eTime, maxAllowedTimeInterval);
    if (StylePrint == 1)
      std::cerr << "iTime=" << iTime << " date=" << strPres
                << " transport=" << eTTS.eTransport << " temp=" << eTTS.eTemp
                << " salt=" << eTTS.eSalt << "\n";
    if (StylePrint == 2)
      std::cerr << eTTS.eTransport << " " << eTTS.eTemp << " " << eTTS.eSalt
                << "\n";
    if (StylePrint == 3)
      std::cerr << eTTS.eTransport << "\n";
    if (StylePrint == 4)
      std::cerr << eTTS.eTemp << "\n";
  }
}

// clang-format off
#endif  // SRC_OCEAN_RIVER_H_
// clang-format on
