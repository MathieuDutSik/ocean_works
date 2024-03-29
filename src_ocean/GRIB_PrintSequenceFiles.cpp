// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Model_grids.h"
int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 5) {
      std::cerr << "GRIB_PrintSequenceFiles [HisPrefix] [ModelName] "
                   "[shortName] [ListFile]\n";
      std::cerr
          << "with HisPrefix     the input place where the files are located\n";
      std::cerr << "with ModelName     the kind of model used\n";
      std::cerr << "with shortName     the name of the variable of interest "
                   "(typically 10u)\n";
      std::cerr
          << " and ListFile      the filename where the data is written\n";
      return -1;
    }
    std::string HisPrefix = argv[1];
    std::string ModelName = argv[2];
    std::string shortName = argv[3];
    std::string ListFile = argv[4];
    //
    std::string GridFile = "unset";  // for grib there is no separate grid
    TripleModelDesc eTriple{ModelName, GridFile, "unset", HisPrefix, {}};
    ArrayHistory eArr = ReadArrayHistory(eTriple);
    //
    std::ofstream os(ListFile);
    int nbTime = eArr.ListTime.size();
    for (int iTime = 0; iTime < nbTime; iTime++) {
      std::vector<GRIB_MessageInfo> eList = eArr.ListListMessages[iTime];
      std::vector<GRIB_MessageInfo> ListSel;
      for (auto &eMesg : eList) {
        if (eMesg.shortName == shortName) {
          ListSel.push_back(eMesg);
        }
      }
      if (ListSel.size() > 1) {
        std::cerr << " We should have only one matching entry\n";
        throw TerminalException{1};
      }
      if (ListSel.size() == 1) {
        std::string FileName = ListSel[0].FileName;
        std::vector<std::string> LStr = STRING_Split(FileName, "/");
        std::string FileNameRed = LStr[LStr.size() - 1];
        os << FileNameRed << "\n";
      }
    }
    std::cerr << "Normal termination of GRIB_PrintSequenceFiles\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRIB_PrintSequenceFiles\n";
    exit(e.eVal);
  }
  runtime(time1);
}
