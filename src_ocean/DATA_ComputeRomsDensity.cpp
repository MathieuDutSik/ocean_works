// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Floats.h"



void print_result(std::vector<double> const& Dens, std::string const& method, std::string const& FileOut) {
  size_t siz = Dens.size();
  if (method == "line") {
    std::ofstream os(FileOut);
    for (size_t i=0; i<siz; i++)
      os << " " << Dens[i];
    os << "\n";
    return;
  }
  std::cerr << "Failed to find a relevant method\n";
  throw TerminalException{1};
}




int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    FullNamelist eFull = NAMELIST_ComputeRomsDensity();
    if (argc != 2) {
      std::cerr << "DATA_ComputeRomsDensity [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    SingleBlock const& BlockPROC = eFull.get_block("PROC");
    std::vector<double> ListTemp = BlockPROC.get_list_double("ListTemp");
    std::vector<double> ListSalt = BlockPROC.get_list_double("ListSalt");
    std::vector<double> ListDep = BlockPROC.get_list_double("ListDep");
    std::string method = BlockPROC.get_string("method");
    std::string FileOut = BlockPROC.get_string("FileOut");
    //
    long siz = ListTemp.size();
    Eigen::Tensor<double,3> Temp(1,1,siz), Salt(1,1,siz), Dep(1,1,siz);
    for (long i=0; i<siz; i++) {
      Temp(0, 0, i) = ListTemp[i];
      Salt(0, 0, i) = ListSalt[i];
      Dep(0, 0, i) = ListDep[i];
    }
    Eigen::Tensor<double, 3> Dens = ComputeDensityAnomaly(Salt, Temp, Dep);
    std::vector<double> V_dens(siz);
    for (long i=0; i<siz; i++)
      V_dens[i] = Dens(0,0,i);
    print_result(V_dens, method, FileOut);
    std::cerr << "Normal termination of DATA_ComputeRomsDensity\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in DATA_ComputeRomsDensity\n";
    exit(e.eVal);
  }
  runtime(time1);
}
