// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Floats.h"



void print_result(std::vector<double> const& Dens, std::string const& method, std::string const& FileOut) {
  size_t siz = Dens.size();
  if (method == "line") {
    std::ofstream os(FileOut);
    for (size_t i=0; i<siz; i++)
      os << " " << Dens[i];
  }
  std::cerr << "Failed to find a relevant method\n";
  throw TerminalException{1};
}




int main(int argc, char *argv[]) {
  srand_random_set();
  SingletonTime time1;
  try {
    FullNamelist eFull = NAMELIST_ComputeRomsDensity();
    if (argc != 2) {
      std::cerr << "DATA_ComputeRomsDensity [file.nml]\n";
      std::cerr << "with file.nml the file describing the plotting routines\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    std::map<std::string, std::vector<double>> const& ListListDoubleValues = eFull.ListBlock.at("PROC").ListListDoubleValues;
    std::vector<double> ListTemp = ListListDoubleValues.at("ListTemp");
    std::vector<double> ListSalt = ListListDoubleValues.at("ListSalt");
    std::vector<double> ListDep = ListListDoubleValues.at("ListDep");
    std::string method = eFull.ListBlock.at("PROC").ListStringValues.at("method");
    std::string FileOut = eFull.ListBlock.at("PROC").ListStringValues.at("FileOut");
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
    std::cerr << "Error in PLOT_float\n";
    exit(e.eVal);
  }
  runtime(time1);
}
