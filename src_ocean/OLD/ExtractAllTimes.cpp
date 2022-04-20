#include "Basic_netcdf.h"
int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "ReadWWM_Namelist\n");
    return -1;
  }
  std::string HisPrefix = argv[1];
  std::string eStringTime = argv[2];
  ArrayHistory eArr = NC_ReadArrayHistory(HisPrefix, eStringTime);
  NETCDF_PrintHistoryArray(std::cout, eArr);
}
