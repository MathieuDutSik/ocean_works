#include "Basic_netcdf.h"
int main(int argc, char *argv[])
{
  if (argc != 2) {
    fprintf(stderr, "Read_ROMS_grid is used as\n");
    fprintf(stderr, "Read_ROMS_grid [GridFile]\n");
    return -1;
  }
  std::string eFileName=argv[1];
  GridArray GrdArr=NC_ReadRomsGridFile(eFileName);
}
