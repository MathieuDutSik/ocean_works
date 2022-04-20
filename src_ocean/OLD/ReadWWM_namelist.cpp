#include "Namelist.h"
int main(int argc, char *argv[]) {
  FullNamelist eFullWWM = NAMELIST_GetStandardWWM();
  if (argc != 2) {
    fprintf(stderr, "ReadWWM_Namelist\n");
    return -1;
  }
  std::string eFileName = argv[1];
  NAMELIST_ReadNamelistFile(eFileName, eFullWWM);
  NAMELIST_WriteNamelistFile(std::cout, eFullWWM);
}
