#ifndef INCLUDE_ONE_DIM_INTERPOLATION_H
#define INCLUDE_ONE_DIM_INTERPOLATION_H

#include "MAT_Matrix.h"

std::vector<double>
OneDimInterpolation_vector(MyVector<double> const &VectVal,
                           MyVector<double> const &VectX,
                           std::vector<double> const &ListX) {
  int len = VectVal.size();
  double minX = VectX(0);
  //  double maxX=VectX(len-1);
  for (int i = 1; i < len; i++) {
    double eDiff = VectVal(i) - VectVal(i - 1);
    if (eDiff < 0) {
      std::cerr << "We should have VectVal increasing from 0 to len-1\n";
      throw TerminalException{1};
    }
  }
  int nbCase = ListX.size();
  std::vector<double> RetVect(nbCase);
  for (int iCase = 0; iCase < nbCase; iCase++) {
    double eVal = 0;
    double eX = ListX[iCase];
    if (eX < minX)
      eVal = VectVal(0);
    else {
      bool IsAssigned = false;
      for (int i = 1; i < len; i++) {
        if (VectX(i - 1) <= eX && eX < VectX(i)) {
          IsAssigned = true;
          double alpha1 = (eX - VectX(i - 1)) / (VectX(i) - VectX(i - 1));
          double alpha2 = (VectX(i) - eX) / (VectX(i) - VectX(i - 1));
          eVal = alpha1 * VectVal(i) + alpha2 * VectVal(i - 1);
        }
      }
      if (!IsAssigned) {
        eVal = VectVal(len - 1);
      }
    }
    RetVect[iCase] = eVal;
  }
  return RetVect;
}

#endif
