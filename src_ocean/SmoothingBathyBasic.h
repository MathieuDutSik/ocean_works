// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_OCEAN_SMOOTHINGBATHYBASIC_H_
#define SRC_OCEAN_SMOOTHINGBATHYBASIC_H_

#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include "Model_grids.h"
#include "Triangulations.h"

MyMatrix<double> GetRoughnessFactor(MyMatrix<double> const &TheBathy,
                                    MyMatrix<int> const &ELE) {
  int nbVert = TheBathy.rows();
  MyMatrix<double> Fret(nbVert, 1);
  MyMatrix<int> Nb(nbVert, 1);
  for (int iVert = 0; iVert < nbVert; iVert++)
    Fret(iVert, 0) = 0;
  int nbEle = ELE.rows();
  for (int iEle = 0; iEle < nbEle; iEle++)
    for (int i = 0; i < 3; i++) {
      int j = NextIdx(3, i);
      int iPt = ELE(iEle, i);
      int jPt = ELE(iEle, j);
      double dep1 = TheBathy(iPt, 0);
      double dep2 = TheBathy(jPt, 0);
      double rFrac = T_abs(dep1 - dep2) / (dep1 + dep2);
      if (rFrac > Fret(iPt))
        Fret(iPt) = rFrac;
      if (rFrac > Fret(jPt))
        Fret(jPt) = rFrac;
    }
  return Fret;
}

// clang-format off
#endif  // SRC_OCEAN_SMOOTHINGBATHYBASIC_H_
// clang-format on
