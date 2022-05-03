#ifndef SRC_OCEAN_WW3_INCLUDES_H_
#define SRC_OCEAN_WW3_INCLUDES_H_

extern "C" {
void write_wavewatch_header_(int *ChoiceFile, int *NX, int *NY, int *GTYPE);
void write_wavewatch_entry_two_field_(char const *FileName, int const *TFN,
                                      int *NX, int *NY, float *U, float *V);
void write_wavewatch_entry_one_field_(char const *FileName, int const *TFN,
                                      int *NX, int *NY, float *F);
}

// clang-format off
#endif  // SRC_OCEAN_WW3_INCLUDES_H_
// clang-format on
