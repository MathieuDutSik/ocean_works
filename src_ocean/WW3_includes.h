#ifndef WW3_EXPORT_FUNCTION_INCLUDE
#define WW3_EXPORT_FUNCTION_INCLUDE

extern "C" {
  void write_wavewatch_header_(int *ChoiceFile, int *NX, int *NY, int *GTYPE);
  void write_wavewatch_entry_two_field_(char const *FileName, int const *TFN, int *NX, int *NY, float *U, float*V);
  void write_wavewatch_entry_one_field_(char const *FileName, int const *TFN, int *NX, int *NY, float *F);
}

#endif
