#ifndef GridGenSterile_h
#define GridGenSterile_h

#include "TROOT.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "OscProb/OscCalc.h"
#include "NueAna/Extrapolation/Background.h"
#include "TRandom.h"
#include "NueAna/MultiBinAna/GridGen.h"

using namespace std;

class GridGenSterile : public GridGen {
public:
  GridGenSterile();
  virtual ~GridGenSterile();

  void SetFrozenDeltaMSq41(Double_t dmsq41 = 1);     // [eV2]
  void SetFrozenTheta34(Double_t theta34 = 0);     // [Pi]
  void SetFrozenDelta13(Double_t delta13 = 0);
  void SetFrozenDelta14(Double_t delta14 = 0);
  void SetFrozenDelta24();
  void SetFrozenDeltaEff(Double_t deltaeff = 0);

  void SetSinSqTheta14(Double_t ssqth14h = 0, Double_t ssqth14l = 0.5, Int_t ssqth14nsteps = 100);
  void SetSinSqTheta24(Double_t ssqth24h = 0, Double_t ssqth24l = 0.5, Int_t ssqth24nsteps = 100);

  void RunMultiBinOscParErrsFourFrozen(string s = "OscParErrDistributions.root");

protected:
  Double_t FrozenDeltaMSq41;
  Double_t FrozenTheta34;
  Double_t FrozenDelta13;
  Double_t FrozenDeltaEff;
  Double_t FrozenDelta14;
  Double_t FrozenDelta24;

  Double_t SinSqTheta14High;
  Double_t SinSqTheta14Low;
  Double_t SinSqTheta24High;
  Double_t SinSqTheta24Low;
  Int_t SinSqTheta14NSteps;
  Int_t SinSqTheta24NSteps;
};

#endif /* ifndef GridGenSterile_h */
