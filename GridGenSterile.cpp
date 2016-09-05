#include "NueAna/MultiBinAna/GridGenSterile.h"

GridGenSterile::GridGenSterile() {
  SetOutputFile();

  // Set up DeltaChiSquare surface grid.
  SetSinSqTheta14();
  SetSinSqTheta24();

  // Number of fake experiments at each grid point of param space.
  SetNExperiments();

  /**
   * Set up frozen values for profiled params in 3+1 sterile model.
   * DeltaMSq41 is NOT a profiled param. However, treating DeltaMSq41 on the same footing as
   * profiled params allows us to easily process the DeltaChiSquare surfaces later.
   */
  SetFrozenTheta34();
  SetFrozenDelta13();
	SetFrozenDeltaEff();
  SetFrozenDelta14();
	SetFrozenDelta24();
  SetFrozenDeltaMSq41();

  /**
   * Set up 3-flavor model's params and their uncertainties.
   */
  SetTheta12();
  SetTheta23();
  SetTheta13();
  SetDeltaMSq12();
  SetAbsValDeltaMSq23();
  SetNormalHierarchy(true);
}

GridGenSterile::GridGenSterile() {
}

void GridGenSterile::SetFrozenTheta34(Double_t theta34) {
  FrozenTheta34 = theta34;
  return;
}

void GridGenSterile::SetFrozenDeltaMSq41(Double_t dmsq41) {
  FrozenDeltaMSq41 = dmsq41;
  return;
}

void GridGenSterile::SetFrozenDelta13(Double_t delta13) {
  FrozenDelta13 = delta13;
  return;
}

void GridGenSterile::SetFrozenDelta14(Double_t delta14) {
  FrozenDelta14 = delta14;
  return;
}

void GridGenSterile::SetFrozenDelta24() {
  FrozenDelta24 = FrozenDeltaEff + FrozenDelta14;
  return;
}

void GridGenSterile::SetFrozenDeltaEff(Double_t deltaeff) {
  FrozenDeltaEff = deltaeff;
  return;
}

void GridGenSterile::SetSinSqTheta14(Double_t ssqt14h, Double_t ssqth14l, Int_t ssqth14nsteps) {
  SinSqTheta14High = ssqt14h;
  SinSqTheta14Low = ssqth14l;
  SinSqTheta14NSteps = ssqth14nsteps;
  return;
}

void GridGenSterile::SetSinSqTheta24(Double_t ssqt2h, Double_t ssqth24l, Int_t ssqth24nsteps) {
  SinSqTheta24High = ssqt24h;
  SinSqTheta24Low = ssqth24l;
  SinSqTheta24NSteps = ssqth24nsteps;
  return;
}

void GridGen::RunMultiBinOscParErrsFourFrozen(string s) {
  /**
   * Check if Extrapolation2D vector is empty.
   */
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  /**
   * Initialize the Extrapolate2D.
   */
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  /**
   * Pass the fixed params to Extrapolate2D.
   */
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->SetOscPar(OscPar::kDeltaM14, FrozenDeltaMSq41);
    Extrap[ie]->SetOscPar(OscPar::kTh34,     FrozenTheta34);
    Extrap[ie]->SetOscPar(OscPar::kDelta,    FrozenDelta13);
    Extrap[ie]->SetOscPar(OscPar::kDelta24,  FrozenDelta24);
    Extrap[ie]->SetOscPar(OscPar::kDelta14,  FrozenDelta14);
    Extrap[ie]->SetOscPar(OscPar::kTh12,     Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh23,     Theta23);
    Extrap[ie]->SetOscPar(OscPar::kTh23,     Theta13);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);
    if (!NormalHier) {
      Extrap[ie]->InvertMassHierarchy();
    }
    Extrap[ie]->OscillatePrediction();
  }

  int nbins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  int nPID = Extrap[0]->GetNPID();
  vector<double> sig, bkgd;
  vector< vector<double> > nc, numucc, bnue, tau, nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;

  /**
   * Variable sig, bkgd are vectors of nbins dimension. They count the number of signal and
   * background events in each bin.
   * Variable oscparerr_offdiag is an upper triangular matrix of nbins-by-nbins dimension.
   * Variable oscparerr is a nbins-dim vector.
   * Variable like nc, numucc, bnue, tau, nue are a vector of nbisn-dim, each componen is a vector of
   * Extrap.size() dimension. Here, we only have one Extrapolate2D so they are just a vector of only
   * one number.
   */
  for (unsigned int i = 0; i < nbins; i++) {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnue.push_back(vector<double>() );
    tau.push_back(vector<double>() );
    nue.push_back(vector<double>() );
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        oscparerr_offdiag.push_back(0);
      }
    }
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      nc[i].push_back(0);
      numucc[i].push_back(0);
      bnue[i].push_back(0);
      tau[i].push_back(0);
      nue[i].push_back(0);
    }
  }

  double sinsqtheta14;
  double sinsqtheta24;
  vector<TTree *> ftree;
  vector< vector<TTree *> > ftree2;
  TTree * ttmp;
  noff = 0;
  for (unsigned int i = 0; i < nbins; i++) {
    ttmp = new TTree(Form("Bin_%i", i), Form("Bin_%i", i));
    ttmp->Branch("SinSqTheta14", &sinsqtheta14, "SinSqTheta14/D");
    ttmp->Branch("SinSqTheta24", &sinsqtheta24, "SinSqTheta24/D");
    ttmp->Branch("Signal", &sig[i], "Signal/D");
    ttmp->Branch("Background", &bkgd[i], "Background/D");
    ttmp->Branch("DNExp_DOscPars", &oscparerr[i], "DNExp_DOscPars/D");
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        ttmp->Branch(Form("Bin_%i_Bin_%i", i, k), &oscparerr_offdiag[noff], Form("Bin_%i_Bin_%i/D", i, k));
        noff++;
      }
    }
    ftree.push_back(ttmp);
    ttmp->Reset();

    // ftree2 is a vector of TTree which has Extrap.size() components.
    ftree2.push_back(vector<TTree *>() );
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      ttmp = new TTree(Form("Bin_%i_Run_%i", i, j), Form("Bin_%i_Run_%i", i, j));
      ttmp->Branch("SinSqTheta14", &sinsqtheta14, "SinSqTheta14/D");
      ttmp->Branch("SinSqTheta24", &sinsqtheta24, "SinSqTheta24/D");
      ttmp->Branch("Signal", &nue[i][j], "Signal/D");
      ttmp->Branch("NC", &nc[i][j], "NC/D");
      ttmp->Branch("NuMuCC", &numucc[i][j], "NuMuCC/D");
      ttmp->Branch("BNueCC", &bnue[i][j], "BNueCC/D");
      ttmp->Branch("NuTauCC", &tau[i][j], "NuTauCC/D");
      ftree2[i].push_back(ttmp);
      ttmp->Reset();
    }
  }

  /**
   * delnexphist is a vector of size nbins of TH1D histograms. nexp and nobs are vectors of size nbins of
   * doubles.
   */
  vector<double> nexp, nobs;
  vector<TH1D *> delnexphist;
  vector<TH1D *> delnexphist_offdiag;
  for (unsigned int i = 0; i < nbins; i++) {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i", i), "", 400, -1, 1));
    nexp.push_back(0);
    nobs.push_back(0);
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i", i, k), "", 400, -1, 1));
      }
    }
  }


  TFile * f = new TFile(gSystem->ExpandPathName(s.c_str()), "RECREATE");

  gRandom->SetSeed(0);

  double SinSq2Theta14Increment = 0;
  if (SinSq2Theta14NSteps > 0) {
    SinSq2Theta14Increment = (SinSq2Theta14High - SinSq2Theta14Low) / SinSq2Theta14NSteps;
  }
  if (SinSq2Theta24NSteps > 0) {
    SinSq2Theta24Increment = (SinSq2Theta24High - SinSq2Theta24Low) / SinSq2Theta24NSteps;
  }

  unsigned int l = 0;
  for (unsigned int ith14 = 0; ith14 < SinSq2Theta14NSteps + 1; ith14++) {
    sinsqtheta14 = (ith14 * SinSq2Theta14Increment + SinSq2Theta14Low);

    for (unsigned int ith24 = 0; ith24 < SinSq2Theta24NSteps + 1; ith24++) {
      sinsqtheta24 = (ith24 * SinSq2Theta24Increment + SinSq2Theta24Low);

      // Get nominal prediction
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetOscPar(OscPar::kDeltaM14, FrozenDeltaMSq41);
        Extrap[ie]->SetOscPar(OscPar::kTh34,     FrozenTheta34);
        Extrap[ie]->SetOscPar(OscPar::kDelta,    FrozenDelta13);
        Extrap[ie]->SetOscPar(OscPar::kDelta24,  FrozenDelta24);
        Extrap[ie]->SetOscPar(OscPar::kDelta14,  FrozenDelta14);
        Extrap[ie]->SetOscPar(OscPar::kTh12,     Theta12);
        Extrap[ie]->SetOscPar(OscPar::kTh23,     Theta23);
        Extrap[ie]->SetOscPar(OscPar::kTh23,     Theta13);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);   // Note that DeltaMSq23 is always positive
        if (!NormalHier) {
          Extrap[ie]->InvertMassHierarchy();
        }
        Extrap[ie]->OscillatePrediction();
      }

      for (unsigned int i = 0; i < nbins; i++) {
        sig[i] = 0;
        bkgd[i] = 0;
        ir = int(i / nPID);
        ip = i % nPID;
        for (ie = 0; ie < Extrap.size(); ie++) {
          bkgd[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1));
          sig[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i + 1));

          nc[i][ie] = Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip + 1, ir + 1);
          numucc[i][ie] = Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip + 1, ir + 1);
          bnue[i][ie] = Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip + 1, ir + 1);
          tau[i][ie] = Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip + 1, ir + 1);
          nue[i][ie] = Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip + 1, ir + 1);
        }
        nexp[i] = sig[i] + bkgd[i];
      }

      // Generate pseudo experiments
      noff = 0;
      for (unsigned int i = 0; i < nbins; i++) {
        delnexphist[i]->Reset();
        delnexphist[i]->SetName(Form("DeltaNexp_%i_%i_Diag_%i", ith14, ith24, i));
        for (unsigned int k = 0; k < nbins; k++) {
          if (k > i) {
            delnexphist_offdiag[noff]->Reset();
            delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_%i_OffDiag_%i_%i", ith14, ith24, i, k));
            noff++;
          }
        }
      }

      double theta12;
      double theta13;
      double theta23;
      double dm21;
      double dm32;
      for (Int_t u = 0; u < NumExpts; u++) {
        theta12 = Theta12 + AsymGaus(dTheta12_dn, dTheta12_up);
        theta13 = Theta13 + AsymGaus(dTheta13_dn, dTheta13_up);
        theta23 = Theta23 + AsymGaus(dTheta23_dn, dTheta23_up);
        dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn, dDeltaMSq12_up);
        dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn, dDeltaMSq23_up);

        for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
          Extrap[ie]->SetOscPar(OscPar::kTh12, theta12);
          Extrap[ie]->SetOscPar(OscPar::kTh13, theta13);
          Extrap[ie]->SetOscPar(OscPar::kTh23, theta23);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM12, dm21);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM23, dm32);
          if (!NormalHier) {
            Extrap[ie]->InvertMassHierarchy();
          }
          Extrap[ie]->OscillatePrediction();
        }

        noff = 0;
        for (unsigned int i = 0; i < nbins; i++) {
          nobs[i] = 0;
          for (ie = 0; ie < Extrap.size(); ie++) {
            nobs[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1));
            nobs[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i + 1));
          }
          delnexphist[i]->Fill((nobs[i] - nexp[i]) / (nexp[i]));
        }
        for (unsigned int i = 0; i < nbins; i++) {
          for (k = 0; k < nbins; k++) {
            if (k > i) {
              delnexphist_offdiag[noff]->Fill((nobs[i] - nexp[i]) * (nobs[k] - nexp[k]) / (nexp[i] * nexp[k]));
              noff++;
            }
          }
        }
      }

      noff = 0;
      for (unsigned int i = 0; i < nbins; i++) {
        oscparerr[i] = delnexphist[i]->GetRMS();
        delnexphist[i]->Write();
        for (unsigned int k = 0; k < nbins; k++) {
          if (k > i) {
            oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
            if (delnexphist_offdiag[noff]->GetMean() < 0) {
              oscparerr_offdiag[noff] = -1. * oscparerr_offdiag[noff];
            }
            delnexphist_offdiag[noff]->Write();
            noff++;
          }
        }
      }
      f->Close();

      gROOT->cd("/");

      for (unsigned int i = 0; i < nbins; i++) {
        for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
          ftree2[i][ie]->Fill();
        }
        ftree[i]->Fill();
      }

      f = new TFile(gSystem->ExpandPathName(s.c_str()), "UPDATE");

      if (l % 100 == 0) {
        cout << 100. * l / ((SinSqTheta14NSteps + 1) * (SinSqTheta24NSteps + 1)) << "% complete" << endl;
      }
      l++;
    }
  }

  f->Close();

  double nPOTNear, nPOTFar;

  TTree * paramtree = new TTree("paramtree", "paramtree");
  paramtree->Branch("nearPOT", &nPOTNear, "nearPOT/D");
  paramtree->Branch("farPOT", &nPOTFar, "farPOT/D");
  paramtree->Branch("Theta12", &Theta12, "Theta12/D");
  paramtree->Branch("Theta23", &Theta23, "Theta23/D");
  paramtree->Branch("Theta12", &Theta13, "Theta12/D");
  paramtree->Branch("DeltaMSq23", &dm32, "DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12", &DeltaMSq12, "DeltaMSq12/D");
  paramtree->Branch("Theta34", &FrozenTheta34, "Theta34/D");
  paramtree->Branch("DeltaMSq41", &FrozenDeltaMSq41, "DeltaMSq41/D");
  paramtree->Branch("Delta13", &FrozenDelta13, "Delta13/D");
	paramtree->Branch("Delta14", &FrozenDelta14, "Delta14/D");
  paramtree->Branch("Delta24", &FrozenDelta24, "Delta24/D");

  dm32 = DeltaMSq23;
  if (!NormalHier) {
    dm32 = -1. * DeltaMSq23;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  for (unsigned int i = 0; i < nbins; i++) {
    ftree[i]->Write();
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      ftree2[i][ie]->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}
