/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2021 The WASA-at-COSY Collaboration
* Aleksander K.                 2019-09
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2021-02
***********************************************/

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TInterpreter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <cassert>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TArrow.h>
#include <vector>
#include <fstream>
#include <iostream>

/////////////////////////////////////FITTING FUNCTIONS//////////////////////////////////////

//POLYNOMIALS function (for background)
Int_t bkgdParam[2];

//first degree polynomial
bkgdParam[0] = 1;
Double_t bkgdPol1(Double_t *x, Double_t *par) {
  Double_t fpol = 0.;
  Double_t c = 1.;
  for (Int_t i=0; i<bkgdParam[0]+1; i++) {
    fpol += par[i]*c;
    c*=x[0];
  };
  return fpol;
}

//second degree polynomial
bkgdParam[1] = 2;
Double_t bkgdPol2(Double_t *x, Double_t *par) {
  Double_t fpol = 0.;
  Double_t c = 1.;
  for (Int_t i = 0; i<bkgdParam[1]+1; i++) {
    fpol += par[i]*c;
    c*=x[0];
  };
  return fpol;
}

////////////////////////////////////////////////////////////////////////////////////////////

void excitationFunction() {

  //Read different ROOT files
  TFile* f[3];
  f[0] = new TFile("input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
  f[1] = new TFile("../Efficiency/output/files/Efficiency.root","READ");
  f[2] = new TFile("../LuminosityDetermination/output/files/LuminosityQuasyFreeReaction.root","READ");

  //DATA
  TH1F *hQ_signal;
  f[0]->cd("Histograms");
  hQ_signal = (TH1F*)gDirectory->Get("DATA_lev2_cut4/hQ_lev2_cut4");    //take histogram

  //Efficiency: calculated as Nacc/Ngen for sum of all Bs and Gamma ranges
  TH1F *hEfficiency;
  f[1]->cd();
  hEfficiency = (TH1F*)gDirectory->Get("hEfficiency_lev2_cut4");   //take histogram

  //Luminosity: L(Q) pd_ppn reaction with shadowing effect taken into account
  TH1F *hLuminosity;
  f[2]->cd();
  hLuminosity = (TH1F*)gDirectory->Get("hLuminosity");   //take histogram

  TF1 *fitFuncLuminosity = new TF1("fitFuncLuminosity", "[0]*x*x*x + [1]*x*x + [2]*x + [3]", -70.,30.);
  //parameters for luminosity fit function
  fitFuncLuminosity->SetParameter(0,2.38324e-05);
  fitFuncLuminosity->SetParameter(1,-7.58369e-04);
  fitFuncLuminosity->SetParameter(2,2.44746e-02);
  fitFuncLuminosity->SetParameter(3,6.55615e+01);

////////////////////////////EXCITATION FUNCTIONS - cross section////////////////////////////

  TH1F *hSignal_normLumEff = new TH1F("hSignal_normLumEff","",40,-70.,30.);
  TH1F *hSignal_normLumEff_copy;

////////////////////////////Efficiency and Luminosity correction////////////////////////////

  for(Int_t l=1; l<41; l++) {

    //Signal
    Double_t N = hQ_signal->GetBinContent(l);
    Double_t N_error = hQ_signal->GetBinError(l);
    Double_t binCenter = hQ_signal->GetBinCenter(l);

    //Efficiency
    Double_t eff = hEfficiency->GetBinContent(l);
    Double_t eff_error = hEfficiency->GetBinError(l);
    //Double_t eff_error = 0.;

    //Luminosity
    //Double_t lum = hLuminosity->GetBinContent(l);
    //Double_t lum_error = 0.;
    Double_t lum = fitFuncLuminosity->Eval(binCenter);
    //Luminosity error from the regresion
    Double_t lum_err = 0.;
    Double_t cl = 0.95;   //cl is calculated at cl = 0.95, there are errors predicted by regression line
    TFitResultPtr r = hLuminosity->Fit(fitFuncLuminosity,"SQ+");
    r->GetConfidenceIntervals(2,1,1,&binCenter,&lum_err,cl);
    Double_t lum_error = lum_err;

    //Statistical error N/(L*E)
    Double_t statError_normLumEff = TMath::Sqrt((N_error*N_error)/(lum*lum*eff*eff)+(N*N*lum_error*lum_error)/(lum*lum*lum*lum*eff*eff)+(N*N*eff_error*eff_error)/(eff*eff*eff*eff*lum*lum));

    hSignal_normLumEff->SetBinContent(l,(N/(lum*eff)));
    hSignal_normLumEff->SetBinError(l,statError_normLumEff);

  }

  //
  hSignal_normLumEff_copy = (TH1F*)hSignal_normLumEff->Clone("hSignal_normLumEff_copy");

  //Fit background
  TF1 *fitBkgdPol[2];
  fitBkgdPol[0] = new TF1("fitBkgdPol1",bkgdPol1,-70,30,bkgdParam[0]+1);
  hSignal_normLumEff->Fit("fitBkgdPol1","","",-70,30);

  fitBkgdPol[1] = new TF1("fitBkgdPol2",bkgdPol2,-70,30,bkgdParam[1]+1);
  hSignal_normLumEff->Fit("fitBkgdPol2","","",-70,30);

/////////////////////////////////////////HISTOGRAMS/////////////////////////////////////////

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(55);

  gStyle->SetTitleFont(62,"XYZ");
  gStyle->SetLabelFont(62,"XYZ");
  gStyle->SetTextFont(62);

  //Number of selected events after application of all selection criteria
  TCanvas* myCanvas00 = new TCanvas;

  Double_t Ymax00 = 1.2*(hQ_signal->GetBinContent(hQ_signal->GetMaximumBin()));

  hQ_signal->SetTitle("");
  hQ_signal->GetXaxis()->SetTitle("excess energy [MeV]");
  hQ_signal->GetXaxis()->SetTitleSize(0.06);
  hQ_signal->GetXaxis()->SetTitleOffset(0.85);
  hQ_signal->GetXaxis()->SetLabelSize(0.05);
  hQ_signal->GetXaxis()->SetRangeUser(-70.,30.);
  hQ_signal->GetYaxis()->SetTitle("events");
  hQ_signal->GetYaxis()->SetTitleSize(0.06);
  hQ_signal->GetYaxis()->SetTitleOffset(1.1);
  hQ_signal->GetYaxis()->SetLabelSize(0.05);
  hQ_signal->GetYaxis()->SetRangeUser(0.,Ymax00);

  hQ_signal->SetLineWidth(1);
  hQ_signal->SetLineColor(1);
  hQ_signal->SetMarkerColor(1);
  hQ_signal->SetMarkerSize(0.75);
  hQ_signal->SetMarkerStyle(23);
  hQ_signal->Draw("E1");

  //legend
  TLegend *myLegend00 = new TLegend(0.525, 0.815, 0.885, 0.885);
  myLegend00->SetFillStyle(0); myLegend00->SetFillColor(0); myLegend00->SetLineColor(0); myLegend00->SetTextSize(0.04);
  myLegend00->AddEntry(hQ_signal, "experimental points", "ep");
  myLegend00->Draw();

  myCanvas00->Print("output/plots/signal_dppi0.png", "png");
  myCanvas00->Print("output/plots/signal_dppi0.eps", "eps");

  //PL
  TCanvas* myCanvas00pl = new TCanvas;

  hQ_signal->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
  hQ_signal->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
  hQ_signal->Draw("E1");

  TLegend *myLegend00pl = new TLegend(0.525, 0.815, 0.885, 0.885);
  myLegend00pl->SetFillStyle(1001); myLegend00pl->SetFillColor(19); myLegend00pl->SetLineColor(1); myLegend00pl->SetTextSize(0.04); myLegend00pl->SetBorderSize(5);
  myLegend00pl->AddEntry(hQ_signal, "dane eksperymentalne", "ep");
  myLegend00pl->Draw();

  myCanvas00pl->Print("output/plots/pl/signal_dppi0_pl.png", "png");
  myCanvas00pl->Print("output/plots/pl/signal_dppi0_pl.eps", "eps");

  ////
  //Normalized events: experimental excitation function
  TCanvas* myCanvas01 = new TCanvas;

  Double_t Ymax01 = 1.2*(hSignal_normLumEff_copy->GetBinContent(hSignal_normLumEff_copy->GetMaximumBin()));

  hSignal_normLumEff_copy->SetTitle("");
  hSignal_normLumEff_copy->GetXaxis()->SetTitle("excess energy [MeV]");
  hSignal_normLumEff_copy->GetXaxis()->SetTitleSize(0.06);
  hSignal_normLumEff_copy->GetXaxis()->SetTitleOffset(0.85);
  hSignal_normLumEff_copy->GetXaxis()->SetLabelSize(0.05);
  hSignal_normLumEff_copy->GetXaxis()->SetRangeUser(-70.,30.);
  hSignal_normLumEff_copy->GetYaxis()->SetTitle("normalized events [nb]");
  hSignal_normLumEff_copy->GetYaxis()->SetTitleSize(0.06);
  hSignal_normLumEff_copy->GetYaxis()->SetTitleOffset(1.0);
  hSignal_normLumEff_copy->GetYaxis()->SetLabelSize(0.05);
  hSignal_normLumEff_copy->GetYaxis()->SetRangeUser(0.,Ymax01);

  hSignal_normLumEff_copy->SetLineWidth(1);
  hSignal_normLumEff_copy->SetLineColor(kMagenta+3);
  hSignal_normLumEff_copy->SetMarkerColor(kMagenta+3);
  hSignal_normLumEff_copy->SetMarkerSize(0.75);
  hSignal_normLumEff_copy->SetMarkerStyle(23);
  hSignal_normLumEff_copy->Draw("E1");

  fitBkgdPol[0]->SetLineColor(kOrange+1);
  fitBkgdPol[0]->SetLineStyle(1);
  fitBkgdPol[0]->SetLineWidth(1);
  fitBkgdPol[0]->Draw("same");
  fitBkgdPol[1]->SetLineColor(kMagenta-4);
  fitBkgdPol[1]->SetLineStyle(1);
  fitBkgdPol[1]->SetLineWidth(1);
  fitBkgdPol[1]->Draw("same");

  //legend
  TLegend *myLegend01 = new TLegend(0.525, 0.720, 0.885, 0.885);
  myLegend01->SetFillStyle(0); myLegend01->SetFillColor(0); myLegend01->SetLineColor(0); myLegend01->SetTextSize(0.04);
  myLegend01->AddEntry(hSignal_normLumEff_copy, "experimental points", "ep");
  myLegend01->AddEntry(fitBkgdPol[0], "polynomial  1" , "l");
  myLegend01->AddEntry(fitBkgdPol[1], "polynomial  2" , "l");
  myLegend01->Draw();

  myCanvas01->Print("output/plots/signal_dppi0_normLumEff.png", "png");
  myCanvas01->Print("output/plots/signal_dppi0_normLumEff.eps", "eps");

  //PL
  TCanvas* myCanvas01pl = new TCanvas;

  hSignal_normLumEff_copy->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
  hSignal_normLumEff_copy->GetYaxis()->SetTitle("zdarzenia znormalizowane, nb");
  hSignal_normLumEff_copy->Draw("E1");

  fitBkgdPol[0]->Draw("same");
  fitBkgdPol[1]->Draw("same");

  TLegend *myLegend01pl = new TLegend(0.525, 0.720, 0.885, 0.885);
  myLegend01pl->SetFillStyle(1001); myLegend01pl->SetFillColor(19); myLegend01pl->SetLineColor(1); myLegend01pl->SetTextSize(0.04); myLegend01pl->SetBorderSize(5);
  myLegend01pl->AddEntry(hSignal_normLumEff_copy, "dane eksperymentalne", "ep");
  myLegend01pl->AddEntry(fitBkgdPol[0], "wielomian  1" , "l");
  myLegend01pl->AddEntry(fitBkgdPol[1], "wielomian  2" , "l");
  myLegend01pl->Draw();

  myCanvas01pl->Print("output/plots/pl/signal_dppi0_normLumEff_pl.png", "png");
  myCanvas01pl->Print("output/plots/pl/signal_dppi0_normLumEff_pl.eps", "eps");

}
