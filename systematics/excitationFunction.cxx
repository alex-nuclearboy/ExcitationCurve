/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2021 The WASA-at-COSY Collaboration
* Aleksander K.                 2019-09
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2021-09
***********************************************/

//Macro for studies of systematic uncertainties

#include <TH1F.h>
#include <TH1D.h>
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
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TArrow.h>
#include <vector>
#include <fstream>
#include <iostream>

/////////////////////////////////////FITTING FUNCTIONS//////////////////////////////////////

//POLYNOMIALS function (for background)
//first degree polynomial
Int_t bkgdParam1 = 1;
Double_t bkgdPol1(Double_t *x, Double_t *par) {
  Double_t fpol = 0.;
  Double_t c = 1.;
  for (Int_t i=0; i<bkgdParam1+1; i++) {
    fpol += par[i]*c;
    c*=x[0];
  };
  return fpol;
}

//second degree polynomial
Int_t bkgdParam2 = 2;
Double_t bkgdPol2(Double_t *x, Double_t *par) {
  Double_t fpol = 0.;
  Double_t c = 1.;
  for (Int_t i = 0; i<bkgdParam2+1; i++) {
    fpol += par[i]*c;
    c*=x[0];
  };
  return fpol;
}

////////////////////////////////////////////////////////////////////////////////////////////

void excitationFunction() {

  //Read different ROOT files
  TFile* f[3];
  f[0] = new TFile("../input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
  f[1] = new TFile("../../Efficiency/output/files/Efficiency.root","READ");
  f[2] = new TFile("../../LuminosityDetermination/output/files/LuminosityQuasyFreeReaction.root","READ");

  //DATA
  TH1D *hQ_signal[12];
  f[0]->cd("Histograms");
  hQ_signal[0] = (TH1D*)gDirectory->Get("DATA_lev2_cut4/hQ_lev2_cut4");    //take main histogram
  //systematics
  hQ_signal[1] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_1_cut4");
  hQ_signal[2] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_2_cut4");
  hQ_signal[3] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_1_1");
  hQ_signal[4] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_1_2");
  hQ_signal[5] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_2_1");
  hQ_signal[6] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_2_2");
  hQ_signal[7] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_3_1");
  hQ_signal[8] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_3_2");
  hQ_signal[9] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_4_1");
  hQ_signal[10] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_4_2");

  //Efficiency: calculated as Nacc/Ngen for sum of all Bs and Gamma ranges
  TH1D *hEfficiency[12];
  f[1]->cd();
  hEfficiency[0] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4");   //take main histogram
  //systematics
  hEfficiency[1] = (TH1D*)gDirectory->Get("hEfficiency_lev2_1_cut4");
  hEfficiency[2] = (TH1D*)gDirectory->Get("hEfficiency_lev2_2_cut4");
  hEfficiency[3] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_1_1");
  hEfficiency[4] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_1_2");
  hEfficiency[5] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_2_1");
  hEfficiency[6] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_2_2");
  hEfficiency[7] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_3_1");
  hEfficiency[8] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_3_2");
  hEfficiency[9] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_4_1");
  hEfficiency[10] = (TH1D*)gDirectory->Get("hEfficiency_lev2_cut4_4_2");

  //Luminosity: L(Q) pd_ppn reaction with shadowing effect taken into account
  TH1D *hLuminosity;
  f[2]->cd();
  hLuminosity = (TH1D*)gDirectory->Get("hLuminosity");   //take histogram

  TF1 *fitFuncLuminosity = new TF1("fitFuncLuminosity", "[0]*x*x*x + [1]*x*x + [2]*x + [3]", -70.,30.);
  //parameters for luminosity fit function
  fitFuncLuminosity->SetParameter(0,2.38324e-05);
  fitFuncLuminosity->SetParameter(1,-7.58369e-04);
  fitFuncLuminosity->SetParameter(2,2.44746e-02);
  fitFuncLuminosity->SetParameter(3,6.55615e+01);

////////////////////////////EXCITATION FUNCTIONS - cross section////////////////////////////

  TH1D *hSignal_normLumEff[12];
  for (Int_t i=0; i<12; i++){
    hSignal_normLumEff[i] = new TH1D(Form("hSignal_normLumEff_%d",i),"",40,-70.,30.);
  }
  TH1D *hSignal_normLumEff_copy;
  TH1D *hSignal_normLumEff_stat = new TH1D("ExcitationFunctionStat","",40,-70,30);
  TH1D *hSignal_normLumEff_syst = new TH1D("ExcitationFunctionSyst","",40,-70,30);

////////////////////////////Efficiency and Luminosity correction////////////////////////////

  Double_t N[12];
  Double_t N_error[12];
  Double_t binCenter[12];

  Double_t eff[12];
  Double_t eff_error[12];

  Double_t lum[12];
  Double_t lum_error[12];
  Double_t lum_err[12];
  Double_t cl = 0.95;   //cl is calculated at cl = 0.95, there are errors predicted by regression line
  TFitResultPtr r[12];

  Double_t statError_normLumEff[12];

  Double_t excit[12],systError[6];
  Double_t systErr_EC = 0.;
  Double_t systErr_IM = 0.;
  Double_t systErr_OA = 0.;
  Double_t systErr_MM = 0.;
  Double_t systErr_DM = 0.;

  ofstream resultsFile;
  resultsFile.open("Results.dat", ios::trunc);
  resultsFile<<Form("sigma  stat  syst   %%")<<endl;

  for(Int_t l=1; l<41; l++) {

    //Signal
    N[0] = hQ_signal[0]->GetBinContent(l);
    N_error[0] = hQ_signal[0]->GetBinError(l);
    binCenter[0] = hQ_signal[0]->GetBinCenter(l);

    //Efficiency
    eff[0] = hEfficiency[0]->GetBinContent(l);
    eff_error[0] = hEfficiency[0]->GetBinError(l);
    //eff_error[0] = 0.;

    //Luminosity
    //lum[0] = hLuminosity->GetBinContent(l);
    //lum_error[0] = 0.;
    lum[0] = fitFuncLuminosity->Eval(binCenter[0]);
    //Luminosity error from the regresion
    r[0] = hLuminosity->Fit(fitFuncLuminosity,"SQ+");
    r[0]->GetConfidenceIntervals(2,1,1,&binCenter[0],&lum_err[0],cl);
    lum_error[0] = lum_err[0];

    //Statistical error N/(L*E)
    statError_normLumEff[0] = TMath::Sqrt((N_error[0]*N_error[0])/(lum[0]*lum[0]*eff[0]*eff[0])+(N[0]*N[0]*lum_error[0]*lum_error[0])/(lum[0]*lum[0]*lum[0]*lum[0]*eff[0]*eff[0])+(N[0]*N[0]*eff_error[0]*eff_error[0])/(eff[0]*eff[0]*eff[0]*eff[0]*lum[0]*lum[0]));

    hSignal_normLumEff[0]->SetBinContent(l,(N[0]/(lum[0]*eff[0])));
    hSignal_normLumEff[0]->SetBinError(l,statError_normLumEff[0]);

    excit[0] = N[0]/(lum[0]*eff[0]);

    for (Int_t k=1; k<11; k++){

      //Signal
      N[k] = hQ_signal[k]->GetBinContent(l);
      N_error[k] = hQ_signal[k]->GetBinError(l);
      binCenter[k] = hQ_signal[k]->GetBinCenter(l);

      //Efficiency
      eff[k] = hEfficiency[k]->GetBinContent(l);
      eff_error[k] = hEfficiency[k]->GetBinError(l);
      //eff_error[k] = 0.;

      //Luminosity
      //lum[k] = hLuminosity->GetBinContent(l);
      //lum_error[k] = 0.;
      lum[k] = fitFuncLuminosity->Eval(binCenter[k]);
      //Luminosity error from the regresion
      r[k] = hLuminosity->Fit(fitFuncLuminosity,"SQ+");
      r[k]->GetConfidenceIntervals(2,1,1,&binCenter[k],&lum_err[k],cl);
      lum_error[k] = lum_err[k];

      //Statistical error N/(L*E)
      statError_normLumEff[k] = TMath::Sqrt((N_error[k]*N_error[k])/(lum[k]*lum[k]*eff[k]*eff[k])+(N[k]*N[k]*lum_error[k]*lum_error[k])/(lum[k]*lum[k]*lum[k]*lum[k]*eff[k]*eff[k])+(N[k]*N[k]*eff_error[k]*eff_error[k])/(eff[k]*eff[k]*eff[k]*eff[k]*lum[k]*lum[k]));

      hSignal_normLumEff[k]->SetBinContent(l,(N[k]/(lum[k]*eff[k])));
      hSignal_normLumEff[k]->SetBinError(l,statError_normLumEff[k]);

      excit[k] = N[k]/(lum[k]*eff[k]);

    }

    systError[1] = (TMath::Abs(excit[1]-excit[0]) + TMath::Abs(excit[2]-excit[0]))/2;
    systError[2] = (TMath::Abs(excit[3]-excit[0]) + TMath::Abs(excit[4]-excit[0]))/2;
    systError[3] = (TMath::Abs(excit[5]-excit[0]) + TMath::Abs(excit[6]-excit[0]))/2;
    systError[4] = (TMath::Abs(excit[7]-excit[0]) + TMath::Abs(excit[8]-excit[0]))/2;
    systError[5] = (TMath::Abs(excit[9]-excit[0]) + TMath::Abs(excit[10]-excit[0]))/2;

    //Double_t err = TMath::Sqrt(4.8*4.8 + 4.*4. + 17.*17.)*excit[0]*0.01;
    Double_t err = TMath::Sqrt(4.8*4.8 + 4.*4.)*excit[0]*0.01;
    systError[0] = TMath::Sqrt(systError[1]*systError[1] + systError[2]*systError[2] + systError[3]*systError[3] + systError[4]*systError[4] + systError[5]*systError[5])+ err;

    hSignal_normLumEff_stat->SetBinContent(l,excit[0]);
    hSignal_normLumEff_stat->SetBinError(l,statError_normLumEff[0]);

    hSignal_normLumEff_syst->SetBinContent(l,excit[0]);
    hSignal_normLumEff_syst->SetBinError(l,systError[0]);

    systErr_EC += TMath::Power(systError[1]/excit[0],2)/40;
    systErr_IM += TMath::Power(systError[2]/excit[0],2)/40;
    systErr_OA += TMath::Power(systError[3]/excit[0],2)/40;
    systErr_MM += TMath::Power(systError[4]/excit[0],2)/40;
    systErr_DM += TMath::Power(systError[5]/excit[0],2)/40;

    resultsFile<<Form("%.2f %.2f %.2f %.2f",excit[0],statError_normLumEff[0],systError[0],((systError[0]/excit[0])*100))<<endl;

  }

  resultsFile.close();

  //Systematic uncertainties
  cout<<"sigma = "<<hSignal_normLumEff[0]->Integral()<<endl;
  cout<<"syst err EC = "<<TMath::Sqrt(systErr_EC)*100<<endl;
  cout<<"syst err IM = "<<TMath::Sqrt(systErr_IM)*100<<endl;
  cout<<"syst err OA = "<<TMath::Sqrt(systErr_OA)*100<<endl;
  cout<<"syst err MM = "<<TMath::Sqrt(systErr_MM)*100<<endl;
  cout<<"syst err DM = "<<TMath::Sqrt(systErr_DM)*100<<endl;
  cout<<"syst err = "<<TMath::Sqrt(systErr_EC + systErr_IM + systErr_OA + systErr_MM + systErr_DM)*100<<endl;

  //
  hSignal_normLumEff_copy = (TH1D*)hSignal_normLumEff[0]->Clone("hSignal_normLumEff_copy");

  //create root file
  TFile* myFile = new TFile("ExcitationFunction_normLumEff.root","RECREATE");
  //write new root file
  myFile->cd();
  hSignal_normLumEff_copy->Write(Form("hNormalizedEvents"));
  for (Int_t i=0; i<11; i++){
    hSignal_normLumEff[i]->Write(Form("hNormalizedEvents_%d",i));
  }
  hSignal_normLumEff_stat->Write("hNormalizedEvents_stat"));
  hSignal_normLumEff_syst->Write("hNormalizedEvents_syst"));
  myFile->Close();

  //Fit background
  TF1 *fitBkgdPol[2];
  fitBkgdPol[0] = new TF1("fitBkgdPol1",bkgdPol1,-70,30,bkgdParam1+1);
  hSignal_normLumEff[0]->Fit("fitBkgdPol1","","",-70,30);

  fitBkgdPol[1] = new TF1("fitBkgdPol2",bkgdPol2,-70,30,bkgdParam2+1);
  hSignal_normLumEff[0]->Fit("fitBkgdPol2","","",-70,30);

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

  Double_t Ymax00 = 1.2*(hQ_signal[0]->GetBinContent(hQ_signal[0]->GetMaximumBin()));

  hQ_signal[0]->SetTitle("");
  hQ_signal[0]->GetXaxis()->SetTitle("excess energy [MeV]");
  hQ_signal[0]->GetXaxis()->SetTitleSize(0.06);
  hQ_signal[0]->GetXaxis()->SetTitleOffset(0.85);
  hQ_signal[0]->GetXaxis()->SetLabelSize(0.05);
  hQ_signal[0]->GetXaxis()->SetRangeUser(-70.,30.);
  hQ_signal[0]->GetYaxis()->SetTitle("events");
  hQ_signal[0]->GetYaxis()->SetTitleSize(0.06);
  hQ_signal[0]->GetYaxis()->SetTitleOffset(1.1);
  hQ_signal[0]->GetYaxis()->SetLabelSize(0.05);
  hQ_signal[0]->GetYaxis()->SetRangeUser(0.,Ymax00);

  hQ_signal[0]->SetLineWidth(1);
  hQ_signal[0]->SetLineColor(1);
  hQ_signal[0]->SetMarkerColor(1);
  hQ_signal[0]->SetMarkerSize(0.75);
  hQ_signal[0]->SetMarkerStyle(23);
  hQ_signal[0]->Draw("E1");

  hQ_signal[1]->SetLineWidth(1);
  hQ_signal[1]->SetLineColor(41);
  //hQ_signal[1]->Draw("same C");

  hQ_signal[2]->SetLineWidth(1);
  hQ_signal[2]->SetLineColor(41);
  //hQ_signal[2]->Draw("same C");

  hQ_signal[3]->SetLineWidth(1);
  hQ_signal[3]->SetLineColor(45);
  //hQ_signal[3]->Draw("same C");

  hQ_signal[4]->SetLineWidth(1);
  hQ_signal[4]->SetLineColor(45);
  //hQ_signal[4]->Draw("same C");

  hQ_signal[5]->SetLineWidth(1);
  hQ_signal[5]->SetLineColor(40);
  //hQ_signal[5]->Draw("same C");

  hQ_signal[6]->SetLineWidth(1);
  hQ_signal[6]->SetLineColor(40);
  //hQ_signal[6]->Draw("same C");

  hQ_signal[7]->SetLineWidth(1);
  hQ_signal[7]->SetLineColor(51);
  //hQ_signal[7]->Draw("same C");

  hQ_signal[8]->SetLineWidth(1);
  hQ_signal[8]->SetLineColor(51);
  //hQ_signal[8]->Draw("same C");

  hQ_signal[9]->SetLineWidth(1);
  hQ_signal[9]->SetLineColor(53);
  //hQ_signal[9]->Draw("same C");

  hQ_signal[10]->SetLineWidth(1);
  hQ_signal[10]->SetLineColor(53);
  //hQ_signal[10]->Draw("same C");

  //legend
  TLegend *myLegend00 = new TLegend(0.550, 0.825, 0.885, 0.885);
  myLegend00->SetFillStyle(0); myLegend00->SetFillColor(0); myLegend00->SetLineColor(0); myLegend00->SetTextSize(0.04);
  myLegend00->AddEntry(hQ_signal[0], "experimental points", "ep");
  //myLegend00->AddEntry(hQ_signal[1], "#DeltaE(PSB)-#DeltaE(SEC)", "l");
  //myLegend00->AddEntry(hQ_signal[3], "invariant mass", "l");
  //myLegend00->AddEntry(hQ_signal[5], "#theta^{CM}_{p#pi}", "l");
  //myLegend00->AddEntry(hQ_signal[7], "missing mass", "l");
  //myLegend00->AddEntry(hQ_signal[9], "deuteron momentum", "l");
  myLegend00->Draw();

  myCanvas00->Print("plots/hAcceptedDATA_dppi0_syst.png", "png");
  myCanvas00->Print("plots/hAcceptedDATA_dppi0_syst.eps", "eps");

  //PL
  TCanvas* myCanvas00pl = new TCanvas;

  hQ_signal[0]->GetXaxis()->SetTitle("\\hbox{energia dostępna [MeV]}");
  hQ_signal[0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
  hQ_signal[0]->Draw("same E1");

  //hQ_signal[1]->Draw("same C");
  //hQ_signal[2]->Draw("same C");
  //hQ_signal[3]->Draw("same C");
  //hQ_signal[4]->Draw("same C");
  //hQ_signal[5]->Draw("same C");
  //hQ_signal[6]->Draw("same C");
  //hQ_signal[7]->Draw("same C");
  //hQ_signal[8]->Draw("same C");
  //hQ_signal[9]->Draw("same C");
  //hQ_signal[10]->Draw("same C");

  TLegend *myLegend00pl = new TLegend(0.525, 0.815, 0.885, 0.885);
  myLegend00pl->SetFillStyle(1001); myLegend00pl->SetFillColor(19); myLegend00pl->SetLineColor(1); myLegend00pl->SetTextSize(0.04); myLegend00pl->SetBorderSize(5);
  myLegend00pl->AddEntry(hQ_signal[0], "dane eksperymentalne", "ep");
  //myLegend00pl->AddEntry(hQ_signal[1], "#DeltaE(PSB)-#DeltaE(SEC)", "l");
  //myLegend00pl->AddEntry(hQ_signal[3], "invariant mass", "l");
  //myLegend00pl->AddEntry(hQ_signal[5], "#theta^{CM}_{p#pi}", "l");
  //myLegend00pl->AddEntry(hQ_signal[7], "missing mass", "l");
  //myLegend00pl->AddEntry(hQ_signal[9], "deuteron momentum", "l");
  myLegend00pl->Draw();

  myCanvas00pl->Print("plots/hAcceptedDATA_dppi0_syst_pl.png", "png");
  myCanvas00pl->Print("plots/hAcceptedDATA_dppi0_syst_pl.eps", "eps");

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

  myCanvas01->Print("plots/hSignal_normLumEff.png", "png");
  myCanvas01->Print("plots/hSignal_normLumEff.eps", "eps");

  //PL
  TCanvas* myCanvas01pl = new TCanvas;

  hSignal_normLumEff_copy->GetXaxis()->SetTitle("\\hbox{energia dostępna [MeV]}");
  hSignal_normLumEff_copy->GetYaxis()->SetTitle("zdarzenia znormalizowane [nb]");
  hSignal_normLumEff_copy->Draw("E1");

  fitBkgdPol[0]->Draw("same");
  fitBkgdPol[1]->Draw("same");

  TLegend *myLegend01pl = new TLegend(0.525, 0.720, 0.885, 0.885);
  myLegend01pl->SetFillStyle(1001); myLegend01pl->SetFillColor(19); myLegend01pl->SetLineColor(1); myLegend01pl->SetTextSize(0.04); myLegend01pl->SetBorderSize(5);
  myLegend01pl->AddEntry(hSignal_normLumEff_copy, "dane eksperymentalne", "ep");
  myLegend01pl->AddEntry(fitBkgdPol[0], "wielomian  1" , "l");
  myLegend01pl->AddEntry(fitBkgdPol[1], "wielomian  2" , "l");
  myLegend01pl->Draw();

  myCanvas01pl->Print("plots/hSignal_normLumEff_pl.png", "png");
  myCanvas01pl->Print("plots/hSignal_normLumEff_pl.eps", "eps");
/*
  //
  TCanvas* myCanvas02 = new TCanvas;

  Double_t Ymax02 = 1.2*(hSignal_normLumEff_copy->GetBinContent(hSignal_normLumEff_copy->GetMaximumBin()));

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
  hSignal_normLumEff_copy->GetYaxis()->SetRangeUser(0.,Ymax02);

  hSignal_normLumEff_copy->SetLineWidth(1);
  hSignal_normLumEff_copy->SetLineColor(kMagenta+3);
  hSignal_normLumEff_copy->SetMarkerColor(kMagenta+3);
  hSignal_normLumEff_copy->SetMarkerSize(0.75);
  hSignal_normLumEff_copy->SetMarkerStyle(23);
  hSignal_normLumEff_copy->Draw("E1");

  hSignal_normLumEff[1]->SetLineWidth(1);
  hSignal_normLumEff[1]->SetLineColor(41);
  hSignal_normLumEff[1]->Draw("same C");

  hSignal_normLumEff[2]->SetLineWidth(1);
  hSignal_normLumEff[2]->SetLineColor(41);
  hSignal_normLumEff[2]->Draw("same C");

  hSignal_normLumEff[3]->SetLineWidth(1);
  hSignal_normLumEff[3]->SetLineColor(45);
  hSignal_normLumEff[3]->Draw("same C");

  hSignal_normLumEff[4]->SetLineWidth(1);
  hSignal_normLumEff[4]->SetLineColor(45);
  hSignal_normLumEff[4]->Draw("same C");

  hSignal_normLumEff[5]->SetLineWidth(1);
  hSignal_normLumEff[5]->SetLineColor(40);
  hSignal_normLumEff[5]->Draw("same C");

  hSignal_normLumEff[6]->SetLineWidth(1);
  hSignal_normLumEff[6]->SetLineColor(40);
  hSignal_normLumEff[6]->Draw("same C");

  hSignal_normLumEff[7]->SetLineWidth(1);
  hSignal_normLumEff[7]->SetLineColor(51);
  hSignal_normLumEff[7]->Draw("same C");

  hSignal_normLumEff[8]->SetLineWidth(1);
  hSignal_normLumEff[8]->SetLineColor(51);
  hSignal_normLumEff[8]->Draw("same C");

  hSignal_normLumEff[9]->SetLineWidth(1);
  hSignal_normLumEff[9]->SetLineColor(53);
  hSignal_normLumEff[9]->Draw("same C");

  hSignal_normLumEff[10]->SetLineWidth(1);
  hSignal_normLumEff[10]->SetLineColor(53);
  hSignal_normLumEff[10]->Draw("same C");

  //legend
  TLegend *myLegend02 = new TLegend(0.525, 0.610, 0.885, 0.885);
  myLegend02->SetFillStyle(0); myLegend02->SetFillColor(0); myLegend02->SetLineColor(0); myLegend02->SetTextSize(0.04);
  myLegend02->AddEntry(hSignal_normLumEff_copy, "experimental points", "ep");
  myLegend02->AddEntry(hSignal_normLumEff[1], "#DeltaE(PSB)-#DeltaE(SEC)", "l");
  myLegend02->AddEntry(hSignal_normLumEff[3], "invariant mass", "l");
  myLegend02->AddEntry(hSignal_normLumEff[5], "#theta^{CM}_{p#pi}", "l");
  myLegend02->AddEntry(hSignal_normLumEff[7], "missing mass", "l");
  myLegend02->AddEntry(hSignal_normLumEff[9], "deuteron momentum", "l");
  myLegend02->Draw();

  myCanvas02->Print("plots/signal_dppi0_normLumEff_diff.png", "png");
  myCanvas02->Print("plots/signal_dppi0_normLumEff_diff.eps", "eps");
*/
  //
  TCanvas* myCanvas03 = new TCanvas;

  //hSignal_normLumEff_syst->SetTitle("");
  hSignal_normLumEff_syst->GetXaxis()->SetTitle("excess energy [MeV]");
  hSignal_normLumEff_syst->GetXaxis()->SetTitleSize(0.06);
  hSignal_normLumEff_syst->GetXaxis()->SetTitleOffset(0.85);
  hSignal_normLumEff_syst->GetXaxis()->SetLabelSize(0.05);
  hSignal_normLumEff_syst->GetXaxis()->SetRangeUser(-70.,30.);
  hSignal_normLumEff_syst->GetYaxis()->SetTitle("normalized events [nb]");
  hSignal_normLumEff_syst->GetYaxis()->SetTitleSize(0.06);
  hSignal_normLumEff_syst->GetYaxis()->SetTitleOffset(1.0);
  hSignal_normLumEff_syst->GetYaxis()->SetLabelSize(0.05);
  hSignal_normLumEff_syst->GetYaxis()->SetRangeUser(0.,Ymax01);

  //hSignal_normLumEff_syst->SetLineWidth(2);
  hSignal_normLumEff_syst->SetLineColor(30);
  hSignal_normLumEff_syst->SetMarkerColor(30);
  //hSignal_normLumEff_syst->SetFillColor(2);
  //hSignal_normLumEff_syst->SetFillStyle(3001);
  hSignal_normLumEff_syst->SetMarkerSize(0.75);
  hSignal_normLumEff_syst->SetMarkerStyle(2);

  hSignal_normLumEff_syst->Draw("E1");

  //hSignal_normLumEff_stat->SetLineWidth(2);
  hSignal_normLumEff_stat->SetLineColor(kMagenta+3);
  hSignal_normLumEff_stat->SetMarkerColor(kMagenta+3);
  hSignal_normLumEff_stat->SetMarkerSize(0.75);
  hSignal_normLumEff_stat->SetMarkerStyle(23);
  hSignal_normLumEff_stat->Draw("same E1");

  fitBkgdPol[0]->Draw("same");
  fitBkgdPol[1]->Draw("same");

  //legend
  TLegend *myLegend03 = new TLegend(0.525, 0.670, 0.885, 0.885);
  myLegend03->SetFillStyle(0); myLegend03->SetFillColor(0); myLegend03->SetLineColor(0); myLegend03->SetTextSize(0.04);
  myLegend03->AddEntry(hSignal_normLumEff_stat, "experimental points", "pe");
  myLegend03->AddEntry(hSignal_normLumEff_syst, "systematic uncertainties", "pe");
  myLegend03->AddEntry(fitBkgdPol[0], "polynomial  1" , "l");
  myLegend03->AddEntry(fitBkgdPol[1], "polynomial  2" , "l");
  myLegend03->Draw();

  myCanvas03->Print("plots/hSignal_normLumEff_syst.png","png");
  myCanvas03->Print("plots/hSignal_normLumEff_syst.eps","eps");

  //PL
  TCanvas* myCanvas03pl = new TCanvas;

  hSignal_normLumEff_syst->GetXaxis()->SetTitle("\\hbox{energia dostępna [MeV]}");
  hSignal_normLumEff_syst->GetYaxis()->SetTitle("zdarzenia znormalizowane [nb]");
  hSignal_normLumEff_syst->Draw("E1");

  hSignal_normLumEff_stat->Draw("same E1");

  fitBkgdPol[0]->Draw("same");
  fitBkgdPol[1]->Draw("same");

  TLegend *myLegend03pl = new TLegend(0.515, 0.665, 0.885, 0.885);
  myLegend03pl->SetFillStyle(1001); myLegend03pl->SetFillColor(19); myLegend03pl->SetLineColor(1); myLegend03pl->SetTextSize(0.04); myLegend03pl->SetBorderSize(5);
  myLegend03pl->AddEntry(hSignal_normLumEff_stat, "dane eksperymentalne", "ep");
  myLegend03pl->AddEntry(hSignal_normLumEff_syst, "\\hbox{niepewności systemat.}", "ep");
  myLegend03pl->AddEntry(fitBkgdPol[0], "wielomian  1" , "l");
  myLegend03pl->AddEntry(fitBkgdPol[1], "wielomian  2" , "l");
  myLegend03pl->Draw();

  myCanvas03pl->Print("plots/signal_dppi0_normLumEff_syst_pl.png", "png");
  myCanvas03pl->Print("plots/signal_dppi0_normLumEff_syst_pl.eps", "eps");

}
