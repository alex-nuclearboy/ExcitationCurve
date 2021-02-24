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

//Breit–Wigner distribution (for signal)
Double_t signalBW(Double_t *x, Double_t *par) {
  Double_t funcBW = 0.;
  funcBW = (par[0]*(par[1]/2)*(par[1]/2))/((x[0]-par[2])*(x[0]-par[2])+(par[1]/2)*(par[1]/2));
  return funcBW;
}

//Total function: BW + polynomial 1
Double_t totalFuncBWpol1(Double_t *x, Double_t *par) {
  //BW
  Double_t funcBW = 0.;
  funcBW = (par[0]*(par[1]/2)*(par[1]/2))/((x[0]-par[2])*(x[0]-par[2])+(par[1]/2)*(par[1]/2));
  //POL1
  Double_t fpol = 0.;
  Double_t c = 1.;
  for(Int_t i=0; i<bkgdParam1+1; i++) {
    fpol += par[i+3]*c;
    c*=x[0];
  };
  return funcBW + fpol;
}

//Total function: BW + polynomial 2
Double_t totalFuncBWpol2(Double_t *x, Double_t *par) {
  //BW
  Double_t funcBW = 0.;
  funcBW = (par[0]*(par[1]/2)*(par[1]/2))/((x[0]-par[2])*(x[0]-par[2])+(par[1]/2)*(par[1]/2));
  //POL2
  Double_t fpol = 0.;
  Double_t c = 1.;
  for (Int_t i=0; i<bkgdParam2+1; i++) {
    fpol += par[i+3]*c;
    c*=x[0];
  };
  return funcBW + fpol;
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
  fitBkgdPol[0] = new TF1("fitBkgdPol1",bkgdPol1,-70,30,bkgdParam1+1);
  hSignal_normLumEff->Fit("fitBkgdPol1","","",-70,30);

  fitBkgdPol[1] = new TF1("fitBkgdPol2",bkgdPol2,-70,30,bkgdParam2+1);
  hSignal_normLumEff->Fit("fitBkgdPol2","","",-70,30);

////////////////////////////////////Fitting with BW + pol////////////////////////////////////

  TF1 *fitBWpol1[51][41],*fitPol1_start[51][41],*fitPol1[51][41],*fitBW1[51][41];
  TF1 *fitBWpol2[51][41],*fitPol2_start[51][41],*fitPol2[51][41],*fitBW2[51][41];
  Double_t a[2][51][41],b[2][51][41],c[2][51][41],d[2][51][41];
  Double_t delta_a[2][51][41],delta_b[2][51][41],delta_c[2][51][41],delta_d[2][51][41];

  //Number of parameters
  Int_t signalParamFit = 3;
  Int_t bkgdParamFit[2];
  bkgdParamFit[0] = bkgdParam1 + 1;
  bkgdParamFit[1] = bkgdParam2 + 1;
  Int_t totalParamFit[2];
  totalParamFit[0] = bkgdParamFit[0] + signalParamFit;
  totalParamFit[1] = bkgdParamFit[1] + signalParamFit;

  Double_t beginFit[2];
  Double_t endFit[2];
  beginFit[0] = -65.; //-58.5
  endFit[0] = 27.5;   //30
  beginFit[1] = -65.; //-60
  endFit[1] = 27.5;   //30

  Double_t chi2_red[2][51][41];

  Double_t amplBW[2];
  Double_t delta_amplBW[2];

  TH1F *hAmplitude[2];
  hAmplitude[0] = new TH1F("hAmplitude_1","",46,4.5,50.5);
  hAmplitude[1] = new TH1F("hAmplitude_2","",46,4.5,50.5);

  Double_t XS_CL90[2][51][41];

  TH1F *hXS_uppLimit[41];
  TH1F *hXS_uppLimit_syst[41];
  for (int k = 0; k < 41; k++) {
      hXS_uppLimit[k]=new TH1F(Form("hXS_uppLimit_Bs%d",k),"",46,4.5,50.5);
      hXS_uppLimit_syst[k]=new TH1F(Form("hXS_uppLimit_syst_Bs%d",k),"",46,4.5,50.5);
  }
  TH2F *hXS_uppLimit_2D = new TH2F("hXS_uppLimit_2D","",41,-40.5,0.5,46,4.5,50.5);;

  ////////
  for (Int_t Bs = 0; Bs < 41; Bs++) {
    for (Int_t Gamma = 5; Gamma < 51; Gamma++) {

      //POLYNOMIAL 1
      //background fit
      fitPol1_start[Gamma][Bs] = new TF1(Form("fitPol1_start_%d_%d",Gamma,Bs),bkgdPol1,-70.,30.,bkgdParamFit[0]);
      //fit
      hSignal_normLumEff->Fit(Form("fitPol1_start_%d_%d",Gamma,Bs),"","",beginFit[0],endFit[0]);

      //BW + pol1
      fitBWpol1[Gamma][Bs] = new TF1(Form("fitBWpol1_%d_%d",Gamma,Bs),totalFuncBWpol1,-70.,30.,totalParamFit[0]);
      fitBWpol1[Gamma][Bs]->SetParameter(0,0); //Amplitude
      fitBWpol1[Gamma][Bs]->FixParameter(1,Gamma);
      fitBWpol1[Gamma][Bs]->FixParameter(2,-Bs);

      for(Int_t k=3; k<totalParamFit[0]; k++) {
        fitBWpol1[Gamma][Bs]->SetParameter(k,fitPol1_start[Gamma][Bs]->GetParameter(k-3));
      }

      //fit
      hSignal_normLumEff->Fit(Form("fitBWpol1_%d_%d",Gamma,Bs),"","",beginFit[0],endFit[0]);

      //signal parameters (from total function)
      fitBW1[Gamma][Bs] = new TF1(Form("fitBW1_%d_%d",Gamma,Bs),signalBW,-70.,30.,signalParamFit);
      for(Int_t j=0; j<signalParamFit; j++) {
        fitBW1[Gamma][Bs]->SetParameter(j,fitBWpol1[Gamma][Bs]->GetParameter(j));
      }

      //background (parameters from total function)
      fitPol1[Gamma][Bs] = new TF1(Form("fitPol1_%d_%d",Gamma,Bs),bkgdPol1,-70.,30.,bkgdParamFit[0]);
      for(Int_t i=0; i<bkgdParamFit[0]; i++) {
        fitPol1[Gamma][Bs]->SetParameter(i,fitBWpol1[Gamma][Bs]->GetParameter(i+3));
      }

      //Chi^2 calculation for fitPol1
      //TFitResultPtr fitResultPol1 = hSignal_normLumEff->Fit(Form("fitBWpol1_%d_%d",Gamma,Bs),"s0");
      //Double_t chi2_pol1 = fitResultPol1->Chi2();
      //chi2_red[0][Gamma][Bs] = chi2_pol1/(40-totalParamFit[0]);

      /*
      Double_t chi2_pol1 = 0.;
      for(Int_t l = 1; l < 41; l++) {
        Double_t N_norm = hSignal_normLumEff->GetBinContent(l);
        Double_t N_norm_err = hSignal_normLumEff->GetBinError(l);
        Double_t binCenter_norm = hSignal_normLumEff->GetBinCenter(l);
        Double_t fittingValue = fitBWpol1[Gamma][Bs]->Eval(binCenter_norm);
        Double_t N_fit = N_norm - fittingValue;
        Double_t N_fit_sqr = N_fit*N_fit;
        Double_t N_fit_err = N_norm_err*N_norm_err;
        chi2_pol1 = chi2_pol1 + N_fit_sqr/N_fit_err;
      }
      chi2_red[0][Gamma][Bs] = chi2_pol1/(40-totalParamFit[0]);
      */

      //POLYNOMIAL 2
      //background fit
      fitPol2_start[Gamma][Bs] = new TF1(Form("fitPol2_start_%d_%d",Gamma,Bs),bkgdPol2,-70.,30.,bkgdParamFit[1]);
      //fit
      hSignal_normLumEff->Fit(Form("fitPol2_start_%d_%d",Gamma,Bs),"","",beginFit[1],endFit[1]);

      //BW + pol1
      fitBWpol2[Gamma][Bs] = new TF1(Form("fitBWpol2_%d_%d",Gamma,Bs),totalFuncBWpol2,-70.,30.,totalParamFit[1]);
      fitBWpol2[Gamma][Bs]->SetParameter(0,0); //Amplitude
      fitBWpol2[Gamma][Bs]->FixParameter(1,Gamma);
      fitBWpol2[Gamma][Bs]->FixParameter(2,-Bs);

      for(Int_t k=3; k<totalParamFit[1]; k++) {
        fitBWpol2[Gamma][Bs]->SetParameter(k,fitPol2_start[Gamma][Bs]->GetParameter(k-3));
      }

      //fit
      hSignal_normLumEff->Fit(Form("fitBWpol2_%d_%d",Gamma,Bs),"","",beginFit[1],endFit[1]);

      //signal parameters (from total function)
      fitBW2[Gamma][Bs] = new TF1(Form("fitBW2_%d_%d",Gamma,Bs),signalBW,-70.,30.,signalParamFit);
      for(Int_t j=0; j<signalParamFit; j++) {
        fitBW2[Gamma][Bs]->SetParameter(j,fitBWpol2[Gamma][Bs]->GetParameter(j));
      }

      //background (parameters from total function)
      fitPol2[Gamma][Bs] = new TF1(Form("fitPol2_%d_%d",Gamma,Bs),bkgdPol2,-70.,30.,bkgdParamFit[1]);
      for(Int_t i=0; i<bkgdParamFit[1]; i++) {
        fitPol2[Gamma][Bs]->SetParameter(i,fitBWpol2[Gamma][Bs]->GetParameter(i+3));
      }

      //Chi^2 calculation for fitPol2
      //TFitResultPtr fitResultPol2 = hSignal_normLumEff->Fit(Form("fitBWpol2_%d_%d",Gamma,Bs),"s0");
      //Double_t chi2_pol2 = fitResultPol2->Chi2();
      //chi2_red[1][Gamma][Bs] = chi2_pol2/(40-totalParamFit[1]);

      /*
      Double_t chi2_pol2 = 0.;
      for(Int_t l = 1; l < 41; l++) {
        Double_t N_norm = hSignal_normLumEff->GetBinContent(l);
        Double_t N_norm_err = hSignal_normLumEff->GetBinError(l);
        Double_t binCenter_norm = hSignal_normLumEff->GetBinCenter(l);
        Double_t fittingValue = fitBWpol2[Gamma][Bs]->Eval(binCenter_norm);
        Double_t N_fit = N_norm - fittingValue;
        Double_t N_fit_sqr = N_fit*N_fit;
        Double_t N_fit_err = N_norm_err*N_norm_err;
        chi2_pol2 = chi2_pol2 + N_fit_sqr/N_fit_err;
      }
      chi2_red[1][Gamma][Bs] = chi2_pol2/(40-totalParamFit[1]);
      */

/////////////////////////////////////////UPPER LIMIT/////////////////////////////////////////

      //parameters: BW + pol1
      a[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParameter(0);  //Amplitude
      b[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParameter(3);
      c[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParameter(4);

      delta_a[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParError(0);
      delta_b[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParError(3);
      delta_c[0][Gamma][Bs] = fitBWpol1[Gamma][Bs]->GetParError(4);

      //parameters: BW + pol2
      a[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParameter(0);  //Amplitude
      b[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParameter(3);
      c[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParameter(4);
      d[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParameter(5);

      delta_a[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParError(0);
      delta_b[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParError(3);
      delta_c[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParError(4);
      delta_d[1][Gamma][Bs] = fitBWpol2[Gamma][Bs]->GetParError(5);

      //
      Double_t scaleFactor = 1.;
      amplBW[0] = a[0][Gamma][Bs]*scaleFactor;
      delta_amplBW[0] = delta_a[0][Gamma][Bs]*scaleFactor;
      amplBW[1] = a[1][Gamma][Bs]*scaleFactor;
      delta_amplBW[1] = delta_a[1][Gamma][Bs]*scaleFactor;

      //XS (CL=90%), k=1.64485
      XS_CL90[0][Gamma][Bs] = 1.64485*delta_amplBW[0]; //pol1
      XS_CL90[1][Gamma][Bs] = 1.64485*delta_amplBW[1]; //pol2

      Double_t XS_CL90_average = 0.5*(XS_CL90[0][Gamma][Bs] + XS_CL90[1][Gamma][Bs]);

      Double_t systErr_fit = 0.5*(TMath::Abs(XS_CL90[0][Gamma][Bs] - XS_CL90[1][Gamma][Bs]));
      Double_t systErr_fit_percent = 100*(systErr_fit/XS_CL90_average);

      //Systematic errors:
      //4.8 - syst lum
      //4.0 - norm lum
      //1.0 - syst Fermi mom distr
      //8.5 - syst from different cuts
      //17.0 - N* distribution

      Double_t systErr_total_percent = TMath::Sqrt(4.8*4.8 + 4.*4. + 1*1 + systErr_fit_percent*systErr_fit_percent + 8.5*8.5 + 17.*17.);
      Double_t systErr_total = 0.01*systErr_total_percent*XS_CL90_average;

      ofstream newFile;
      newFile.open("output/upperLimit.dat", ios::app);
      newFile<<Form("%d %d %g %g %g %g %g",Bs,Gamma,XS_CL90[0][Gamma][Bs],XS_CL90[1][Gamma][Bs],XS_CL90_average,systErr_total,systErr_total_percent)<<endl;
      newFile.close();

      Int_t binNumber = Gamma - 4;

      hAmplitude[0]->SetBinContent(binNumber,amplBW[0]);
      hAmplitude[0]->SetBinError(binNumber,delta_amplBW[0]);
      hAmplitude[1]->SetBinContent(binNumber,amplBW[1]);
      hAmplitude[1]->SetBinError(binNumber,delta_amplBW[1]);

      hXS_uppLimit[Bs]->SetBinContent(binNumber,XS_CL90_average);
      hXS_uppLimit[Bs]->SetBinError(binNumber,0.);

      hXS_uppLimit_syst[Bs]->SetBinContent(binNumber,systErr_total);
      hXS_uppLimit_syst[Bs]->SetBinError(binNumber,0.);

      hXS_uppLimit_2D->Fill(-Bs,Gamma,XS_CL90_average);

    }

  }

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
  //myCanvas00->Print("output/plots/signal_dppi0.eps", "eps");

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
  //myCanvas01->Print("output/plots/signal_dppi0_normLumEff.eps", "eps");

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

  //FITTING
  TCanvas* myCanvas02a = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {
    for(Int_t Gamma=5; Gamma<55; Gamma=Gamma+5) {

      //hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV, #chi^{2} = %.2f",-Bs,Gamma,Chi2_red1[Gamma][Bs]));
      hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV",-Bs,Gamma));
      hSignal_normLumEff_copy->GetXaxis()->SetTitle("excess energy [MeV]");
      hSignal_normLumEff_copy->GetXaxis()->SetTitleSize(0.06);
      hSignal_normLumEff_copy->GetXaxis()->SetTitleOffset(0.85);
      hSignal_normLumEff_copy->GetXaxis()->SetLabelSize(0.05);
      hSignal_normLumEff_copy->GetXaxis()->SetRangeUser(-65.,30.);
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

      fitBWpol1[Gamma][Bs]->SetLineWidth(1);
      fitBWpol1[Gamma][Bs]->SetLineColor(kCyan-3);
      fitBWpol1[Gamma][Bs]->SetLineStyle(1);
      fitPol1[Gamma][Bs]->SetLineWidth(1);
      fitPol1[Gamma][Bs]->SetLineColor(kOrange+1);
      fitPol1[Gamma][Bs]->SetLineStyle(1);

      fitPol1[Gamma][Bs]->Draw("same");
      fitBWpol1[Gamma][Bs]->Draw("same");

      //legend
      TLegend *myLegend02a = new TLegend(0.525, 0.720, 0.885, 0.885);
      myLegend02a->SetFillStyle(0); myLegend02a->SetFillColor(0); myLegend02a->SetLineColor(0); myLegend02a->SetTextSize(0.04);
      myLegend02a->AddEntry(hSignal_normLumEff_copy, "experimental points", "ep");
      myLegend02a->AddEntry(fitBWpol1[Gamma][Bs], "BW + polynomial 1", "l");
      myLegend02a->AddEntry(fitPol1[Gamma][Bs], "polynomial 1" , "l");
      myLegend02a->Draw();

      myCanvas02a->Print(Form("output/plots/signal_dppi0_normLumEff_fitBWpol1_Bs%d_G%d.png",Bs,Gamma),"png");
    }
  }

  TCanvas* myCanvas02b = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {
    for(Int_t Gamma=5; Gamma<55; Gamma=Gamma+5) {

      //hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV, #chi^{2} = %.2f",-Bs,Gamma,Chi2_red2[Gamma][Bs]));
      hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV",-Bs,Gamma));
      hSignal_normLumEff_copy->Draw("E1");

      fitBWpol2[Gamma][Bs]->SetLineWidth(1);
      fitBWpol2[Gamma][Bs]->SetLineColor(kCyan-3);
      fitBWpol2[Gamma][Bs]->SetLineStyle(1);

      fitPol2[Gamma][Bs]->SetLineWidth(1);
      fitPol2[Gamma][Bs]->SetLineColor(kMagenta-4);
      fitPol2[Gamma][Bs]->SetLineStyle(1);

      fitPol2[Gamma][Bs]->Draw("same");
      fitBWpol2[Gamma][Bs]->Draw("same");


      //legend
      TLegend *myLegend02b = new TLegend(0.525, 0.720, 0.885, 0.885);
      myLegend02b->SetFillStyle(0); myLegend02b->SetFillColor(0); myLegend02b->SetLineColor(0); myLegend02b->SetTextSize(0.04);
      myLegend02b->AddEntry(hSignal_normLumEff_copy, "experimental points", "ep");
      myLegend02b->AddEntry(fitBWpol2[Gamma][Bs], "BW + polynomial 2", "l");
      myLegend02b->AddEntry(fitPol2[Gamma][Bs], "polynomial 2" , "l");
      myLegend02b->Draw();

      myCanvas02b->Print(Form("output/plots/signal_dppi0_normLumEff_fitBWpol2_Bs%d_G%d.png",Bs,Gamma),"png");

    }
  }

  TCanvas* myCanvas02a_pl = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {
    for(Int_t Gamma=5; Gamma<55; Gamma=Gamma+5) {

      //hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV, #chi^{2} = %.2f",-Bs,Gamma,Chi2_red1[Gamma][Bs]));
      hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV",-Bs,Gamma));
      hSignal_normLumEff_copy->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
      hSignal_normLumEff_copy->GetYaxis()->SetTitle("zdarzenia znormalizowane, nb");
      hSignal_normLumEff_copy->Draw("E1");

      fitPol1[Gamma][Bs]->Draw("same");
      fitBWpol1[Gamma][Bs]->Draw("same");

      //legend
      TLegend *myLegend02a_pl = new TLegend(0.525, 0.720, 0.885, 0.885);
      myLegend02a_pl->SetFillStyle(1001); myLegend02a_pl->SetFillColor(19); myLegend02a_pl->SetLineColor(1); myLegend02a_pl->SetTextSize(0.04); myLegend02a_pl->SetBorderSize(5);
      myLegend02a_pl->AddEntry( hSignal_normLumEff_copy, "dane eksperymentalne", "ep");
      myLegend02a_pl->AddEntry( fitBWpol1[Gamma][Bs], "BW + wielomian 1", "l");
      myLegend02a_pl->AddEntry( fitPol1[Gamma][Bs], "wielomian 1" , "l");
      myLegend02a_pl->Draw();

      myCanvas02a_pl->Print(Form("output/plots/pl/signal_dppi0_normLumEff_fitBWpol1_Bs%d_G%d_pl.png",Bs,Gamma),"png");

      if(Bs == 30 && Gamma == 15){
        myCanvas02a_pl->Print("output/plots/pl/signal_dppi0_normLumEff_fitBWpol1_Bs30_G15_pl.eps","eps");
      }

    }
  }

  TCanvas* myCanvas02b_pl = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {
    for(Int_t Gamma=5; Gamma<55; Gamma=Gamma+5) {

      //hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV, #chi^{2} = %.2f",-Bs,Gamma,Chi2_red2[Gamma][Bs]));
      hSignal_normLumEff_copy->SetTitle(Form("B_{s} = %d MeV, #Gamma = %d MeV",-Bs,Gamma));
      hSignal_normLumEff_copy->Draw("E1");

      fitPol2[Gamma][Bs]->Draw("same");
      fitBWpol2[Gamma][Bs]->Draw("same");

      //legend
      TLegend *myLegend02b_pl = new TLegend(0.525, 0.720, 0.885, 0.885);
      myLegend02b_pl->SetFillStyle(1001); myLegend02b_pl->SetFillColor(19); myLegend02b_pl->SetLineColor(1); myLegend02b_pl->SetTextSize(0.04); myLegend02b_pl->SetBorderSize(5);
      myLegend02b_pl->AddEntry( hSignal_normLumEff_copy, "dane eksperymentalne", "ep");
      myLegend02b_pl->AddEntry( fitBWpol2[Gamma][Bs], "BW + wielomian 2", "l");
      myLegend02b_pl->AddEntry( fitPol2[Gamma][Bs], "wielomian 2" , "l");
      myLegend02b_pl->Draw();

      myCanvas02b_pl->Print(Form("output/plots/pl/signal_dppi0_normLumEff_fitBWpol2_Bs%d_G%d_pl.png",Bs,Gamma),"png");

      if(Bs == 5 && Gamma == 40){
        myCanvas02b_pl->Print("output/plots/pl/signal_dppi0_normLumEff_fitBWpol2_Bs5_G40_pl.eps","eps");
      }

    }
  }

  //Amplitude
  TCanvas* myCanvas03 = new TCanvas;

  for (Int_t i=0; i<2; i++) {
    hAmplitude[i]->GetXaxis()->SetTitle("#Gamma [MeV]");
    hAmplitude[i]->GetXaxis()->SetTitleSize(0.06);
    hAmplitude[i]->GetXaxis()->SetTitleOffset(0.85);
    hAmplitude[i]->GetXaxis()->SetLabelSize(0.05);
    hAmplitude[i]->GetYaxis()->SetTitle("Amplitude [nb]");
    hAmplitude[i]->GetYaxis()->SetTitleSize(0.06);
    hAmplitude[i]->GetYaxis()->SetTitleOffset(0.85);
    hAmplitude[i]->GetYaxis()->SetLabelSize(0.05);

    hAmplitude[i]->SetLineWidth(2);
    hAmplitude[i]->SetLineColor(1);
    hAmplitude[i]->SetMarkerColor(2);
    hAmplitude[i]->SetMarkerSize(1.0);
    hAmplitude[i]->SetMarkerStyle(24);

    hAmplitude[i]->Draw("");

    myCanvas03->Print(Form("output/plots/hAmplitude_%d.png",i),"png");

  }

  //UPPER LIMIT
  TCanvas* myCanvas04 = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {

    hXS_uppLimit[Bs]->SetTitle(Form("B_{s} = %d MeV", -Bs));
    hXS_uppLimit[Bs]->GetYaxis()->SetRangeUser(0,35);
    hXS_uppLimit[Bs]->GetXaxis()->SetTitle("#Gamma [MeV]");
    hXS_uppLimit[Bs]->GetXaxis()->SetTitleSize(0.06);
    hXS_uppLimit[Bs]->GetXaxis()->SetTitleOffset(0.95);
    hXS_uppLimit[Bs]->GetXaxis()->SetLabelSize(0.05);
    hXS_uppLimit[Bs]->GetYaxis()->SetTitle("#sigma_{upp}^{CL=90%} [nb]");
    hXS_uppLimit[Bs]->GetYaxis()->SetTitleSize(0.06);
    hXS_uppLimit[Bs]->GetYaxis()->SetTitleOffset(0.8);
    hXS_uppLimit[Bs]->GetYaxis()->SetLabelSize(0.05);

    hXS_uppLimit[Bs]->SetLineWidth(2);
    hXS_uppLimit[Bs]->SetLineColor(2);
    hXS_uppLimit[Bs]->SetFillColor(2);
    hXS_uppLimit[Bs]->SetFillStyle(3354);
    hXS_uppLimit_syst[Bs]->SetLineWidth(1);
    hXS_uppLimit_syst[Bs]->SetLineColor(1);
    hXS_uppLimit_syst[Bs]->SetFillColor(4);
    hXS_uppLimit_syst[Bs]->SetFillStyle(3344);
    //hXS_uppLimit[Bs]->SetMarkerColor(2);
    //hXS_uppLimit[Bs]->SetMarkerStyle(20);
    hXS_uppLimit[Bs]->DrawCopy("LF2");
    hXS_uppLimit_syst[Bs]->DrawCopy("same LF2");

    TPaveText *text04 = new TPaveText(15.,25.,15.,25.,"excldd");
    text04->SetTextFont(62); text04->SetTextSize(0.05);
    text04->SetTextAlign(12);
    text04->SetShadowColor(0); text04->SetFillColor(0);
    text04->SetBorderSize(0);
    text04->AddText("Excluded");
    text04->Draw();

    TLegend *myLegend04 = new TLegend(0.420, 0.755, 0.885, 0.885);
    myLegend04->SetFillStyle(0); myLegend04->SetFillColor(0); myLegend04->SetLineColor(0); myLegend04->SetTextSize(0.04);
    myLegend04->AddEntry(hXS_uppLimit[Bs], "upper limit", "l");
    myLegend04->AddEntry(hXS_uppLimit_syst[Bs], "systematic uncertainties", "f");
    myLegend04->Draw();

    myCanvas04->Print(Form("output/plots/upperLimit_dppi0_Bs%d.png",Bs),"png");
    //myCanvas04->Print(Form("output/plots/upperLimit_dppi0_Bs%d.eps",Bs),"eps");

  }

  TCanvas* myCanvas04_pl = new TCanvas;

  for(Int_t Bs=0; Bs<45; Bs=Bs+5) {

    hXS_uppLimit[Bs]->SetTitle(Form("B_{s} = %d MeV", -Bs));
    hXS_uppLimit[Bs]->GetXaxis()->SetTitle("#Gamma, MeV");
    hXS_uppLimit[Bs]->GetYaxis()->SetTitle("#sigma_{granica}^{CL=90%}, nb");
    hXS_uppLimit[Bs]->DrawCopy("LF2");
    hXS_uppLimit_syst[Bs]->DrawCopy("same LF2");

    TLegend *myLegend04_pl = new TLegend(0.420, 0.755, 0.885, 0.885);
    myLegend04_pl->SetFillStyle(1001); myLegend04_pl->SetFillColor(19); myLegend04_pl->SetLineColor(1); myLegend04_pl->SetTextSize(0.04); myLegend04_pl->SetBorderSize(5);
    myLegend04_pl->AddEntry(hXS_uppLimit[Bs], "\\hbox{górna granica}", "l");
    myLegend04_pl->AddEntry(hXS_uppLimit_syst[Bs], "\\hbox{niepewności systematyczne}", "f");
    myLegend04_pl->Draw();

    myCanvas04_pl->Print(Form("output/plots/pl/upperLimit_dppi0_Bs%d_pl.png",Bs),"png");
    myCanvas04_pl->Print(Form("output/plots/pl/upperLimit_dppi0_Bs%d_pl.eps",Bs),"eps");

  }

  gStyle->SetPadRightMargin(0.16);

  TCanvas* myCanvas05 = new TCanvas;

  hXS_uppLimit_2D->SetTitle("");
  hXS_uppLimit_2D->GetXaxis()->SetRangeUser(-50.,10.);
  hXS_uppLimit_2D->GetYaxis()->SetRangeUser(0.,55.);
  hXS_uppLimit_2D->GetXaxis()->SetTitle("B_{s} [MeV]");
  hXS_uppLimit_2D->GetXaxis()->SetTitleSize(0.06);
  hXS_uppLimit_2D->GetXaxis()->SetTitleOffset(0.95);
  hXS_uppLimit_2D->GetXaxis()->SetLabelSize(0.05);
  hXS_uppLimit_2D->GetYaxis()->SetTitle("#Gamma [MeV]");
  hXS_uppLimit_2D->GetYaxis()->SetTitleSize(0.06);
  hXS_uppLimit_2D->GetYaxis()->SetTitleOffset(0.8);
  hXS_uppLimit_2D->GetYaxis()->SetLabelSize(0.05);
  hXS_uppLimit_2D->GetZaxis()->SetTitle("#sigma_{upp}^{CL=90%} [nb]");
  hXS_uppLimit_2D->GetZaxis()->SetTitleOffset(0.8);
  hXS_uppLimit_2D->GetZaxis()->SetTitleSize(0.06);
  hXS_uppLimit_2D->GetZaxis()->SetLabelSize(0.05);

  hXS_uppLimit_2D->SetMinimum(11.5);
  //hXS_uppLimit_2D->SetMaximum(30.5);

  hXS_uppLimit_2D->Draw("colz");

  myCanvas05->Print("output/plots/upperLimit_dppi0_2D.png","png");
  //myCanvas04->Print("output/plots/upperLimit_dppi0_2D.eps","eps");

  //PL
  TCanvas* myCanvas05_pl = new TCanvas;

  hXS_uppLimit_2D->GetXaxis()->SetTitle("B_{s}, MeV");
  hXS_uppLimit_2D->GetYaxis()->SetTitle("#Gamma, MeV");
  hXS_uppLimit_2D->GetZaxis()->SetTitle("#sigma_{granica}^{CL=90%}, nb");
  hXS_uppLimit_2D->Draw("colz");

  myCanvas05_pl->Print("output/plots/pl/upperLimit_dppi0_2D_pl.png","png");
  myCanvas05_pl->Print("output/plots/pl/upperLimit_dppi0_2D_pl.eps","eps");

}
