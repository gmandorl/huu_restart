#include <stdio.h>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>
#include "TROOT.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TGaxis.h"
// #include <boost/format.hpp>
#include <array>
#include <algorithm>
#include <iterator>
#include <sstream>

#define SWAP(A, B) { Float_t t = A; A = B; B = t; }
#define SWAP2(A, B) { B.swap(A); }

void bubblesort(Float_t *a, std::string *str, int n)
{
  int i, j;
  for (i = n - 1; i >= 0; i--)
  {
    for (j = 0; j < i; j++)
    {
      if (a[j] < a[j + 1]){
        SWAP( a[j], a[j + 1] );
		SWAP2( str[j], str[j+1]);
		}
    }
  }
}


Float_t Discr(TH1F *h1, TH1F *h2){  
  Int_t h1s = h1->GetNbinsX();
  Int_t h2s = h2->GetNbinsX();
  if (h1s != h2s) {
  	printf("h1 size %d   h2 size %d \n",h1s,h2s);
  	return -1;
  }
  if (h1->Integral()!=0 ) h1->Scale(1./h1->Integral());
  if (h2->Integral()!=0 ) h2->Scale(1./h2->Integral());
  Float_t adiff = 0;
  for (Int_t i=0;i<h1s;i++){
	adiff +=  TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i));
	}
  return adiff/2;
}


Float_t add_underFlow_overFlow(TH1F *h){ 
    
    h->SetBinContent(1, h->GetBinContent(0) + h->GetBinContent(1));
    h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX()) + h->GetBinContent(h->GetNbinsX() + 1)); 
    
    h->SetBinError(1, sqrt(h->GetBinError(0)*h->GetBinError(0) + h->GetBinError(1)*h->GetBinError(1)));
    h->SetBinError(h->GetNbinsX(), sqrt(h->GetBinError(h->GetNbinsX())*h->GetBinError(h->GetNbinsX()) + h->GetBinError(h->GetNbinsX() + 1)*h->GetBinError(h->GetNbinsX() + 1))); 
    
}
    
    
    
    
    
void draw2Dhistos(TH2F * h2Dhistos_DY, TH2F * h2Dhistos_VBF, std::string figName) {
    
    h2Dhistos_DY->Scale(1./h2Dhistos_DY->Integral());
    h2Dhistos_VBF->Scale(1./h2Dhistos_VBF->Integral());
    
    h2Dhistos_DY->SetTitle("DY");
    h2Dhistos_VBF->SetTitle("H->mumu");
    
//     for(int ny = 1; ny < h2Dhistos_DY->GetNbinsY()+1;ny++) {
//         float integral_DY = 0;
//         float integral_VBF = 0;
//         for(int nx = 1; nx < h2Dhistos_DY->GetNbinsX()+1;nx++) { 
//             integral_DY += h2Dhistos_DY->GetBinContent(nx, ny);
//             integral_VBF += h2Dhistos_VBF->GetBinContent(nx, ny);
//         }
//         for(int nx = 1; nx < h2Dhistos_DY->GetNbinsX()+1;nx++) { 
//             h2Dhistos_DY->SetBinContent(nx, ny, h2Dhistos_DY->GetBinContent(nx, ny)/integral_DY);
//             h2Dhistos_VBF->SetBinContent(nx, ny, h2Dhistos_VBF->GetBinContent(nx, ny)/integral_DY);
//         }            
//             
//     }
    
         
    TCanvas * c = new TCanvas("c", "c", 1200, 600);
    c->Divide(2,1);
    gStyle->SetPaintTextFormat("2.2f");
    gPad->SetLogz();
    c->cd(1);
    gStyle->SetOptStat(0000);  
    gPad->SetLogz();
    gPad->SetRightMargin(.2);
    h2Dhistos_DY->Draw("colz text");
//     h2Dhistos_DY->Draw("CONT1Z");
    c->cd(2);
    gStyle->SetOptStat(0000);  
    gPad->SetLogz();
    gPad->SetRightMargin(.2);
//     h2Dhistos_VBF->Draw("CONT1Z");
    h2Dhistos_VBF->Draw("colz text");
    c->Print(("plotsDirectory/plot2D/"+figName+".png").c_str());
    
}
    
    
    
    
    
using namespace std;

int main(int argc, char* argv[]){

    
std::cout << "I can see here" << std::endl;

gROOT->ProcessLine(".x setTDRStyle.C");


int set_type=atoi(argv[1]); // 0 - analysis, 1 - control region , top
int era=atoi(argv[2]); // 2016, 2017 or 2018
int Zanalysis=atoi(argv[3]); // 2016, 2017 or 2018


// const int nfiles  = 10; //9;  //11;
// const int nfiles  = 11;

// const int nfiles  = 15;
const int nfiles  = 8;


// const int nfiles  = 14;



// TString leg_names[nfiles] = {"Data", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets", "EW lljj","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
// TString leg_names_2000[nfiles] = {"Data", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets", "EW lljj","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};



TString leg_names[nfiles] = {"Data", "Other bkg","Top", "EW lljj interf", "EW lljj", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
TString leg_names_2000[nfiles] = {"Data", "Other bkg","Top", "EW lljj interf", "EW lljj", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};






// TString leg_names[nfiles] = {"Data","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets", "EW lljj","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
// TString leg_names_2000[nfiles] = {"Data","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets", "EW lljj","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};



// TString leg_names[nfiles] = {"Data", "EW lljj", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
// TString leg_names_2000[nfiles] = {"Data", "EW lljj", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};


// TString leg_names[nfiles] = {"Data", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
// TString leg_names_2000[nfiles] = {"Data", "WW DPS", "WZjj", "WWJJToLNuLNu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "DYtautau", "DY + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};

// TString leg_names[nfiles] = {"Data","W(l#nu) + jets","WW + jets","Single Top", "t#bar{t}", "Z + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};
// TString leg_names_2000[nfiles] = {"Data","W(l#nu) + jets","WW + jets","Single Top", "t#bar{t}", "Z + jets","VBF H #rightarrow ll x20","GluGlu H #rightarrow ll x20"};


//---------------------------------------------------------------------------------------------------------------------- OK
/// change EWK to VBF bla
//TString leg_names[nfiles] = {"Data","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "Z + jets (AMC)","VBF Z #rightarrow ll"};
TString set_names[2] = {"Dimuon","Dielectron"}; 
//if (set_type==0)leg_names[0] = "Data (DoubleB)";
//if (set_type==1)leg_names[0] = "Data (SingleB)";

// if (Zanalysis) set_names[0]="Z Control Region";

/// change EWK to VBF bla-----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//TString file_names[nfiles] =        {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};

TString DY_name = "DYInclusivetoLL";
TString DY_file_name = "DYInclusivetoLL_M";

DY_name = "DYJetstoLL_amc";
DY_file_name = "DYJetstoLL_amc";

// DY_name = "DYJetstoLL_amc_Filter105";
// DY_file_name = "DYJetstoLL_amc_Filter105";

// DY_name = "DYJetsToLL_M-105To160-madgraphMLM";
// DY_file_name = "DYJetsToLL_M-105To160-madgraphMLM";



//DY_name = "DYJetstoLL_amc_Pt";
//DY_file_name = "DYJetstoLL_amc_Pt";

//DY_name = "DYJetstoLL_Final_amc";
//DY_file_name = "DYJetstoLL_Final_amc";


// TString file_names[nfiles] =        {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","VBF_HToMuMu", "GluGlu_HToMuMu"};

// TString file_names[nfiles] =        {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};


// TString file_names[nfiles] =        {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};


// TString file_names[nfiles] =        {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "EWK_LLJJ", "VBF_HToMuMu", "GluGlu_HToMuMu"};


// TString file_names[nfiles] =        {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};




// TString file_names[nfiles] =        {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon", "EWK_LLJJ", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};



TString file_names[nfiles] =         {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_QCDup[nfiles] =   {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_QCDdown[nfiles] = {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_JESup[nfiles] =   {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_JESdown[nfiles] = {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};





// ------- Z ordered sample ----------

// TString file_names_Z[nfiles] =        {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "EWK_LLJJ", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup_Z[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "EWK_LLJJ", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown_Z[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "EWK_LLJJ", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup_Z[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "EWK_LLJJ", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown_Z[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu_EWK_noTop","WJetsToLNu","WW","ZZ","WZ","ST","TT", "EWK_LLJJ", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};




TString file_names_Z[nfiles] =         {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_QCDup_Z[nfiles] =   {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_QCDdown_Z[nfiles] = {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_JESup_Z[nfiles] =   {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
TString file_names_JESdown_Z[nfiles] = {"SingleMuon", "Other","Top","EWK_LLJJ_INT","EWK_LLJJ",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};










// TString file_names[nfiles] =        {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon", "WWTo2L2Nu_DoubleScattering", "WLLJJ_WToLNu_EWK", "WWJJToLNuLNu","WJetsToLNu","WW","ZZ","WZ","ST","TT", "DYJetstoTauTau_amc",DY_file_name, "VBF_HToMuMu", "GluGlu_HToMuMu"};



// TString file_names[nfiles] =        {"SingleMuon","WJetsToLNu","WW","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_QCDdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
// TString file_names_JESdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};




int bg_begin;
int qcd_begin=18;
 bg_begin=1;

 


// int FILLCOLOR[nfiles] = {1,kMagenta+3, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int LINECOLOR[nfiles] = {1,kMagenta+3, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,2};
// int LINEWIDTH[nfiles] = {1,1,1,1,1,1,1,1,1,1,1,1,1,3,3};
// int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001};

int FILLCOLOR[nfiles] = {1, kGreen+1,kTeal+10,kMagenta,kMagenta+3,kOrange,kRed-4,kMagenta+1};
int LINECOLOR[nfiles] = {1, kGreen+1,kTeal+10,kMagenta,kMagenta+3,kOrange,kRed-4,kMagenta+1};
int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,2};
int LINEWIDTH[nfiles] = {1,1,1,1,1,1,3,3};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001};

// int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,1,1,1,1,2};
// int LINEWIDTH[nfiles] = {1,1,1,1,1,1,1,1,1,1,3,3};
// int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001};


// int FILLCOLOR[nfiles] = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int LINECOLOR[nfiles] = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3,kAzure+1,kBlue-4,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,1,1,1,1,1,1,2};
// int LINEWIDTH[nfiles] = {1,1,1,1,1,1,1,1,1,1,1,1,3,3};
// int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001};


	
//int order[nfiles] = {0,1,2,3,4,5,6,7,8,9,10};// number in file_names array as they should appear in the legend
// int order[nfiles] = {0,1,2,3,4,5,6,7,8,9};// number in file_names array as they should appear in the legend
// int order[nfiles] = {0,1,2,3,4,5,6,7,8,9,10,11};// number in file_names array as they should appear in the legend
// int order[nfiles] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};// number in file_names array as they should appear in the legend
// int order[nfiles] = {0,1,2,3,4,5,6,7};// number in file_names array as they should appear in the legen

// int order[nfiles] = {0,12,1,2,3,4,5,6,7,8,9,10,11,13,14};// number in file_names array as they should appear in the legend
int order[nfiles] = {0,1,2,3,4,5,6,7};// number in file_names array as they should appear in the legend



// int FILLCOLOR_Z[nfiles] = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kMagenta+3,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int LINECOLOR_Z[nfiles] = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kMagenta+3,kOrange+2,kOrange-2,kRed-4,kMagenta+1};
// int order_Z[nfiles] = {0,1,2,3,4,5,6,7,8,9,12,10,11,13,14};// number in file_names array as they should appear in the legend
// 
// if (Zanalysis) {
//         for (int n=0; n<nfiles; n++) {
//             FILLCOLOR[n] = FILLCOLOR_Z[n];
//             LINECOLOR[n] = LINECOLOR_Z[n];
//             order[n] = order_Z[n];
//             
//             file_names[n]           = file_names_Z[n];
//             file_names_QCDup[n]     = file_names_QCDup_Z[n];
//             file_names_QCDdown[n]   = file_names_QCDdown_Z[n];
//             file_names_JESup[n]     = file_names_JESup_Z[n];
//             file_names_JESdown[n]   = file_names_JESdown_Z[n];
//         }
// }
    
    // if (Zanalysis) {
//     FILLCOLOR = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange+2,kOrange-2,kMagenta+3,kRed-4,kMagenta+1};
//     LINECOLOR = {1, kPink, kPink-4, kPink-9,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange+2,kOrange-2,kMagenta+3,kRed-4,kMagenta+1};
//     order = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};// number in file_names array as they should appear in the legend
// }



int order_legend[nfiles]; 
for (int i=0;i<nfiles;i++){
	order_legend[order[i]]=i;
}


TString data_name[2] = {"SingleMuon","SingleElectron"};

TString set[3]={"_mu","_el"}; 

for (int i=0;i<nfiles;i++){
	if (i==0) file_names[i] = data_name[set_type];
        file_names[i].Prepend("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/histoFileDir/"); 
//         file_names[i].Prepend("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/HistoFileDir/"); 

	file_names[i].Append(set[set_type]);
	file_names_QCDup[i] =   file_names[i];
	file_names_QCDdown[i] = file_names[i];
	file_names_JESup[i] =   file_names[i];
	file_names_JESdown[i] = file_names[i];


//         if (i==0) file_names[i].Append("_QCDScalenom_JESnom_v25_reskim");
//     	if (i!=0) {
// 		file_names[i].Append("_QCDScalenom_JESnom_v25_reskim");
// 		file_names_QCDup[i].Append("_QCDScaleup_JESnom_v25_reskim");
// 		file_names_QCDdown[i].Append("_QCDScaledown_JESnom_v25_reskim");
// 		file_names_JESup[i].Append("_QCDScalenom_JESup_v25_reskim");
// 		file_names_JESdown[i].Append("_QCDScalenom_JESdown_v25_reskim");
// 	}

	
//         if (i==0) file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_v25_reskim");
//     	if (i!=0) {
// 		file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_v25_reskim");
// 		file_names_QCDup[i].Append("_QCDScaleup_JESnom_JERnom_PUnom_v25_reskim");
// 		file_names_QCDdown[i].Append("_QCDScaledown_JESnom_JERnom_PUnom_v25_reskim");
// 		file_names_JESup[i].Append("_QCDScalenom_JESup_JERnom_PUnom_v25_reskim");
// 		file_names_JESdown[i].Append("_QCDScalenom_JESdown_JERnom_PUnom_v25_reskim");
// 	}
	
        if (era == 2017) {
            if (i==0) file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2017");
            if (i!=0) {
                    file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2017");
                    file_names_QCDup[i].Append("_QCDScaleup_JESnom_JERnom_PUnom_nano_2017");
                    file_names_QCDdown[i].Append("_QCDScaledown_JESnom_JERnom_PUnom_nano_2017");
                    file_names_JESup[i].Append("_QCDScalenom_JESup_JERnom_PUnom_nano_2017");
                    file_names_JESdown[i].Append("_QCDScalenom_JESdown_JERnom_PUnom_nano_2017");
            }
        }
	

	if (era == 2016) {
            if (i==0) file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2016");
            if (i!=0) {
                    file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2016");
                    file_names_QCDup[i].Append("_QCDScaleup_JESnom_JERnom_PUnom_nano_2016");
                    file_names_QCDdown[i].Append("_QCDScaledown_JESnom_JERnom_PUnom_nano_2016");
                    file_names_JESup[i].Append("_QCDScalenom_JESup_JERnom_PUnom_nano_2016");
                    file_names_JESdown[i].Append("_QCDScalenom_JESdown_JERnom_PUnom_nano_2016");
            }
        }

        
	if (era == 2018) {
            if (i==0) file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2018");
            if (i!=0) {
                    file_names[i].Append("_QCDScalenom_JESnom_JERnom_PUnom_nano_2018");
                    file_names_QCDup[i].Append("_QCDScaleup_JESnom_JERnom_PUnom_nano_2018");
                    file_names_QCDdown[i].Append("_QCDScaledown_JESnom_JERnom_PUnom_nano_2018");
                    file_names_JESup[i].Append("_QCDScalenom_JESup_JERnom_PUnom_nano_2018");
                    file_names_JESdown[i].Append("_QCDScalenom_JESdown_JERnom_PUnom_nano_2018");
            }
        }
	
	
	
	
	
//	if (i==0) file_names[i].Append("_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//	if (i!=0) {
//		file_names[i].Append("_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//		file_names_QCDup[i].Append("_QCDScaleup_JESnom_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//		file_names_QCDdown[i].Append("_QCDScaledown_JESnom_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//		file_names_JESup[i].Append("_QCDScalenom_JESup_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//		file_names_JESdown[i].Append("_QCDScalenom_JESdown_v25_bdt_alldata4_qglweightsnorm_reminiaod");
//	}


	file_names[i].Append(".root");
	file_names_QCDup[i].Append(".root");
	file_names_QCDdown[i].Append(".root");
	file_names_JESup[i].Append(".root");
	file_names_JESdown[i].Append(".root");
}

//TString dir_name= "plots_mdg_ht_bdt_alldata_v25_reaod";
//TString dir_name= "plots_amc_bdt_alldata_v25_reaod";
//TString dir_name= "plots_amc_bdt_axis2_v25_reaod";
//TString dir_name= "plots_amc_bdt_alldata4_qglnormChrisEl_v25_reaod";
//---------------------------------------------------------------------------------------------------------------------- OK
//create directory for output plots and change the name here
//TString dir_name= "plots_amc_bdt_alldata4_qglnorm_v25_reaod";
TString dir_name= "plotsDirectory";
//TString dir_name= "plots_mdg_ht_bdt_axis2_v25_reaod";
//TString dir_name= "plots_amc_bdt_v25";
//TString dir_name= "plots_amc_herwig_bdt_v25";
// dir_name.Append(set[set_type]+"/");
dir_name.Append("/");
//TString dir_name = "plots_amc/";
Float_t lumi = 35900*0.93;
// Float_t lumi = 41000;
if (era==2018) lumi = 59970;
if (era==2017) lumi = 41530;



TLegend *leg = new TLegend(0.77,0.55,0.92,0.9); //without writing about SF
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetTextSize(0.025);
TLegend *leg_2000 = new TLegend(0.77,0.65,0.92,0.9); //with  writing about SF
leg_2000 ->SetFillColor(0);
leg_2000 ->SetBorderSize(0);
leg_2000 ->SetTextFont(42);
leg_2000 ->SetTextSize(0.025);

//const int nhistos = 100; //40//52
//TString hist_names[nhistos]={"hMqq", "hEtaQQ","hHTsoft","hSoft_n2","hSoft_n5","hSoft_n10","hHTsoftEWK","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hHTsoftEWK_bdt","hSoft_n2EWK_bdt", "hSoft_n5EWK_bdt","hSoft_n10EWK_bdt","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult", "hmet",   "hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt","hV_mass", "hqgl", "hqgl2", "hZll_mass", "hZll_pt", "hZll_phi", "hZll_eta",  "hrho", "hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", "hDeltaRelQQ", "hRptHard", "hEtaQQSum", "hPhiZQ1", "hZll_y", "hZll_ystar", "hZll_zstar", "hMqq_log","hJet3_pt","hPhiQQ","hJets12_pt","hJets12_pt_log","hJet1q_pt_log","hJet2q_pt_log","hbdt","hbdt_atanh","hJet3_pt_bdt", "hAdJetHT_bdt","hHT","hlheHT_log","hNAdJets", "hNAdJets_bdt","hNAdJets_bdt2", "hJet3_pt_bdt2", "hAdJetHT_bdt2","hNAdJets_mjj1", "hJet3_pt_mjj1", "hAdJetHT_mjj1","hNAdJets_mjj2", "hJet3_pt_mjj2", "hAdJetHT_mjj2", "hHTsoftEWK_bdt2","hSoft_n2EWK_bdt2","hSoft_n5EWK_bdt2","hSoft_n10EWK_bdt2", "hHTsoftEWK_mjj1","hSoft_n2EWK_mjj1","hSoft_n5EWK_mjj1","hSoft_n10EWK_mjj1","hHTsoftEWK_mjj2","hSoft_n2EWK_mjj2","hSoft_n5EWK_mjj2","hSoft_n10EWK_mjj2", "hJet1q_eta_bdt","hJet1q_eta_bdt2","hJet2q_eta_bdt","hJet2q_eta_bdt2","hbdt_atanh2","hsoftleadTrackPt", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta" , "hJet3_pt_new","hJet3_eta_bdt","hJet3_eta_bdt2"};


/*
const int nhistos =105; //79 //40//52
TString hist_names[nhistos]={"hMqq", "hZll_mass", "hSelectionCuts","hdeltaM", "hdeltaMRel", "hJet1q_pt","hEtaQQ","hHTsoftEWK","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hPVs", "hJet1q_eta","hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2","hJet2q_mult", "hVtype", "hVtypeSim", "hmet","hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt", "hqgl", "hqgl2", "hZll_pt", "hZll_phi", "hZll_eta","hrho","hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", "hDeltaRelQQ", "hRptHard", "hEtaQQSum", "hPhiZQ1", "hZll_y", "hZll_ystar", "hZll_zstar", "hMqq_log", "hlheV_pt","hJet3_pt","hPhiQQ", "hJet1q_pt_log", "hJet2q_pt_log","hJets12_pt","hJets12_pt_log", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta", "hJet3_pt_new", "hMaxJetBTagCSV","hVirtual_Pt1","hVirtual_Pt2", "hMaxJetBTagCMVA","hthetastar_W1","hthetastar_W2", "hMaxSecondJetBTagCMVA", "hMaxSecondJetBTagCSV","hVirtual_eta1","hVirtual_eta2","hNAdJets","hgen_mass", "hlheNj","hlheNpNLO", "hHT","hlheHT_log","hweights_weighted","hweights", "hThetaStarJet","hThetaPlanes","hThetaStar","hDiffmass","hThetaStarAbs","hTheta_HiggsJ1","hTheta_HiggsJ2","hthetastar_W2toHW1","hthetastar_W1toHW2","hthetastar_HtoWW", "hTotalEnergy","hTotalEnergylog","hVirtual_phi1","hParton_M1","hParton_M2","hVirtual_phi2","hWWmass", "hEnergy_fraction_Parton1","hPz","hPzAbs" ,"hVirtual_Wmass1","hVirtual_Wmass2", "hInvariant_Masslog","hInvariant_Mass" , "hEnergy_fraction_Parton2","hBDT_VBF","hBDT_VBF_atanh"};
*/



// const int nhistos =96 ; //79 //40//52
// TString hist_names[nhistos]={"hMqq", "hZll_mass", "hEtaQQ","hHTsoftEWK","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult", "hmet",   "hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt", "hqgl", "hqgl2", "hZll_pt", "hZll_phi", "hZll_eta",  "hrho", "hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", "hDeltaRelQQ", "hRptHard", "hEtaQQSum", "hPhiZQ1", "hZll_y", "hZll_ystar", "hZll_zstar", "hMqq_log", "hlheV_pt","hPhiQQ", "hJets12_pt","hJets12_pt_log", "hJet1q_pt_log", "hJet2q_pt_log", "hHT","hlheHT_log", "hlheNj", "hNAdJets", "hJet3_pt", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta", "hJet3_pt_new", "hMaxJetBTagCSV","hVirtual_Pt1","hVirtual_Pt2", "hMaxJetBTagCMVA","hthetastar_W1","hthetastar_W2", "hMaxSecondJetBTagCMVA", "hMaxSecondJetBTagCSV","hVirtual_eta1","hVirtual_eta2", "hThetaStarJet", "hThetaPlanes", "hThetaStar","hDiffmass", "hThetaStarAbs","hTheta_HiggsJ1","hTheta_HiggsJ2","hthetastar_W2toHW1","hthetastar_W1toHW2","hthetastar_HtoWW", "hTotalEnergy","hTotalEnergylog","hVirtual_phi1","hParton_M1","hParton_M2","hVirtual_phi2","hWWmass", "hEnergy_fraction_Parton1","hPz","hPzAbs" ,"hVirtual_Wmass1","hVirtual_Wmass2", "hInvariant_Masslog","hInvariant_Mass" , "hEnergy_fraction_Parton2", "hBDT_VBF" , "hBDT_VBF_atanh"};



// const int nhistos =118 ; //79 //40//52
// TString hist_names[nhistos]={"hMqq", "hZll_mass", "hEtaQQ","hHTsoftEWK","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult", "hmet",   "hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt", "hqgl", "hqgl2", "hqglAtanh", "hqgl2Atanh", "hZll_pt", "hZll_phi", "hZll_eta",  "hrho", "hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", "hDeltaRelQQ", "hRpt", "hRptAtanh", "hEtaQQSum", "hPhiZQ1", "hZll_y", "hZll_ystar", "hZll_zstar", "hMqq_log", "hlheV_pt","hPhiQQ", "hJets12_pt","hJets12_pt_log", "hJet1q_pt_log", "hJet2q_pt_log", "hHT","hlheHT_log", "hlheNj", "hNAdJets", "hJet3_pt", "hJet3_pt_log", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta", "hJet3_pt_new", "hMaxJetBTagCSV", "hBDiscriminator_CSV", "hBDiscriminator_CMVA","hVirtual_Pt1","hVirtual_Pt2","hVirtual_Pt1_log", "hVirtual_Pt2_log", "hgen_mass", "hgenJetMass",  "hMaxJetBTagCMVA","hthetastar_W1","hthetastar_W2", "hMaxSecondJetBTagCMVA", "hMaxSecondJetBTagCSV","hVirtual_eta1","hVirtual_eta2", "hThetaStarJet", "hThetaPlanes", "hThetaStar", "hThetaStarJetAtanh", "hThetaPlanesAtanh","hDiffmass", "hThetaStarAbs","hTheta_HiggsJ1","hTheta_HiggsJ2","hthetastar_W2toHW1","hthetastar_W1toHW2","hthetastar_HtoWW", "hmumujj_pt", "hmumujj_pt", "hmumujj_ptLog", "hEnergy_fraction_Parton1_log", "hdeltaR1", "hdeltaR2",  "hEnergy_fraction_Parton2_log","hdeltaM", "hdeltaMRel", "hTotalEnergy","hTotalEnergylog","hVirtual_phi1","hVirtual_phi2","hWWmass", "hEnergy_fraction_Parton1","hPz","hPzAbs", "hPzAbsLog","hVirtual_Wmass1","hVirtual_Wmass2", "hVirtual_Wmass1_log", "hVirtual_Wmass2_log",  "hInvariant_MassLog","hInvariant_Mass" , "hEnergy_fraction_Parton2", /*"hNNplusBDT_atanh",*/ /*"hNNoutput" , "hNNoutput_atanh",*/ "hBDT_VBF" , "hBDT_VBF_atanh"};


// ULTIMA VERSIONE 
// const int nhistos =140 ; //79 //40//52
// TString hist_names[nhistos]={"hMqq", "hZll_mass", "hHll_mass_precise", "hHll_mass_unprecise", "hEtaQQ","hHTsoftEWK","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hPVs", "hJet1q_pt", "hJet1q_eta", "hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult", "hmet",   "hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt", "hqgl", "hqgl2", "hqglAtanh", "hqgl2Atanh", "hZll_pt", "hZll_phi", "hZll_eta",  "hrho", "hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", "hDeltaRelQQ", "hRpt", "hRptAtanh", "hEtaQQSum", "hPhiZQ1", "hZll_y", "hZll_ystar", "hZll_zstar", "hMqq_log", "hlheV_pt","hPhiQQ", "hJets12_pt","hJets12_pt_log", "hJet1q_pt_log", "hJet2q_pt_log", "hHT","hlheHT_log", "hlheNj", "hNAdJets", "hJet3_pt", "hJet3_pt_log", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta", "hJet3_pt_new", "hMaxJetBTagCSV", "hBDiscriminator_CSV", "hBDiscriminator_CMVA","hVirtual_Pt1","hVirtual_Pt2","hVirtual_Pt1_log", "hVirtual_Pt2_log", "hgen_mass", "hgenJetMass",  "hMaxJetBTagCMVA","hthetastar_W1","hthetastar_W2", "hMaxSecondJetBTagCMVA", "hMaxSecondJetBTagCSV","hVirtual_eta1","hVirtual_eta2", "hThetaStarJet", "hThetaPlanes", "hThetaStar", "hThetaStarJetAtanh", "hThetaPlanesAtanh","hDiffmass", "hThetaStarAbs","hTheta_HiggsJ1","hTheta_HiggsJ2","hthetastar_W2toHW1","hthetastar_W1toHW2","hthetastar_HtoWW", "hmumujj_pt", "hmumujj_pt", "hmumujj_ptLog", "hEnergy_fraction_Parton1_log", "hdeltaR1", "hdeltaR2",  "hEnergy_fraction_Parton2_log","hdeltaM", "hdeltaMRel", "hTotalEnergy","hTotalEnergylog","hVirtual_phi1","hVirtual_phi2","hWWmass", "hEnergy_fraction_Parton1","hPz","hPzAbs", "hPzAbsLog","hVirtual_Wmass1","hVirtual_Wmass2", "hVirtual_Wmass1_log", "hVirtual_Wmass2_log",  "hInvariant_MassLog","hInvariant_Mass" , "hEnergy_fraction_Parton2", "hgenJetMass_bothMatching", "hgenJetMass_minus1", "hgenJetMass_oneIsNotMatching", "hgenJetMass_noMatching", "hmatchLeading", "hmatchSubleading", "hminAbsEta", "hminAbsEta_cut06", "hminAbsEta_cut08", "hminAbsEta_GEN", "hgenJetMass_cut06", "hgenJetMass_cut08", "hgenJetMass_cut10", "hgenJetMass_cut12", "hMqq_cut06", "hMqq_cut08", "hMqq_cut10", "hMqq_cut12", "hminAbsEta_cut10", "hminAbsEta_cut12", /*"hNNplusBDT_atanh", "hNNoutput" , "hNNoutput_atanh",*/ "hBDT_VBF" , "hBDT_VBF_atanh"};





const int nhistos =149; //79 //40//52
TString hist_names[nhistos]={/*"hBDT_VBF_atanh", */"hMqq", "hSelectionCuts", "hMqq_cut06", "hMqq_cut08", /*"hMqq_cut10", "hMqq_cut12",*/ "hZll_mass", "hZll_mass_complete"/*, "hHll_mass_precise", "hHll_mass_unprecise"*/ , "hEtaQQ","hHTsoftEWK", "hPtSoftJets","hSoft_n2EWK","hSoft_n5EWK","hSoft_n10EWK","hPVs", "hJet1q_neHEF", "hJet1q_chHEF", "hJet1q_neEmEF", "hJet1q_chEmEF", "hJet2q_neHEF", "hJet2q_chHEF", "hJet2q_neEmEF", "hJet2q_chEmEF" , "hJetUnder50_neHEF", "hJetUnder50_chHEF", "hJetUnder50_neEmEF", "hJetUnder50_chEmEF", "hJetOver50_neHEF", "hJetOver50_chHEF", "hJetOver50_neEmEF", "hJetOver50_chEmEF", "hJetUnder50_eta", "hJetOver50_eta", "hJet1q_pt", "hJet1q_eta", "hJet1q_phi", "hJet2q_phi", "hJet1q_ptd", "hJet1q_axis2", "hJet1q_mult", "hJet2q_pt", "hJet2q_eta", "hJet2q_ptd", "hJet2q_axis2", "hJet2q_mult", "hmet",   "hJet1q_leadTrackPt", "hJet2q_leadTrackPt", "hqq_pt", "hqgl", "hqgl2", "hqglAtanh", "hqgl2Atanh", "hZll_pt", "hZll_phi", "hZll_eta",  "hrho", "hlepton1_pt", "hlepton2_pt", "hlepton1_eta", "hlepton2_eta","hlepton1_iso03", "hlepton2_iso03", /*"hGenLepton_matching1", "hGenLepton_matching2",*/ "hnGenLep",  "hDeltaRelQQ", "hRpt", "hRptAtanh", "hEtaQQSum", "hPhiZQ1", "hPhiZQ2", "hEtaHQ1", "hEtaHQ2", "hZll_y", "hZll_ystar", "hZll_zstar", "hZll_zstar_log", "hMqq_log", "hlheV_pt","hPhiQQ", "hJets12_pt","hJets12_pt_log", "hJet1q_pt_log", "hJet2q_pt_log", "hHT","hlheHT_log"/*, "hlheHT_log_BDTgt08","hlheHT", "hlheHT_BDTgt08"*/, "hlheNj", "hNAdJets", "hJet3_pt", "hJet3_pt_log", "hsoftleadTrackEta", "hAdJetHT", "hJet3_eta", "hJet3_pt_new", "hMaxJetBTagCSV", "hBDiscriminator_CSV", "hBDiscriminator_CMVA","hVirtual_Pt1","hVirtual_Pt2","hVirtual_Pt1_log", "hVirtual_Pt2_log", "hgen_mass", /*"hgenJetMass",*/  "hMaxJetBTagCMVA","hthetastar_W1","hthetastar_W2", "hMaxSecondJetBTagCMVA", "hMaxSecondJetBTagCSV","hVirtual_eta1","hVirtual_eta2", "hThetaStarJet", "hThetaPlanes", "hThetaStar", "hThetaStarJetAtanh", "hThetaPlanesAtanh","hDiffmass",/* "hIsolatedElectrons", "hpuId_jetWithoutGenJet", "hpuId_AllJets",*/  "hThetaStarAbs","hTheta_HiggsJ1","hTheta_HiggsJ2","hthetastar_W2toHW1","hthetastar_W1toHW2","hthetastar_HtoWW", /*"hmaxAbsEta_BDTmFixgt08", "hmaxAbsEta_BDTmFixgt1", "hmaxAbsEta",*/ "hminAbsEta", "hmumujj_pt", "hmumujj_pt", "hmumujj_ptLog", "hEnergy_fraction_Parton1_log", "hdeltaR1", "hdeltaR2",  "hEnergy_fraction_Parton2_log","hdeltaM", "hdeltaMRel", "hTotalEnergy","hTotalEnergylog","hVirtual_phi1","hVirtual_phi2","hWWmass", "hEnergy_fraction_Parton1","hPz","hPzAbs", "hPzAbsLog","hVirtual_Wmass1","hVirtual_Wmass2", "hVirtual_Wmass1_log", "hVirtual_Wmass2_log",  "hInvariant_MassLog","hInvariant_Mass" , "hEnergy_fraction_Parton2",  /*"hNNplusBDT_atanh", "hNNoutput_atanh_mFix", "hNNoutput_atanh_mFix_zFix","hNNoutput_atanh_mFix_controlRegion", "hNNoutput_atanh_mFix_controlRegionDown", "hNNoutput_atanh_mFix_controlRegionUp",  /*"hNNoutput" , "hNNoutput_atanh","hBDT_VBF_atanh_m125_Rpt0", "hBDT_VBF_atanh_m125_WM1Log5", "hBDT_VBF_atanh_m125_zStar0", "hBDT_VBF_atanh_m125_ptll300", "hBDT_VBF_atanh_m125_softN50", "hBDT_VBF_atanh_m125_MqqLog75", "hBDT_VBF_atanh_m125_mumujjPt",  "hBDT_VBF_atanh_m125_DEtajj", "hBDT_VBF_atanh_m125_Theta2", "hBDT_VBF_atanh_m125_onlyMqq", "hBDT_VBF_atanh_m125_onlyRpt", "hBDT_VBF_atanh_m125_onlyWM1", "hBDT_VBF_atanh_m125_onlyzStar", "hBDT_VBF_atanh_m125_onlyptll", "hBDT_VBF_atanh_m125_onlysoftN5", "hBDT_VBF_atanh_m125_onlymumujjPt", "hBDT_VBF_atanh_m125_onlyDEtajj", "hBDT_VBF_atanh_m125_onlyTheta2", "hBDT_VBF_atanh_m125ForAll",*/ "hBDT_VBF_atanh_m125ControlRegion", "hBDT_VBF_atanh_m125ControlRegionUp", "hBDT_VBF_atanh_m125ControlRegionDown", "hBDT_VBF" , "hBDT_VBF_atanh"};





// const int nhistos =23; //79 //40//52
// TString hist_names[nhistos]={"hMqq", "hmatchLeading", "hmatchSubleading", "hminAbsEta", "hminAbsEta_cut06", "hminAbsEta_cut08", "hminAbsEta_GEN", "hZll_mass", "hdeltaR1", "hdeltaR2", "hgenJetMass_cut06", "hgenJetMass_cut08", "hgenJetMass_cut10", "hgenJetMass", "hgenJetMass_cut12", "hMqq_cut06", "hMqq_cut08", "hMqq_cut10", "hMqq_cut12", "hminAbsEta_cut10", "hminAbsEta_cut12",  "hBDT_VBF" , "hBDT_VBF_atanh"};

std::vector<std::string> variablesName_in_2D_plot =  {"singleMuTrigger", "doubleMuTrigger"}; // {"Zll_mass", "deltaM"};
// std::vector<std::string> variablesName_in_2D_plot = {"Zll_mass", "deltaM", "Xparton1Log", "Xparton2Log", "RpT", "zStar"};
// std::vector<std::string> variablesName_in_2D_plot = {"Zll_mass", "deltaM", "deltaR1", "deltaR2", "RpT", "zStar", "BDToutput", "genJetMassLeading", "hgenJetMassMatched", "Mqq", "indexFirstJet", "indexSecondJet"};
std::vector<std::string> hist2D_names;
std::vector<TH2F*> h2Dhistos_DY;
std::vector<TH2F*> h2Dhistos_VBF;
    

// hist2D_names.push_back("histo2D_Zll_mass_deltaM");    
// hist2D_names.push_back("histo2D_Xparton1Log_Xparton2Log");        

for (int i = 0; i < variablesName_in_2D_plot.size(); i++) {
    for (int j = i+1; j < variablesName_in_2D_plot.size(); j++) {

        hist2D_names.push_back(("histo2D_"+variablesName_in_2D_plot[i]+"_"+variablesName_in_2D_plot[j]).c_str());        
        
    }
}
   
   
// hist2D_names.push_back("hNN_Vs_BDT_atanh");     
        

std::vector<TH1F*> hisoSensitivity_binByBin_VECTOR;


std::array<int,200> LOGY_array = {};


// ,"hpdgId"

TString hist_names_sum[nhistos]={};
TString sum_histos_names[nhistos]={};
std::string hist_names_sort[nhistos];


for (int i=0;i<nhistos;i++){
	hist_names_sort[i] = hist_names[i];
	hist_names_sum[i] = hist_names[i];
	hist_names_sum[i].Append("_sum");
	sum_histos_names[i] = hist_names[i];
	sum_histos_names[i].Append("_sum0");
}


TString stacks_names[nhistos];
for (int i=0;i<nhistos;i++){
	stacks_names[i] = hist_names[i];
	stacks_names[i].Prepend("s");
}


TString output_names[nhistos];
for (int i=0;i<nhistos;i++){
	output_names[i] = hist_names[i];
	output_names[i].Prepend(dir_name);
	output_names[i].Append("_v25.png");
}





TH1F *data_histos[nhistos];
TH1F *data_histos2[nhistos];
TH1F *data_histosQCDup[nhistos];
TH1F *data_histosQCDlo[nhistos];
TH1F *data_histosJESup[nhistos];
TH1F *data_histosJESlo[nhistos];
TH1F *data_histosTrig[nhistos];
TH1F *signal_histos[nhistos];
TH1F *signal_histos2[nhistos];
TH1F *gluonFu_histos[nhistos];
TH1F *tthbb_histos[nhistos];
TH1F *tthnbb_histos[nhistos];
TH1F *gf_histos[nhistos];
TH1F *sum_histos[nhistos];
TH1F *sum_histosUp[nhistos];
TH1F *sum_histosDown[nhistos];
TH1F *histos_forUnc[nhistos];
TH1F *histos_for_legened[nhistos];
TH1F *discr_histos[nhistos];//Mqq,delta eta, delta phi, qgl, btag //12,13,14,21,22
TH1F *hBkgVis[nhistos];
TH1F *hBkgUncUp[nhistos];
TH1F *hBkgUncLo[nhistos];
TH1F *hBkgUncUpTrig[nhistos];
TH1F *hBkgUncLoTrig[nhistos];
TH1F *hBkgQCDUp[nhistos];
TH1F *hBkgQCDLo[nhistos];
TH1F *hBkgJESUp[nhistos];
TH1F *hBkgJESLo[nhistos];



int files=0; 
THStack *stacks[nhistos];
TH1F * stackHisto[nhistos];
for (int i=0;i<nhistos;++i){
	stacks[i] = new THStack(stacks_names[i],"");
}

Double_t totalBG=0.;
Double_t totalQCD=0.;
Double_t totalMC=0.;
Double_t totalData=0.;
Double_t totalDataQCD=0.;
ofstream out_efficiency;
ofstream out_discrimination;
out_efficiency.open(dir_name+"efficiency.txt"); 

//Float_t MC_data[2] = {1., 1.};//ht  without qcd

//Float_t MC_data[2] = {0.871,0.915};//{0.78,1.};//////////amc  (after muon corrections)
//Float_t MC_dataup[2] = {0.92,0.97};//{0.832,1.};//////////amc
//Float_t MC_datalow[2] = {0.78,0.808};//{0.701,1.};//////////amc
//Float_t MC_datajesup[2] = {0.871,0.855};//{0.733,1.};//////////amc
//Float_t MC_datajeslow[2] = {0.954,1}; //{0.861,1.};//////////amc
//Float_t MC_data[2] = {1.04,1.};////////////mdg lo
//Float_t MC_dataup[2] = {1.16,1.};//////////mdg up
//Float_t MC_datalow[2] = {0.838,1.};//////////mdg lo
//Float_t MC_data[2] =  {1.14, 1.19}; // {0.893, 1.19};//ht
//Float_t MC_dataup[2] = {1.11, 1.16}; //{0.872,1.16};//////////ht up
//Float_t MC_datalow[2] =  {1.07,1.13}; // {0.844,1.13};//////////ht lo
//Float_t MC_datajesup[2] = {1.06,1.12}; //{0.839,1.12};//////////ht lo  
//Float_t MC_datajeslow[2] = {1.24,1.3}; //{0.976,1.3};/////////ht lo 
//Float_t MC_data[2] = {1.14, 1.19};//ht  without qcd
//Float_t MC_dataup[2] = {1.11,1.16};//////////ht up   wo wcd
//Float_t MC_datalow[2] = {1.07,1.13};//////////ht lo   wo qcd
//Float_t MC_datajesup[2] = {1.14,1.19};//////////ht lo   wo qcd
//Float_t MC_datajeslow[2] = {1.24,1.3};/////////ht lo wo qcd
//Float_t MC_data[2] = {1.13, 1.19};//ht  Pt sum
//Float_t MC_dataup[2] = {1.11,1.16};//////////ht  PT sum
//Float_t MC_datalow[2] = {1.07,1.13};//////////ht PT sum
//Float_t MC_datajesup[2] = {1.06,1.12};//////////ht PT sum
//Float_t MC_datajeslow[2] = {1.24,1.3};/////////ht PT sum 
//Float_t MC_data[2] = {1.25,1.19};//ht  Pt sumEtaQQ
//Float_t MC_dataup[2] = {1.19,1.17};//////////ht  PT sumEtaQQ
//Float_t MC_datalow[2] = {1.18,1.13};//////////ht PT sumEtaQQ
//Float_t MC_datajesup[2] = {1.25,1.11};//////////ht PT sumEtaQQ
//Float_t MC_datajeslow[2] = {1.35,1.31};/////////ht PT sumEtaQQ
//
//Float_t MC_data[2] = {1.14, 1.19};//ht  without qcd
//Float_t MC_dataup[2] = {1.11,1.16};//////////ht up   wo wcd
//Float_t MC_datalow[2] = {1.07,1.12};//////////ht lo   wo qcd
//Float_t MC_datajesup[2] = {1.13,1.18};//////////ht lo   wo qcd
//Float_t MC_datajeslow[2] = {1.24,1.29};/////////ht lo wo qcd





// DY bg ratio is calculated here
float ratio[5];
float ratioError[5];
TString file_names_mc[nfiles];
for (int counter_ratio=0; counter_ratio < 5;counter_ratio++) {

    float dataInt = 0;
    float MCint = 0;
    float DYint =0;
    float MCsign =0;
    
    double dataError = 0;
    double MCError = 0;
    double DYError =0;
    
    
    //Set file_names_mc
    file_names_mc[0] = file_names[0];
    for (int i=1;i<nfiles;i++) {
        if (counter_ratio==0) file_names_mc[i] = file_names[i];
        if (counter_ratio==1) file_names_mc[i] = file_names_QCDup[i];
        if (counter_ratio==2) file_names_mc[i] = file_names_QCDdown[i];
        if (counter_ratio==3) file_names_mc[i] = file_names_JESup[i];
        if (counter_ratio==4) file_names_mc[i] = file_names_JESdown[i];
        if (((counter_ratio==1) || (counter_ratio==2))&& (i==5)) {
            file_names_mc[i]=file_names[i];
            cout << " No QCD up e down \n " << file_names[i] << endl;
        }
    }

    files=0;
    for (int fileIterator = 0; fileIterator < nfiles; fileIterator++) {

        TFile *file_initial_mc;
        std::cout << file_names[fileIterator] << std::endl;

        file_initial_mc = TFile::Open(file_names_mc[fileIterator]);
        string file_name_tag = file_names_mc[fileIterator].Data();

        TH1F *histos_mc[100];
        int hist=0;

        histos_mc[hist] = (TH1F*)file_initial_mc->Get(hist_names[hist])->Clone("mc");
           
        if (fileIterator==0) dataInt = histos_mc[hist]->IntegralAndError(0,histos_mc[hist]->GetNbinsX()+1, dataError); 
        if (fileIterator!=0) {
            histos_mc[hist]->Scale(lumi);
            
/*            if (file_name_tag.find("HToMuMu")==std::string::npos) MCint+=histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1); 
            else MCsign +=histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1);*/ 
            
//             std::cout << "INTEGRALI    " << file_name_tag << "  \t " << histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1) << std::endl;

            if (file_name_tag.find(DY_name)!=std::string::npos) {
                DYint= histos_mc[hist]->IntegralAndError(0,histos_mc[hist]->GetNbinsX()+1, DYError);
                
            }
            else MCint+=histos_mc[hist]->IntegralAndError(0,histos_mc[hist]->GetNbinsX()+1, MCError); 

            
        }

    }

    float tmp_ratio = (dataInt-MCint)/DYint;
    ratio[counter_ratio] = tmp_ratio;
    ratioError[counter_ratio] = DYError/DYint*ratio[counter_ratio];
//     ratio[counter_ratio] = 1.;
    if(counter_ratio == 0)  out_efficiency << "ratio  "  << counter_ratio << " : " << ratio[0] <<   " \t\t DYint : " << DYint <<   " \t\t MCint : " << MCint  <<   " \t\t dataInt : " << dataInt << endl;

}




TH1F *histos_check[nhistos];
TH1F *histos_JESup_check[nhistos];




files=0;
for (int fileIterator = 0; fileIterator < nfiles; fileIterator++) {

	TFile *file_initial;
        std::cout << file_names[fileIterator] << std::endl;
  	file_initial = TFile::Open(file_names[fileIterator]);
	string file_name_tag = file_names[fileIterator].Data();
	string leg_name_tag = leg_names[fileIterator].Data();

        TH1F *histos[nhistos];

	for (int hist=0;hist<nhistos;++hist){
        std::cout << hist << " \t "  << hist_names[hist] << std::endl;
	//	if (hist_names[hist].CompareTo("hlepton1_eta")==0) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(4);
	//	if (hist_names[hist].CompareTo("hlepton2_eta")==0) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt")==0) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt_atanh")==0) ((TH1F*)file_initial->Get(hist_names[hist]))->Rebin(10);
		histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("h");
                if(file_name_tag.find("DYJetstoLL")!=std::string::npos)   histos_check[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("h");
                    
                double bkgIntegralErrorForBUG = 0;
//                 if (hist<3) out_efficiency<< "search for bug 0 " << leg_names[order[fileIterator]]<<"\t\t\t"<< std::setprecision(8)<<histos[hist]->IntegralAndError(1,histos[hist]->GetNbinsX(), bkgIntegralErrorForBUG) << " +- " << bkgIntegralErrorForBUG  <<endl;


                if (fileIterator==0) {
                    data_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("data");
                    if (!((hist_names[hist].CompareTo("hJetOver50_eta")==0)|| (hist_names[hist].CompareTo("hJetUnder50_eta")==0))) add_underFlow_overFlow(data_histos[hist]);
                    
                }
//                 if (hist<3) out_efficiency<< "search for bug 1 " << leg_names[order[fileIterator]]<<"\t\t\t"<< std::setprecision(8)<<histos[hist]->IntegralAndError(1,histos[hist]->GetNbinsX(), bkgIntegralErrorForBUG) << " +- " << bkgIntegralErrorForBUG  <<endl;
		if (fileIterator>0)histos[hist]->Scale(lumi); 
//                 if (hist<3) out_efficiency<< "search for bug 2 " << leg_names[order[fileIterator]]<<"\t\t\t"<< std::setprecision(8)<<histos[hist]->IntegralAndError(1,histos[hist]->GetNbinsX(), bkgIntegralErrorForBUG) << " +- " << bkgIntegralErrorForBUG  <<endl;
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos[hist]->Scale(ratio[0]);  
//                 if ((file_name_tag.find("HToMuMu")==std::string::npos) && (file_name_tag.find("SingleMuon")==std::string::npos))  histos[hist]->Scale(ratio[0]);  
                
                


                
                
                if (!((hist_names[hist].CompareTo("hJetOver50_eta")==0) || (hist_names[hist].CompareTo("hJetUnder50_eta")==0))) add_underFlow_overFlow(histos[hist]);

//                 histos[hist]->SetBinContent(1, histos[hist]->GetBinContent(0) + histos[hist]->GetBinContent(1));
//                 histos[hist]->SetBinContent(histos[hist]->GetNbinsX(), histos[hist]->GetBinContent(histos[hist]->GetNbinsX()) + histos[hist]->GetBinContent(histos[hist]->GetNbinsX() + 1));
    
    
// 		if ((hist_names[hist].CompareTo("hNAdJets_mjj1")==0) && ((file_name_tag.find("DYJetstoLL")!=std::string::npos)|| (file_name_tag.find("EWK_LLJJ")!=std::string::npos)) ) 
// 			cout<< " mjj > 1500 ,Integral " <<leg_name_tag<<"   "<< histos[hist]->Integral()<<endl; 
// 		if ((hist_names[hist].CompareTo("hNAdJets_mjj2")==0) && ((file_name_tag.find("DYJetstoLL")!=std::string::npos)|| (file_name_tag.find("EWK_LLJJ")!=std::string::npos)) ) 
// 			cout<< " mjj > 2500 ,Integral " <<leg_name_tag<<"   "<< histos[hist]->Integral()<<endl; 
		if ((hist_names[hist].CompareTo("hNAdJets_bdt")==0) && ((file_name_tag.find(DY_name)!=std::string::npos)|| (file_name_tag.find("EWK_LLJJ")!=std::string::npos)) ) 
			cout<< " bdt > 0.92 ,Integral " <<leg_name_tag<<"   "<< histos[hist]->Integral()<<endl; 
		if ((hist_names[hist].CompareTo("hNAdJets_bdt2")==0) && ((file_name_tag.find(DY_name)!=std::string::npos)|| (file_name_tag.find("EWK_LLJJ")!=std::string::npos)) ) 
			cout<< " bdt > 0.84 ,Integral " <<leg_name_tag<<"   "<< histos[hist]->Integral()<<endl; 
////////////////////////////
////////////////////////////
////////////////////////////
// Signal now is included in the sum total MC


		if (fileIterator==1) 	sum_histos[hist] = (TH1F*)histos[hist]->Clone(sum_histos_names[hist]);
		if (fileIterator>1)	sum_histos[hist]->Add(histos[hist]); 


////////////////////////////
////////////////////////////
////////////////////////////
		if (fileIterator>0) histos[hist]->Sumw2(kFALSE);    //  COMMENT THIS FOR ERRORS
// 		if (fileIterator>0) histos[hist]->Sumw2();
//		if (hist==1) cout<<fileIterator<<"   "<<histos[1]->Integral() <<endl;


                float signalMultipyFactor = 20.;
//                 if ((hist == nhistos-3) || (hist == nhistos-3))  signalMultipyFactor = 20.;
	//	if (fileIterator==1) {
		if (fileIterator==nfiles-2) {
			signal_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"newhist");
			signal_histos[hist]->Scale(lumi);

			signal_histos[hist]->Scale(signalMultipyFactor);
			signal_histos[hist]->Sumw2(kFALSE);
// 			signal_histos[hist]->Sumw2();
			signal_histos[hist]->SetLineColor(LINECOLOR[fileIterator]);
			signal_histos[hist]->SetLineWidth(LINEWIDTH[fileIterator]);
// 			signal_histos[hist]->SetLineColor(kRed+2);
			signal_histos[hist]->SetLineStyle(LINESTYLE[fileIterator]);
			signal_histos2[hist]=(TH1F*)signal_histos[hist]->Clone("signalHist2");
                        
                        add_underFlow_overFlow(signal_histos[hist]);
                        add_underFlow_overFlow(signal_histos2[hist]);
//            for (int n=0; n < signal_histos[hist]->GetNbinsX()+2; ++n) {
//	            signal_histos[hist]->SetBinError(n, signal_histos[hist]->GetBinError(n));
//	            signal_histos2[hist]->SetBinError(n, signal_histos2[hist]->GetBinError(n));
//            }
		}
  

                if (fileIterator==nfiles-1) {
                        gluonFu_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"newhist");
                        gluonFu_histos[hist]->Scale(lumi);

			gluonFu_histos[hist]->Scale(signalMultipyFactor);
			gluonFu_histos[hist]->Sumw2(kFALSE);
// 			gluonFu_histos[hist]->Sumw2();
			gluonFu_histos[hist]->SetLineColor(LINECOLOR[fileIterator]);
                        gluonFu_histos[hist]->SetLineWidth(LINEWIDTH[fileIterator]);
                        gluonFu_histos[hist]->SetLineStyle(2);

// 			gluonFu_histos[hist]->SetLineColor(kMagenta);
			gluonFu_histos[hist]->SetLineStyle(LINESTYLE[fileIterator]);
//			gluonFu_histos[hist]=(TH1F*)signal_histos[hist]->Clone("signalHist2");
//            for (int n=0; n < gluonFu_histos[hist]->GetNbinsX()+2; ++n)
//	            gluonFu_histos[hist]->SetBinError(n, gluonFu_histos[hist]->GetBinError(n));
                        add_underFlow_overFlow(gluonFu_histos[hist]);
                        
		}

//        for (int n=0; n < histos[hist]->GetNbinsX()+2; ++n)
//		    histos[hist]->SetBinError(n, histos[hist]->GetBinError(n)/2.);
		histos[hist]->SetLineColor(LINECOLOR[fileIterator]);
		histos[hist]->SetLineStyle(LINESTYLE[fileIterator]);
		histos[hist]->SetLineWidth(LINEWIDTH[fileIterator]);
		histos[hist]->SetFillStyle(FILLSTYLE[fileIterator]);
		if ((fileIterator!=0)) histos[hist]->SetFillColor(FILLCOLOR[fileIterator]);

                
		if (fileIterator==0) {
			histos[hist]->SetMarkerStyle(20);
			data_histos[hist]->SetLineColor(1);
			data_histos[hist]->SetMarkerStyle(10);
			data_histos[hist]->SetMarkerSize(.8);
		}
	 	//if (files>=bg_begin) stacks[hist]->Add(histos[hist]);
	 	if (fileIterator>=1) {
                    stacks[hist]->Add(histos[hist]);
                    if (fileIterator==1)stackHisto[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("stackHisto");
                    else stackHisto[hist]->Add(histos[hist]);
                }
		if (hist==0) histos_for_legened[fileIterator] = (TH1F*)histos[0]->Clone("newd");
		if (fileIterator==bg_begin)	discr_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("discr");
		if ((fileIterator>bg_begin)&&(fileIterator!=(nfiles-2))){
			discr_histos[hist]->Add(histos[hist]); 
                }

                

            }

            
            std::cout << "Reading 2D histos" << std::endl;
//             for(int n=0; n<hist2D_names.size(); n++) {
//                 if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  {
//                     TH2F* h =  (TH2F*)file_initial->Get(hist2D_names[n].c_str())->Clone(("h2D"+hist2D_names[n]).c_str());
//                     h->Scale(ratio[0]);
//                     h2Dhistos_DY.push_back(h); 
//                     
//                 }
//                 if (file_name_tag.find("VBF_HToMuMu")!=std::string::npos)  {
//                     TH2F* h =  (TH2F*)file_initial->Get(hist2D_names[n].c_str())->Clone(("h2D"+hist2D_names[n]).c_str());
//                     h2Dhistos_VBF.push_back(h); 
//                 }
//             }
            std::cout << "finished reading 2D histos" << std::endl;
            
            double bkgIntegralError = 0;
            if (fileIterator>=bg_begin) totalBG+=histos[0]->Integral(1,histos[0]->GetNbinsX());   //from 1 to histos[0]->GetNbinsX() because I added overflow and underflow
            if (fileIterator>0) totalMC+=histos[0]->Integral(1,histos[0]->GetNbinsX());
            if (fileIterator==0) {totalData+=histos[0]->Integral(1,histos[0]->GetNbinsX()); }
            if (fileIterator==0) out_efficiency<<"Sample  \t\t\t yield(per "<< lumi<<" pb^-1)"<<endl;
//             if (fileIterator==0) out_efficiency<<leg_names[order[fileIterator]]<<"\t \t \t"<< std::setprecision(8)<<histos[0]->Integral(1,histos[0]->GetNbinsX()) <<endl;
//             else out_efficiency<<leg_names[order[fileIterator]]<<"\t\t\t  "<<std::setprecision(8)<<histos[0]->Integral(1,histos[0]->GetNbinsX())<<endl;
            float bkgIntegral = histos[0]->IntegralAndError(1,histos[0]->GetNbinsX(), bkgIntegralError);
            out_efficiency<<leg_names[order[fileIterator]]<<"\t\t\t"<< std::setprecision(8)<<histos[0]->IntegralAndError(1,histos[0]->GetNbinsX(), bkgIntegralError) << " +- " << bkgIntegralError  <<endl;

            if (fileIterator==nfiles-1) out_efficiency<<"Total BG"<<"\t \t \t  "<<std::setprecision(8)<<totalBG<<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Total MC"<<"\t \t \t  "<<std::setprecision(8)<<totalMC<<endl;
//            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t \t \t  "<<std::setprecision(3)<<totalData/totalMC<<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t \t \t \t  "<<std::setprecision(8)<<ratio[0]<< " +- " << ratioError[0] <<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t QCDup   \t \t  "<<std::setprecision(8)<<ratio[1]<< " +- " << ratioError[1] <<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t QCDdown \t \t  "<<std::setprecision(8)<<ratio[2]<< " +- " << ratioError[2] <<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t JESup   \t \t  "<<std::setprecision(8)<<ratio[3]<< " +- " << ratioError[3] <<endl;
            if (fileIterator==nfiles-1) out_efficiency<<"Data/MC"<<"\t JESdown \t \t  "<<std::setprecision(8)<<ratio[4]<< " +- " << ratioError[4] <<endl;
//            if (fileIterator==nfiles-1) out_efficiency<<"DY MC scaled with "<<"\t \t \t  "<<std::setprecision(13)<<ratio[0]<<endl;
                if(file_name_tag.find("DYJetstoLL")!=std::string::npos  ) cout << "ratio4  " << fileIterator << "  "  << ratio[0]  << " +- " << ratioError[0] <<   " \t\t integral : " << histos[0]->Integral(0,histos[0]->GetNbinsX()+1) << endl;


}



out_efficiency.close();



std::cout << "drawing 2D histos" << std::endl;

// for(int n=0; n<hist2D_names.size(); n++) {
// 
// draw2Dhistos(h2Dhistos_DY[n], h2Dhistos_VBF[n], hist2D_names[n]);
// }
std::cout << "finished drawing 2D histos" << std::endl;


/////////////////////////////////////


//----------------------------------------------------------------------------------------------------------------------
///JES and QCD part ->no need

files=1;
for (int fileIterator = 1; fileIterator < nfiles; fileIterator++) {

	TFile *file_initial_up;		
	string file_name_tag = file_names_QCDup[fileIterator].Data();
// // // // // //   	if (file_name_tag.find("ST_")!=std::string::npos) {  // this lines are useless
// // // // // // 		file_names_QCDup[fileIterator] = file_names[fileIterator];
// // // // // // 		string file_name_tag = file_names_QCDup[fileIterator].Data();
// // // // // // 	}
   file_initial_up = TFile::Open(file_names_QCDup[fileIterator]);
//	cout<<file_names_QCDup[fileIterator] <<endl;
	TH1F *histos_QCDup[nhistos];
        
        std::cout << "1 \t " <<  file_names[fileIterator] << std::endl;
        
        
	for (int hist=0;hist<nhistos;++hist){

// 		if (hist_names[hist].CompareTo("hbdt")==0) ((TH1F*)file_initial_up->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt_atanh")==0) ((TH1F*)file_initial_up->Get(hist_names[hist]))->Rebin(10);
		histos_QCDup[hist] = (TH1F*)file_initial_up->Get(hist_names[hist])->Clone("hup");
		histos_QCDup[hist]->Scale(lumi);  //top
//		if (hist==0) cout<<histos_QCDup[hist]->Integral() <<endl;
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos_QCDup[hist]->Scale(ratio[1]);  
//                 if (file_name_tag.find("HToMuMu")==std::string::npos)  histos_QCDup[hist]->Scale(ratio[1]);  

// 		if (file_name_tag.find("EWK_LLJJ")!=std::string::npos)	histos_QCDup[hist]->Scale(0.97957); //new signal x-sec
		if (fileIterator==1) 	hBkgQCDUp[hist] = (TH1F*)histos_QCDup[hist]->Clone(sum_histos_names[hist]+"QCDup");
		if (fileIterator>1)	hBkgQCDUp[hist]->Add(histos_QCDup[hist]);
		hBkgQCDUp[hist]->SetLineColor(kRed);
		hBkgQCDUp[hist]->SetLineStyle(2);
	}
	TFile *file_initial_down;
	file_name_tag = file_names_QCDdown[fileIterator].Data();
  	if (file_name_tag.find("ST_")!=std::string::npos) {
		file_names_QCDdown[fileIterator] = file_names[fileIterator];
		string file_name_tag = file_names_QCDdown[fileIterator].Data();
	}
	std::cout << "2 \t " <<  file_names[fileIterator] << std::endl;
//	cout<<file_names_QCDdown[fileIterator] <<endl;
  	file_initial_down = TFile::Open(file_names_QCDdown[fileIterator]);
	TH1F *histos_QCDdown[nhistos];
	for (int hist=0;hist<nhistos;++hist){
// 		if (hist_names[hist].CompareTo("hbdt")==0) ((TH1F*)file_initial_down->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt_atanh")==0) ((TH1F*)file_initial_down->Get(hist_names[hist]))->Rebin(10);
		histos_QCDdown[hist] = (TH1F*)file_initial_down->Get(hist_names[hist])->Clone("hdown");
		histos_QCDdown[hist]->Scale(lumi);  //top
	//	if (hist==0) cout<<histos_QCDdown[hist]->Integral() <<endl;
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos_QCDdown[hist]->Scale(ratio[2]);  
//                 if (file_name_tag.find("HToMuMu")==std::string::npos)  histos_QCDdown[hist]->Scale(ratio[2]);  

// 		if (file_name_tag.find("EWK_LLJJ")!=std::string::npos)	histos_QCDdown[hist]->Scale(0.97957); //new signal x-sec
		if (fileIterator==1) 	hBkgQCDLo[hist] = (TH1F*)histos_QCDdown[hist]->Clone(sum_histos_names[hist]+"QCDdown");
		if (fileIterator>1)	hBkgQCDLo[hist]->Add(histos_QCDdown[hist]);
		hBkgQCDLo[hist]->SetLineColor(kBlue);
		hBkgQCDLo[hist]->SetLineStyle(3);
	}

}

files=1;
for (int fileIterator = 1; fileIterator < nfiles; fileIterator++) {
    
	TFile *file_initial_up;		
  	file_initial_up = TFile::Open(file_names_JESup[fileIterator]);
	string file_name_tag = file_names_JESup[fileIterator].Data();
        cout << "file_names_JESup[fileIterator] " << file_names_JESup[fileIterator] << endl;
	TH1F *histos_JESup[nhistos];
	for (int hist=0;hist<nhistos;++hist){
// 		if (hist_names[hist].CompareTo("hbdt")==0) ((TH1F*)file_initial_up->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt_atanh")==0) ((TH1F*)file_initial_up->Get(hist_names[hist]))->Rebin(10);
		histos_JESup[hist] = (TH1F*)file_initial_up->Get(hist_names[hist])->Clone("hup");
                if(file_name_tag.find("DYJetstoLL")!=std::string::npos)
                    histos_JESup_check[hist] = (TH1F*)file_initial_up->Get(hist_names[hist])->Clone("data_JESup_check");
		histos_JESup[hist]->Scale(lumi);  //top
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos_JESup[hist]->Scale(ratio[3]);  
//                 if (file_name_tag.find("HToMuMu")!=std::string::npos)  histos_JESup[hist]->Scale(ratio[3]);
                
// 		if (file_name_tag.find("EWK_LLJJ")==std::string::npos)	histos_JESup[hist]->Scale(0.97957); //new signal x-sec
		if (fileIterator==1) 	hBkgJESUp[hist] = (TH1F*)histos_JESup[hist]->Clone(sum_histos_names[hist]+"JESup");
		if (fileIterator>1)	hBkgJESUp[hist]->Add(histos_JESup[hist]);
		hBkgJESUp[hist]->SetLineColor(kRed);
		hBkgJESUp[hist]->SetLineStyle(2);
                

                    
                    
	}
	TFile *file_initial_down;
  	file_initial_down = TFile::Open(file_names_JESdown[fileIterator]);
	file_name_tag = file_names_JESdown[fileIterator].Data();
	TH1F *histos_JESdown[nhistos];
	for (int hist=0;hist<nhistos;++hist){
// 		if (hist_names[hist].CompareTo("hbdt")==0) ((TH1F*)file_initial_down->Get(hist_names[hist]))->Rebin(4);
// 		if (hist_names[hist].CompareTo("hbdt_atanh")==0) ((TH1F*)file_initial_down->Get(hist_names[hist]))->Rebin(10);
		histos_JESdown[hist] = (TH1F*)file_initial_down->Get(hist_names[hist])->Clone("hdown");
		histos_JESdown[hist]->Scale(lumi);  //top
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos_JESdown[hist]->Scale(ratio[4]); 
//                 if (file_name_tag.find("HToMuMu")==std::string::npos)  histos_JESdown[hist]->Scale(ratio[4]);
                
// 		if (file_name_tag.find("EWK_LLJJ")!=std::string::npos)	histos_JESdown[hist]->Scale(0.97957); //new signal x-sec
		if (fileIterator==1) 	hBkgJESLo[hist] = (TH1F*)histos_JESdown[hist]->Clone(sum_histos_names[hist]+"JESdown");
		if (fileIterator>1)	hBkgJESLo[hist]->Add(histos_JESdown[hist]);
		hBkgJESLo[hist]->SetLineColor(kBlue);
		hBkgJESLo[hist]->SetLineStyle(3);
	}

}



// //////////////TROVA BACO///////////////////////////////////
// 
// TCanvas * c_check = new TCanvas("c_check", "", 800, 800);
// c_check->cd();
// // histos_JESup_check[2]->Draw();
// // histos_check[2]->Draw("same");
// // histos_check[2]->SetLineColor(2);
// // c_check->Print("check.png");
// // 
// // for(int i = 0 ; i < histos_JESup_check[2]->GetNbinsX(); i++ )
// // histos_JESup_check[2]->SetBinContent(i, histos_check[2]->GetBinContent(i)>0 ? histos_JESup_check[2]->GetBinContent(i)/histos_check[2]->GetBinContent(i) - 1. : 0.38);
// // 
// // histos_JESup_check[2]->GetXaxis()->SetRangeUser(0,200);
// // histos_JESup_check[2]->GetYaxis()->SetRangeUser(-0.4,0.4);
// // histos_JESup_check[2]->Draw();
// // c_check->Print("check_ratio.png");
// // 
// 
//  
// hBkgJESUp[2]->Draw();
// sum_histos[2]->Draw("same");
// sum_histos[2]->SetLineColor(2);
// c_check->Print("check.png");
// 
// for(int i = 0 ; i < hBkgJESUp[2]->GetNbinsX(); i++ )
// hBkgJESUp[2]->SetBinContent(i, sum_histos[2]->GetBinContent(i)>0 ? hBkgJESUp[2]->GetBinContent(i)/sum_histos[2]->GetBinContent(i) - 1. : 0.38);
// 
// hBkgJESUp[2]->GetXaxis()->SetRangeUser(0,200);
// hBkgJESUp[2]->GetYaxis()->SetRangeUser(-0.4,0.4);
// hBkgJESUp[2]->Draw();
// c_check->Print("check_ratio.png");
// 
// 
// 
// 
// //////////////END TROVA BACO///////////////////////////////////


/////////////////////////////////////---------------------------------------------------------------------------------------------


for (int hist=0;hist<nhistos;hist++){
	hBkgUncUp[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncUp");
	hBkgUncLo[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncLo");
  	hBkgVis[hist]   = (TH1F*)sum_histos[hist]->Clone("hbkgVis");
  	for(int i=0;i<hBkgUncUp[hist]->GetNbinsX();i++) {
  		float e = 0.0;
    	if (sum_histos[hist]->GetBinContent(i+1) != 0) {
      	e = sum_histos[hist]->GetBinError(i+1)/sum_histos[hist]->GetBinContent(i+1);
   	}
    	hBkgUncUp[hist]->SetBinContent(i+1,e);
    	hBkgUncLo[hist]->SetBinContent(i+1,-e);
	}
  	hBkgVis[hist]->SetMarkerSize(0);
 	hBkgVis[hist]->SetFillColor(kBlack);
 	hBkgVis[hist]->SetFillStyle(3004);
 	hBkgUncUp[hist]->SetLineColor(kBlack);
 	hBkgUncUp[hist]->SetLineWidth(1);
 	hBkgUncUp[hist]->SetFillColor(kBlack);
 	hBkgUncUp[hist]->SetFillStyle(3004);
	hBkgUncLo[hist]->SetLineColor(kBlack);
 	hBkgUncLo[hist]->SetLineWidth(1);
  	hBkgUncLo[hist]->SetFillColor(kBlack);
	hBkgUncLo[hist]->SetFillStyle(3004);
}

Float_t TSF[2] = {1.,1.};
Float_t kfactors[2] = {0.,0.};
//kfactors[set_type] = 1./(MC_data[set_type] * TSF[set_type]);
kfactors[set_type] = 1./(ratio[0] * TSF[set_type]);
//cout<<"kfactor = "<<kfactors[set_type]<<endl;
//cout<<"Data/MC = "<< MC_data[set_type]<<endl;


for (int i=0;i<nfiles-2;i++){
	if (i==0) { 
//             leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"P");
//             leg_2000->AddEntry(histos_for_legened[order_legend[i]],leg_names_2000[i],"P");
            
            leg->AddEntry(histos_for_legened[i],leg_names[i],"P");
            leg_2000->AddEntry(histos_for_legened[i],leg_names_2000[i],"P");
        }
	if (i>=bg_begin) {
            leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"F");
            leg_2000->AddEntry(histos_for_legened[order_legend[i]],leg_names_2000[i],"F");
        }
}
	leg->AddEntry(histos_for_legened[nfiles-2],leg_names[nfiles-2],"L");
	leg->AddEntry(histos_for_legened[nfiles-1],leg_names[nfiles-1],"L");
        leg_2000->AddEntry(histos_for_legened[nfiles-2],leg_names_2000[nfiles-2],"L");
	leg_2000->AddEntry(histos_for_legened[nfiles-1],leg_names_2000[nfiles-1],"L");
        
leg->AddEntry(hBkgUncUp[0],"MC stat. unc.","F");
leg_2000->AddEntry(hBkgUncUp[0],"MC stat. unc.","F");




// COMPUTE d VALUES FOR ALL THE HISTOGRAMS
Float_t discriminators[nhistos];
for (int d=0;d<nhistos;d++){
	discriminators[d] = Discr(discr_histos[d],signal_histos2[d]);
}

bubblesort(discriminators, hist_names_sort,nhistos);

//out_discrimination.open("Aftertriggercorr2/"+trigger[set_type]+dir_name+"discrimination.txt");
out_discrimination.open(dir_name+"discrimination.txt");
for (int d=0;d<nhistos;d++){
	if (d==0) out_discrimination<<"Variable &\t d"<<endl;
	out_discrimination<<"$"<<hist_names_sort[d]<<"$"<<" & \t "<< std::setprecision(2)<< discriminators[d]<<endl;
}
out_discrimination.close();

TString lumiString = " fb^{-1} (13 TeV)";
if (era==2016) lumiString = "35.9" + lumiString;
if (era==2017) lumiString = "41.53" + lumiString;
if (era==2018) lumiString = "59.97" + lumiString;
TLatex* tex = new TLatex(0.75,0.95,lumiString);
tex->SetNDC();
tex->SetTextAlign(35);
tex->SetTextFont(42);
tex->SetTextSize(0.035);
tex->SetLineWidth(2);
TLatex *tex1 = new TLatex(0.17,0.95,"CMS");
tex1->SetNDC();
tex1->SetTextAlign(20);
tex1->SetTextFont(61);
tex1->SetTextSize(0.04);
tex1->SetLineWidth(2);
TLatex* tex2 = new TLatex(0.27,0.89,"Work in progress");
tex2->SetNDC();
tex2->SetTextAlign(20);
tex2->SetTextFont(52);
tex2->SetTextSize(0.035);
tex2->SetLineWidth(2);
TString temp_str;
temp_str.Form("%2.2f",kfactors[set_type]);
//		temp_str.Form("%2.2f",(1./MC_data[set_type]));
TString k_factor_str;
TLatex* tex_set = new TLatex(0.667,0.86,set_names[set_type]);
tex_set->SetNDC();
tex_set->SetTextAlign(20);
tex_set->SetTextFont(42);
tex_set->SetTextSize(0.03);
tex_set->SetLineWidth(2);
TLatex* tex_k = new TLatex(0.63,0.89,k_factor_str);
tex_k->SetNDC();
tex_k->SetTextAlign(20);
tex_k->SetTextFont(42);
tex_k->SetTextSize(0.03);
tex_k->SetLineWidth(2);
	float left2 = gStyle->GetPadLeftMargin();
	float right2 = gStyle->GetPadRightMargin();
	float top2 = gStyle->GetPadTopMargin();
	float bottom2 = gStyle->GetPadBottomMargin();
	TPaveText pCMSset(0.57,1.-top2*2.,0.67,0.92,"NDC");
	pCMSset.SetTextFont(42);
	pCMSset.SetTextSize(top2*0.75);
	pCMSset.SetTextAlign(12);
	pCMSset.SetFillStyle(-1);
	pCMSset.SetBorderSize(0);
	if (!Zanalysis) pCMSset.AddText(set_names[set_type]);
	else pCMSset.AddText("Z region");

	

        
        
        
        
//------------------------------------- DRAW THE HISTOGRAMS -----------------------------------------------------------
        

for (int i=0;i<nhistos;i++){
//	for (int i=0;i<4;i++){
		//temp_str.Form("%2.2f",Discr(discr_histos[i],signal_histos2[i]));
		Float_t d_value = Discr((TH1F*)data_histos[i]->Clone("data_clone_"+hist_names_sum[i]),signal_histos2[i]);
		temp_str.Form("%2.2f",d_value);
		TString disc_value = temp_str.Prepend(" d = ");
		TLatex *disc_value_text = new TLatex(0.62,0.83,disc_value);
      disc_value_text->SetNDC();
      disc_value_text->SetTextAlign(20);
      disc_value_text->SetTextFont(42);
      disc_value_text->SetTextSize(0.03);
      disc_value_text->SetLineWidth(2);
		
		temp_str.Form("%2d",i);
		TString can_name="c1";
		can_name.Append(temp_str);
		TCanvas *c1 = new TCanvas(can_name,"",900,750);
		c1->cd();
		gPad->SetLogy(0);
		c1->SetBottomMargin(.3);
		c1->SetRightMargin(.25);
	
	//	bool LOGY=false;
	///	if (hist_names[i].CompareTo("hMqq")==0) LOGY=true; 
	//	if (hist_names[i].CompareTo("hqq_pt")==0) LOGY=true; 
		bool LOGY=true;

		Double_t xmin = signal_histos[i]->GetBinCenter(0);
		Double_t xmax = signal_histos[i]->GetBinCenter(signal_histos[i]->GetNbinsX())+signal_histos[i]->GetBinWidth(signal_histos[i]->GetNbinsX());
		if (hist_names[i].CompareTo("hPVs")==0) {
// 			xmax=30;
			LOGY=false;
		}
                if (hist_names[i].CompareTo("hZll_mass")==0) {
                    LOGY=false;
		}
		if ((hist_names[i].CompareTo("hThetaStar")==0) || (hist_names[i].CompareTo("hThetaStarAbs")==0)) {
			LOGY=false;
		}
		string tmp_hist_name = hist_names[i].Data();
		if ((hist_names[i].CompareTo("hAdJetHT_bdt")==0) || (hist_names[i].CompareTo("hNAdJetHT_bdt")==0) ||  (hist_names[i].CompareTo("hJet3_pt_bdt")==0)|| (hist_names[i].CompareTo("hJet3_eta_bdt")==0) ||(tmp_hist_name.find("EWK_bdt")!=std::string::npos) ||(tmp_hist_name.find("_bdt2")!=std::string::npos ) ||(tmp_hist_name.find("_mjj1")!=std::string::npos) ||(tmp_hist_name.find("_mjj2")!=std::string::npos ) ) {
			LOGY=false;
		}
		
		if ((hist_names[i].CompareTo("hJet1q_neEmEF")==0) || (hist_names[i].CompareTo("hJet1q_chEmEF")==0) ||  (hist_names[i].CompareTo("hJet1q_neHEF")==0)|| (hist_names[i].CompareTo("hJet1q_chHEF")==0) || (hist_names[i].CompareTo("hJet2q_neEmEF")==0) || (hist_names[i].CompareTo("hJet2q_chEmEF")==0) ||  (hist_names[i].CompareTo("hJet2q_neHEF")==0)|| (hist_names[i].CompareTo("hJet2q_chHEF")==0)) {
			LOGY=false;
		}
		if ((hist_names[i].CompareTo("hJetUnder50_neHEF")==0) || (hist_names[i].CompareTo("hJetUnder50_chHEF")==0) ||  (hist_names[i].CompareTo("hJetUnder50_neEmEF")==0)|| (hist_names[i].CompareTo("hJetUnder50_chEmEF")==0) || (hist_names[i].CompareTo("hJetOver50_chHEF")==0) || (hist_names[i].CompareTo("hJetOver50_neHEF")==0) ||  (hist_names[i].CompareTo("hJetOver50_neEmEF")==0)|| (hist_names[i].CompareTo("hJetOver50_chEmEF")==0)/*|| (hist_names[i].CompareTo("hJetOver50_eta")==0)|| (hist_names[i].CompareTo("hJetUnder50_eta")==0)*/) {
			LOGY=false;
		}
		
		if (hist_names[i].CompareTo("hMqq_log")==0) {
			xmin=4.5;
			xmax=9;
		}
		if (hist_names[i].CompareTo("hAdJetHT_bdt")==0) {
			xmin=-50;
			xmax=450;
		}
		if (hist_names[i].CompareTo("hbdt_atanh")==0) {
			xmin=0;
			xmax=3.5;
		}
//		if (hist_names[i].CompareTo("hAdJetHT")==0) {
//			xmax=200;
//		}

		TH1F *frame = new TH1F("frame","",1,xmin,xmax);
		TGaxis::SetExponentOffset(-0.07,0,"xy");
		frame->Reset();
		frame->SetMinimum(0.1);
                
                frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*1.2,sum_histos[i]->GetMaximum()*1.2) );
//                 frame->SetMaximum(1);
		if (LOGY==true) {
			gPad->SetLogy();	
                        frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*100,sum_histos[i]->GetMaximum()*100) );
// // // // //                         frame->SetMaximum(sum_histos[i]->GetMaximum()*1.1 );
			if (hist_names[i].CompareTo("hHT")==0) frame->SetMaximum(sum_histos[i]->GetMaximum()*100 );
		}
		TGaxis::SetMaxDigits(4);
                frame->GetXaxis()->SetTitleOffset(0.91);
                frame->SetStats(0);
		frame->GetYaxis()->SetNdivisions(505);
	 	frame->GetXaxis()->SetLabelSize(0.0);
		char name[1000];
		if (tmp_hist_name.find("GeV")==std::string::npos) {
			if (data_histos[i]->GetBinWidth(1)>1) sprintf(name,"Events / %1.0f",data_histos[i]->GetBinWidth(1));
			else sprintf(name,"Events / %1.2f",data_histos[i]->GetBinWidth(1));
		} else {
      	sprintf(name,"Events / %1.0f %s",data_histos[i]->GetBinWidth(1),"GeV");
		}
		frame->GetYaxis()->SetTitle(name);

      frame->Draw();
		tex->Draw();
		tex1->Draw();
		tex2->Draw();
                pCMSset.Draw("same");
	//	tex_k->Draw();
//		tex_set->Draw();
                
                
		if (!Zanalysis) if ((d_value>0.01)&&(d_value<1.)) disc_value_text->Draw();
                stacks[i]->Draw("same");	
		signal_histos[i]->Draw("same hist");
		gluonFu_histos[i]->Draw("same hist");
                

                
//                 TString Name_Hist_String(data_histos[i]->GetName());
//                 std::cout << Name_Hist_String << std::endl;


                float totalSensitivity=0.;
                if  (hist_names[i].CompareTo("hZll_mass")==0) {
                    if (Zanalysis) data_histos[i]->Draw("Psame");
                    else {
                        TH1F *h_below_115 = (TH1F*)data_histos[i]->Clone("Mll_below_115");
                        TH1F *h_above_130 = (TH1F*)data_histos[i]->Clone("Mll_above_130");
                        
                        h_below_115->GetXaxis()->SetRangeUser(110,120);
                        h_above_130->GetXaxis()->SetRangeUser(130,155);
                        std::cout << "HERE---------------------------------------------------------------------------------------------------------------------------" << std::endl;


                        
                        h_below_115->Draw("Psame");
                        h_above_130->Draw("Psame");
                    }
                }
                else {
                    if ((hist_names[i].CompareTo("hNNoutput")!=0) && (hist_names[i].CompareTo("hNNplusBDT_atanh")!=0) && (hist_names[i].CompareTo("hNNoutput_atanh")!=0) && (hist_names[i].CompareTo("hBDT_VBF")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125ForAll")!=0) && (hist_names[i].CompareTo("hZll_mass")!=0) && (hist_names[i].CompareTo("hHll_mass_precise")!=0) && (hist_names[i].CompareTo("hHll_mass_unprecise")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_Rpt0")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_WM1Log5")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_zStar0")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_ptll300")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_softN50")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_MqqLog75")!=0)) data_histos[i]->Draw("Psame");
                    
                    
                    
                    
                    if (hist_names[i].CompareTo("hBDT_VBF")==0) {
                        if (!Zanalysis) data_histos[i]->GetXaxis()->SetRangeUser(-1.,0.2);
                        data_histos[i]->Draw("Psame");
                    }

                    if (hist_names[i].CompareTo("hBDT_VBF_atanh")==0) {
                       if (!Zanalysis) data_histos[i]->GetXaxis()->SetRangeUser(0.,0.5);
                        data_histos[i]->Draw("Psame");
                    }
                    
                    if (hist_names[i].CompareTo("hBDT_VBF_atanh_m125ForAll")==0) {
//                         data_histos[i]->GetXaxis()->SetRangeUser(0.,1.2);
                        data_histos[i]->Draw("Psame");
                    }
                    
                    if (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_zStar0")==0) {
//                         data_histos[i]->GetXaxis()->SetRangeUser(0.,0.9);
                        data_histos[i]->Draw("Psame");
                    }
                    
                    if ((hist_names[i].CompareTo("hBDT_VBF_atanh_m125_Rpt0")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_WM1Log5")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_ptll300")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_softN50")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_MqqLog75")==0))  {
//                         data_histos[i]->GetXaxis()->SetRangeUser(0.,1.2);
                        data_histos[i]->Draw("Psame");
                    }
                    
                    if (hist_names[i].CompareTo("hNNoutput")==0) {
                        data_histos[i]->GetXaxis()->SetRangeUser(0.,0.5);
                        data_histos[i]->Draw("Psame");
                    }

                    if (hist_names[i].CompareTo("hNNoutput_atanh")==0) {
                        data_histos[i]->GetXaxis()->SetRangeUser(0.,2.);
                        data_histos[i]->Draw("Psame");
                    }


                    if (hist_names[i].CompareTo("hNNplusBDT_atanh")==0) {
                        data_histos[i]->GetXaxis()->SetRangeUser(0.,4.);
                        data_histos[i]->Draw("Psame");
                    }
                    
                }

//                 data_histos[nhistos-2]->GetXaxis()->SetRangeUser(-1.,0.5);
//                 data_histos[nhistos-1]->GetXaxis()->SetRangeUser(0.,0.5);

//                 if (i>=nhistos-2) {
//                         for (int n = 1; n <= signal_histos[i]->GetXaxis()->GetNbins(); ++n) {
//                             float histoBin = signal_histos[i]->GetBinContent(n)/20.;
//                             float staskPos = stackHisto[i]->GetBinContent(n);
//                             if (staskPos < 0.) staskPos = 0.;
//                             totalSensitivity += (histoBin > 0.00001) ? histoBin*histoBin/(histoBin + staskPos) : 0.;
//                             std::cout <<"totalSensitivity " << n << " " << totalSensitivity <<" \t histoBin  " << histoBin  <<" \t stackHisto[i]->GetBinContent(n)  " << stackHisto[i]->GetBinContent(n) << " \t " << staskPos << std::endl;
// //                            totalSensitivity += sensitivity*sensitivity;
//                         }
                        
                        
                        
                if ((hist_names[i].CompareTo("hNNoutput")==0) || (hist_names[i].CompareTo("hNNplusBDT_atanh")==0) || (hist_names[i].CompareTo("hNNoutput_atanh")==0) || (hist_names[i].CompareTo("hBDT_VBF")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh")==0)) {
//                 if ((hist_names[i].CompareTo("hNNoutput_atanh")==0) ) {
                    TH1F * hisoSensitivity_binByBin = new TH1F (hist_names[i], "", signal_histos[i]->GetNbinsX(), 0, signal_histos[i]->GetNbinsX());
                    hisoSensitivity_binByBin->GetXaxis()->SetTitle(signal_histos[i]->GetXaxis()->GetTitle());
                    
                    float totalSensitivitySquared=0.;   
                    float totalSensitivitySquaredErrorSquared = 0.;
                        for (int n = 1; n <= signal_histos[i]->GetXaxis()->GetNbins(); ++n) {
                            float histoBin = signal_histos[i]->GetBinContent(n)/20.;  //  (signal_histos[i]->GetBinContent(n) + gluonFu_histos[i]->GetBinContent(n))/20.;
                            float histoBinGlu = gluonFu_histos[i]->GetBinContent(n)/20.;
                            float staskPos = stackHisto[i]->GetBinContent(n) - signal_histos[i]->GetBinContent(n)/20. - gluonFu_histos[i]->GetBinContent(n)/20.;
                            
                            float histoBinVBFError  = signal_histos[i]->GetBinError(n)/20.;
                            float histoBinGluError  = gluonFu_histos[i]->GetBinError(n)/20.;
                            float histoBinSigError  = sqrt(histoBinVBFError*histoBinVBFError + histoBinGluError*histoBinGluError);
                            float histoSTACK_Error  = stackHisto[i]->GetBinError(n);
                            float histoBinBKGError  = sqrt(histoSTACK_Error*histoSTACK_Error - histoBinSigError*histoBinSigError);
                            

                            
                            
                            float Signal_RelativeError = histoBinSigError/(histoBin + histoBinGlu);
                            float BKG_RelativeError = histoBinBKGError/staskPos;
                            float SensitivityRelativeError = (histoBin > 0.00001 && staskPos > 0.00001) ? (4*Signal_RelativeError*Signal_RelativeError + BKG_RelativeError*BKG_RelativeError) : 0.;
                            float SensitivitySquaredErrorSquared = SensitivityRelativeError*histoBin/staskPos*histoBin/staskPos;
                            totalSensitivitySquaredErrorSquared += SensitivitySquaredErrorSquared;
                        
                            
//                             std::cout <<"bin " << n << " \t\t signal VBF " << histoBin  <<" \t +-  " << histoBinVBFError <<  std::endl;
//                             std::cout <<"\t\t signal GLU " << histoBinGlu  <<" \t +-  " << histoBinGluError <<  std::endl;
//                             std::cout <<"\t\t  VBF + GLU " << histoBin + histoBinGlu  <<" \t +-  " << histoBinSigError <<  std::endl;
//                             std::cout <<"\t\t stack      " << staskPos + histoBin + histoBinGlu <<" \t +-  " << histoSTACK_Error <<  std::endl;
//                             std::cout <<"\t\t BKG        " << staskPos  <<" \t +-  " << histoBinBKGError <<  std::endl;
//                             std::cout <<"\t\t relErrors  " <<  Signal_RelativeError <<" \t  " << BKG_RelativeError  << " \t  " << SensitivityRelativeError <<  std::endl;
                            
                            
//                             TFile *file_initial;
//                             file_initial = TFile::Open(file_names[7]);
//                             TH1F * h_check = (TH1F *) file_initial->Get("hBDT_VBF_atanh");
//                             staskPos = h_check->GetBinContent(n)*lumi* ratio[0];
                            
                            if (staskPos < 0.) {staskPos = 0.; std::cout << "ATTENTION: staskPos < 0.       -----------------------------------------------------------" <<  std::endl;}
                            
                            float binSensitivitySquared = (histoBin > 0.00001 && staskPos > 0.00001) ? histoBin*histoBin/(staskPos) : 0.;
                            totalSensitivitySquared += binSensitivitySquared;
                            hisoSensitivity_binByBin->SetBinContent(n, sqrt(binSensitivitySquared));
                            std::cout <<"bin " << n << " \t signal          " << histoBin  <<" \t bkg    " << staskPos << "  \t " << "bin sensitivity  " << binSensitivitySquared << " \t +- " << SensitivitySquaredErrorSquared << std::endl;
                            std::cout <<"\t Sensitivity^2   " <<   totalSensitivitySquared<< " \t +-  " << totalSensitivitySquaredErrorSquared  <<  std::endl;

                        }
                        
                        
                        hisoSensitivity_binByBin_VECTOR.push_back(hisoSensitivity_binByBin);
                        
//                        data_histos[nhistos-1]->GetXaxis()->SetRangeUser(-1.,0.5);
                        data_histos[i]->Draw("Psame");


                        std::ostringstream sigmaString;
//                         sigmaString<<std::setprecision(3)<<totalSensitivity;
                        sigmaString<<std::setprecision(2)<<sqrt(totalSensitivitySquared) ;//<< " +-" << sqrt(totalSensitivitySquaredErrorSquared/2./sqrt(totalSensitivitySquared));
                        TLatex* texSig = new TLatex(0.50,0.85,("Sensitivity = " + sigmaString.str() ).c_str());
                        texSig->SetNDC();
                        texSig->SetTextAlign(35);
                        texSig->SetTextFont(42);
                        texSig->SetTextSize(0.04);
                        texSig->SetLineWidth(2);

                        if (!Zanalysis) texSig->Draw();
                    }


		hBkgVis[i]->Draw("same E2");

		leg->Draw("same");
                    
                gPad->RedrawAxis();
	
  		TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);
 		pad2->SetTopMargin(0.73);
 		pad2->SetRightMargin(0.25);
 		pad2->SetFillColor(0);
  		pad2->SetFillStyle(0);
  		pad2->Draw();
  		pad2->cd(0);
  		gPad->SetGridy();

		TH1F *frame2 = new TH1F("frame2","",1,xmin,xmax);
                frame2->SetMinimum(-.5);
                frame2->SetMaximum(.5);
                frame2->SetMinimum(-0.9);
                frame2->SetMaximum(0.9);
                frame2->SetMinimum(-0.45);
                frame2->SetMaximum(0.45);
        
      frame2->SetStats(0);
      frame2->SetTitleFont(42,"x");
		frame2->SetTitleFont(42,"y");
      frame2->SetTitleSize(0.13, "XYZ");
		frame2->GetYaxis()->SetNdivisions(505);
 		frame2->GetYaxis()->SetTickLength(0.06);
  		frame2->GetYaxis()->SetTitleSize(0.04);
  		frame2->GetYaxis()->SetTitleOffset(1.5);
 		frame2->GetYaxis()->SetLabelSize(0.03);
  		frame2->GetYaxis()->CenterTitle(kTRUE);
  		frame2->GetXaxis()->SetTitleSize(0.05);
  		frame2->GetXaxis()->SetLabelSize(0.04);
		frame2->SetXTitle(signal_histos[i]->GetXaxis()->GetTitle());
		if (hist_names[i].CompareTo("hNAdJets_bdt2")==0 ) 
	   	frame2->SetXTitle("N of jets, BDT > 0.84");
// 		if (hist_names[i].CompareTo("hNAdJets_mjj1")==0 ) 
// 	   	frame2->SetXTitle("N of jets, m(qq) > 1500");
// 		if (hist_names[i].CompareTo("hNAdJets_mjj2")==0 ) 
// 	   	frame2->SetXTitle("N of jets, m(qq) > 2500");
		frame2->SetYTitle("Data / MC - 1");
		frame2->Draw();	

 		Double_t aa[2] = {xmin,xmax};
                Double_t bb[2] = {0,0};
                TGraph *cons = new TGraph(2,aa,bb);
                cons->SetLineStyle(2);
		cons->Draw("Lsame");
		
		data_histos2[i] = (TH1F*)data_histos[i]->Clone("new");
		data_histosQCDup[i] = (TH1F*)hBkgQCDUp[i]->Clone("newQCDup");
		data_histosQCDup[i]->SetLineColor(kBlue);
		data_histosQCDup[i]->SetLineStyle(2);
		data_histosQCDlo[i] = (TH1F*)hBkgQCDLo[i]->Clone("newQCDlow");
		data_histosQCDlo[i]->SetLineStyle(2);
		data_histosQCDlo[i]->SetLineColor(kBlue);
// 		cout<<data_histosQCDlo[i]->Integral()<< "    "<<data_histosQCDup[i]->Integral()<<endl;
		data_histosJESup[i] = (TH1F*)hBkgJESUp[i]->Clone("newJESup");
		data_histosJESup[i]->SetLineColor(kRed);
		data_histosJESup[i]->SetLineStyle(2);
		data_histosJESlo[i] = (TH1F*)hBkgJESLo[i]->Clone("newJESlow");
		data_histosJESlo[i]->SetLineStyle(2);
		data_histosJESlo[i]->SetLineColor(kRed);
		float chi2_data_mc = data_histos2[i]->Chi2Test(sum_histos[i],"UWCHI2/NDF");
		float ks = data_histos2[i]->KolmogorovTest(sum_histos[i]);
		data_histos2[i]->Add(sum_histos[i],-1);
		data_histos2[i]->Divide(sum_histos[i]);

		for (int j=0;j<data_histos2[i]->GetNbinsX();j++){
			if (sum_histos[i]->GetBinContent(j+1)!=0) data_histos2[i]->SetBinError(j+1, TMath::Sqrt(data_histos[i]->GetBinContent(j+1))/sum_histos[i]->GetBinContent(j+1)) ;
			else 
                        if (data_histos[i]->GetBinContent(j+1)!=0) data_histos2[i]->SetBinError(j+1, TMath::Sqrt(data_histos[i]->GetBinContent(j+1))/data_histos[i]->GetBinContent(j+1)) ;
					else data_histos2[i]->SetBinError(j+1,0.); 
		}

// 		 if (i!=nhistos-1) data_histos2[i]->Draw("PEsame");
//                  else {
//                         data_histos2[nhistos-1]->GetXaxis()->SetRangeUser(0.,2.4);
//                         data_histos2[nhistos-1]->Draw("Psame");
//                  }
//                 data_histos2[i]->Draw("PEsame");

		if ((hist_names[i].CompareTo("hNNoutput")!=0) && (hist_names[i].CompareTo("hNNplusBDT_atanh")!=0) && (hist_names[i].CompareTo("hNNoutput_atanh")!=0) && (hist_names[i].CompareTo("hBDT_VBF")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125ForAll")!=0) && (hist_names[i].CompareTo("hZll_mass")!=0)  && (hist_names[i].CompareTo("hHll_mass_unprecise")!=0) && (hist_names[i].CompareTo("hHll_mass_precise")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_Rpt0")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_WM1Log5")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_zStar0")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_ptll300")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_softN50")!=0) && (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_MqqLog75")!=0)) 
                 data_histos2[i]->Draw("PEsame");

        if (hist_names[i].CompareTo("hBDT_VBF")==0) {
             if (!Zanalysis) data_histos2[i]->GetXaxis()->SetRangeUser(-1.,0.2);
             data_histos2[i]->Draw("Psame");
        }

        if (hist_names[i].CompareTo("hBDT_VBF_atanh")==0) {
             if (!Zanalysis) data_histos2[i]->GetXaxis()->SetRangeUser(0.,0.5);
             data_histos2[i]->Draw("Psame");
        }
        
        if (hist_names[i].CompareTo("hBDT_VBF_atanh_m125ForAll")==0) {
//              data_histos2[i]->GetXaxis()->SetRangeUser(0.,1.2);
             data_histos2[i]->Draw("Psame");
        }
        
        
        if (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_zStar0")==0) {
//             data_histos2[i]->GetXaxis()->SetRangeUser(0.,0.9);
            data_histos2[i]->Draw("Psame");
        }
                    
                    

        if ((hist_names[i].CompareTo("hBDT_VBF_atanh_m125_Rpt0")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_WM1Log5")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_zStar0")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_ptll300")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_softN50")==0) || (hist_names[i].CompareTo("hBDT_VBF_atanh_m125_MqqLog75")==0))  {
//             data_histos2[i]->GetXaxis()->SetRangeUser(0.,1.2);
            data_histos2[i]->Draw("Psame");
        }
                    
        
        if (hist_names[i].CompareTo("hNNoutput")==0) {
             data_histos2[i]->GetXaxis()->SetRangeUser(0.,0.5);
             data_histos2[i]->Draw("Psame");
        }

        if (hist_names[i].CompareTo("hNNoutput_atanh")==0) {
             data_histos2[i]->GetXaxis()->SetRangeUser(0.,2.);
             data_histos2[i]->Draw("Psame");
        }
        
        
        if (hist_names[i].CompareTo("hNNplusBDT_atanh")==0) {
             data_histos2[i]->GetXaxis()->SetRangeUser(0.,4.);
             data_histos2[i]->Draw("Psame");
        }
        
        

		if (hist_names[i].CompareTo("hZll_mass")==0) {
// 		if (hist_names[i].CompareTo("NOhZllMASS")==0) {
            data_histos2[i]->GetXaxis()->SetRangeUser(100.,115);
            TH1F *h_below_115_2 = (TH1F*)data_histos2[i]->Clone("Mll_below_115");
            TH1F *h_above_130_2 = (TH1F*)data_histos2[i]->Clone("Mll_above_130");

            h_below_115_2->GetXaxis()->SetRangeUser(110,120);
            h_above_130_2->GetXaxis()->SetRangeUser(130,155);

            h_below_115_2->Draw("Psame");
            h_above_130_2->Draw("Psame");
        }


		data_histosQCDup[i]->Add(sum_histos[i],-1);
		data_histosQCDup[i]->Divide(sum_histos[i]);
		data_histosQCDup[i]->Draw("HISTsame");
		data_histosQCDlo[i]->Add(sum_histos[i],-1);
		data_histosQCDlo[i]->Divide(sum_histos[i]);
		data_histosQCDlo[i]->Draw("HISTsame");
		data_histosJESup[i]->Add(sum_histos[i],-1);
		data_histosJESup[i]->Divide(sum_histos[i]);
		data_histosJESup[i]->Draw("HISTsame");
		data_histosJESlo[i]->Add(sum_histos[i],-1);
		data_histosJESlo[i]->Divide(sum_histos[i]);
		data_histosJESlo[i]->Draw("HISTsame");
		hBkgUncUp[i]->Draw("HIST same");
		hBkgUncLo[i]->Draw("HIST same");
		
		temp_str.Form("%2.2f",chi2_data_mc);
		TString chi2_str = temp_str.Prepend("#chi^{2} = ");
		TLatex *chi2_latex = new TLatex(0.21,0.23,temp_str);
 	  	chi2_latex->SetNDC();
  	 	chi2_latex->SetTextAlign(20);
  		chi2_latex->SetTextFont(42);
   	chi2_latex->SetTextSize(0.025);
   	chi2_latex->SetLineWidth(2);
		chi2_latex->Draw();
		temp_str.Form("%2.2f",ks);
		TString ks_str = temp_str.Prepend(", KS = ");
		TLatex *ks_latex = new TLatex(0.31,0.23,temp_str);
 	  	ks_latex->SetNDC();
  	 	ks_latex->SetTextAlign(20);
  		ks_latex->SetTextFont(42);
   	ks_latex->SetTextSize(0.025);
   	ks_latex->SetLineWidth(2);
		ks_latex->Draw();
		TLegend *leg_jes = new TLegend(0.16,0.14,0.3,0.17); //without writing about SF
		leg_jes->SetFillStyle(0);
		leg_jes->SetBorderSize(0);
		leg_jes->SetTextFont(42);
		leg_jes->SetTextSize(0.02);
		leg_jes->AddEntry(data_histosQCDup[i],"QCD scale up/down","L");
		leg_jes->AddEntry(data_histosJESup[i],"JES up/down","L");
		leg_jes->Draw();
		


		pad2->RedrawAxis();
		c1->Print(output_names[i]);
		c1->Delete();
	}



// TCanvas *canv = new TCanvas("canv","",900,750);	
// for (int n = 0; n < hisoSensitivity_binByBin_VECTOR.size(); n++) {
//     
// 
//     canv->cd();
//     gStyle->SetOptStat(0000);
//     hisoSensitivity_binByBin_VECTOR[n]->SetLineColor(4);
//     hisoSensitivity_binByBin_VECTOR[n]->SetLineWidth(2);
//     hisoSensitivity_binByBin_VECTOR[n]->Draw();
//     std::string nameSig = hisoSensitivity_binByBin_VECTOR[n]->GetName();
//     canv->Print(("plotsDirectory/"+nameSig+"_BinByBinSensitivity_"+std::to_string(n)+".png").c_str());
// }

	
	
return 0;
}
