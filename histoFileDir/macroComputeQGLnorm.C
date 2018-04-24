#include <ctype.h>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm> 
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
#include <TProfile.h>
#include "math.h"




void macroComputeQGLnorm() {
    
    
    const int fileNumber = 36;
    TString nomiDeiFile[fileNumber] = {"DY0JetsToLL_M", "DY1JetsToLL_M", "DY2JetsToLL_M", "DY3JetsToLL_M", "DY4JetsToLL_M", "DYJetsToLL_M-105To160-madgraphMLM", "DYJetsToLL_M", "DYJetsToLL_Pt-0To50_amc", "DYJetsToLL_Pt-100To250_amc", "DYJetsToLL_Pt-250To400_amc", "DYJetsToLL_Pt-400To650_amc", "DYJetsToLL_Pt-50To100_amc", "DYJetstoLL_HT100_200", "DYJetstoLL_HT1200_2500", "DYJetstoLL_HT200_400", "DYJetstoLL_HT2500_Inf", "DYJetstoLL_HT400_600", "DYJetstoLL_HT600_800", "DYJetstoLL_HT800_1200", "DYJetstoLL_amc_0J", "DYJetstoLL_amc_1J", "DYJetstoLL_amc_2J", "DYJetstoLL_amc_Filter105", "DYJetstoLL_amc_M-50", "GluGlu_HToMuMu", "ST_s-channel", "ST_t-channel_antitop_4f_inclusiveDecays", "ST_t-channel_top_4f_inclusiveDecays", "ST_tW_antitop", "ST_tW_top", "TT", "VBF_HToMuMu", "WJetsToLNu", "WW", "WZ", "ZZ"};
    
    
//         const int fileNumber = 8;
//     TString nomiDeiFile[fileNumber] = {"ST_tW_antitop", "ST_tW_top", "TT", "VBF_HToMuMu", "WJetsToLNu", "WW", "WZ", "ZZ"};


    ofstream out_txt;
    out_txt.open("qgl_normalizations.txt"); 
    

    for (int curr_file=0; curr_file < fileNumber; curr_file++) {
        TFile * f_qglCorrected      = TFile::Open(nomiDeiFile[curr_file]+"_mu_QCDScalenom_JESnom_v25_reskim.root","read");
        TFile * f_NoqglCorrected    = TFile::Open(nomiDeiFile[curr_file]+"_mu_QCDScalenom_JESnom_v25_reskim_noQGLcorrection.root","read");
        
        
        TString histname = "hmet";
        TH1F * h_qglCorrected       = (TH1F*) f_qglCorrected->Get(histname);
        TH1F * h_NoqglCorrected     = (TH1F*) f_NoqglCorrected->Get(histname);
        
        
        float I_qglCorrected        = h_qglCorrected->Integral(0, h_qglCorrected->GetNbinsX()+1);
        float I_NoqglCorrected      = h_NoqglCorrected->Integral(0, h_qglCorrected->GetNbinsX()+1);
        
        std::cout << "qgl_norm[\"" << nomiDeiFile[curr_file] << "\"]=" << I_NoqglCorrected/I_qglCorrected << ";" << std::endl;
        
        out_txt << "qgl_norm[\"" << nomiDeiFile[curr_file] << "\"]=" << I_NoqglCorrected/I_qglCorrected << ";" << std::endl;
        
    }
    
    
    out_txt.close();
    
}





