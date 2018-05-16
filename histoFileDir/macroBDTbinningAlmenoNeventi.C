#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"
#include "TROOT.h"




int  FindBinDown(TH1F *hBDT_VBF_atanh_findBinning, int binLimitUp, float minNumberOfEventPerBin, int MinNumberOfBin_inBinning) {
        
    int binLimitDown = 0.;
        for(int n = binLimitUp-MinNumberOfBin_inBinning; n > 0; n--) {
            if (hBDT_VBF_atanh_findBinning->Integral(n+1, binLimitUp) >= minNumberOfEventPerBin ) {
                binLimitDown = n;
                break;
            }
             
        }
        
    return binLimitDown;
}












void macroBDTbinningAlmenoNeventi () {
    std::cout << "I can see here" << std::endl;

    float lumi = 35900.;
    float xMax = 5;
    
    ofstream out_bdt_binning;
    out_bdt_binning.open("bdt_binning.txt"); 

    TFile* fileDY= new TFile ("DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScalenom_JESnom_v25_reskim.root");
//     TH1F *hBDT_VBF_atanh_findBinning= (TH1F*) fileDY->Get("hBDT_VBF_atanh_findBinning");
    TH1F *hBDT_VBF_atanh_findBinning= (TH1F*) fileDY->Get("hNNoutput_atanh_findBinning");
    hBDT_VBF_atanh_findBinning->Scale(lumi);


    
    
    
    
    
    
    
    if(1) { //THIS CODE COMPUTE THE BINS WITH AT LEAST "minNumberOfEventPerBin" EVENTS AND AT LEAST "binMinWidth" BIG
        
        float minNumberOfEventPerBin = 0.5;
        float binMinWidth = 0.15;
        
        std::vector<float> binning_BDT;
        
        float Nbins_binning = hBDT_VBF_atanh_findBinning->GetNbinsX();
        int MinNumberOfBin_inBinning = binMinWidth/xMax*Nbins_binning;
            
        int binLimitUp   = Nbins_binning+1;
        int binLimitDown = Nbins_binning;

    
        int binLimitDown_Nth= Nbins_binning;
        int binLimitDown_Nplus1th= Nbins_binning;


        std::cout << "prima del while     " << binLimitDown  << " \t " << Nbins_binning <<  std::endl;
        while(binLimitDown>0) {
            binning_BDT.push_back((1.*binLimitDown*xMax)/Nbins_binning);
            binLimitUp = binLimitDown;
            binLimitDown_Nth        = FindBinDown(hBDT_VBF_atanh_findBinning, binLimitUp, minNumberOfEventPerBin, MinNumberOfBin_inBinning);
            binLimitDown_Nplus1th   = FindBinDown(hBDT_VBF_atanh_findBinning, binLimitUp, minNumberOfEventPerBin+1, MinNumberOfBin_inBinning);
            binLimitDown = (binLimitDown_Nth+binLimitDown_Nplus1th)/2;
            binLimitDown = binLimitDown_Nth;
            
//             if(binning_BDT.size()>1) std::cout << ((int) (binning_BDT[binning_BDT.size()-1]*Nbins_binning/xMax))<< " - " << ((int) (binning_BDT[binning_BDT.size()-2]*Nbins_binning/xMax)) << " \t " << binning_BDT[binning_BDT.size()-1] << " - " << binning_BDT[binning_BDT.size()-2] << " \t" << hBDT_VBF_atanh_findBinning->Integral(((int) (binning_BDT[binning_BDT.size()-1]*Nbins_binning/xMax)), ((int) (binning_BDT[binning_BDT.size()-2]*Nbins_binning/xMax))) << std::endl;
                        
        }
                        
                        
        binning_BDT.push_back(0.);

        
        for(int n = 0; n < binning_BDT.size()-1; n++) {
            
             std::cout << binning_BDT[n+1] << " - " << binning_BDT[n] << " \t" << hBDT_VBF_atanh_findBinning->Integral(((int) (binning_BDT[n+1]*Nbins_binning/xMax)), ((int) (binning_BDT[n]*Nbins_binning/xMax))) << std::endl;
            
        }
        
        std::cout << "binning_BDT.size()    " << binning_BDT.size() << std::endl;
        ;
        std::cout << "{";
        for(int n = binning_BDT.size()-1; n >= 0; n--) std::cout << binning_BDT[n] << ", \t";
        std::cout  << "};" << std::endl;
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    if(0) {  // THIS CODE MAKE THE LAST BIN START FROM THE FIRST EMPTY BIN. ALL THE OTHER BINS ARE EQUAL
        
        int firstBinVuoto = 1;
        while (hBDT_VBF_atanh_findBinning->GetBinContent(firstBinVuoto) > 0)
            firstBinVuoto++;
            
            
        float sum = 0;
        int bin_piu_a_sinistra_di_quelli_uniti = firstBinVuoto;
        while (sum <= 0) {
            bin_piu_a_sinistra_di_quelli_uniti--;
            sum = hBDT_VBF_atanh_findBinning->Integral(bin_piu_a_sinistra_di_quelli_uniti, hBDT_VBF_atanh_findBinning->GetNbinsX());
        }
        
        for (int n = 0; n < hBDT_VBF_atanh_findBinning->GetNbinsX(); n++) {
            std::cout << n << " \t " << hBDT_VBF_atanh_findBinning->GetBinContent(n) << std::endl;

        }
                
        std::cout << sum << " \t " << hBDT_VBF_atanh_findBinning->GetBinContent(bin_piu_a_sinistra_di_quelli_uniti) << std::endl;

        std::cout << "const int binNUmberBDT = " <<  bin_piu_a_sinistra_di_quelli_uniti+1 << std::endl;
        std::cout << "float BDT_bin[binNUmberBDT] = { ";
        for (int n = 0; n < bin_piu_a_sinistra_di_quelli_uniti; n++) {
            
            std::cout << n*xMax/hBDT_VBF_atanh_findBinning->GetNbinsX() << ", " ;
            
        }
        std::cout << "4}" << std::endl;

    }
    
}



        
        
        
        
