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




float findMean (TH1F * h) {
    float mean = 0;
    for (int n=1; n <= h->GetNbinsX(); n++) mean += h->GetBinContent(n);
    return mean/h->GetNbinsX();
}

float findRMS (TH1F * h) {
    float RMS = 0;
    float mean = findMean(h);
    for (int n=1; n <= h->GetNbinsX(); n++) RMS += (h->GetBinContent(n) - mean) * (h->GetBinContent(n) - mean);
    RMS = sqrt( RMS/(h->GetNbinsX()-1) );
    return RMS;
}

void write_PDF_variation(TH1F * hNNoutput_atanh, TH1F * hNNoutput_atanh_PDFvariation, int n, TFile* fileHist) {
    
    std::string name = hNNoutput_atanh->GetName();
    TH1F * hNNoutput_atanh_Up   = (TH1F*) hNNoutput_atanh->Clone((name+"_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Up").c_str());
    TH1F * hNNoutput_atanh_Down = (TH1F*) hNNoutput_atanh->Clone((name+"_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Down").c_str());
    
    
    float RMS = findRMS(hNNoutput_atanh_PDFvariation);
//     std::cout << "bin value     " << hNNoutput_atanh_Up->GetBinContent(n)<< " \t mean   " << findMean(hNNoutput_atanh_PDFvariation)  << " \t RMS   " << RMS << std::endl;
    hNNoutput_atanh_Up->SetBinContent  (n, hNNoutput_atanh_Up->GetBinContent(n)   + RMS);
    hNNoutput_atanh_Down->SetBinContent(n, hNNoutput_atanh_Down->GetBinContent(n) - RMS);
    
    fileHist->cd();
    hNNoutput_atanh_Up->Write();              // they are automaticaly saved becaused they are clone of something already saved
    hNNoutput_atanh_Down->Write();
}




void macroAddNNdiscriminator () {
    
    std::cout << "I can see here" << std::endl;


    
//     std::string pathTree= "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/mvaTree/";
    std::string pathTree= "/scratch/mandorli/Hmumu/AddNNinformationToMVAnTuples/CMSSW_9_4_4/src/leonardoCode/ROOT/";
    std::string pathHistoFile= "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/histoFileDir/";

//     const int nfiles = 15;
//     std::string file_names[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST_tW_antitop", "ST_tW_top", "ST_s-channel", "ST_t-channel_antitop_4f_inclusiveDecays", "ST_t-channel_top_4f_inclusiveDecays","TT","DYJetsToLL_M-105To160-madgraphMLM","DYJetsToLL_M-105To160-amcatnloFXFX","VBF_HToMuMu", "GluGlu_HToMuMu"};

//     const int nfiles = 15;
//     std::string file_names[nfiles] = {"SingleMuon","WToLNu_2J","WW","ZZ","WZ","ST_tW_antitop", "ST_tW_top", "ST_s-channel", "ST_t-channel_antitop_4f_inclusiveDecays", "ST_t-channel_top_4f_inclusiveDecays","TT","DYJetsToLL_M-105To160-madgraphMLM","DYJetsToLL_M-105To160-amcatnloFXFX","VBF_HToMuMu", "GluGlu_HToMuMu"};
    
    
    const int nfiles = 18;  
    std::string file_names[nfiles] = {"SingleMuon","WToLNu_2J","W3JetsToLNu", "W4JetsToLNu","WW","ZZ","WZ","ST_tW_antitop", "ST_tW_top", "ST_t-channel_antitop_4f_inclusiveDecays", "ST_t-channel_top_4f_inclusiveDecays","TT", "TTToSemilepton", "TTTo2L2Nu","DYJetsToLL_M-105To160-madgraphMLM","DYJetsToLL_M-105To160-amcatnloFXFX","VBF_HToMuMu", "GluGlu_HToMuMu"};
    

    
    
//     const int nfiles = 10; //Mqq > 1000
//     std::string file_names[nfiles] = {"SingleMuon","WToLNu_2J","WW","ST_tW_antitop", "ST_t-channel_antitop_4f_inclusiveDecays","TT","DYJetsToLL_M-105To160-madgraphMLM","DYJetsToLL_M-105To160-amcatnloFXFX","VBF_HToMuMu", "GluGlu_HToMuMu"};
    
    
    
    
//     const int nfiles = 1;
//     std::string file_names[nfiles] = {"VBF_HToMuMu"};
    
    
    std::string qcdScale[9] = {"nom", "up", "down", "nom", "nom", "nom", "nom", "nom", "nom"};
    std::string jesScale[9] = {"nom", "nom", "nom", "up", "down", "nom", "nom", "nom", "nom"};
    std::string jerScale[9] = {"nom", "nom", "nom", "nom", "nom", "up", "down", "nom", "nom"};
    std::string puScale[9] =  {"nom", "nom", "nom", "nom", "nom", "nom", "nom", "up", "down"};
    

        
//         const int binNUmberNN = 31;   // NN in Agense thesis
//         float NN_bin[binNUmberNN] = {0,     0.149,  0.299,  0.449,  0.599,  0.749,  0.899,  1.049,  1.199,  1.349,  1.499,  1.649,  1.799,  1.949,  2.099,  2.249,  2.399, 2.549,  2.699,  2.849,  2.999,  3.149,  3.299,  3.449,  3.599,  3.749,  3.9115,         4.0615,         4.2115,         4.4385,       5};
        
        
        
        const int binNUmberNN = 31;   // NN with nano
        float NN_bin[binNUmberNN] = {0,     0.1025,         0.2525,         0.4025,         0.5525,         0.7025,         0.8525,         1.0025,         1.1525,         1.3025,         1.4525,       1.6025,         1.7525,         1.9025,         2.0525,         2.2025,         2.3525,         2.5025,         2.6525,         2.8025,         2.9525,       3.1025,         3.2525,         3.4025,         3.5525,         3.7025,         3.8525,         4.0025,         4.1525,         4.3025,         5};
        
//         const int binNUmberNN = 26;   // NN withou mll
//         float NN_bin[binNUmberNN] = {0,     0.144,  0.294,  0.444,  0.594,  0.744,  0.894,  1.044,  1.194,  1.344,  1.494,  1.644,  1.794,  1.944,  2.094,  2.244,  2.394,  2.544,  2.694,    2.844,  2.994,  3.144,  3.2975,         3.678,  4.126,  5};
            
        
        
//         const int binNUmberNN = 31;   // NN with Mqq > 200 and Mqq > 250 and Mqq > 300 and Mqq > 400 and Mqq > 600
//         float NN_bin[binNUmberNN] = {0,     0.149,  0.299,  0.449,  0.599,  0.749,  0.899,  1.049,  1.199,  1.349,  1.499,  1.649,  1.799,  1.949,  2.099,  2.249,  2.399, 2.549,  2.699,  2.849,  2.999,  3.149,  3.299,  3.449,  3.599,  3.749,  3.9115,         4.0615,         4.2115,         4.4385, 5};    
            
        
//         const int binNUmberNN = 31;   // NN with Mqq > 800 
//         float NN_bin[binNUmberNN] = {0,     0.1335,         0.2835,         0.4335,         0.5835,         0.7335,         0.8835,         1.0335,         1.1835,         1.3335,  1.4835,  1.6335,         1.7835,         1.9335,         2.0835,         2.2335,         2.3835,         2.5335,         2.6835,         2.8335,  2.9835,  3.1335,         3.2835,         3.4335,         3.5835,         3.7335,         3.9115,         4.0615,         4.2115,         4.4385, 5};
            
 
/*        const int binNUmberNN = 31;   // NN with Mqq > 200 and Mqq > 250 and Mqq > 300 and Mqq > 400 and Mqq > 600
        float NN_bin[binNUmberNN] = {0,     0.002,  0.152,  0.302,  0.452,  0.602,  0.752,  0.902,  1.052,  1.202,  1.352,  1.502,  1.652,  1.802,  1.952,  2.102,  2.252,  2.402,  2.552,    2.702,  2.852,  3.002,  3.152,  3.302,  3.452,  3.602,  3.752,  4.046,  4.2115,         4.4385,         5};  */ 
            


        
        
        for (int uncIdx = 0; uncIdx < 9; ++uncIdx) {
//        for (int jesIdx = 0; jesIdx < 5; ++jesIdx) { 
        
        
            for(int n = 0; n < nfiles; n++) {
                
                if (file_names[n].compare("SingleMuon")==0 && uncIdx != 0) continue;
                
                
                std::string fileWithTree = file_names[n];
                std::string fileWithHist = file_names[n];
                
                fileWithTree =  pathTree + "main_tmva_tree_" + fileWithTree + "_v25mu_QCDScale"+qcdScale[uncIdx]+"_JES"+jesScale[uncIdx]+"_JER"+jerScale[uncIdx]+"_PU"+puScale[uncIdx]+".root";
//                 fileWithHist =  pathHistoFile + fileWithHist + "_mu_QCDScale"+qcdScale[uncIdx]+"_JES"+jesScale[uncIdx]+"_JER"+jerScale[uncIdx]+"_PU"+puScale[uncIdx]+"_v25_reskim.root";
                fileWithHist =  pathHistoFile + fileWithHist + "_mu_QCDScale"+qcdScale[uncIdx]+"_JES"+jesScale[uncIdx]+"_JER"+jerScale[uncIdx]+"_PU"+puScale[uncIdx]+"_nano_2016.root";
                
                std::cout << "root -l " << fileWithTree << std::endl;
                TFile* fileTree= new TFile (fileWithTree.c_str(), "read");
                
            
            
                TH1F *hNNoutput=new TH1F("hNNoutput","",25,0,1);
                hNNoutput->GetXaxis()->SetTitle("NN output");        
//                 TH1F *hNNoutput_atanh=new TH1F("hNNoutput_atanh","",25,0,5);
                TH1F *hNNoutput_atanh=new TH1F("hNNoutput_atanh","",binNUmberNN-1, NN_bin);
                hNNoutput_atanh->GetXaxis()->SetTitle("tanh^{-1}(NN output)");   
                TH1F *hNNoutput_atanh_findBinning=new TH1F("hNNoutput_atanh_findBinning","",10000,0,5);
                hNNoutput_atanh_findBinning->GetXaxis()->SetTitle("tanh^{-1}(NN output)");   
                
                
                TH1F *hNNplusBDT_atanh=new TH1F("hNNplusBDT_atanh","",binNUmberNN-1, 0, NN_bin[binNUmberNN-1]+4);
                hNNplusBDT_atanh->GetXaxis()->SetTitle("tanh^{-1}(NN output) + 3 x tanh^{-1}((BDT output + 1)/2 )");  
                
                TH2F *hNN_Vs_BDT_atanh=new TH2F("hNN_Vs_BDT_atanh","",20, 0, 4, 32, 0, 1.6);
                hNN_Vs_BDT_atanh->GetXaxis()->SetTitle("tanh^{-1}(NN output)");
                hNN_Vs_BDT_atanh->GetYaxis()->SetTitle("tanh^{-1}((BDT output + 1)/2 )");
                
                
                const int number_LHE_weights_pdf = 102;
                
                std::vector<TH1F*> hNNoutput_atanh_PDFvariation;
                for (int n = 0; n <= binNUmberNN; n++) {
                    TH1F *hNNoutput_atanh_PDFvariation_toAdd = new TH1F(("hNNoutput_atanh_PDFvariation_toAdd"+std::to_string(n)).c_str(),"",number_LHE_weights_pdf, 0, number_LHE_weights_pdf);
                    hNNoutput_atanh_PDFvariation_toAdd->GetXaxis()->SetTitle(("tanh^{-1}(( NN output) variation in bin "+std::to_string(n)).c_str());
                    hNNoutput_atanh_PDFvariation.push_back(hNNoutput_atanh_PDFvariation_toAdd);
                }
    
    
                   
                                
                                
                                
                                
                                
                                
                                
                
                float NNout = 0;
                float genweight = 1;
                float BDToutput = 0;
                
                TTree * tree = (TTree*) fileTree->Get("tree");
                std::cout << fileTree <<  " \t "<< tree <<  " \t "<<  tree->GetEntries() << std::endl;
                tree->SetBranchAddress("NNout",&NNout); 
                tree->SetBranchAddress("BDToutput",&BDToutput); 
                tree->SetBranchAddress("genweight",&genweight); 

                float genweightVECTOR[number_LHE_weights_pdf] ;
                for (int n = 0; n < number_LHE_weights_pdf; n++) tree->SetBranchAddress(("genweightVECTOR_"+std::to_string(n)).c_str(),&genweightVECTOR[n]); 
                
                int nentries = tree->GetEntries();
                std::cout << fileTree <<  " \t "<< tree <<  " \t "<< nentries << std::endl;
                
                
                
                
                
                bool savePDFunc = false;
                if (uncIdx == 0 && file_names[n].compare("VBF_HToMuMu")==0) savePDFunc = true;
//                 std::cout << "savePDFunc     " << savePDFunc << " \t " << (uncIdx == 0) << " \t " << (file_names[n].compare("VBF_HToMuMu")==0)  << std::endl;
                for (int entry=0; entry<nentries;++entry){
            
                    tree->GetEntry(entry);
            
                    hNNoutput->Fill(NNout, genweight);
                    hNNoutput_atanh->Fill(atanh(NNout), genweight);
                    hNNoutput_atanh_findBinning->Fill(atanh(NNout), genweight);
                    
                    
                    hNNplusBDT_atanh->Fill(atanh(NNout) + 3*atanh((BDToutput+1.)/2.), genweight);
                    hNN_Vs_BDT_atanh->Fill(atanh(NNout) , atanh((BDToutput+1.)/2.), genweight);
                    
                    if (savePDFunc) {   
                        int binNumber = hNNoutput_atanh->FindBin(atanh(NNout));
                        for (int n = 0; n <= number_LHE_weights_pdf; n++) {
                            hNNoutput_atanh_PDFvariation[binNumber]->Fill(n,genweight*genweightVECTOR[n]);
                        }
                    }
                }
     

                    
            
                
                TFile* fileHist= new TFile (fileWithHist.c_str(), "update");
                fileHist->cd();
                
                hNNoutput->Write();
                hNNoutput_atanh->Write();
                hNNplusBDT_atanh->Write();
//                 hNN_Vs_BDT_atanh->Write();
                
                if (savePDFunc) {
                for (int n = 0; n < hNNoutput_atanh_PDFvariation.size(); n++) {
                        hNNoutput_atanh_PDFvariation[n]->Write();
                        write_PDF_variation(hNNoutput_atanh, hNNoutput_atanh_PDFvariation[n], n, fileHist);
                    }
                }  
                    
                
                hNNoutput_atanh_findBinning->Write();
                
                fileHist->Close();
            
        
            
                
            }
//         }
    }
    
    
    std::cout << " FINE! " << std::endl;
    
    
}



        
        
        
        
