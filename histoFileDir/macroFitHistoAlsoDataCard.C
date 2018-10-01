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

#include "dataCard.h"


void addErrorToHist (TH1F * histo, float sigma, int n) {
        histo->SetBinContent(n, histo->GetBinContent(n) + histo->GetBinError(n)*sigma);
}


void writePDFErrors (TString BDTversionName, TString fileTag, TFile * file_FitHisto, TFile * file, float lumi) {
    
    
    
    for (int n = 0; n <= ((TH1F*) file->Get(BDTversionName))->GetNbinsX(); ++n) {
        TString nameUp  = BDTversionName;
        TString nameDown  = BDTversionName;
        nameUp  = nameUp.Append("_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Up");
        nameDown  = nameDown.Append("_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Down");
        
        TH1F * histo_PDF_Up = (TH1F*) file->Get(nameUp)->Clone();
        histo_PDF_Up->Scale(lumi);  
        TH1F * histo_PDF_Down = (TH1F*) file->Get(nameDown)->Clone();
        histo_PDF_Down->Scale(lumi);
    
        file_FitHisto->cd();
        histo_PDF_Up->Write();
        histo_PDF_Down->Write();
    
    }
}





void writeStatisticalErrors (TString BDTversionName, TString fileTag, TFile * file_FitHisto, TFile * file, float lumi) {
    

            
            
    
        TH1F * histo_BDT_nom = (TH1F*) file->Get(BDTversionName)->Clone();
        
        for (int n = 0; n < histo_BDT_nom->GetNbinsX()+1; ++n) {
            TH1F * histoStat_BDT_Up = (TH1F*) file->Get(BDTversionName)->Clone();
            TH1F * histoStat_BDT_Down = (TH1F*) file->Get(BDTversionName)->Clone();

            std::string binNumber = std::to_string(n);
            
            histoStat_BDT_Up->SetName((BDTversionName+"_mu_"+fileTag+"_Stat"+fileTag+binNumber+"_Up").Data());
            histoStat_BDT_Down->SetName((BDTversionName+"_mu_"+fileTag+"_Stat"+fileTag+binNumber+"_Down").Data());
            
//             histoStat_BDT_Up->SetName((BDTversionName+"_mu_"+fileTag+"_Stat"+binNumber+"_Up").Data());
//             histoStat_BDT_Down->SetName((BDTversionName+"_mu_"+fileTag+"_Stat"+binNumber+"_Down").Data());
            

            
            addErrorToHist(histoStat_BDT_Up,1,n);
            addErrorToHist(histoStat_BDT_Down,-1,n);
            
            
            histoStat_BDT_Up->Scale(lumi);
            histoStat_BDT_Down->Scale(lumi);
        
                
            file_FitHisto->cd();
                        
            histoStat_BDT_Up->Write();
            histoStat_BDT_Down->Write();
        }
}




float compute_DY_normalization(int eventData, std::vector<float> normalization_vector) {

    float sumMC = 0;
    for (int n=0;n<normalization_vector.size();n++) {
        std::cout << n << " \t " << normalization_vector[n] << std::endl;
        sumMC = sumMC + normalization_vector[n];
    }
    sumMC = sumMC - normalization_vector[normalization_vector.size()-1];                
    return eventData - sumMC;
}



void macroFitHistoAlsoDataCard () {
    std::cout << "I can see here" << std::endl;

    float lumi = 35900.;
    float FIXED_lumi = 35900.;
//     FIXED_lumi = 35800 + 44980;   // 2016 + 1017
//     FIXED_lumi = 120000;          // Expected    2016 + 1017 + 2018
//     FIXED_lumi = 300000;          // 300 fb -1  
//     FIXED_lumi = 3000000;          // 3000 fb -1  
    
    float scaleLumi = FIXED_lumi/lumi;
    
    std::string discriminator = "BDT";
    discriminator = "NN";
    
    
    

//////////////////// TIGHT PRESELECTION 15 MAY 2018  /////////////////////////////////////////////////////////////////////////
        float k_factor_nom      = 0.96821916;
        float k_factor_QCDup    = 1.0551188     / k_factor_nom;    // I have to divide by k_factor_nom because DY lumi is multiplied by k_factor_nom
        float k_factor_QCDdown  = 0.90215605    / k_factor_nom;
        float k_factor_JESup    = 0.81042933    / k_factor_nom;
        float k_factor_JESdown  = 1.1614661     / k_factor_nom;




        
    const int Nsyst = 4;
    int Nsyst_NoConst = Nsyst;
    //TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","JES","JER","QGL"};
    TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","QGL"};
    std::vector<std::vector<std::string>> uncVariation;
    uncVariation.push_back({"nom", "up", "down", "nom", "nom", "nom", "nom", "nom", "nom"});
    uncVariation.push_back({"nom", "nom", "nom", "up", "down", "nom", "nom", "nom", "nom"});
    uncVariation.push_back({"nom", "nom", "nom", "nom", "nom", "up", "down", "nom", "nom"});
    uncVariation.push_back({"nom", "nom", "nom", "nom", "nom", "nom", "nom", "up", "down"});


    
    TString DY_name = "DYJetstoLL_amc";
    TString DY_file_name = "DYJetstoLL_amc";



    
//     const int nfiles  = 9;
    // TString file_names[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
//     TString file_names[nfiles] = {"WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
//     TString file_names_QCDup[nfiles];
//     TString file_names_QCDdown[nfiles];
//     TString file_names_JESup[nfiles];
//     TString file_names_JESdown[nfiles];
//     TString file_names_JERup[nfiles];
//     TString file_names_JERdown[nfiles];
//         TString file_names_JERup[nfiles];
//     TString file_names_JERdown[nfiles];
//     
// //     TString file_data_name = "WW";
    

    std::vector<std::vector<float>> eventMC;
    std::vector<std::vector<TString>> file_names;
//     TString file_tag[nfiles]= {"WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
//     std::vector<TString> file_tag= {"VBF_HToMuMu","WJetsToLNu","WW","ZZ","WZ","ST","TT","GluGlu_HToMuMu",DY_file_name};  // DY is the last one because k-factor has to be computed first
    std::vector<TString> file_tag= {"VBF_HToMuMu","WW"/*,"ZZ","WZ"*/,"ST","TT","GluGlu_HToMuMu",DY_file_name};  // DY is the last one because k-factor has to be computed first
    int nfiles  = file_tag.size();

    
    
    for (int j=0;j<uncVariation[0].size();j++) {
        file_names.push_back(file_tag);
        eventMC.push_back({0.});
        for (int i=1;i<nfiles;i++) eventMC[j].push_back({0.});
    }
        
        
    for (int i=0;i<nfiles;i++){
        for (int j=0;j<uncVariation[0].size();j++){
            
            
//              file_names[i][j] = file_tag[i];
             file_names[j][i].Append("_mu");
             file_names[j][i].Append("_QCDScale"+uncVariation[0][j]+"_JES"+uncVariation[1][j]+"_JER"+uncVariation[2][j]+"_PU"+uncVariation[3][j]+"_v25_reskim.root");
             
             
//         file_tag[i] = file_names[i];
//  
//         file_names[i].Append("_mu");
//         file_names_QCDup[i]   = file_names[i];
//         file_names_QCDdown[i] = file_names[i];
//         file_names_JESup[i]   = file_names[i];
//         file_names_JESdown[i] = file_names[i];
//         file_names_JERup[i]   = file_names[i];
//         file_names_JERdown[i] = file_names[i];
//         
//         file_names[i].Append("_QCDScalenom_JESnom_v25_reskim.root");
//         file_names_QCDup[i].Append("_QCDScaleup_JESnom_v25_reskim.root");
//         file_names_QCDdown[i].Append("_QCDScaledown_JESnom_v25_reskim.root");
//         file_names_JESup[i].Append("_QCDScalenom_JESup_v25_reskim.root");
//         file_names_JESdown[i].Append("_QCDScalenom_JESdown_v25_reskim.root");
//         file_names_JERup[i].Append("_QCDScaleup_JESup_v25_reskim.root");
//         file_names_JERdown[i].Append("_QCDScaledown_JESdown_v25_reskim.root");
        }
    }
    
      
      
      
    TString histTitle;
    TString hTitle;
    TString hTitle_atanh;
    if (discriminator.compare("BDT") == 0) {
        histTitle       = "hBDT_VBF_";
        hTitle          = "hBDT_VBF";
        hTitle_atanh    = "hBDT_VBF_atanh";
    }

    if (discriminator.compare("NN") == 0) {
        histTitle = "hNNoutput_";
        hTitle          = "hNNoutput";
        hTitle_atanh    = "hNNoutput_atanh";
    }

///////////////////////// DATA ///////////////////////////////////////
    TString file_data_name = "SingleMuon_mu_QCDScalenom_JESnom_JERnom_PUnom_v25_reskim.root";
    std::cout <<  "Opening    " << file_data_name << std::endl;
    TFile * file_data   = new TFile (file_data_name);
    TH1F * histo_data_obs  = (TH1F*) file_data->Get(hTitle_atanh.Data());
    histo_data_obs->Scale(scaleLumi);
    histo_data_obs->SetName((hTitle_atanh+"_mu_data_obs").Data());
    float eventData = histo_data_obs->Integral(0,histo_data_obs->GetNbinsX()+1);
    int binNumber = histo_data_obs->GetNbinsX();
/////////////////////// END DATA /////////////////////////////////////
    
    
    
    file_tag[file_tag.size()-1] = "DYJetsToLL";



    TFile * file_FitHisto= new TFile ((discriminator+"Systematics.root").c_str(), "recreate");
    file_FitHisto->cd();


    

    
    
    
    for (int j=0;j<uncVariation[0].size();j++){
        for (int i=0;i<nfiles;i++){
        
        std::cout <<  j << "\t" << i << std::endl;
        
        TFile * file         = new TFile (file_names[j][i]);
//         TFile * file_QCDup   = new TFile (file_names_QCDup[i]);
//         TFile * file_QCDdown = new TFile (file_names_QCDdown[i]);
//         TFile * file_JESup   = new TFile (file_names_JESup[i]);
//         TFile * file_JESdown = new TFile (file_names_JESdown[i]);
//         TFile * file_JERup   = new TFile (file_names_JERup[i]);
//         TFile * file_JERdown = new TFile (file_names_JERdown[i]);

        
        
//         if(file_tag[i].CompareTo("DYJetsToLL")==0) {lumi = lumi * k_factor_nom; cout << "filetag test:  " << file_tag[i] << endl;}  // DY_highMass_NLO
//         else lumi = FIXED_lumi;
        lumi = FIXED_lumi;
            
        
//-------------------------------------------------------------------------------AGGIUNGERE QUESTI ERRORI-----------------------------------------------------------------------------

        if (j==0) {

            
            
            TH1F * histo_BDT_nom = (TH1F*) file->Get(hTitle.Data())->Clone();
            TH1F * histo_atanhBDT_nom = (TH1F*) file->Get(hTitle_atanh.Data())->Clone();
            histo_BDT_nom->SetName((histTitle+"mu_"+file_tag[i]).Data());
            histo_atanhBDT_nom->SetName((histTitle+"atanh_mu_"+file_tag[i]).Data());

    

            TH1F * histoSys_BDT_Up_      = (TH1F*) file->Get(hTitle.Data())->Clone((histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data());
            TH1F * histoSys_atanhBDT_Up_ = (TH1F*) file->Get(hTitle_atanh.Data())->Clone((histTitle+"atanh_mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data()); 
            

            std::cout << file_names[j][i] << " \t " << lumi << std::endl;
            std::cout << file_names[j][i] << " \t " << histo_BDT_nom->Integral(0,histo_BDT_nom->GetNbinsX()+1)*FIXED_lumi << std::endl;
            
            
            histoSys_BDT_Up_->Scale(lumi);
            histoSys_atanhBDT_Up_->Scale(lumi);
            histo_BDT_nom->Scale(lumi);
            histo_atanhBDT_nom->Scale(lumi);
            
            eventMC[j][i]= histo_BDT_nom->Integral(0,histo_BDT_nom->GetNbinsX()+1);    
            

        
            if (i == nfiles-1) {
                std::vector<float> normalization_vector;
                for (int n = 0; n < nfiles; n++) normalization_vector.push_back(eventMC[j][n]);
                float DY_realNormalization = compute_DY_normalization(eventData, normalization_vector);
                
                float kFact = DY_realNormalization/normalization_vector[normalization_vector.size()-1];
                std::cout << "k-factor \t " << eventData << "\t" << DY_realNormalization <<  "\t" << normalization_vector[normalization_vector.size()-1]  << " \t " << kFact << std::endl;
                
                
                histoSys_BDT_Up_->Scale(kFact);
                histoSys_atanhBDT_Up_->Scale(kFact);
                histo_BDT_nom->Scale(kFact);
                histo_atanhBDT_nom->Scale(kFact);
                
                lumi = kFact * FIXED_lumi;
                eventMC[j][i]=kFact * eventMC[j][i];
            }
            
            
            file_FitHisto->cd();
            
            TString histTitleLinear = hTitle;
            writeStatisticalErrors (histTitleLinear, file_tag[i], file_FitHisto, file, lumi);
            histTitleLinear = hTitle_atanh;
            writeStatisticalErrors (histTitleLinear, file_tag[i], file_FitHisto, file, lumi);
            
            std::cout << " Prima di writePDFErrors " << std::endl;
            if (file_tag[i].CompareTo("VBF_HToMuMu")==0) writePDFErrors (hTitle_atanh, "VBF_HToMuMu", file_FitHisto, file, lumi);
            std::cout << " Dopo writePDFErrors "  << std::endl;
            
            histoSys_BDT_Up_->Write();
            histoSys_atanhBDT_Up_->Write();
            histo_BDT_nom->Write();
            histo_atanhBDT_nom->Write();  

            
        }
        else {
        
            
            TH1F * histoSys_get           = (TH1F*) file->Get(hTitle.Data());
            TH1F * histoSys_atanh_get     = (TH1F*) file->Get(hTitle_atanh.Data());
            
            TH1F * histoSys               = (TH1F*) histoSys_get->Clone(histoSys_get->GetName());
            TH1F * histoSys_atanh         = (TH1F*) histoSys_atanh_get->Clone(histoSys_atanh_get->GetName());
            
            
            if (j==1) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Up").Data());
            if (j==1) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Up").Data());
            if (j==2) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Down").Data());
            if (j==2) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Down").Data());
            
            if (j==3) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Up").Data());
            if (j==3) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Up").Data());
            if (j==4) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Down").Data());
            if (j==4) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Down").Data());
            
            if (j==5) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Up").Data());
            if (j==5) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Up").Data());
            if (j==6) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Down").Data());
            if (j==6) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Down").Data());
            
            if (j==7) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_PU_Up").Data());
            if (j==7) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_PU_Up").Data());
            if (j==8) histoSys->SetName((histTitle+"mu_"+file_tag[i]+"_PU_Down").Data());
            if (j==8) histoSys_atanh->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_PU_Down").Data());
            

            
                  

            
            
            histoSys->Scale(lumi);
            histoSys_atanh->Scale(lumi);
            eventMC[j][i]= histoSys->Integral(0,histoSys->GetNbinsX()+1);  
            
/*            
            if (i == nfiles-1) {
                float sumMC = 0;
                for (int n=0;n<nfiles;n++) {
                    std::cout << j << " \t "<< n << " \t " << eventMC[j][n] << std::endl;
                    sumMC = sumMC + eventMC[j][n];
                }
                sumMC = sumMC - eventMC[j][nfiles-1];
                float kFact = (eventData - sumMC)/eventMC[j][nfiles-1];
                std::cout << "k-factor \t " << eventData << "\t" << sumMC <<  "\t" << eventMC[j][nfiles-1] << " \t " << (eventData - sumMC) << " \t " << kFact << std::endl;
                
                histoSys->Scale(kFact);
                histoSys_atanh->Scale(kFact);
            }*/
            
            
            
            if (i == nfiles-1) {
                std::vector<float> normalization_vector;
                for (int n = 0; n < nfiles; n++) normalization_vector.push_back(eventMC[j][n]);
                float DY_realNormalization = compute_DY_normalization(eventData, normalization_vector);
                
                float kFact = DY_realNormalization/normalization_vector[normalization_vector.size()-1];
                std::cout << "k-factor \t " << eventData << "\t" << DY_realNormalization <<  "\t" << normalization_vector[normalization_vector.size()-1]  << " \t " << kFact << std::endl;
                
                histoSys->Scale(kFact);
                histoSys_atanh->Scale(kFact);
                
                eventMC[j][i]=kFact * eventMC[j][i];
            } 
            
            
            
            
            
            file_FitHisto->cd();
            
            
            
            histoSys->Write();
            histoSys_atanh->Write();
        
        
        
        
//         TH1F * histoSys_BDT_QCDup_get           = (TH1F*) file_QCDup->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_QCDup_get      = (TH1F*) file_QCDup->Get(hTitle_atanh.Data());
//         TH1F * histoSys_BDT_QCDdown_get         = (TH1F*) file_QCDdown->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_QCDdown_get    = (TH1F*) file_QCDdown->Get(hTitle_atanh.Data());
//         TH1F * histoSys_BDT_JESup_get           = (TH1F*) file_JESup->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_JESup_get      = (TH1F*) file_JESup->Get(hTitle_atanh.Data());
//         TH1F * histoSys_BDT_JESdown_get         = (TH1F*) file_JESdown->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_JESdown_get    = (TH1F*) file_JESdown->Get(hTitle_atanh.Data());
// //         TH1F * histoSys_BDT_JERup_get           = (TH1F*) file_JERup->Get(hTitle.Data());
// //         TH1F * histoSys_atanhBDT_JERup_get      = (TH1F*) file_JERup->Get(hTitle_atanh.Data());
// //         TH1F * histoSys_BDT_JERdown_get         = (TH1F*) file_JERdown->Get(hTitle.Data());
// //         TH1F * histoSys_atanhBDT_JERdown_get    = (TH1F*) file_JERdown->Get(hTitle_atanh.Data());


                
//         TH1F * histoSys_BDT_QCDup           = (TH1F*) histoSys_BDT_QCDup_get->Clone(histoSys_BDT_QCDup_get->GetName());
//         TH1F * histoSys_atanhBDT_QCDup      = (TH1F*) histoSys_atanhBDT_QCDup_get->Clone(histoSys_atanhBDT_QCDup_get->GetName());
//         TH1F * histoSys_BDT_QCDdown         = (TH1F*) histoSys_BDT_QCDdown_get->Clone(histoSys_BDT_QCDdown_get->GetName());
//         TH1F * histoSys_atanhBDT_QCDdown    = (TH1F*) histoSys_atanhBDT_QCDdown_get->Clone(histoSys_atanhBDT_QCDdown_get->GetName());
//         TH1F * histoSys_BDT_JESup           = (TH1F*) histoSys_BDT_JESup_get->Clone(histoSys_BDT_JESup_get->GetName());
//         TH1F * histoSys_atanhBDT_JESup      = (TH1F*) histoSys_atanhBDT_JESup_get->Clone(histoSys_atanhBDT_JESup_get->GetName());
//         TH1F * histoSys_BDT_JESdown         = (TH1F*) histoSys_BDT_JESdown_get->Clone(histoSys_BDT_JESdown_get->GetName());
//         TH1F * histoSys_atanhBDT_JESdown    = (TH1F*) histoSys_atanhBDT_JESdown_get->Clone(histoSys_atanhBDT_JESdown_get->GetName());
// //         TH1F * histoSys_BDT_JERup           = (TH1F*) histoSys_BDT_JERup_get->Clone(histoSys_BDT_JERup_get->GetName());
// //         TH1F * histoSys_atanhBDT_JERup      = (TH1F*) histoSys_atanhBDT_JERup_get->Clone(histoSys_atanhBDT_JERup_get->GetName());
// //         TH1F * histoSys_BDT_JERdown         = (TH1F*) histoSys_BDT_JERdown_get->Clone(histoSys_BDT_JERdown_get->GetName());
// //         TH1F * histoSys_atanhBDT_JERdown    = (TH1F*) histoSys_atanhBDT_JERdown_get->Clone(histoSys_atanhBDT_JERdown_get->GetName());



                
//         histoSys_BDT_QCDup->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Up").Data());
//         histoSys_atanhBDT_QCDup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Up").Data());
//         histoSys_BDT_QCDdown->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Down").Data());
//         histoSys_atanhBDT_QCDdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Down").Data());
//         histoSys_BDT_JESup->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Up").Data());
//         histoSys_atanhBDT_JESup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Up").Data());
//         histoSys_BDT_JESdown->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Down").Data());
//         histoSys_atanhBDT_JESdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Down").Data());
// //         histoSys_BDT_JERup->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Up").Data());
// //         histoSys_atanhBDT_JERup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Up").Data());
// //         histoSys_BDT_JERdown->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Down").Data());
// //         histoSys_atanhBDT_JERdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Down").Data());

      
//         std::cout << histoSys_BDT_JESup << "      "+histTitle+"atanh_mu_"+file_tag[i]+"_JER_Down \t" << histoSys_BDT_JESup->GetEntries() << std::endl;
        
    
        
//         histoSys_BDT_QCDup->Scale(lumi * k_factor_QCDup);
//         histoSys_atanhBDT_QCDup->Scale(lumi * k_factor_QCDup);
//         histoSys_BDT_QCDdown->Scale(lumi * k_factor_QCDdown);
//         histoSys_atanhBDT_QCDdown->Scale(lumi * k_factor_QCDdown);
//         histoSys_BDT_JESup->Scale(lumi * k_factor_JESup);
//         histoSys_atanhBDT_JESup->Scale(lumi * k_factor_JESup);
//         histoSys_BDT_JESdown->Scale(lumi * k_factor_JESdown);
//         histoSys_atanhBDT_JESdown->Scale(lumi * k_factor_JESdown);
// //         histoSys_BDT_JERup->Scale(lumi);
// //         histoSys_atanhBDT_JERup->Scale(lumi);
// //         histoSys_BDT_JERdown->Scale(lumi);
// //         histoSys_atanhBDT_JERdown->Scale(lumi);
                
//         histo_BDT_nom->Scale(lumi);;
//         histo_atanhBDT_nom->Scale(lumi); 
        
        
//         histoSys_BDT_Up_->Scale(lumi);
//         histoSys_atanhBDT_Up_->Scale(lumi);



        
        
        
//         file_FitHisto->cd();
// 
//         
//         histoSys_BDT_QCDup->Write();
//         histoSys_atanhBDT_QCDup->Write();
//         histoSys_BDT_QCDdown->Write();
//         histoSys_atanhBDT_QCDdown->Write();
//         histoSys_BDT_JESup->Write();
//         histoSys_atanhBDT_JESup->Write();
//         histoSys_BDT_JESdown->Write();
//         histoSys_atanhBDT_JESdown->Write();
// //         histoSys_BDT_JERup->Write();
// //         histoSys_atanhBDT_JERup->Write();
// //         histoSys_BDT_JERdown->Write();
// //         histoSys_atanhBDT_JERdown->Write();

            }
        
        }
    }
    
    


    file_FitHisto->cd();
    histo_data_obs->Write();
    
    
    
    
    file_FitHisto->Close();
    
    
    ofstream eventMC_txt;
    eventMC_txt.open("eventMC.txt"); 
    for (int j=0;j<uncVariation[0].size()-1;j++){
        for (int i=0;i<nfiles;i++){
            eventMC_txt<< file_names[j][i] << "  \t\t\t "<< eventMC[j][i] <<endl;    
            std::cout << j <<  " \t"  << i << "  \t\t\t "<< eventMC[j][i] <<endl;       
        }
        eventMC_txt <<endl; 
        std::cout << endl;

    }
    eventMC_txt.close();
    
    dataCard(eventMC, file_tag, hTitle_atanh, discriminator, eventData, binNumber);
    
}






        
        
        
        
