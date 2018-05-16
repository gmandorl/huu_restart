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




void addErrorToHist (TH1F * histo, float sigma, int n) {
        histo->SetBinContent(n, histo->GetBinContent(n) + histo->GetBinError(n)*sigma);
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


void macroFitHisto () {
    std::cout << "I can see here" << std::endl;

    float lumi = 35900.;
    
    ////////////// DY highMass NLO ////////////////////////////
//     float k_factor_nom    = 0.98280948;
// 
//         
// //     float k_factor_QCDup    = 1.0399461;
// //     float k_factor_QCDdown  = 0.90371829;
// //     float k_factor_JESup    = 0.78923333;
// //     float k_factor_JESdown  = 1.1755885;
//     
//     float k_factor_QCDup    = 1.0399461     / k_factor_nom;    // I have to divide by k_factor_nom because DY lumi is multiplied by k_factor_nom
//     float k_factor_QCDdown  = 0.90371829    / k_factor_nom;
//     float k_factor_JESup    = 0.78923333    / k_factor_nom;
//     float k_factor_JESdown  = 1.1755885     / k_factor_nom;
    
    
    
    //      NN 
//         float k_factor_nom      = 0.97476554;
//         float k_factor_QCDup    = 1.0604774     / k_factor_nom;    // I have to divide by k_factor_nom because DY lumi is multiplied by k_factor_nom
//         float k_factor_QCDdown  = 0.91036564    / k_factor_nom;
//         float k_factor_JESup    = 0.81818968    / k_factor_nom;
//         float k_factor_JESdown  = 1.1656951     / k_factor_nom;
    
    
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

    TString DY_name = "DYInclusivetoLL";
    TString DY_file_name = "DYInclusivetoLL_M";

    DY_name = "DYJetstoLL_amc";
    DY_file_name = "DYJetstoLL_amc";

    //DY_name = "DYJetstoLL_amc_Pt";
    //DY_file_name = "DYJetstoLL_amc_Pt";


    // DY_name = "DYJetstoLL_Final_amc";
    // DY_file_name = "DYJetstoLL_Final_amc";


    const int nfiles  = 9;
    // TString file_names[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
    TString file_names[nfiles] = {"WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
    TString file_names_QCDup[nfiles];
    TString file_names_QCDdown[nfiles];
    TString file_names_JESup[nfiles];
    TString file_names_JESdown[nfiles];
    TString file_names_JERup[nfiles];
    TString file_names_JERdown[nfiles];
    TString file_data_name = "SingleMuon";
//     TString file_data_name = "WW";
    
    TString file_tag[nfiles];

    //TString file_names_QCDup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
    //TString file_names_QCDdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
    //TString file_names_JESup[nfiles] =  {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};
    //TString file_names_JESdown[nfiles] ={"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT",DY_file_name,"VBF_HToMuMu", "GluGlu_HToMuMu"};


    for (int i=0;i<nfiles;i++){
        file_tag[i] = file_names[i];
        
        
 
        file_names[i].Append("_mu");
        file_names_QCDup[i]   = file_names[i];
        file_names_QCDdown[i] = file_names[i];
        file_names_JESup[i]   = file_names[i];
        file_names_JESdown[i] = file_names[i];
        file_names_JERup[i]   = file_names[i];
        file_names_JERdown[i] = file_names[i];
        
        file_names[i].Append("_QCDScalenom_JESnom_v25_reskim.root");
        file_names_QCDup[i].Append("_QCDScaleup_JESnom_v25_reskim.root");
        file_names_QCDdown[i].Append("_QCDScaledown_JESnom_v25_reskim.root");
        file_names_JESup[i].Append("_QCDScalenom_JESup_v25_reskim.root");
        file_names_JESdown[i].Append("_QCDScalenom_JESdown_v25_reskim.root");
        file_names_JERup[i].Append("_QCDScaleup_JESup_v25_reskim.root");
        file_names_JERdown[i].Append("_QCDScaledown_JESdown_v25_reskim.root");
        
    }
    
    
    file_data_name.Append("_mu");
    file_data_name.Append("_QCDScalenom_JESnom_v25_reskim.root");

    file_tag[6] = "DYJetsToLL";



//    TFile * fileDY_Jet   = new TFile (fileDY_name_Jet.c_str());
//    TFile * fileDY_Pt    = new TFile (fileDY_name_Pt.c_str());
    TFile * file_FitHisto= new TFile ("BDTSystematics.root", "recreate");
    file_FitHisto->cd();


    
//         TString histTitle = "hBDT_VBF_";
//         TString hTitle          = "hBDT_VBF";
//         TString hTitle_atanh    = "hBDT_VBF_atanh";
        
        TString histTitle = "hNNoutput_";
        TString hTitle          = "hNNoutput";
        TString hTitle_atanh    = "hNNoutput_atanh";

    for (int i=0;i<nfiles;i++){
//         if (i==0) Nsyst_NoConst = 1;
//         else Nsyst_NoConst = Nsyst;


        
        TFile * file         = new TFile (file_names[i]);
        TFile * file_QCDup   = new TFile (file_names_QCDup[i]);
        TFile * file_QCDdown = new TFile (file_names_QCDdown[i]);
        TFile * file_JESup   = new TFile (file_names_JESup[i]);
        TFile * file_JESdown = new TFile (file_names_JESdown[i]);
//         TFile * file_JERup   = new TFile (file_names_JERup[i]);
//         TFile * file_JERdown = new TFile (file_names_JERdown[i]);

        
        
        if(file_tag[i].CompareTo("DYJetsToLL")==0) {lumi = lumi * k_factor_nom; cout << "filetag test:  " << file_tag[i] << endl;}  // DY_highMass_NLO
//         if(file_tag[i].CompareTo(DY_file_name)==0) {lumi = lumi * 0.997; cout << "filetag test:  " << file_tag[i] << endl;}       // DY jet-multiplicity NLO
        else lumi = 35900.;
            
            
        
//-------------------------------------------------------------------------------AGGIUNGERE QUESTI ERRORI-----------------------------------------------------------------------------

        
        for (int j=1;j<Nsyst_NoConst;++j){

    
            cout << "filename = " << file_names[i] << "   \t" <<  histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[j]+"_Up"<< endl;
 
            TH1F * histoSys_BDT_Up          = (TH1F*) file->Get(hTitle.Data())->Clone((histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[j]+"_Up").Data());
            TH1F * histoSys_atanhBDT_Up     = (TH1F*) file->Get(hTitle_atanh.Data())->Clone((histTitle+"atanh_mu_"+file_tag[i]+"_"+uncertainty_name[j]+"_Up").Data());
            TH1F * histoSys_BDT_Down        = (TH1F*) file->Get(hTitle.Data())->Clone((histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[j]+"_Down").Data());
            TH1F * histoSys_atanhBDT_Down   = (TH1F*) file->Get(hTitle_atanh.Data())->Clone((histTitle+"atanh_mu_"+file_tag[i]+"_"+uncertainty_name[j]+"_Down").Data());
            

                
                
            file_FitHisto->cd();

            histoSys_BDT_Up->Scale(lumi);
            histoSys_atanhBDT_Up->Scale(lumi);
            histoSys_BDT_Down->Scale(lumi);
            histoSys_atanhBDT_Down->Scale(lumi);
            
            histoSys_BDT_Up->Write();
            histoSys_atanhBDT_Up->Write();
            histoSys_BDT_Down->Write();
            histoSys_atanhBDT_Down->Write();

        }
        
        
        
        
        TH1F * histo_BDT_nom = (TH1F*) file->Get(hTitle.Data())->Clone();
        TH1F * histo_atanhBDT_nom = (TH1F*) file->Get(hTitle_atanh.Data())->Clone();
        histo_BDT_nom->SetName((histTitle+"mu_"+file_tag[i]).Data());
        histo_atanhBDT_nom->SetName((histTitle+"atanh_mu_"+file_tag[i]).Data());
        
        TString histTitleLinear = hTitle;
        writeStatisticalErrors (histTitleLinear, file_tag[i], file_FitHisto, file, lumi);
        histTitleLinear = hTitle_atanh;
        writeStatisticalErrors (histTitleLinear, file_tag[i], file_FitHisto, file, lumi);
        

        
        
        
        
/*        TH1F * histoSys_BDT_Up_ = (TH1F*) file->Get()->Clone((histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data());
        TH1F * histoSys_atanhBDT_Up_ = (TH1F*) file->Get()->Clone((histTitle+"atanh_mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data());  */      

        TH1F * histoSys_BDT_Up_      = (TH1F*) file->Get(hTitle.Data())->Clone((histTitle+"mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data());
        TH1F * histoSys_atanhBDT_Up_ = (TH1F*) file->Get(hTitle_atanh.Data())->Clone((histTitle+"atanh_mu_"+file_tag[i]+"_"+uncertainty_name[0]+"_Up").Data());       

        
        
        TH1F * histoSys_BDT_QCDup_get           = (TH1F*) file_QCDup->Get(hTitle.Data());
        TH1F * histoSys_atanhBDT_QCDup_get      = (TH1F*) file_QCDup->Get(hTitle_atanh.Data());
        TH1F * histoSys_BDT_QCDdown_get         = (TH1F*) file_QCDdown->Get(hTitle.Data());
        TH1F * histoSys_atanhBDT_QCDdown_get    = (TH1F*) file_QCDdown->Get(hTitle_atanh.Data());
        TH1F * histoSys_BDT_JESup_get           = (TH1F*) file_JESup->Get(hTitle.Data());
        TH1F * histoSys_atanhBDT_JESup_get      = (TH1F*) file_JESup->Get(hTitle_atanh.Data());
        TH1F * histoSys_BDT_JESdown_get         = (TH1F*) file_JESdown->Get(hTitle.Data());
        TH1F * histoSys_atanhBDT_JESdown_get    = (TH1F*) file_JESdown->Get(hTitle_atanh.Data());
//         TH1F * histoSys_BDT_JERup_get           = (TH1F*) file_JERup->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_JERup_get      = (TH1F*) file_JERup->Get(hTitle_atanh.Data());
//         TH1F * histoSys_BDT_JERdown_get         = (TH1F*) file_JERdown->Get(hTitle.Data());
//         TH1F * histoSys_atanhBDT_JERdown_get    = (TH1F*) file_JERdown->Get(hTitle_atanh.Data());


                
        TH1F * histoSys_BDT_QCDup           = (TH1F*) histoSys_BDT_QCDup_get->Clone(histoSys_BDT_QCDup_get->GetName());
        TH1F * histoSys_atanhBDT_QCDup      = (TH1F*) histoSys_atanhBDT_QCDup_get->Clone(histoSys_atanhBDT_QCDup_get->GetName());
        TH1F * histoSys_BDT_QCDdown         = (TH1F*) histoSys_BDT_QCDdown_get->Clone(histoSys_BDT_QCDdown_get->GetName());
        TH1F * histoSys_atanhBDT_QCDdown    = (TH1F*) histoSys_atanhBDT_QCDdown_get->Clone(histoSys_atanhBDT_QCDdown_get->GetName());
        TH1F * histoSys_BDT_JESup           = (TH1F*) histoSys_BDT_JESup_get->Clone(histoSys_BDT_JESup_get->GetName());
        TH1F * histoSys_atanhBDT_JESup      = (TH1F*) histoSys_atanhBDT_JESup_get->Clone(histoSys_atanhBDT_JESup_get->GetName());
        TH1F * histoSys_BDT_JESdown         = (TH1F*) histoSys_BDT_JESdown_get->Clone(histoSys_BDT_JESdown_get->GetName());
        TH1F * histoSys_atanhBDT_JESdown    = (TH1F*) histoSys_atanhBDT_JESdown_get->Clone(histoSys_atanhBDT_JESdown_get->GetName());
//         TH1F * histoSys_BDT_JERup           = (TH1F*) histoSys_BDT_JERup_get->Clone(histoSys_BDT_JERup_get->GetName());
//         TH1F * histoSys_atanhBDT_JERup      = (TH1F*) histoSys_atanhBDT_JERup_get->Clone(histoSys_atanhBDT_JERup_get->GetName());
//         TH1F * histoSys_BDT_JERdown         = (TH1F*) histoSys_BDT_JERdown_get->Clone(histoSys_BDT_JERdown_get->GetName());
//         TH1F * histoSys_atanhBDT_JERdown    = (TH1F*) histoSys_atanhBDT_JERdown_get->Clone(histoSys_atanhBDT_JERdown_get->GetName());



                
        histoSys_BDT_QCDup->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Up").Data());
        histoSys_atanhBDT_QCDup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Up").Data());
        histoSys_BDT_QCDdown->SetName((histTitle+"mu_"+file_tag[i]+"_QCD_Down").Data());
        histoSys_atanhBDT_QCDdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_QCD_Down").Data());
        histoSys_BDT_JESup->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Up").Data());
        histoSys_atanhBDT_JESup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Up").Data());
        histoSys_BDT_JESdown->SetName((histTitle+"mu_"+file_tag[i]+"_JES_Down").Data());
        histoSys_atanhBDT_JESdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JES_Down").Data());
//         histoSys_BDT_JERup->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Up").Data());
//         histoSys_atanhBDT_JERup->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Up").Data());
//         histoSys_BDT_JERdown->SetName((histTitle+"mu_"+file_tag[i]+"_JER_Down").Data());
//         histoSys_atanhBDT_JERdown->SetName((histTitle+"atanh_mu_"+file_tag[i]+"_JER_Down").Data());

      
        std::cout << histoSys_BDT_JESup << "      "+histTitle+"atanh_mu_"+file_tag[i]+"_JER_Down \t" << histoSys_BDT_JESup->GetEntries() << std::endl;
        
    
        
        histoSys_BDT_QCDup->Scale(lumi * k_factor_QCDup);
        histoSys_atanhBDT_QCDup->Scale(lumi * k_factor_QCDup);
        histoSys_BDT_QCDdown->Scale(lumi * k_factor_QCDdown);
        histoSys_atanhBDT_QCDdown->Scale(lumi * k_factor_QCDdown);
        histoSys_BDT_JESup->Scale(lumi * k_factor_JESup);
        histoSys_atanhBDT_JESup->Scale(lumi * k_factor_JESup);
        histoSys_BDT_JESdown->Scale(lumi * k_factor_JESdown);
        histoSys_atanhBDT_JESdown->Scale(lumi * k_factor_JESdown);
//         histoSys_BDT_JERup->Scale(lumi);
//         histoSys_atanhBDT_JERup->Scale(lumi);
//         histoSys_BDT_JERdown->Scale(lumi);
//         histoSys_atanhBDT_JERdown->Scale(lumi);
                
        histo_BDT_nom->Scale(lumi);;
        histo_atanhBDT_nom->Scale(lumi); 
        
        
//         histoSys_BDT_Up_->Scale(lumi);
//         histoSys_atanhBDT_Up_->Scale(lumi);



        
        
        
        file_FitHisto->cd();

        
        histoSys_BDT_QCDup->Write();
        histoSys_atanhBDT_QCDup->Write();
        histoSys_BDT_QCDdown->Write();
        histoSys_atanhBDT_QCDdown->Write();
        histoSys_BDT_JESup->Write();
        histoSys_atanhBDT_JESup->Write();
        histoSys_BDT_JESdown->Write();
        histoSys_atanhBDT_JESdown->Write();
//         histoSys_BDT_JERup->Write();
//         histoSys_atanhBDT_JERup->Write();
//         histoSys_BDT_JERdown->Write();
//         histoSys_atanhBDT_JERdown->Write();

        
        histoSys_BDT_Up_->Write();
        histoSys_atanhBDT_Up_->Write();


        histo_BDT_nom->Write();
        histo_atanhBDT_nom->Write();           
    }
    
    
    TFile * file_data   = new TFile (file_data_name);
    TH1F * histo_data_obs  = (TH1F*) file_data->Get(hTitle_atanh.Data());
    histo_data_obs->SetName((hTitle_atanh+"_mu_data_obs").Data());

    file_FitHisto->cd();
    histo_data_obs->Write();
    
    
    
    
    file_FitHisto->Close();
}






        
        
        
        
