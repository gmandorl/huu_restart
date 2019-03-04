#include <ctype.h>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1D.h"
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
#include <TEfficiency.h>
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
#include <random>
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.h"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.h"
#include "EWcorr.C"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/rochcor2016.h"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/rochcor2016.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/RoccoR.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/RoccoR.h"
// #include "/afs/cern.ch/user/g/gimandor/private/CMSSW_8_0_24/src/giulioMandorli/2016/RoccoR.cc"
#include "2016/RoccoR.cc"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

// #include "parse_json.hh"
// #include "LightweightNeuralNetwork.hh"
// #include "NNet_comparison.h"


#include "VBFFilter.h"



#define getName(var)  #var


// #include "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/ratioDYFilter.h"


//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
//#include "NNet_comparison.h"

//#include "/afs/pi.infn.it/user/mandorli/mandorli/Hmumu/CMSSW_8_0_25/src/mva/giulioMandorli/LWTNN/interface/LightweightNeuralNetwork.h"
//#include "/afs/pi.infn.it/user/mandorli/mandorli/Hmumu/CMSSW_8_0_25/src/mva/giulioMandorli/LWTNN/interface/parse_json.h"  


Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t erf2( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]))+ (1.-par[0]);
}


TString Correct_fileTag (TString file_tag_hBDT_name) {
    TString fileTag_names[10] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetsToLL","VBF_HToMuMu", "GluGlu_HToMuMu"};
    TString fileTag_toReturn = file_tag_hBDT_name;

    if ((file_tag_hBDT_name.CompareTo("ST_tW_top") == 0) || (file_tag_hBDT_name.CompareTo("ST_tW_antitop") == 0) || (file_tag_hBDT_name.CompareTo("ST_s-channel") == 0) || (file_tag_hBDT_name.CompareTo("ST_t-channel_top_4f_inclusiveDecays") == 0) || (file_tag_hBDT_name.CompareTo("ST_t-channel_antitop_4f_inclusiveDecays") == 0)) fileTag_toReturn = fileTag_names[5];

    if ((file_tag_hBDT_name.CompareTo("DYJetstoLL_amc_0J") == 0) || (file_tag_hBDT_name.CompareTo("DYJetstoLL_amc_1J") == 0) || (file_tag_hBDT_name.CompareTo("DYJetstoLL_amc_2J") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-0To50_amc") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-50To100_amc") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-100To250_amc") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-250To400_amc") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-400To650_amc") == 0) || (file_tag_hBDT_name.CompareTo("DYJetsToLL_Pt-650ToInf_amc") == 0)) fileTag_toReturn = fileTag_names[7];


    return fileTag_toReturn;
}


#define SWAP2(A, B) { TLorentzVector t = A; A = B; B = t; }
void SortByEta(std::vector<TLorentzVector> &jets){
  int i, j;
	int n=jets.size();
  for (i = n - 1; i >= 0; i--){
    for (j = 0; j < i; j++){
      if (jets[j].Eta() < jets[j + 1].Eta() ){
        SWAP2( jets[j], jets[j + 1] );
		}
    }
	}
}

float getScaleFactor(TH2F *scaleMap, double pt, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
	 int biny;
    if (abs==0) biny = scaleMap->GetYaxis()->FindBin(eta);
    else biny = scaleMap->GetYaxis()->FindBin(TMath::Abs(eta));
  //  std::cout<<binx<<": ,"<<biny<<std::endl;
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) && (biny != 0) && (biny != scaleMap->GetNbinsY()+1)) {
        sfactor = scaleMap->GetBinContent(binx, biny);
        sf_err = scaleMap->GetBinError(binx, biny);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}
float getScaleFactor1D(TH1F *scaleMap, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
	 int binx;
    if (abs==0) binx = scaleMap->GetXaxis()->FindBin(eta);
    else binx = scaleMap->GetXaxis()->FindBin(TMath::Abs(eta));
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) ) {
        sfactor = scaleMap->GetBinContent(binx);
        sf_err = scaleMap->GetBinError(binx);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}



double SumEntries(TChain * tree_Runs, std::string countVariableName, bool readArray, int LHEelement) {
    double sum = 0;
    std::cout << "QCD VARIATION  " << LHEelement << "  ----------------------" << std::endl;
//     std::cout << "tree_Runs  " << tree_Runs << "   "  << readArray << std::endl;
    if (readArray) {
        Double_t x[9]; 
        double genEventSumw; 
        tree_Runs->SetBranchAddress("genEventSumw",&genEventSumw);
        tree_Runs->SetBranchAddress(countVariableName.c_str(),x);
        for (int n_file = 0; n_file <  tree_Runs->GetEntries(); n_file++) {
            tree_Runs->GetEntry(n_file);
//             std::cout << "variation n_file  " << x[LHEelement] * genEventSumw << "  \t sum   "  << sum  << std::endl;
            sum += x[LHEelement]/x[4] * genEventSumw;
        }  
        if (sum < 1.) return 1.;
        else return sum;//tree_Runs->GetEntries(); 
    }
    
    else {
        double x = 0;
        tree_Runs->SetBranchAddress(countVariableName.c_str(),&x);
//         std::cout << "tree_Runs  " << tree_Runs << "   "  << readArray << std::endl;
        for (int n_file = 0; n_file <  tree_Runs->GetEntries(); n_file++) {
            
            tree_Runs->GetEntry(n_file);
//             std::cout << "n_file  " << x << "\t\t sum   " << sum << std::endl;
            sum += x;
        }
        return sum;
    }
    
}    


                
float findMean (TH1D * h) {
    float mean = 0;
    for (int n=1; n <= h->GetNbinsX(); n++) mean += h->GetBinContent(n);
    return mean/h->GetNbinsX();
}

float findRMS (TH1D * h) {
    float RMS = 0;
    float mean = findMean(h);
    for (int n=1; n <= h->GetNbinsX(); n++) RMS += (h->GetBinContent(n) - mean) * (h->GetBinContent(n) - mean);
    RMS = sqrt( RMS/(h->GetNbinsX()-1) );
    return RMS;
}




void write_PDF_variation(TH1D * hBDT_VBF_atanh, TH1D * hBDT_VBF_atanh_PDFvariation, int n) {
    
    std::string name = hBDT_VBF_atanh->GetName();
    TH1D * hBDT_VBF_atanh_Up   = (TH1D*) hBDT_VBF_atanh->Clone((name+"_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Up").c_str());
    TH1D * hBDT_VBF_atanh_Down = (TH1D*) hBDT_VBF_atanh->Clone((name+"_mu_VBF_HToMuMu_PDF"+std::to_string(n)+"_Down").c_str());
    
    
    float RMS = findRMS(hBDT_VBF_atanh_PDFvariation);
//     std::cout << "bin value     " << hBDT_VBF_atanh_copy->GetBinContent(n)<< " \t mean   " << findMean(hBDT_VBF_atanh_PDFvariation)  << " \t RMS   " << RMS << std::endl;
    hBDT_VBF_atanh_Up->SetBinContent  (n, hBDT_VBF_atanh_Up->GetBinContent(n)   + RMS);
    hBDT_VBF_atanh_Down->SetBinContent(n, hBDT_VBF_atanh_Down->GetBinContent(n) - RMS);
    
//     hBDT_VBF_atanh_Up->Write();              // they are automaticaly saved becaused they are clone of something already saved
//     hBDT_VBF_atanh_Down->Write();
}

                        
                        
                        

float computeThetaStar (TLorentzVector p1, TLorentzVector p2) {
    
        TLorentzVector pSum = p1 + p2;
        TVector3 BoostVector_toRestFrame = pSum.BoostVector();
        TLorentzVector p1_newSys = p1;

        p1_newSys.Boost(-BoostVector_toRestFrame);
        TVector3 Sum_direction  = pSum.Vect();
        TVector3 p1_new_direction = p1_newSys.Vect();
        return Sum_direction.Dot(p1_new_direction)/Sum_direction.Mag()/p1_new_direction.Mag();
        
}


typedef std::map<double, int> JetList;
const int njets = 120;




typedef struct {
    Float_t eta[njets];
    Float_t pt[njets];
    Float_t JEC_corr[njets];
    Float_t JEC_corr_up[njets];
    Float_t JEC_corr_down[njets];
    Float_t JER_corr[njets];
    Float_t JER_corr_up[njets];
    Float_t JER_corr_down[njets];
    Float_t phi[njets];
	Float_t mass[njets];
	Float_t btagCSV[njets];
    Float_t btagCMVA[njets];
	Int_t nsoft;
	Float_t soft_pt[njets];
	Float_t soft_eta[njets];
	Float_t soft_mass[njets];
	Float_t qgl[njets];
	Int_t nsoft2;
	Int_t nsoft5;
	Int_t nsoft10;
//     Int_t nsoftActivityEWKJets;
    UInt_t  nsoftActivityEWKJets;
    Int_t EWKnsoft2;
	Int_t EWKnsoft5;
	Int_t EWKnsoft10;
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Float_t lepton1_pt;
	Float_t lepton2_pt;
	Float_t jets12;
	Int_t partonFlavour[njets];
    Float_t softLeadingJet_pt;
	Float_t EWKHTsoft;
	Float_t EWKsoft_pt[njets];
	Float_t EWKsoft_eta[njets];
	Float_t EWKsoft_phi[njets];
	Float_t EWKsoft_mass[njets];
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
        
        Float_t VBFselected[njets];
        int muonIdx1[njets];
        int muonIdx2[njets];
        int electronIdx1[njets];
        int electronIdx2[njets];
        int Jet_genJetIdx[njets];
        
        int jetIdx1;
        int jetIdx2;
        
        int photonIdx1[njets];
        int photonIdx2[njets];

        
         Float_t pt_JERup[njets];   
         Float_t mass_JERup[njets];
         Float_t pt_JESup[njets];   
         Float_t mass_JESup[njets];
         Float_t pt_JERdown[njets];   
         Float_t mass_JERdown[njets];
         Float_t pt_JESdown[njets];   
         Float_t mass_JESdown[njets];
         
         Float_t chEmEF[njets];
         Float_t neEmEF[njets];      
         Float_t chHEF[njets]; 
         Float_t neHEF[njets]; 
         
} Jets;


typedef struct {
        
    
        Float_t ll_mass;
        Float_t ll_mass_biasUp;
        Float_t ll_mass_biasDown;
        Float_t ll_mass_resolution;
            
//         Float_t Zll_pt;
   	Float_t Mqq; 
        Float_t MqqLog;
	Float_t DeltaEtaQQ;
  	Float_t q1_eta;
    Float_t met_pt;
    Float_t EWKHTsoft;
    Float_t btagCMVA;
    Float_t btagCMVA_leading;
    Float_t btagCSV;
    Float_t btagCMVA_second;
    Float_t btagCSV_second;
	Float_t softLeadingJet_pt;
    Float_t cosThetaStar;
	Float_t cosThetaStarAbs;
    Float_t cosThetaPlane;
	Float_t absCosThetaStarJet;
                       
    Float_t lepton1_pt;
	Float_t lepton2_pt;
	Float_t jets12;    
	Float_t qq_pt;
	Float_t Jet3_eta;
	Float_t Jet3_pt;
	Float_t axis2_jet1;
	Float_t axis2_jet2;
	Float_t Jet3q_pt; 
	Float_t Jet2q_pt; 
	Float_t Jet2q_leadTrackPt;
	Float_t Jet1q_pt;
        Float_t Jet1q_eta;
        Float_t Jet2q_eta;
        
        Float_t Jet1q_chEmEF;
        Float_t Jet2q_chEmEF;
        Float_t Jet1q_neEmEF;
        Float_t Jet2q_neEmEF;

        Float_t Jet2q_chHEF;
        Float_t Jet1q_chHEF;
        Float_t Jet2q_neHEF;
        Float_t Jet1q_neHEF;

        
        Float_t Jet1q_puId;
        Float_t Jet2q_puId;
        
        
	Float_t Jet1q_leadTrackPt;
	Float_t ll_ystar;
	Float_t ll_zstar;
	Float_t Rpt;
	Float_t ll_pt;
	Float_t ll_eta;
        Float_t ll_pt_Log;
	Float_t qgl_1q;
	Float_t qgl_2q;
        Float_t qgl_1qAtanh;
        Float_t qgl_2qAtanh;
        
        Float_t Jet1q_ptLog;
        Float_t Jet2q_ptLog;
        Float_t Jet3q_ptLog;
        
        Float_t Jet1q_pt_jerUp;
        Float_t Jet1q_pt_jerDown;
        Float_t Jet2q_pt_jerUp;
        Float_t Jet2q_pt_jerDown;
        Float_t genJetPt_match1;
        Float_t genJetPt_match2;
        
        Int_t genJetIdx0;
        Int_t genJetIdx1;
        Int_t genJetIdx2;
        Int_t softActivityEWK_njets2;
        Int_t softActivityEWK_njets5;
	Int_t softActivityEWK_njets10;
        
        Int_t genJetIdx[30];
        Float_t Jet_pt[30];
        
        Int_t genJetIdx0_genJetMassWithoutLepton;
        Int_t genJetIdx1_genJetMassWithoutLepton;
        Int_t genJetIdx2_genJetMassWithoutLepton;
        
        Float_t randomVariable;
        Float_t xSection;
	
	Float_t Inv_mass;
        Float_t Invariant_MassLog;
        Float_t X_parton1;
        Float_t X_parton2;
        Float_t W_mass_virtual1;
        Float_t W_mass_virtual2;
        Float_t W_Pt_virtual1;
        Float_t W_Pt_virtual2;
        Float_t W_eta_virtual1;
        Float_t W_eta_virtual2;
        Float_t W_phi_virtual1;
        Float_t W_phi_virtual2;
        Float_t thetastarW1;
        Float_t thetastarW2;
        Float_t thetastarW2toHW1;
        Float_t thetastarW1toHW2;
        Float_t thetastarHtoWW;
        Float_t theta1;
        Float_t theta2;
        
        Float_t WWmass;
        Float_t impulsoZ;
        Float_t energytot;
        Float_t energytotLog;
        Float_t diffMassWWH;
        Float_t deltaMRel;
        Float_t deltaM;
        Float_t normalizedDistance_from_mH;
        Float_t mumujj_pt;
        
        
                
        Float_t quark1_pt; 
        Float_t quark1_eta;  
        Float_t quark1_phi;
        Float_t quark1_mass;       
        Float_t quark2_pt;
        Float_t quark2_eta;
        Float_t quark2_phi;  
        Float_t quark2_mass;
        
        Float_t genMuon1_pt;
        Float_t genMuon1_eta;
        Float_t genMuon1_phi;
        Float_t genMuon1_mass;      
        Float_t genMuon2_pt;
        Float_t genMuon2_eta;
        Float_t genMuon2_phi;
        Float_t genMuon2_mass;
        
        Float_t PhiZQ1;
        Float_t PhiZQ2;

        Float_t cosThetaPlaneAtanh;
        Float_t absCosThetaStarJetAtanh;
        Float_t X_parton1Log;
        Float_t X_parton2Log;
        Float_t W_mass_virtual1Log;
        Float_t W_mass_virtual2Log;
        Float_t W_Pt_virtual1Log;
        Float_t W_Pt_virtual2Log;
        
        Float_t minAbsEta;
        Float_t maxAbsEta;
        Float_t qqMass_skim;
        Float_t MuonMass_skim;
        Float_t genJetMass_leading;

        Float_t JetPuId_maxAbsEta;
        Int_t countJet25;
        
        Float_t Muon1_relIso04;
        Float_t Muon2_relIso04;
        Float_t LHE_weights_scale_wUp;
        Float_t LHE_weights_scale_wDown;
        
        Float_t BDToutput;
        
        float weightMVA;
        float genweight;	
        float genweight_noQGLcorrection;
        std::vector<float> genweightVECTOR;

        
        
}TMVAstruct;




using namespace std;


int main(int argc, char* argv[]){

std::cout << "I can see here" << std::endl;
//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");

bool andrewSelection = false;
// andrewSelection = true;


bool fromNANO = true;

        UChar_t lheNpNLO;     // if from nanoAOD
        UInt_t nselLeptons; // if from nanoAOD
        UChar_t lheNj_f;   // if from nanoAOD
        UInt_t nJets;         // if from nanoAOD
        Bool_t HLT_IsoMu27;      // if from nanoAOD
	Bool_t HLT_IsoTkMu27;      // if from nanoAOD
	Bool_t HLT_IsoTkMu24;      // if from nanoAOD
	Bool_t HLT_IsoMu24;      // if from nanoAOD
        UInt_t nGenJet;     // if from nanoAOD



        Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;      // if from nanoAOD
        Bool_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;      // if from nanoAOD
        Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;      // if from nanoAOD
        Bool_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;      // if from nanoAOD

        



// fromNANO = false;  
// 
//    
//     Float_t lheNpNLO;     // if not from nanoAOD
//      Int_t nJets;       // if not from nanoAOD
//      Float_t lheNj_f;   // if not from nanoAOD
//      Int_t nselLeptons; // if not from nanoAOD
//      Int_t HLT_IsoMu27;       // if not from nanoAOD
//      Int_t HLT_IsoTkMu27;       // if not from nanoAOD
//      Int_t HLT_IsoTkMu24;       // if not from nanoAOD
//      Int_t HLT_IsoMu24;       // if not from nanoAOD
//      Int_t nGenJet;    // if not from nanoAOD
      
      

     



// TrackTag nnResult;

TString file_name = std::string(argv[1]);
TString file_tag = std::string(argv[2]);
TString region = std::string(argv[3]); 
int data = atoi(argv[4]);
int applyQCDScaleWeight = atoi(argv[5]);
TString QCDScaleWeight_str = std::string(argv[6]);
int applyJESWeight = atoi(argv[7]);
TString JESWeight_str = std::string(argv[8]);
TString JERWeight_str = std::string(argv[9]);
TString PUWeight_str = std::string(argv[10]);
TString heppyVersion = std::string(argv[11]);
TString postfix = std::string(argv[12]);
TString output = std::string(argv[13]);

// TString heppyVersion = std::string(argv[9]);
// TString postfix = std::string(argv[10]);
// TString output = std::string(argv[11]);

bool MVAtree_to_fill = false;
if ( ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) && ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) && ((JERWeight_str.CompareTo("none")==0)||(JERWeight_str.CompareTo("nom")==0)) && ((PUWeight_str.CompareTo("none")==0)||(PUWeight_str.CompareTo("nom")==0))) MVAtree_to_fill = true;
MVAtree_to_fill = true;


if (data == 1 && !(((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) && ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) && ((JERWeight_str.CompareTo("none")==0)||(JERWeight_str.CompareTo("nom")==0)) && ((PUWeight_str.CompareTo("none")==0)||(PUWeight_str.CompareTo("nom")==0)))) return 0;

ofstream qgl_normalization_file;
qgl_normalization_file.open("qgl_normalization_file.txt", std::ofstream::out | std::ofstream::app); 

std::map <TString, float> xsec;
std::map <TString, float> qgl_norm;
xsec["SingleMuon"] = 1.;
xsec["SingleMuonB"] = 1.;
xsec["SingleMuonB1"] = 1.;
xsec["SingleMuonC"] = 1.;
xsec["SingleMuonD"] = 1.;
xsec["SingleMuonE"] = 1.;
xsec["SingleMuonF"] = 1.;
xsec["SingleMuonG"] = 1.;
xsec["SingleMuonH"] = 1.;
xsec["SingleMuonH1"] = 1.;
xsec["SingleElectron"] =  1.;
xsec["SingleElectronB"] =  1.;
xsec["SingleElectronC"] =  1.;
xsec["SingleElectronD"] =  1.;
xsec["SingleElectronE"] =  1.;
xsec["SingleElectronF"] =  1.;
xsec["SingleElectronG"] =  1.;
xsec["DYJetstoLL"] =  5765.4;
xsec["DYJetstoLL_amc"] =  5765.4;
xsec["DYJetstoLL_HT100"] =  5765.4;
//xsec["DYJetstoLL_HT100_200"] = 173.96106;
//xsec["DYJetstoLL_HT200_400"] = 48.27802 ;
//xsec["DYJetstoLL_HT400_600"] =6.68755 ;
//xsec["DYJetstoLL_HT600_Inf"] = 2.588804;
//xsec["DYJetstoLL_HT100_200"] = 147.40 ; 
//xsec["DYJetstoLL_HT200_400"] = 40.99 ; 
//xsec["DYJetstoLL_HT400_600"] = 5.678 ; 
//xsec["DYJetstoLL_HT600_Inf"] = 2.198; 
xsec["DYJetstoLL_HT100_200"] = 181.302; 
xsec["DYJetstoLL_HT200_400"] =50.4177  ; 
xsec["DYJetstoLL_HT400_600"] =6.98394; 
xsec["DYJetstoLL_HT600_Inf"] =2.70354 ;
xsec["DYJetstoLL_HT600_800"] = 1.6814;
xsec["DYJetstoLL_HT800_1200"] = 0.7754;
xsec["DYJetstoLL_HT1200_2500"] = 0.186;
xsec["DYJetstoLL_HT2500_Inf"] = 0.00438495;





xsec["DYJetsToLL_Pt-0To50_amc"] = 5451; 
xsec["DYJetsToLL_Pt-50To100_amc"] = 354.8; 
xsec["DYJetsToLL_Pt-100To250_amc"] = 81.22; 
xsec["DYJetsToLL_Pt-250To400_amc"] = 2.991 ; 
xsec["DYJetsToLL_Pt-400To650_amc"] = 0.3882 ; 
xsec["DYJetsToLL_Pt-650ToInf_amc"] = 0.03737 ;



xsec["DY1JetsToLL_M-50_LHEZpT_50-150"] = 316.6; 
xsec["DY1JetsToLL_M-50_LHEZpT_150-250"] = 9.543; 
xsec["DY1JetsToLL_M-50_LHEZpT_250-400"] = 1.098; 
xsec["DY1JetsToLL_M-50_LHEZpT_400-inf"] = 0.1193; 


xsec["DY2JetsToLL_M-50_LHEZpT_50-150"] = 169.6; 
xsec["DY2JetsToLL_M-50_LHEZpT_150-250"] = 15.65; 
xsec["DY2JetsToLL_M-50_LHEZpT_250-400"] = 2.737; 
xsec["DY2JetsToLL_M-50_LHEZpT_400-inf"] = 0.4477; 


 

xsec["DYJetsToLL_M-100to200"] = 226.6;
xsec["DYJetsToLL_M-200to400"] = 7.77;
xsec["DYJetsToLL_M-400to500"] = 0.4065;
xsec["DYJetsToLL_M-500to700"] = 0.2334	;
xsec["DYJetsToLL_M-700to800"] = 0.03614;
xsec["DYJetsToLL_M-1500to2000"] = 0.00218;




xsec["DYJetsToLL_Pt-50To100"] = 354.8;
xsec["DYJetsToLL_Pt-100To250"] = 81.22;
xsec["DYJetsToLL_Pt-250To400"] = 2.991;
xsec["DYJetsToLL_Pt-400To650"] = 0.3882;
xsec["DYJetsToLL_Pt-650ToInf"] = 0.03737;



// xsec["DYJetstoLL_Pt-100_amc"] = 5765.4; 
// xsec["DYJetstoLL_Pt-100To250_amc"] = 83.12; 
// xsec["DYJetstoLL_Pt-250To400_amc"] =3.047 ; 
// xsec["DYJetstoLL_Pt-400To650_amc"] = 0.3921 ; 
// xsec["DYJetstoLL_Pt-650ToInf_amc"] = 0.03636 ;


xsec["DYJetsToLL_M-105To160-madgraphMLM"] = 41.25;

xsec["DYJetstoLL_amc_M-50"] = 5765.40;
xsec["DYJetstoLL_amc_Inclusive"] = 4585.27;
xsec["DYJetstoLL_amc_Filter105"] = 41.81; 
xsec["DYJetsToLL_M-105To160-amcatnloFXFX"] = 41.81;//42.73;  //42.95;//43.19;   //4.319e+01 +- 1.837e-01 pb
xsec["DYJetstoLL_amc_0J"] = 4578.455;  // 4620.52; //   4756.180;    //4585.27; //4732.;  normalize to 5765.4 pb, 1.032 
xsec["DYJetstoLL_amc_1J"] = 851.764;  //885.216;   //853.198;//880.5; 
xsec["DYJetstoLL_amc_2J"] = 335.180; //338.26; //341.492;  //325.194;//335.6; 


xsec["DYJetsToTauTau_ForcedMuDecay"] = 5765.40 / 3. * 0.1739 * 0.1739 ;  

xsec["DYJetsToLL_M-50_VBFFilter-amcatnloFXFX"] = 41.81;
xsec["DYJetsToLL_M-50_VBFFilter-madgraphMLM"] = 41.25;
xsec["DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX"] = 41.81*0.0385983;//Mqq>400 pt2>25
xsec["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"] = 41.25*0.0307832/**0.0322284*/;


xsec["DYJetstoLL_madgraph"] = 1.162*4963.; 
xsec["DYJetsToLL_M"] = 1.162*4963.; 
xsec["DY1JetsToLL_M"] = 1.162*1012.; 
xsec["DY2JetsToLL_M"] = 1.162*334.7;
xsec["DY3JetsToLL_M"] = 1.162*95.02;
xsec["DY4JetsToLL_M"] = 1.162*54.52;
xsec["DY0JetsToLL_M"] = xsec["DYJetsToLL_M"] - xsec["DY1JetsToLL_M"] - xsec["DY2JetsToLL_M"] - xsec["DY3JetsToLL_M"] - xsec["DY4JetsToLL_M"]; 

xsec["VBF_HToMuMu"] = 0.0008204;//------------------------------------------------------------------------------------------------
xsec["GluGlu_HToMuMu"] = 0.01053;//------------------------------------------------------------------------------------------------

xsec["TT"] =809.;
xsec["WW"] =118.7;
xsec["WZ"] = 47.13;
xsec["ZZ"] =16.523;
xsec["EWK_LLJJ"]=1.664;
xsec["EWK_LLJJ_herwig"]=1.664;
xsec["EKW_LLJJ_pythia8"]=1.664;
xsec["interference"]=0.128;
xsec["EWK_LLJJ_INT"]=0.128;

xsec["QCD_HT100to200"] = 27990000;
xsec["QCD_HT200to300"] = 1712000 ;
xsec["QCD_HT300to500"] = 347700;
xsec["QCD_HT500to700"] = 32100 ;
xsec["QCD_HT700to1000"] = 6831;
xsec["QCD_HT1000to1500"] = 1207 ;
xsec["QCD_HT1500to2000"] =  119.9;
xsec["QCD_HT2000toInf"] = 25.24;

xsec["ST_tW"] = 71.7 ;			//inclusive decays
xsec["ST_tW_top"] = 35.85  ;			//inclusive decays
xsec["ST_tW_antitop"] = 35.85  ;			//inclusive decays

xsec["ST_s-channel"] = 3.36; //leptonic decays
xsec["ST_t-channel_top_4f_inclusiveDecays"] = 136.02;
xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 80.95;
xsec["ST_t-channel_top_4f_leptonDecays"] = 44.33;  //leptonDecays  , multiplied with BR 0.325
xsec["ST_t-channel_antitop_4f_leptonDecays"] = 26.38;//leptonDecays ,multiplied with BR 0.325


//////////////////////  TOP X SEC = 0  //////////////////////////////////////////////////////////////////////
// xsec["ST_tW"] = 0. ;			//inclusive decays
// xsec["ST_tW_top"] = 0.  ;			//inclusive decays
// xsec["ST_tW_antitop"] = 0.  ;			//inclusive decays
// xsec["ST_s-channel"] = 0.; //leptonic decays
// xsec["ST_t-channel_top_4f_inclusiveDecays"] = 0.;
// xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 0.;
// xsec["ST_t-channel_top_4f_leptonDecays"] = 0.;  //leptonDecays  , multiplied with BR 0.325
// xsec["ST_t-channel_antitop_4f_leptonDecays"] = 0.;//leptonDecays ,multiplied with BR 0.325
// xsec["TT"] =0.;//809.;


xsec["WJetsToLNu_amc"]  = 61526.7; //not going to use these ones

xsec["WJetsToLNu"] =  61526.7;
xsec["WToLNu_0J"]  =  50131.98;
xsec["WToLNu_1J"]  =  8426.09;
xsec["WToLNu_2J"]  =  3172.96;


xsec["W4JetsToLNu"]  =  50131.98;
xsec["W3JetsToLNu"]  =  8426.09;
xsec["W2JetsToLNu"]  =  3172.96;
xsec["W1JetsToLNu"]  =  3172.96;


xsec["WJetsToLNu_HT-100To200"]  =  1395.0;
xsec["WJetsToLNu_HT-200To400"]  =  407.9;
xsec["WJetsToLNu_HT-400To600"]  =  57.48;
xsec["WJetsToLNu_HT-600To800"]  =  12.87;
xsec["WJetsToLNu_HT-1200To2500"]  =  1.074;
xsec["WJetsToLNu_HT-2500ToInf"]  =  0.008001;



xsec["WpWpJJ_EWK"]  =  0.027;
xsec["WWTo2L2Nu_DoubleScattering"]  =  1.62;
xsec["WLLJJ_WToLNu_EWK"]  =  0.0176;
xsec["WWJJToLNuLNu_EWK_noTop"]  =  0.3452;




//////////////////////////////////////// BR (W -> lnu)       = (10.71 +- 0.16) + (10.63 +- 0.15) = 21.34 +- 0.11
///////////// W BRANCHING RATIO //////// BR (W -> tau nu)    = 11.38 +- 0.21
//////////////////////////////////////// BR (W -> had)       = 67.41 +- 0.27

//////////////////////////////////////// BR (Z -> ll)        = (3.363 +- 0.004) + (3.366 +- 0.007) = 6.729 +- 0.16
///////////// Z BRANCHING RATIO //////// BR (Z -> tau tau)   = 3.370 +- 0.008
//////////////////////////////////////// BR (Z -> inv)       = 20.00 +- 0.06
//////////////////////////////////////// BR (Z -> had)       = 69.91 +- 0.06

xsec["WWTo1L1Nu2Q"] = 118.7 * 21.34 * 67.41 / 10000.;
xsec["WZTo1L1Nu2Q"] = 47.13 * 21.34 * 69.91 / 10000.;
xsec["ZZTo2L2Q"]    = 16.523 * 2. * 6.729 * 69.91 / 10000.;



//////////////////////////////////////// BR (t -> W b)        = 1.
///////////// W BRANCHING RATIO //////// 1 - BR (W -> had)    = 32.59 +- 0.27
//////////////////////////////////////// BR (W -> lnu)        = (10.71 +- 0.16) + (10.63 +- 0.15) + (11.38 +- 0.21) = 32.72  +- 0.23
//////////////////////////////////////// 1 - BR (W -> lnu)    = 67.28 +- 0.23

xsec["TTJets_DiLept"]   = 85.65;
xsec["TTTo2L2Nu"]       = 809. * 32.72 * 32.72 / 10000.;
xsec["TTToSemilepton"]  = 809. * 2. * 32.72 * 67.28  / 10000.;
xsec["TTToHadronic"]    = 809. - xsec["TTTo2L2Nu"] - xsec["TTToSemilepton"];


//k factors 1.21 are not included// 
// xsec["WJetsToLNu_HT100To200"] = 1345 ;
// xsec["WJetsToLNu_HT200To400"] = 359.7  ;
// xsec["WJetsToLNu_HT400To600"] = 48.91;
// xsec["WJetsToLNu_HT600To800"] =12.05;
// xsec["WJetsToLNu_HT800To1200"] = 5.501;
// xsec["WJetsToLNu_HT1200To2500"] = 1.329;
// xsec["WJetsToLNu_HT2500ToInf"] = 0.03216;

xsec["WJetsToLNu_HT100"]  = 61526.7;
xsec["WJetsToLNu_HT100To200"] = 1627.45 ;
xsec["WJetsToLNu_HT200To400"] = 435.236  ;
xsec["WJetsToLNu_HT400To600"] = 59.18109;
xsec["WJetsToLNu_HT600To800"] =14.5805;
xsec["WJetsToLNu_HT800To1200"] = 6.656210;
xsec["WJetsToLNu_HT1200To2500"] = 1.608089;
xsec["WJetsToLNu_HT2500ToInf"] = 0.0389135;

xsec["TTZToLLNuNu"] = 0.2529;
xsec["tZq_ll"]=0.0758;


if (region.CompareTo("el")==0) {
qgl_norm["EWK_LLJJ"]=0.938595977;
qgl_norm["EWK_LLJJ_herwig"]=1;
qgl_norm["TT"]=1.05998369;
qgl_norm["WW"]=0.981059169;
qgl_norm["WZ"]=0.956107275;
qgl_norm["ZZ"]=0.970928133;
qgl_norm["ST_tW_antitop"]=0.999659674;
qgl_norm["ST_tW_top"]=0.980208626;
qgl_norm["ST_s-channel"]=0.99449528;
qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.998267176;
qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=0.856539545;
qgl_norm["WJetsToLNu"]=1;
qgl_norm["DYJetstoLL_amc_0J"]=0.996136594;
qgl_norm["DYJetstoLL_amc_1J"]=0.956949008;
qgl_norm["DYJetstoLL_amc_2J"]=0.952277759;
qgl_norm["VBF_HToMuMu"]=1.;
qgl_norm["GluGlu_HToMuMu"]=1.;
}

if (region.CompareTo("mu")==0) {
// qgl_norm["EWK_LLJJ_herwig"]=1;
// qgl_norm["EKW_LLJJ_pythia8"]=1;
// qgl_norm["EWK_LLJJ_INT"]=1;
// qgl_norm["WpWpJJ_EWK"]  =  1;
// qgl_norm["WWTo2L2Nu_DoubleScattering"]  =  1;
// qgl_norm["WLLJJ_WToLNu_EWK"]  =  1;
// qgl_norm["WWJJToLNuLNu_EWK_noTop"]  =  1;

// // // // qgl_norm["EWK_LLJJ"]=0.939774091;
// // // // qgl_norm["TT"]=1.069615388;
// // // // qgl_norm["WW"]=0.927930632;
// // // // qgl_norm["WZ"]=0.967820083;
// // // // qgl_norm["ZZ"]=0.964110717;
// // // // qgl_norm["ST_tW_antitop"]=1.012466766;
// // // // qgl_norm["ST_tW_top"]=0.990673372;
// // // // qgl_norm["ST_s-channel"]=0.911006075;
// // // // qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.986731357;
// // // // qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=0.994925429;
// // // // qgl_norm["WJetsToLNu"]=1;
// // // // qgl_norm["DYJetstoLL_amc_M-50"]=1;      
// // // // qgl_norm["DYJetstoLL_amc_0J"]=1.002031915;
// // // // qgl_norm["DYJetstoLL_amc_1J"]=0.966710372;
// // // // qgl_norm["DYJetstoLL_amc_2J"]=0.954117783;
// // // // 
// // // // 
// // // // qgl_norm["DYJetsToLL_Pt-0To50_amc"] = 1.; 
// // // // qgl_norm["DYJetsToLL_Pt-50To100_amc"] = 1.; 
// // // // qgl_norm["DYJetsToLL_Pt-100To250_amc"] = 1.; 
// // // // qgl_norm["DYJetsToLL_Pt-250To400_amc"] = 1. ; 
// // // // qgl_norm["DYJetsToLL_Pt-400To650_amc"] = 1. ; 
// // // // qgl_norm["DYJetsToLL_Pt-650ToInf_amc"] = 1. ;
// // // // 
// // // // 
// // // // qgl_norm["DYJetsToLL_M-105To160-madgraphMLM"]=1.;
// // // // 
// // // // qgl_norm["DYJetstoLL_madgraph"]=1.;
// // // // qgl_norm["DYJetstoLL_amc_Filter105"]=1.;
// // // // qgl_norm["DYJetstoLL_amc_Inclusive"]=1.;
// // // // qgl_norm["DYJetsToLL_M"]=1.;
// // // // qgl_norm["DY1JetsToLL_M"]=1.;
// // // // qgl_norm["DY2JetsToLL_M"]=1.;
// // // // qgl_norm["DY3JetsToLL_M"]=1.;
// // // // qgl_norm["DY4JetsToLL_M"]=1.;
// // // // 
// // // // qgl_norm["VBF_HToMuMu"]=1.;
// // // // qgl_norm["GluGlu_HToMuMu"]=1.;



qgl_norm["DY1JetsToLL_M-50_LHEZpT_50-150"] = 1.; 
qgl_norm["DY1JetsToLL_M-50_LHEZpT_150-250"] = 1.; 
qgl_norm["DY1JetsToLL_M-50_LHEZpT_250-400"] = 1.; 
qgl_norm["DY1JetsToLL_M-50_LHEZpT_400-inf"] = 1.; 


qgl_norm["DY2JetsToLL_M-50_LHEZpT_50-150"] = 1.; 
qgl_norm["DY2JetsToLL_M-50_LHEZpT_150-250"] = 1.; 
qgl_norm["DY2JetsToLL_M-50_LHEZpT_250-400"] = 1.; 
qgl_norm["DY2JetsToLL_M-50_LHEZpT_400-inf"] = 1.; 


qgl_norm["DYJetsToLL_M-100to200"] = 226.6;
qgl_norm["DYJetsToLL_M-200to400"] = 7.77;
qgl_norm["DYJetsToLL_M-400to500"] = 0.4065;
qgl_norm["DYJetsToLL_M-500to700"] = 0.2334	;
qgl_norm["DYJetsToLL_M-700to800"] = 0.03614;
qgl_norm["DYJetsToLL_M-1500to2000"] = 0.00218;




qgl_norm["DYJetsToLL_Pt-50To100"] = 354.8;
qgl_norm["DYJetsToLL_Pt-100To250"] = 81.22;
qgl_norm["DYJetsToLL_Pt-250To400"] = 2.991;
qgl_norm["DYJetsToLL_Pt-400To650"] = 0.3882;
qgl_norm["DYJetsToLL_Pt-650ToInf"] = 0.03737;



qgl_norm["DYJetsToLL_M"]=0.972332;
qgl_norm["DY0JetsToLL_M"]=0.996343;
qgl_norm["DY1JetsToLL_M"]=0.970908;
qgl_norm["DY2JetsToLL_M"]=0.963111;
qgl_norm["DY3JetsToLL_M"]=0.953453;
qgl_norm["DY4JetsToLL_M"]=0.981916;
// qgl_norm["DYJetsToLL_M-105To160-madgraphMLM"]=0.967472;

qgl_norm["DYJetsToLL_Pt-0To50_amc"]=0.971096;
qgl_norm["DYJetsToLL_Pt-100To250_amc"]=0.953761;
qgl_norm["DYJetsToLL_Pt-250To400_amc"]=0.949793;
qgl_norm["DYJetsToLL_Pt-400To650_amc"]=0.939615;
qgl_norm["DYJetsToLL_Pt-50To100_amc"]=0.95937;

qgl_norm["DYJetstoLL_HT100_200"]=0.961927;
qgl_norm["DYJetstoLL_HT200_400"]=0.964326;
qgl_norm["DYJetstoLL_HT400_600"]=0.97434;
qgl_norm["DYJetstoLL_HT600_800"]=0.940367;
qgl_norm["DYJetstoLL_HT800_1200"]=0.94855;
qgl_norm["DYJetstoLL_HT1200_2500"]=0.937126;
qgl_norm["DYJetstoLL_HT2500_Inf"]=0.93129;

// qgl_norm["DYJetstoLL_amc_0J"]=0.99609;
// qgl_norm["DYJetstoLL_amc_1J"]=0.966311;
// qgl_norm["DYJetstoLL_amc_2J"]=0.961447;
// qgl_norm["DYJetstoLL_amc_Filter105"]=0.964836;
// qgl_norm["DYJetsToLL_M-105To160-amcatnloFXFX"]=0.964836;
qgl_norm["DYJetstoLL_amc_M-50"]=0.965673;

// qgl_norm["DYJetsToLL_M-50_VBFFilter-amcatnloFXFX"] = 1.;
// qgl_norm["DYJetsToLL_M-50_VBFFilter-madgraphMLM"] = 1.;
// qgl_norm["DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX"] = 1.;
// qgl_norm["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"] = 1.;

qgl_norm["DYJetsToTauTau_ForcedMuDecay"] = 1.; 


// qgl_norm["GluGlu_HToMuMu"]=0.987766;
// qgl_norm["VBF_HToMuMu"]=0.974783;

// qgl_norm["ST_s-channel"]=1.45921;
// qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=1.01835;
// qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.999555;
// qgl_norm["ST_tW_antitop"]=0.993577;
// qgl_norm["ST_tW_top"]=1.00096;
// qgl_norm["TT"]=0.994732;

qgl_norm["WJetsToLNu"]=1;
// qgl_norm["WW"]=1.00899;
// qgl_norm["WZ"]=0.941571;
// qgl_norm["ZZ"]=0.930012;



// qgl_norm["WToLNu_0J"]  =  1.;
// qgl_norm["WToLNu_1J"]  =  1.;
// qgl_norm["WToLNu_2J"]  =  1.;


qgl_norm["W4JetsToLNu"]  =  1.;
qgl_norm["W3JetsToLNu"]  =  1.;
qgl_norm["W2JetsToLNu"]  =  1.;
qgl_norm["W1JetsToLNu"]  =  1.;


qgl_norm["WWTo1L1Nu2Q"] = 1.;
qgl_norm["WZTo1L1Nu2Q"] = 0.941571;
qgl_norm["ZZTo2L2Q"]    = 0.930012;

// qgl_norm["TTTo2L2Nu"]       = 1.06;
// qgl_norm["TTToSemilepton"]  = 1.06;
qgl_norm["TTToHadronic"]  = 1.06;





qgl_norm["DYJetstoLL_amc_0J"]  = 1.;
qgl_norm["DYJetstoLL_amc_1J"]  = 0.976177;
qgl_norm["DYJetstoLL_amc_2J"]  = 1.17075;
qgl_norm["EWK_LLJJ_INT"]  = 0.896597;
qgl_norm["EKW_LLJJ_pythia8"]  = 0.934949;
qgl_norm["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"]  = 0.929512;


qgl_norm["WToLNu_0J"]  = 1.;
qgl_norm["WToLNu_1J"]  = 1.;
qgl_norm["WToLNu_2J"]  = 0.821382;

qgl_norm["ST_tW_top"]  = 0.982312;
qgl_norm["ST_tW_antitop"]  = 0.939564;
qgl_norm["ST_s-channel"]  = 1.;
qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]  = 0.97891;
qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]  = 1.04428;


qgl_norm["DYJetsToLL_M-105To160-madgraphMLM"]  = 0.94263;
qgl_norm["DYJetsToLL_M-105To160-amcatnloFXFX"]  = 0.942454;
qgl_norm["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"]  = 0.929779; // This is the one the code compute: therefore it is wrong



qgl_norm["TT"]  = 0.950831;
qgl_norm["TTToSemilepton"]  = 0.947242;
qgl_norm["TTTo2L2Nu"]  = 0.974071;
qgl_norm["TTJets_DiLept"] = 1000;

qgl_norm["GluGlu_HToMuMu"]  = 0.947775;
qgl_norm["VBF_HToMuMu"]  = 0.961363;


qgl_norm["WpWpJJ_EWK"]  = 1.;
qgl_norm["WWTo2L2Nu_DoubleScattering"]  = 0.928416;
qgl_norm["WLLJJ_WToLNu_EWK"]  = 0.93014;
qgl_norm["WWJJToLNuLNu_EWK_noTop"]  = 0.942922;
qgl_norm["WW"]  = 0.881484;
qgl_norm["WZ"]  = 0.944434;
qgl_norm["ZZ"]  = 0.88568;


}


 int counter=0;


int whichQCDScaleWeight;
if ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) whichQCDScaleWeight=0;
if (QCDScaleWeight_str.CompareTo("up")==0) whichQCDScaleWeight=1;
if (QCDScaleWeight_str.CompareTo("down")==0) whichQCDScaleWeight=2;
int whichJESWeight;
if ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) whichJESWeight=0;
if (JESWeight_str.CompareTo("up")==0) whichJESWeight=1;
if (JESWeight_str.CompareTo("down")==0) whichJESWeight=2;

int whichJERWeight = 0;
if (JERWeight_str.CompareTo("up")==0) whichJERWeight=1;
if (JERWeight_str.CompareTo("down")==0) whichJERWeight=2;

int whichPUWeight = 0;
if (PUWeight_str.CompareTo("up")==0) whichPUWeight=1;
if (PUWeight_str.CompareTo("down")==0) whichPUWeight=2;


// if ((JESWeight_str.CompareTo("up")==0) && (QCDScaleWeight_str.CompareTo("up")==0)) whichJERWeight=1;
// if ((JESWeight_str.CompareTo("down")==0) && (QCDScaleWeight_str.CompareTo("down")==0)) whichJERWeight=2;

if ( !((file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0))) whichQCDScaleWeight = 0;

float gen_pos=0; 
float gen_neg=0; 
float gen_pos_weight=0; 
float gen_neg_weight=0; 

	int MVAcount = 0;
	int MVAcountMAX = 500000;
        
	Float_t presel=0;
	Float_t presel_vtype[10] = {0,0,0,0,0,0,0,0,0};
	Float_t pos_puweight=0;
	Float_t all_puweight=0.;
	Float_t puweight = 1;
      Float_t puweightUp = 1;
	Float_t puweightDown = 1;
	Float_t PU=1.;
	Float_t genweight;
	Float_t bTagWeight;
	double genweight0;
        Float_t genweight_noQGLcorrection;
	float  trigWeight_tree;
// 	Int_t global_counter = 0;
	Int_t HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200;
	Int_t HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460;
	Int_t HLT_IsoMu22;
	Int_t HLT_IsoTkMu22;

        
        float qqMass_skim = 0;
        float MuonMass_skim = 0;
        
	Int_t HLT_Ele27_eta2p1;
	TFile *file_initial;
	TChain *tree_initial;

	Int_t nvLeptons;

	const int brLeptons=13;
	Float_t vLeptons_pt[30], vLeptons_eta[30], vLeptons_phi[30], vLeptons_mass[30], vLeptons_SF_IdCutLoose[30], vLeptons_SF_IdCutTight[30], vLeptons_SF_IsoLoose[30], vLeptons_SF_IsoTight[30],vLeptons_SF_trk_eta[30], vLeptons_SF_HLT_RunD4p2[30],vLeptons_SF_HLT_RunD4p3[30], vLeptons_relIso03[30], vLeptons_eleSieie[30], vLeptons_eleHoE[30], vLeptons_eleDEta[30],vLeptons_eleDPhi[30], vLeptons_eleEcalClusterIso[30], vLeptons_eleHcalClusterIso[30],vLeptons_dr03TkSumPt[30]  ;
	Int_t vLeptons_charge[30], vLeptons_pdgId[30],vLeptons_trackerLayers[30] ; 

	Float_t selLeptons_pt[30], selLeptons_eta[30], selLeptons_phi[30], selLeptons_mass[30], selLeptons_SF_IdCutLoose[30], selLeptons_SF_IdCutTight[30], selLeptons_SF_IsoLoose[30], selLeptons_SF_IsoTight[30],selLeptons_SF_trk_eta[30], selLeptons_SF_HLT_RunD4p2[30],selLeptons_SF_HLT_RunD4p3[30], selLeptons_relIso04[30], selLeptons_relIso03[30], selLeptons_eleSieie[30], selLeptons_eleHoE[30], selLeptons_eleDEta[30],selLeptons_eleDPhi[30], selLeptons_eleEcalClusterIso[30], selLeptons_eleHcalClusterIso[30],selLeptons_dr03TkSumPt[30], selLeptons_dz[30], selLeptons_dxy[30] ;

        Float_t Electron_pt[10], Electron_eta[10], Electron_phi[10], Electron_mass[10], Electron_pfRelIso03_all[10], Electron_dz[10], Electron_dxy[10];
        UInt_t  nElectron=0;
        Int_t Electron_charge[10],Electron_pdgId[10];

        Float_t Photon_eta[10], Photon_pt[10],Photon_pfRelIso03_all[10];
        
        Bool_t Photon_mvaID_WP90[10], Electron_mvaSpring16GP_WP90[10];
        
        
//     UInt_t nGenPart;
// 	Float_t GenPart_pt[30];
// 	Float_t GenPart_eta[30];
// 	Float_t GenPart_phi[30];
// 	Float_t GenPart_mass[30];
// 	Int_t GenPart_pdgId[30];

//         std::vector<Float_t> GenPart_pt;
// 	std::vector<Float_t> GenPart_eta;
// 	std::vector<Float_t> GenPart_phi;
// 	std::vector<Float_t> GenPart_mass;
// 	std::vector<Int_t> GenPart_pdgId;
        
    Int_t nGenLep;
	Float_t GenLep_pt[30];
	Float_t GenLep_eta[30];
	Float_t GenLep_phi[30];
	Float_t GenLep_mass[30];
	Int_t GenLep_pdgId[30];
	
        
        UInt_t nGenPart;
	Float_t GenPart_pt[300];
	Float_t GenPart_eta[300];
	Float_t GenPart_phi[300];
	Float_t GenPart_mass[300];
	Int_t GenPart_pdgId[300];
        Int_t GenPart_statusFlags[300];

        
    Int_t nGenLepRecovered;
	Float_t GenLepRecovered_pt[30];
	Float_t GenLepRecovered_eta[30];
	Float_t GenLepRecovered_phi[30];
	Float_t GenLepRecovered_mass[30];
	Int_t GenLepRecovered_pdgId[30];
	
        

        

        std::vector<float> genJetsWithoutLeptonsP4_pt;
        std::vector<float> genJetsWithoutLeptonsP4_eta;
        std::vector<float> genJetsWithoutLeptonsP4_phi;
        std::vector<float> genJetsWithoutLeptonsP4_mass;
        int ngenJetsWithoutLeptonsP4;
    

        Float_t GenJet_pt[30];
        Float_t GenJet_eta[30];
        Float_t GenJet_phi[30];
        Float_t GenJet_mass[30];

        
	Int_t selLeptons_charge[30], selLeptons_pdgId[30], selLeptons_looseIdPOG[30]/*, selLeptons_mediumIdPOG[30], selLeptons_tightIdPOG[30]*/, selLeptons_trackerLayers[30],  selLeptons_eleMVAIdSppring16GenPurp[30]; 
        UChar_t Muon_genPartFlav[30]; 
        Bool_t selLeptons_softIdPOG[30], selLeptons_mediumIdPOG[30], selLeptons_tightIdPOG[30];
        
	TString str_leptons[brLeptons] = {"vLeptons_pt", "vLeptons_eta", "vLeptons_phi", "vLeptons_mass", "vLeptons_charge", "vLeptons_pdgId", "vLeptons_SF_IdCutLoose", "vLeptons_SF_IdCutTight", "vLeptons_SF_IsoLoose","vLeptons_SF_IsoTight","vLeptons_SF_trk_eta","vLeptons_SF_HLT_RunD4p2","vLeptons_SF_HLT_RunD4p3"};


	
////////////////////////////


        TFile* file_id_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunBCDEF.root");
        TH2F* id_mu_bf = (TH2F*)file_id_mu_bf->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
        TFile* file_id_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunGH.root");
        TH2F* id_mu_aft = (TH2F*)file_id_mu_aft->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");

        TFile* file_trig_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunBCDEF.root");
        TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
        TFile* file_trig_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunGH.root");
        TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");

        TFile* file_iso_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_LooseISO_LooseID_pt_eta_RunBCDEF.root");
        TH2F* iso_mu_bf = (TH2F*)file_iso_mu_bf->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
        TFile* file_iso_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_LooseISO_LooseID_pt_eta_RunGH.root");
        TH2F* iso_mu_aft = (TH2F*)file_iso_mu_aft->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");


        TFile* file_tracker_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_ScaleFactor_tracker_80x.root");
        TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");
        TFile* file_trig_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_Tight27AfterIDISO.root");
        TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_Tight27AfterIDISO");
        TFile* file_id_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_EIDISO_ZH.root");
        TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_EIDISO_ZH");

        TFile* file_track_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_Muons_trk_SF_RunBCDEF.root");
        TH1F* track_mu_bf = (TH1F*)file_track_mu_bf->Get("TriggerEffMap_Graph");
        TFile* file_track_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_Muons_trk_SF_RunGH.root");
        TH1F* track_mu_aft = (TH1F*)file_track_mu_aft->Get("TriggerEffMap_Graph");

        
        
        
        
        
        
        
        
        
//         TFile* f_param= TFile::Open("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/out.root");
        TFile* f_param= TFile::Open("out.root");
        TH2F *hmuon= (TH2F*)f_param->Get("PtErrParametrization");
        
        
        TFile* f_forwardJetWeight= TFile::Open("forwardJetPtCorrection/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root");
//         TFile* f_forwardJetWeight= TFile::Open("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/forwardJetPtCorrection/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root");
        TEfficiency *h_forwardJetWeight= (TEfficiency*)f_forwardJetWeight->Get("prefireEfficiencyMap");
//         TFile* f_IsoEG30= TFile::Open("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/forwardJetPtCorrection/Map_Jet_L1IsoEG30eff_bxm1_looseJet_SingleMuon_Run2016B-H.root");
        TFile* f_IsoEG30= TFile::Open("forwardJetPtCorrection/Map_Jet_L1IsoEG30eff_bxm1_looseJet_SingleMuon_Run2016B-H.root");
        TEfficiency *h_IsoEG30Weight= (TEfficiency*)f_IsoEG30->Get("prefireEfficiencyMap");
        
//         cout<< "fparam "<<f_param<<endl;
//         cout<< "hmuon  "<<hmuon<<endl;
// 	RoccoR  *rc = new RoccoR("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/rcdata.2016.v3/");
        RoccoR  *rc = new RoccoR("2016/rcdata.2016.v3/");

	
	file_initial = TFile::Open(file_name);
	
	tree_initial = (TChain*)file_initial->Get("tree");
	if (fromNANO) tree_initial = (TChain*)file_initial->Get("Events");
// 	Int_t events_generated;
	double events_generated;

        Float_t events_generated_muF_QCDUp=1;
	Float_t events_generated_muF_QCDDown=1;
	Float_t events_generated_muR_QCDUp=1;
	Float_t events_generated_muR_QCDDown=1;
	Float_t events_generated_muFR_QCDUp=1;
	Float_t events_generated_muFR_QCDDown=1;
	TH1F *countPos;
	TH1F *countNeg;
	TH1F *countLHEScale;
	TH1F *countLHEPdf;
	TH1F *countWeighted;
	if ((data!=1)){
            if (fromNANO) {
                TChain* tree_Runs = (TChain*)file_initial->Get("Runs");

                
                events_generated = SumEntries(tree_Runs, "genEventSumw", false, -1);
//                 events_generated = SumEntries(tree_Runs, "genEventCount", false, -1);
                
//                 if (whichQCDScaleWeight==1) events_generated = events_generated * SumEntries(tree_Runs, "LHEScaleSumw", true, 8);    // whichQCDScaleWeight==1    is up variation
//                 if (whichQCDScaleWeight==2) events_generated = events_generated * SumEntries(tree_Runs, "LHEScaleSumw", true, 0);    // whichQCDScaleWeight==1    is down variation
                
                if (whichQCDScaleWeight==1) events_generated = SumEntries(tree_Runs, "LHEScaleSumw", true, 8);    // whichQCDScaleWeight==1    is up variation
                if (whichQCDScaleWeight==2) events_generated = SumEntries(tree_Runs, "LHEScaleSumw", true, 0);    // whichQCDScaleWeight==1    is down variation
                

                std::cout << "events_generated    " << whichQCDScaleWeight << "     " << events_generated << std::endl;
                

                
            }
             else {   
                countPos = (TH1F*)file_initial->Get("CountPosWeight");
                countNeg = (TH1F*)file_initial->Get("CountNegWeight");
                countWeighted = (TH1F*)file_initial->Get("CountWeighted");
                countLHEScale = (TH1F*)file_initial->Get("CountWeightedLHEWeightScale");
                countLHEPdf=	(TH1F*)file_initial->Get("CountWeightedLHEWeightPdf");
                //	events_generated = countWeighted->GetBinContent(1);
                //	if (whichQCDScaleWeight==0) events_generated = countPos->GetBinContent(1) - countNeg->GetBinContent(1);
                if (whichQCDScaleWeight==0) 	events_generated = countWeighted->GetBinContent(1);
                else {
                    events_generated = countLHEScale->GetBinContent( countLHEScale->FindBin( 3 + whichQCDScaleWeight) );
                    if (events_generated==0) events_generated =  countPos->GetBinContent(1) - countNeg->GetBinContent(1);
                }
                events_generated_muF_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(0) );
                events_generated_muF_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(1) );
                events_generated_muR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(2) );
                events_generated_muR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(3) );
                events_generated_muFR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(4) );
                events_generated_muFR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(5) );
            }
        }
    else events_generated = 1;

//	if (file_tag.CompareTo("EWK_LLJJ")==0)  events_generated=events_generated/2.00703;
    Jets Jet;
    Float_t v_type;
    Float_t v_typeSim;
    Float_t wrong_type=0.;


	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	
	Float_t lheHT;	

	
        int lheNj;
	Float_t lheV_pt;	

        
        

	Float_t bdt;
	int passSel;
	int passSel_JESR[4];
	Float_t bdt_JESR[4];

	Float_t met_pt;
	Float_t met_phi;
        
        Float_t met_pt_JERup;
        Float_t met_phi_JERup;      
        Float_t met_pt_JESup;
        Float_t met_phi_JESup; 
        Float_t met_pt_JERdown;
        Float_t met_phi_JERdown;      
        Float_t met_pt_JESdown;
        Float_t met_phi_JESdown; 
        
   

        
        
	Jets GenHiggsSisters;

	int pos_weight_presel=0;
 	Int_t selLeptons_tightId[20];
	Float_t  selLeptons_chargedHadRelIso03[20], selLeptons_pfRelIso03[20];
	Float_t vLeptons_dz[20], vLeptons_edz[20];
	Float_t cosThetaStar = 0;
	Float_t cosThetaStarJet = 0;
        Float_t absCosThetaStarJet = 0;
	Float_t cosThetaPlane = 0;


	Int_t nGenVbosons;
	Float_t GenVbosons_pt[1];
	Int_t GenVbosons_pdgId[1];
	Float_t VtypeSim; 
        
        const int number_LHE_weights_pdf = 103;
Float_t LHE_weights_pdf_wgt[number_LHE_weights_pdf];
Float_t LHE_weights_scale_wgt[10];
	
	float V_mass;
	ULong64_t evt;
	UInt_t run;
        UInt_t luminosityBlock;
        
        int runToWrite = 0;
        int lumiToWrite = 0;
	ULong64_t event; 
//         ratioDYFilter ratiosClass;
//         ratiosClass.setNumber();

        
    float BDT_VBF = 0;
    
    
    
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,0.01);
    
    
    
//     bool fromNANO = true;
    if (fromNANO) {//if skim from nanoAOD

//    tree_initial->SetBranchAddress("nGenPart",&nGenPart);
//    tree_initial->SetBranchAddress("GenPart_pt",&GenPart_pt);
//    tree_initial->SetBranchAddress("GenPart_eta",&GenPart_eta);
//    tree_initial->SetBranchAddress("GenPart_phi",&GenPart_phi);
//    tree_initial->SetBranchAddress("GenPart_mass",&GenPart_mass);
//    tree_initial->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
        

   
//    tree_initial->SetBranchAddress("Vtype",&v_type);
//    tree_initial->SetBranchAddress("V_mass",&V_mass);
//    tree_initial->SetBranchAddress("rho",&rho);
    tree_initial->SetBranchAddress("nJet",&nJets);
//     tree_initial->SetBranchAddress("Jet_pt_nom",Jet.pt);
//    tree_initial->SetBranchAddress("Jet_corr_JECUp",Jet.JEC_corr_up);
//    tree_initial->SetBranchAddress("Jet_corr_JECDown",Jet.JEC_corr_down);
//    tree_initial->SetBranchAddress("Jet_corr",Jet.JEC_corr);
//    tree_initial->SetBranchAddress("Jet_corr_JERUp",Jet.JER_corr_up);
//    tree_initial->SetBranchAddress("Jet_corr_JERDown",Jet.JER_corr_down);
//    tree_initial->SetBranchAddress("Jet_corr_JER",Jet.JER_corr);

          
    
        if (data==1) {
            tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
            tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
        }
        else {
            tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
            tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
//             tree_initial->SetBranchAddress("Jet_pt_nom",Jet.pt);
//             tree_initial->SetBranchAddress("Jet_mass_nom",Jet.mass);
            tree_initial->SetBranchAddress("Jet_pt_jerUp",Jet.pt_JERup);
            tree_initial->SetBranchAddress("Jet_mass_jerUp",Jet.mass_JERup);
            tree_initial->SetBranchAddress("Jet_pt_jesTotalUp",Jet.pt_JESup);
            tree_initial->SetBranchAddress("Jet_mass_jesTotalUp",Jet.mass_JESup);
            tree_initial->SetBranchAddress("Jet_pt_jerDown",Jet.pt_JERdown);
            tree_initial->SetBranchAddress("Jet_mass_jerDown",Jet.mass_JERdown);
            tree_initial->SetBranchAddress("Jet_pt_jesTotalDown",Jet.pt_JESdown);
            tree_initial->SetBranchAddress("Jet_mass_jesTotalDown",Jet.mass_JESdown);
        }
            
     
tree_initial->SetBranchAddress("jetIdx1",&Jet.jetIdx1);
tree_initial->SetBranchAddress("jetIdx2",&Jet.jetIdx2);

        tree_initial->SetBranchAddress("Jet_chEmEF",Jet.chEmEF);
        tree_initial->SetBranchAddress("Jet_neEmEF",Jet.neEmEF);
        tree_initial->SetBranchAddress("Jet_chHEF",Jet.chHEF);
        tree_initial->SetBranchAddress("Jet_neHEF",Jet.neHEF);



        tree_initial->SetBranchAddress("Jet_eta",Jet.eta);
        tree_initial->SetBranchAddress("Jet_phi",Jet.phi);
	tree_initial->SetBranchAddress("Jet_btagCSVV2",Jet.btagCSV);
	tree_initial->SetBranchAddress("Jet_btagCMVA",Jet.btagCMVA);
//	tree_initial->SetBranchAddress("Jet_blike_VBF",Jet.blike_VBF);
	tree_initial->SetBranchAddress("Jet_jetId",Jet.id);	
	tree_initial->SetBranchAddress("Jet_puId",Jet.puId);
// 	tree_initial->SetBranchAddress("Jet_leadTrackPt",Jet.leadTrackPt);
	tree_initial->SetBranchAddress("Jet_partonFlavour",Jet.partonFlavour);

        tree_initial->SetBranchAddress("Jet_qgl",Jet.qgl);

	tree_initial->SetBranchAddress("Jet_ptd",Jet.ptd);
// 	tree_initial->SetBranchAddress("Jet_axis2",Jet.axis2);
	tree_initial->SetBranchAddress("Jet_area",Jet.axis2);
// 	tree_initial->SetBranchAddress("Jet_mult",Jet.mult);
        tree_initial->SetBranchAddress("Jet_nConstituents",Jet.mult);
        tree_initial->SetBranchAddress("Jet_VBFselected",Jet.VBFselected);
        tree_initial->SetBranchAddress("Jet_muonIdx1",Jet.muonIdx1);
        tree_initial->SetBranchAddress("Jet_muonIdx2",Jet.muonIdx2);
        tree_initial->SetBranchAddress("Jet_electronIdx1",Jet.electronIdx1);
        tree_initial->SetBranchAddress("Jet_electronIdx2",Jet.electronIdx2);
        
        tree_initial->SetBranchAddress("Jet_genJetIdx",Jet.Jet_genJetIdx);
        
        if (data==1) {
            tree_initial->SetBranchAddress("MET_pt",&met_pt);
            tree_initial->SetBranchAddress("MET_phi",&met_phi);
        }
        else {
//             tree_initial->SetBranchAddress("MET_pt",&met_pt);
//             tree_initial->SetBranchAddress("MET_phi",&met_phi);
            tree_initial->SetBranchAddress("MET_pt_nom",&met_pt);
            tree_initial->SetBranchAddress("MET_phi_nom",&met_phi);
            tree_initial->SetBranchAddress("MET_pt_jerUp",&met_pt_JERup);
            tree_initial->SetBranchAddress("MET_phi_jerUp",&met_phi_JERup);
            tree_initial->SetBranchAddress("MET_pt_jesTotalUp",&met_pt_JESup);
            tree_initial->SetBranchAddress("MET_phi_jesTotalUp",&met_phi_JESup);
            tree_initial->SetBranchAddress("MET_pt_jerDown",&met_pt_JERdown);
            tree_initial->SetBranchAddress("MET_phi_jerDown",&met_phi_JERdown);
            tree_initial->SetBranchAddress("MET_pt_jesTotalDown",&met_pt_JESdown);
            tree_initial->SetBranchAddress("MET_phi_jesTotalDown",&met_phi_JESdown); 
        }


        tree_initial->SetBranchAddress("MuonMass",&MuonMass_skim);
        tree_initial->SetBranchAddress("qqMass",&qqMass_skim);
        
        
        
//         tree_initial->SetBranchAddress("SoftActivityJet_pt",Jet.soft_pt);
//         tree_initial->SetBranchAddress("SoftActivityJet_eta",Jet.soft_eta);
//         tree_initial->SetBranchAddress("SoftActivityJet_mass",Jet.soft_mass);
        tree_initial->SetBranchAddress("SoftActivityJetHT",&Jet.HTsoft);
//         tree_initial->SetBranchAddress("SoftActivityJetNjets2",&Jet.nsoft2);
//         tree_initial->SetBranchAddress("SoftActivityJetNjets5",&Jet.nsoft5);
//         tree_initial->SetBranchAddress("SoftActivityJetNjets10",&Jet.nsoft10);
//	tree_initial->SetBranchAddress("softActivityEWK_HT",&Jet.EWKHTsoft);
//	tree_initial->SetBranchAddress("softActivityEWK_njets2",&Jet.EWKnsoft2);
//	tree_initial->SetBranchAddress("softActivityEWK_njets5",&Jet.EWKnsoft5);
//	tree_initial->SetBranchAddress("softActivityEWK_njets10",&Jet.EWKnsoft10);
        tree_initial->SetBranchAddress("SoftActivityJetNjets2",&Jet.EWKnsoft2);
        tree_initial->SetBranchAddress("SoftActivityJetNjets5",&Jet.EWKnsoft5);
        tree_initial->SetBranchAddress("SoftActivityJetNjets10",&Jet.EWKnsoft10);
        
	tree_initial->SetBranchAddress("nSoftActivityJet",&Jet.nsoftActivityEWKJets);
	tree_initial->SetBranchAddress("SoftActivityJet_pt",Jet.EWKsoft_pt);
	tree_initial->SetBranchAddress("SoftActivityJet_eta",Jet.EWKsoft_eta);
	tree_initial->SetBranchAddress("SoftActivityJet_phi",Jet.EWKsoft_phi);
//	tree_initial->SetBranchAddress("softActivityEWKJets_mass",Jet.EWKsoft_mass);

	tree_initial->SetBranchAddress("genWeight",&genweight);
//	tree_initial->SetBranchAddress("bTagWeight",&bTagWeight);
	tree_initial->SetBranchAddress("puWeight",&puweight);
	tree_initial->SetBranchAddress("puWeightUp",&puweightUp);
	tree_initial->SetBranchAddress("puWeightDown",&puweightDown);
	//tree_initial->SetBranchAddress("PhiZQ1",&PhiZQ1);
	tree_initial->SetBranchAddress("nPVs",&nPVs);

// 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",&HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200);
// 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",&HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460);
// 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu22_v",&HLT_IsoMu22);
// 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu22_v",&HLT_IsoTkMu22);
 	tree_initial->SetBranchAddress("HLT_IsoMu24",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_IsoTkMu24",&HLT_IsoTkMu24);
 	tree_initial->SetBranchAddress("HLT_IsoMu27",&HLT_IsoMu27);
//  	tree_initial->SetBranchAddress("HLT_IsoTkMu27",&HLT_IsoTkMu27);
        
        tree_initial->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
        tree_initial->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);   
        tree_initial->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
        tree_initial->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);   
        
        
      

// 	tree_initial->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v",&HLT_Ele27_eta2p1);
//	tree_initial->SetBranchAddress("Jet_pt_regVBF",Jet.pt_regVBF);
	tree_initial->SetBranchAddress("json",&JSON);
    
//	tree_initial->SetBranchAddress("GenHiggsSisters_pt",GenHiggsSisters.pt);
//      tree_initial->SetBranchAddress("GenHiggsSisters_eta",GenHiggsSisters.eta);
//      tree_initial->SetBranchAddress("GenHiggsSisters_phi",GenHiggsSisters.phi);
//	tree_initial->SetBranchAddress("GenHiggsSisters_mass",GenHiggsSisters.mass);
//	tree_initial->SetBranchAddress("nGenVbosons",&nGenVbosons);
//	tree_initial->SetBranchAddress("GenVbosons_pt",GenVbosons_pt);
//	tree_initial->SetBranchAddress("GenVbosons_pdgId",GenVbosons_pdgId);
//	tree_initial->SetBranchAddress("VtypeSim",&VtypeSim);

        tree_initial->SetBranchAddress("nGenJet",&nGenJet);
        tree_initial->SetBranchAddress("GenJet_pt",GenJet_pt);
        tree_initial->SetBranchAddress("GenJet_eta",GenJet_eta);
        tree_initial->SetBranchAddress("GenJet_phi",GenJet_phi);
        tree_initial->SetBranchAddress("GenJet_mass",GenJet_mass);
                    
        tree_initial->SetBranchAddress("nGenPart",&nGenPart);
        tree_initial->SetBranchAddress("GenPart_pt",GenPart_pt);
        tree_initial->SetBranchAddress("GenPart_eta",GenPart_eta);
        tree_initial->SetBranchAddress("GenPart_phi",GenPart_phi);
        tree_initial->SetBranchAddress("GenPart_mass",GenPart_mass);
        tree_initial->SetBranchAddress("GenPart_pdgId",GenPart_pdgId);
        tree_initial->SetBranchAddress("GenPart_statusFlags",GenPart_statusFlags);

        
        
        
// 	tree_initial->SetBranchAddress("nvLeptons",&nvLeptons);
// 	tree_initial->SetBranchAddress("vLeptons_pt",vLeptons_pt);
// 	tree_initial->SetBranchAddress("vLeptons_eta",vLeptons_eta);
// 	tree_initial->SetBranchAddress("vLeptons_phi",vLeptons_phi);
// 	tree_initial->SetBranchAddress("vLeptons_mass",vLeptons_mass);
// 	tree_initial->SetBranchAddress("vLeptons_charge",vLeptons_charge);
// 	tree_initial->SetBranchAddress("vLeptons_pdgId",vLeptons_pdgId);
// 	tree_initial->SetBranchAddress("vLeptons_SF_IdCutLoose",vLeptons_SF_IdCutLoose);
// 	tree_initial->SetBranchAddress("vLeptons_SF_IdCutTight",vLeptons_SF_IdCutTight);
// 	tree_initial->SetBranchAddress("vLeptons_SF_IsoLoose",vLeptons_SF_IsoLoose);
// 	tree_initial->SetBranchAddress("vLeptons_SF_IsoTight",vLeptons_SF_IsoTight);
// 	tree_initial->SetBranchAddress("vLeptons_SF_trk_eta",vLeptons_SF_trk_eta);
// 	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p2",vLeptons_SF_HLT_RunD4p2);
// 	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p3",vLeptons_SF_HLT_RunD4p3);
// 	tree_initial->SetBranchAddress("vLeptons_trackerLayers", vLeptons_trackerLayers);





        tree_initial->SetBranchAddress("nMuon",&nselLeptons);
        if (data==1) tree_initial->SetBranchAddress("Muon_pt",selLeptons_pt);
        else         tree_initial->SetBranchAddress("Muon_pt_corrected",selLeptons_pt);
//         else         tree_initial->SetBranchAddress("Muon_pt",selLeptons_pt);
        tree_initial->SetBranchAddress("Muon_eta",selLeptons_eta);
        tree_initial->SetBranchAddress("Muon_phi",selLeptons_phi);
        tree_initial->SetBranchAddress("Muon_mass",selLeptons_mass);
        tree_initial->SetBranchAddress("Muon_charge",selLeptons_charge);
        tree_initial->SetBranchAddress("Muon_pdgId",selLeptons_pdgId);
        tree_initial->SetBranchAddress("Muon_genPartFlav",Muon_genPartFlav);

//         tree_initial->SetBranchAddress("selLeptons_looseIdPOG",selLeptons_looseIdPOG);   # loose cut on muons is already implemented
        tree_initial->SetBranchAddress("Muon_softId",selLeptons_softIdPOG);   
        tree_initial->SetBranchAddress("Muon_mediumId",selLeptons_mediumIdPOG);   
        tree_initial->SetBranchAddress("Muon_tightId",selLeptons_tightIdPOG);   
        tree_initial->SetBranchAddress("Muon_dz",selLeptons_dz);          // abs(Muon_dz) < 0.2 is not implemented in nanoAOD
        tree_initial->SetBranchAddress("Muon_dxy",selLeptons_dxy);        // Muon_dxy < 0.05 is not implemented in nanoAOD
        tree_initial->SetBranchAddress("Muon_pfRelIso04_all",selLeptons_relIso04);
        tree_initial->SetBranchAddress("Muon_pfRelIso03_all",selLeptons_relIso03);
        tree_initial->SetBranchAddress("Muon_nTrackerLayers", selLeptons_trackerLayers);
//      tree_initial->SetBranchAddress("selLeptons_eleMVAIdSppring16GenPurp",selLeptons_eleMVAIdSppring16GenPurp);
// 	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
// 	tree_initial->SetBranchAddress("selLeptons_eleHoE",selLeptons_eleHoE);
// 	tree_initial->SetBranchAddress("selLeptons_eleDEta",selLeptons_eleDEta);
// 	tree_initial->SetBranchAddress("selLeptons_eleDPhi",selLeptons_eleDPhi);
// 	tree_initial->SetBranchAddress("selLeptons_eleEcalClusterIso",selLeptons_eleEcalClusterIso);
// 	tree_initial->SetBranchAddress("selLeptons_eleHcalClusterIso",selLeptons_eleHcalClusterIso);
// 	tree_initial->SetBranchAddress("selLeptons_dr03TkSumPt",selLeptons_dr03TkSumPt);
	
        
        
        tree_initial->SetBranchAddress("nElectron",&nElectron);
        tree_initial->SetBranchAddress("Electron_pt",Electron_pt);
        tree_initial->SetBranchAddress("Electron_eta",Electron_eta);
        tree_initial->SetBranchAddress("Electron_phi",Electron_phi);
        tree_initial->SetBranchAddress("Electron_mass",Electron_mass);
        tree_initial->SetBranchAddress("Electron_charge",Electron_charge);
        tree_initial->SetBranchAddress("Electron_pdgId",Electron_pdgId);   
        tree_initial->SetBranchAddress("Electron_pfRelIso03_all",Electron_pfRelIso03_all);   
        tree_initial->SetBranchAddress("Electron_dz",Electron_dz);   
        tree_initial->SetBranchAddress("Electron_dxy",Electron_dxy);   
        tree_initial->SetBranchAddress("Electron_mvaSpring16GP_WP90",Electron_mvaSpring16GP_WP90);
        
        
        tree_initial->SetBranchAddress("Jet_photonIdx1",Jet.photonIdx1);   
        tree_initial->SetBranchAddress("Jet_photonIdx2",Jet.photonIdx2);   

        tree_initial->SetBranchAddress("Photon_eta",Photon_eta);   
        tree_initial->SetBranchAddress("Photon_pt",Photon_pt);   
        tree_initial->SetBranchAddress("Photon_mvaID_WP90",Photon_mvaID_WP90);   
        tree_initial->SetBranchAddress("Photon_pfRelIso03_all",Photon_pfRelIso03_all);   
        
           


                
                
	tree_initial->SetBranchAddress("event",&evt);
	tree_initial->SetBranchAddress("run",&run);
	tree_initial->SetBranchAddress("luminosityBlock",&luminosityBlock);

        tree_initial->SetBranchAddress("LHEPdfWeight",LHE_weights_pdf_wgt);
        tree_initial->SetBranchAddress("LHEScaleWeight",LHE_weights_scale_wgt);
        tree_initial->SetBranchAddress("LHE_HT",&lheHT);
        tree_initial->SetBranchAddress("LHE_Njets",&lheNj_f);
        tree_initial->SetBranchAddress("LHE_Vpt",&lheV_pt);
        tree_initial->SetBranchAddress("LHE_NpNLO",&lheNpNLO);       
        
    }
    else {

    tree_initial->SetBranchAddress("nGenLep",&nGenLep);
    tree_initial->SetBranchAddress("GenLep_pt",GenLep_pt);
    tree_initial->SetBranchAddress("GenLep_eta",GenLep_eta);
    tree_initial->SetBranchAddress("GenLep_phi",GenLep_phi);
    tree_initial->SetBranchAddress("GenLep_mass",GenLep_mass);
    tree_initial->SetBranchAddress("GenLep_pdgId",GenLep_pdgId);
    tree_initial->SetBranchAddress("nGenLepRecovered",&nGenLepRecovered);
    tree_initial->SetBranchAddress("GenLepRecovered_pt",GenLepRecovered_pt);
    tree_initial->SetBranchAddress("GenLepRecovered_eta",GenLepRecovered_eta);
    tree_initial->SetBranchAddress("GenLepRecovered_phi",GenLepRecovered_phi);
    tree_initial->SetBranchAddress("GenLepRecovered_mass",GenLepRecovered_mass);
    tree_initial->SetBranchAddress("GenLepRecovered_pdgId",GenLepRecovered_pdgId);
    tree_initial->SetBranchAddress("Vtype",&v_type);
//     tree_initial->SetBranchAddress("VtypeSim",&v_typeSim);
    tree_initial->SetBranchAddress("V_mass",&V_mass);
    tree_initial->SetBranchAddress("rho",&rho);
    tree_initial->SetBranchAddress("nJet",&nJets);
    tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
    tree_initial->SetBranchAddress("Jet_corr_JECUp",Jet.JEC_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JECDown",Jet.JEC_corr_down);
    tree_initial->SetBranchAddress("Jet_corr",Jet.JEC_corr);
    tree_initial->SetBranchAddress("Jet_corr_JERUp",Jet.JER_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JERDown",Jet.JER_corr_down);
    tree_initial->SetBranchAddress("Jet_corr_JER",Jet.JER_corr);
    tree_initial->SetBranchAddress("Jet_eta",Jet.eta);
    tree_initial->SetBranchAddress("Jet_phi",Jet.phi);
	tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
	tree_initial->SetBranchAddress("Jet_btagCSV",Jet.btagCSV);
	tree_initial->SetBranchAddress("Jet_btagCMVA",Jet.btagCMVA);
	tree_initial->SetBranchAddress("Jet_blike_VBF",Jet.blike_VBF);
	tree_initial->SetBranchAddress("Jet_id",Jet.id);	
	tree_initial->SetBranchAddress("Jet_puId",Jet.puId);
 	tree_initial->SetBranchAddress("Jet_leadTrackPt",Jet.leadTrackPt);
 	tree_initial->SetBranchAddress("Jet_partonFlavour",Jet.partonFlavour);
	
	tree_initial->SetBranchAddress("met_pt",&met_pt);
	tree_initial->SetBranchAddress("met_phi",&met_phi);
	
	tree_initial->SetBranchAddress("softActivityJets_pt",Jet.soft_pt);
	tree_initial->SetBranchAddress("softActivityJets_eta",Jet.soft_eta);
	tree_initial->SetBranchAddress("softActivityJets_mass",Jet.soft_mass);
	tree_initial->SetBranchAddress("softActivity_HT",&Jet.HTsoft);
	tree_initial->SetBranchAddress("softActivity_njets2",&Jet.nsoft2);
	tree_initial->SetBranchAddress("softActivity_njets5",&Jet.nsoft5);
	tree_initial->SetBranchAddress("softActivity_njets10",&Jet.nsoft10);
	tree_initial->SetBranchAddress("softActivityEWK_HT",&Jet.EWKHTsoft);
	tree_initial->SetBranchAddress("softActivityEWK_njets2",&Jet.EWKnsoft2);
	tree_initial->SetBranchAddress("softActivityEWK_njets5",&Jet.EWKnsoft5);
	tree_initial->SetBranchAddress("softActivityEWK_njets10",&Jet.EWKnsoft10);
	tree_initial->SetBranchAddress("nsoftActivityEWKJets",&Jet.nsoftActivityEWKJets);
	tree_initial->SetBranchAddress("softActivityEWKJets_pt",Jet.EWKsoft_pt);
	tree_initial->SetBranchAddress("softActivityEWKJets_eta",Jet.EWKsoft_eta);
	tree_initial->SetBranchAddress("softActivityEWKJets_phi",Jet.EWKsoft_phi);
	tree_initial->SetBranchAddress("softActivityEWKJets_mass",Jet.EWKsoft_mass);
    tree_initial->SetBranchAddress("Jet_qgl",Jet.qgl);
	tree_initial->SetBranchAddress("genWeight",&genweight);
	tree_initial->SetBranchAddress("bTagWeight",&bTagWeight);
	tree_initial->SetBranchAddress("puWeight",&puweight);
	tree_initial->SetBranchAddress("puWeightUp",&puweightUp);
	tree_initial->SetBranchAddress("puWeightDown",&puweightDown);
	//tree_initial->SetBranchAddress("PhiZQ1",&PhiZQ1);
	tree_initial->SetBranchAddress("nPVs",&nPVs);
	tree_initial->SetBranchAddress("Jet_ptd",Jet.ptd);
	tree_initial->SetBranchAddress("Jet_axis2",Jet.axis2);
	tree_initial->SetBranchAddress("Jet_mult",Jet.mult);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",&HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",&HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu22_v",&HLT_IsoMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu22_v",&HLT_IsoTkMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v",&HLT_IsoMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu24_v",&HLT_IsoTkMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v",&HLT_IsoTkMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v",&HLT_Ele27_eta2p1);
	tree_initial->SetBranchAddress("Jet_pt_regVBF",Jet.pt_regVBF);
	tree_initial->SetBranchAddress("json",&JSON);
    
	tree_initial->SetBranchAddress("GenHiggsSisters_pt",GenHiggsSisters.pt);
        tree_initial->SetBranchAddress("GenHiggsSisters_eta",GenHiggsSisters.eta);
        tree_initial->SetBranchAddress("GenHiggsSisters_phi",GenHiggsSisters.phi);
	tree_initial->SetBranchAddress("GenHiggsSisters_mass",GenHiggsSisters.mass);
	tree_initial->SetBranchAddress("nGenVbosons",&nGenVbosons);
	tree_initial->SetBranchAddress("GenVbosons_pt",GenVbosons_pt);
	tree_initial->SetBranchAddress("GenVbosons_pdgId",GenVbosons_pdgId);
	tree_initial->SetBranchAddress("VtypeSim",&VtypeSim);

        tree_initial->SetBranchAddress("nGenJet",&nGenJet);
        tree_initial->SetBranchAddress("GenJet_pt",GenJet_pt);
        tree_initial->SetBranchAddress("GenJet_eta",GenJet_eta);
        tree_initial->SetBranchAddress("GenJet_phi",GenJet_phi);
        tree_initial->SetBranchAddress("GenJet_mass",GenJet_mass);
                
        

        
        
        
	tree_initial->SetBranchAddress("nvLeptons",&nvLeptons);
	tree_initial->SetBranchAddress("vLeptons_pt",vLeptons_pt);
	tree_initial->SetBranchAddress("vLeptons_eta",vLeptons_eta);
	tree_initial->SetBranchAddress("vLeptons_phi",vLeptons_phi);
	tree_initial->SetBranchAddress("vLeptons_mass",vLeptons_mass);
	tree_initial->SetBranchAddress("vLeptons_charge",vLeptons_charge);
	tree_initial->SetBranchAddress("vLeptons_pdgId",vLeptons_pdgId);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutLoose",vLeptons_SF_IdCutLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutTight",vLeptons_SF_IdCutTight);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoLoose",vLeptons_SF_IsoLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoTight",vLeptons_SF_IsoTight);
	tree_initial->SetBranchAddress("vLeptons_SF_trk_eta",vLeptons_SF_trk_eta);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p2",vLeptons_SF_HLT_RunD4p2);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p3",vLeptons_SF_HLT_RunD4p3);
	tree_initial->SetBranchAddress("vLeptons_trackerLayers", vLeptons_trackerLayers);





	tree_initial->SetBranchAddress("nselLeptons",&nselLeptons);
	tree_initial->SetBranchAddress("selLeptons_pt",selLeptons_pt);
	tree_initial->SetBranchAddress("selLeptons_eta",selLeptons_eta);
	tree_initial->SetBranchAddress("selLeptons_phi",selLeptons_phi);
	tree_initial->SetBranchAddress("selLeptons_mass",selLeptons_mass);
	tree_initial->SetBranchAddress("selLeptons_charge",selLeptons_charge);
	tree_initial->SetBranchAddress("selLeptons_pdgId",selLeptons_pdgId);
	tree_initial->SetBranchAddress("selLeptons_looseIdPOG",selLeptons_looseIdPOG);
	tree_initial->SetBranchAddress("selLeptons_relIso04",selLeptons_relIso04);
	tree_initial->SetBranchAddress("selLeptons_relIso03",selLeptons_relIso03);
    tree_initial->SetBranchAddress("selLeptons_eleMVAIdSppring16GenPurp",selLeptons_eleMVAIdSppring16GenPurp); 
	tree_initial->SetBranchAddress("selLeptons_trackerLayers", selLeptons_trackerLayers);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleHoE",selLeptons_eleHoE);
	tree_initial->SetBranchAddress("selLeptons_eleDEta",selLeptons_eleDEta);
	tree_initial->SetBranchAddress("selLeptons_eleDPhi",selLeptons_eleDPhi);
	tree_initial->SetBranchAddress("selLeptons_eleEcalClusterIso",selLeptons_eleEcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_eleHcalClusterIso",selLeptons_eleHcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_dr03TkSumPt",selLeptons_dr03TkSumPt);
	
	tree_initial->SetBranchAddress("evt",&evt);
	tree_initial->SetBranchAddress("run",&run);
	tree_initial->SetBranchAddress("lumi",&luminosityBlock);


	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);
	tree_initial->SetBranchAddress("lheHT",&lheHT);
	tree_initial->SetBranchAddress("lheNj",&lheNj_f);
	tree_initial->SetBranchAddress("lheV_pt",&lheV_pt);
 	tree_initial->SetBranchAddress("hlheNpNLO",&lheNpNLO);       
        


    }



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

	


    
    TH1D * hIsolatedElectrons = new TH1D("hIsolatedElectrons","", 4, 0.,4.);
    hIsolatedElectrons->GetXaxis()->SetTitle("#iso Electrons"); 
    
    
    TH1D *hdeltaR1 = new TH1D("hdeltaR1","", 30, 0.,5.);
    hdeltaR1->GetXaxis()->SetTitle("#Delta R (lead reco jet, genJet)");
    TH1D *hdeltaR2 = new TH1D("hdeltaR2","", 30, 0.,5.);
    hdeltaR2->GetXaxis()->SetTitle("#Delta R (sublead reco jet, genJet)");
    
    
    TH1D *hHiggsSister1_Leading_eta = new TH1D("hHiggsSister1_Leading_eta","", 20, 0.,1.);
    hHiggsSister1_Leading_eta->GetXaxis()->SetTitle("SelectionCut");
    TH1D *hHiggsSister2_Subleading_eta = new TH1D("hHiggsSister2_Subleading_eta","", 20, 0.,1.);
    hHiggsSister2_Subleading_eta->GetXaxis()->SetTitle("SelectionCut");
    
        TH1D *hHiggsSister1_Leading_phi = new TH1D("hHiggsSister1_Leading_phi","", 20, 0.,3.5);
    hHiggsSister1_Leading_phi->GetXaxis()->SetTitle("SelectionCut");
    TH1D *hHiggsSister2_Subleading_phi = new TH1D("hHiggsSister2_Subleading_phi","", 20, 0.,3.5);
    hHiggsSister2_Subleading_phi->GetXaxis()->SetTitle("SelectionCut");
    
    TH1D *hHiggsSister1_Leading_R = new TH1D("hHiggsSister1_Leading_R","", 20, 0.,1.);
    hHiggsSister1_Leading_R->GetXaxis()->SetTitle("#DeltaR(j1, q1");
    TH1D *hHiggsSister2_Subleading_R = new TH1D("hHiggsSister2_Subleading_R","", 20, 0.,1.);
    hHiggsSister2_Subleading_R->GetXaxis()->SetTitle("#DeltaR(j2, q2");
	
    TH1D *hHiggsSister1_HiggsSister2_R = new TH1D("hHiggsSister1_HiggsSister2_R","", 20, 0.,10);
    hHiggsSister1_HiggsSister2_R->GetXaxis()->SetTitle(" #DeltaR(q1, q2)");
	
    TH1D *hSelectionCuts = new TH1D("hSelectionCuts","", 20, -0.5,19.5);
    hSelectionCuts->GetXaxis()->SetTitle("SelectionCut");
    	
    TH1D *hSelectionCuts2 = new TH1D("hSelectionCuts2","", 20, -0.5,19.5);
    hSelectionCuts2->GetXaxis()->SetTitle("SelectionCut double");
    
    TH1D *hlheNpNLO = new TH1D("hlheNpNLO","", 6, -0.5,5.5);
    hlheNpNLO->GetXaxis()->SetTitle("hlheNpNLO"); 

    TH1D *hVtype = new TH1D("hVtype","", 7,-1.,6.);
    hVtype->GetXaxis()->SetTitle("vtype");

    TH1D *hVtypeSim = new TH1D("hVtypeSim","", 21,-10.5,10.5);
    hVtypeSim->GetXaxis()->SetTitle("vtypeSim");
    
    TH1D *hmumujj_pt = new TH1D("hmumujj_pt","", 100,0.,500);
    hmumujj_pt->GetXaxis()->SetTitle("pT(#mu#mu jj)");
    TH1D *hmumujj_ptLog = new TH1D("hmumujj_ptLog","", 20,0.,10);
    hmumujj_ptLog->GetXaxis()->SetTitle("log(pT(#mu#mu jj))");
    
    
	TH1D *hMqq = new TH1D("hMqq","",60,0.,3000.);
// 	TH1D *hMqq = new TH1D("hMqq","",50,0.,1000.);
	hMqq->GetXaxis()->SetTitle("m(jj) [GeV]");
        TH1D *hMqq_cut06 = new TH1D("hMqq_cut06","",60,0.,3000.);
	hMqq_cut06->GetXaxis()->SetTitle("m(jj) [GeV]");
        TH1D *hMqq_cut08 = new TH1D("hMqq_cut08","",60,0.,3000.);
	hMqq_cut08->GetXaxis()->SetTitle("m(jj) [GeV]");
        TH1D *hMqq_cut10 = new TH1D("hMqq_cut10","",60,0.,3000.);
	hMqq_cut10->GetXaxis()->SetTitle("m(jj) [GeV]");
        TH1D *hMqq_cut12 = new TH1D("hMqq_cut12","",60,0.,3000.);
	hMqq_cut12->GetXaxis()->SetTitle("m(jj) [GeV]");
        
        
        
	TH1D *hMqq_log = new TH1D("hMqq_log","",40.,5.,9.);
	hMqq_log->GetXaxis()->SetTitle("ln(m(jj)) [GeV]");
	TH1D *hqq_pt = new TH1D("hqq_pt","",40.,0.,500.);
	hqq_pt->GetXaxis()->SetTitle("p_{T}(jj) [GeV]");
	TH1D *hlepton1_pt = new TH1D("hlepton1_pt","",40.,0.,400.);
	hlepton1_pt->GetXaxis()->SetTitle("leading lepton p_{T} [GeV]");
	TH1D *hlepton2_pt = new TH1D("hlepton2_pt","",30.,0.,300.);
	hlepton2_pt->GetXaxis()->SetTitle("subleading lepton p_{T} [GeV]");
	TH1D *hlepton1_eta = new TH1D("hlepton1_eta","",80.,-3.,3.);
	hlepton1_eta->GetXaxis()->SetTitle("leading lepton #eta");
	TH1D *hlepton2_eta = new TH1D("hlepton2_eta","",80.,-3.,3.);
	hlepton2_eta->GetXaxis()->SetTitle("subleading lepton #eta");
  

	TH1D *hlepton1_iso03 = new TH1D("hlepton1_iso03","",40.,0.,0.3);
	hlepton1_iso03->GetXaxis()->SetTitle("leading lepton iso");
	
	TH1D *hlepton2_iso03= new TH1D("hlepton2_iso03","",40.,0.,0.3);
	hlepton2_iso03->GetXaxis()->SetTitle("subleading lepton iso");


    
    
 
	TH1D *hEtaQQ = new TH1D("hEtaQQ","",54,0.,9.);
	hEtaQQ->GetXaxis()->SetTitle("|#Delta#eta_{jj}|");
	
	TH1D *hbdt = new TH1D("hbdt","",100,-1.,1.);
	hbdt->GetXaxis()->SetTitle("BDT output");
	TH1D *hbdt_atanh = new TH1D("hbdt_atanh","",500,0.,5.);
	hbdt_atanh->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	float bining[50];
	bining[0]=0.;
	for (int i=1;i<31;i++)
		bining[i]=bining[i-1]+0.1;
	bining[31] = 3.5;
			
	TH1D *hbdt_atanh2 = new TH1D("hbdt_atanh2","",31,bining);
	hbdt_atanh2->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	TH1D *hPhiQQ = new TH1D("hPhiQQ","",32,0.,3.2);
	hPhiQQ->GetXaxis()->SetTitle("|#Delta#phi_{jj}|");
    
	
	TH1D *hEtaSoftJets = new TH1D("hEtaSoftJets","",12,-3.,3.);
	hEtaSoftJets->GetXaxis()->SetTitle("|#eta^{soft}|");
    
	
	TH1D *hMassSoftJets = new TH1D("hMassSoftJets","",10,0.,100.);
	hMassSoftJets->GetXaxis()->SetTitle("m^{soft}");


	TH1D *hSoft_n2 = new TH1D("hSoft_n2","",25,0.,25.);
	hSoft_n2->GetXaxis()->SetTitle("N soft jets, p_{T} > 2 GeV");
	TH1D *hSoft_n5 = new TH1D("hSoft_n5","",10,0.,10.);
	hSoft_n5->GetXaxis()->SetTitle("N soft jets, p_{T} > 5 GeV");
	TH1D *hSoft_n10 = new TH1D("hSoft_n10","",6,0.,6.);
	hSoft_n10->GetXaxis()->SetTitle("N soft jets, p_{T} > 10 GeV");

	TH1D *hHTsoftEWK = new TH1D("hHTsoftEWK","",30,0.,300.);
	hHTsoftEWK->GetXaxis()->SetTitle("EWK H_{T}^{soft} (GeV)" );
	TH1D *hSoft_n2EWK = new TH1D("hSoft_n2EWK","",25,0.,25.);
	hSoft_n2EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV");
	TH1D *hSoft_n5EWK = new TH1D("hSoft_n5EWK","",10,0.,10.);
	hSoft_n5EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV");
	TH1D *hSoft_n10EWK = new TH1D("hSoft_n10EWK","",6,0.,6.);
	hSoft_n10EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV");
	
	TH1D *hHTsoftEWK_bdt = new TH1D("hHTsoftEWK_bdt","",30,0.,300.);
	hHTsoftEWK_bdt->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.92 (GeV)" );
	TH1D *hSoft_n2EWK_bdt = new TH1D("hSoft_n2EWK_bdt","",25,0.,25.);
	hSoft_n2EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.92");
	TH1D *hSoft_n5EWK_bdt = new TH1D("hSoft_n5EWK_bdt","",10,0.,10.);
	hSoft_n5EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.92");
	TH1D *hSoft_n10EWK_bdt = new TH1D("hSoft_n10EWK_bdt","",6,0.,6.);
	hSoft_n10EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.92");


	TH1D *hHTsoftEWK_bdt2 = new TH1D("hHTsoftEWK_bdt2","",30,0.,300.);
	hHTsoftEWK_bdt2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.84 (GeV)" );
	TH1D *hSoft_n2EWK_bdt2 = new TH1D("hSoft_n2EWK_bdt2","",25,0.,25.);
	hSoft_n2EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.84");
	TH1D *hSoft_n5EWK_bdt2 = new TH1D("hSoft_n5EWK_bdt2","",10,0.,10.);
	hSoft_n5EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.84");
	TH1D *hSoft_n10EWK_bdt2 = new TH1D("hSoft_n10EWK_bdt2","",6,0.,6.);
	hSoft_n10EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.84");

// 	TH1D *hHTsoftEWK_mjj1 = new TH1D("hHTsoftEWK_mjj1","",30,0.,300.);
// 	hHTsoftEWK_mjj1->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 1500 (GeV)" );
// 	TH1D *hSoft_n2EWK_mjj1 = new TH1D("hSoft_n2EWK_mjj1","",25,0.,25.);
// 	hSoft_n2EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 1500");
// 	TH1D *hSoft_n5EWK_mjj1 = new TH1D("hSoft_n5EWK_mjj1","",10,0.,10.);
// 	hSoft_n5EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 1500");
// 	TH1D *hSoft_n10EWK_mjj1 = new TH1D("hSoft_n10EWK_mjj1","",6,0.,6.);
// 	hSoft_n10EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 1500");
// 	
// 	TH1D *hHTsoftEWK_mjj2 = new TH1D("hHTsoftEWK_mjj2","",30,0.,300.);
// 	hHTsoftEWK_mjj2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 2500 (GeV)" );
// 	TH1D *hSoft_n2EWK_mjj2 = new TH1D("hSoft_n2EWK_mjj2","",25,0.,25.);
// 	hSoft_n2EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 2500");
// 	TH1D *hSoft_n5EWK_mjj2 = new TH1D("hSoft_n5EWK_mjj2","",10,0.,10.);
// 	hSoft_n5EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 2500");
// 	TH1D *hSoft_n10EWK_mjj2 = new TH1D("hSoft_n10EWK_mjj2","",6,0.,6.);
// 	hSoft_n10EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 2500");



        TH1D* hqgl = new TH1D("hqgl","",20.,0.,1.);
        hqgl->GetXaxis()->SetTitle("QGL 1^{st} q-jet");

        TH1D* hqgl2 = new TH1D("hqgl2","",20.,0.,1.);
        hqgl2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");
        
        
        TH1D* hqgl_noQGLweight = new TH1D("hqgl_noQGLweight","",20.,0.,1.);
        hqgl_noQGLweight->GetXaxis()->SetTitle("QGL 1^{st} q-jet");
        TH1D* hqgl_noQGLweight2 = new TH1D("hqgl_noQGLweight2","",20.,0.,1.);
        hqgl_noQGLweight2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");
        
        TH1D* hqglAtanh = new TH1D("hqglAtanh","",42.,-7.,7.);
        hqglAtanh->GetXaxis()->SetTitle("tanh^{-1}(QGL 1^{st} q-jet)");

        TH1D* hqgl2Atanh = new TH1D("hqgl2Atanh","",42.,-7.,7.);
        hqgl2Atanh->GetXaxis()->SetTitle("tanh^{-1}(QGL 2^{nd} q-jet)");
	
	TH1D *hPtSoftJets = new TH1D("hPtSoftJets","",30,0.,300);
	hPtSoftJets->GetXaxis()->SetTitle("p_{T}^{soft} (GeV)");
        TH1D *hPtSoftJets2 = new TH1D("hPtSoftJets2", "", 20, 0., 200.);
        hPtSoftJets2->GetXaxis()->SetTitle("2nd Soft Jet p_{T} (GeV)");
        TH1D *hPtSoftJets3 = new TH1D("hPtSoftJets3", "", 50, 0., 200.);
        hPtSoftJets3->GetXaxis()->SetTitle("3rd Soft Jet p_{T} (GeV)");
	
	TH1D *hcosOqqbb = new TH1D("hcosOqqbb","",100,-1.,1.);
	hcosOqqbb->GetXaxis()->SetTitle("cos(#theta_{bb_qq})");
	TH1D *hEtaQB1 = new TH1D("hEtaQB1","",160.,-8.,8.);
	hEtaQB1->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}");
	TH1D *hEtaQB2 = new TH1D("hEtaQB2","",160.,-8.,8.);
	hEtaQB2->GetXaxis()->SetTitle("#Delta#eta_{qb}^{backward}");
	TH1D *hPhiQB1 = new TH1D("hPhiQB1","",32,0.,3.2);
	hPhiQB1->GetXaxis()->SetTitle("#Delta#phi_{qb}^{forward}");
	TH1D *hPhiQB2 = new TH1D("hPhiQB2","",32,0.,3.2);
	hPhiQB2->GetXaxis()->SetTitle("#Delta#phi_{qb}^{backward}");
	TH1D *hx1 = new TH1D("hx1","",100.,0.,1.);
	hx1->GetXaxis()->SetTitle("x_{1}");
	TH1D *hx2 = new TH1D("hx2","",100.,0.,1.);
	hx2->GetXaxis()->SetTitle("x_{2}");
	TH1D *hVB1_mass = new TH1D("hVB1_mass","",100,0.,1000.);
	hVB1_mass->GetXaxis()->SetTitle("M_{W'_{1}} (GeV)");
	TH1D *hVB2_mass = new TH1D("hVB2_mass","",100.,0.,1000.);
	hVB2_mass->GetXaxis()->SetTitle("M_{W'_{2}} (GeV)");

	TH1D* hEtot = new TH1D("hEtot","",150.,0.,6000.);
	hEtot->GetXaxis()->SetTitle("E^{tot} (GeV)");
	TH1D* hPxtot= new TH1D("hPxtot","",100,-500.,500.);
	hPxtot->GetXaxis()->SetTitle("P_{x}^{tot} (GeV)");
	TH1D* hPytot= new TH1D("hPytot","",100,-500.,500.);
	hPytot->GetXaxis()->SetTitle("P_{y}^{tot} (GeV)");
	TH1D* hPztot= new TH1D("hPztot","",100,-5000.,5000);
	hPztot->GetXaxis()->SetTitle("P_{z}^{tot} (GeV)");

	
	TH1D *hPtqqll = new TH1D("hPtqqll","",50.,0.,500.);
	hPtqqll->GetXaxis()->SetTitle("p_{T} of qqll system (GeV)");
	TH1D *hPhiqqll = new TH1D("hPhiqqll","",32,-3.2,3.2);
	hPhiqqll->GetXaxis()->SetTitle("-#phi of qqll system");
	TH1D *hEtaqqll = new TH1D("hEtaqqll","",160,0,8);
	hEtaqqll->GetXaxis()->SetTitle("#eta of qqll system");

	TH1D *hnPVs = new TH1D("hPVs","",50,0,50);
	hnPVs->GetXaxis()->SetTitle("nPVs");
	
        TH1D *hNjet25 = new TH1D("hNjet25","",10,0,10);
	hNjet25->GetXaxis()->SetTitle("number of jets with pT>25");
	
        
        
    TH1D *hdeltaM=new TH1D("hdeltaM","",20,0,20);
    hdeltaM->GetXaxis()->SetTitle("#Delta M(#mu#mu)");
    TH1D *hdeltaMRel=new TH1D("hdeltaMRel","",20,0,0.2);
    hdeltaMRel->GetXaxis()->SetTitle("#Delta M(#mu#mu)/M(#mu#mu)");
    TH1D *hnormalizedDistance_from_mH=new TH1D("hdnormalizedDistance_from_mH","",40,0,3);
    hnormalizedDistance_from_mH->GetXaxis()->SetTitle("hnormalizedDistance_from_mH");
    
    
        TH1D* hZll_mass_biasUp = new TH1D("hZll_mass_biasUp","",35,110.,145.);
        hZll_mass_biasUp->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");
        TH1D* hZll_mass_biasDown = new TH1D("hZll_mass_biasDown","",35,110.,145.);
        hZll_mass_biasDown->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");
        TH1D* hZll_mass_resolution = new TH1D("hZll_mass_resolution","",35,110.,145.);
        hZll_mass_resolution->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");


        
        
	TH1D* hZll_mass = new TH1D("hZll_mass","",45,110.,155.);
// 	TH1D* hZll_mass = new TH1D("hZll_mass","",75,50.,200.);
	hZll_mass->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");
	TH1D* hZll_pt = new TH1D("hZll_pt","",40,0.,400.);
	hZll_pt->GetXaxis()->SetTitle("p_{T}(#mu#mu) (GeV)");
	TH1D* hZll_eta = new TH1D("hZll_eta","",20,-5,5.);
	hZll_eta->GetXaxis()->SetTitle("#eta(#mu#mu)");
	TH1D* hZll_phi = new TH1D("hZll_phi","",32,-3.2,3.2);
	hZll_phi->GetXaxis()->SetTitle("#phi(#mu#mu)");

        TH1D* hHll_mass_unprecise = new TH1D("hHll_mass_unprecise","",35,110.,145.);
	hHll_mass_unprecise->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");
        TH1D* hHll_mass_precise = new TH1D("hHll_mass_precise","",35,110.,145.);
	hHll_mass_precise->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");

        
	TH1D *hJet1q_pt = new TH1D("hJet1q_pt","",30,0,500.);
	hJet1q_pt->GetXaxis()->SetTitle("p_{T} 1^{st} q-jet");
	TH1D *hJet1q_eta = new TH1D("hJet1q_eta","",20,-5,5);
	hJet1q_eta->GetXaxis()->SetTitle("#eta 1^{st} q-jet");
	TH1D *hJet1q_ptd = new TH1D("hJet1q_ptd","",50,0,1);
	hJet1q_ptd->GetXaxis()->SetTitle("ptd 1^{st} q-jet");
	TH1D *hJet1q_axis2= new TH1D("hJet1q_axis2","",32,0.,0.16);
	hJet1q_axis2->GetXaxis()->SetTitle("#sigma_{2} 1^{st} q-jet");
	TH1D *hJet1q_mult= new TH1D("hJet1q_mult","",30,0,30.);
	hJet1q_mult->GetXaxis()->SetTitle("N 1^{st} q-jet");
	TH1D *hJet1q_leadTrackPt= new TH1D("hJet1q_leadTrackPt","",20,0,100.);
	hJet1q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 1^{st} q-jet");
	TH1D *hJet1q_leadTrackEta= new TH1D("hJet1q_leadTrackEta","",20,-5,5.);
	hJet1q_leadTrackEta->GetXaxis()->SetTitle("leading track #eta 1^{st} q-jet");
	TH1D *hJets12_pt = new TH1D("hJets12_pt","",30,0,600.);
	hJets12_pt->GetXaxis()->SetTitle("|p_{T}(j_{1}) + p_{T}(j_{2})| (GeV)");
	TH1D *hJets12_pt_log = new TH1D("hJets12_pt_log","",60,3,9);
	hJets12_pt_log->GetXaxis()->SetTitle("ln|p_{T}(j_{1}) + p_{T}(j_{2})| (GeV)");
	TH1D *hJet1q_phi = new TH1D("hJet1q_phi","",32,-3.2,3.2);
	hJet1q_phi->GetXaxis()->SetTitle("#phi 1^{st} q-jet");
	TH1D *hJet2q_phi = new TH1D("hJet2q_phi","",32,-3.2,3.2);
	hJet2q_phi->GetXaxis()->SetTitle("#phi 2^{nd} q-jet");


        TH1D *hJet_genJetIdx = new TH1D("hJet_genJetIdx","",40,-25,15);
        hJet_genJetIdx->GetXaxis()->SetTitle("genJet idx1; 10 - genJet idx2");
        
        
        
	TH1D *hJet2q_pt = new TH1D("hJet2q_pt","",30,0,300);
	hJet2q_pt->GetXaxis()->SetTitle("p_{T} 2^{nd} q-jet");
	TH1D *hJet2q_eta = new TH1D("hJet2q_eta","",20,-5,5);
	hJet2q_eta->GetXaxis()->SetTitle("#eta 2^{nd} q-jet");
	TH1D *hJet2q_ptd = new TH1D("hJet2q_ptd","",50,0,1);
	hJet2q_ptd->GetXaxis()->SetTitle("ptd 2^{nd} q-jet");
	TH1D *hJet2q_axis2= new TH1D("hJet2q_axis2","",32,0.,0.16);
	hJet2q_axis2->GetXaxis()->SetTitle("#sigma_{2} 2^{nd} q-jet");

	
	TH1D *hJet2q_pt_log = new TH1D("hJet2q_pt_log","",60,2,7);
	hJet2q_pt_log->GetXaxis()->SetTitle("log p_{T} 2^{nd} q-jet");
	TH1D *hJet1q_pt_log = new TH1D("hJet1q_pt_log","",60,2,7);
	hJet1q_pt_log->GetXaxis()->SetTitle("log p_{T} 1^{st} q-jet");
	


	TH1D *hJet2q_mult= new TH1D("hJet2q_mult","",30,0,30.);
	hJet2q_mult->GetXaxis()->SetTitle("N 2^{nd} q-jet");
	TH1D *hJet2q_leadTrackPt= new TH1D("hJet2q_leadTrackPt","",20,0,100.);
	hJet2q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 2^{nd} q-jet");
	
	TH1D *hJet3_pt = new TH1D("hJet3_pt","",18,20,200);
	hJet3_pt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
        TH1D *hJet3_pt_log = new TH1D("hJet3_pt_log","",60,2,7);
	hJet3_pt_log->GetXaxis()->SetTitle("log(p_{T} 3^{rd} jet)");
        
	TH1D *hJet3_pt_new = new TH1D("hJet3_pt_new","",13,0,195);
	hJet3_pt_new->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
	TH1D *hJet3_eta = new TH1D("hJet3_eta","",20,-5,5);
	hJet3_eta->GetXaxis()->SetTitle("#eta 3^{rd} jet");
	TH1D *hJet3_eta_bdt = new TH1D("hJet3_eta_bdt","",20,-5,5);
	hJet3_eta_bdt->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.92");
	TH1D *hJet3_eta_bdt2 = new TH1D("hJet3_eta_bdt2","",20,-5,5);
	hJet3_eta_bdt2->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.84");
        
        TH1D *hJet3_etaRatio = new TH1D("hJet3_etaRatio","",20,-5,5);
	hJet3_etaRatio->GetXaxis()->SetTitle("#etaRatio 3^{rd} jet");
	
	TH1D *hsoftleadTrackPt= new TH1D("hsoftleadTrackPt","",40,0,200.);
	hsoftleadTrackPt->GetXaxis()->SetTitle("leading track p_{T}");
	TH1D *hsoftleadTrackEta= new TH1D("hsoftleadTrackEta","",20,-5,5.);
	hsoftleadTrackEta->GetXaxis()->SetTitle("leading s track #eta ");
	
	TH1D *hAdJetHT = new TH1D("hAdJetHT","",30,0,300);
        hAdJetHT->GetXaxis()->SetTitle("additional jets H_{T} (GeV)");

        TH1D *hmet = new TH1D("hmet","",40,0.,250.);
	hmet->GetXaxis()->SetTitle("MET p_{T} (GeV)");
	TH1D *hrho = new TH1D("hrho","",40,0.,40.);
	hrho->GetXaxis()->SetTitle("rho");
	TH1D *hHT = new TH1D("hHT","",50,0.,1000.);
	hHT->GetXaxis()->SetTitle("lhe H_{T} (GeV)" );
	TH1D *hlheHT_log = new TH1D("hlheHT_log","",50,0.,10.);
	hlheHT_log->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );
        TH1D *hlheHT_log_BDTgt1 = new TH1D("hlheHT_log_BDTgt08","",50,0.,10.);
	hlheHT_log_BDTgt1->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );
        TH1D *hlheHT = new TH1D("hlheHT","",20,0.,400.);
	hlheHT->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );
        TH1D *hlheHT_BDTgt1 = new TH1D("hlheHT_BDTgt08","",20,0.,400.);
	hlheHT_BDTgt1->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );
        
	TH1D *hlheNj = new TH1D("hlheNj","",6,0.,6);
	hlheNj->GetXaxis()->SetTitle("lhe N jets" );

        
        
	TH1D *hDeltaRelQQ = new TH1D("hDeltaRelQQ","",25.,0.,1.);
	hDeltaRelQQ->GetXaxis()->SetTitle("#Delta_{rel}(jj)");
// 	TH1D *hRpt = new TH1D("hRpt","",25.,0.,0.4);
	TH1D *hRpt = new TH1D("hRpt","",50.,0.,1);
	hRpt->GetXaxis()->SetTitle("R(p_{T})");
        TH1D *hRptAtanh = new TH1D("hRptAtanh","",35.,-3.5,3.5);
	hRptAtanh->GetXaxis()->SetTitle("tanh^{-1}(R(p_{T}))");
        
	TH1D *hEtaQQSum = new TH1D("hEtaQQSum","",90,0.,9.);
	hEtaQQSum->GetXaxis()->SetTitle("|#eta(j_{1})| + |#eta(j_{2})| ");
	TH1D *hPhiZQ1 = new TH1D("hPhiZQ1","",32,0.,3.2);
	hPhiZQ1->GetXaxis()->SetTitle("|#Delta#phi(#mu#mu,j_{1})|");
        TH1D *hPhiZQ2 = new TH1D("hPhiZQ2","",32,0.,3.2);
	hPhiZQ2->GetXaxis()->SetTitle("|#Delta#phi(#mu#mu,j_{2})|");
        
        TH1D *hEtaHQ1 = new TH1D("hEtaHQ1","",25,0.,10.);
	hEtaHQ1->GetXaxis()->SetTitle("|#Delta#eta(#mu#mu,j_{1})|");
        TH1D *hEtaHQ2 = new TH1D("hEtaHQ2","",25,0.,10.);
	hEtaHQ2->GetXaxis()->SetTitle("|#Delta#eta(#mu#mu,j_{2})|");
        
	TH1D* hZll_y = new TH1D("hZll_y","",20,-4,4.);
	hZll_y->GetXaxis()->SetTitle("y(#mu#mu)");
	TH1D* hZll_ystar = new TH1D("hZll_ystar","",20,-6,6.);
	hZll_ystar->GetXaxis()->SetTitle("y*(#mu#mu)");
	TH1D* hZll_zstar_log = new TH1D("hZll_zstar_log","",20,-8,3.);
	hZll_zstar_log->GetXaxis()->SetTitle("log(z*(#mu#mu))");
        TH1D* hZll_zstar = new TH1D("hZll_zstar","",40,0,3);
	hZll_zstar->GetXaxis()->SetTitle("|z*(#mu#mu)|");
	TH1D* hlheV_pt = new TH1D("hlheV_pt","",40,0.,400.);
	hlheV_pt->GetXaxis()->SetTitle("lheV_pt (GeV)");

        TH1D* hzstar = new TH1D("hzstar","",40,0,2);
	hzstar->GetXaxis()->SetTitle("|z*(#mu#mu)| con 1/2");
        
	TH1D *hJet3_pt_bdt = new TH1D("hJet3_pt_bdt","",13,0,195);
	hJet3_pt_bdt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.92 (GeV)");
	TH1D *hAdJetHT_bdt = new TH1D("hAdJetHT_bdt","",30,0,450);
	hAdJetHT_bdt->GetXaxis()->SetTitle("additional jets HT, BDT>0.92 (GeV)");
	TH1D *hNAdJets_bdt = new TH1D("hNAdJets_bdt","",10,0,10);
	hNAdJets_bdt->GetXaxis()->SetTitle("N of jets, BDT>0.92");
	TH1D *hNAdJets = new TH1D("hNAdJets","",10,0,10);
	hNAdJets->GetXaxis()->SetTitle("N of jets");

	TH1D *hJet3_pt_bdt2 = new TH1D("hJet3_pt_bdt2","",13,0,195);
	hJet3_pt_bdt2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.84 (GeV)");
	TH1D *hAdJetHT_bdt2 = new TH1D("hAdJetHT_bdt2","",30,0,450);
	hAdJetHT_bdt2->GetXaxis()->SetTitle("additional jets HT, BDT>0.84 (GeV)");
	TH1D *hNAdJets_bdt2 = new TH1D("hNAdJets_bdt2","",10,0,10);
	hNAdJets_bdt2->GetXaxis()->SetTitle("N of jets, BDT>0.84");
// 	TH1D *hJet3_pt_mjj1 = new TH1D("hJet3_pt_mjj1","",13,0,195);
// 	hJet3_pt_mjj1->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 1500 (GeV)");
// 	TH1D *hAdJetHT_mjj1 = new TH1D("hAdJetHT_mjj1","",30,0,450);
// 	hAdJetHT_mjj1->GetXaxis()->SetTitle("additional jets HT, m(qq) > 1500 (GeV)");
// 	TH1D *hNAdJets_mjj1 = new TH1D("hNAdJets_mjj1","",10,0,10);
// 	hNAdJets_mjj1->GetXaxis()->SetTitle("N of jets, m(qq) > 1500");
// 	TH1D *hJet3_pt_mjj2 = new TH1D("hJet3_pt_mjj2","",13,0,195);
// 	hJet3_pt_mjj2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 2500 (GeV)");
// 	TH1D *hAdJetHT_mjj2 = new TH1D("hAdJetHT_mjj2","",30,0,450);
// 	hAdJetHT_mjj2->GetXaxis()->SetTitle("additional jets HT, m(qq) > 2500 (GeV)");
// 	TH1D *hNAdJets_mjj2 = new TH1D("hNAdJets_mjj2","",10,0,10);
// 	hNAdJets_mjj2->GetXaxis()->SetTitle("N of jets, m(qq) > 2500 ");


	TH1D *hJet1q_eta_bdt = new TH1D("hJet1q_eta_bdt","",20,-5,5);
	hJet1q_eta_bdt->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.92");
	TH1D *hJet2q_eta_bdt = new TH1D("hJet2q_eta_bdt","",20,-5,5);
	hJet2q_eta_bdt->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.92");
	
	TH1D *hJet1q_eta_bdt2 = new TH1D("hJet1q_eta_bdt2","",20,-5,5);
	hJet1q_eta_bdt2->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.84");
	TH1D *hJet2q_eta_bdt2 = new TH1D("hJet2q_eta_bdt2","",20,-5,5);
	hJet2q_eta_bdt2->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.84");
	
	TH1D *hpdgId= new TH1D("hpdgId","",40,-20,20);
	hpdgId->GetXaxis()->SetTitle("type of particle");
    
    TH1D *hweights= new TH1D("hweights","",300,-1.5,1.5);
    hweights->GetXaxis()->SetTitle("sign(weight)");
	 
    TH1D *hweights_weighted= new TH1D("hweights_weighted","",300,-1.5,1.5);
    hweights_weighted->GetXaxis()->SetTitle("sign(weight)weighted");
    
	TH1D *hveto_jet3pt_nom = new TH1D("hveto_jet3pt_nom","",9,0,270);
	TH1D *hveto_jet3pt_denom = new TH1D("hveto_jet3pt_denom","",9,0,270);
	TH1D *hveto_ht_nom = new TH1D("hveto_ht_nom","",14,0,420);
	TH1D *hveto_ht_denom = new TH1D("hveto_ht_denom","",14,0,420);
	TH1D *hveto_softht_nom = new TH1D("hveto_softht_nom","",8,0,320);
	TH1D *hveto_softht_denom = new TH1D("hveto_softht_denom","",8,0,320);
	TH1D *hveto_softpt_nom = new TH1D("hveto_softpt_nom","",6,0,180);
	TH1D *hveto_softpt_denom = new TH1D("hveto_softpt_denom","",6,0,180);


        
//        float pointWhereBinsChange = 2.2;
//        int binNumber_smallBins = 22;
//        const int binNUmberBDT = 30;
//        float BDT_bin[binNUmberBDT];
//        for (int n = 0; n < binNUmberBDT; ++n ) {
//            if (n<binNumber_smallBins) BDT_bin[n] = n*(pointWhereBinsChange/binNumber_smallBins);
//            else BDT_bin[n] = binNumber_smallBins*(pointWhereBinsChange/binNumber_smallBins) + (n-binNumber_smallBins)*(2.*pointWhereBinsChange/binNumber_smallBins);
//        }



//         const int binNUmberBDT = 28;
//         float BDT_bin[binNUmberBDT];
//         for (int n = 0; n < binNUmberBDT-2; ++n ) BDT_bin[n] = n/10.;
//         BDT_bin[26] = 2.8;
//         BDT_bin[27] = 4.0;
        
//         const int binNUmberBDT = 31;
//         float BDT_bin[binNUmberBDT] = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 4};
        
//         const int binNUmberBDT = 10;
//         float BDT_bin[binNUmberBDT] = {0.,        0.1507,         0.3003,         0.4499,         0.5995,         0.7491,         0.8987,         1.0483,         1.1979,         2.5};
        
//         const int binNUmberBDT = 15;
//         float BDT_bin[binNUmberBDT] ={0,     0.0728,         0.1728,         0.2728,         0.3728,         0.4728,         0.5728,         0.6728,         0.7728,       0.8728,  0.9728,         1.0728,         1.1728,         1.2728,        1.5};

        
//         const int binNUmberBDT = 34;
//         float BDT_bin[binNUmberBDT] = {0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 2.};

        
//         const int binNUmberBDT = 28;
//          float BDT_bin[binNUmberBDT] = {0,     0.4854,         0.6354,         0.7854,         0.9354,         1.0854,         1.2354,         1.3854,         1.5354,       1.6854,  1.8354,         1.9854,         2.1354,         2.2854,         2.4354,         2.5854,         2.7354,         2.8854,       3.0354,  3.1854,         3.3354,         3.4854,         3.6354,         3.7854,         3.9354,         4.0854,         4.2354,       6}; //MLP
         
         
//          const int binNUmberBDT = 32;
//         float BDT_bin[binNUmberBDT] = {0,     0.1692,         0.4692,         0.7692,         1.0692,         1.3692,         1.6692,         1.9692,         2.2692,         2.5692,  2.8692,  3.1692,         3.4692,         3.7692,         4.0692,         4.3692,         4.6692,         4.9692,         5.2692,         5.5692,  5.8692,  6.1692,         6.4692,         6.7692,         7.0692,         7.3692,         7.6692,         7.9692,         8.2716,         8.5716,  9.1092,  10};   //Sum Of BDT and NN
        

        
        
//                   const int binNUmberBDT = 10;
//         float BDT_bin[binNUmberBDT] = {0, 0.2628,  0.4128,         0.5628,         0.7128,         0.8628,         1.0128,         1.1628,         1.3128,         1.5};



//  const int binNUmberBDT = 10;
//         float BDT_bin[binNUmberBDT] = {0,     0.18855,         0.33855,        0.48855,        0.63855,        0.78855,        0.93855,        1.08855,        1.23855,        1.5};
 
//         const int binNUmberBDT = 14;
//         float BDT_bin[binNUmberBDT] = {0,     0.0732,         0.1731,         0.273,  0.3729,         0.4728,         0.5727,         0.6726,         0.7725,         0.8724,        0.9723,         1.0722,         1.1721,         1.5};
        
//         const int binNUmberBDT = 10;
//         float BDT_bin[binNUmberBDT] = {0, 0.17175,        0.32175,        0.47175,        0.62175,        0.77175,        0.92175,      1.07175,         1.22175,        1.5};
        
//         const int binNUmberBDT = 10;
//         float BDT_bin[binNUmberBDT] = {0,         0.16485,        0.31485,        0.46485,        0.61485,        0.76485,        0.91485,        1.06485,     1.21485,         2.5};
        
        
//         const int binNUmberBDT = 16;   // Mqq > 250. This is the binning I user in Agese thesis
//         float BDT_bin[binNUmberBDT] = {0,     0.0225,         0.1224,         0.2223,         0.3222,         0.4221,         0.522,  0.6219,         0.7218,         0.8217,        0.9216,         1.0215,         1.1214,         1.2213,         1.3842,         2.5};
        
        
//         const int binNUmberBDT = 16;// Mqq > 200.   and Mqq > 250.   and Mqq > 300.   and Mqq>400 and Mqq > 600
//         float BDT_bin[binNUmberBDT] = {0,     0.0225,         0.1224,         0.2223,    0.3222,    0.4221,     0.522,  0.6219,      0.7218,         0.8217,         0.9216,   1.0215,         1.1214,         1.2213,    1.3842,  1.6};
                                   
        
//         const int binNUmberBDT = 16;// nanoAOD files with old training
//         float BDT_bin[binNUmberBDT] = {0,     0.063,  0.163,  0.263,  0.363,  0.463,  0.563,  0.663,  0.763,  0.863,  0.963,  1.063,  1.163,  1.263,  1.363,  1.6};
        
//         const int binNUmberBDT = 15;// nanoAOD files with old training
//         float BDT_bin[binNUmberBDT] = {0,     0.0725,  0.1725,   0.2725,   0.3725,  0.4725,   0.5725,   0.6725,  0.7725,  0.8725,  0.9725,  1.0725,  1.1725,  1.2725, 1.6};
                
//         const int binNUmberBDT = 15;// nanoAOD files with old training with EKW_LLJJ and EKW_LLJJ_INT
//         float BDT_bin[binNUmberBDT] = {0,     0.0955,  0.1955, 0.2955, 0.3955,   0.4955,  0.5955, 0.6955, 0.7955, 0.8955, 0.9955, 1.0955,  1.1955,  1.2955,   1.6 };
        
        
        
//         const int binNUmberBDT = 22;// old BDT with nano and new samples, last bins are splitted
//         float BDT_bin[binNUmberBDT] = {0,     0.0955,  0.1955, 0.2955, 0.3955,   0.4955,  0.5955, 0.6955, 0.7455, 0.7955, 0.8455, 0.8955, 0.9455, 0.9955, 1.0499, 1.0955, 1.1499,  1.1955, 1.2499,  1.2955, 1.30,  1.6 };
        
        

        
        const int binNUmberBDT_splitted = 22;// first BDT in october Training, last bins are splitted
        float BDT_bin_splitted[binNUmberBDT_splitted] = {0,  0.017,  0.117,  0.217,  0.317,  0.417,  0.517,  0.617,  0.717,  0.817,  0.867,  0.917,  0.967,  1.017,  1.067,  1.117,  1.167,  1.217,  1.267,  1.319,  1.4,  1.6 };
//         const int binNUmberBDT = 16;
//         float BDT_bin[binNUmberBDT] = {0,     0.017,  0.117,  0.217,  0.317,  0.417,  0.517,  0.617,  0.717,  0.817,  0.917,  1.017,  1.117,  1.217,  1.319,  1.6 };
            
        
//         const int binNUmberBDT = 15;// first BDT in october Training, binning for DY_LO + VBF Z 
//         float BDT_bin[binNUmberBDT] = {0, 0.0205,  0.1205,  0.2205,  0.3205,  0.4205, 0.5205, 0.6205, 0.7205,  0.8205, 0.9205,  1.0205, 1.1205,  1.2205,   1.6};
//         {0,         0.22125,        0.37125,        0.52125,        0.67125,        0.82125,        0.97125,        1.12125,     1.27125,         2.5};
//         float BDT_bin[binNUmberBDT] = {0,  0.1,  0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9,    1,      1.1,    1.2,    1.3,    1.4,    1.581,  1.8};

//         const int binNUmberBDT = 16;// first BDT in october Training, binning for DY_LO + VBF Z 
//         float BDT_bin[binNUmberBDT] = {0,  0.037,  0.137,  0.237,  0.337,  0.437,  0.537,  0.637,  0.737,  0.837,  0.937,  1.037,  1.137,  1.237,  1.337,  1.6};
            
            
/*        const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704
        float BDT_bin[binNUmberBDT] = {0,     0.048,  0.148,  0.248,  0.348,  0.448,  0.548,  0.648,  0.748,  0.848,  0.948,  1.048,  1.2095, 1.3195,   1.6};  */   
//         const int binNUmberBDT = 25;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.05)
//         float BDT_bin[binNUmberBDT] = {0, 0.048, 0.098,  0.148,  0.198,  0.248,  0.298,  0.348,  0.398,  0.448,  0.498,  0.548,  0.598,  0.648,  0.698,  0.748,  0.798,  0.848,  0.898,  0.948,   0.998,  1.048,  1.2095,  1.3195,  1.6};     
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.1) Mqq>400
//         float BDT_bin[binNUmberBDT] = {0,     0.048,  0.148,  0.248,  0.348,  0.448,  0.548,  0.648,  0.748,  0.848,  0.948,  1.048,  1.2095, 1.3195, 1.6};
            
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704  no qgl_2qAtanh   (binning  width = 0.1) Mqq>400
//         float BDT_bin[binNUmberBDT] = {0,     0.026,  0.126,  0.226,  0.326,  0.426,  0.526,  0.626,  0.726,  0.826,  0.926,  1.026,  1.126,  1.317, 1.6};

        
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.1) Mqq>400
//         float BDT_bin[binNUmberBDT] = { 0, 0.1435, 0.2435, 0.3435, 0.4435,   0.5435,  0.6435, 0.7435, 0.8435, 0.9435,  1.0435, 1.1435, 1.2435,  1.34525, 1.6};
        
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.1) Mqq>400
//         float BDT_bin[binNUmberBDT] = {0., 0.12475,  0.22475, 0.32475,  0.42475, 0.52475, 0.62475,  0.72475,  0.82475,  0.92475,  1.02475,  1.12475,  1.22475,  1.32475,  1.6};
        
        
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.1) Mqq>400    nuovo skimming
//         float BDT_bin[binNUmberBDT] = {0,     0.0415,  0.1415,  0.2415,  0.3415,  0.4415,   0.5415,   0.6415, 0.7415,  0.8415,    0.9415,     1.0415,     1.1415,  1.288,  1.6};
        
        
        
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_NLO + DY_VBFFilter_NLO + VBF Z,       BDT 1704     (binning  width = 0.1) Mqq>400     skimming NanoAODv3
//         float BDT_bin[binNUmberBDT] = {0, 0.0095,  0.1095,   0.2095,  0.3095,  0.4095,  0.5095,  0.6095,  0.7095, 0.8095,   0.9095,  1.0095, 1.1095,   1.2095,  1.6};
        
        
//         const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + DY_VBFFilter_NLO + VBF Z, computed with macroBDT... on madgraphMLM and....       BDT 1704     (binning  width = 0.1) Mqq>400     skimming NanoAODv3
//         float BDT_bin[binNUmberBDT] = {0,     0.0945,   0.1945,  0.2945,   0.3945,   0.4945,  0.5945,   0.6945,   0.7945,   0.8945,   0.9945, 1.0945,   1.1945,   1.3385, 1.6};
            
        const int binNUmberBDT = 14;//  BDT with DY_VBFFilter, binning for DY_LO + DY_VBFFilter_NLO + VBF Z, computed with macroBDT... on madgraphMLM and....  selection with puId > 6 or |eta|<2.5     BDT 1704     (binning  width = 0.1) Mqq>400     skimming NanoAODv3
        float BDT_bin[binNUmberBDT] = {0,     0.0945,   0.1945,  0.2945,   0.3945,  0.4945,  0.5945, 0.6945,  0.7945, 0.8945, 0.9945,  1.0945, 1.1945, 1.6};
        
            
//         const int binNUmberBDT = 14;//  BDT with DY_VBFFilter, binning for DY_NLO + VBF Z,       BDT 1704     (binning  width = 0.1) Mqq>400     skimming NanoAODv3
//         float BDT_bin[binNUmberBDT] = {{0,     0.0905,         0.1905,         0.2905,  0.3905,   0.4905,    0.5905,  0.6905,   0.7905,   0.8905,    0.9905,    1.0905,     1.1905, 1.6};
        
        
        
//         BINNING TO VALIDATION WITHOUT VBF Z 
//         const int binNUmberBDT = 14;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1704     (binning  width = 0.1) Mqq>400
//         float BDT_bin[binNUmberBDT] ={0,     0.0165,         0.1165,         0.2165,         0.3165,         0.4165,         0.5165,         0.6165,         0.7165,         0.8165,         0.9165,         1.0165,     1.305, 1.6 };
        
        
/*            const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_LO + VBF Z     BDT 1782     (binning  width = 0.1)
            float BDT_bin[binNUmberBDT] = {0, 0.141, 0.241, 0.341, 0.441,  0.541,  0.641,  0.741,  0.841,  0.941,  1.041,  1.141,  1.241,  1.341,  1.6};  */   
        
        
            
//             const int binNUmberBDT = 16;//  BDT with DY_VBFFilter, binning for DY_NLO + VBF Z +   VBF Z INT   BDT 1782     (binning  width = 0.1)     pt subleading > 30 in selection
//             float BDT_bin[binNUmberBDT] = {0,     0.035,  0.135,  0.235,  0.335,  0.435,  0.535,  0.635,  0.735,  0.835,  0.935,  1.035,  1.135,  1.244,  1.384,  1.6};
                
//             const int binNUmberBDT = 15;//  BDT with DY_VBFFilter, binning for DY_NLO + VBF Z +   VBF Z INT   BDT 1782     (binning  width = 0.1)     pt subleading > 25 in selection
//             float BDT_bin[binNUmberBDT] =  {0, 0.0645,   0.1645,   0.2645,  0.3645,  0.4645,  0.5645,   0.6645,   0.7645,  0.8645,  0.9645, 1.0645, 1.1645,  1.3855,  1.6};




//         TEST TRAIN % VARIABLES WITHOUT SOFTN5
//         const int binNUmberBDT = 17;//       BDT 1588 + zStar
//         float BDT_bin[binNUmberBDT] = {0., 0.1245, 0.2245, 0.3245, 0.4245, 0.5245, 0.6245, 0.7245, 0.8245,  0.9245,  1.0245,   1.1245, 1.2245,    1.3245,  1.4245,  1.5245, 1.8};
        
        
        
//         const int binNUmberBDT = 15; // Mqq > 1000
//         float BDT_bin[binNUmberBDT] = {0,     0.108,  0.2079,         0.3078,         0.4077,         0.5076,         0.6075,         0.7074,         0.8073,         0.9072,         1.0071,   1.107,  1.2069,         1.3842,     2.};

        
        
     TH1D *hBDT_VBF_atanh_findBinning = new TH1D("hBDT_VBF_atanh_findBinning","",10000, 0, 5);
     hBDT_VBF_atanh_findBinning->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
     TH1D * hBDT_VBF_atanh_FITTING = new TH1D("hBDT_VBF_atanh_FITTING","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_FITTING->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
       
       
    TH1D *hBDT_VBF_atanh = new TH1D("hBDT_VBF_atanh","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
   

    
    
    TH1D *hBDT_VBF_atanh_m125ForAll = new TH1D("hBDT_VBF_atanh_m125ForAll","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125ForAll->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
     
    
    TH1D *hBDT_VBF_atanh_m125ControlRegion = new TH1D("hBDT_VBF_atanh_m125ControlRegion","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125ControlRegion->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125ControlRegionDown = new TH1D("hBDT_VBF_atanh_m125ControlRegionDown","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125ControlRegionDown->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125ControlRegionUp = new TH1D("hBDT_VBF_atanh_m125ControlRegionUp","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125ControlRegionUp->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    
    
    
    TH1D *hBDT_VBF_atanh_m125_MqqLog75 = new TH1D("hBDT_VBF_atanh_m125_MqqLog75","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_MqqLog75->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_Rpt0 = new TH1D("hBDT_VBF_atanh_m125_Rpt0","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_Rpt0->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_WM1Log5 = new TH1D("hBDT_VBF_atanh_m125_WM1Log5","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_WM1Log5->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_zStar0 = new TH1D("hBDT_VBF_atanh_m125_zStar0","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_zStar0->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_ptll300 = new TH1D("hBDT_VBF_atanh_m125_ptll300","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_ptll300->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_softN50 = new TH1D("hBDT_VBF_atanh_m125_softN50","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_softN50->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_mumujjPt = new TH1D("hBDT_VBF_atanh_m125_mumujjPt","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_mumujjPt->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_DEtajj = new TH1D("hBDT_VBF_atanh_m125_DEtajj","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_DEtajj->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_Theta2 = new TH1D("hBDT_VBF_atanh_m125_Theta2","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_Theta2->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");

    
    
    TH1D *hBDT_VBF_atanh_m125_onlyMqq = new TH1D("hBDT_VBF_atanh_m125_onlyMqq","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyMqq->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyRpt = new TH1D("hBDT_VBF_atanh_m125_onlyRpt","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyRpt->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyWM1 = new TH1D("hBDT_VBF_atanh_m125_onlyWM1","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyWM1->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyzStar = new TH1D("hBDT_VBF_atanh_m125_onlyzStar","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyzStar->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyptll = new TH1D("hBDT_VBF_atanh_m125_onlyptll","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyptll->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlysoftN5 = new TH1D("hBDT_VBF_atanh_m125_onlysoftN5","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlysoftN5->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlymumujjPt = new TH1D("hBDT_VBF_atanh_m125_onlymumujjPt","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlymumujjPt->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyDEtajj = new TH1D("hBDT_VBF_atanh_m125_onlyDEtajj","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyDEtajj->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_m125_onlyTheta2 = new TH1D("hBDT_VBF_atanh_m125_onlyTheta2","",binNUmberBDT_splitted-1, BDT_bin_splitted);
    hBDT_VBF_atanh_m125_onlyTheta2->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    
    

            
    TH1D *hBDT_VBF_atanh_biasUp = new TH1D("hBDT_VBF_atanh_biasUp","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_biasUp->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_biasDown = new TH1D("hBDT_VBF_atanh_biasDown","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_biasDown->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_resolutionUp = new TH1D("hBDT_VBF_atanh_resolutionUp","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_resolutionUp->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_resolutionDown = new TH1D("hBDT_VBF_atanh_resolutionDown","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_resolutionDown->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_QGLUp = new TH1D("hBDT_VBF_atanh_QGLUp","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_QGLUp->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_QGLDown = new TH1D("hBDT_VBF_atanh_QGLDown","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_QGLDown->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
        
        
    TH1D *hBDT_VBF_atanh_genMqqGt350 = new TH1D("hBDT_VBF_atanh_genMqqGt350","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_genMqqGt350->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_genMqqLs350 = new TH1D("hBDT_VBF_atanh_genMqqLs350","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_genMqqLs350->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_genMqq0 = new TH1D("hBDT_VBF_atanh_genMqq0","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_genMqq0->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_negativeFraction = new TH1D("hBDT_VBF_atanh_negativeFraction","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_negativeFraction->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_notWeighted = new TH1D("hBDT_VBF_atanh_notWeighted","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_notWeighted->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    
    
    TH1D *hBDT_VBF_atanh_cumulativeLHENpNLO_0 = new TH1D("hBDT_VBF_atanh_cumulativeLHENpNLO_0","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_cumulativeLHENpNLO_0->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_cumulativeLHENpNLO_1 = new TH1D("hBDT_VBF_atanh_cumulativeLHENpNLO_1","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_cumulativeLHENpNLO_1->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    TH1D *hBDT_VBF_atanh_cumulativeLHENpNLO_2 = new TH1D("hBDT_VBF_atanh_cumulativeLHENpNLO_2","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh_cumulativeLHENpNLO_2->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    
    
    
    
    std::vector<TH1D*> hBDT_VBF_atanh_PDFvariation;
    for (int n = 0; n <= binNUmberBDT; n++) {
        TH1D *hBDT_VBF_atanh_PDFvariation_toAdd = new TH1D(("hBDT_VBF_atanh_PDFvariation_toAdd"+std::to_string(n)).c_str(),"",number_LHE_weights_pdf, 0, number_LHE_weights_pdf);
        hBDT_VBF_atanh_PDFvariation_toAdd->GetXaxis()->SetTitle(("tanh^{-1}(( BDT output + 1.)/2.) variation in bin "+std::to_string(n)).c_str());
        hBDT_VBF_atanh_PDFvariation.push_back(hBDT_VBF_atanh_PDFvariation_toAdd);

    }
    
    
    int NTemporaneo = 19;
    float BinningTemporaneo[NTemporaneo] = {-1., -0.8365, -0.7365, -0.6365, -0.5365, -0.4365, -0.3365, -0.2365, -0.1365, -0.0365,     0.0635,         0.1635,         0.2635,         0.3635,         0.4635,         0.5635,         0.6635,         0.7635, 1};
    int bin_new_NumberBDT = 40;
    float upperLimitBDT = 1.;//bin_new_NumberBDT*0.2;
//     TH1D *hBDT_VBF = new TH1D("hBDT_VBF","",bin_new_NumberBDT,-1.,upperLimitBDT);
    TH1D *hBDT_VBF = new TH1D("hBDT_VBF","",NTemporaneo-1.,BinningTemporaneo);
    hBDT_VBF->GetXaxis()->SetTitle(" BDT output  ");
    
    TH1D *hThetaPlanes = new TH1D("hThetaPlanes","",60,-0.1,1.1);
    hThetaPlanes->GetXaxis()->SetTitle("cos(#theta_{(#mu#mu)(jj)})"); 
    TH1D *hThetaPlanesAtanh = new TH1D("hThetaPlanesAtanh","",40,-6.,6.);
    hThetaPlanesAtanh->GetXaxis()->SetTitle("tanh^{-1}(cos(#theta_{(#mu#mu)(jj)}))"); 
    
    TH1D *hThetaStarJet = new TH1D("hThetaStarJet","",60,-0.1,1.1);
    hThetaStarJet->GetXaxis()->SetTitle("cos(#theta_{jj}*)"); 
    TH1D *hThetaStarJetAtanh = new TH1D("hThetaStarJetAtanh","",40,-6.,6.);
    hThetaStarJetAtanh->GetXaxis()->SetTitle("tanh^{-1}(cos(#theta_{jj}*))"); 
    
    TH1D *hThetaStar = new TH1D("hThetaStar","",50,-1.1,1.1);
    hThetaStar->GetXaxis()->SetTitle("cos(#theta*)");
    
    TH1D *hThetaStarAbs = new TH1D("hThetaStarAbs","",50,0,1.1);
    hThetaStarAbs->GetXaxis()->SetTitle("|cos(#theta*)|");


    TH1D *hBDiscriminator_CSV = new TH1D("hBDiscriminator_CSV","",60,-0.3,1.);
    hBDiscriminator_CSV->GetXaxis()->SetTitle("max CSV");
    
    
    
    TH1D *hMaxJetBTagCSV = new TH1D("hMaxJetBTagCSV","",60,-0.3,1.);
    hMaxJetBTagCSV->GetXaxis()->SetTitle("max CSV");
    TH1D *hMaxSecondJetBTagCSV = new TH1D("hMaxSecondJetBTagCSV","",60,-0.3,1.);
    hMaxSecondJetBTagCSV->GetXaxis()->SetTitle("second max CSV");

    TH1D *hMaxJetBTagCMVA = new TH1D("hMaxJetBTagCMVA","",50,-1,1.);
    hMaxJetBTagCMVA->GetXaxis()->SetTitle("max CMVA ");
    TH1D *hMaxSecondJetBTagCMVA = new TH1D("hMaxSecondJetBTagCMVA","",50,-1.,1.);
    hMaxSecondJetBTagCMVA->GetXaxis()->SetTitle("second max CMVA ");
    
    
    TH1D *hBDiscriminator_CMVA = new TH1D("hBDiscriminator_CMVA","",50,-1,1.);
    hBDiscriminator_CMVA->GetXaxis()->SetTitle("tanh^{-1}( max CMVA )");
    
    TH1D * hmaxAbsEta = new TH1D("hmaxAbsEta","",50,0,5.);
    hmaxAbsEta->GetXaxis()->SetTitle("max(|#eta(j1)|, |#eta(j2)|)");
    TH1D * hmaxAbsEta_BDTmFixgt1 = new TH1D("hmaxAbsEta_BDTmFixgt1","",50,0,5.);
    hmaxAbsEta_BDTmFixgt1->GetXaxis()->SetTitle("max(|#eta(j1)|, |#eta(j2)|)");
    TH1D * hmaxAbsEta_BDTmFixgt08 = new TH1D("hmaxAbsEta_BDTmFixgt08","",50,0,5.);
    hmaxAbsEta_BDTmFixgt08->GetXaxis()->SetTitle("max(|#eta(j1)|, |#eta(j2)|)");
    
    TH1D * hminAbsEta = new TH1D("hminAbsEta","",40,0,4.);
    hminAbsEta->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    TH1D * hminAbsEta_cut06 = new TH1D("hminAbsEta_cut06","",40,0,4.);
    hminAbsEta_cut06->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1D * hminAbsEta_cut08 = new TH1D("hminAbsEta_cut08","",40,0,4.);
    hminAbsEta_cut08->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1D * hminAbsEta_cut10 = new TH1D("hminAbsEta_cut10","",40,0,4.);
    hminAbsEta_cut10->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1D * hminAbsEta_cut12 = new TH1D("hminAbsEta_cut12","",40,0,4.);
    hminAbsEta_cut12->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    
    TH1D * hminAbsEta_GEN = new TH1D("hminAbsEta_GEN","",40,0,4.);
    hminAbsEta_GEN->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    ///Agnese
    
    TH1D *hInvariant_Mass=new TH1D("hInvariant_Mass","",60,0,6000);
    hInvariant_Mass->GetXaxis()->SetTitle("M(jj#mu#mu) [GeV]");
    TH1D *hInvariant_MassLog=new TH1D("hInvariant_MassLog","",50,5,10);
    hInvariant_MassLog->GetXaxis()->SetTitle("log(M(jj#mu#mu)) [GeV]");
    
    TH1D *hTotalEnergy = new TH1D("hTotalEnergy","",50,0,6000);
    hTotalEnergy->GetXaxis()->SetTitle("E(jj#mu#mu) [GeV]");
    TH1D *hTotalEnergylog = new TH1D("hTotalEnergylog","",50,4,10);
    hTotalEnergylog->GetXaxis()->SetTitle("E(jj#mu#mu) [GeV]");
    
    TH1D *hPz= new TH1D("hPz","",50,-5000,5000);
    hPz->GetXaxis()->SetTitle("Pz(jj#mu#mu) [GeV]");
    TH1D *hPzAbs= new TH1D("hPzAbs","",50,0,5000);
    hPzAbs->GetXaxis()->SetTitle("|Pz(jj#mu#mu)| [GeV]");
    TH1D *hPzAbsLog= new TH1D("hPzAbsLog","",50,0,10);
    hPzAbsLog->GetXaxis()->SetTitle("log(|Pz(jj#mu#mu)|)");
   
    
    
    TH1D *hEnergy_fraction_Parton1=new TH1D("hEnergy_fraction_Parton1","",100,0,1);
    hEnergy_fraction_Parton1->GetXaxis()->SetTitle("X(parton1)");
    TH1D *hEnergy_fraction_Parton2=new TH1D("hEnergy_fraction_Parton2","",100,0,1);
    hEnergy_fraction_Parton2->GetXaxis()->SetTitle("X(parton2)");
    TH1D *hEnergy_fraction_Parton1_log=new TH1D("hEnergy_fraction_Parton1_log","",100,-6,0);
    hEnergy_fraction_Parton1_log->GetXaxis()->SetTitle("log(X(parton1))");
    TH1D *hEnergy_fraction_Parton2_log=new TH1D("hEnergy_fraction_Parton2_log","",100,-6,0);
    hEnergy_fraction_Parton2_log->GetXaxis()->SetTitle("log(X(parton2))");
    
    TH1D *hVirtual_Wmass1=new TH1D("hVirtual_Wmass1","",70,-600,100);
    hVirtual_Wmass1->GetXaxis()->SetTitle("M(V_{1}) [GeV]");
    TH1D *hVirtual_Wmass2=new TH1D("hVirtual_Wmass2","",70,-600,100);
    hVirtual_Wmass2->GetXaxis()->SetTitle("M(V_{2}) [GeV]");
    
    TH1D *hVirtual_Wmass1_log=new TH1D("hVirtual_Wmass1_log","",60,2,8);
    hVirtual_Wmass1_log->GetXaxis()->SetTitle("log( -M(V_{1} )");
    TH1D *hVirtual_Wmass2_log=new TH1D("hVirtual_Wmass2_log","",60,2,8);
    hVirtual_Wmass2_log->GetXaxis()->SetTitle("log( -M(V_{2} )");   
    TH1D *hWWmass=new TH1D("hWWmass","",55,-200,350);
    hWWmass->GetXaxis()->SetTitle("M(V+V) [GeV]");
    TH1D *hDiffmass=new TH1D("hDiffmass","",55,-350,200);
    hDiffmass->GetXaxis()->SetTitle("M(V+V)-M_ll [GeV]");
    
    TH1D *hVirtual_Pt1=new TH1D("hVirtual_Pt1","",80,0,400);
    hVirtual_Pt1->GetXaxis()->SetTitle("p_{T}(V_{1}) [GeV]");
    TH1D *hVirtual_Pt2=new TH1D("hVirtual_Pt2","",80,0,400);
    hVirtual_Pt2->GetXaxis()->SetTitle("p_{T}(V_{2}) [GeV]");
    
    TH1D *hVirtual_Pt1_log=new TH1D("hVirtual_Pt1_log","",80,2,8);
    hVirtual_Pt1_log->GetXaxis()->SetTitle("log(p_{T}(V_{1}))");
    TH1D *hVirtual_Pt2_log=new TH1D("hVirtual_Pt2_log","",80,2, 8);
    hVirtual_Pt2_log->GetXaxis()->SetTitle("log(p_{T}(V_{2}))");
    
    
    TH1D *hVirtual_eta1=new TH1D("hVirtual_eta1","",25,0.,5);
    hVirtual_eta1->GetXaxis()->SetTitle("#eta (V_{1})");
    TH1D *hVirtual_eta2=new TH1D("hVirtual_eta2","",25,0.,5);
    hVirtual_eta2->GetXaxis()->SetTitle("#eta (V_{2})");
    
    TH1D *hVirtual_phi1=new TH1D("hVirtual_phi1","",50,-5,5);
    hVirtual_phi1->GetXaxis()->SetTitle("#phi (V_{1})");
    TH1D *hVirtual_phi2=new TH1D("hVirtual_phi2","",50,-5,5);
    hVirtual_phi2->GetXaxis()->SetTitle("#phi (V_{2})");
    
    TH1D *hParton_M1=new TH1D("hParton_M1","",50,-5,5);
    hParton_M1->GetXaxis()->SetTitle("M(parton1) [GeV]");
    TH1D *hParton_M2=new TH1D("hParton_M2","",50,-5,5);
    hParton_M2->GetXaxis()->SetTitle("M(parton2) [GeV]");
    
    TH1D *hTheta_HiggsJ1=new TH1D ("hTheta_HiggsJ1","",50,-1.1,1.1);
    hTheta_HiggsJ1->GetXaxis()->SetTitle("cos(#theta)(h,q1)");
    
    TH1D *hTheta_HiggsJ2=new TH1D ("hTheta_HiggsJ2","",50,-1.1,1.1);
    hTheta_HiggsJ2->GetXaxis()->SetTitle("cos(#theta)(h,q2)");
    
    TH1D *hthetastar_W1=new TH1D("hthetastar_W1","",50,-1.1,1.1);
    hthetastar_W1->GetXaxis()->SetTitle("cos(#theta*)(p1,q1)");
    TH1D *hthetastar_W2=new TH1D("hthetastar_W2","",50,-1.1,1.1);
    hthetastar_W2->GetXaxis()->SetTitle("cos(#theta*)(p2,q2)");
    
    TH1D *hthetastar_W2toHW1= new TH1D("hthetastar_W2toHW1","",50,-1.1,1.1);
    hthetastar_W2toHW1->GetXaxis()->SetTitle("cos(#theta*)(h,w1)");
    TH1D *hthetastar_W1toHW2=new TH1D("hthetastar_W1toHW2","",50,-1.1,1.1);
    hthetastar_W1toHW2->GetXaxis()->SetTitle("cos(#theta*)(h,w2)");
    TH1D *hthetastar_HtoWW=new TH1D("hthetastar_HtoWW","",50,-1.1,1.1);
    hthetastar_HtoWW->GetXaxis()->SetTitle("cos(#theta*)(w1,w2)");
     
    TH2F *hbinning_MeeJet= new TH2F ("hbinning_MeeJet","",3,0,3,3,0,3);
    hbinning_MeeJet->GetXaxis()->SetTitle("M_ee binning");
    hbinning_MeeJet->GetYaxis()->SetTitle("lhe_NJ"); 
    
    TH2F *hbinning_MmumuJet= new TH2F ("hbinning_MmumuJet","",3,0,3,3,0,3);
    hbinning_MmumuJet->GetXaxis()->SetTitle("M_mumu binning");
    hbinning_MmumuJet->GetYaxis()->SetTitle("lhe_NJ"); 
    
    TH2F *hbinning_MttJet= new TH2F ("hbinning_MttJet","",3,0,3,3,0,3);
    hbinning_MttJet->GetXaxis()->SetTitle("M_tt binning");
    hbinning_MttJet->GetYaxis()->SetTitle("lhe_NJ");     
    
    TH1D *hgen_mass=new TH1D ("hgen_mass","",100,80,180);  
    hgen_mass->GetXaxis()->SetTitle("GenLeptons Mass [GeV]");
    
    TH1D *hnGenLep=new TH1D ("hnGenLep","",5,0,5);  
    hnGenLep->GetXaxis()->SetTitle("number of genLeptons");
    
    TH1D *hGenLepton_matching1=new TH1D ("hGenLepton_matching1","",50,0,1);  
    hGenLepton_matching1->GetXaxis()->SetTitle("number of genLeptons");
    TH1D *hGenLepton_matching2=new TH1D ("hGenLepton_matching2","",50,0,1);  
    hGenLepton_matching2->GetXaxis()->SetTitle("number of genLeptons");
    
    TH1D *hnGenJet=new TH1D ("hnGenJet","",10,0,10);  
    hnGenJet->GetXaxis()->SetTitle("number of genLeptons");
    TH1D *hnGenJetsWithoutLeptons=new TH1D ("hnGenJetsWithoutLeptons","",10,0,10);  
    hnGenJetsWithoutLeptons->GetXaxis()->SetTitle("number of genLeptons");
     
    TH1D *hpuId_jetWithoutGenJet=new TH1D ("hpuId_jetWithoutGenJet","",10,0,10);  
    hpuId_jetWithoutGenJet->GetXaxis()->SetTitle("puId of leading jets with genJetIdx=-1");
    TH1D *hpuId_AllJets=new TH1D ("hpuId_AllJets","",10,0,10);  
    hpuId_AllJets->GetXaxis()->SetTitle("puId of leading jets");
    
    TH1D *hrecoJetPt_genJetPt1=new TH1D ("hrecoJetPt_genJetPt1","",50,-1,1);  
    hrecoJetPt_genJetPt1->GetXaxis()->SetTitle("number of genLeptons");
    TH1D *hrecoJetPt_genJetPt2=new TH1D ("hrecoJetPt_genJetPt2","",50,-1,1);  
    hrecoJetPt_genJetPt2->GetXaxis()->SetTitle("number of genLeptons");
    
    TH1D *hgenJetMass=new TH1D ("hgenJetMass","",150,0,1500);  
    hgenJetMass->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1D *hgenJetMass_cut06=new TH1D ("hgenJetMass_cut06","",100,0,3000);  
    hgenJetMass_cut06->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1D *hgenJetMass_cut08=new TH1D ("hgenJetMass_cut08","",100,0,3000);  
    hgenJetMass_cut08->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1D *hgenJetMass_cut10=new TH1D ("hgenJetMass_cut10","",100,0,3000);  
    hgenJetMass_cut10->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1D *hgenJetMass_cut12=new TH1D ("hgenJetMass_cut12","",100,0,3000);  
    hgenJetMass_cut12->GetXaxis()->SetTitle("GenJet Mass [GeV]");

    
    TH1D *hgenJetMass_matched=new TH1D ("hgenJetMass_matched","",100,0,500);  
    hgenJetMass_matched->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    
	
	
	TProfile *hprof_htsoft_pu  = new TProfile("hprof_htsoft_pu","",50,0.,50,0.,300.);
	hprof_htsoft_pu->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu->GetYaxis()->SetTitle("<EWK H_{T}^{soft}> (GeV)");

	TProfile *hprof_htsoft_pu_bdt  = new TProfile("hprof_htsoft_pu_bdt","",50,0.,50,0.,300.);
	hprof_htsoft_pu_bdt->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_bdt->GetYaxis()->SetTitle("<EWK H_{T}^{soft}>, BDT > 0.92 (GeV)");

	TProfile *hprof_htsoft_pu_rms  = new TProfile("hprof_htsoft_pu_rms","",50,0.,50,0.,300.,"s");
	hprof_htsoft_pu_rms->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_rms->GetYaxis()->SetTitle("<EWK H_{T}^{soft}> (GeV)");

	TProfile *hprof_htsoft_pu_rms_bdt  = new TProfile("hprof_htsoft_pu_rms_bdt","",50,0.,50,0.,300.,"s");
	hprof_htsoft_pu_rms_bdt->GetXaxis()->SetTitle("number of PVs");
	hprof_htsoft_pu_rms_bdt->GetYaxis()->SetTitle("<EWK H_{T}^{soft}>, BDT > 0.92 (GeV)");


        
        
        
        
        TH2F * hMll_deltaM = new TH2F ("hMll_deltaM", "hMll_deltaM",12,115.,135., 5,0,0.1);
        hMll_deltaM->GetXaxis()->SetTitle("m(#mu#mu) (GeV)");
        hMll_deltaM->GetYaxis()->SetTitle("#DeltaM (GeV)");
        
        TH2F * histo_TrueIndices_MVAIndices = new TH2F ("histo_TrueIndices_MVAIndices", "histo_TrueIndices_MVAIndices",2,0.,2, 2,0,2);
        histo_TrueIndices_MVAIndices->GetXaxis()->SetTitle("MVA indices are correct");
        histo_TrueIndices_MVAIndices->GetYaxis()->SetTitle("leading and subleading are correct");

        TH2F * histo_TrueIndices_MVAIndices_3Jets = new TH2F ("histo_TrueIndices_MVAIndices_3Jets", "histo_TrueIndices_MVAIndices in events with 3 jets",2,0.,2, 2,0,2);
        histo_TrueIndices_MVAIndices_3Jets->GetXaxis()->SetTitle("MVA indices are correct");
        histo_TrueIndices_MVAIndices_3Jets->GetYaxis()->SetTitle("leading and subleading are correct");
        

        std::vector<TH2F*> histo2D_vector;
//         std::vector<std::string> variablesName_in_2D_plot = {"Zll_mass", "deltaM", "deltaR1", "deltaR2", "RpT", "zStar", "BDToutput", "genJetMassLeading", "hgenJetMassMatched", "Mqq" };
        std::vector<std::string> variablesName_in_2D_plot = {/*"singleMuTrigger", "doubleMuTrigger",*/ "pTmumu", "deltaPhiQQ", "higgsSisterJetIndex1", "higgsSisterJetIndex2"};

//         std::vector<float> limitDown  = {115,   0.,     -6,     -6,     0.,     -8,     0,      0,      0};
//         std::vector<float> limitUp    = {135,   0.1,    0.,     0.,     0.5,    3.,     1.6,    1500,    1500};
//         std::vector<int> binNumber  =   {12,    5,      6,      6,      10,     5,      16,     100,     100};
        
//         std::vector<float> limitDown  = {115,   0.,     0,     0,     0.,     -8,     0,      0,      0,      0};
//         std::vector<float> limitUp    = {135,   0.1,    6,     6,     0.5,    3.,     1.6,    1500,    1500,    1500};
//         std::vector<int> binNumber  =   {12,    5,      20,      20,      10,     5,      16,     100,     100,     100};
        
        
        
        std::vector<float> limitDown  = {0, 0, -1, -1};
        std::vector<float> limitUp    = {200, 3.2, 6, 6};
        std::vector<int> binNumber  =   {20,  32, 7, 7};
        
        
        for (int i = 0; i < variablesName_in_2D_plot.size(); i++) {
            for (int j = i+1; j < variablesName_in_2D_plot.size(); j++) {

                TH2F * histo2D = new TH2F (("histo2D_"+variablesName_in_2D_plot[i]+"_"+variablesName_in_2D_plot[j]).c_str(), "", binNumber[i], limitDown[i], limitUp[i], binNumber[j], limitDown[j], limitUp[j]);
                histo2D->GetXaxis()->SetTitle(variablesName_in_2D_plot[i].c_str());
                histo2D->GetYaxis()->SetTitle(variablesName_in_2D_plot[j].c_str());
                histo2D_vector.push_back(histo2D);
                
                
            }
        }
        
        

//        const int numArray= 109;  //64+8 
//         const inDt numArray= 111;  //64+8 
//         TH1D* histArray[numArray] = { hMqq, hEtaQQ,hHTsoft,hSoft_n2,hSoft_n5,hSoft_n10,hHTsoftEWK,hSoft_n2EWK,hSoft_n5EWK,hSoft_n10EWK,hHTsoftEWK_bdt,hSoft_n2EWK_bdt,hSoft_n5EWK_bdt, hSoft_n10EWK_bdt,hnPVs, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hmet,   hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt,hV_mass, hqgl, hqgl2, hZll_mass, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRptHard, hEtaQQSum, hPhiZQ1, hZll_y, hZll_ystar, hZll_zstar, hMqq_log, hlheV_pt, hJet3_pt, hlheHT_log, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh,hbdt_atanh2 , hlepton1_iso03, hlepton2_iso03, hveto_jet3pt_nom, hveto_jet3pt_denom, hveto_ht_nom, hveto_ht_denom, hveto_softht_nom, hveto_softht_denom, hveto_softpt_nom, hveto_softpt_denom, hJet2q_phi, hJet1q_pffhi, hNAdJets, hNAdJets_bdt, hJet3_pt_bdt, hAdJetHT_bdt, hNAdJets_bdt2, hJet3_pt_bdt2, hAdJetHT_bdt2,hNAdJets_mjj1, hJet3_pt_mjj1, hAdJetHT_mjj1,hNAdJets_mjj2, hJet3_pt_mjj2, hAdJetHT_mjj2, hHTsoftEWK_bdt2,hSoft_n2EWK_bdt2,hSoft_n5EWK_bdt2, hSoft_n10EWK_bdt2,hHTsoftEWK_mjj1, hSoft_n2EWK_mjj1,hSoft_n5EWK_mjj1,hSoft_n10EWK_mjj1, hHTsoftEWK_mjj2,hSoft_n2EWK_mjj2,hSoft_n5EWK_mjj2,hSoft_n10EWK_mjj2 ,hJet1q_eta_bdt, hJet1q_eta_bdt2, hJet2q_eta_bdt, hJet2q_eta_bdt2, hsoftleadTrackPt, hsoftleadTrackEta, hAdJetHT, hJet3_eta , hJet3_pt_new , hJet3_eta_bdt, hJet3_eta_bdt2, hThetaStar, hMaxJetBTag};


        TH1D* histo_without_qglCorrection = new TH1D ("histo_without_qglCorrection", "", 1, 0, 2); 
        
        const int numArray= 252;//64+8 
        TH1D* histArray[numArray] = { hMqq,hMqq_cut06, hMqq_cut08, hMqq_cut10, hMqq_cut12, hgenJetMass_cut06, hgenJetMass_cut08, hgenJetMass_cut10, hgenJetMass_cut12, hEtaQQ,hSoft_n2,hSoft_n5,hSoft_n10,hHTsoftEWK,hSoft_n2EWK,hSoft_n5EWK,hSoft_n10EWK,hHTsoftEWK_bdt,hSoft_n2EWK_bdt,hSoft_n5EWK_bdt, hSoft_n10EWK_bdt,hnPVs, hNjet25, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hVtype, hVtypeSim, hmet, hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt, hqgl, hqgl_noQGLweight, hqgl_noQGLweight2, hqgl2,hqglAtanh, hqgl2Atanh, hHll_mass_precise, hHll_mass_unprecise, hZll_mass, hZll_mass_biasUp, hZll_mass_biasDown, hZll_mass_resolution, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRpt, hRptAtanh, hEtaQQSum, hPhiZQ1, hPhiZQ2, hEtaHQ1, hEtaHQ2, hZll_y, hZll_ystar, hZll_zstar, hzstar, hZll_zstar_log, hMqq_log, hlheV_pt, hlheNpNLO, hJet3_pt, hJet3_pt_log, hlheHT_log, hlheHT_log_BDTgt1, hlheHT, hlheHT_BDTgt1, hlheNj, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh,hbdt_atanh2 , hlepton1_iso03, hlepton2_iso03, hveto_jet3pt_nom, hveto_jet3pt_denom, hveto_ht_nom, hveto_ht_denom, hveto_softht_nom, hveto_softht_denom, hveto_softpt_nom, hveto_softpt_denom, hJet2q_phi, hJet1q_phi, hJet_genJetIdx, hNAdJets, hNAdJets_bdt, hJet3_pt_bdt, hAdJetHT_bdt, hNAdJets_bdt2, hJet3_pt_bdt2, hAdJetHT_bdt2, hHTsoftEWK_bdt2,hSoft_n2EWK_bdt2,hSoft_n5EWK_bdt2, hPtSoftJets, hJet1q_eta_bdt, hJet1q_eta_bdt2, hJet2q_eta_bdt, hJet2q_eta_bdt2, hsoftleadTrackPt, hsoftleadTrackEta, hAdJetHT, hJet3_eta, hJet3_etaRatio, hJet3_pt_new , hJet3_eta_bdt, hJet3_eta_bdt2, hBDT_VBF, hBDT_VBF_atanh, hBDT_VBF_atanh_m125ForAll, hBDT_VBF_atanh_m125ControlRegion, hBDT_VBF_atanh_m125ControlRegionDown, hBDT_VBF_atanh_m125ControlRegionUp, hBDT_VBF_atanh_m125_Rpt0, hBDT_VBF_atanh_m125_WM1Log5, hBDT_VBF_atanh_m125_zStar0, hBDT_VBF_atanh_m125_ptll300, hBDT_VBF_atanh_m125_softN50, hBDT_VBF_atanh_m125_MqqLog75, hBDT_VBF_atanh_m125_mumujjPt,  hBDT_VBF_atanh_m125_DEtajj, hBDT_VBF_atanh_m125_Theta2, hBDT_VBF_atanh_m125_onlyMqq, hBDT_VBF_atanh_m125_onlyRpt, hBDT_VBF_atanh_m125_onlyWM1, hBDT_VBF_atanh_m125_onlyzStar, hBDT_VBF_atanh_m125_onlyptll, hBDT_VBF_atanh_m125_onlysoftN5, hBDT_VBF_atanh_m125_onlymumujjPt, hBDT_VBF_atanh_m125_onlyDEtajj, hBDT_VBF_atanh_m125_onlyTheta2, hBDT_VBF_atanh_biasUp, hBDT_VBF_atanh_biasDown, hBDT_VBF_atanh_resolutionDown, hBDT_VBF_atanh_resolutionUp, hBDT_VBF_atanh_QGLUp, hBDT_VBF_atanh_QGLDown, hBDT_VBF_atanh_genMqqGt350, hBDT_VBF_atanh_genMqqLs350, hBDT_VBF_atanh_genMqq0, hBDT_VBF_atanh_negativeFraction, hBDT_VBF_atanh_findBinning, hBDT_VBF_atanh_FITTING, hBDT_VBF_atanh_cumulativeLHENpNLO_0, hBDT_VBF_atanh_cumulativeLHENpNLO_1, hBDT_VBF_atanh_cumulativeLHENpNLO_2, hThetaStarJet, hThetaPlanes, hThetaStarJetAtanh, hThetaPlanesAtanh, hThetaStar, hThetaStarAbs, hMaxJetBTagCSV,hBDiscriminator_CSV, hBDiscriminator_CMVA, hweights_weighted,hweights,hdeltaMRel,hdeltaM, hnormalizedDistance_from_mH, hmaxAbsEta_BDTmFixgt08, hmaxAbsEta_BDTmFixgt1,hmaxAbsEta, hminAbsEta, hminAbsEta_GEN, hminAbsEta_cut06, hminAbsEta_cut08, hminAbsEta_cut10, hminAbsEta_cut12, hIsolatedElectrons, hHiggsSister1_Leading_eta,hHiggsSister2_Subleading_eta,hHiggsSister1_Leading_phi,hHiggsSister2_Subleading_phi,hHiggsSister1_Leading_R,hHiggsSister2_Subleading_R, hHiggsSister1_HiggsSister2_R, hdeltaR1, hdeltaR2 ,  hMaxJetBTagCMVA,hTotalEnergy,hTotalEnergylog,hWWmass,hDiffmass,hpdgId,hgen_mass,hnGenLep, hnGenJet, hnGenJetsWithoutLeptons, hpuId_jetWithoutGenJet, hpuId_AllJets, hGenLepton_matching1, hGenLepton_matching2, hrecoJetPt_genJetPt1, hrecoJetPt_genJetPt2, hgenJetMass, hgenJetMass_matched, hEnergy_fraction_Parton2_log,hEnergy_fraction_Parton1_log, hmumujj_pt, hmumujj_ptLog, hEnergy_fraction_Parton1,hPz,hPzAbs, hPzAbsLog, hInvariant_MassLog, hInvariant_Mass, hthetastar_W2toHW1, hthetastar_W1toHW2,hthetastar_HtoWW, hEnergy_fraction_Parton2, hVirtual_Wmass1,hVirtual_Wmass2,hVirtual_Wmass1_log,hVirtual_Wmass2_log,hVirtual_Pt1,hVirtual_Pt2, hVirtual_Pt1_log, hVirtual_Pt2_log, hVirtual_eta1,hVirtual_eta2,hTheta_HiggsJ1,hTheta_HiggsJ2,hthetastar_W1,hthetastar_W2, hVirtual_phi1, hVirtual_phi2, hParton_M1,hParton_M2,hMaxSecondJetBTagCSV, hMaxSecondJetBTagCMVA, hSelectionCuts};
     
        
    
        
        
        
//
//         std::cout << "Length of array = " << (sizeof(histArray)/sizeof(*histArray)) << std::end;
        
        for (int i=0;i<numArray;i++){
            histArray[i]->Sumw2();
        }

        hprof_htsoft_pu->Sumw2();
        hprof_htsoft_pu_bdt->Sumw2();
        hprof_htsoft_pu_rms->Sumw2();
        hprof_htsoft_pu_rms_bdt->Sumw2();


		
		TString cut_flow_names[30] = {"triggers","2jets events","q1_pt>50","q2_pt>30","Mqq>200","leptons_pt<20","(mll-mz)<15"};
		Float_t cut_flow[30] = {0,0,0,0,0,0,0};

	float qq_matching = 0;
	float qq_matching_all = 0;
	

	int nentries = tree_initial->GetEntries() ;
	

//	TF1 *func_lheHT = new TF1("func_lheHT","([0]+[1]*x+[2]*x*x+[3]*x*x*x)*TMath::Exp(-1*[4]*x)",60,4000);
//	func_lheHT->FixParameter(0,  -1.71063e+00);
//	func_lheHT->FixParameter(1, 6.90159e-02 );
///	func_lheHT->FixParameter(2, -2.83168e-04);
//	func_lheHT->FixParameter(3, 4.69007e-07);
//	func_lheHT->FixParameter(4, 7.07950e-03 );

	TF1* func_lheHT = new TF1("func_lheHT","pol6",4.2,7.8);
	func_lheHT->FixParameter(0,89.0139);
	func_lheHT->FixParameter(1,-275.535);
	func_lheHT->FixParameter(2,195.308);
	func_lheHT->FixParameter(3,-61.2467);
	func_lheHT->FixParameter(4,9.8217);
	func_lheHT->FixParameter(5,-0.791744);
	func_lheHT->FixParameter(6,0.0255211);
//	TF1* func_EtaQQ = new TF1("func_EtaQQ","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,16.);
//	func_EtaQQ->FixParameter(0,1.33558e+00);
//	func_EtaQQ->FixParameter(1,-2.60070e-01 );
//	func_EtaQQ->FixParameter(2,5.87759e-02);
//	func_EtaQQ->FixParameter(3,-4.64109e-03 );

	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);
	TF1* func_EtaQQ = new TF1("func_EtaQQ","pol7",0,10);
	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);

	func_Mqq->FixParameter(0,-1378.361758);
	func_Mqq->FixParameter(1,1327.383178);
	func_Mqq->FixParameter(2,-529.261759);
	func_Mqq->FixParameter(3,111.854413);
	func_Mqq->FixParameter(4,-13.208941);
	func_Mqq->FixParameter(5,0.826140);
	func_Mqq->FixParameter(6,-0.021376);


	TF1* func_qgl_q = new TF1("func_qgl_q","pol3",0.,1.);
	func_qgl_q->FixParameter(0,0.981581);
	func_qgl_q->FixParameter(1,-0.255505);
	func_qgl_q->FixParameter(2,0.929524);
	func_qgl_q->FixParameter(3,-0.666978);
	TF1* func_qgl_g = new TF1("func_qgl_g","pol7",0.,1.);
	func_qgl_g->FixParameter(0,0.612992);
	func_qgl_g->FixParameter(1,6.27);
	func_qgl_g->FixParameter(2,-34.3663);
	func_qgl_g->FixParameter(3,92.8668);
	func_qgl_g->FixParameter(4,-99.927);
	func_qgl_g->FixParameter(5,-21.1421);
	func_qgl_g->FixParameter(6, 113.218);
	func_qgl_g->FixParameter(7,-55.7067);
//pythia8 quark (|eta|<2.0, pT inclusive, pythia ): -0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581
//pythia8 gluon (|eta|<2.0, pT inclusive, pythia ): -55.7067*x^7 + 113.218*x^6 -21.1421*x^5 -99.927*x^4 + 92.8668*x^3 -34.3663*x^2 + 6.27*x + 0.612992

	TF1* interference_func = new TF1("interference_func","pol7",0,10);
	interference_func->FixParameter(0,-3236.73);
	interference_func->FixParameter(1,3158.71);
	interference_func->FixParameter(2,-1314.93);
	interference_func->FixParameter(3,302.849);
	interference_func->FixParameter(4,-41.6913);
	interference_func->FixParameter(5,3.4312);
	interference_func->FixParameter(6,-0.156337);
	interference_func->FixParameter(7,0.00304253);

//	rochcor2016 *rmcor = new rochcor2016();




    bool plotOutput = true;
    bool plotOutput_NN = false;
    bool histo_Preparation_To_Fit = false;

//     plotOutput = false;
//     plotOutput_NN = true;
     //histo_Preparation_To_Fit = true;

//     if (plotOutput) plotOutput_NN = false;
    if ( ((QCDScaleWeight_str.CompareTo("none")!=0)&&(QCDScaleWeight_str.CompareTo("nom")!=0)) || ((JESWeight_str.CompareTo("none")!=0)&&(JESWeight_str.CompareTo("nom")!=0)) ) histo_Preparation_To_Fit = false;



	const int Nsyst = 6;
//	TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","LHE_weights_scale_muF","LHE_weights_scale_muR","JES","JER","MDG_NLO_corr","int_shape","QGL"};
TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","JES","JER","QGL"};
	TH1D *hbdtUp[20];
	TH1D *hbdt_atanhUp[20];
	TH1D *hbdt_atanhUp_pdf[20];
	TH1D *hbdtDown[20];
	TH1D *hbdt_atanhDown[20];
	TH1D *hbdt_atanhDown_pdf[20];
    TString file_tag_hBDT_name = file_tag;
    file_tag_hBDT_name = Correct_fileTag (file_tag_hBDT_name);
	for (int i=0;i<Nsyst;i++){
		std::string histTitleUp;
		std::string histTitleDown;

		histTitleUp = ((std::string) hBDT_VBF->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] + "_Up";
		hbdtUp[i] = (TH1D*) hBDT_VBF->Clone(histTitleUp.c_str());
		histTitleUp = ((std::string) hBDT_VBF_atanh->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] +"_Up";
		hbdt_atanhUp[i] = (TH1D*) hBDT_VBF_atanh->Clone(histTitleUp.c_str());


//        hbdtUp[i]->Sumw2();
//        hbdt_atanhUp[i]->Sumw2();

		if (i!=0) {
			histTitleDown = ((std::string) hBDT_VBF->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] + "_Down";
			hbdtDown[i] = (TH1D*) hBDT_VBF->Clone(histTitleDown.c_str());
			histTitleDown = ((std::string) hBDT_VBF_atanh->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] +"_Down";
			hbdt_atanhDown[i] = (TH1D*) hBDT_VBF_atanh->Clone(histTitleDown.c_str());
//            hbdtDown[i]->Sumw2();
//            hbdt_atanhDown[i]->Sumw2();
		}
	}



float scaleWeightsUp[20];
float scaleWeightsDown[20];
string file_tag_str = file_tag.Data();
TString uncNotAppl[20] = {"","data","WW_WZ_ZZ_ST_tW","WW_WZ_ZZ_ST_tW","WW_WZ_ZZ_ST_tW","","","",""};

int Nsyst_NoConst = Nsyst;
if (data==1) Nsyst_NoConst = 1;

//for (int current_syst=0;current_syst<Nsyst;current_syst++){
////for (int current_syst=0;current_syst<1;current_syst++){
////	if ((current_syst!=0) && (current_syst!=9)) continue;

//	if ((data==1)&&(current_syst!=0)) continue;
//	if (( file_tag.CompareTo("EWKinterference")==0 )&&(current_syst!=0)) continue;
//	string uncNotAppl_str = uncNotAppl[current_syst].Data();
//	if (current_syst!=0) if (uncNotAppl_str.find(file_tag_str)!=std::string::npos) continue;

//	bdt = 0;
//	passSel=0;
//	genweight = 0 ;
//	cout<<current_syst<<endl;
//}


        TFile fileMVA("mvaTree/main_tmva_tree_"+file_tag+"_"+postfix+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_JER"+JERWeight_str+"_PU"+PUWeight_str+".root","recreate");
//         TFile fileMVAJet("jetMVAtree/"+file_tag+"_"+postfix+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_JER"+JERWeight_str+"_PU"+PUWeight_str+".root","recreate");

                
        TMVAstruct TMVA;
        for (int n = 0; n <= number_LHE_weights_pdf; n++) TMVA.genweightVECTOR.push_back(0);
            
        TTree *treeMVA = new TTree("TMVAtree","TMVAtree");

        
        treeMVA->Branch("ll_mass",&TMVA.ll_mass,"ll_mass/F");
        treeMVA->Branch("ll_mass_biasUp",&TMVA.ll_mass_biasUp,"ll_mass_biasUp/F");
        treeMVA->Branch("ll_mass_biasDown",&TMVA.ll_mass_biasDown,"ll_mass_biasDown/F");
        treeMVA->Branch("ll_mass_resolution",&TMVA.ll_mass_resolution,"ll_mass_resolution/F");
        treeMVA->Branch("diffMassWWH",&TMVA.diffMassWWH,"diffMassWWH/F");
	treeMVA->Branch("ll_pt",&TMVA.ll_pt,"ll_pt/F");
	treeMVA->Branch("ll_eta",&TMVA.ll_eta,"ll_eta/F");
	treeMVA->Branch("ll_ystar",&TMVA.ll_ystar,"ll_ystar/F");
        treeMVA->Branch("ll_zstar",&TMVA.ll_zstar,"ll_zstar/F");   
        treeMVA->Branch("ll_pt_Log",&TMVA.ll_pt_Log,"ll_pt_Log/F");
        
        treeMVA->Branch("q1_eta",&TMVA.q1_eta,"q1_eta/F");
        treeMVA->Branch("met_pt",&TMVA.met_pt,"met/F");
        treeMVA->Branch("EWKHTsoft",&TMVA.EWKHTsoft,"EWKHTsoft/F");
        treeMVA->Branch("btagCMVA",&TMVA.btagCMVA,"btagCMVA/F");
        treeMVA->Branch("btagCSV",&TMVA.btagCSV,"btagCSV/F");
        treeMVA->Branch("btagCMVA_second",&TMVA.btagCMVA_second,"btagCMVA_second/F");
        treeMVA->Branch("btagCSV_second",&TMVA.btagCSV_second,"btagCSV_second/F");
        treeMVA->Branch("softLeadingJet_pt",&TMVA.softLeadingJet_pt,"softLeadingJet_pt/F");
        treeMVA->Branch("cosThetaStar",&TMVA.cosThetaStar,"cosThetaStar/F");
        treeMVA->Branch("cosThetaStarAbs",&TMVA.cosThetaStarAbs,"cosThetaStarAbs/F");        
        treeMVA->Branch("absCosThetaStarJet",&TMVA.absCosThetaStarJet,"absCosThetaStarJet/F");        
        treeMVA->Branch("cosThetaPlane",&TMVA.cosThetaPlane,"cosThetaPlane/F");        
        //treeMVA->Branch("PhiZQ1", &PhiZQ1,"PhiZQ1/F");
        
	treeMVA->Branch("Mqq",&TMVA.Mqq,"Mqq/F");
        treeMVA->Branch("MqqLog",&TMVA.MqqLog,"MqqLog/F");
	treeMVA->Branch("evt",&evt,"evt/l");
	treeMVA->Branch("run",&runToWrite,"run/I");
	treeMVA->Branch("luminosityBlock",&lumiToWrite,"luminosityBlock/I");
	treeMVA->Branch("DeltaEtaQQ",&TMVA.DeltaEtaQQ,"DeltaEtaQQ/F");
	treeMVA->Branch("Jet2q_pt",&TMVA.Jet2q_pt,"Jet2q_pt/F");
	treeMVA->Branch("Jet1q_pt",&TMVA.Jet1q_pt,"Jet1q_pt/F");
        treeMVA->Branch("Jet2q_eta",&TMVA.Jet2q_eta,"Jet2q_eta/F");
        treeMVA->Branch("Jet1q_eta",&TMVA.Jet1q_eta,"Jet1q_eta/F");
        
        treeMVA->Branch("Jet2q_chEmEF",&TMVA.Jet2q_chEmEF,"Jet2q_chEmEF/F");
        treeMVA->Branch("Jet1q_chEmEF",&TMVA.Jet1q_chEmEF,"Jet1q_chEmEF/F");
        treeMVA->Branch("Jet2q_neEmEF",&TMVA.Jet2q_neEmEF,"Jet2q_neEmEF/F");
        treeMVA->Branch("Jet1q_neEmEF",&TMVA.Jet1q_neEmEF,"Jet1q_neEmEF/F");
        
        treeMVA->Branch("Jet2q_chHEF",&TMVA.Jet2q_chHEF,"Jet2q_chHEF/F");
        treeMVA->Branch("Jet1q_chHEF",&TMVA.Jet1q_chHEF,"Jet1q_chHEF/F");
        treeMVA->Branch("Jet2q_neHEF",&TMVA.Jet2q_neHEF,"Jet2q_neHEF/F");
        treeMVA->Branch("Jet1q_neHEF",&TMVA.Jet1q_neHEF,"Jet1q_neHEF/F");
        
        
        
        treeMVA->Branch("Jet1q_puId",&TMVA.Jet1q_puId,"Jet1q_puId/F");
        treeMVA->Branch("Jet2q_puId",&TMVA.Jet2q_puId,"Jet2q_puId/F");
     
        treeMVA->Branch("Jet1q_pt_jerUp",&TMVA.Jet1q_pt_jerUp,"Jet1q_pt_jerUp/F");
        treeMVA->Branch("Jet1q_pt_jerDown",&TMVA.Jet1q_pt_jerDown,"Jet1q_pt_jerDown/F");
        treeMVA->Branch("Jet2q_pt_jerUp",&TMVA.Jet2q_pt_jerUp,"Jet2q_pt_jerUp/F");
        treeMVA->Branch("Jet2q_pt_jerDown",&TMVA.Jet2q_pt_jerDown,"Jet2q_pt_jerDown/F");
        
        
        treeMVA->Branch("genJetPt_match1",&TMVA.genJetPt_match1,"genJetPt_match1/F");
        treeMVA->Branch("genJetPt_match2",&TMVA.genJetPt_match2,"genJetPt_match2/F");

        treeMVA->Branch("genJetIdx0",&TMVA.genJetIdx0,"genJetIdx0/I");
        treeMVA->Branch("genJetIdx1",&TMVA.genJetIdx1,"genJetIdx1/I");
        treeMVA->Branch("genJetIdx2",&TMVA.genJetIdx2,"genJetIdx2/I");
        treeMVA->Branch("genJetIdx",&TMVA.genJetIdx,"genJetIdx[30]/I");
        
        treeMVA->Branch("genJetIdx0_genJetMassWithoutLepton",&TMVA.genJetIdx0_genJetMassWithoutLepton,"genJetIdx0_genJetMassWithoutLepton/I");
        treeMVA->Branch("genJetIdx1_genJetMassWithoutLepton",&TMVA.genJetIdx1_genJetMassWithoutLepton,"genJetIdx1_genJetMassWithoutLepton/I");
        treeMVA->Branch("genJetIdx2_genJetMassWithoutLepton",&TMVA.genJetIdx2_genJetMassWithoutLepton,"genJetIdx2_genJetMassWithoutLepton/I");
        
        treeMVA->Branch("Jet1q_ptLog",&TMVA.Jet1q_ptLog,"Jet1q_ptLog/F");
        treeMVA->Branch("Jet2q_ptLog",&TMVA.Jet2q_ptLog,"Jet2q_ptLog/F");
        treeMVA->Branch("Jet3q_ptLog",&TMVA.Jet3q_ptLog,"Jet3q_ptLog/F");
        
        treeMVA->Branch("Jet_pt",&TMVA.Jet_pt,"Jet_pt[30]/F");

                
//	treeMVA->Branch("Jet1q_leadTrackPt",&TMVA.Jet1q_leadTrackPt,"Jet1q_leadTrackPt/F");
//	treeMVA->Branch("Jet2q_leadTrackPt",&TMVA.Jet2q_leadTrackPt,"Jet2q_leadTrackPt/F");
//	treeMVA->Branch("axis2_jet1",&TMVA.axis2_jet1,"axis2_jet1/F");
//	treeMVA->Branch("axis2_jet2",&TMVA.axis2_jet2,"axis2_jet2/F");
	treeMVA->Branch("qq_pt",&TMVA.qq_pt,"qq_pt/F");
	treeMVA->Branch("qgl_1q",&TMVA.qgl_1q,"qgl_1q/F");
	treeMVA->Branch("qgl_2q",&TMVA.qgl_2q,"qgl_2q/F");

	treeMVA->Branch("Jet3_pt",&TMVA.Jet3_pt,"Jet3_pt/F");
	treeMVA->Branch("Rpt",&TMVA.Rpt,"Rpt/F");
    treeMVA->Branch("Jet3_eta",&TMVA.Jet3_eta,"Jet_eta/F");
    treeMVA->Branch("mumujj_pt",&TMVA.mumujj_pt,"mumujj_pt/F");
        
    
    
    



        
    
        
        treeMVA->Branch("softActivityEWK_njets2",&TMVA.softActivityEWK_njets2,"softActivityEWK_njets2/I");
        treeMVA->Branch("softActivityEWK_njets5",&TMVA.softActivityEWK_njets5,"softActivityEWK_njets5/I");
        treeMVA->Branch("softActivityEWK_njets10",&TMVA.softActivityEWK_njets10,"softActivityEWK_njets10/I");
        
        
        treeMVA->Branch("randomVariable",&TMVA.randomVariable,"randomVariable/F");
        treeMVA->Branch("xSection",&TMVA.xSection,"CrossSectionOfTheSample/F");

 

        treeMVA->Branch("Inv_mass",&TMVA.Inv_mass,"Inv_mass/F");
        treeMVA->Branch("Invariant_MassLog",&TMVA.Invariant_MassLog,"Invariant_MassLog/F");
        treeMVA->Branch("X_parton1",&TMVA.X_parton1,"X_parton1/F");
	treeMVA->Branch("X_parton2",&TMVA.X_parton2,"X_parton2/F");
	treeMVA->Branch("W_mass_virtual1",&TMVA.W_mass_virtual1,"W_mass_virtual1/F");
	treeMVA->Branch("W_mass_virtual2",&TMVA.W_mass_virtual2,"W_mass_virtual2/F");
	treeMVA->Branch("W_Pt_virtual1",&TMVA.W_Pt_virtual1,"W_Pt_virtual1/F");
 
        
        treeMVA->Branch("W_Pt_virtual2",&TMVA.W_Pt_virtual2,"W_Pt_virtual2/F");
        treeMVA->Branch("W_eta_virtual1",&TMVA.W_eta_virtual1,"W_eta_virtual1/F");
	treeMVA->Branch("W_eta_virtual2",&TMVA.W_eta_virtual2,"W_eta_virtual2/F");
	treeMVA->Branch("W_phi_virtual1",&TMVA.W_phi_virtual1,"W_phi_virtual1/F");  
	treeMVA->Branch("lepton1_pt",&TMVA.lepton1_pt,"lepton1_pt/F");
	treeMVA->Branch("lepton2_pt",&TMVA.lepton2_pt,"lepton2_pt/F");  
	treeMVA->Branch("jets12",&TMVA.jets12,"jets12/F");
	treeMVA->Branch("W_phi_virtual2",&TMVA.W_phi_virtual2,"W_phi_virtual2/F");
        treeMVA->Branch("thetastarW1",&TMVA.thetastarW1,"thetastarW1/F");
	treeMVA->Branch("thetastarW2",&TMVA.thetastarW2,"thetastarW2/F");
	treeMVA->Branch("thetastarW2toHW1",&TMVA.thetastarW2toHW1,"thetastarW2toHW1/F");
	treeMVA->Branch("thetastarW1toHW2",&TMVA.thetastarW1toHW2,"thetastarW1toHW2/F");
	treeMVA->Branch("thetastarHtoWW",&TMVA.thetastarHtoWW,"thetastarHtoWW/F");  
	
	treeMVA->Branch("theta1",&TMVA.theta1,"theta1/F");
        treeMVA->Branch("theta2",&TMVA.theta2,"theta2/F");

	treeMVA->Branch("WWmass",&TMVA.WWmass,"WWmass/F");
	treeMVA->Branch("impulsoZ",&TMVA.impulsoZ,"impulsoZ/F");  
    treeMVA->Branch("energytot",&TMVA.energytot,"energytot/F");
	treeMVA->Branch("energytotLog",&TMVA.energytotLog,"energytotLog/F");  
    treeMVA->Branch("deltaMRel", &TMVA.deltaMRel,"deltaMRel/F");
    treeMVA->Branch("deltaM",  &TMVA.deltaM,"deltaM/F"); 
    treeMVA->Branch("normalizedDistance_from_mH",  &TMVA.normalizedDistance_from_mH,"normalizedDistance_from_mH/F"); 
    
    
                       
        treeMVA->Branch("qgl_1qAtanh",&TMVA.qgl_1qAtanh,"qgl_1qAtanh/F");
        treeMVA->Branch("qgl_2qAtanh",&TMVA.qgl_2qAtanh,"qgl_2qAtanh/F");
        
        treeMVA->Branch("cosThetaPlaneAtanh",&TMVA.cosThetaPlaneAtanh,"cosThetaPlaneAtanh/F");
        treeMVA->Branch("absCosThetaStarJetAtanh",&TMVA.absCosThetaStarJetAtanh,"absCosThetaStarJetAtanh/F");
        treeMVA->Branch("X_parton1Log",&TMVA.X_parton1Log,"X_parton1Log/F");
        treeMVA->Branch("X_parton2Log",&TMVA.X_parton2Log,"X_parton2Log/F");
        treeMVA->Branch("W_mass_virtual1Log",&TMVA.W_mass_virtual1Log,"W_mass_virtual1Log/F");
        treeMVA->Branch("W_mass_virtual2Log",&TMVA.W_mass_virtual2Log,"W_mass_virtual2Log/F");
        treeMVA->Branch("W_Pt_virtual1Log",&TMVA.W_Pt_virtual1Log,"W_Pt_virtual1Log/F");
        treeMVA->Branch("W_Pt_virtual2Log",&TMVA.W_Pt_virtual2Log,"W_Pt_virtual2Log/F");
        
        treeMVA->Branch("maxAbsEta",&TMVA.maxAbsEta,"maxAbsEta/F");
        treeMVA->Branch("minAbsEta",&TMVA.minAbsEta,"minAbsEta/F");
        treeMVA->Branch("qqMass_skim",&TMVA.qqMass_skim,"qqMass_skim/F");
        treeMVA->Branch("MuonMass_skim",&TMVA.MuonMass_skim,"MuonMass_skim/F");  
        treeMVA->Branch("genJetMass_leading",&TMVA.genJetMass_leading,"genJetMass_leading/F");  
//         treeMVA->Branch("genJetMass_leading",&TMVA.genJetMass_leading,"genJetMass_leading/F");  

        treeMVA->Branch("JetPuId_maxAbsEta",&TMVA.JetPuId_maxAbsEta,"JetPuId_maxAbsEta/F");
        treeMVA->Branch("Muon1_relIso04",&TMVA.Muon1_relIso04,"Muon1_relIso04/F");
        treeMVA->Branch("Muon2_relIso04",&TMVA.Muon2_relIso04,"Muon2_relIso04/F");
         
         treeMVA->Branch("countJet25",&TMVA.countJet25,"countJet25/I");
        
        treeMVA->Branch("nGenJet",&nGenJet,"nGenJet/I");  
        treeMVA->Branch("GenJet_pt",GenJet_pt,"GenJet_pt[30]/F");  
        treeMVA->Branch("GenJet_eta",GenJet_eta,"GenJet_eta[30]/F");  
        treeMVA->Branch("GenJet_phi",GenJet_phi,"GenJet_phi[30]/F");  
        treeMVA->Branch("GenJet_mass",GenJet_mass,"GenJet_mass[30]/F");  

        
        treeMVA->Branch("quark1_pt",&TMVA.quark1_pt,"quark1_pt/F");  
        treeMVA->Branch("quark1_eta",&TMVA.quark1_eta,"quark1_eta/F");  
        treeMVA->Branch("quark1_phi",&TMVA.quark1_phi,"quark1_phi/F");  
        treeMVA->Branch("quark1_mass",&TMVA.quark1_mass,"quark1_mass/F");        
        treeMVA->Branch("quark2_pt",&TMVA.quark2_pt,"quark2_pt/F");  
        treeMVA->Branch("quark2_eta",&TMVA.quark2_eta,"quark2_eta/F");  
        treeMVA->Branch("quark2_phi",&TMVA.quark2_phi,"quark2_phi/F");  
        treeMVA->Branch("quark2_mass",&TMVA.quark2_mass,"quark2_mass/F");
        
        treeMVA->Branch("genMuon1_pt",&TMVA.genMuon1_pt,"genMuon1_pt/F");  
        treeMVA->Branch("genMuon1_eta",&TMVA.genMuon1_eta,"genMuon1_eta/F");  
        treeMVA->Branch("genMuon1_phi",&TMVA.genMuon1_phi,"genMuon1_phi/F");  
        treeMVA->Branch("genMuon1_mass",&TMVA.genMuon1_mass,"genMuon1_mass/F");        
        treeMVA->Branch("genMuon2_pt",&TMVA.genMuon2_pt,"genMuon2_pt/F");  
        treeMVA->Branch("genMuon2_eta",&TMVA.genMuon2_eta,"genMuon2_eta/F");  
        treeMVA->Branch("genMuon2_phi",&TMVA.genMuon2_phi,"genMuon2_phi/F");  
        treeMVA->Branch("genMuon2_mass",&TMVA.genMuon2_mass,"genMuon2_mass/F");
        
        

        
        
        treeMVA->Branch("ngenJetsWithoutLeptonsP4",&ngenJetsWithoutLeptonsP4,"ngenJetsWithoutLeptonsP4/I");  
        treeMVA->Branch("genJetsWithoutLeptonsP4_pt",&genJetsWithoutLeptonsP4_pt);  
        treeMVA->Branch("genJetsWithoutLeptonsP4_eta",&genJetsWithoutLeptonsP4_eta);  
        treeMVA->Branch("genJetsWithoutLeptonsP4_phi",&genJetsWithoutLeptonsP4_phi);  
        treeMVA->Branch("genJetsWithoutLeptonsP4_mass",&genJetsWithoutLeptonsP4_mass);  
        
//         treeMVA->Branch("btagCMVA_leading",&TMVA.btagCMVA_leading,"btagCMVA_leading/F");

        treeMVA->Branch("BDToutput",&TMVA.BDToutput,"BDToutput/F");
        
        treeMVA->Branch("weightMVA",  &TMVA.weightMVA,"weightMVA/F"); 
        treeMVA->Branch("genweight",  &TMVA.genweight,"genweight/F"); 
        treeMVA->Branch("genweight_noQGLcorrection",  &TMVA.genweight_noQGLcorrection,"genweight_noQGLcorrection/F"); 
        
        
        treeMVA->Branch("LHE_weights_scale_wUp",&TMVA.LHE_weights_scale_wUp,"LHE_weights_scale_wUp/F");
        treeMVA->Branch("LHE_weights_scale_wDown",&TMVA.LHE_weights_scale_wDown,"LHE_weights_scale_wDown/F");
        
        
        for (int n = 0; n < number_LHE_weights_pdf; n++)
        treeMVA->Branch(("genweightVECTOR_"+std::to_string(n)).c_str(),  &TMVA.genweightVECTOR[n],"genweightVECTOR_n/F"); 

        //if  ((file_tag.CompareTo("TT")==0)) MVAcountMAX = 3300;  //3900;
        if  ((file_tag.CompareTo("DYJetstoLL_madgraph")==0)) MVAcountMAX = 100000;
        if  ((file_tag.CompareTo("DY0JetsToLL_M")==0)) 675;//   323 events        // MVAcountMAX = 1500;    //   713 events
        if  ((file_tag.CompareTo("DY1JetsToLL_M")==0)) 2313; // MVAcountMAX = 5700;
        if  ((file_tag.CompareTo("DY2JetsToLL_M")==0)) 2565; // MVAcountMAX = 4600;
        if  ((file_tag.CompareTo("DY3JetsToLL_M")==0)) 1591; // MVAcountMAX = 2550;
        if  ((file_tag.CompareTo("DY4JetsToLL_M")==0)) 1197; // MVAcountMAX = 1970;
//         if  ((file_tag.CompareTo("VBF_HToMuMu")==0)) MVAcountMAX = 200000;  //30000;//

        if  ((file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0)) MVAcountMAX = 1000000;  
        if  ((file_tag.CompareTo("DYJetsToLL_M-105To160-madgraphMLM")==0)) MVAcountMAX = 1000000;  


//         if  ((file_tag.CompareTo("TT")==0)) MVAcountMAX = 8568;  
        if  ((file_tag.CompareTo("EKW_LLJJ_INT")==0)) MVAcountMAX = 46;  
//         if  ((file_tag.CompareTo("TTTo2L2Nu")==0)) MVAcountMAX = 11889.6;  
//         if  ((file_tag.CompareTo("DYJetstoLL_madgraph")==0)) MVAcountMAX = 100000;
//         if  ((file_tag.CompareTo("DY0JetsToLL_M")==0)) MVAcountMAX = 1500;    //   713 events
//         if  ((file_tag.CompareTo("DY1JetsToLL_M")==0)) MVAcountMAX = 5700;
//         if  ((file_tag.CompareTo("DY2JetsToLL_M")==0)) MVAcountMAX = 4600;
//         if  ((file_tag.CompareTo("DY3JetsToLL_M")==0)) MVAcountMAX = 2550;
//         if  ((file_tag.CompareTo("DY4JetsToLL_M")==0)) MVAcountMAX = 1970;
//         if  ((file_tag.CompareTo("VBF_HToMuMu")==0)) MVAcountMAX = 30000;  



    TMVA::Reader *readerJetMVA = new TMVA::Reader("Silent");
    TMVA::Reader *reader = new TMVA::Reader("Silent");
    float temp_softActivityEWK_njets5 = 0.;

//    bool plotOutput = false;
//    bool plotOutput_NN = false;
//    bool histo_Preparation_To_Fit = false;


//     plotOutput = true;
////     plotOutput_NN = true;
////     histo_Preparation_To_Fit = true;
//    if (plotOutput) plotOutput_NN = false;
//    if ( ((QCDScaleWeight_str.CompareTo("none")!=0)&&(QCDScaleWeight_str.CompareTo("nom")!=0)) || ((JESWeight_str.CompareTo("none")!=0)&&(JESWeight_str.CompareTo("nom")!=0)) ) histo_Preparation_To_Fit = false;

    
    
    
    
        TFile fileMVAJet("jetMVAtree/"+file_tag+"_"+postfix+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_JER"+JERWeight_str+"_PU"+PUWeight_str+".root","recreate");
        fileMVAJet.cd();


        TTree *treeMVAjet = new TTree("TMVAtree","TMVAtree");

        int JET_MVA_fromVBF = 0;
        treeMVAjet->Branch("fromVBF",&JET_MVA_fromVBF,"fromVBF/I");
        int JET_MVA_index = 0;
        treeMVAjet->Branch("index",&JET_MVA_index,"index/I");
        float JET_MVA_p = 0;
        treeMVAjet->Branch("p",&JET_MVA_p,"p/F");
        float JET_MVA_pt = 0;
        treeMVAjet->Branch("pt",&JET_MVA_pt,"pt/F");
        float JET_MVA_eta = 0;
        treeMVAjet->Branch("eta",&JET_MVA_eta,"eta/F");
        float JET_MVA_phi = 0;
        treeMVAjet->Branch("phi",&JET_MVA_phi,"phi/F");
        float JET_MVA_mass = 0;
        treeMVAjet->Branch("mass",&JET_MVA_mass,"mass/F");
        int JET_MVA_id = 0;
        treeMVAjet->Branch("id",&JET_MVA_id,"id/I");
        int JET_MVA_puId = 0;
        treeMVAjet->Branch("puId",&JET_MVA_puId,"puId/I");
        int JET_MVA_Njets = 0;
        treeMVAjet->Branch("Njets",&JET_MVA_Njets,"Njets/I");
        
        float JET_MVA_qgl = 0;
        treeMVAjet->Branch("qgl",&JET_MVA_qgl,"qgl/F");
        float JET_MVA_area = 0;
        treeMVAjet->Branch("area",&JET_MVA_area,"area/F");
        float JET_MVA_mult = 0;
        treeMVAjet->Branch("mult",&JET_MVA_mult,"mult/F");
        
        float JET_MVA_pt1 = 0;
        treeMVAjet->Branch("pt1",&JET_MVA_pt1,"pt1/F");
        float JET_MVA_eta1 = 0;
        treeMVAjet->Branch("eta1",&JET_MVA_eta1,"eta1/F");
        float JET_MVA_phi1 = 0;
        treeMVAjet->Branch("phi1",&JET_MVA_phi1,"phi1/F");
        float JET_MVA_pt2 = 0;
        treeMVAjet->Branch("pt2",&JET_MVA_pt2,"pt2/F");
        float JET_MVA_eta2 = 0;
        treeMVAjet->Branch("eta2",&JET_MVA_eta2,"eta2/F");
        float JET_MVA_phi2 = 0;
        treeMVAjet->Branch("phi2",&JET_MVA_phi2,"phi2/F");
        float JET_MVA_pt3 = 0;
        treeMVAjet->Branch("pt3",&JET_MVA_pt2,"pt3/F");
        float JET_MVA_eta3 = 0;
        treeMVAjet->Branch("eta3",&JET_MVA_eta2,"eta3/F");
        float JET_MVA_phi3 = 0;
        treeMVAjet->Branch("phi3",&JET_MVA_phi2,"phi3/F");
        
        
        float JET_MVA_higgh_pt = 0;
        treeMVAjet->Branch("higgh_pt",&JET_MVA_higgh_pt,"higgh_pt/F");
        float JET_MVA_higgh_eta = 0;
        treeMVAjet->Branch("higgh_eta",&JET_MVA_higgh_eta,"higgh_eta/F");
        float JET_MVA_higgh_phi = 0;
        treeMVAjet->Branch("higgh_phi",&JET_MVA_higgh_phi,"higgh_phi/F");
        float JET_MVA_higgh_mass = 0;
        treeMVAjet->Branch("higgh_mass",&JET_MVA_higgh_mass,"higgh_mass/F");
        
        int JET_MVA_etaIdx = 0;
        treeMVAjet->Branch("etaIdx",&JET_MVA_etaIdx,"etaIdx/I");
        int JET_MVA_NjetsEtaMinus = 0;
        treeMVAjet->Branch("NjetsEtaMinus",&JET_MVA_NjetsEtaMinus,"NjetsEtaMinus/I");
        int JET_MVA_NjetsEtaPlus = 0;
        treeMVAjet->Branch("NjetsEtaPlus",&JET_MVA_NjetsEtaPlus,"NjetsEtaPlus/I");
        
        
        float JET_MVA_index_float = JET_MVA_index;
        float JET_MVA_puId_float = JET_MVA_puId;
        float JET_MVA_Njets_float = JET_MVA_Njets;
        float JET_MVA_etaIdx_float = JET_MVA_etaIdx;
        float JET_MVA_NjetsEtaMinus_float = JET_MVA_NjetsEtaMinus;
        float JET_MVA_NjetsEtaPlus_float = JET_MVA_NjetsEtaPlus;
        
        readerJetMVA->AddVariable("pt",&JET_MVA_pt); 
        readerJetMVA->AddVariable("eta",&JET_MVA_eta); 
        readerJetMVA->AddVariable("index",&JET_MVA_index_float); 
        readerJetMVA->AddVariable("area",&JET_MVA_area); 
        readerJetMVA->AddVariable("puId",&JET_MVA_puId_float); 
        readerJetMVA->AddVariable("Njets",&JET_MVA_Njets_float); 
        readerJetMVA->AddVariable("p",&JET_MVA_p); 
        
        
        
        readerJetMVA->AddVariable("etaIdx",&JET_MVA_etaIdx_float); 
        readerJetMVA->AddVariable("NjetsEtaMinus",&JET_MVA_NjetsEtaMinus_float); 
        readerJetMVA->AddVariable("NjetsEtaPlus",&JET_MVA_NjetsEtaPlus_float); 
        
        readerJetMVA->AddVariable("pt1",&JET_MVA_pt1); 
        readerJetMVA->AddVariable("eta1",&JET_MVA_eta1); 
        readerJetMVA->AddVariable("pt2",&JET_MVA_pt2); 
        readerJetMVA->AddVariable("eta2",&JET_MVA_eta2); 
//         readerJetMVA->AddVariable("pt3",&JET_MVA_pt3); 
//         readerJetMVA->AddVariable("eta3",&JET_MVA_eta3); 
        
    
    
//         readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_index_puId_Njets_p.weights.xml");
//         readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_index_puId_Njets_p_pt1_eta1_pt2_eta2_pt3_eta3.weights.xml");
//         readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_index_puId_Njets_p.weights_etaIdx_NjetsEtaMinus_NjetsEtaPlus.xml");
//         readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_puId_Njets_p.weights.xml");
//         readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_index.weights.xml");
readerJetMVA->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/jetBDT/firstTest/TMVAClassification_BDTG_pt_eta_index_area_puId_Njets_p_etaIdx_NjetsEtaMinus_NjetsEtaPlus_pt1_eta1_pt2_eta2_pt3_eta3.weights.xml");

    if (plotOutput) {
        
// // // // // // // // //    VARIABILI GIUGNO
// // // // // // // // //     reader->AddVariable("ll_mass",&TMVA.ll_mass);    
// // // // // // // // //     reader->AddVariable("Mqq",&TMVA.Mqq);
// // // // // // // // //     
// // // // // // // // //     reader->AddVariable("RptHard",&TMVA.RptHard);
// // // // // // // // //     reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
// // // // // // // // // //    
// // // // // // // // //     reader->AddVariable("ll_pt",&TMVA.ll_pt);
// // // // // // // // //     reader->AddVariable("ll_eta",&TMVA.ll_eta);
// // // // // // // // //     reader->AddVariable("Jet2q_pt",&TMVA.Jet2q_pt);
// // // // // // // // //      reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
    
        
        
        
        
/////////////////////////VARIABILI APRILE BTDG/////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("Rpt",&TMVA.Rpt);
// //         reader->AddVariable("W_mass_virtual1Log",&TMVA.W_mass_virtual1Log); 
//         reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1Log); 
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
//         reader->AddVariable("ll_pt",&TMVA.ll_pt);
// //         reader->AddVariable("W_mass_virtual2Log",&TMVA.W_mass_virtual2Log);
//         reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2Log);
        
        
//         {"ll_mass","MqqLog", "Rpt","W_mass_virtual1Log","ll_zstar","softActivityEWK_njets5", "ll_pt","W_mass_virtual2Log"}
        
        
        
        
        /////////////////////////VARIABILI  MAGGIO /////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("Rpt",&TMVA.Rpt);
//         reader->AddVariable("W_mass_virtual1Log",&TMVA.W_mass_virtual1Log); 
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
//         reader->AddVariable("ll_pt",&TMVA.ll_pt);
        
                
        /////////////////////////VARIABILI  FILIPPO /////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("Mqq",&TMVA.Mqq);
//         reader->AddVariable("Rpt",&TMVA.Rpt);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);
//         reader->AddVariable("qq_pt",&TMVA.qq_pt);
//         reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1);
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
        
//////////////////////////  VARIABILI TRAINING NANO 2016/////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("mumujj_pt",&TMVA.mumujj_pt);
//         reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//         reader->AddVariable("ll_pt",&TMVA.ll_pt);
//         reader->AddVariable("theta2",&TMVA.theta2); 
        
//////////////////////////  VARIABILI TRAINING NANO 2016  WITH DY_VBFFILTER/////////////////////////////////////////////////////////
        reader->AddVariable("ll_mass",&TMVA.ll_mass);   
        reader->AddVariable("MqqLog",&TMVA.MqqLog);
        reader->AddVariable("mumujj_pt",&TMVA.mumujj_pt);
        reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
        reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
        reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
        reader->AddVariable("ll_pt",&TMVA.ll_pt);
        reader->AddVariable("theta2",&TMVA.theta2); 
        
        reader->AddVariable("impulsoZ",&TMVA.impulsoZ);
        reader->AddVariable("maxAbsEta",&TMVA.maxAbsEta);
//         reader->AddVariable("qgl_2qAtanh",&TMVA.qgl_2qAtanh); 
        
   
        
//////////////////////////  VARIABILI TRAINING WITHOUT SOFTN5/////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("mumujj_pt",&TMVA.mumujj_pt);
//         reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);  
        
        
        
/////////////////////////VARIABILI APRILE MLP/////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("Rpt",&TMVA.Rpt);
//         reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1); 
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
//         reader->AddVariable("ll_pt",&TMVA.ll_pt);
//         reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2); 
//         reader->AddVariable("qgl_1qAtanh",&TMVA.qgl_1qAtanh); 
//         
        
        
        
/////////////////////////VARIABILI APRILE BTDG W/O WMASSES /////////////////////////////////////////////////////////
//         reader->AddVariable("ll_mass",&TMVA.ll_mass);   
//         reader->AddVariable("MqqLog",&TMVA.MqqLog);
//         reader->AddVariable("Rpt",&TMVA.Rpt);
//         reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//         reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
//         reader->AddVariable("ll_pt",&TMVA.ll_pt);
//         reader->AddVariable("qgl_1qAtanh",&TMVA.qgl_1qAtanh); 
//         reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1); 
//         reader->AddSpectator("ll_mass",&TMVA.ll_mass);
//         reader->AddSpectator("deltaMRel",&TMVA.deltaMRel);
        
        
// // // // //     reader->AddVariable("ll_mass",&TMVA.ll_mass);    
// // // // //     reader->AddVariable("Mqq",&TMVA.Mqq);
// // // // //     
// // // // //     reader->AddVariable("RptHard",&TMVA.Rpt);
// // // // //     reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
// // // // //     
// // // // // //     reader->AddVariable("softActivityEWK_njets5",&TMVA.softActivityEWK_njets5);   // this is int
// // // // //     reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
// // // // // 
// // // // //     
// // // // // //     reader->AddVariable("ll_pt",&TMVA.ll_pt);
// // // // // //     
// // // // // //     
// // // // // //     reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2); 
// // // // // //     reader->AddVariable("W_Pt_virtual1",&TMVA.W_Pt_virtual1);
// // // // //     reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1); 


    
    
//    reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
//     reader->AddVariable("thetastarHtoWW",&TMVA.thetastarHtoWW);    
//     reader->AddVariable("W_Pt_virtual1",&TMVA.W_Pt_virtual1);
   //reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1); 
   // reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2); 
   // 
       
//       reader->AddVariable("ll_mass",&TMVA.ll_mass);
//        reader->AddVariable("Mqq",&TMVA.Mqq);
//    reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
//    reader->AddVariable("q1_eta",&TMVA.q1_eta);
 //   reader->AddVariable("ll_pt",&TMVA.ll_pt);
 
   // reader->AddVariable("met_pt",&TMVA.met_pt);
   
//    reader->AddVariable("qq_pt",&TMVA.qq_pt);
//       reader->AddVariable("RptHard",&TMVA.RptHard);
//       reader->AddVariable("thetastarW1",&TMVA.thetastarW1);
       //reader->AddVariable("theta1",&TMVA.theta1);
//       reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1);
//      reader->AddVariable("W_Pt_virtual1",&TMVA.W_Pt_virtual1);
   
    //reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);
//    reader->AddVariable("softActivityEWK_njets5",&TMVA.softActivityEWK_njets5);
//    reader->AddVariable("btagCMVA",&TMVA.btagCMVA);
//    reader->AddVariable("btagCMVA_second",&TMVA.btagCMVA_second);
//    reader->AddVariable("btagCSV_second",&TMVA.btagCSV_second);
//    reader->AddVariable("btagCSV",&TMVA.btagCSV);
//    reader->AddVariable("qgl_1q",&TMVA.qgl_1q);
//    reader->AddVariable("qgl_2q",&TMVA.qgl_2q);
//    reader->AddVariable("cosThetaStar",&TMVA.cosThetaStar);
//    reader->AddVariable("cosThetaStarAbs",&TMVA.cosThetaStarAbs);
//    reader->AddVariable("cosThetaStarJet",&TMVA.cosThetaStarJet);
//    reader->AddVariable("cosThetaPlane",&TMVA.cosThetaPlane);
//    reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
//    reader->AddVariable("ll_ystar",&TMVA.ll_ystar);
    
//    reader->AddVariable("Jet1q_pt",&TMVA.Jet1q_pt);
   // reader->AddVariable("Jet2q_pt",&TMVA.Jet2q_pt);
//    reader->AddVariable("ll_eta",&TMVA.ll_eta);
    //
//     reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
     // reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
    //reader->AddVariable("W_eta_virtual1",&TMVA.W_eta_virtual1);
    //reader->AddVariable("E_parton1",&TMVA.E_parton1);
     //reader->AddVariable("ll_pt",&TMVA.ll_pt);
    //reader->AddVariable("energytot",&TMVA.energytot);
   //reader->AddVariable("deltaMRel",&TMVA.deltaMRel);
   
//    reader->BookMVA("BDTG", "/gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkimOutputMVA/Classification_BDTG.weights.xml");
//     reader->BookMVA("BDTG", "BDTClassification/Classification_BDTG.weights.xml");
//     reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/training/CMSSW_8_0_28/src/training/TMVA-v4.2.0/test/weights/TMVAClassification_BDTG_nomoremuaxis2jet2q_v25/Classification_BDTG.weights.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/Classification_BDTG.weights_filerootJune_ll_mass_Mqq_DeltaEtaQQ_ll_eta_RptHard_EWKHTsoft_ll_pt_Jet2q_pt.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/trainingForRecoveringJune/Classification_BDTG.weights_JuneOption_ll_mass_Mqq_RptHard_DeltaEtaQQ_llPt_llEta_Jet2qPt_EWKHTsoft.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/januaryTraining/Classification_BDTG.weights_80000Event_200Tree_4Deep_mll_Mqq_RptHard_llZstar_softN5_llPt_Wmass2_Wpt1.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/marchTraining/Classification_BDTG.weights_80000Event_100Tree_2Deep_mll_Mqq_RptHard_llZstar_softN5_Wmass1.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt_Wmass2.xml");


//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_llZstar_softN5_llPt_qgl1qAtanh.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_llZstar_softN5_llPt_qgl1qAtanh_Wmass1.xml");

//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/mayTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt_Wmass2.xml");

        
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/mayTraining/Classification_BDTG_1377SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt.xml");
        
//         TRAINING TESI AGNESE
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/mayTraining/Classification_BDTG_1477SeedGen.weights_30000Event_100Tree_2Deep_SigCirca1_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt.xml");
  
// //         TRAINING NANO 2016
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/octoberTraining/Classification_BDTG_1657SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2.xml");

// //         TRAINING NANO 2016 WITH DY_VBFFILTER
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/octoberTraining/Classification_BDTG_1704SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2_mumujjPz_maxAbsEta_qgl2Atanh.xml");
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/octoberTraining/Classification_BDTG_1704SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2_mumujjPz_maxAbsEta.xml");
            reader->BookMVA("BDTG", "BDTxml/octoberTraining/Classification_BDTG_1704SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2_mumujjPz_maxAbsEta.xml");
//                reader->BookMVA("BDTG", "BDTxml/decemberTraining/Classification_BDTG_1782SeedGen.weights_60000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_softN5_llZstar_llPt_cosMuMuJ2_mumujjPz_maxAbsEta.xml");

// //             TRAINING WITHOUT NSOFT5
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/novemberTraining/Classification_BDTG_1561SeedGen.weights_30000Event_100Tree_2Deep_mll_MqqLog_mumujjPt_DEtajj_llZstar.xml");
            
//         TRAINING FILIPPO
//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/filippoTraining/TMVAClassification_BDTG.weights.xml");
        
        
//         reader->BookMVA("MLP", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_MLP.weights_35000Event_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt_Wmass2_qgl1qAtanh.xml");

    }  


    
   // /afs/cern.ch/user/a/abonavit/private/tesi/training/CMSSW_8_0_28/src/training/TMVA-v4.2.0/test/weights/TMVAClassification_BDTG_nomoremuaxis2jet2q_v25
//"ll_mass", "Mqq","W_Pt_virtual2", "RptHard","thetastarHtoWW","DeltaEtaQQ","W_eta_virtual1","E_parton1","EWKHTsoft","theta1"}



    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    int count5 = 0;
    int count6 = 0;
    int count7 = 0;
    int count8 = 0;
    int count9 = 0;
    int count10 = 0;
    int count11 = 0;
    int count12 = 0;
    int count13 = 0;
    int count14 = 0;

double x = 0;
int numb = 0;

// nentries = 15;
// nentries = 10000;
 for (int entry=0; entry<nentries;++entry){
//   for (int entry=50000; entry<54000;++entry){


//       if ( (entry != 613748) && (entry != 613766) && (entry != 613772)) continue;
//  if ( entry < 21400) continue;
//      std::cout  << "entry   " << entry  << std::endl;
     
     
        for (int n = 0; n < number_LHE_weights_pdf; n++) LHE_weights_pdf_wgt[n] = 1;   // has to be before tree_initial->GetEntry(entry)
        for (int n = 0; n < number_LHE_weights_pdf; n++) TMVA.genweightVECTOR[n] = 0;
//         std::cout << "HERE 1   "<< std::endl;

        if (entry%10000 == 0) {
            if (file_tag.CompareTo("DYJetsToLL_M-105To160-amcatnloFXFX")==0) std::cout << "DY-amcatnloFXFX      ";
            std::cout << "Writing " << entry << "th event" << std::endl;
        }

        tree_initial->GetEntry(entry);
        ++count1;
        
        
//         std::cout << "HERE 0   "<< std::endl; 
//       if ( (evt != 18451352) && (evt != 122582287) ) continue;
//         if (fromNANO && qqMass_skim < 240 ) continue;
            
//       if ( (evt != 328339) && (evt != 81) && (evt != 50643) && (evt != 397525) ) continue;
//         if ( evt != 285 ) continue;
        
        if (nJets > njets ) std::cout << "The problem is: nJets > njets   "<< std::endl;
    
// if (VtypeSim  == -1) cout << "v_typeSim 0 " << v_typeSim << "  " << VtypeSim << endl;
//   cout << entry << "  \t  " <<evt << endl;

//         cout << "\n" << entry << "  \t  " <<evt <<  "  \t 0 ";
//         if ( (evt != 266642512) && (evt != 795769461) && (evt != 57276101) && (evt != 442855523) ) continue;
//         if ( (evt != 51401127) && (evt != 1109295104) && (evt != 387406063) && (evt != 554841123) && (evt != 13393783) && (evt != 1467957513) && (evt != 1376213582) ) continue;
// if (  (evt != 187606427)  &&  (evt != 552232638)  &&  (evt != 1963639737)  &&  (evt != 1095358725)  &&  (evt != 575736556)  &&  (evt != 1320157056)  &&  (evt != 4179710311)  &&  (evt != 1152993139)  &&  (evt != 162577672)  &&  (evt != 283191324)  &&  (evt != 269552083)  &&  (evt != 1397928629)  &&  (evt != 1103027880)  &&  (evt != 186313561) ) continue;// 
//         if ( evt != 529282994) continue;
//         cout << "\n" << entry << "  \t  " <<evt <<  std::endl;
//         
//         std::cout << "HERE 2   " << run << " \t " << luminosityBlock <<  std::endl;
//         std::cout << "HERE 3   " << Jet.pt[0] << std::endl;
        
//         if (  (evt != 575736556)  &&  (evt != 1164776276)  &&  (evt != 2715131981)  &&  (evt != 187606427)  &&  (evt != 1376213582)  &&  (evt != 1467957513)  &&  (evt != 13393783)  &&  (evt != 554841123)  &&  (evt != 387406063)  &&  (evt != 1109295104)  &&  (evt != 51401127)  &&  (evt != 1891982335)  &&  (evt != 1400206699)  &&  (evt != 154868855)  &&  (evt != 902708164)  &&  (evt != 103365171)  &&  (evt != 1963760263)  &&  (evt != 270053610)  &&  (evt != 2751806208)  &&  (evt != 91630457)  &&  (evt != 1652237041)  &&  (evt != 1244994407)  &&  (evt != 530597996)  &&  (evt != 3623607603)  &&  (evt != 171007692)  &&  (evt != 1596511456)  &&  (evt != 20369896)  &&  (evt != 985055631)  &&  (evt != 5826527)  &&  (evt != 430522442)  &&  (evt != 171933071)  &&  (evt != 401106596)  &&  (evt != 2305371621)  &&  (evt != 1851782961)  &&  (evt != 166965309)  &&  (evt != 187172132)  &&  (evt != 534376316)  &&  (evt != 106567679)  &&  (evt != 1343118334)  &&  (evt != 59981121)  &&  (evt != 725969178)  &&  (evt != 1186197122)  &&  (evt != 1994628265)  &&  (evt != 323019882)  &&  (evt != 551365065)  &&  (evt != 2771650456)  &&  (evt != 247817064)  &&  (evt != 1380434242)  &&  (evt != 299783414)  &&  (evt != 66255556)  &&  (evt != 419784069)  &&  (evt != 177379782)  &&  (evt != 887488264)  &&  (evt != 236363582)  &&  (evt != 1785667667)  &&  (evt != 1917126934)  &&  (evt != 580564799) ) continue;
        
//         genJetsWithoutLeptonsP4.clear();
        genJetsWithoutLeptonsP4_pt.clear();
        genJetsWithoutLeptonsP4_phi.clear();
        genJetsWithoutLeptonsP4_eta.clear();
        genJetsWithoutLeptonsP4_mass.clear();
        ngenJetsWithoutLeptonsP4 = 0;
        
        
        TMVA.btagCMVA_leading = 0;
        TMVA.btagCMVA = 0;
        TMVA.btagCSV=0;
        TMVA.btagCMVA_second = 0;
        TMVA.btagCSV_second=0;
        TMVA.softLeadingJet_pt=0;
        TMVA.cosThetaStar=0;
        TMVA.cosThetaStarAbs=0;
        
        
        TMVA.deltaMRel=0;
        TMVA.normalizedDistance_from_mH=0;
        TMVA.lepton1_pt=0;
        TMVA.lepton2_pt=0;
        TMVA.ll_mass = 0;
        TMVA.ll_mass_biasUp = 0;
        TMVA.ll_mass_biasDown = 0;
        TMVA.ll_mass_resolution = 0;
        TMVA.Mqq=0;
	TMVA.qq_pt=0;
	TMVA.DeltaEtaQQ=0;
	TMVA.Jet3_pt=0;
	TMVA.Jet3_eta=0;
	TMVA.axis2_jet1=0;
	TMVA.axis2_jet2=0;
	TMVA.Jet2q_pt=0; 
	TMVA.Jet2q_leadTrackPt=0;
	TMVA.Jet1q_pt=0;
        TMVA.Jet1q_eta=0;
        TMVA.Jet2q_eta=0;
        
        TMVA.Jet1q_chEmEF=0;
        TMVA.Jet2q_chEmEF=0;
        TMVA.Jet1q_neEmEF=0;
        TMVA.Jet2q_neEmEF=0;
        
        
	TMVA.Jet1q_leadTrackPt=0;
	TMVA.ll_ystar=0;
	TMVA.ll_zstar=0;
	TMVA.Rpt=0;
	TMVA.ll_pt=0;
        TMVA.ll_pt_Log=0;
        TMVA.ll_eta=0;
	TMVA.qgl_1q=0;
	TMVA.qgl_2q=0;
        TMVA.qgl_1qAtanh=0;
	TMVA.qgl_2qAtanh==0;
        
        TMVA.Jet1q_puId = -1;
        TMVA.Jet2q_puId = -1;
        
        TMVA.Jet1q_pt_jerUp = 0;
        TMVA.Jet1q_pt_jerDown = 0;
        TMVA.Jet2q_pt_jerUp = 0;
        TMVA.Jet2q_pt_jerDown = 0;
        
        for (int n=0; n<30; n++) {
            TMVA.Jet_pt[n] = 0;
            TMVA.genJetIdx[n] = -1;
        }
        
        TMVA.genJetPt_match1 = 0;
        TMVA.genJetPt_match2 = 0;
        
        TMVA.genJetIdx0 = -1;
        TMVA.genJetIdx1 = -1;
        TMVA.genJetIdx2 = -1;
        TMVA.genJetIdx0_genJetMassWithoutLepton = -1;
        TMVA.genJetIdx1_genJetMassWithoutLepton = -1;
        TMVA.genJetIdx2_genJetMassWithoutLepton = -1;
        
        TMVA.minAbsEta=10;
        TMVA.maxAbsEta=0;
        TMVA.countJet25=0;
        
//         std::cout << "HERE 0   "<< std::endl;
        TMVA.cosThetaPlaneAtanh =0;
        TMVA.absCosThetaStarJetAtanh =0;
        TMVA.X_parton1Log=0;
        TMVA.X_parton2Log=0;
        TMVA.W_mass_virtual1Log=0;
        TMVA.W_mass_virtual2Log=0;
        TMVA.W_Pt_virtual1Log=0;
        TMVA.W_Pt_virtual2Log=0;
               
        
        TMVA.quark1_pt=0;
        TMVA.quark1_eta=0;
        TMVA.quark1_phi=0;
        TMVA.quark1_mass=0;
            
        TMVA.quark2_pt=0;
        TMVA.quark2_eta=0;
        TMVA.quark2_phi=0;
        TMVA.quark2_mass=0;
        
        TMVA.genMuon1_pt=0;
        TMVA.genMuon1_eta=0;
        TMVA.genMuon1_phi=0;
        TMVA.genMuon1_mass=0;        
        TMVA.genMuon2_pt=0;
        TMVA.genMuon2_eta=0;
        TMVA.genMuon2_phi=0;
        TMVA.genMuon2_mass=0;
        
        
        TMVA.qqMass_skim=0;
        TMVA.MuonMass_skim=0;
        TMVA.genJetMass_leading=0;
        
        TMVA.LHE_weights_scale_wUp = 1;           
        TMVA.LHE_weights_scale_wDown = 1;           
    
        TMVA.softActivityEWK_njets2=0;
        TMVA.softActivityEWK_njets5=0;
	TMVA.softActivityEWK_njets10=0;
        TMVA.randomVariable=0;
        TMVA.xSection=0;
        TMVA.weightMVA=0;
        TMVA.genweight=0;
        TMVA.genweight_noQGLcorrection=0;
        TMVA.deltaM=0;
        
        TMVA.BDToutput=0;
        
        int nIsolatedElectrons=0;
//         std::cout << "HERE 4   "<< std::endl;
        
        float maxBTagCSV = -0.2;
        float maxBTagCMVA = -1.2;
        float maxSecondBTagCSV = -0.2;
        float maxSecondBTagCMVA = -1.2;
        
        float maxBTagCMVA_leading = -1.2;
        float maxSecondBTagCMVA_leading = -1.2;
        
        lumiToWrite = luminosityBlock;
        runToWrite = run;
        
//           std::cout << "HERE 5   "<< std::endl;
        
        
//     if (fromNANO) {
// 
//         nGenLepRecovered = 0;
//         nGenLep = 0;
//         for (int n=0; n < nGenPart; n++) {
//             if (abs(GenPart_pdgId[n]) == 13) {
//                 GenLep_pt[nGenLep]    = GenPart_pt[n];
//                 GenLep_eta[nGenLep]   = GenPart_eta[n];
//                 GenLep_phi[nGenLep]   = GenPart_phi[n];
//                 GenLep_mass[nGenLep]  = GenPart_mass[n];  
//                 GenLep_pdgId[nGenLep] = GenPart_pdgId[n];
//                 nGenLep++;
//             }
//         }
//    
//     }
          
          
            
        lheNj = (int) lheNj_f;
        
        
        
// CORRECTIONS DUE TO L1 TRIGGER         THEY ARE BEFORE JES AND JER PT CORRECTIONS
    float forwardJetCorrection = 1.;
    if ((data!=1)) {
        if (fromNANO) {
            for (int i=0;i<nJets;i++){
                float jetWeight = 1.;
                if (abs(Jet.eta[i]) > 2 && abs(Jet.eta[i]) < 3.5) {
//                     std::cout << Jet.pt[i] << " \t " << abs(Jet.eta[i])  << " \t " << h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),Jet.pt[i]) << " \t " <<  h_forwardJetWeight->GetEfficiency(h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),Jet.pt[i]))<< std::endl;
                    double ptToUse =  Jet.pt[i];
                    float w = 0;
                    do {
                        w = h_forwardJetWeight->GetEfficiency(h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),ptToUse));
//                         if (w == 0 && ptToUse > 100.) std::cout << ptToUse << " \t " << abs(Jet.eta[i])  << " \t " << h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),ptToUse) << " \t " <<  h_forwardJetWeight->GetEfficiency(h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),ptToUse))<< std::endl;
                        ptToUse = ptToUse * 0.9;
                    }
                    while (w == 0 && ptToUse > 100.);
                        
                    jetWeight = 1. - w;
                    forwardJetCorrection = forwardJetCorrection * jetWeight;
                }
            }
            
            
            for (int i=0;i<nJets;i++){
                float jetWeight = 1.;
                if (abs(Jet.eta[i]) < 3.5) {
                    jetWeight = 1. - h_forwardJetWeight->GetEfficiency(h_forwardJetWeight->FindFixBin(abs(Jet.eta[i]),Jet.pt[i]));
                    forwardJetCorrection = forwardJetCorrection * jetWeight;
                }
            }
            
            
        }
    }
               
               genweight = genweight * forwardJetCorrection;   // <---------------------------------------------------- PREFIRING TO APPLY! DO not comment this line
    
        
               
//  std::cout << "HERE 1  "<< std::endl;
    if ((data!=1)) {
            
        if (fromNANO) {
            for (int i=0;i<nJets;i++){
                if (whichJESWeight==1)  {Jet.pt[i]= Jet.pt_JESup[i];    Jet.mass[i]= Jet.mass_JESup[i];}	
                if (whichJESWeight==2)  {Jet.pt[i]= Jet.pt_JESdown[i];  Jet.mass[i]= Jet.mass_JESdown[i];}	
                if (whichJERWeight==1)  {Jet.pt[i]= Jet.pt_JERup[i];    Jet.mass[i]= Jet.mass_JERup[i];}	
                if (whichJERWeight==2)  {Jet.pt[i]= Jet.pt_JERdown[i];  Jet.mass[i]= Jet.mass_JERdown[i];}	
            }
            if (whichJESWeight==1)  {met_pt= met_pt_JESup;    met_phi= met_phi_JESup;}	
            if (whichJESWeight==2)  {met_pt= met_pt_JESdown;  met_phi= met_phi_JESdown;}	
            if (whichJERWeight==1)  {met_pt= met_pt_JERup;    met_phi= met_phi_JERup;}	
            if (whichJERWeight==2)  {met_pt= met_pt_JERdown;  met_phi= met_phi_JERdown;}	
        }
        else {        
            for (int i=0;i<nJets;i++){
                if (whichJESWeight==1) if(Jet.JEC_corr_up[i] > 0 && Jet.JEC_corr[i] > 0)    {Jet.pt[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i];    Jet.mass[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i];}	
                if (whichJESWeight==2) if(Jet.JEC_corr_down[i] > 0 && Jet.JEC_corr[i] > 0)  {Jet.pt[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];  Jet.mass[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];}	
                if (whichJERWeight==1) if(Jet.JER_corr_up[i] > 0 && Jet.JER_corr[i] > 0)    {Jet.pt[i]*= Jet.JER_corr_up[i]/Jet.JER_corr[i];    Jet.mass[i]*= Jet.JER_corr_up[i]/Jet.JER_corr[i];}	
                if (whichJERWeight==2) if(Jet.JER_corr_down[i] > 0 && Jet.JER_corr[i] > 0)  {Jet.pt[i]*= Jet.JER_corr_down[i]/Jet.JER_corr[i];  Jet.mass[i]*= Jet.JER_corr_down[i]/Jet.JER_corr[i];}	
            }
        }
    }
        

        if (!fromNANO && JSON!=1) {
            continue;
        }

	//	if ((file_tag.CompareTo("EWK_LLJJ")==0) &&(evt%2!=0)) continue;

// 		if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
	//	if (region.CompareTo("el")==0) if (!(v_type==1)) continue;


	
	
    


    
    
    
//     std::cout <<  forwardJetCorrection << std::endl << std::endl ;
        
	
		if (data==1) PU=1.;
		else {
                    if ((PUWeight_str.CompareTo("none")==0)||(PUWeight_str.CompareTo("nom")==0))    PU=puweight;
                    if (PUWeight_str.CompareTo("up")==0)                                            PU=puweightUp;
                    if (PUWeight_str.CompareTo("down")==0)                                          PU=puweightDown;
                }
		genweight0 = genweight/TMath::Abs(genweight);
                genweight0 = genweight;
		if (fromNANO) genweight=genweight*PU;   
                else genweight=genweight/TMath::Abs(genweight)*PU;   
                
                
                
// 		cout << "PU  " << PU << " \t   genweight     " << genweight << " \t xsec[file_tag]   " << xsec[file_tag] << " \t events_generated     " << events_generated <<  endl; 

//		genweight/=events_generated/xsec[file_tag]; 
                if (fromNANO) {
                    
                    
                    if  ((data!=1) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==1)  {
                        genweight*= ( LHE_weights_scale_wgt[8] > 0.001 ? LHE_weights_scale_wgt[8] : 1.);
                    }
                    if  ((data!=1) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0) ) if (whichQCDScaleWeight==2)  {
                        genweight*= ( LHE_weights_scale_wgt[0] > 0.001 ? LHE_weights_scale_wgt[0] : 1.);
                    }
                }
//                 else {
//                     if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==1)  genweight*=LHE_weights_scale_wgt[4];
//                     if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==2)       genweight*=LHE_weights_scale_wgt[5];
//                     
//                 }
                
                TMVA.LHE_weights_scale_wUp = ( LHE_weights_scale_wgt[8] > 0.001 ? LHE_weights_scale_wgt[8] : 1.);
                TMVA.LHE_weights_scale_wDown = ( LHE_weights_scale_wgt[0] > 0.001 ? LHE_weights_scale_wgt[0] : 1.);
                
                

                

		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];


 
		if  ((file_tag.CompareTo("DYJetstoLL_HT100")==0)) if (lheHT>100) continue;  
//		if  ((file_tag.CompareTo("DYJetstoLL_Pt-100_amc")==0)) if (lheV_pt>100) continue;  
		if  ((file_tag.CompareTo("WJetsToLNu_HT100")==0)) if (lheHT>100) continue;  


		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector qq;
		int good_jets = 0;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		float jet3_pt = 0;
		float jet3_eta=0;


		/////////////////////////////////////////////////////////////////////
		/////////////////////////PRESELECTION////////////////////////////////
		////////////////////////////////////////////////////////////////////
		cut_flow[0]+=genweight;
                

                

                
        Int_t pdgId_mu1;
        Int_t pdgId_mu2;  
        TLorentzVector Genlepton1;
        TLorentzVector Genlepton2;
        TLorentzVector GenLeptons;
          
        
        if(nGenLep==2){
            pdgId_mu1=GenLep_pdgId[0];
            pdgId_mu2=GenLep_pdgId[1];
            Genlepton1.SetPtEtaPhiM(GenLep_pt[0],GenLep_eta[0], GenLep_phi[0], GenLep_mass[0]);
            Genlepton2.SetPtEtaPhiM(GenLep_pt[1],GenLep_eta[1], GenLep_phi[1], GenLep_mass[1]);
    	}
    			
    	if(nGenLepRecovered==2){	
            pdgId_mu1=GenLepRecovered_pdgId[0];
            pdgId_mu2=GenLepRecovered_pdgId[1];
            Genlepton1.SetPtEtaPhiM(GenLepRecovered_pt[0],GenLepRecovered_eta[0], GenLepRecovered_phi[0], GenLepRecovered_mass[0] );
            Genlepton2.SetPtEtaPhiM(GenLepRecovered_pt[1],GenLepRecovered_eta[1], GenLepRecovered_phi[1], GenLepRecovered_mass[1] );
    	}
    	

    	

    	
    	
        GenLeptons=Genlepton1+Genlepton2;
        Float_t gen_mass = GenLeptons.M();
		Float_t gen_pt = GenLeptons.Pt();
		//Float_t gen_eta = GenLeptons.Eta();
		Float_t gen_phi = GenLeptons.Phi();
                

                
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////  FIND HIGGS SISTERS  ////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
                
        for (int n=0; n<nGenPart; n++) { 
            if (GenPart_pdgId[n]==25 and n<nGenPart-2) { 
                int leadingIdx = 3;
                int subleadingIdx = 4;
                if (GenPart_pt[4]>GenPart_pt[3]) {
                    leadingIdx = 4;
                    subleadingIdx = 3;
                }
                TMVA.quark1_pt   = GenPart_pt[leadingIdx];
                TMVA.quark1_eta  = GenPart_eta[leadingIdx];
                TMVA.quark1_phi  = GenPart_phi[leadingIdx];
                TMVA.quark1_mass = GenPart_mass[leadingIdx];
                TMVA.quark2_pt   = GenPart_pt[subleadingIdx];
                TMVA.quark2_eta  = GenPart_eta[subleadingIdx];
                TMVA.quark2_phi  = GenPart_phi[subleadingIdx];
                TMVA.quark2_mass = GenPart_mass[subleadingIdx];
            
            }
        }
                
                

                
//   std::cout << "HERE 2   "<< std::endl; 
                
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////// FILTER GENJET INVARIANT MASS ////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<TLorentzVector> genLeptons;
        std::vector<TLorentzVector> genJet;
        int ngenLepIdx = 0;
        if (data!=1) {
            
            if (fromNANO) {
//                 std::cout << "nGenPart \t " << nGenPart<<  std::endl;
                
                for (int n=0; n<nGenPart; n++) {
                    if (file_tag.CompareTo("EKW_LLJJ_pythia8")==0) continue;
                    if ( ( (abs(GenPart_pdgId[n])==11) || (abs(GenPart_pdgId[n])==13) || (abs(GenPart_pdgId[n])==15)) && GenPart_statusFlags[n]%256>=128 ) {
                        TLorentzVector genLepton_tmp;
                        genLepton_tmp.SetPtEtaPhiM(GenPart_pt[n],GenPart_eta[n], GenPart_phi[n], GenPart_mass[n]);
                        genLeptons.push_back(genLepton_tmp);
                        
//                         std::cout << GenPart_pdgId[n] << " \t " << GenPart_pt[n] << "     \t " << GenPart_eta[n] << "   \t " << GenPart_mass[n] << " \t " << ngenLepIdx <<  std::endl;
                        
                        if(ngenLepIdx==1){
                            pdgId_mu2=GenPart_pdgId[n];
                            pdgId_mu2=GenPart_statusFlags[n];
                            Genlepton2 = genLepton_tmp;
                            ngenLepIdx++;
                        }
                        
                        if(ngenLepIdx==0){
                            pdgId_mu1=GenPart_pdgId[n];
                            pdgId_mu1=GenPart_statusFlags[n];
                            Genlepton1 = genLepton_tmp;
                            ngenLepIdx++;
                        }
                        
                    }
                }
                
            }
            else {
                for (int n=0; n<nGenLep; n++) {
                    TLorentzVector genLepton_tmp;
                    genLepton_tmp.SetPtEtaPhiM(GenLep_pt[n],GenLep_eta[n], GenLep_phi[n], GenLep_mass[n]);
                    genLeptons.push_back(genLepton_tmp);
                }

                for (int n=0; n<nGenLepRecovered; n++) {
                    TLorentzVector genLepton_tmp;
                    genLepton_tmp.SetPtEtaPhiM(GenLepRecovered_pt[n],GenLepRecovered_eta[n], GenLepRecovered_phi[n], GenLepRecovered_mass[n]);
                    genLeptons.push_back(genLepton_tmp);
                }
            }
            
            
            if (genLeptons.size()>0) {
                TMVA.genMuon1_pt = genLeptons[0].Pt();
                TMVA.genMuon1_eta = genLeptons[0].Eta();
                TMVA.genMuon1_phi = genLeptons[0].Phi();
                TMVA.genMuon1_mass = genLeptons[0].M();
            }
            if (genLeptons.size()>1) {
                TMVA.genMuon2_pt = genLeptons[1].Pt();
                TMVA.genMuon2_eta = genLeptons[1].Eta();
                TMVA.genMuon2_phi = genLeptons[1].Phi();
                TMVA.genMuon2_mass = genLeptons[1].M();
            }

            
            for (int n=0; n<nGenJet; n++) {
                TLorentzVector genJet_tmp;
                genJet_tmp.SetPtEtaPhiM(GenJet_pt[n],GenJet_eta[n], GenJet_phi[n], GenJet_mass[n]);
                genJet.push_back(genJet_tmp);
            }
        }
        
//         float genJetMass = VBFFilter(genLeptons, genJet);
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        


                        
                        
        // This redefinition are here for nanoAOD skim
        GenLeptons=Genlepton1+Genlepton2;
        gen_mass = GenLeptons.M();
        gen_pt = GenLeptons.Pt();
//         gen_eta = GenLeptons.Eta();
        gen_phi = GenLeptons.Phi();


        // std::cout << "HERE 4   " << nJets << std::endl;
        
        for (int i=0;i<nJets;i++){

//             std::cout << entry << " \t  "  << i << " \t  "  << (Jet.id[i]>0) << " \t  "  << (Jet.puId[i]>0) << " \t  "  << (abs(Jet.eta[i])<4.7) << " \t  "  << !((abs(Jet.eta[i])>2.5)&&(Jet.puId[i]<7)) << " \t  "  << " \t  "  << (Jet.muonIdx1[i]>-1) <<  " \t " <<  (selLeptons_relIso04[max(Jet.muonIdx1[i],0)]*(Jet.muonIdx1[i]>-1) <0.25) <<  " \t " <<  (abs(selLeptons_dz[max(Jet.muonIdx1[i],0)])*(Jet.muonIdx1[i]>-1) < 0.2) <<  " \t " << (abs(selLeptons_dxy[max(Jet.muonIdx1[i],0)])*(Jet.muonIdx1[i]>-1) < 0.05) << " \t " << (Jet.electronIdx1[i]>-1) <<  " \t " <<  (Electron_pfRelIso03_all[max(Jet.electronIdx1[i],0)]*(Jet.electronIdx1[i]>-1) <0.25) <<  " \t " << (abs(Electron_dz[max(Jet.electronIdx1[i],0)])*(Jet.electronIdx1[i]>-1) < 0.2) <<  " \t " << (abs(Electron_dxy[max(Jet.electronIdx1[i],0)])*(Jet.electronIdx1[i]>-1) < 0.05)  << std::endl;

//             if (i+1==nJets) std::cout <<std::endl;
            
//             std::cout << "HERE 4   " << nJets << " \t " << Jet.pt[i] << std::endl;
            TLorentzVector jet0;
//             if (!((Jet.id[i]>2)&&(Jet.puId[i]>0))) continue;     
            if (!((Jet.id[i]>0)&&(Jet.puId[i]>0))) continue;
//             if (!((Jet.id[i]>0)&&(Jet.puId[i]>4))) continue;
//             if (!((Jet.id[i]>0)&&(Jet.puId[i]>6))) continue;
            
            
            if ((abs(Jet.eta[i])>2.5)&&(Jet.puId[i]<7)) continue;
//             if ((abs(Jet.eta[i])>2.5)&&(Jet.puId[i]<6)) continue;
            
            if (abs(Jet.eta[i])>4.7) continue;
             
            if (fromNANO) {
                if (!andrewSelection) {
                if (Jet.muonIdx1[i]>-1 && selLeptons_relIso04[Jet.muonIdx1[i]] <0.25 && abs(selLeptons_dz[Jet.muonIdx1[i]]) < 0.2 && abs(selLeptons_dxy[Jet.muonIdx1[i]]) < 0.05) continue;
                if (Jet.muonIdx2[i]>-1 && selLeptons_relIso04[Jet.muonIdx2[i]] <0.25 && abs(selLeptons_dz[Jet.muonIdx2[i]]) < 0.2 && abs(selLeptons_dxy[Jet.muonIdx2[i]]) < 0.05) continue;
                if (Jet.electronIdx1[i]>-1 && Electron_pfRelIso03_all[Jet.electronIdx1[i]] <0.25 && abs(Electron_dz[Jet.electronIdx1[i]]) < 0.2 && abs(Electron_dxy[Jet.electronIdx1[i]]) < 0.05) continue;
                if (Jet.electronIdx2[i]>-1 && Electron_pfRelIso03_all[Jet.electronIdx2[i]] <0.25 && abs(Electron_dz[Jet.electronIdx2[i]]) < 0.2 && abs(Electron_dxy[Jet.electronIdx2[i]]) < 0.05) continue;
                }
                else {
                if(abs(Jet.eta[i])>4.7) continue;
                
                if (Jet.photonIdx1[i]>-1 && Photon_pt[Jet.photonIdx1[i]] > 15 && abs(Photon_eta[Jet.photonIdx1[i]])<2.5 && Photon_mvaID_WP90[Jet.photonIdx1[i]] && Photon_pfRelIso03_all[Jet.photonIdx1[i]]<0.25) continue;
                if (Jet.photonIdx2[i]>-1 && Photon_pt[Jet.photonIdx2[i]] > 15 && abs(Photon_eta[Jet.photonIdx2[i]])<2.5 && Photon_mvaID_WP90[Jet.photonIdx2[i]] && Photon_pfRelIso03_all[Jet.photonIdx2[i]]<0.25) continue;

                if (Jet.muonIdx1[i]>-1 && selLeptons_relIso04[Jet.muonIdx1[i]] <0.25 && selLeptons_softIdPOG[Jet.muonIdx1[i]] && selLeptons_pt[Jet.muonIdx1[i]]>10 && abs(selLeptons_eta[Jet.muonIdx1[i]])<2.4 ) continue;
                if (Jet.muonIdx2[i]>-1 && selLeptons_relIso04[Jet.muonIdx2[i]] <0.25 && selLeptons_softIdPOG[Jet.muonIdx2[i]] && selLeptons_pt[Jet.muonIdx2[i]]>10 && abs(selLeptons_eta[Jet.muonIdx2[i]])<2.4 ) continue;

                if (Jet.electronIdx1[i]>-1 && Electron_pfRelIso03_all[Jet.electronIdx1[i]] <0.25 && Electron_mvaSpring16GP_WP90[Jet.electronIdx1[i]] && Electron_pt[Jet.electronIdx1[i]]>10 && abs(Electron_eta[Jet.electronIdx1[i]])<2.5 ) continue;
                if (Jet.electronIdx2[i]>-1 && Electron_pfRelIso03_all[Jet.electronIdx2[i]] <0.25 && Electron_mvaSpring16GP_WP90[Jet.electronIdx2[i]] && Electron_pt[Jet.electronIdx2[i]]>10 && abs(Electron_eta[Jet.electronIdx2[i]])<2.5 ) continue;
                }
            }
                
            jet0.SetPtEtaPhiM(Jet.pt[i],Jet.eta[i],Jet.phi[i],Jet.mass[i]);

            if(Jet.pt[i]>20) {
                
                if (maxBTagCSV < Jet.btagCSV[i]) {
                    maxSecondBTagCSV = maxBTagCSV;
                    maxBTagCSV = Jet.btagCSV[i];
                }
                if (maxBTagCMVA < Jet.btagCMVA[i]) {
                    maxSecondBTagCMVA = maxBTagCMVA;
                    maxBTagCMVA = Jet.btagCMVA[i];
                }

                if (maxBTagCSV > Jet.btagCSV[i] && maxSecondBTagCSV < Jet.btagCSV[i]) maxSecondBTagCSV = Jet.btagCSV[i];
                if (maxBTagCMVA > Jet.btagCMVA[i] && maxSecondBTagCMVA < Jet.btagCMVA[i]) maxSecondBTagCMVA = Jet.btagCMVA[i];
            }
            if (maxBTagCMVA_leading < -1) maxBTagCMVA_leading = Jet.btagCMVA[i];
            if (maxBTagCMVA_leading > -1 && maxSecondBTagCMVA_leading < -1) maxSecondBTagCMVA_leading = Jet.btagCMVA[i];
//             
//             if (maxBTagCSV < 0) maxBTagCSV = Jet.btagCSV[i];
//             if (maxBTagCSV > 0 && maxSecondBTagCSV < 0) maxSecondBTagCSV = Jet.btagCSV[i];
                
            bool condition = false;
            if (good_jets < 2) condition = true;
            else if (jet0.Eta() < max(jets_pv[0].Eta(), jets_pv[1].Eta()) && jet0.Eta() > min(jets_pv[0].Eta(), jets_pv[1].Eta())) condition = true;  
            condition = true;  

//             if (good_jets == 0 && (Jet.jetIdx1 != i) ) {
//                 std::cout << entry << " \tjetIdx1   "  << Jet.jetIdx1 << " \t " << i  << " \t " << nJets << " \t " << selLeptons_pt[max(Jet.muonIdx1[i],0)]*(Jet.muonIdx1[i]>-1) <<  " \t " << Jet.muonIdx1[i] <<  " \t " << Jet.pt[i]  << std::endl;
//             }
//             if (good_jets == 1 && (Jet.jetIdx2 != i) ) {
//                 std::cout << entry << " \tjetIdx2   "  << Jet.jetIdx2 << " \t " << i  << " \t " << nJets << " \t "  << selLeptons_pt[max(Jet.muonIdx1[i],0)]*(Jet.muonIdx1[i]>-1) <<  " \t " << Jet.muonIdx1[i] <<  " \t " << Jet.pt[i]  << std::endl;
//             }
            
            if (condition) {
//                 std::cout << "HERE 4.4   " << nJets << " \t " << Jet.pt[i] << std::endl;
//                 if (good_jets == 0 && Jet.jetIdx1 != i) std::cout << " \tjetIdx1   "  << Jet.jetIdx1 << " \t " << i  << " \t "  << selLeptons_relIso04[Jet.muonIdx1[i]] <<  " \t " << Jet.muonIdx1[i] <<  " \t " << Jet.pt[i]  << std::endl;
//                 if (good_jets == 1 && Jet.jetIdx2 != i) std::cout << " \tjetIdx2   "  << Jet.jetIdx2 << " \t " << i  << " \t "  << selLeptons_relIso04[Jet.muonIdx2[i]] <<  " \t " << Jet.muonIdx2[i] <<  " \t " << Jet.pt[i]  << std::endl;
                jets_pv.push_back(jet0);
                jets_indices.push_back(i);
                good_jets++;
            }
        }

        if (good_jets==1) {
                jets_pv.push_back(jets_pv[0]);
                jets_indices.push_back(jets_indices[0]);
                good_jets++;
        }


            x += genweight0*xsec[file_tag]/events_generated;
            numb += 1;
        
            
            
            
            
            
//             std::vector<float> jet_MVA_value;
//             fileMVA.cd();
//                 for (int i=0;i<good_jets;i++) {
//                     
//                     if(jets_pv[i].Pt()<25) {jet_MVA_value.push_back(-1); continue;}
//                         
// //                     if ((indexJetHiggsSister1 == i) || (indexJetHiggsSister2 == i)) JET_MVA_fromVBF = 1;
//                     else JET_MVA_fromVBF = 0;
//                     JET_MVA_index = i;
//                     JET_MVA_pt = jets_pv[i].Pt();
//                     JET_MVA_p = jets_pv[i].P();
//                     JET_MVA_eta = jets_pv[i].Eta();
//                     JET_MVA_phi = jets_pv[i].Phi();
//                     JET_MVA_mass = jets_pv[i].M();
//                     JET_MVA_id = Jet.id[jets_indices[i]];
//                     JET_MVA_puId = Jet.puId[jets_indices[i]];
//                     JET_MVA_Njets = TMVA.countJet25;
//                     JET_MVA_qgl = Jet.qgl[jets_indices[i]];
//                     JET_MVA_area = Jet.axis2[jets_indices[i]];
//                     JET_MVA_mult = Jet.mult[jets_indices[i]];
//                     
//                     JET_MVA_pt1 = 0;
//                     JET_MVA_phi1 = 0;
//                     JET_MVA_eta1 = 0;
//                     JET_MVA_pt2 = 0;
//                     JET_MVA_eta2 = 0;
//                     JET_MVA_phi2 = 0;
//                     JET_MVA_pt3 = 0;
//                     JET_MVA_eta3 = 0;
//                     JET_MVA_phi3 = 0;
//                     int idxToWrite = 1;
//                     for (int j=0;j<good_jets;j++) {
//                         if ((jets_pv[j].Pt()<25) || (j==i)) continue;
//                         
//                         if (idxToWrite==1) {JET_MVA_pt1 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta1 = jets_pv[j].Eta(); JET_MVA_phi1 = jets_pv[j].Phi();}
//                         if (idxToWrite==2) {JET_MVA_pt2 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta2 = jets_pv[j].Eta(); JET_MVA_phi2 = jets_pv[j].Phi();}
//                         if (idxToWrite==3) {JET_MVA_pt3 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta3 = jets_pv[j].Eta(); JET_MVA_phi3 = jets_pv[j].Phi();}
//                         idxToWrite++;
//                     }
// 
//         
// //                     JET_MVA_higgh_pt = Zll_pt;
// //                     JET_MVA_higgh_eta = Zll_eta;
// //                     JET_MVA_higgh_phi = Zll_pt;
// //                     JET_MVA_higgh_mass = Zll_mass;
//                     
//                     
//                     JET_MVA_etaIdx = -JET_MVA_NjetsEtaMinus;
//                     for (int j=0;j<good_jets;j++)        
//                         if (jets_pv[j].Pt()>25 && JET_MVA_eta>jets_pv[j].Eta() && i!=j) JET_MVA_etaIdx++;
//         
// //                     if (JET_MVA_eta>Zll_eta) JET_MVA_etaIdx++;
//                     if (JET_MVA_etaIdx>=0) JET_MVA_etaIdx++;
//         
//         
// //                     if (TMVA.countJet25>2) std::cout  << i << "    " << Zll_eta << " \t "  << JET_MVA_Njets - JET_MVA_NjetsEtaMinus - JET_MVA_NjetsEtaPlus << " \t "  << JET_MVA_NjetsEtaMinus << " \t "  << JET_MVA_NjetsEtaPlus << " \t "  << JET_MVA_eta << " \t "  << JET_MVA_eta1 << " \t "  << JET_MVA_eta2 << " \t "  << JET_MVA_eta3 << " \t "  << JET_MVA_etaIdx << " \t "  << std::endl;
// 
//                     
//                     JET_MVA_index_float = JET_MVA_index;
//                     JET_MVA_puId_float = JET_MVA_puId;
//                     JET_MVA_Njets_float = JET_MVA_Njets;
//                     JET_MVA_etaIdx_float = JET_MVA_etaIdx;
//                     JET_MVA_NjetsEtaMinus_float = JET_MVA_NjetsEtaMinus;
//                     JET_MVA_NjetsEtaPlus_float = JET_MVA_NjetsEtaPlus;
//         
//                     jet_MVA_value.push_back(readerJetMVA->EvaluateMVA("BDTG"));
// //                     if (TMVA.countJet25>2) std::cout << "AppenaValutato " << i << " \t " << readerJetMVA->EvaluateMVA("BDTG") << std::endl;
// //                     jet_MVA_value.push_back(0);
//                 
// //                     if (TMVA.countJet25>2) treeMVAjet->Fill(); 
//             }
            
            
            
/*        int higherJetMVA = 0;
        int secondJetMVA = 1;
        float higherMVA = -0.99;
        for (int i=0;i<good_jets;i++) {
            if (jet_MVA_value[i] > higherMVA) {
                higherMVA = jet_MVA_value[i];
                higherJetMVA = i;
            }
        }
        
        higherMVA = -0.99;
        for (int i=0;i<good_jets;i++) {
            if (jet_MVA_value[i] > higherMVA && i!=higherJetMVA) {
                higherMVA = jet_MVA_value[i];
                secondJetMVA = i;
            }
        }    */             
            
            

/*        if(good_jets>0) jets_indices[0] = higherJetMVA;
        if(good_jets>1) jets_indices[1] = secondJetMVA; */         
            
            
//             if (entry > nentries -  2)  std::cout << "HERE   "  << x  << " \t " << x*35900  << " \t " << numb << " \t " << std::endl;
            
            
//         for (int n = 0; n < jets_pv.size(); n++)  std::cout << n <<  " \t   "  << jets_pv[n].Pt() << " \t " << jets_pv[n].Eta()  << " \t "  << jets_pv[n].Phi() <<  " \t " << jets_pv[n].M()   << std::endl;
 
//    std::cout << "HERE 5   "  <<  genweight << " \t " << LHE_weights_scale_wgt[8]  << " \t " << xsec[file_tag] << " \t " << events_generated  << std::endl;
        
//         hSelectionCuts->Fill(0., genweight*xsec[file_tag]/events_generated);
//         hSelectionCuts2->Fill(0., 1.64871e-08);
//         if (good_jets<2) continue;
//         else  ++count2;    
        
        
//         if (!(jets_indices[0]==Jet.jetIdx1 && jets_indices[1]==Jet.jetIdx2) ) std::cout << evt << " \t " <<  jets_indices[0] << " \t " << jets_indices[1]  << " \t " << Jet.jetIdx1 << " \t " << Jet.jetIdx2  <<  std::endl;


        if (maxSecondBTagCMVA_leading > maxBTagCMVA_leading) {float tmpBTag = maxSecondBTagCMVA_leading; maxSecondBTagCMVA_leading = maxBTagCMVA_leading; maxBTagCMVA_leading = tmpBTag;}
//         if (maxSecondBTagCSV > maxBTagCSV) {float tmpBTag = maxSecondBTagCSV; maxSecondBTagCSV = maxBTagCSV; maxBTagCSV = tmpBTag;}


//         if(Jet.qgl[jets_indices[0]]<-0.5) std::cout << good_jets << " \t " << entry << " \t jets_indices[0]:  " << jets_indices[0] << " \t " << Jet.qgl[jets_indices[0]] << std::endl;
//         if(Jet.qgl[jets_indices[1]]<-0.5) std::cout << good_jets << " \t " << entry << " \t jets_indices[1]:  " << jets_indices[1] << " \t " << Jet.qgl[jets_indices[1]] << std::endl;

        if(good_jets>0 &&  Jet.qgl[jets_indices[0]]<-0.5) Jet.qgl[jets_indices[0]] = 0.0000001;
        if(good_jets>1 &&  Jet.qgl[jets_indices[1]]<-0.5) Jet.qgl[jets_indices[1]] = 0.0000001;

        
// std::cout << evt << " \t "  << " \t " << jets_indices[0] << " \t "  << jets_indices[1] << " \t " << (Jet.pt_JERdown[jets_indices[0]]-Jet.pt[jets_indices[0]])/(Jet.pt[jets_indices[0]]-GenJet_pt[Jet.Jet_genJetIdx[jets_indices[0]]])  << " \t " << (Jet.pt_JERdown[jets_indices[0]]-Jet.pt[jets_indices[0]])/(Jet.pt[jets_indices[0]]-GenJet_pt[jets_indices[0]]) << std::endl;
        
        if (maxBTagCMVA>-0.853213 && maxBTagCMVA <-0.853211) maxBTagCMVA = -0.999;
        if (maxSecondBTagCMVA>-0.853213 && maxSecondBTagCMVA <-0.853211) maxSecondBTagCMVA = -0.999;
        if (maxBTagCMVA_leading>-0.853213 && maxBTagCMVA_leading <-0.853211) maxBTagCMVA_leading = -0.999;
        Int_t EWKnsoft2_toSubtract = 0;
        Int_t EWKnsoft5_toSubtract = 0;
        Int_t EWKnsoft10_toSubtract = 0;
        float EWKHTsoft_toSubtract = 0.;

        Int_t Nsoft2 = 0;
        Int_t Nsoft5 = 0;
        Int_t Nsoft10 = 0;
        float HTsoft = 0;
        
        // //      SOFT EWK 3rd JET IN THE RAPIDITY GAP
//          for (int i=0;i<Jet.nsoftActivityEWKJets;i++){
//              TLorentzVector jet0;
//              bool found_softLeadingJet = false;
//              jet0.SetPtEtaPhiM(Jet.EWKsoft_pt[i],Jet.EWKsoft_eta[i],Jet.EWKsoft_phi[i],Jet.EWKsoft_mass[i]);
//              if (jet0.Eta() > max(jets_pv[0].Eta(), jets_pv[1].Eta()) || jet0.Eta() < min(jets_pv[0].Eta(), jets_pv[1].Eta()) || 
//                   jet0.DeltaR(Qjet1) < 0.4 || jet0.DeltaR(Qjet2) < 0.4 || jet0.DeltaR(lepton1) < 0.4 || jet0.DeltaR(lepton2) < 0.4 ) {
//                  if (Jet.EWKsoft_pt[i] > 1)  EWKHTsoft_toSubtract += Jet.EWKsoft_pt[i];
//                  if (Jet.EWKsoft_pt[i] > 2)  EWKnsoft2_toSubtract++;
//                  if (Jet.EWKsoft_pt[i] > 5)  EWKnsoft5_toSubtract++;
//                  if (Jet.EWKsoft_pt[i] > 10)  EWKnsoft10_toSubtract++;
//              }
//              else 
//                  if (!found_softLeadingJet) {
//                      Jet.softLeadingJet_pt = Jet.EWKsoft_pt[i];
//                      found_softLeadingJet = true;
//                  }
//          }

        
        
       
            
            
        
        
//         UNCOMMENT THIS -------------------------------------------------------------------- THIS part of code is already written below
//         Qjet1 = jets_pv[0];
//         Qjet2 = jets_pv[1];
//         
//         if (Qjet1.Pt() < Qjet2.Pt() ) {
//             TLorentzVector Qjet1_tmp = Qjet1;
//             Qjet1 = Qjet2;
//             Qjet2 = Qjet1_tmp;
//             int jet1Idx_tmp = jets_indices[0];
//             jets_indices[0] = jets_indices[1];
//             jets_indices[1] = jet1Idx_tmp;
//         }
//         
// 		if (good_jets>=3) {
// 			jet3_pt=jets_pv[2].Pt();
// 			jet3_eta=jets_pv[2].Eta();
// 		}
// 		qq=Qjet1+Qjet2;
// 		Float_t Mqq = qq.M();
// 		Float_t qq_pt = qq.Pt();
// 		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
// 		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));

//////////////////leptons////////////////
		TLorentzVector lepton1;
		TLorentzVector lepton2;
		TLorentzVector Zll;

            

//         if (abs(Mqq-qqMass_skim)>20) std::cout << evt << " \tjetIdx   "  << Jet.jetIdx1 << " \t " << Qjet1.Pt()  << " \t "  << Jet.jetIdx2 << " \t " << Qjet2.Pt() << " \t " << qqMass_skim  << " \t "  <<Mqq  << std::endl;  
//                   
//         hSelectionCuts->Fill(0., /*genweight*/1.*xsec[file_tag]/events_generated);

//         if(qqMass_skim > 400 && Mqq < 400)  std::cout << evt << " \tjetIdx   "  << Jet.jetIdx1 << " \t " << Qjet1.Pt()  << " \t "  << Jet.jetIdx2 << " \t " << Qjet2.Pt() << " \t " << qqMass_skim  << " \t "  <<Mqq  << std::endl;  
        
//         if (Mqq<180) continue;
//         if (Mqq<200) continue;
// //         if (Mqq<250) continue;      //----------------------------------------------------------------------------------------------------
//         if (!andrewSelection) if (Mqq<400) continue; 
//                 
        
//         hSelectionCuts->Fill(1, /*genweight*/1.*xsec[file_tag]/events_generated);


                
		int idx_1stLepton = 0;
		int idx_2ndLepton = 0;
		int count_l=0;


                if (region.CompareTo("mu")==0) {
		    count_l=0;
		    idx_1stLepton = 0;
		    idx_2ndLepton = 0;
		    for (int i=0; i<nselLeptons;i++ ){
                        if (fromNANO) {

                            if ( ( (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0) ) && abs((int) Muon_genPartFlav[i]) != 15 ) continue;
                            selLeptons_looseIdPOG[i] = 1;
                            if ( !(abs(selLeptons_dz[i]) < 0.2 && abs(selLeptons_dxy[i]) < 0.05) ) {
//                                 std::cout << "selLeptons_dz  " << selLeptons_dz[i] << "  \t  " << selLeptons_dxy[i] << std::endl;
                                continue;
                            }
                        }

                        if (!andrewSelection) if (!((selLeptons_looseIdPOG[i]>0) && (selLeptons_relIso04[i]<0.25) && (TMath::Abs(selLeptons_pdgId[i])==13 )) ) continue;
			    if (andrewSelection)  if (!((selLeptons_mediumIdPOG[i]>0) && (selLeptons_relIso04[i]<0.25) && (TMath::Abs(selLeptons_pdgId[i])==13 )) ) continue;
                            

                            if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;

                            if (count_l==1) {
				    idx_2ndLepton=i;
				    lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				    count_l++;
				    break;
			    }
// 			    if ( (count_l==0) && (selLeptons_relIso04[i]<0.15) && (selLeptons_tightIdPOG[i]>0)) {
			    if (count_l==0) {
                                if (andrewSelection) {if ((selLeptons_relIso04[i]>0.15) || (!selLeptons_tightIdPOG[i])) continue;}
				    idx_1stLepton=i;
				    lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				    count_l++;
			    }
		    }
		}



		
        
        if (fromNANO) {
            for (int i=0; i<nElectron;i++ ){
                if ( abs(selLeptons_dz[i]) < 0.2 && selLeptons_dxy[i] < 0.05 && Electron_pfRelIso03_all[i] < 0.25 && Electron_pt[i] > 15 ) nIsolatedElectrons++;
            }
        }
            
            
		
//        TLorentzVector Hll = lepton1+lepton2;
//        if  ( Hll.M()<100 )  continue;
//        else  ++count14;
//        std::cout << "Hll.M()   " << Hll.M()<< std::endl;



// 		lepton1.SetPtEtaPhiM(vLeptons_pt[0], vLeptons_eta[0], vLeptons_phi[0], vLeptons_mass[0]);	
// 		int idx_2ndLepton = 0;
// 		for (int i=1; i<nvLeptons;i++ ){
// 			if (vLeptons_charge[0]*vLeptons_charge[i] < 0) {
// 				idx_2ndLepton=i;
// 				break;
// 			}
// 		}
// 		lepton2.SetPtEtaPhiM(vLeptons_pt[idx_2ndLepton], vLeptons_eta[idx_2ndLepton], vLeptons_phi[idx_2ndLepton], vLeptons_mass[idx_2ndLepton]);
// 		float qter1 = 1.0;
// 		float qter2 = 1.0;
// 		float mu_correction1 = 1.0;
// 		float mu_correction2 = 1.0;
// 		if (data!=1) 	if (region.CompareTo("mu")==0) {
// 			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
// 			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
// 			}
// 		if (data==1) 	if (region.CompareTo("mu")==0){
// 			rmcor->momcor_data(lepton1, vLeptons_charge[0],  0, qter1);
// 			rmcor->momcor_data(lepton2, vLeptons_charge[idx_2ndLepton], 0, qter2);
// 			}

///////////////muon corrections 2016 calibration////////////////
// RoccoR  *rc = new RoccoR("/scratch/mandorli/Hmumu/CMSSW_8_0_25/src/giulioMandorli/2016/rcdata.2016.v3/");



        if (!fromNANO) {
            if (data!=1) {
		if (region.CompareTo("mu")==0 && !fromNANO) {
			double dataSF1, dataSF2 ;
			double mcSF1, mcSF2;
			if (data==1) {
				dataSF1 = rc->kScaleDT(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), 0, 0);
				dataSF2 = rc->kScaleDT(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*dataSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*dataSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
			if (data!=1) {
				double u1 = gRandom->Rndm();
				double u2 = gRandom->Rndm();
//                                 u1 = 0.5;
//                                 u2 = 0.5;
				mcSF1 = rc->kScaleAndSmearMC(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(),  selLeptons_trackerLayers[idx_1stLepton], u1, u2, 0, 0);
				mcSF2 = rc->kScaleAndSmearMC(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(),  selLeptons_trackerLayers[idx_2ndLepton], u1, u2, 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*mcSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*mcSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
		}
            }
        }
        

        //////////////////////////////////////////////////////////////// --------------- TriggerRequest
// 			if  (region.CompareTo("mu")==0) if (!((HLT_IsoMu24==1) || (HLT_IsoTkMu24==1)  )) continue; 
//             else  ++count4;
// 			if  (region.CompareTo("el")==0) if (!(HLT_Ele27_eta2p1 == 1)) continue;
                        
            
//         TFile* file_id_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunBCDEF.root");
//         TH2F* id_mu_bf = (TH2F*)file_id_mu_bf->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
//         TFile* file_id_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunGH.root");
//         TH2F* id_mu_aft = (TH2F*)file_id_mu_aft->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
// 
//         TFile* file_trig_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunBCDEF.root");
//         TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
//         TFile* file_trig_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunGH.root");
//         TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
// 
//         TFile* file_iso_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_LooseISO_LooseID_pt_eta_RunBCDEF.root");
//         TH2F* iso_mu_bf = (TH2F*)file_iso_mu_bf->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
//         TFile* file_iso_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_LooseISO_LooseID_pt_eta_RunGH.root");
//         TH2F* iso_mu_aft = (TH2F*)file_iso_mu_aft->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
// 
// 
//         TFile* file_tracker_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_ScaleFactor_tracker_80x.root");
//         TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");
//         TFile* file_trig_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_Tight27AfterIDISO.root");
//         TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_Tight27AfterIDISO");
//         TFile* file_id_el = TFile::Open("2016/triggerCorrection/TriggerEffMap_EIDISO_ZH.root");
//         TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_EIDISO_ZH");
// 
//         TFile* file_track_mu_bf = TFile::Open("2016/triggerCorrection/TriggerEffMap_Muons_trk_SF_RunBCDEF.root");
//         TH1F* track_mu_bf = (TH1F*)file_track_mu_bf->Get("TriggerEffMap_Graph");
//         TFile* file_track_mu_aft = TFile::Open("2016/triggerCorrection/TriggerEffMap_Muons_trk_SF_RunGH.root");
//         TH1F* track_mu_aft = (TH1F*)file_track_mu_aft->Get("TriggerEffMap_Graph");

        double TrkEff = 1.;
                if (data!=1) {
			if (region.CompareTo("mu")==0) {
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =20.1/36.4*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
	
				float eff1_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				float eff1_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ; 
				abs=0; 
				float eff1_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton1.Eta(), SF_mu_bf_err1,abs );  	
				float eff2_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton2.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton2.Eta(), SF_mu_bf_err1,abs );  	

				genweight*= eff1*eff1_id*eff2_id*eff1_iso*eff2_iso*eff1_tracker*eff2_tracker; 	
				TrkEff*= eff1*eff1_id*eff2_id*eff1_iso*eff2_iso*eff1_tracker*eff2_tracker; 	
			}

			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1 = getScaleFactor(trig_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs );  	
				float eff1_id =getScaleFactor(id_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_id =getScaleFactor(id_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				float eff1_tracker =getScaleFactor(tracker_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_tracker =getScaleFactor(tracker_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				genweight*=eff1* eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
// 				genweight0*=eff1* eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}

		}

		// 		cout << "  genweight 2    " << genweight  <<  endl; 
              
		string file_tag_str = file_tag.Data();  //ALREADY DECLARED: CHECK!
// 		if  (file_tag_str.find("DYJetstoLL")!=std::string::npos)  genweight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 


                
        double wwww = genweight0*PU*forwardJetCorrection*xsec[file_tag]/events_generated;
//         double wwww = genweight0*xsec[file_tag]/events_generated;
        
        float weightCorection = wwww;
        weightCorection = weightCorection - genweight0/events_generated*xsec[file_tag];
        
//         wwww = genweight0/events_generated;
//         if (forwardJetCorrection>0.1) wwww = wwww/forwardJetCorrection;
        hSelectionCuts->Fill(0., wwww);
        hSelectionCuts2->Fill(0., genweight0/events_generated*xsec[file_tag]);
        if (count_l>0) hSelectionCuts2->Fill(1., genweight0/events_generated*xsec[file_tag]);
        if (count_l>1) hSelectionCuts2->Fill(2., genweight0/events_generated*xsec[file_tag]);
        if (good_jets<2) continue;
//         if (count_l>1) hSelectionCuts2->Fill(3., genweight0/events_generated*xsec[file_tag]);
        else  ++count2;    
                
        

                
        Qjet1 = jets_pv[0];
        Qjet2 = jets_pv[1];
        

        
        
        if (Qjet1.Pt() < Qjet2.Pt() ) {
            TLorentzVector Qjet1_tmp = Qjet1;
            Qjet1 = Qjet2;
            Qjet2 = Qjet1_tmp;
            int jet1Idx_tmp = jets_indices[0];
            jets_indices[0] = jets_indices[1];
            jets_indices[1] = jet1Idx_tmp;
        }
        
		if (good_jets>=3) {
			jet3_pt=jets_pv[2].Pt();
			jet3_eta=jets_pv[2].Eta();
		}
		qq=Qjet1+Qjet2;
		Float_t Mqq = qq.M();
		Float_t qq_pt = qq.Pt();
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));
                
                
                
                
                
                
                Float_t Mqq_log = TMath::Log(Mqq);	

	//	if ((Mqq_log< 8.3227 )&&(Mqq_log>5.101)) if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0))  genweight*=func_Mqq->Eval(Mqq_log);		
		if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0))  genweight*=func_Mqq->Eval(Mqq_log);		



		float qgl_weight=1.;
		int apply_qgl=0;
		if (!( (data==1)|| (Jet.partonFlavour[jets_indices[0]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[0]])>=2) || (Jet.qgl[jets_indices[0]] < 0) ) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) < 4 ) qgl_weight=func_qgl_q->Eval(Jet.qgl[jets_indices[0]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) ==21 ) qgl_weight=func_qgl_g->Eval(Jet.qgl[jets_indices[0]]);
		}
		if (qgl_weight!=1.) apply_qgl+=1;
// // 	//	cout<<qgl_weight<<endl;
        genweight*=qgl_weight;     

        float qgl_weight_bothJet = qgl_weight;
        
// 		cout << "  genweight 3    " << genweight  <<  endl; 

                            

        
        
        qgl_weight=1.;
        if (!( (data==1)|| (Jet.partonFlavour[jets_indices[1]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[1]])>=2) || (Jet.qgl[jets_indices[1]] < 0)) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) < 4 ) qgl_weight=func_qgl_q->Eval(Jet.qgl[jets_indices[1]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) ==21 ) qgl_weight=func_qgl_g->Eval(Jet.qgl[jets_indices[1]]);
		}
		if (qgl_weight!=1.) apply_qgl+=1;
// 		if (data!=1) genweight*=qgl_norm[file_tag];

		
		genweight*=qgl_weight;
                qgl_weight_bothJet = qgl_weight_bothJet * qgl_weight;
                
                
                
//                 std::cout << "qgl_weight_bothJet  " << qgl_weight_bothJet << " \t "  <<  qgl_weight << " \t "  <<  std::endl;

                                        
// 	        cout  << "4 \t\t " << genweight << "  *  " << xsec[file_tag] << "/" << events_generated << endl;
// 		cout << "  genweight 4    " << genweight  <<  endl; 

//                 cout << "3   genweight     " << genweight <<  endl; 

		
//         std::cout << entry << "\t lepton1  " << lepton1.Pt() << " \t  " << lepton1.Phi() << " \t  " << lepton1.Eta() << " \t  " << lepton1.M() << " \t  "  << std::endl;
//         std::cout << entry << "\t lepton2  " << lepton2.Pt() << " \t  " << lepton2.Phi() << " \t  " << lepton2.Eta() << " \t  " << lepton2.M() << " \t  "  << std::endl;

                
//     ------------------------------------------- PROBLEM IN SOFT ACTIVITY -----------------------------------------------            
//             std::cout << evt << " \t ------------------------------------------------------------- "  << std::endl;   
                //      SOFT EWK 3rd JET IN THE RAPIDITY GAP
//          for (int i=0;i<Jet.nsoftActivityEWKJets;i++){
//              TLorentzVector jet0;
//              bool found_softLeadingJet = false;
//              jet0.SetPtEtaPhiM(Jet.EWKsoft_pt[i],Jet.EWKsoft_eta[i],Jet.EWKsoft_phi[i],0.);
// //              std::cout << jet0.Pt() << " \t " <<  jet0.Eta() << " \t " << jets_pv[0].Pt() << " \t " <<  jets_pv[0].Eta() << " \t " << jets_pv[1].Pt() << " \t " <<  jets_pv[1].Eta() << " \t " << lepton1.Pt() << " \t " <<  lepton1.Eta() << " \t " << lepton2.Pt() << " \t " <<  lepton2.Eta() <<   std::endl;
// //              jet0.SetPtEtaPhiM(Jet.EWKsoft_pt[i],Jet.EWKsoft_eta[i],Jet.EWKsoft_phi[i],Jet.EWKsoft_mass[i]);
//              if (jet0.Eta() > max(jets_pv[0].Eta(), jets_pv[1].Eta()) || jet0.Eta() < min(jets_pv[0].Eta(), jets_pv[1].Eta()) || 
//                   jet0.DeltaR(Qjet1) < 0.1 || jet0.DeltaR(Qjet2) < 0.1 || jet0.DeltaR(lepton1) < 0.1 || jet0.DeltaR(lepton2) < 0.1 ) {
// //                  std::cout << (jet0.Eta() >  max(jets_pv[0].Eta() ,jets_pv[1].Eta())) << " \t " << (jet0.Eta() < min(jets_pv[0].Eta(), jets_pv[1].Eta())) << " \t " <<  (jet0.DeltaR(Qjet1) < 0.3)  << " \t " << (jet0.DeltaR(Qjet2) < 0.3) << " \t " << (jet0.DeltaR(lepton1) < 0.3) << " \t " << (jet0.DeltaR(lepton2) < 0.3) << " \t " << Jet.EWKsoft_pt[i] <<   std::endl;
//                  if (Jet.EWKsoft_pt[i] > 1)  EWKHTsoft_toSubtract += Jet.EWKsoft_pt[i];
//                  if (Jet.EWKsoft_pt[i] > 2)  EWKnsoft2_toSubtract++;
//                  if (Jet.EWKsoft_pt[i] > 5)  EWKnsoft5_toSubtract++;
//                  if (Jet.EWKsoft_pt[i] > 10)  EWKnsoft10_toSubtract++;
//              }
//              else 
//                  if (!found_softLeadingJet) {
//                      Jet.softLeadingJet_pt = Jet.EWKsoft_pt[i];
//                      found_softLeadingJet = true;
//                  }
//          }
//          
// //          std::cout << entry << " \t " << EWKHTsoft_toSubtract << " \t " << Jet.EWKnsoft2 << " \t " << EWKnsoft2_toSubtract << " \t " << Jet.EWKnsoft5 << " \t "  <<  EWKnsoft5_toSubtract << " \t " << Jet.EWKnsoft10 << " \t " << EWKnsoft10_toSubtract  << std::endl;
// 
//          if (EWKnsoft2_toSubtract > Jet.EWKnsoft2) EWKnsoft2_toSubtract = Jet.EWKnsoft2;
//          if (EWKnsoft5_toSubtract > Jet.EWKnsoft5) EWKnsoft5_toSubtract = Jet.EWKnsoft5;
//          if (EWKnsoft10_toSubtract > Jet.EWKnsoft10) EWKnsoft10_toSubtract = Jet.EWKnsoft10;
         
//     ----------------------------------------- END PROBLEM IN SOFT ACTIVITY -----------------------------------------------  
         


                       //      SOFT EWK 3rd JET IN THE RAPIDITY GAP
         for (int i=0;i<Jet.nsoftActivityEWKJets;i++){
             TLorentzVector jet0;
             bool found_softLeadingJet = false;
             jet0.SetPtEtaPhiM(Jet.EWKsoft_pt[i],Jet.EWKsoft_eta[i],Jet.EWKsoft_phi[i],0.);
             if (jet0.Eta() < max(jets_pv[0].Eta(), jets_pv[1].Eta()) && jet0.Eta() > min(jets_pv[0].Eta(), jets_pv[1].Eta()) && 
                  jet0.DeltaR(Qjet1) > 0.2 && jet0.DeltaR(Qjet2) > 0.2 && jet0.DeltaR(lepton1) > 0.2 && jet0.DeltaR(lepton2) > 0.2 ) {

                 if (Jet.EWKsoft_pt[i] > 1)  HTsoft += Jet.EWKsoft_pt[i];
                 if (Jet.EWKsoft_pt[i] > 2)  Nsoft2++;
                 if (Jet.EWKsoft_pt[i] > 5)  Nsoft5++;
                 if (Jet.EWKsoft_pt[i] > 10)  Nsoft10++;
             }
             else 
                 if (!found_softLeadingJet) {
                     Jet.softLeadingJet_pt = Jet.EWKsoft_pt[i];
                     found_softLeadingJet = true;
                 }
         }
         
         
		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();
		
//         std::cout << "Zll_mass  " << Zll_mass << std::endl;


//             float Zll_mass_biasUp = 0;
//             float Zll_mass_biasDown = 0;
//             float Zll_mass_resolution = 0;
//                         
//             TLorentzVector lepton1_tmp;
//             TLorentzVector lepton2_tmp;  
//                   
//             lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*1.001,lepton1.Eta(),lepton1.Phi(),lepton1.M());
//             lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*1.001,lepton2.Eta(),lepton2.Phi(),lepton2.M());
//             Zll_mass_biasUp = (lepton1_tmp + lepton2_tmp).M();
//      
//                               
//             lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*0.999,lepton1.Eta(),lepton1.Phi(),lepton1.M());
//             lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*0.999,lepton2.Eta(),lepton2.Phi(),lepton2.M());
//             Zll_mass_biasDown = (lepton1_tmp + lepton2_tmp).M();
//             
//                         
//             lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*(1.+distribution(generator)),lepton1.Eta(),lepton1.Phi(),lepton1.M());
//             lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*(1.+distribution(generator)),lepton2.Eta(),lepton2.Phi(),lepton2.M());
//             Zll_mass_resolution = (lepton1_tmp + lepton2_tmp).M();
            
            
            
//         if (genweight<0.)  continue;
            
//             if(genweight > 0.) genweight = 1.;
//             else genweight = -1.;

        float weightHistoCut = 1;
        weightHistoCut = genweight;
        weightHistoCut = genweight/qgl_weight_bothJet/events_generated*xsec[file_tag]; 
        weightHistoCut = genweight0/events_generated*xsec[file_tag]; 
//         wwww = genweight0*PU*forwardJetCorrection*xsec[file_tag]/events_generated;
        weightHistoCut = wwww;       
        genweight/=events_generated/xsec[file_tag]; 
        genweight_noQGLcorrection = genweight/qgl_weight_bothJet;        
// genweight=1;

//         float weightCorection = weightHistoCut;
//         weightCorection = weightCorection - genweight0/events_generated*xsec[file_tag];
//         hSelectionCuts->Fill(0., weightCorection);
                
//         hSelectionCuts->Fill(0., weightHistoCut);

        // /*/
        if (Mqq<180) continue;
        if (Mqq<200) continue;
//         if (Mqq<250) continue;      //----------------------------------------------------------------------------------------------------
        if (!andrewSelection && Mqq<400) continue; 
//                 
        
        
        hSelectionCuts->Fill(1, weightHistoCut);
        

        
        
//         if (Qjet1.Pt() < 50) continue;
        if (!andrewSelection) if (Qjet1.Pt() < 35) continue;
        if (Qjet1.Pt() < 30) continue;
        if (Qjet1.Pt() < 20) continue;
        else  ++count5;
        hSelectionCuts->Fill(2, weightHistoCut);
        
        if (andrewSelection) if (Qjet2.Pt() < 30) continue;
//         if (Qjet2.Pt() < 30) continue;
        if (Qjet2.Pt() < 25) continue;
        if (Qjet2.Pt() < 20) continue;
//      if (Qjet2.Pt() < 22) continue;
        else  ++count6;
        hSelectionCuts->Fill(3, weightHistoCut);

        
//         if (count_l<2 && Zll_mass > 115 ) {
//         if (count_l<2) {
//             std::cout << entry << "\t lepton1  " << lepton1.Pt() << " \t  " << lepton1.Phi() << " \t  " << lepton1.Eta() << " \t  " << lepton1.M() << " \t  "  << std::endl;
//             std::cout << entry << "\t lepton2  " << lepton2.Pt() << " \t  " << lepton2.Phi() << " \t  " << lepton2.Eta() << " \t  " << lepton2.M() << " \t  "  << std::endl;
//         }
        
        weightHistoCut = weightHistoCut * TrkEff;
        if (count_l<2)  continue; 
        hSelectionCuts->Fill(4, weightHistoCut);


        
        
        if (Zll_mass < 110 ) continue;
        if (!andrewSelection) if (Zll_mass < 115 ) continue;      //---------------------------------------------------------------------------------------------------- comment for ControlRegion
        else  ++count12;
        hSelectionCuts->Fill(5, weightHistoCut);
                		
        if (Zll_mass > 155 ) continue;
//         if (Zll_mass > 150 ) continue;
//         if (!andrewSelection) if (Zll_mass > 145 ) continue;
        if (!andrewSelection) if (Zll_mass > 135 ) continue;      //---------------------------------------------------------------------------------------------------- comment for ControlRegion
        hSelectionCuts->Fill(6, weightHistoCut);
        
        
        
//         std::cout << "HERE 8   "  << std::endl;
        
        
//          std::cout << entry << "\t HLT_IsoMu24  " << HLT_IsoMu24 << " \t HLT_IsoTkMu24 " << HLT_IsoTkMu24  << std::endl;
                
        if  (region.CompareTo("mu")==0) if (!((HLT_IsoMu24==1) || (HLT_IsoTkMu24==1)  )) continue; 
        else  ++count4;
        hSelectionCuts->Fill(7, weightHistoCut);  //Trigger only      
                
                
        if (good_jets<2) continue;
        else  ++count2;    
        hSelectionCuts->Fill(8, weightHistoCut);
        
        
        if (qqDeltaEta < 1.) continue;
        if (!andrewSelection) if (qqDeltaEta < 2.5) continue;
        hSelectionCuts->Fill(9, weightHistoCut);
        
        
  //  /*/    
        if (count_l<2)  continue;
        else  ++count3;///*/ 
        hSelectionCuts->Fill(10, weightHistoCut);
        if ((selLeptons_charge[idx_1stLepton]*selLeptons_charge[idx_2ndLepton]) >0) continue;                                   // redundant
        hSelectionCuts->Fill(11, weightHistoCut);
        


        
		cut_flow[4]+=genweight;
                if (!andrewSelection) {
                    if (region.CompareTo("mu")==0) {
			if (lepton1.Pt()<30) continue;
                        else  ++count8;
                        hSelectionCuts->Fill(12, weightHistoCut);
			if (lepton2.Pt()<10) continue;
                        else  ++count9;
                        hSelectionCuts->Fill(13, weightHistoCut);
                    }
		}	
	 	if (region.CompareTo("el")==0) {
			if (lepton1.Pt()<30) continue;	
			if (lepton2.Pt()<20) continue;
		}	
		cut_flow[5]+=genweight;
	 	if (region.CompareTo("el")==0) {
			if (TMath::Abs(lepton1.Eta())>2.1) continue;	
			if (TMath::Abs(lepton2.Eta())>2.1) continue;
		}	
	 	if (!andrewSelection) {
                    if (region.CompareTo("mu")==0) {
			if (TMath::Abs(lepton1.Eta())>2.4) continue;
                        else  ++count10;
                        hSelectionCuts->Fill(14, weightHistoCut);
			if (TMath::Abs(lepton2.Eta())>2.4) continue;
                        else  ++count11;
                        hSelectionCuts->Fill(15, weightHistoCut);
                    }
		}


		
                

//                 if (maxBTagCMVA > 0.95) continue;
                if (andrewSelection) {if (maxBTagCSV > 0.8) continue; }
                else if (maxBTagCMVA > 0.8) continue;     //----------------------------------------------------------------------------------------------------
                    //----------------------------------------------------------------------------------------------------
//                 else  ++count13;
                hSelectionCuts->Fill(16, weightHistoCut);///*/



                
                
                presel+=genweight;

		presel_vtype[(int)(v_type+1)]+=genweight;

		
		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;              
		float Zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;        
		float ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity())/2. ;
		float zstar = Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() ) ;
                Zll_ystar = ystar;             
		Zll_zstar = zstar;        
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
                Float_t PhiZQ2 = TMath::Abs(Zll.DeltaPhi(Qjet2));
                Float_t EtaHQ1 = TMath::Abs(Zll.Eta() - Qjet1.Eta());
                Float_t EtaHQ2 = TMath::Abs(Zll.Eta() - Qjet2.Eta());
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t Rpt = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 


//ADDITIONAL CUTS
// 		if (!andrewSelection) if (abs(Zll_zstar) > 2.5 ) continue;
// // //                 if (qqDeltaEta>8) continue;       
//                 hSelectionCuts->Fill(17, weightHistoCut);
// //                 if (!andrewSelection) if (Rpt>0.4) continue;      //----------------------------------------------------------------------------------------------------
                
                weightHistoCut = weightHistoCut * qgl_weight_bothJet * 0.959045;
                
                hSelectionCuts->Fill(18, weightHistoCut);
                
// std::cout << "HERE 6   "<< std::endl; 
//                   cout << entry << "  \t  " <<evt << endl;
//                 cout << "\r" << entry << "  \t  " <<evt <<  "  \t 1 ";
//                  cout << "\n" << entry << "  \t  " <<evt <<  "  \t 1 ";
                  
//          if (abs((genweight-genweight0*forwardJetCorrection*PU*TrkEff*qgl_weight_bothJet*xsec[file_tag]/events_generated)/genweight)>0.1)  cout  << "DDD  " <<  genweight << "  \t  " <<genweight0*forwardJetCorrection*PU*TrkEff*qgl_weight_bothJet*xsec[file_tag]/events_generated << std::endl;
        
        
//                 float genJetMass = VBFFilter(genLeptons, genJet);
            float genJetMass_leading = 0;
            std::vector<TLorentzVector> genJetsWithoutLeptonsP4= VBFFilter(genLeptons, genJet);
            std::vector<TLorentzVector> recoLepton;
            recoLepton.push_back(lepton1);
            recoLepton.push_back(lepton2);
//             std::vector<TLorentzVector> genJetsWithoutLeptonsP4= VBFFilter(recoLepton, genJet);
            if (genJetsWithoutLeptonsP4.size() > 1) genJetMass_leading = (genJetsWithoutLeptonsP4[0] + genJetsWithoutLeptonsP4[1]).M();
            else if(file_tag.CompareTo("DYJetsToLL_M-105To160_VBFFilter-madgraphMLM")==0) continue;
            
//             if ((file_tag.CompareTo("DYJetsToLL_M-105To160-madgraphMLM")==0) && genJetMass_leading>350) {
// // //                 cout << "INSIDE" << file_tag.CompareTo("DYJetsToLL_M-105To160-madgraphMLM") << "  \t  " <<genJetMass_leading << "----------------------------------------------------------" << endl;
//                 continue;
//             }

//             if ((file_tag.CompareTo("DYJetsToLL_M-105To160-amcatnloFXFX")==0) && genJetMass_leading>350) continue;
//             if ((file_tag.CompareTo("DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX")==0) && genJetMass_leading>350) continue;
//             if ((file_tag.CompareTo("DYJetsToLL_M-105To160-amcatnloFXFX")==0) && genJetMass_leading<350) continue;
            if ((file_tag.CompareTo("DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX")==0) && genJetMass_leading<350) continue;
    
            for (int n=0; n < genJet.size(); n++) {
                if (genJetsWithoutLeptonsP4.size() > 0)  {if (genJet[n].DeltaR(genJetsWithoutLeptonsP4[0]) < 0.01) TMVA.genJetIdx0_genJetMassWithoutLepton = n;}
                if (genJetsWithoutLeptonsP4.size() > 1)  {if (genJet[n].DeltaR(genJetsWithoutLeptonsP4[1]) < 0.01) TMVA.genJetIdx1_genJetMassWithoutLepton = n;}
                if (genJetsWithoutLeptonsP4.size() > 2)  {if (genJet[n].DeltaR(genJetsWithoutLeptonsP4[2]) < 0.01) TMVA.genJetIdx2_genJetMassWithoutLepton = n;}
            }
            
            
            Double_t genLepton_deltaR1 = 10;
            Double_t genLepton_deltaR2 = 10;
            int genLepThatMatchToLeadingLepton = 0;
            int genLepThatMatchToSubleadingLepton = 0;
            Double_t tmp_R=0;
            for (int n=0; n<genLeptons.size(); n++) { 
                tmp_R = lepton1.DeltaR(genLeptons[n]);
                if (tmp_R < genLepton_deltaR1) {
                    genLepton_deltaR1 = tmp_R;
                    genLepThatMatchToLeadingLepton = n;
                }
            }
                
            for (int n=0; n<genLeptons.size(); n++) { 
                if (n!=genLepThatMatchToLeadingLepton ) {
                    genLepton_deltaR2 = min(genLepton_deltaR2, lepton2.DeltaR(genLeptons[n]));                
                    genLepThatMatchToSubleadingLepton = n;
                }
            }
            
            
//          READING GENJET_PT MATCHED WITH LEADING RECO JETS   
// std::cout << "HERE 7   "<< std::endl;
//             Double_t GenJet_deltaR1 = 10;
//             Double_t GenJet_deltaR2 = 10;
//             int genLepThatMatchToLeadingJet = 0;
            float jet1_deltaGenPT = 1;
            float jet2_deltaGenPT = 1;
             if ((data!=1) && (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0))  {
//                 if (file_tag.CompareTo("EKW_LLJJ_pythia8")==0) {
//                     jet1_deltaGenPT = 0;
//                     jet2_deltaGenPT = 0;
//                 }
//                 else {
                TMVA.genJetPt_match1 = genJet[jets_indices[0]].Pt();
                TMVA.genJetPt_match2 = genJet[jets_indices[1]].Pt();
                jet1_deltaGenPT = (Qjet1.Pt() - TMVA.genJetPt_match1) / Qjet1.Pt();
                jet2_deltaGenPT = (Qjet2.Pt() - TMVA.genJetPt_match2) / Qjet2.Pt();
//                 }
             }
//              std::cout << "HERE 8   " << genLeptons.size()<< " \t " << genLepThatMatchToLeadingLepton << " \t " << genLepThatMatchToSubleadingLepton << std::endl;
//             Double_t tmp_jetR=0;
//             for (int n=0; n<genJet.size(); n++) { 
//                 tmp_jetR = Qjet1.DeltaR(genJet[n]);
//                 if (tmp_jetR < GenJet_deltaR1) {
//                     GenJet_deltaR1 = tmp_jetR;
//                     genLepThatMatchToLeadingJet = n;
//                     if (GenJet_deltaR1<0.4) jet1_deltaGenPT = (Qjet1.Pt() - genJet[n].Pt()) / Qjet1.Pt();
//                 }
//             }
//                 
//             for (int n=0; n<genJet.size(); n++) { 
//                 if (n!=genLepThatMatchToLeadingJet ) {
//                     tmp_jetR = Qjet2.DeltaR(genJet[n]);
//                     if (tmp_jetR < GenJet_deltaR2) {
//                         GenJet_deltaR2 = tmp_jetR;                
//                         if (GenJet_deltaR2<0.4) jet2_deltaGenPT = (Qjet2.Pt() - genJet[n].Pt()) / Qjet2.Pt();
//                     }
//                 }
//             }
            
            
            
            float Zll_mass_biasUp = 0;
            float Zll_mass_biasDown = 0;
            float Zll_mass_resolutionUp = 0;
            float Zll_mass_resolutionDown = 0;
                        
            TLorentzVector lepton1_tmp;
            TLorentzVector lepton2_tmp;  
             
            lepton1_tmp.SetPtEtaPhiM(lepton1.Pt(),lepton1.Eta(),lepton1.Phi(),lepton1.M());
            lepton2_tmp.SetPtEtaPhiM(lepton2.Pt(),lepton2.Eta(),lepton2.Phi(),lepton2.M());
                
            if ((data!=1) /*&& (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0)*/) {
//                 lepton1_tmp.SetPtEtaPhiM((lepton1.Pt()-genLeptons[genLepThatMatchToLeadingLepton].Pt())*1.1+genLeptons[genLepThatMatchToLeadingLepton].Pt(),lepton1.Eta(),lepton1.Phi(),lepton1.M());
//                 lepton2_tmp.SetPtEtaPhiM((lepton2.Pt()-genLeptons[genLepThatMatchToSubleadingLepton].Pt())*1.1+genLeptons[genLepThatMatchToSubleadingLepton].Pt(),lepton2.Eta(),lepton2.Phi(),lepton2.M());
                lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*1.002,lepton1.Eta(),lepton1.Phi(),lepton1.M());
                lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*1.002,lepton2.Eta(),lepton2.Phi(),lepton2.M());
            }
            Zll_mass_biasUp = (lepton1_tmp + lepton2_tmp).M();
            

            if ((data!=1) /*&& (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0)*/) {                  
//                 lepton1_tmp.SetPtEtaPhiM((lepton1.Pt()-genLeptons[genLepThatMatchToLeadingLepton].Pt())*0.9+genLeptons[genLepThatMatchToLeadingLepton].Pt(),lepton1.Eta(),lepton1.Phi(),lepton1.M());
//                 lepton2_tmp.SetPtEtaPhiM((lepton2.Pt()-genLeptons[genLepThatMatchToSubleadingLepton].Pt())*0.9+genLeptons[genLepThatMatchToSubleadingLepton].Pt(),lepton2.Eta(),lepton2.Phi(),lepton2.M());
                lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*0.998,lepton1.Eta(),lepton1.Phi(),lepton1.M());
                lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*0.998,lepton2.Eta(),lepton2.Phi(),lepton2.M());
            }
            Zll_mass_biasDown = (lepton1_tmp + lepton2_tmp).M();
                        
//             if ((data!=1) && (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0)) {
//                 lepton1_tmp.SetPtEtaPhiM(lepton1.Pt()*(1.+distribution(generator)),lepton1.Eta(),lepton1.Phi(),lepton1.M());
//                 lepton2_tmp.SetPtEtaPhiM(lepton2.Pt()*(1.+distribution(generator)),lepton2.Eta(),lepton2.Phi(),lepton2.M());
//             }
                        
            
            if ((data!=1) && (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0) && (genLeptons.size()>1)) {                  
                lepton1_tmp.SetPtEtaPhiM((lepton1.Pt()-genLeptons[genLepThatMatchToLeadingLepton].Pt())*1.1+genLeptons[genLepThatMatchToLeadingLepton].Pt(),lepton1.Eta(),lepton1.Phi(),lepton1.M());
                lepton2_tmp.SetPtEtaPhiM((lepton2.Pt()-genLeptons[genLepThatMatchToSubleadingLepton].Pt())*1.1+genLeptons[genLepThatMatchToSubleadingLepton].Pt(),lepton2.Eta(),lepton2.Phi(),lepton2.M());
            }
            Zll_mass_resolutionUp = (lepton1_tmp + lepton2_tmp).M();
            
            if ((data!=1) && (file_tag.CompareTo("EKW_LLJJ_pythia8")!=0) && (genLeptons.size()>1)) {                  
                lepton1_tmp.SetPtEtaPhiM((lepton1.Pt()-genLeptons[genLepThatMatchToLeadingLepton].Pt())*0.9+genLeptons[genLepThatMatchToLeadingLepton].Pt(),lepton1.Eta(),lepton1.Phi(),lepton1.M());
                lepton2_tmp.SetPtEtaPhiM((lepton2.Pt()-genLeptons[genLepThatMatchToSubleadingLepton].Pt())*0.9+genLeptons[genLepThatMatchToSubleadingLepton].Pt(),lepton2.Eta(),lepton2.Phi(),lepton2.M());
            }
            Zll_mass_resolutionDown = (lepton1_tmp + lepton2_tmp).M();

            
            
//             std::cout << "HERE 7   "<< std::endl; 
            
            ngenJetsWithoutLeptonsP4 = genJetsWithoutLeptonsP4.size();
            for (int n=0; n<ngenJetsWithoutLeptonsP4;n++){ 
                genJetsWithoutLeptonsP4_pt.push_back(genJetsWithoutLeptonsP4[n].Pt());
                genJetsWithoutLeptonsP4_phi.push_back(genJetsWithoutLeptonsP4[n].Phi());
                genJetsWithoutLeptonsP4_eta.push_back(genJetsWithoutLeptonsP4[n].Eta());
                genJetsWithoutLeptonsP4_mass.push_back(genJetsWithoutLeptonsP4[n].M());                
            }
            
            
            float minAbsEta_GEN = 10;
            if (genJetsWithoutLeptonsP4.size() > 1) minAbsEta_GEN = min(abs(genJetsWithoutLeptonsP4[0].Eta()), abs(genJetsWithoutLeptonsP4[1].Eta()));
            
            float maxAbsEta = 0;
            float minAbsEta = 10;
            minAbsEta = min(abs(Qjet1.Eta()), abs(Qjet2.Eta()));
            maxAbsEta = max(abs(Qjet1.Eta()), abs(Qjet2.Eta()));
            TMVA.JetPuId_maxAbsEta = Jet.puId[jets_indices[0]];
            if (abs(Qjet1.Eta()) < abs(Qjet2.Eta()))  Jet.puId[jets_indices[2]];
            
            TMVA.Jet1q_puId = Jet.puId[jets_indices[0]];           
            TMVA.Jet2q_puId = Jet.puId[jets_indices[1]];
                        
                        
//             if (genJetMass_leading > 500 ) continue;
//             std::cout << std::endl << Qjet1.Pt() << "  " << Qjet1.Eta() << " \t " << Qjet2.Pt() << "  " << Qjet2.Eta() << " \t\t " << genJetsWithoutLeptonsP4[0].Pt() << "  " << genJetsWithoutLeptonsP4[0].Eta() << " \t " << genJetsWithoutLeptonsP4[1].Pt() << "  " << genJetsWithoutLeptonsP4[1].Eta() << " \t\t " << Mqq << "  " << genJetMass_leading << std::endl << std::endl;
            

            float deltaR1 = 100;
            float deltaR2 = 100;
            float deltaR_tmp = 100;
            int indexFisrtJet = 0;
            int indexSecondJet = 0;
            
            for (int n=0; n < genJetsWithoutLeptonsP4.size(); n++) {
                
                deltaR_tmp = genJetsWithoutLeptonsP4[n].DeltaR(Qjet1);
                if (deltaR1 > deltaR_tmp) {
                    deltaR1 = deltaR_tmp;
                    indexFisrtJet = n;
                }
}

deltaR_tmp = 100;
           for (int n=0; n < genJetsWithoutLeptonsP4.size(); n++) {
                if (indexFisrtJet==n) continue;
                deltaR_tmp = genJetsWithoutLeptonsP4[n].DeltaR(Qjet2);
                if (deltaR2 > deltaR_tmp) {
                    deltaR2 = deltaR_tmp;
                    indexSecondJet = n;
                }
            } 
            
            float genJetMass_matched = 0;
            if (genJetsWithoutLeptonsP4.size() > 1 && deltaR1<0.3 && deltaR2<0.3) genJetMass_matched = (genJetsWithoutLeptonsP4[indexFisrtJet] + genJetsWithoutLeptonsP4[indexSecondJet]).M();
//             if (Mqq - genJetMass_leading < 400) continue;
            

            
//         TVector3 BoostVector_toMuonRestFrame = Zll.BoostVector();
//         TLorentzVector positiveMuon_newSys;
//         TLorentzVector positiveMuon_OldSys;
//         if (selLeptons_charge[idx_1stLepton] > 0) positiveMuon_newSys = lepton1;
//         else positiveMuon_newSys = lepton2;
//         positiveMuon_OldSys = positiveMuon_newSys;
//         positiveMuon_newSys.Boost(-BoostVector_toMuonRestFrame);
//         TVector3 H_direction  = Zll.Vect();
//         TVector3 mu_direction = positiveMuon_newSys.Vect();
//         cosThetaStar = H_direction.Dot(mu_direction)/H_direction.Mag()/mu_direction.Mag();
                
        if (selLeptons_charge[idx_1stLepton] > 0) cosThetaStar = computeThetaStar (lepton1, lepton2) ;
        else cosThetaStar = computeThetaStar (lepton2, lepton1)   ;   
                

        cosThetaStarJet = computeThetaStar (Qjet1, Qjet2) ;
        cosThetaStarJet = cosThetaStarJet;
        absCosThetaStarJet = abs(cosThetaStarJet);

        
        TVector3 mu_1 = lepton1.Vect();
        TVector3 mu_2 = lepton2.Vect();
        TVector3 j_1 = Qjet1.Vect();
        TVector3 j_2 = Qjet2.Vect();


        cosThetaPlane = (mu_1.Cross(mu_2)).Dot(j_1.Cross(j_2))/mu_1.Mag()/mu_2.Mag()/j_1.Mag()/j_2.Mag();
        cosThetaPlane = abs(cosThetaPlane);
        
        float Inv_mass=0;
        float mumujj_pt=0;
        float Invariant_MassLog=0;
        float impulsoZ=0;
        float energytot=0;
        float energytotLog=0;
        
        TLorentzVector Sum;

        TLorentzVector W1_p4,W2_p4;
        Sum= lepton1+lepton2+Qjet1+Qjet2;  //Sum=Higgs+Qjet1+Qjet2
        mumujj_pt=Sum.Pt();
        Inv_mass=Sum.M();
        Invariant_MassLog=TMath::Log(Inv_mass);
        impulsoZ=Sum.Pz();
        energytot=Sum.E();
        energytotLog=TMath::Log(Sum.E());
        
        TLorentzVector Hll;
        Hll.SetPtEtaPhiM(Zll_pt,Zll_eta,Zll_phi,Zll_mass);
        
        int M_ll_cutbin=1;
        if(gen_mass<105) M_ll_cutbin=0;
        if(gen_mass>160) M_ll_cutbin=2;

        TVector3 hll3D= Hll.Vect();
        
        float theta1=0; 
        float theta2=0;
        float X_parton1=0;
        float X_parton2=0;
        float W_mass_virtual1=0;
        float W_mass_virtual2=0;
        float W_Pt_virtual1=0;
        float W_Pt_virtual2=0;
        float W_eta_virtual1=0;
        float W_eta_virtual2=0;
        float W_phi_virtual1=0;
        float W_phi_virtual2=0;
        float Parton_mass1=0;
        float Parton_mass2=0;
        float thetastarW1=0;
        float thetastarW2=0;
        float thetastarW2toHW1=0;
        float thetastarW1toHW2=0;
        float thetastarHtoWW=0;
        float WWmass=0;

       //// theta=1;
        theta1 = hll3D.Dot(j_1)/hll3D.Mag()/j_1.Mag();
        theta2 = hll3D.Dot(j_2)/hll3D.Mag()/j_2.Mag();

        
        X_parton1=0.5*(Sum.E()+Sum.Pz())/6500;
        X_parton2=0.5*(Sum.E()-Sum.Pz())/6500;
        
        TLorentzVector parton1_p4;
        TLorentzVector parton2_p4;
        
        parton1_p4.SetPxPyPzE(0, 0, X_parton1*6500, X_parton1*6500);
        parton2_p4.SetPxPyPzE(0, 0, -X_parton2*6500, X_parton2*6500);
        
//         if(Qjet1.Eta()>Qjet2.Eta()){
//         W1_p4= parton1_p4-Qjet1;
//         W2_p4=parton2_p4-Qjet2;
//         }else{
//         W1_p4=parton1_p4-Qjet2;
//         W2_p4=parton2_p4-Qjet1;
//         }
        
///////// Ordinati nel pt del jet che li produce ////////////////
        if(Qjet1.Eta()>Qjet2.Eta()){
        W1_p4= parton1_p4-Qjet1;
        W2_p4=parton2_p4-Qjet2;
        }else{
        W1_p4=parton2_p4-Qjet1;
        W2_p4=parton1_p4-Qjet2;
        }
/////////////////////////////////////////////////////////////////

        WWmass=(W1_p4+W2_p4).M();
                
        thetastarW1=computeThetaStar(Qjet1,parton1_p4);
        thetastarW2=computeThetaStar(Qjet2,parton2_p4);
        thetastarW2toHW1= computeThetaStar(W1_p4, Hll-W1_p4);
        thetastarW1toHW2=computeThetaStar(W2_p4, Hll- W2_p4);
        thetastarHtoWW=computeThetaStar(W1_p4,W2_p4);
        //if(!(thetastarHtoWW>-1&&thetastarHtoWW<1)) cout<< thetastarHtoWW <<" "<< (W1_p4+W2_p4).M()<< endl;
        
        W_mass_virtual1=W1_p4.M();
        W_mass_virtual2=W2_p4.M();
        
        W_Pt_virtual1=W1_p4.Pt();
        W_Pt_virtual2=W2_p4.Pt();
        W_eta_virtual1=W1_p4.Eta();
        W_eta_virtual2=W2_p4.Eta();
        W_phi_virtual1=W1_p4.Phi();
        W_phi_virtual2=W2_p4.Phi();
        
        Parton_mass1=parton1_p4.M();
        Parton_mass2=parton2_p4.M();
        
        Double_t MuonPtErr1=0;
        Double_t MuonPtErr2=0;
        MuonPtErr1=hmuon->GetBinContent(hmuon->FindBin(log(lepton1.Pt()),abs(lepton1.Eta())));
        MuonPtErr2=hmuon->GetBinContent(hmuon->FindBin(log(lepton2.Pt()),abs(lepton2.Eta())));
       
       Float_t deltaMRel=0;
       Float_t deltaM=0;  
       Float_t normalizedDistance_from_mH=0;
       deltaMRel=TMath::Sqrt(0.5*(pow(MuonPtErr1,2)+pow(MuonPtErr2,2)));
       deltaM=Zll_mass*deltaMRel;
       if (deltaM>0.00001) normalizedDistance_from_mH = abs(Zll_mass - 125)/deltaM;
       else normalizedDistance_from_mH = 5.;
       

       
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// END PRESELECTION/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			counter++;
//        float TMVA.weightMVA = genweight;
                        
//             if(genweight > 0.) genweight = 1.;
//             else genweight = -1.;

            
// 		genweight/=events_generated/xsec[file_tag]; //  <------spostato prima della preselection

//         cout << file_tag << " \t\t " << xsec[file_tag] << "/" << events_generated << endl;
// 	        cout  << "2 \t\t " << genweight << "  *  " << xsec[file_tag] << "/" << events_generated << endl;

                        
            float atanhCMVA = atanh(maxBTagCMVA);
//            atanhCMVA = maxBTagCMVA;
//            if (atanhCMVA < -5.)  atanhCMVA = -5;
           

            
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// LEONARDO VARIABLE /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    
            float HiggsSister1_Leading_eta = 0;
            float HiggsSister2_Subleading_eta = 0;
            float HiggsSister1_Leading_phi = 0;
            float HiggsSister2_Subleading_phi = 0;
            float HiggsSister1_Leading_R = 0;
            float HiggsSister2_Subleading_R = 0;


            
            
            HiggsSister1_Leading_eta    = abs(jets_pv[0].Eta() - GenHiggsSisters.eta[0]);
            HiggsSister2_Subleading_eta = abs(jets_pv[1].Eta() - GenHiggsSisters.eta[0]);
            
            HiggsSister1_Leading_phi    = abs(jets_pv[0].Phi() - GenHiggsSisters.phi[0]);
            HiggsSister2_Subleading_phi = abs(jets_pv[1].Phi() - GenHiggsSisters.phi[0]);

            HiggsSister1_Leading_R      = sqrt(HiggsSister1_Leading_eta*HiggsSister1_Leading_eta + HiggsSister1_Leading_phi*HiggsSister1_Leading_phi);
            HiggsSister2_Subleading_R   = sqrt(HiggsSister2_Subleading_eta*HiggsSister2_Subleading_eta + HiggsSister2_Subleading_phi*HiggsSister2_Subleading_phi);

            if (HiggsSister1_Leading_R > HiggsSister2_Subleading_R) {   // Leading_jet match higgsSister[0]
                HiggsSister2_Subleading_eta = abs(jets_pv[1].Eta() - GenHiggsSisters.eta[1]);
                HiggsSister2_Subleading_phi = abs(jets_pv[1].Phi() - GenHiggsSisters.phi[1]);
                HiggsSister2_Subleading_R   = sqrt(HiggsSister2_Subleading_eta*HiggsSister2_Subleading_eta + HiggsSister2_Subleading_phi*HiggsSister2_Subleading_phi);
            }
            else {  // Leading_jet match higgsSister[1]
                HiggsSister1_Leading_eta    = abs(jets_pv[0].Eta() - GenHiggsSisters.eta[1]);
                HiggsSister1_Leading_phi    = abs(jets_pv[0].Phi() - GenHiggsSisters.phi[1]);
                HiggsSister1_Leading_R      = sqrt(HiggsSister1_Leading_eta*HiggsSister1_Leading_eta + HiggsSister1_Leading_phi*HiggsSister1_Leading_phi);
            }
            
            

            hHiggsSister1_Leading_eta->Fill(HiggsSister1_Leading_eta, genweight);
            hHiggsSister2_Subleading_eta->Fill(HiggsSister2_Subleading_eta, genweight);
            
            hHiggsSister1_Leading_phi->Fill(HiggsSister1_Leading_phi, genweight);
            hHiggsSister2_Subleading_phi->Fill(HiggsSister2_Subleading_phi, genweight);
            


            TLorentzVector  HiggsSister1;
            TLorentzVector HiggsSister2;
            HiggsSister1.SetPtEtaPhiM(GenPart_pt[3], GenPart_eta[3], GenPart_phi[3], GenPart_mass[3]);
            HiggsSister2.SetPtEtaPhiM(GenPart_pt[4], GenPart_eta[4], GenPart_phi[4], GenPart_mass[4]);
            

            
            float r1 = 10;
            float r2 = 10;
            int indexJetHiggsSister1 = -1;
            int indexJetHiggsSister2 = -1;
            float minDeltaR1 = 15;
            float minDeltaR2 = 15;
            if (file_tag.CompareTo("VBF_HToMuMu")==0) {
            for (int i=0;i<good_jets;i++) {
//                 if(jets_pv[i].Pt()<25) continue;
                
                r1 = HiggsSister1.DeltaR(jets_pv[i]);
                r2 = HiggsSister2.DeltaR(jets_pv[i]);

                if (minDeltaR1>r1) {
                    minDeltaR1=r1;
                    indexJetHiggsSister1 = i;
                }
                if (minDeltaR2>r2) {
                    minDeltaR2=r2;
                    indexJetHiggsSister2 = i;
                }
            }
            
            
            
                       
//             for (int i=0;i<good_jets;i++) {
//                 if(jets_pv[i].Pt()<25) continue;
//                 
//                 r2 = HiggsSister2.DeltaR(jets_pv[i]);
// 
// 
//                 if (minDeltaR2>r2 && i!=indexJetHiggsSister1) {
//                     minDeltaR2=r2;
//                     indexJetHiggsSister2 = i;
//                 }
//             }
            
            
            
            

            if (indexJetHiggsSister1>indexJetHiggsSister2) {
                float rtmp = indexJetHiggsSister1;
                indexJetHiggsSister1 = indexJetHiggsSister2;
                indexJetHiggsSister2 = rtmp;
            }
                
            float minAngle= 0.99;
            hHiggsSister1_Leading_R->Fill(min(minDeltaR1,minAngle), genweight);
            hHiggsSister2_Subleading_R->Fill(min(minDeltaR2,minAngle), genweight);
            if (indexJetHiggsSister1 == indexJetHiggsSister2)  hHiggsSister1_HiggsSister2_R->Fill(HiggsSister1.DeltaR(HiggsSister2), genweight);

            
            TMVA.countJet25 = 0;
            for (int i=0;i<good_jets;i++) {
//                 std::cout << good_jets << " \t " << jets_pv[i].Pt() << std::endl;
                if (jets_pv[i].Pt() >= 25 ) TMVA.countJet25++;
            }

            if (TMVA.countJet25>1) hSelectionCuts2->Fill(3., genweight0/events_generated*xsec[file_tag]);  
            
            JET_MVA_NjetsEtaMinus = 0;
            JET_MVA_NjetsEtaPlus = 0;
            for (int i=0;i<good_jets;i++) {
                if(jets_pv[i].Pt()<25) continue;
                if (Zll_eta>jets_pv[i].Eta()) JET_MVA_NjetsEtaMinus++;
                else JET_MVA_NjetsEtaPlus++;                
            }
                   
                                
            std::vector<float> jet_MVA_value;
            fileMVA.cd();
                for (int i=0;i<good_jets;i++) {
                    
                    if(jets_pv[i].Pt()<25) {jet_MVA_value.push_back(-1); continue;}
                        
                    if ((indexJetHiggsSister1 == i) || (indexJetHiggsSister2 == i)) JET_MVA_fromVBF = 1;
                    else JET_MVA_fromVBF = 0;
                    JET_MVA_index = i;
                    JET_MVA_pt = jets_pv[i].Pt();
                    JET_MVA_p = jets_pv[i].P();
                    JET_MVA_eta = jets_pv[i].Eta();
                    JET_MVA_phi = jets_pv[i].Phi();
                    JET_MVA_mass = jets_pv[i].M();
                    JET_MVA_id = Jet.id[jets_indices[i]];
                    JET_MVA_puId = Jet.puId[jets_indices[i]];
                    JET_MVA_Njets = TMVA.countJet25;
                    JET_MVA_qgl = Jet.qgl[jets_indices[i]];
                    JET_MVA_area = Jet.axis2[jets_indices[i]];
                    JET_MVA_mult = Jet.mult[jets_indices[i]];
                    
                    JET_MVA_pt1 = 0;
                    JET_MVA_phi1 = 0;
                    JET_MVA_eta1 = 0;
                    JET_MVA_pt2 = 0;
                    JET_MVA_eta2 = 0;
                    JET_MVA_phi2 = 0;
                    JET_MVA_pt3 = 0;
                    JET_MVA_eta3 = 0;
                    JET_MVA_phi3 = 0;
                    int idxToWrite = 1;
                    for (int j=0;j<good_jets;j++) {
                        if ((jets_pv[j].Pt()<25) || (j==i)) continue;
                        
                        if (idxToWrite==1) {JET_MVA_pt1 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta1 = jets_pv[j].Eta(); JET_MVA_phi1 = jets_pv[j].Phi();}
                        if (idxToWrite==2) {JET_MVA_pt2 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta2 = jets_pv[j].Eta(); JET_MVA_phi2 = jets_pv[j].Phi();}
                        if (idxToWrite==3) {JET_MVA_pt3 = jets_pv[j].Pt()/JET_MVA_pt; JET_MVA_eta3 = jets_pv[j].Eta(); JET_MVA_phi3 = jets_pv[j].Phi();}
                        idxToWrite++;
                    }

        
                    JET_MVA_higgh_pt = Zll_pt;
                    JET_MVA_higgh_eta = Zll_eta;
                    JET_MVA_higgh_phi = Zll_pt;
                    JET_MVA_higgh_mass = Zll_mass;
                    
                    
                    JET_MVA_etaIdx = -JET_MVA_NjetsEtaMinus;
                    for (int j=0;j<good_jets;j++)        
                        if (jets_pv[j].Pt()>25 && JET_MVA_eta>jets_pv[j].Eta() && i!=j) JET_MVA_etaIdx++;
        
//                     if (JET_MVA_eta>Zll_eta) JET_MVA_etaIdx++;
                    if (JET_MVA_etaIdx>=0) JET_MVA_etaIdx++;
        
        
//                     if (TMVA.countJet25>2) std::cout  << i << "    " << Zll_eta << " \t "  << JET_MVA_Njets - JET_MVA_NjetsEtaMinus - JET_MVA_NjetsEtaPlus << " \t "  << JET_MVA_NjetsEtaMinus << " \t "  << JET_MVA_NjetsEtaPlus << " \t "  << JET_MVA_eta << " \t "  << JET_MVA_eta1 << " \t "  << JET_MVA_eta2 << " \t "  << JET_MVA_eta3 << " \t "  << JET_MVA_etaIdx << " \t "  << std::endl;

                    
                    JET_MVA_index_float = JET_MVA_index;
                    JET_MVA_puId_float = JET_MVA_puId;
                    JET_MVA_Njets_float = JET_MVA_Njets;
                    JET_MVA_etaIdx_float = JET_MVA_etaIdx;
                    JET_MVA_NjetsEtaMinus_float = JET_MVA_NjetsEtaMinus;
                    JET_MVA_NjetsEtaPlus_float = JET_MVA_NjetsEtaPlus;
        
                    jet_MVA_value.push_back(readerJetMVA->EvaluateMVA("BDTG"));
//                     if (TMVA.countJet25>2) std::cout << "AppenaValutato " << i << " \t " << readerJetMVA->EvaluateMVA("BDTG") << std::endl;
//                     jet_MVA_value.push_back(0);
                
//                     if (TMVA.countJet25>2) treeMVAjet->Fill(); 
                    if (TMVA.countJet25==2) treeMVAjet->Fill(); 
                    
//                     if (TMVA.countJet25==2 && JET_MVA_fromVBF==0) std::cout << entry << " \t " << i << " \t " << indexJetHiggsSister2 << " \t " << indexJetHiggsSister1 << " \t " << JET_MVA_Njets << " \t " << std::endl; 
                        
                        
            }
        

//         for (int i=0;i<good_jets;i++)std::cout << "Qui fuori " << i << " \t " << jet_MVA_value[i] << std::endl;
            
        int higherJetMVA = 0;
        int secondJetMVA = 1;
        float higherMVA = -0.99;
        for (int i=0;i<good_jets;i++) {
            if (jet_MVA_value[i] > higherMVA) {
                higherMVA = jet_MVA_value[i];
                higherJetMVA = i;
            }
        }
        
        higherMVA = -0.99;
        for (int i=0;i<good_jets;i++) {
            if (jet_MVA_value[i] > higherMVA && i!=higherJetMVA) {
                higherMVA = jet_MVA_value[i];
                secondJetMVA = i;
            }
        }
        
        float x_indices = ((secondJetMVA==indexJetHiggsSister2 && higherJetMVA==indexJetHiggsSister1) || (secondJetMVA==indexJetHiggsSister1 && higherJetMVA==indexJetHiggsSister2)) + 0.5;
        float y_indices = (1==indexJetHiggsSister2 && 0==indexJetHiggsSister1) + 0.5;
        
        histo_TrueIndices_MVAIndices->Fill(x_indices, y_indices, genweight);
        if (TMVA.countJet25>2) histo_TrueIndices_MVAIndices_3Jets->Fill(x_indices, y_indices, genweight);
        
//           if (TMVA.countJet25>2 && jet_MVA_value[2] > jet_MVA_value[2])
//             if (TMVA.countJet25>2) std::cout  << higherJetMVA << " \t "  << secondJetMVA << " \t "  << indexJetHiggsSister1 << " \t "  << indexJetHiggsSister2 << " \t "  << x_indices << " \t "  << y_indices << " \t "  << std::endl;
//             std::cout << higherJetMVA << " \t " << secondJetMVA << std::endl;
//             if (TMVA.countJet25>2) std::cout  << " ------------------------------------------------------- "  << std::endl;
            
            
            }
            
            
            
            
            
        fileMVA.cd();
        
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// END LEONARDO VARIABLE /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            
////////////////////////////TO NORMALIZE FILTER105////////////////////////////////////////////////
//             float ratioForComparisonIn_lheNj_Mll_eventNumber = 1.;
//            if(file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0){
//                 if(lheNj==0)  ratioForComparisonIn_lheNj_Mll_eventNumber = xsec["DYJetstoLL_amc_0J"]/xsec["DYJetstoLL_amc_Filter105"];
//                 if(lheNj==1)  ratioForComparisonIn_lheNj_Mll_eventNumber = xsec["DYJetstoLL_amc_1J"]/xsec["DYJetstoLL_amc_Filter105"];//genweight=0;
//                 if(lheNj==2)  ratioForComparisonIn_lheNj_Mll_eventNumber = xsec["DYJetstoLL_amc_2J"]/xsec["DYJetstoLL_amc_Filter105"];
//                 genweight=genweight*ratioForComparisonIn_lheNj_Mll_eventNumber;
//                 
//            }
///////////////////////////////////////////////////////////////////////////////////////////////////





           //Agnese
/////////////////////////////////////////////////////////////////////////////////////////////////////////
           	    hlheHT_log->Fill(TMath::Log(lheHT) ,genweight);
                    hlheHT->Fill(lheHT ,genweight);
                    hlheNj->Fill(lheNpNLO ,genweight);
		        hgen_mass->Fill(gen_mass,genweight);
                        hnGenJet->Fill(nGenJet,genweight);
                        hnGenJetsWithoutLeptons->Fill(genJetsWithoutLeptonsP4.size(),genweight);
                        if (fromNANO) hnGenLep->Fill(ngenLepIdx,genweight);
                        else hnGenLep->Fill(nGenLep + nGenLepRecovered,genweight);
                        hGenLepton_matching1->Fill(genLepton_deltaR1,genweight);
                        hGenLepton_matching2->Fill(genLepton_deltaR2,genweight);
                        
                        hrecoJetPt_genJetPt1->Fill(jet1_deltaGenPT,genweight);
                        hrecoJetPt_genJetPt2->Fill(jet2_deltaGenPT,genweight);
                        
                        hgenJetMass->Fill(genJetMass_leading,genweight);
                        hgenJetMass_matched->Fill(genJetMass_matched,genweight);
                        
				if(abs(GenLep_pdgId[0])==11 && abs(GenLep_pdgId[1])==11 )
				hbinning_MeeJet->Fill(M_ll_cutbin+0.5,lheNj+0.5, genweight);
				
				if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ) {
				hbinning_MmumuJet->Fill(M_ll_cutbin+0.5,lheNj+0.5, genweight);
				
				
				}
				
				if(abs(GenLep_pdgId[0])==15 && abs(GenLep_pdgId[1])==15 )
				hbinning_MttJet->Fill(M_ll_cutbin+0.5,lheNj+0.5, genweight);
				
				hpdgId->Fill(GenLep_pdgId[0],genweight);
				hpdgId->Fill(GenLep_pdgId[1],genweight);
				
////////////////////////////////////////////////////////////////////////////////////////

//            if(file_tag.CompareTo("DYJetstoLL_amc_0J")==0){
//             if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
//                 if((M_ll_cutbin==1)&&(lheNj==0))
//            genweight=genweight*(0.042/4.779);
//            }
//            }
//            if(file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0){
//             if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
//                 if((M_ll_cutbin==1)&&(lheNj==0))
//            genweight=genweight*(0.042/4.779);
//            else genweight=0;
//            }
//            else genweight=0;
//            }


//---------------------STITCHING------------------------------------------------------------



           float ratio = 1.;
//             if(file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0){
//                 if(abs(GenLep_pdgId[0])==11 && abs(GenLep_pdgId[1])==11 ){
//                         if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==2))  ratio = 1.;
//                     }
//                 if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
//                         if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 0.00290597/0.00290597;
//                         if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 0.0172459/0.0171639;
//                         if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 0.0196249/0.0192204;
//                         if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 0.0429222/4.77991;
//                         if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 0.189544/4.31358;
//                         if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 0.211394/1.93457;
//                         if((M_ll_cutbin==2)&&(lheNj==0))  ratio =  0.000469661/0.000469661;
//                         if((M_ll_cutbin==2)&&(lheNj==1))  ratio =  0.00111306/0.00111306;
//                         if((M_ll_cutbin==2)&&(lheNj==2))  ratio =  0.000978254/0.000918971;
//                     }
//                 if(abs(GenLep_pdgId[0])==15 && abs(GenLep_pdgId[1])==15 ){
//                         if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==0))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==1))  ratio = 1.;
//                         if((M_ll_cutbin==2)&&(lheNj==2))  ratio = 1.;
//                     }
//             }

//             if ((file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0) || (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0) ) {
//                 if (VtypeSim < 0.5) {
//                     
//                     ratio = ratiosClass.provideRatio(file_tag, abs(GenLep_pdgId[0]),M_ll_cutbin ,lheNj);
//                     cout << "VtypeSim " << VtypeSim << "ratio " << ratio << endl;
//                 }
//             }


//         bool stitching = (nGenLep + nGenLepRecovered) == 2;
//         
//         if ((file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0) || (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0) ) {
//             
//             if(stitching && abs(GenLep_pdgId[0]) == 13 && abs(GenLep_pdgId[1]) == 13) {
//                 if (file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0) ratio = 1;
//                 else ratio = 0;
//             }
//             else {
//                 if (file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0) ratio = 0;
//                 else ratio = 1;   
//             }
//             
//         }
                
                
                

            
//            genweight = genweight * ratio;

//---------------------END STITCHING------------------------------------------------------------



               if ( (MVAtree_to_fill && MVAcount < MVAcountMAX) || (plotOutput) || (plotOutput_NN)) {
//                 if ( MVAtree_to_fill && MVAcount < MVAcountMAX) {
//                     std::cout << "tree entries    " << treeMVA->GetEntries() << std::endl;
                        TMVA.Inv_mass=Inv_mass;
                        TMVA.Invariant_MassLog=Invariant_MassLog;
                        TMVA.X_parton1=X_parton1;
                        TMVA.X_parton2=X_parton2;
                        TMVA.W_mass_virtual1=abs(W_mass_virtual1);
                        TMVA.W_mass_virtual2=abs(W_mass_virtual2);
                        TMVA.W_Pt_virtual1=W_Pt_virtual1;
                        TMVA.W_Pt_virtual2=W_Pt_virtual2;
                        TMVA.W_eta_virtual1=W_eta_virtual1;
                        TMVA.W_eta_virtual2=W_eta_virtual2;
                        TMVA.W_phi_virtual1=W_phi_virtual1;
                        TMVA.W_phi_virtual2=W_phi_virtual2;
                        TMVA.thetastarW1=thetastarW1;
                        TMVA.thetastarW2=thetastarW2;
                        TMVA.thetastarW2toHW1=thetastarW2toHW1;
                        TMVA.thetastarW1toHW2=thetastarW1toHW2;
                        TMVA.thetastarHtoWW=thetastarHtoWW;
                        TMVA.theta1=theta1;
                        TMVA.theta2=theta2;

                        TMVA.WWmass=WWmass;
                        TMVA.impulsoZ=log(abs(impulsoZ));
                        TMVA.energytot=energytot;
                        TMVA.energytotLog=energytotLog;
                        TMVA.mumujj_pt = log(mumujj_pt);
                        TMVA.lepton1_pt=lepton1.Pt();
                        TMVA.lepton2_pt=lepton2.Pt();
                        TMVA.jets12=Qjet1.Pt() + Qjet2.Pt();
                        
                        TMVA.deltaMRel=deltaMRel;
                    TMVA.deltaM=deltaM;
                    TMVA.normalizedDistance_from_mH=normalizedDistance_from_mH;
                    TMVA.ll_mass = Zll_mass;
                    TMVA.ll_mass_biasUp = Zll_mass_biasUp;
                    TMVA.ll_mass_biasDown = Zll_mass_biasDown;
                    TMVA.ll_mass_resolution = Zll_mass_resolutionUp;
                    TMVA.diffMassWWH=WWmass-Zll_mass;

                    TMVA.q1_eta = Qjet1.Eta();
                    TMVA.met_pt = met_pt;
                    TMVA.EWKHTsoft = Jet.EWKHTsoft - EWKHTsoft_toSubtract;
                    TMVA.EWKHTsoft = HTsoft;
                    TMVA.btagCMVA = maxBTagCMVA;
                    TMVA.btagCMVA_leading = maxBTagCMVA_leading;
                    TMVA.btagCSV= maxBTagCSV;
                    TMVA.btagCMVA_second = maxSecondBTagCMVA;
                    TMVA.btagCSV_second= maxSecondBTagCSV;
                    TMVA.softLeadingJet_pt = Jet.softLeadingJet_pt;

                    TMVA.cosThetaStar=cosThetaStar;
                    TMVA.cosThetaStarAbs=abs(cosThetaStar);
                    TMVA.cosThetaPlane = cosThetaPlane;
                    TMVA.absCosThetaStarJet = absCosThetaStarJet;

                    TMVA.Muon1_relIso04 = selLeptons_relIso04[idx_1stLepton];
                    TMVA.Muon2_relIso04 = selLeptons_relIso04[idx_2ndLepton];
                    
                                        
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    TMVA.qgl_1qAtanh = atanh(1.999995*(Jet.qgl[jets_indices[0]]-0.5));
                    TMVA.qgl_2qAtanh = atanh(1.999995*(Jet.qgl[jets_indices[1]]-0.5));
                    TMVA.cosThetaPlaneAtanh = atanh(2.*cosThetaPlane-1.);
                    TMVA.absCosThetaStarJetAtanh = atanh(1.999995*(absCosThetaStarJet-0.5));
                    TMVA.X_parton1Log=log(X_parton1);
                    TMVA.X_parton2Log=log(X_parton2);
                    TMVA.W_mass_virtual1Log=log(abs(W_mass_virtual1));
                    TMVA.W_mass_virtual2Log=log(abs(W_mass_virtual2));
                    TMVA.W_Pt_virtual1Log=log(W_Pt_virtual1);
                    TMVA.W_Pt_virtual2Log=log(W_Pt_virtual2);                    
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    
                    TMVA.MuonMass_skim=MuonMass_skim;
                    TMVA.qqMass_skim=qqMass_skim;
                    TMVA.genJetMass_leading=genJetMass_leading;
                    
                    TMVA.Mqq = Mqq;
                    TMVA.MqqLog = log(Mqq);
                    
                    TMVA.Jet3_eta=jet3_eta;
                    
                    TMVA.Jet1q_pt = Qjet1.Pt();
                    TMVA.Jet2q_pt = Qjet2.Pt();
                    TMVA.Jet1q_eta = Qjet1.Eta();
                    TMVA.Jet2q_eta = Qjet2.Eta();
                    TMVA.Jet3q_pt = jet3_pt;
                    TMVA.Jet3q_ptLog = log(jet3_pt);
                    TMVA.Jet2q_ptLog = log(Qjet2.Pt());
                    TMVA.Jet1q_ptLog = log(Qjet1.Pt());
                    TMVA.PhiZQ1 =PhiZQ1;
                    TMVA.PhiZQ2 =PhiZQ2;

//                     std::cout << std::setprecision(4) << TMVA.Jet1q_eta << TMVA.genJetPt_match1 << " \t "  << TMVA.Jet1q_pt << " \t " << Jet.pt_JERup[jets_indices[0]] << " \t "  << Jet.pt_JERdown[jets_indices[0]] << std::endl;
                    
//                     std::cout << evt << " \t " << GenJet_pt[Jet.Jet_genJetIdx[jets_indices[0]]] << " \t " << TMVA.Jet1q_pt << " \t " << Jet.pt_JERup[jets_indices[0]] << " \t "  << Jet.pt_JERdown[jets_indices[0]] << " \t " << jets_indices[0] << " \t "  << jets_indices[1] << " \t " << (Jet.pt_JERdown[jets_indices[0]]-TMVA.Jet1q_pt)/(TMVA.Jet1q_pt-TMVA.genJetPt_match1) << " \t " << (Jet.pt_JERdown[jets_indices[0]]-TMVA.Jet1q_pt)/(TMVA.Jet1q_pt-GenJet_pt[Jet.Jet_genJetIdx[jets_indices[0]]]) << std::endl;
                    
//                     printf ("eta %4.2f \t  %4.2f \t  %4.2f \t  %4.2f \t  %4.2f \t", TMVA.Jet1q_eta, TMVA.genJetPt_match1, TMVA.Jet1q_pt, Jet.pt_JERup[jets_indices[0]],Jet.pt_JERdown[jets_indices[0]]);
//                     printf (" %4.5f\n", (Jet.pt_JERdown[jets_indices[0]]-TMVA.Jet1q_pt)/(TMVA.Jet1q_pt-TMVA.genJetPt_match1));
                    
                    TMVA.Jet1q_pt_jerUp = Jet.pt_JERup[jets_indices[0]];
                    TMVA.Jet1q_pt_jerDown = Jet.pt_JERdown[jets_indices[0]];
                    TMVA.Jet2q_pt_jerUp = Jet.pt_JERup[jets_indices[1]];
                    TMVA.Jet2q_pt_jerDown = Jet.pt_JERdown[jets_indices[1]];
        
                    TMVA.Jet1q_chEmEF = Jet.chEmEF[jets_indices[0]];
                    TMVA.Jet2q_chEmEF = Jet.chEmEF[jets_indices[1]];
                    TMVA.Jet1q_neEmEF = Jet.neEmEF[jets_indices[0]];
                    TMVA.Jet2q_neEmEF = Jet.neEmEF[jets_indices[1]];
                    
                    TMVA.Jet1q_chHEF = Jet.chHEF[jets_indices[0]];
                    TMVA.Jet2q_chHEF = Jet.chHEF[jets_indices[1]];
                    TMVA.Jet1q_neHEF = Jet.neHEF[jets_indices[0]];
                    TMVA.Jet2q_neHEF = Jet.neHEF[jets_indices[1]];
                    
                    
                    TMVA.ll_zstar =  (TMath::Abs(Zll_zstar) > 0.000001) ? log(TMath::Abs(Zll_zstar)) : -13.8155105579642736;
                    TMVA.ll_ystar =  Zll_ystar;
    // 		if (TMath::Exp((-1)*Jet_axis2[jets_indices[0]]) < 0.2)  TMVA.axis2_jet1 = TMath::Exp((-1)*Jet_axis2[jets_indices[0]]);
    // 		if (TMath::Exp((-1)*Jet_axis2[jets_indices[1]]) < 0.2) TMVA.axis2_jet2 = TMath::Exp((-1)*Jet_axis2[jets_indices[1]]);
    // 		TMVA.Jet1q_leadTrackPt = Jet_leadTrackPt[jets_indices[0]];
                    TMVA.Rpt = Rpt;
                    TMVA.DeltaEtaQQ = qqDeltaEta;
    // 		TMVA.Jet2q_leadTrackPt = Jet_leadTrackPt[jets_indices[1]];
                    TMVA.qq_pt = log(qq_pt);
                    TMVA.ll_pt = Zll_pt;
                    TMVA.ll_pt_Log = log(Zll_pt);
                    TMVA.ll_eta = Zll_eta;
                    
                    TMVA.Jet3_pt = jet3_pt;
                    TMVA.qgl_1q = Jet.qgl[jets_indices[0]];
                    TMVA.qgl_2q = Jet.qgl[jets_indices[1]];
                    
                    TMVA.minAbsEta = minAbsEta;
                    TMVA.maxAbsEta = maxAbsEta;
                    
                    TMVA.genJetIdx0 = Jet.Jet_genJetIdx[jets_indices[0]];
                    TMVA.genJetIdx1 = Jet.Jet_genJetIdx[jets_indices[1]];
                    if(jets_indices.size()>2) TMVA.genJetIdx2 = Jet.Jet_genJetIdx[jets_indices[2]];
                    for (int n=0; n<jets_indices.size() && n<30 ;n++) {
                        TMVA.genJetIdx[n] = Jet.Jet_genJetIdx[jets_indices[n]];
                        TMVA.Jet_pt[n] = Jet.pt[jets_indices[n]];
                    }
                    
                    if(Jet.qgl[jets_indices[0]]<-0.5) std::cout << entry << " \t jets_indices[0]:  " << jets_indices[0] << " \t " << Jet.qgl[jets_indices[0]] << std::endl;
                    if(Jet.qgl[jets_indices[1]]<-0.5) std::cout << entry << " \t jets_indices[1]:  " << jets_indices[1] << " \t " << Jet.qgl[jets_indices[1]] << std::endl;
                    
                    
                    TMVA.softActivityEWK_njets2 = Nsoft2; // Jet.EWKnsoft2 - EWKnsoft2_toSubtract;
                    TMVA.softActivityEWK_njets5 = Nsoft5; // Jet.EWKnsoft5 - EWKnsoft5_toSubtract;
                    TMVA.softActivityEWK_njets10 = Nsoft10; // Jet.EWKnsoft10 - EWKnsoft10_toSubtract;
//                     std::cout << entry << " TMVA     \t Jet.EWKnsoft5   " << Jet.EWKnsoft5 << " \t EWKnsoft5_toSubtract   " << EWKnsoft5_toSubtract<< " \t TMVA.softActivityEWK_njets5   "  << TMVA.softActivityEWK_njets5 << std::endl;
//                    TMVA.weightMVA = genweight; 
                   TMVA.genweight = genweight; 
                   TMVA.genweight_noQGLcorrection = genweight_noQGLcorrection; 
//                    TMVA.weightMVA = weightMVA;
                   TMVA.weightMVA = 1.;
                   if (file_tag.CompareTo("TT")==0)                   TMVA.weightMVA = 24.22/2.527;
                   if (file_tag.CompareTo("EKW_LLJJ_INT")==0)         TMVA.weightMVA = 24.22/17.13;
                   if (file_tag.CompareTo("EKW_LLJJ_pythia8")==0)     TMVA.weightMVA = 24.22/17.13;
                   if (file_tag.CompareTo("TTTo2L2Nu")==0)            TMVA.weightMVA = 24.22/26.0359;
                   if (file_tag.CompareTo("TTToSemilepton")==0)       TMVA.weightMVA = 24.22/9.2955;
                   TMVA.weightMVA = genweight; 
                   if (file_tag.CompareTo("VBF_HToMuMu")==0)  TMVA.weightMVA = genweight * 891.900; //this factor set the same normalization to bkg DYNLO
//                    if (file_tag.CompareTo("DYJetsToLL_M-105To160-madgraphMLM")==0)  TMVA.weightMVA = genweight * 1.28866; 
//                    if (file_tag.CompareTo("DYJetsToLL_M-105To160_VBFFilter-madgraphMLM")==0)  TMVA.weightMVA = genweight * 1.28866; 
                   
                   for (int n = 0; n <= number_LHE_weights_pdf; n++) TMVA.genweightVECTOR[n] = LHE_weights_pdf_wgt[n];

                           
//                     TMVA.randomVariable = gRandom->Rndm();
                    TMVA.xSection = xsec[file_tag];
                    MVAcount +=1;
//                     if ( MVAtree_to_fill && MVAcount < MVAcountMAX) treeMVA->Fill();   


                    
                    

                    temp_softActivityEWK_njets5 = (float) TMVA.softActivityEWK_njets5;
                    if (plotOutput)  {
//                     if (0)  {
//                        std::cout << "plotOutput  "  << std::endl;
//                         TMVA.ll_mass = 125;
                        BDT_VBF = reader->EvaluateMVA("BDTG");
//                         BDT_VBF = reader->EvaluateMVA("MLP");
                     if ( data!=1 && (atanh((BDT_VBF+1.)/2.)> 0.8)) hlheHT_log_BDTgt1->Fill(TMath::Log(lheHT) ,genweight);
                     if ( data!=1 && (atanh((BDT_VBF+1.)/2.)> 0.8)) hlheHT_BDTgt1->Fill(lheHT ,genweight);
//                         if(atanh((BDT_VBF+1.)/2.)<0.9) continue;
                        
                        TMVA.BDToutput = BDT_VBF;
                        
                        hBDT_VBF->Fill(BDT_VBF-0.01,genweight);
                        hBDT_VBF_atanh->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        hBDT_VBF_atanh_FITTING->Fill(0.);
                        hBDT_VBF_atanh_findBinning->Fill(atanh((BDT_VBF+1.)/2.),genweight);

                        hBDT_VBF_atanh_QGLUp->Fill(atanh((BDT_VBF+1.)/2.),genweight*qgl_weight_bothJet);
                        hBDT_VBF_atanh_QGLDown->Fill(atanh((BDT_VBF+1.)/2.),genweight/qgl_weight_bothJet);
                        
                        if (genJetMass_leading>350) hBDT_VBF_atanh_genMqqGt350->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        else hBDT_VBF_atanh_genMqqLs350->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        if (genJetMass_leading<1) hBDT_VBF_atanh_genMqq0->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        
                        if (genweight<0) hBDT_VBF_atanh_negativeFraction->Fill(atanh((BDT_VBF+1.)/2.));
                        hBDT_VBF_atanh_notWeighted->Fill(atanh((BDT_VBF+1.)/2.));
                        
                        
                        if (atanh((BDT_VBF+1.)/2.) > 0.6 ) hMqq_cut06->Fill(Mqq,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 0.8 ) hMqq_cut08->Fill(Mqq,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.0 ) hMqq_cut10->Fill(Mqq,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.1721 ) hMqq_cut12->Fill(Mqq,genweight);
                        
                        
                        if (atanh((BDT_VBF+1.)/2.) > 0.6 ) hgenJetMass_cut06->Fill(genJetMass_leading,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 0.8 ) hgenJetMass_cut08->Fill(genJetMass_leading,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.0 ) hgenJetMass_cut10->Fill(genJetMass_leading,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.2 ) hgenJetMass_cut12->Fill(genJetMass_leading,genweight);
                        
                        
                        if (atanh((BDT_VBF+1.)/2.) > 0.6 ) hminAbsEta_cut06->Fill(minAbsEta,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 0.8 ) hminAbsEta_cut08->Fill(minAbsEta,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.0 ) hminAbsEta_cut10->Fill(minAbsEta,genweight);
                        if (atanh((BDT_VBF+1.)/2.) > 1.2 ) hminAbsEta_cut12->Fill(minAbsEta,genweight);
                        
                        
                        if (lheNpNLO == 0 ) hBDT_VBF_atanh_cumulativeLHENpNLO_0->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        if (lheNpNLO == 1 ) hBDT_VBF_atanh_cumulativeLHENpNLO_1->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        if (lheNpNLO == 2 ) hBDT_VBF_atanh_cumulativeLHENpNLO_2->Fill(atanh((BDT_VBF+1.)/2.),genweight);                        
                        
                    }
                    

                    

                    
                    if ( MVAtree_to_fill && MVAcount < MVAcountMAX)   treeMVA->Fill();   
         
                    
                    
                    
                    
                    
                    


                    
                    
//                     lwt::ValueMap nnout;
//                     float vars[18] = { TMVA.Mqq, TMVA.DeltaEtaQQ, TMVA.q1_eta, TMVA.ll_pt, TMVA.ll_mass, TMVA.met_pt, TMVA.EWKHTsoft, TMVA.qq_pt, TMVA.RptHard, (float) TMVA.softActivityEWK_njets5, TMVA.btagCMVA, TMVA.btagCSV, TMVA.qgl_1q, TMVA.qgl_2q, TMVA.cosThetaStar, TMVA.cosThetaStarAbs};
                    float vars[18] = { TMVA.Mqq, TMVA.DeltaEtaQQ, TMVA.q1_eta, TMVA.ll_pt, TMVA.ll_mass, TMVA.met_pt, TMVA.EWKHTsoft, TMVA.qq_pt, TMVA.Rpt, (float) TMVA.softActivityEWK_njets5, TMVA.btagCSV, TMVA.qgl_1q, TMVA.qgl_2q, TMVA.ll_zstar, TMVA.ll_ystar, TMVA.Jet1q_pt, TMVA.Jet2q_pt, TMVA.ll_eta};
                    
                    

//                    if (plotOutput_NN && false)  {
// 
//                        
// //                         BDT_VBF = nnResult.result(vars) - 0.0000001;
// //                         BDT_VBF = BDT_VBF + 4.*reader->EvaluateMVA("BDTG");
// //                         std::cout << "BDT_VBF:   " << BDT_VBF << std::endl;
// 
//                         
//                         float  atanhBDT = atanh(BDT_VBF) + 4*atanh((reader->EvaluateMVA("BDTG")+1.)/2.);
// //                         hBDT_VBF->Fill(atanh((BDT_VBF+1.)/2.),genweight);
//                         hBDT_VBF->Fill(atanhBDT,genweight);
// //                         hBDT_VBF_atanh->Fill(atanh((BDT_VBF+1.)/2.),genweight);
//                         hBDT_VBF_atanh->Fill(atanhBDT,genweight);
//                         hBDT_VBF_atanh_findBinning->Fill(atanhBDT,genweight);
//                     }
        
        
        
                }
                
                
                
             

    if ((file_tag.CompareTo("VBF_HToMuMu")==0) &&  ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) && ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) && ((JERWeight_str.CompareTo("none")==0)||(JERWeight_str.CompareTo("nom")==0)) && ((PUWeight_str.CompareTo("none")==0)||(PUWeight_str.CompareTo("nom")==0))) {
        int binNumber = hBDT_VBF_atanh->FindBin(atanh((TMVA.BDToutput+1.)/2.));
//         for (int n = 1; n <= number_LHE_weights_pdf; n++) hBDT_VBF_atanh_PDFvariation[binNumber]->Fill(atanh((TMVA.BDToutput+1.)/2.),genweight);
        for (int n = 0; n <= number_LHE_weights_pdf; n++) {
            hBDT_VBF_atanh_PDFvariation[binNumber]->Fill(n,genweight*LHE_weights_pdf_wgt[n]);
        }
    }

    float BDT_VBF_atanh_m125 = 0;
    
    if (plotOutput)  {

       

        
        TMVA.ll_mass = Zll_mass_biasUp;
        float BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_biasUp->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
                
        TMVA.ll_mass = Zll_mass_biasDown;
//         TMVA.ll_mass = Zll_mass;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_biasDown->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
                
        TMVA.ll_mass = Zll_mass_resolutionUp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_resolutionUp->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        
        TMVA.ll_mass = Zll_mass_resolutionDown;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_resolutionDown->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        
        
                
        TMVA.ll_mass = 125;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        BDT_VBF_atanh_m125 = atanh((BDT_VBF_m125ForAll+1.)/2.);
//         if (atanh((BDT_VBF_m125ForAll+1.)/2.)< 1.) continue;
        hBDT_VBF_atanh_m125ForAll->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        
//         if ((Zll_mass < 120) || (Zll_mass > 130)) hBDT_VBF_atanh_m125ControlRegion->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
//         if (Zll_mass < 120) hBDT_VBF_atanh_m125ControlRegionDown->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
//         if (Zll_mass > 130) hBDT_VBF_atanh_m125ControlRegionUp->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        
        if ((Zll_mass < 115) || (Zll_mass > 135)) hBDT_VBF_atanh_m125ControlRegion->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        if (Zll_mass < 115) hBDT_VBF_atanh_m125ControlRegionDown->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        if (Zll_mass > 135) hBDT_VBF_atanh_m125ControlRegionUp->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        
        
//         TMVA.ll_mass = Zll_mass;
//         if ( (MVAtree_to_fill && MVAcount < MVAcountMAX) || (plotOutput) || (plotOutput_NN)) if ( MVAtree_to_fill && MVAcount < MVAcountMAX) treeMVA->Fill();  
//         TMVA.ll_mass = 125;
  
  
        float tmpVar = TMVA.MqqLog;
        TMVA.MqqLog = 7.5;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_MqqLog75->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.MqqLog = tmpVar;
        
        tmpVar = TMVA.Rpt;
        TMVA.Rpt = 0;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_Rpt0->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.Rpt = tmpVar;
        
        tmpVar = TMVA.W_mass_virtual1Log;
        TMVA.W_mass_virtual1Log = 5.;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_WM1Log5->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.W_mass_virtual1Log = tmpVar;
        
        tmpVar = TMVA.ll_zstar;
        TMVA.ll_zstar = 0;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_zStar0->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.ll_zstar = tmpVar;
        
        tmpVar = temp_softActivityEWK_njets5;
        temp_softActivityEWK_njets5 = 0;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_softN50->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        temp_softActivityEWK_njets5 = tmpVar;
        
        
        tmpVar = TMVA.ll_pt;
        TMVA.ll_pt = 300.;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_ptll300->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.ll_pt = tmpVar;
        
                
        tmpVar = TMVA.mumujj_pt;
        TMVA.mumujj_pt = 30.;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_mumujjPt->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.mumujj_pt = tmpVar;
        
                
        tmpVar = TMVA.theta2;
        TMVA.theta2 = 0.7;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_Theta2->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.theta2 = tmpVar;
        
                 
        tmpVar = TMVA.DeltaEtaQQ;
        TMVA.DeltaEtaQQ = 4.5;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_DEtajj->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.DeltaEtaQQ = tmpVar;
        
        
        
        
        float MqqLog_tmp                         = TMVA.MqqLog;
        float Rpt_tmp                            = TMVA.Rpt;
        float W_mass_virtual1Log_tmp             = TMVA.W_mass_virtual1Log;
        float ll_zstar_tmp                       = TMVA.ll_zstar;
        float temp_softActivityEWK_njets5_tmp    = temp_softActivityEWK_njets5;
        float ll_pt_tmp                          = TMVA.ll_pt;
        float mumujj_pt_tmp                      = TMVA.mumujj_pt;
        float theta2_tmp                         = TMVA.theta2;
        float DeltaEtaQQ_tmp                     = TMVA.DeltaEtaQQ;
        
        
        TMVA.MqqLog = 7.5;
        TMVA.Rpt = 0;
        TMVA.W_mass_virtual1Log = 5.;
        TMVA.ll_zstar = 0;
        temp_softActivityEWK_njets5 = 0;
        TMVA.ll_pt = 300.;
        TMVA.mumujj_pt = 30.;
        TMVA.theta2 = 0.7;
        TMVA.DeltaEtaQQ = 4.5;
        
        
        TMVA.MqqLog = MqqLog_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyMqq->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.MqqLog = 7.5;
        
        TMVA.Rpt = Rpt_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyRpt->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.Rpt = 0.;
        
        TMVA.W_mass_virtual1Log = W_mass_virtual1Log_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyWM1->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.W_mass_virtual1Log = 5.;
        
        TMVA.ll_zstar = ll_zstar_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyzStar->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.ll_zstar = 0.;
        
        temp_softActivityEWK_njets5 = temp_softActivityEWK_njets5_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlysoftN5->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        temp_softActivityEWK_njets5 = 0.;
        
        TMVA.ll_pt = ll_pt_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyptll->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.ll_pt = 300.;
        
        
        TMVA.mumujj_pt = mumujj_pt_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlymumujjPt->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.mumujj_pt = 30.;
        
        TMVA.theta2 = theta2_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyTheta2->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.theta2 = 0.7;
        
        TMVA.DeltaEtaQQ = DeltaEtaQQ_tmp;
        BDT_VBF_m125ForAll = reader->EvaluateMVA("BDTG");
        hBDT_VBF_atanh_m125_onlyDEtajj->Fill(atanh((BDT_VBF_m125ForAll+1.)/2.),genweight);
        TMVA.DeltaEtaQQ = 4.5;
        
        
    }
    

        
    
            
/////////////////////////////////////////////// BTAG FILIPPO ////////////////////////////////////////////////////////////////////////////////
            hBDiscriminator_CSV->Fill(Jet.btagCSV[jets_indices[0]], genweight);
            hBDiscriminator_CSV->Fill(Jet.btagCSV[jets_indices[1]], genweight);
            hBDiscriminator_CMVA->Fill(Jet.btagCMVA[jets_indices[0]], genweight);
            hBDiscriminator_CMVA->Fill(Jet.btagCMVA[jets_indices[1]], genweight);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       
            
            histo_without_qglCorrection->Fill(1.,genweight/qgl_weight_bothJet);
            
            
            hIsolatedElectrons->Fill(nIsolatedElectrons,genweight);
            
            if ((Zll_mass>130) || (Zll_mass<120)) {
                hmaxAbsEta->Fill(maxAbsEta,genweight);
                if (BDT_VBF_atanh_m125 > 0.8) hmaxAbsEta_BDTmFixgt08->Fill(maxAbsEta,genweight);
                if (BDT_VBF_atanh_m125 > 1.) hmaxAbsEta_BDTmFixgt1->Fill(maxAbsEta,genweight);
            }
            hminAbsEta->Fill(minAbsEta,genweight);
            hminAbsEta_GEN->Fill(minAbsEta_GEN,genweight);
            
       
           hlheNpNLO->Fill(lheNpNLO,genweight);
            hdeltaR1->Fill(deltaR1, genweight);
            hdeltaR2->Fill(deltaR2, genweight);
           
            hdeltaM->Fill(deltaM,genweight);
            hdeltaMRel->Fill(deltaMRel, genweight);
            hnormalizedDistance_from_mH->Fill(normalizedDistance_from_mH, genweight);

           hInvariant_Mass->Fill(Inv_mass,genweight);
           hInvariant_MassLog->Fill(Invariant_MassLog,genweight);
           
           hPz->Fill(Sum.Pz(),genweight);
           hPzAbs->Fill(TMath::Abs(Sum.Pz()),genweight);
           hPzAbsLog->Fill(log(TMath::Abs(Sum.Pz())),genweight);
           
           hTotalEnergy->Fill(Sum.E(),genweight);
           hTotalEnergylog->Fill(TMath::Log(Sum.E()),genweight);
           hmumujj_ptLog->Fill(log(mumujj_pt),genweight);
           hmumujj_pt->Fill(mumujj_pt,genweight);
           
           hEnergy_fraction_Parton1->Fill(X_parton1,genweight);
           hEnergy_fraction_Parton2->Fill(X_parton2,genweight);
           hEnergy_fraction_Parton1_log->Fill(log(X_parton1),genweight);
           hEnergy_fraction_Parton2_log->Fill(log(X_parton2),genweight);
           hVirtual_Wmass1->Fill(W_mass_virtual1,genweight);
           hVirtual_Wmass2->Fill(W_mass_virtual2,genweight);
           hVirtual_Wmass1_log->Fill(log(abs(W_mass_virtual1)),genweight);
           hVirtual_Wmass2_log->Fill(log(abs(W_mass_virtual2)),genweight);        
           
           hVirtual_Pt1->Fill(W_Pt_virtual1,genweight);
           hVirtual_Pt2->Fill(W_Pt_virtual2,genweight);
           hVirtual_Pt1_log->Fill(log(W_Pt_virtual1),genweight);
           hVirtual_Pt2_log->Fill(log(W_Pt_virtual2),genweight);
           hVirtual_eta1->Fill(abs(W_eta_virtual1),genweight);
           hVirtual_eta2->Fill(abs(W_eta_virtual2),genweight);
           hVirtual_phi1->Fill(W_phi_virtual1,genweight);
           hVirtual_phi2->Fill(W_phi_virtual2,genweight);
           
           hthetastar_W1->Fill(thetastarW1,genweight);
           hthetastar_W2->Fill(thetastarW2, genweight);
           hthetastar_W2toHW1->Fill(thetastarW2toHW1,genweight);
           hthetastar_W1toHW2->Fill(thetastarW1toHW2,genweight);
           hthetastar_HtoWW->Fill(thetastarHtoWW,genweight);
           
           hTheta_HiggsJ1->Fill(theta1,genweight);
           hTheta_HiggsJ2->Fill(theta2,genweight);
           hParton_M1->Fill(Parton_mass1,genweight);
           hParton_M2->Fill(Parton_mass2,genweight);
           
           hWWmass->Fill(WWmass,genweight);
           hDiffmass->Fill(WWmass-Zll_mass,genweight);
            
            hThetaPlanes->Fill(cosThetaPlane,genweight);
            hThetaStarJet->Fill(absCosThetaStarJet,genweight);
            hThetaPlanesAtanh->Fill(atanh(2.*cosThetaPlane-1),genweight);
            hThetaStarJetAtanh->Fill(atanh(2.*absCosThetaStarJet-1),genweight);   
            
            hThetaStar->Fill(cosThetaStar,genweight);   
            hThetaStarAbs->Fill(abs(cosThetaStar),genweight);   
            hMaxJetBTagCSV->Fill(maxBTagCSV,genweight);
//             hMaxJetBTagCMVA->Fill(atanhCMVA,genweight);
            hMaxJetBTagCMVA->Fill(maxBTagCMVA,genweight);
            hMaxSecondJetBTagCSV->Fill(maxSecondBTagCSV,genweight);
            hMaxSecondJetBTagCMVA->Fill(maxSecondBTagCMVA,genweight);

            hVtype->Fill(v_type,genweight);
            hVtypeSim->Fill(VtypeSim,genweight);
		   	hMqq->Fill(Mqq,genweight);
//                         if ((atanh(TMVA.BDToutput+1)/2) > 0.6 ) hMqq_cut06->Fill(Mqq,genweight);
//                         if ((atanh(TMVA.BDToutput+1)/2) > 0.8 ) hMqq_cut08->Fill(Mqq,genweight);
//                         if ((atanh(TMVA.BDToutput+1)/2) > 1.0 ) hMqq_cut10->Fill(Mqq,genweight);
//                         if ((atanh(TMVA.BDToutput+1)/2) > 1.2 ) hMqq_cut12->Fill(Mqq,genweight);
		   	hMqq_log->Fill(TMath::Log(Mqq),genweight);
			   hqq_pt->Fill(qq_pt,genweight);
			   hEtaQQ->Fill(qqDeltaEta,genweight);
		 	   hPhiQQ->Fill(qqDeltaPhi,genweight);
                           
                         if (deltaMRel<0.025) hHll_mass_precise->Fill(Zll_mass,genweight);
                         else  hHll_mass_unprecise->Fill(Zll_mass,genweight);
                        hZll_mass_biasUp->Fill(Zll_mass_biasUp,genweight);
		   	hZll_mass_biasDown->Fill(Zll_mass_biasDown,genweight);
                        hZll_mass_resolution->Fill(Zll_mass_resolutionUp,genweight);
		   	hZll_mass->Fill(Zll_mass,genweight);
		   	hZll_pt->Fill(Zll_pt,genweight);
		   	hZll_phi->Fill(Zll_phi,genweight);
		   	hZll_eta->Fill(Zll_eta,genweight);
// 			   hHTsoft->Fill(Jet.HTsoft,genweight);
// 			   hSoft_n2->Fill(Jet.nsoft2, genweight);
// 			   hSoft_n5->Fill(Jet.nsoft5, genweight);
// 			   hSoft_n10->Fill(Jet.nsoft10, genweight);
			   hHTsoftEWK->Fill(HTsoft, genweight); //(Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
			   hSoft_n2EWK->Fill(Nsoft2, genweight); // (Jet.EWKnsoft2 - EWKnsoft2_toSubtract, genweight);
			   hSoft_n5EWK->Fill(Nsoft5, genweight); // (Jet.EWKnsoft5 - EWKnsoft5_toSubtract, genweight);
			   hSoft_n10EWK->Fill(Nsoft10, genweight); // (Jet.EWKnsoft10 - EWKnsoft10_toSubtract, genweight);
                hPtSoftJets->Fill(Jet.softLeadingJet_pt, genweight);
				hnPVs->Fill(nPVs,genweight);
				hqgl->Fill(Jet.qgl[jets_indices[0]],genweight);
				hqgl2->Fill(Jet.qgl[jets_indices[1]],genweight);
                                hqgl_noQGLweight->Fill(Jet.qgl[jets_indices[0]],genweight);
                                hqgl_noQGLweight2->Fill(Jet.qgl[jets_indices[1]],genweight/qgl_weight_bothJet);
                                hqglAtanh->Fill(atanh(2.*Jet.qgl[jets_indices[0]]-1.),genweight);
                                hqgl2Atanh->Fill(atanh(2.*Jet.qgl[jets_indices[1]]-1.),genweight/qgl_weight_bothJet);
                        
				hJet1q_pt->Fill(Qjet1.Pt(),genweight);
				hJet1q_eta->Fill(Qjet1.Eta(),genweight);
				hJet1q_ptd->Fill(Jet.ptd[jets_indices[0]],genweight);
				hJet1q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[0]]),genweight);
				hJet1q_mult->Fill(Jet.mult[jets_indices[0]],genweight);
				hJet1q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[0]],genweight);
				hJet1q_phi->Fill(Qjet1.Phi(),genweight);
				hJet2q_pt->Fill(Qjet2.Pt(),genweight);
				hJet2q_eta->Fill(Qjet2.Eta(),genweight);
				hJet2q_ptd->Fill(Jet.ptd[jets_indices[1]],genweight);
				hJet2q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[1]]),genweight);
				hJet2q_mult->Fill(Jet.mult[jets_indices[1]],genweight);
				hJet2q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[1]],genweight);
				hJet2q_phi->Fill(Qjet2.Phi(),genweight);
                                
                                hJet_genJetIdx->Fill(Jet.Jet_genJetIdx[jets_indices[0]],genweight);
                                hJet_genJetIdx->Fill(-10 - Jet.Jet_genJetIdx[jets_indices[1]],genweight);
                                
				hmet->Fill(met_pt,genweight);
				hrho->Fill(rho,genweight);
				hlepton1_pt->Fill(lepton1.Pt(),genweight);
				hlepton2_pt->Fill(lepton2.Pt(),genweight);
				hlepton1_eta->Fill(lepton1.Eta(),genweight);
				hlepton2_eta->Fill(lepton2.Eta(),genweight);
				hlepton1_iso03->Fill(selLeptons_relIso04[idx_1stLepton],genweight);
				hlepton2_iso03->Fill(selLeptons_relIso04[idx_2ndLepton],genweight);
				hHT->Fill(lheHT ,genweight);
				
				
				hweights_weighted->Fill(genweight>0 ? 1 : -1,abs(genweight));
				hweights->Fill(genweight>0 ? 1 : -1);
				

				hDeltaRelQQ->Fill(DeltaRelQQ,genweight);
				hRpt->Fill(Rpt,genweight);
				hRptAtanh->Fill(atanh(2.*Rpt-1.),genweight);
				hEtaQQSum->Fill(DeltaEtaQQSum,genweight);
				hPhiZQ1->Fill(PhiZQ1,genweight);
                                hPhiZQ2->Fill(PhiZQ2,genweight);
                                hEtaHQ1->Fill(EtaHQ1,genweight);
                                hEtaHQ2->Fill(EtaHQ2,genweight);
				hZll_y->Fill(Zll.Rapidity(),genweight);
				hZll_ystar->Fill(Zll_ystar   ,genweight);
				hZll_zstar_log->Fill(log(TMath::Abs(Zll_zstar)),genweight);
				hZll_zstar->Fill(TMath::Abs(Zll_zstar),genweight);
                                hzstar->Fill(TMath::Abs(zstar),genweight);
				hlheV_pt->Fill(lheV_pt  ,genweight);
				hJet3_pt->Fill(jet3_pt ,genweight);	
				if (good_jets>=3) hJet3_pt_log->Fill(log(jets_pv[2].Pt()) ,genweight);	
				if (good_jets>=3) hJet3_eta->Fill(jet3_eta ,genweight);	
                                if (good_jets>=3) hJet3_etaRatio->Fill((2.*jet3_eta + jets_pv[0].Eta() + jets_pv[1].Eta()/abs(jets_pv   [0].Eta() - jets_pv[1].Eta())) ,genweight);	
				if (good_jets>=3) hJet3_pt_new->Fill(jets_pv[2].Pt(),genweight);
				if (good_jets==2) hJet3_pt_new->Fill(10.,genweight);
				float AdJetHT = 0;
				if (good_jets>=3)
					for (int i=2;i<good_jets;i++) {
						if (jets_pv[i].Pt() > 15 ) AdJetHT+=jets_pv[i].Pt();
                                        }
                                        

                                hNjet25->Fill(TMVA.countJet25, genweight);
                                
                                if (Jet.Jet_genJetIdx[jets_indices[0]] == -1 ) hpuId_jetWithoutGenJet->Fill(Jet.puId[jets_indices[0]],genweight);
                                if (Jet.Jet_genJetIdx[jets_indices[1]] == -1 ) hpuId_jetWithoutGenJet->Fill(Jet.puId[jets_indices[1]],genweight);
                                hpuId_AllJets->Fill(Jet.puId[jets_indices[0]],genweight);
                                hpuId_AllJets->Fill(Jet.puId[jets_indices[1]],genweight);
                                
                                
				if (good_jets==2) hAdJetHT->Fill(0.,genweight);
				if (good_jets>=3) hAdJetHT->Fill(AdJetHT,genweight);
				hsoftleadTrackPt->Fill(Jet.EWKsoft_pt[0],genweight);	
				hsoftleadTrackEta->Fill(Jet.EWKsoft_eta[0],genweight);	
			
				hJets12_pt->Fill((jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJets12_pt_log->Fill(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJet1q_pt_log->Fill(TMath::Log(Qjet1.Pt()),genweight);
				hJet2q_pt_log->Fill(TMath::Log(Qjet2.Pt()),genweight);
				hNAdJets->Fill(good_jets, genweight);	

				if  (file_tag.CompareTo("interference")!=0) {
					hbdt->Fill(bdt,genweight);
					hbdt_atanh->Fill(TMath::ATanH((bdt+1)/2),genweight);
					hbdt_atanh2->Fill(TMath::ATanH((bdt+1)/2),genweight);
				}
				else {
					float interference_weight= interference_func->Eval(TMath::Log(Mqq))  ;
					hbdt->Fill(bdt,genweight*interference_weight);
					hbdt_atanh->Fill(TMath::ATanH((bdt+1)/2),genweight*interference_weight);
					hbdt_atanh2->Fill(TMath::ATanH((bdt+1)/2),genweight*interference_weight);
				}
		 

				hprof_htsoft_pu->Fill(nPVs,Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
				hprof_htsoft_pu_rms->Fill(nPVs,Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);

				if (bdt>0.92) {
					hprof_htsoft_pu_bdt->Fill(nPVs,Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
					hprof_htsoft_pu_rms_bdt->Fill(nPVs,Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
				
				//	hveto_jet3pt_denom->Fill(,genweight)
					for (int i=0;i<hveto_jet3pt_denom->GetNbinsX();i++)
						hveto_jet3pt_denom->Fill(hveto_jet3pt_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_ht_denom->GetNbinsX();i++)
						hveto_ht_denom->Fill(hveto_ht_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softht_denom->GetNbinsX();i++)
						hveto_softht_denom->Fill(hveto_softht_denom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softpt_denom->GetNbinsX();i++)
						hveto_softpt_denom->Fill(hveto_softpt_denom->GetBinCenter(i+1),genweight);
			//		cout<<hveto_jet3pt_denom->GetBinCenter(1)<<" , "<<hveto_jet3pt_denom->GetBinCenter(2)<<endl;
			
					hNAdJets_bdt->Fill(good_jets, genweight);	
			        hHTsoftEWK_bdt->Fill(Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
			 	 	hSoft_n2EWK_bdt->Fill(Jet.EWKnsoft2 - EWKnsoft2_toSubtract, genweight);
			  		hSoft_n5EWK_bdt->Fill(Jet.EWKnsoft5 - EWKnsoft5_toSubtract, genweight);
			  		hSoft_n10EWK_bdt->Fill(Jet.EWKnsoft10 - EWKnsoft10_toSubtract, genweight);
					if (good_jets>=3) {
						hJet3_pt_bdt->Fill(jets_pv[2].Pt(),genweight);
						hJet3_eta_bdt->Fill(jets_pv[2].Eta(),genweight);
					}
					if (good_jets==2) hJet3_pt_bdt->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_bdt->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_bdt->Fill(AdJetHT,genweight);	

					hJet1q_eta_bdt->Fill(Qjet1.Eta(),genweight);
					hJet2q_eta_bdt->Fill(Qjet2.Eta(),genweight);
			
					if (good_jets>2) {
						for (int i=0;i<hveto_jet3pt_nom->GetNbinsX();i++)
							if (jets_pv[2].Pt()>hveto_jet3pt_nom->GetBinCenter(i+1)) hveto_jet3pt_nom->Fill(hveto_jet3pt_nom->GetBinCenter(i+1),genweight);
						for (int i=0;i<hveto_ht_nom->GetNbinsX();i++)
							if (AdJetHT>hveto_ht_nom->GetBinCenter(i+1)) hveto_ht_nom->Fill(hveto_ht_nom->GetBinCenter(i+1),genweight);
					}
					for (int i=0;i<hveto_softht_nom->GetNbinsX();i++)
						if (Jet.EWKHTsoft - EWKHTsoft_toSubtract>hveto_softht_nom->GetBinCenter(i+1)) hveto_softht_nom->Fill(hveto_softht_nom->GetBinCenter(i+1),genweight);
					for (int i=0;i<hveto_softpt_nom->GetNbinsX();i++)
						if (Jet.EWKsoft_pt[0] > hveto_softpt_nom->GetBinCenter(i+1)) hveto_softpt_nom->Fill(hveto_softpt_nom->GetBinCenter(i+1),genweight);
				}
				if (bdt>0.84) {
					hNAdJets_bdt2->Fill(good_jets, genweight);	
			   	hHTsoftEWK_bdt2->Fill(Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
			 	 	hSoft_n2EWK_bdt2->Fill(Jet.EWKnsoft2 - EWKnsoft2_toSubtract, genweight);
			  		hSoft_n5EWK_bdt2->Fill(Jet.EWKnsoft5 - EWKnsoft5_toSubtract, genweight);
			  		hSoft_n10EWK_bdt2->Fill(Jet.EWKnsoft10 - EWKnsoft10_toSubtract, genweight);
					if (good_jets>=3) {
						hJet3_pt_bdt2->Fill(jets_pv[2].Pt(),genweight);
						hJet3_eta_bdt2->Fill(jets_pv[2].Eta(),genweight);
					}
					if (good_jets==2) hJet3_pt_bdt2->Fill(10.,genweight);
					if (good_jets==2) hAdJetHT_bdt2->Fill(10.,genweight);
					if (good_jets>=3) hAdJetHT_bdt2->Fill(AdJetHT,genweight);	
					hJet1q_eta_bdt2->Fill(Qjet1.Eta(),genweight);
					hJet2q_eta_bdt2->Fill(Qjet2.Eta(),genweight);
				}
// 				if (Mqq > 1500) {
// 					hNAdJets_mjj1->Fill(good_jets, genweight);	
// 			   	hHTsoftEWK_mjj1->Fill(Jet.EWKHTsoft,genweight);
// 			 	 	hSoft_n2EWK_mjj1->Fill(Jet.EWKnsoft2, genweight);
// 			  		hSoft_n5EWK_mjj1->Fill(Jet.EWKnsoft5, genweight);
// 			  		hSoft_n10EWK_mjj1->Fill(Jet.EWKnsoft10, genweight);
// 					if (good_jets>=3) hJet3_pt_mjj1->Fill(jets_pv[2].Pt(),genweight);
// 					if (good_jets==2) hJet3_pt_mjj1->Fill(10.,genweight);
// 					if (good_jets==2) hAdJetHT_mjj1->Fill(10.,genweight);
// 					if (good_jets>=3) hAdJetHT_mjj1->Fill(AdJetHT,genweight);	
// 				}
// 				if (Mqq > 2500) {
// 					hNAdJets_mjj2->Fill(good_jets, genweight);	
// 			   	hHTsoftEWK_mjj2->Fill(Jet.EWKHTsoft,genweight);
// 			 	 	hSoft_n2EWK_mjj2->Fill(Jet.EWKnsoft2, genweight);
// 			  		hSoft_n5EWK_mjj2->Fill(Jet.EWKnsoft5, genweight);
// 			  		hSoft_n10EWK_mjj2->Fill(Jet.EWKnsoft10, genweight);
// 					if (good_jets>=3) hJet3_pt_mjj2->Fill(jets_pv[2].Pt(),genweight);
// 					if (good_jets==2) hJet3_pt_mjj2->Fill(10.,genweight);
// 					if (good_jets==2) hAdJetHT_mjj2->Fill(10.,genweight);
// 					if (good_jets>=3) hAdJetHT_mjj2->Fill(AdJetHT,genweight);	
// 				}


float singleMu = 0;
float doubleMu = 0;


if (HLT_IsoMu24 || HLT_IsoTkMu24 ) singleMu = 1;
if (postfix.CompareTo("2017")==0) {
//     std::cout << "postfix  " << postfix << std::endl;
    if (HLT_IsoMu27 ) singleMu = 1;
    else singleMu = 0;
}
if (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ ) doubleMu = 1;
// std::cout << "HLT_IsoMu27  " << HLT_IsoMu27 << "  \t   " << singleMu << "  \t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL  " << HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL  << "  \t   " << doubleMu << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////// 2D HISTOS ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        float JSidx1 = indexJetHiggsSister1+0.5;
        float JSidx2 = indexJetHiggsSister2+0.5;

//         JSidx1 = higherJetMVA+0.5;
//         JSidx2 = secondJetMVA+0.5;
        
//         std::vector<float> variables_in_2D_plot = {Zll_mass, deltaM, deltaR1, deltaR2, Rpt, log(Zll_zstar), atanh((TMVA.BDToutput+1)/2), genJetMass_leading, genJetMass_matched, Mqq};
        std::vector<float> variables_in_2D_plot = {/*singleMu, doubleMu, */Zll_pt, qqDeltaPhi, JSidx1, JSidx2};

        
//         for (int i = 0; i < variables_in_2D_plot.size(); i++) variablesName_in_2D_plot.push_back(getName(variables_in_2D_plot[i]));
        
        int idx = 0;
        for (int i = 0; i < variables_in_2D_plot.size(); i++) {
            for (int j = i+1; j < variables_in_2D_plot.size(); j++) {
//                 if (Jet.Jet_genJetIdx[jets_indices[0]] == -1 && Jet.Jet_genJetIdx[jets_indices[1]] == -1)histo2D_vector[idx]->Fill( variables_in_2D_plot[i], variables_in_2D_plot[j], genweight);;
                if (TMVA.countJet25>2) histo2D_vector[idx]->Fill( variables_in_2D_plot[i], variables_in_2D_plot[j], genweight);;
                idx++;
            }
        }
        
    
//         hMll_deltaM->Fill( Zll_mass, deltaM, genweight);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////// END 2D HISTOS ///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




		
		
		
		
		if (genweight>0) pos_weight_presel++;
		float mcweight=genweight*events_generated/xsec[file_tag]; 

		if (genweight>0) gen_pos_weight+=mcweight;
		if (genweight<0) gen_neg_weight+=mcweight;
		if (genweight>0) gen_pos+=genweight0;
		if (genweight<0) gen_neg+=genweight0;
		


				
// 			global_counter++;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// FILL HISTO SYSTEMATICS///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if(histo_Preparation_To_Fit) {
            scaleWeightsUp[0] = puweight/events_generated;
            scaleWeightsDown[0] = puweight/events_generated;

            scaleWeightsUp[1]   = puweightUp/events_generated;
            scaleWeightsDown[1] = puweightDown/events_generated;
            scaleWeightsUp[2]   = puweight/events_generated_muFR_QCDUp*LHE_weights_scale_wgt[4];
            scaleWeightsDown[2] = puweight/events_generated_muFR_QCDDown*LHE_weights_scale_wgt[5];
//            scaleWeightsUp[3] = puweight/events_generated_muF_QCDUp*LHE_weights_scale_wgt[0];
//            scaleWeightsDown[3] = puweight/events_generated_muF_QCDDown*LHE_weights_scale_wgt[1];
//            scaleWeightsUp[4] = puweight/events_generated_muR_QCDUp*LHE_weights_scale_wgt[2];
//            scaleWeightsDown[4] = puweight/events_generated_muR_QCDDown*LHE_weights_scale_wgt[3];

            scaleWeightsUp[3] =  puweightUp/events_generated;
            scaleWeightsDown[3] =  puweightUp/events_generated;
            scaleWeightsUp[4] =  puweightUp/events_generated;
            scaleWeightsDown[4] =  puweightUp/events_generated;

            scaleWeightsUp[5] = puweight/events_generated;
            scaleWeightsDown[5] = puweight/events_generated;
            scaleWeightsUp[6] = puweight/events_generated;
            scaleWeightsDown[6] = puweight/events_generated;

            scaleWeightsUp[7] = puweight/events_generated;
            scaleWeightsDown[7] = puweight/events_generated;
            scaleWeightsUp[8] = puweight/events_generated;
            scaleWeightsDown[8] = puweight/events_generated;
            scaleWeightsUp[9] = puweight/events_generated;
            scaleWeightsDown[9] = puweight/events_generated;

            scaleWeightsUp[0]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[1]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[2]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[3]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[4]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[5]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[6]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[7]*= qgl_weight;//nominal is with puweight
            scaleWeightsUp[8]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[0]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[1]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[2]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[3]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[4]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[5]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[6]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[7]*= qgl_weight;//nominal is with puweight
            scaleWeightsDown[8]*= qgl_weight;//nominal is with puweight


            hbdtUp[0]->Fill(BDT_VBF,genweight*scaleWeightsUp[0]);
            hbdt_atanhUp[0]->Fill(TMath::ATanH((BDT_VBF+1)/2),genweight*scaleWeightsUp[0]);
            for (int current_syst=1;current_syst<Nsyst_NoConst;current_syst++){
				hbdtUp[current_syst]->Fill(BDT_VBF,genweight*scaleWeightsUp[current_syst]);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((BDT_VBF+1)/2),genweight*scaleWeightsUp[current_syst]);
				hbdtDown[current_syst]->Fill(BDT_VBF,genweight*scaleWeightsDown[current_syst]);
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((BDT_VBF+1)/2),genweight*scaleWeightsDown[current_syst]);
            }
        }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// END HISTO SYSTEMATICS////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// std::cout << "HERE FINE LOOP   "<< std::endl;

        }
// ------------------------------------------------------------------------------------END EventLoop----------------------------

//         hSelectionCuts->Scale(xsec[file_tag]/events_generated);



//     fileMVA->cd();


	treeMVA->Write();

	fileMVA.Close();
        
        
        fileMVAJet.cd();
        treeMVAjet->Write();
        fileMVAJet.Close();
        
        

        
        

        


                
                
//        for (int i=0;i<numArray;i++){
//        for (int n=0; n < histArray[i]->GetNbinsX()+2; ++n)
//		    histArray[i]->SetBinError(n, histArray[i]->GetBinError(n)/2.);
//        }

//std::cout << "count1: " << count1 << std::endl;
//std::cout << "count2: " << count2 << std::endl;
//std::cout << "count3: " << count3 << std::endl;
//std::cout << "count4: " << count4 << std::endl;
//std::cout << "count5: " << count5 << std::endl;
//std::cout << "count6: " << count6 << std::endl;
//std::cout << "count7: " << count7 << std::endl;
//std::cout << "count8: " << count8 << std::endl;
//std::cout << "count9: " << count9 << std::endl;
//std::cout << "count10: " << count10 << std::endl;
//std::cout << "count11: " << count11 << std::endl;
//std::cout << "count12: " << count12 << std::endl;
//std::cout << "count13: " << count13 << std::endl;
//std::cout << "count14: " << count14 << std::endl;

///////////////////////////////////// Writing the output file  ////////////////////////////

		//cout << "Number of events that passed the perselection: " << counter<<endl;
		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_JER"+JERWeight_str+"_PU"+PUWeight_str+"_"+heppyVersion+"_"+postfix+".root","recreate");
// 		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+"_noQGLcorrection.root","recreate");


                

                
        histo_without_qglCorrection->Write();
        hSelectionCuts2->Write();
        float qgl_normalization = 1.;
        if(histArray[0]->Integral(0,histArray[0]->GetNbinsX()+1)>0) qgl_normalization = histo_without_qglCorrection->Integral()   /    histArray[0]->Integral(0,histArray[0]->GetNbinsX()+1);
        
    if (hBDT_VBF_atanh_QGLUp->Integral(0,hBDT_VBF_atanh_QGLUp->GetNbinsX()) > 0 )    hBDT_VBF_atanh_QGLUp->Scale(  hBDT_VBF_atanh->Integral(0,hBDT_VBF_atanh->GetNbinsX()+1) /  hBDT_VBF_atanh_QGLUp->Integral(0,hBDT_VBF_atanh_QGLUp->GetNbinsX()+1));
    if (hBDT_VBF_atanh_QGLDown->Integral(0,hBDT_VBF_atanh_QGLDown->GetNbinsX())> 0 ) hBDT_VBF_atanh_QGLDown->Scale(  hBDT_VBF_atanh->Integral(0,hBDT_VBF_atanh->GetNbinsX()+1) / hBDT_VBF_atanh_QGLDown->Integral(0,hBDT_VBF_atanh_QGLDown->GetNbinsX()+1));
    
    
        /*if(QCDScaleWeight_str.CompareTo("nom")==0 && JESWeight_str.CompareTo("nom")==0 && JERWeight_str.CompareTo("nom")==0 && PUWeight_str.CompareTo("nom")==0)*/ qgl_normalization_file << file_tag << " \t\t " << qgl_normalization  << " \t\t " << QCDScaleWeight_str << JESWeight_str << JERWeight_str << PUWeight_str << std::endl;
        
    hBDT_VBF_atanh_negativeFraction->Divide(hBDT_VBF_atanh_notWeighted);
    hBDT_VBF_atanh_negativeFraction->Scale(1./qgl_normalization);
    
            for (int i=0;i<numArray;++i){

//                 histArray[i]->Scale( qgl_normalization );

                
                histArray[i]->SetLineWidth(2);
                histArray[i]->GetYaxis()->SetTitle("N_{events}");
                histArray[i]->GetYaxis()->SetTitleFont(42);
                histArray[i]->GetYaxis()->SetTitleSize(0.060);

                histArray[i]->GetYaxis()->SetTitleOffset(0.8);
                histArray[i]->SetLineColor(kBlue);
                histArray[i]->Draw();
                histArray[i]->Write();
       		}


       		for (int n = 0; n < histo2D_vector.size(); n++) histo2D_vector[n]->Write();
                histo_TrueIndices_MVAIndices->Write();
                histo_TrueIndices_MVAIndices_3Jets->Write();
                
                if ((file_tag.CompareTo("VBF_HToMuMu")==0) &&  ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) && ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) && ((JERWeight_str.CompareTo("none")==0)||(JERWeight_str.CompareTo("nom")==0)) && ((PUWeight_str.CompareTo("none")==0)||(PUWeight_str.CompareTo("nom")==0))) {
                    for (int n = 0; n < hBDT_VBF_atanh_PDFvariation.size(); n++) {
                        hBDT_VBF_atanh_PDFvariation[n]->Write();
                        write_PDF_variation(hBDT_VBF_atanh, hBDT_VBF_atanh_PDFvariation[n], n);
                    }
                }
                
                
                     
            TH1D * hBDT_VBF_atanh_CumulativeLHENpNLO_0 = (TH1D*) hBDT_VBF_atanh_cumulativeLHENpNLO_0->GetCumulative(false);
            TH1D * hBDT_VBF_atanh_CumulativeLHENpNLO_1 = (TH1D*) hBDT_VBF_atanh_cumulativeLHENpNLO_1->GetCumulative(false);
            TH1D * hBDT_VBF_atanh_CumulativeLHENpNLO_2 = (TH1D*) hBDT_VBF_atanh_cumulativeLHENpNLO_2->GetCumulative(false);
                
            
//             hBDT_VBF_atanh_CumulativeLHENpNLO_0->Write();
//             hBDT_VBF_atanh_CumulativeLHENpNLO_1->Write();
//             hBDT_VBF_atanh_CumulativeLHENpNLO_2->Write();
                
            hbinning_MeeJet->Write();
            hbinning_MmumuJet->Write();
            hbinning_MttJet->Write();

			hprof_htsoft_pu->SetLineWidth(2); 
			hprof_htsoft_pu->SetLineColor(kBlue); 
			hprof_htsoft_pu->Draw();
			hprof_htsoft_pu->Write();
			hprof_htsoft_pu_rms->Draw();
			hprof_htsoft_pu_rms->Write();
			hprof_htsoft_pu_bdt->SetLineWidth(2); 
			hprof_htsoft_pu_bdt->SetLineColor(kBlue); 
			hprof_htsoft_pu_bdt->Draw();
			hprof_htsoft_pu_bdt->Write();
			hprof_htsoft_pu_rms_bdt->Draw();
			hprof_htsoft_pu_rms_bdt->Write();







        if(histo_Preparation_To_Fit) {
            for (int current_syst=0;current_syst<Nsyst_NoConst;current_syst++){

                hbdtUp[current_syst]->Draw();
                hbdtUp[current_syst]->Write();
                if (current_syst!=0) hbdtDown[current_syst]->Draw();
                if (current_syst!=0) hbdtDown[current_syst]->Write();

                hbdt_atanhUp[current_syst]->Draw();
                hbdt_atanhUp[current_syst]->Write();
                if (current_syst!=0)  hbdt_atanhDown[current_syst]->Draw();
                if (current_syst!=0)  hbdt_atanhDown[current_syst]->Write();
            }
        }

    		file.Write();

    		file.Close();


	 ofstream out(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_JER"+JERWeight_str+"_PU"+PUWeight_str+"_"+heppyVersion+"_"+postfix+".txt");
	out<< "positive pure selected = "<<gen_pos<<"  , positive weighted selected =  "<<gen_pos_weight<<" , negative pure selected = "<<gen_neg<< ", negative weighted selected = "<<gen_neg_weight<< ", all evetns in the begining = "<<events_generated<<" , xsec = "<<xsec[file_tag]<<endl;
	out<<"positive weight in so many events : "<<  pos_weight_presel<<endl;

        
	for (int i=0;i<7;i++)
		out<<cut_flow_names[i]<<"\t";
	out<<endl;
	for (int i=0;i<7;i++){
		cut_flow[i]=cut_flow[i];
		out<<cut_flow[i]<<"\t";
	}
	out<<endl;

// return 0;

}
