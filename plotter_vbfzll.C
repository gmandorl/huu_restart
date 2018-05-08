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

#include "parse_json.hh"
#include "LightweightNeuralNetwork.hh"
#include "NNet_comparison.h"


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
const int njets = 30;




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
    Int_t nsoftActivityEWKJets;
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

} Jets;


typedef struct {
        
    
        Float_t ll_mass;
//         Float_t Zll_pt;
   	Float_t Mqq; 
        Float_t MqqLog;
	Float_t DeltaEtaQQ;
  	Float_t q1_eta;
    Float_t met_pt;
    Float_t EWKHTsoft;
    Float_t btagCMVA;
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
	Float_t Jet1q_leadTrackPt;
	Float_t ll_ystar;
	Float_t ll_zstar;
	Float_t Rpt;
	Float_t ll_pt;
	Float_t ll_eta;
	Float_t qgl_1q;
	Float_t qgl_2q;
        Float_t qgl_1qAtanh;
        Float_t qgl_2qAtanh;
        
        Float_t Jet1q_ptLog;
        Float_t Jet2q_ptLog;
        Float_t Jet3q_ptLog;
        
        Int_t softActivityEWK_njets2;
        Int_t softActivityEWK_njets5;
	Int_t softActivityEWK_njets10;
        
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
        
        
        Float_t cosThetaPlaneAtanh;
        Float_t absCosThetaStarJetAtanh;
        Float_t X_parton1Log;
        Float_t X_parton2Log;
        Float_t W_mass_virtual1Log;
        Float_t W_mass_virtual2Log;
        Float_t W_Pt_virtual1Log;
        Float_t W_Pt_virtual2Log;
        
        
        Float_t BDToutput;
        
        float weightMVA;
        float genweight;	

        
}TMVAstruct;




using namespace std;


int main(int argc, char* argv[]){

//std::cout << "I can see here" << std::endl;
//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
//gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");

TrackTag nnResult;

TString file_name = std::string(argv[1]);
TString file_tag = std::string(argv[2]);
TString region = std::string(argv[3]); 
int data = atoi(argv[4]);
int applyQCDScaleWeight = atoi(argv[5]);
TString QCDScaleWeight_str = std::string(argv[6]);
int applyJESWeight = atoi(argv[7]);
TString JESWeight_str = std::string(argv[8]);
TString heppyVersion = std::string(argv[9]);
TString postfix = std::string(argv[10]);
TString output = std::string(argv[11]);

bool MVAtree_to_fill = false;
if ( ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) && ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) ) MVAtree_to_fill = true;
MVAtree_to_fill = true;

std::map <TString, float> xsec;
std::map <TString, float> qgl_norm;
xsec["SingleMuon"] = 1.;
xsec["SingleMuonB"] = 1.;
xsec["SingleMuonC"] = 1.;
xsec["SingleMuonD"] = 1.;
xsec["SingleMuonE"] = 1.;
xsec["SingleMuonF"] = 1.;
xsec["SingleMuonG"] = 1.;
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



xsec["DYJetsToLL_M-50_VBFFilter-amcatnloFXFX"] = 41.81;
xsec["DYJetsToLL_M-50_VBFFilter-madgraphMLM"] = 41.25;
xsec["DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX"] = 41.81;
xsec["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"] = 41.25;


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
xsec["interference"]=1.664;

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

xsec["TTTo2L2Nu"]       = 809. * 32.72 * 32.72 / 10000.;
xsec["TTToSemilepton"]  = 809. * 2. * 32.72 * 67.28  / 10000.;



//k factors 1.21 are not included// 
/*xsec["WJetsToLNu_HT100To200"] = 1345 ;
xsec["WJetsToLNu_HT200To400"] = 359.7  ;
xsec["WJetsToLNu_HT400To600"] = 48.91;
xsec["WJetsToLNu_HT600To800"] =12.05;
xsec["WJetsToLNu_HT800To1200"] = 5.501;
xsec["WJetsToLNu_HT1200To2500"] = 1.329;
xsec["WJetsToLNu_HT2500ToInf"] = 0.03216;
*/
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
// qgl_norm["EWK_LLJJ"]=0.939774091;
// qgl_norm["EWK_LLJJ_herwig"]=1;
// qgl_norm["TT"]=1.069615388;
// qgl_norm["WW"]=0.927930632;
// qgl_norm["WZ"]=0.967820083;
// qgl_norm["ZZ"]=0.964110717;
// qgl_norm["ST_tW_antitop"]=1.012466766;
// qgl_norm["ST_tW_top"]=0.990673372;
// qgl_norm["ST_s-channel"]=0.911006075;
// qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.986731357;
// qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=0.994925429;
// qgl_norm["WJetsToLNu"]=1;
// qgl_norm["DYJetstoLL_amc_M-50"]=1;      
// qgl_norm["DYJetstoLL_amc_0J"]=1.002031915;
// qgl_norm["DYJetstoLL_amc_1J"]=0.966710372;
// qgl_norm["DYJetstoLL_amc_2J"]=0.954117783;
// 
// 
// qgl_norm["DYJetsToLL_Pt-0To50_amc"] = 1.; 
// qgl_norm["DYJetsToLL_Pt-50To100_amc"] = 1.; 
// qgl_norm["DYJetsToLL_Pt-100To250_amc"] = 1.; 
// qgl_norm["DYJetsToLL_Pt-250To400_amc"] = 1. ; 
// qgl_norm["DYJetsToLL_Pt-400To650_amc"] = 1. ; 
// qgl_norm["DYJetsToLL_Pt-650ToInf_amc"] = 1. ;
// 
// 
// qgl_norm["DYJetsToLL_M-105To160-madgraphMLM"]=1.;
// 
// qgl_norm["DYJetstoLL_madgraph"]=1.;
// qgl_norm["DYJetstoLL_amc_Filter105"]=1.;
// qgl_norm["DYJetstoLL_amc_Inclusive"]=1.;
// qgl_norm["DYJetsToLL_M"]=1.;
// qgl_norm["DY1JetsToLL_M"]=1.;
// qgl_norm["DY2JetsToLL_M"]=1.;
// qgl_norm["DY3JetsToLL_M"]=1.;
// qgl_norm["DY4JetsToLL_M"]=1.;
// 
// qgl_norm["VBF_HToMuMu"]=1.;
// qgl_norm["GluGlu_HToMuMu"]=1.;






qgl_norm["DYJetsToLL_M"]=0.972332;
qgl_norm["DY0JetsToLL_M"]=0.996343;
qgl_norm["DY1JetsToLL_M"]=0.970908;
qgl_norm["DY2JetsToLL_M"]=0.963111;
qgl_norm["DY3JetsToLL_M"]=0.953453;
qgl_norm["DY4JetsToLL_M"]=0.981916;
qgl_norm["DYJetsToLL_M-105To160-madgraphMLM"]=0.967472;

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

qgl_norm["DYJetstoLL_amc_0J"]=0.99609;
qgl_norm["DYJetstoLL_amc_1J"]=0.966311;
qgl_norm["DYJetstoLL_amc_2J"]=0.961447;
qgl_norm["DYJetstoLL_amc_Filter105"]=0.964836;
qgl_norm["DYJetsToLL_M-105To160-amcatnloFXFX"]=0.964836;
qgl_norm["DYJetstoLL_amc_M-50"]=0.965673;

qgl_norm["DYJetsToLL_M-50_VBFFilter-amcatnloFXFX"] = 1.;
qgl_norm["DYJetsToLL_M-50_VBFFilter-madgraphMLM"] = 1.;
qgl_norm["DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX"] = 1.;
qgl_norm["DYJetsToLL_M-105To160_VBFFilter-madgraphMLM"] = 1.;


qgl_norm["GluGlu_HToMuMu"]=0.987766;
qgl_norm["VBF_HToMuMu"]=0.974783;

qgl_norm["ST_s-channel"]=1.45921;
qgl_norm["ST_t-channel_antitop_4f_inclusiveDecays"]=1.01835;
qgl_norm["ST_t-channel_top_4f_inclusiveDecays"]=0.999555;
qgl_norm["ST_tW_antitop"]=0.993577;
qgl_norm["ST_tW_top"]=1.00096;
qgl_norm["TT"]=0.994732;

qgl_norm["WJetsToLNu"]=1;
qgl_norm["WW"]=1.00899;
qgl_norm["WZ"]=0.941571;
qgl_norm["ZZ"]=0.930012;



qgl_norm["WToLNu_0J"]  =  1.;
qgl_norm["WToLNu_1J"]  =  1.;
qgl_norm["WToLNu_2J"]  =  1.;

qgl_norm["WWTo1L1Nu2Q"] = 1.;
qgl_norm["WZTo1L1Nu2Q"] = 0.941571;
qgl_norm["ZZTo2L2Q"]    = 0.930012;

qgl_norm["TTTo2L2Nu"]       = 1.06;
qgl_norm["TTToSemilepton"]  = 1.06;



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
if ((JESWeight_str.CompareTo("up")==0) && (QCDScaleWeight_str.CompareTo("up")==0)) whichJERWeight=1;
if ((JESWeight_str.CompareTo("down")==0) && (QCDScaleWeight_str.CompareTo("down")==0)) whichJERWeight=2;


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
	Float_t puweight;
      Float_t puweightUp;
	Float_t puweightDown;
	Float_t PU=1.;
	Float_t genweight;
	Float_t bTagWeight;
	Float_t genweight0;
	float  trigWeight_tree;
// 	Int_t global_counter = 0;
	Int_t HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200;
	Int_t HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460;
	Int_t HLT_IsoMu22;
	Int_t HLT_IsoTkMu22;
	Int_t HLT_IsoMu27;
	Int_t HLT_IsoTkMu27;
	Int_t HLT_IsoTkMu24;
	Int_t HLT_IsoMu24;
	Int_t HLT_Ele27_eta2p1;
	TFile *file_initial;
	TChain *tree_initial;

	Int_t nvLeptons, nselLeptons;
	const int brLeptons=13;
	Float_t vLeptons_pt[30], vLeptons_eta[30], vLeptons_phi[30], vLeptons_mass[30], vLeptons_SF_IdCutLoose[30], vLeptons_SF_IdCutTight[30], vLeptons_SF_IsoLoose[30], vLeptons_SF_IsoTight[30],vLeptons_SF_trk_eta[30], vLeptons_SF_HLT_RunD4p2[30],vLeptons_SF_HLT_RunD4p3[30], vLeptons_relIso03[30], vLeptons_eleSieie[30], vLeptons_eleHoE[30], vLeptons_eleDEta[30],vLeptons_eleDPhi[30], vLeptons_eleEcalClusterIso[30], vLeptons_eleHcalClusterIso[30],vLeptons_dr03TkSumPt[30]  ;
	Int_t vLeptons_charge[30], vLeptons_pdgId[30],vLeptons_trackerLayers[30] ; 

	Float_t selLeptons_pt[30], selLeptons_eta[30], selLeptons_phi[30], selLeptons_mass[30], selLeptons_SF_IdCutLoose[30], selLeptons_SF_IdCutTight[30], selLeptons_SF_IsoLoose[30], selLeptons_SF_IsoTight[30],selLeptons_SF_trk_eta[30], selLeptons_SF_HLT_RunD4p2[30],selLeptons_SF_HLT_RunD4p3[30], selLeptons_relIso04[30], selLeptons_relIso03[30], selLeptons_eleSieie[30], selLeptons_eleHoE[30], selLeptons_eleDEta[30],selLeptons_eleDPhi[30], selLeptons_eleEcalClusterIso[30], selLeptons_eleHcalClusterIso[30],selLeptons_dr03TkSumPt[30] ;


    Int_t nGenLep;
	Float_t GenLep_pt[30];
	Float_t GenLep_eta[30];
	Float_t GenLep_phi[30];
	Float_t GenLep_mass[30];
	Int_t GenLep_pdgId[30];
	
    Int_t nGenLepRecovered;
	Float_t GenLepRecovered_pt[30];
	Float_t GenLepRecovered_eta[30];
	Float_t GenLepRecovered_phi[30];
	Float_t GenLepRecovered_mass[30];
	Int_t GenLepRecovered_pdgId[30];
	
        

    Int_t nGenJet;
        Float_t GenJet_pt[30];
        Float_t GenJet_eta[30];
        Float_t GenJet_phi[30];
        Float_t GenJet_mass[30];

        
	Int_t selLeptons_charge[30], selLeptons_pdgId[30], selLeptons_looseIdPOG[30], selLeptons_trackerLayers[30],  selLeptons_eleMVAIdSppring16GenPurp[30]; 

	TString str_leptons[brLeptons] = {"vLeptons_pt", "vLeptons_eta", "vLeptons_phi", "vLeptons_mass", "vLeptons_charge", "vLeptons_pdgId", "vLeptons_SF_IdCutLoose", "vLeptons_SF_IdCutTight", "vLeptons_SF_IsoLoose","vLeptons_SF_IsoTight","vLeptons_SF_trk_eta","vLeptons_SF_HLT_RunD4p2","vLeptons_SF_HLT_RunD4p3"};


	
////////////////////////////
//	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
//	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
//	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix.root");
//	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix");
//	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix.root");
///	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix");
//	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
//	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
/*	TFile* file_id_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunBCDEF.root");
	TH2F* id_mu_bf = (TH2F*)file_id_mu_bf->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
	TFile* file_id_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunGH.root");
	TH2F* id_mu_aft = (TH2F*)file_id_mu_aft->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");

	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunBCDEF.root");
	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunGH.root");
	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");

	TFile* file_iso_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunBCDEF.root");
	TH2F* iso_mu_bf = (TH2F*)file_iso_mu_bf->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
	TFile* file_iso_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunGH.root");
	TH2F* iso_mu_aft = (TH2F*)file_iso_mu_aft->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");

//	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_MVAIDWP80_80x.root");
//	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_ScaleFactor_MVAIDWP80_80x");
	TFile* file_tracker_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_tracker_80x.root");
	TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");
	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Tight27AfterIDISO.root");
	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_Tight27AfterIDISO");
	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_EIDISO_ZH.root");
	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_EIDISO_ZH");

	TFile* file_track_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunBCDEF.root");
	TH1F* track_mu_bf = (TH1F*)file_track_mu_bf->Get("TriggerEffMap_Graph");
	TFile* file_track_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunGH.root");
	TH1F* track_mu_aft = (TH1F*)file_track_mu_aft->Get("TriggerEffMap_Graph");
*/
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

        TFile* f_param= TFile::Open("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/out.root");
        TH2F *hmuon= (TH2F*)f_param->Get("PtErrParametrization");
        
        cout<< "fparam "<<f_param<<endl;
        cout<< "hmuon  "<<hmuon<<endl;
// 	RoccoR  *rc = new RoccoR("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/rcdata.2016.v3/");
        RoccoR  *rc = new RoccoR("2016/rcdata.2016.v3/");


	
	file_initial = TFile::Open(file_name);
	
	tree_initial = (TChain*)file_initial->Get("tree");
	Int_t events_generated;
	Float_t events_generated_muF_QCDUp;
	Float_t events_generated_muF_QCDDown;
	Float_t events_generated_muR_QCDUp;
	Float_t events_generated_muR_QCDDown;
	Float_t events_generated_muFR_QCDUp;
	Float_t events_generated_muFR_QCDDown;
	TH1F *countPos;
	TH1F *countNeg;
	TH1F *countLHEScale;
	TH1F *countLHEPdf;
	TH1F *countWeighted;
	if ((data!=1)){
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
    else events_generated = 1;

//	if (file_tag.CompareTo("EWK_LLJJ")==0)  events_generated=events_generated/2.00703;
    Jets Jet;
    Float_t v_type;
    Float_t v_typeSim;
    Float_t wrong_type=0.;
    Int_t nJets;
	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	
	Float_t lheHT;	
	Float_t lheNj_f;
        int lheNj;
	Float_t lheV_pt;	
        Float_t lheNpNLO;

	Float_t bdt;
	int passSel;
	int passSel_JESR[4];
	Float_t bdt_JESR[4];

	Float_t met_pt;
	Float_t met_phi;

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

Float_t LHE_weights_pdf_wgt[103];
Float_t LHE_weights_scale_wgt[10];
	
	float V_mass;
	ULong64_t evt;

	ULong64_t event; 
//         ratioDYFilter ratiosClass;
//         ratiosClass.setNumber();


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
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	tree_initial->SetBranchAddress("selLeptons_eleHoE",selLeptons_eleHoE);
	tree_initial->SetBranchAddress("selLeptons_eleDEta",selLeptons_eleDEta);
	tree_initial->SetBranchAddress("selLeptons_eleDPhi",selLeptons_eleDPhi);
	tree_initial->SetBranchAddress("selLeptons_eleEcalClusterIso",selLeptons_eleEcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_eleHcalClusterIso",selLeptons_eleHcalClusterIso);
	tree_initial->SetBranchAddress("selLeptons_dr03TkSumPt",selLeptons_dr03TkSumPt);
	
	tree_initial->SetBranchAddress("evt",&evt);

    float BDT_VBF = 0;
//	tree_initial->SetBranchAddress("selLeptons_pt",&BDT_VBF);
//	tree_initial->SetBranchAddress("BDT_VBF",&BDT_VBF);

//	for (int i=0;i<brLeptons;i++){
//		tree_initial->SetBranchAddress(str_leptons[i],);
//	}	

	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);
	tree_initial->SetBranchAddress("lheHT",&lheHT);
	tree_initial->SetBranchAddress("lheNj",&lheNj_f);
	tree_initial->SetBranchAddress("lheV_pt",&lheV_pt);
 	tree_initial->SetBranchAddress("hlheNpNLO",&lheNpNLO);       
        
//	tree_initial->SetBranchAddress("BDT_VBF",&bdt);



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

	
	
	
	

    
    
    
    TH1F *hdeltaR1 = new TH1F("hdeltaR1","", 30, 0.,5.);
    hdeltaR1->GetXaxis()->SetTitle("#Delta R (lead reco jet, genJet)");
    TH1F *hdeltaR2 = new TH1F("hdeltaR2","", 30, 0.,5.);
    hdeltaR2->GetXaxis()->SetTitle("#Delta R (sublead reco jet, genJet)");
    
    
    TH1F *hHiggsSister1_Leading_eta = new TH1F("hHiggsSister1_Leading_eta","", 20, 0.,1.);
    hHiggsSister1_Leading_eta->GetXaxis()->SetTitle("SelectionCut");
    TH1F *hHiggsSister2_Subleading_eta = new TH1F("hHiggsSister2_Subleading_eta","", 20, 0.,1.);
    hHiggsSister2_Subleading_eta->GetXaxis()->SetTitle("SelectionCut");
    
        TH1F *hHiggsSister1_Leading_phi = new TH1F("hHiggsSister1_Leading_phi","", 20, 0.,3.5);
    hHiggsSister1_Leading_phi->GetXaxis()->SetTitle("SelectionCut");
    TH1F *hHiggsSister2_Subleading_phi = new TH1F("hHiggsSister2_Subleading_phi","", 20, 0.,3.5);
    hHiggsSister2_Subleading_phi->GetXaxis()->SetTitle("SelectionCut");
    
    TH1F *hHiggsSister1_Leading_R = new TH1F("hHiggsSister1_Leading_R","", 20, 0.,1.);
    hHiggsSister1_Leading_R->GetXaxis()->SetTitle("SelectionCut");
    TH1F *hHiggsSister2_Subleading_R = new TH1F("hHiggsSister2_Subleading_R","", 20, 0.,1.);
    hHiggsSister2_Subleading_R->GetXaxis()->SetTitle("SelectionCut");
	
	
    TH1F *hSelectionCuts = new TH1F("hSelectionCuts","", 20, -0.5,19.5);
    hSelectionCuts->GetXaxis()->SetTitle("SelectionCut");
    
    TH1F *hlheNpNLO = new TH1F("hlheNpNLO","", 6, -0.5,5.5);
    hlheNpNLO->GetXaxis()->SetTitle("SelectionCut"); 

    TH1F *hVtype = new TH1F("hVtype","", 7,-1.,6.);
    hVtype->GetXaxis()->SetTitle("vtype");

    TH1F *hVtypeSim = new TH1F("hVtypeSim","", 21,-10.5,10.5);
    hVtypeSim->GetXaxis()->SetTitle("vtypeSim");
    
    TH1F *hmumujj_pt = new TH1F("hmumujj_pt","", 100,0.,500);
    hmumujj_pt->GetXaxis()->SetTitle("pT(#mu#mu jj)");
    TH1F *hmumujj_ptLog = new TH1F("hmumujj_ptLog","", 20,0.,10);
    hmumujj_ptLog->GetXaxis()->SetTitle("log(pT(#mu#mu jj))");
    
    
	TH1F *hMqq = new TH1F("hMqq","",60,0.,3000.);
	hMqq->GetXaxis()->SetTitle("m(qq) (GeV)");
        TH1F *hMqq_cut06 = new TH1F("hMqq_cut06","",60,0.,3000.);
	hMqq_cut06->GetXaxis()->SetTitle("m(qq) (GeV)");
        TH1F *hMqq_cut08 = new TH1F("hMqq_cut08","",60,0.,3000.);
	hMqq_cut08->GetXaxis()->SetTitle("m(qq) (GeV)");
        TH1F *hMqq_cut10 = new TH1F("hMqq_cut10","",60,0.,3000.);
	hMqq_cut10->GetXaxis()->SetTitle("m(qq) (GeV)");
        TH1F *hMqq_cut12 = new TH1F("hMqq_cut12","",60,0.,3000.);
	hMqq_cut12->GetXaxis()->SetTitle("m(qq) (GeV)");
        
        
        
	TH1F *hMqq_log = new TH1F("hMqq_log","",150.,0.,15.);
	hMqq_log->GetXaxis()->SetTitle("ln(m(qq)) (GeV)");
	TH1F *hqq_pt = new TH1F("hqq_pt","",40.,0.,500.);
	hqq_pt->GetXaxis()->SetTitle("p_{T}(qq) (GeV)");
	TH1F *hlepton1_pt = new TH1F("hlepton1_pt","",40.,0.,400.);
	hlepton1_pt->GetXaxis()->SetTitle("leading lepton p_{T} (GeV)");
	TH1F *hlepton2_pt = new TH1F("hlepton2_pt","",30.,0.,300.);
	hlepton2_pt->GetXaxis()->SetTitle("subleading lepton p_{T} (GeV)");
	TH1F *hlepton1_eta = new TH1F("hlepton1_eta","",80.,-4.,4.);
	hlepton1_eta->GetXaxis()->SetTitle("leading lepton #eta");
	TH1F *hlepton2_eta = new TH1F("hlepton2_eta","",80.,-4.,4.);
	hlepton2_eta->GetXaxis()->SetTitle("subleading lepton #eta");
  

	TH1F *hlepton1_iso03 = new TH1F("hlepton1_iso03","",40.,0.,0.3);
	hlepton1_iso03->GetXaxis()->SetTitle("leading lepton iso");
	
	TH1F *hlepton2_iso03= new TH1F("hlepton2_iso03","",40.,0.,0.3);
	hlepton2_iso03->GetXaxis()->SetTitle("subleading lepton iso");


    
    
 
	TH1F *hEtaQQ = new TH1F("hEtaQQ","",90,0.,9.);
	hEtaQQ->GetXaxis()->SetTitle("|#Delta#eta_{qq}|");
	
	TH1F *hbdt = new TH1F("hbdt","",100,-1.,1.);
	hbdt->GetXaxis()->SetTitle("BDT output");
	TH1F *hbdt_atanh = new TH1F("hbdt_atanh","",500,0.,5.);
	hbdt_atanh->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	float bining[50];
	bining[0]=0.;
	for (int i=1;i<31;i++)
		bining[i]=bining[i-1]+0.1;
	bining[31] = 3.5;
			
	TH1F *hbdt_atanh2 = new TH1F("hbdt_atanh2","",31,bining);
	hbdt_atanh2->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	TH1F *hPhiQQ = new TH1F("hPhiQQ","",32,0.,3.2);
	hPhiQQ->GetXaxis()->SetTitle("|#Delta#phi_{qq}|");
    
	
	TH1F *hEtaSoftJets = new TH1F("hEtaSoftJets","",12,-3.,3.);
	hEtaSoftJets->GetXaxis()->SetTitle("|#eta^{soft}|");
    
	
	TH1F *hMassSoftJets = new TH1F("hMassSoftJets","",10,0.,100.);
	hMassSoftJets->GetXaxis()->SetTitle("m^{soft}");


	TH1F *hSoft_n2 = new TH1F("hSoft_n2","",25,0.,25.);
	hSoft_n2->GetXaxis()->SetTitle("N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5 = new TH1F("hSoft_n5","",10,0.,10.);
	hSoft_n5->GetXaxis()->SetTitle("N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10 = new TH1F("hSoft_n10","",6,0.,6.);
	hSoft_n10->GetXaxis()->SetTitle("N soft jets, p_{T} > 10 GeV");

	TH1F *hHTsoftEWK = new TH1F("hHTsoftEWK","",30,0.,300.);
	hHTsoftEWK->GetXaxis()->SetTitle("EWK H_{T}^{soft} (GeV)" );
	TH1F *hSoft_n2EWK = new TH1F("hSoft_n2EWK","",25,0.,25.);
	hSoft_n2EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5EWK = new TH1F("hSoft_n5EWK","",10,0.,10.);
	hSoft_n5EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10EWK = new TH1F("hSoft_n10EWK","",6,0.,6.);
	hSoft_n10EWK->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV");
	
	TH1F *hHTsoftEWK_bdt = new TH1F("hHTsoftEWK_bdt","",30,0.,300.);
	hHTsoftEWK_bdt->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.92 (GeV)" );
	TH1F *hSoft_n2EWK_bdt = new TH1F("hSoft_n2EWK_bdt","",25,0.,25.);
	hSoft_n2EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.92");
	TH1F *hSoft_n5EWK_bdt = new TH1F("hSoft_n5EWK_bdt","",10,0.,10.);
	hSoft_n5EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.92");
	TH1F *hSoft_n10EWK_bdt = new TH1F("hSoft_n10EWK_bdt","",6,0.,6.);
	hSoft_n10EWK_bdt->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.92");


	TH1F *hHTsoftEWK_bdt2 = new TH1F("hHTsoftEWK_bdt2","",30,0.,300.);
	hHTsoftEWK_bdt2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , BDT>0.84 (GeV)" );
	TH1F *hSoft_n2EWK_bdt2 = new TH1F("hSoft_n2EWK_bdt2","",25,0.,25.);
	hSoft_n2EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , BDT>0.84");
	TH1F *hSoft_n5EWK_bdt2 = new TH1F("hSoft_n5EWK_bdt2","",10,0.,10.);
	hSoft_n5EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , BDT>0.84");
	TH1F *hSoft_n10EWK_bdt2 = new TH1F("hSoft_n10EWK_bdt2","",6,0.,6.);
	hSoft_n10EWK_bdt2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , BDT>0.84");

// 	TH1F *hHTsoftEWK_mjj1 = new TH1F("hHTsoftEWK_mjj1","",30,0.,300.);
// 	hHTsoftEWK_mjj1->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 1500 (GeV)" );
// 	TH1F *hSoft_n2EWK_mjj1 = new TH1F("hSoft_n2EWK_mjj1","",25,0.,25.);
// 	hSoft_n2EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 1500");
// 	TH1F *hSoft_n5EWK_mjj1 = new TH1F("hSoft_n5EWK_mjj1","",10,0.,10.);
// 	hSoft_n5EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 1500");
// 	TH1F *hSoft_n10EWK_mjj1 = new TH1F("hSoft_n10EWK_mjj1","",6,0.,6.);
// 	hSoft_n10EWK_mjj1->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 1500");
// 	
// 	TH1F *hHTsoftEWK_mjj2 = new TH1F("hHTsoftEWK_mjj2","",30,0.,300.);
// 	hHTsoftEWK_mjj2->GetXaxis()->SetTitle("EWK H_{T}^{soft} , m(qq) > 2500 (GeV)" );
// 	TH1F *hSoft_n2EWK_mjj2 = new TH1F("hSoft_n2EWK_mjj2","",25,0.,25.);
// 	hSoft_n2EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 2 GeV , m(qq) > 2500");
// 	TH1F *hSoft_n5EWK_mjj2 = new TH1F("hSoft_n5EWK_mjj2","",10,0.,10.);
// 	hSoft_n5EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 5 GeV , m(qq) > 2500");
// 	TH1F *hSoft_n10EWK_mjj2 = new TH1F("hSoft_n10EWK_mjj2","",6,0.,6.);
// 	hSoft_n10EWK_mjj2->GetXaxis()->SetTitle("EWK N soft jets, p_{T} > 10 GeV , m(qq) > 2500");



        TH1F* hqgl = new TH1F("hqgl","",20.,0.,1.);
        hqgl->GetXaxis()->SetTitle("QGL 1^{st} q-jet");

        TH1F* hqgl2 = new TH1F("hqgl2","",20.,0.,1.);
        hqgl2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");
        
        TH1F* hqglAtanh = new TH1F("hqglAtanh","",42.,-7.,7.);
        hqglAtanh->GetXaxis()->SetTitle("tanh^{-1}(QGL 1^{st} q-jet)");

        TH1F* hqgl2Atanh = new TH1F("hqgl2Atanh","",42.,-7.,7.);
        hqgl2Atanh->GetXaxis()->SetTitle("tanh^{-1}(QGL 2^{nd} q-jet)");
	
	TH1F *hPtSoftJets = new TH1F("hPtSoftJets","",30,0.,300);
	hPtSoftJets->GetXaxis()->SetTitle("p_{T}^{soft} (GeV)");
        TH1F *hPtSoftJets2 = new TH1F("hPtSoftJets2", "", 20, 0., 200.);
        hPtSoftJets2->GetXaxis()->SetTitle("2nd Soft Jet p_{T} (GeV)");
        TH1F *hPtSoftJets3 = new TH1F("hPtSoftJets3", "", 20, 0., 200.);
        hPtSoftJets3->GetXaxis()->SetTitle("3rd Soft Jet p_{T} (GeV)");
	
	TH1F *hcosOqqbb = new TH1F("hcosOqqbb","",100,-1.,1.);
	hcosOqqbb->GetXaxis()->SetTitle("cos(#theta_{bb_qq})");
	TH1F *hEtaQB1 = new TH1F("hEtaQB1","",160.,-8.,8.);
	hEtaQB1->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}");
	TH1F *hEtaQB2 = new TH1F("hEtaQB2","",160.,-8.,8.);
	hEtaQB2->GetXaxis()->SetTitle("#Delta#eta_{qb}^{backward}");
	TH1F *hPhiQB1 = new TH1F("hPhiQB1","",32,0.,3.2);
	hPhiQB1->GetXaxis()->SetTitle("#Delta#phi_{qb}^{forward}");
	TH1F *hPhiQB2 = new TH1F("hPhiQB2","",32,0.,3.2);
	hPhiQB2->GetXaxis()->SetTitle("#Delta#phi_{qb}^{backward}");
	TH1F *hx1 = new TH1F("hx1","",100.,0.,1.);
	hx1->GetXaxis()->SetTitle("x_{1}");
	TH1F *hx2 = new TH1F("hx2","",100.,0.,1.);
	hx2->GetXaxis()->SetTitle("x_{2}");
	TH1F *hVB1_mass = new TH1F("hVB1_mass","",100,0.,1000.);
	hVB1_mass->GetXaxis()->SetTitle("M_{W'_{1}} (GeV)");
	TH1F *hVB2_mass = new TH1F("hVB2_mass","",100.,0.,1000.);
	hVB2_mass->GetXaxis()->SetTitle("M_{W'_{2}} (GeV)");

	TH1F* hEtot = new TH1F("hEtot","",150.,0.,6000.);
	hEtot->GetXaxis()->SetTitle("E^{tot} (GeV)");
	TH1F* hPxtot= new TH1F("hPxtot","",100,-500.,500.);
	hPxtot->GetXaxis()->SetTitle("P_{x}^{tot} (GeV)");
	TH1F* hPytot= new TH1F("hPytot","",100,-500.,500.);
	hPytot->GetXaxis()->SetTitle("P_{y}^{tot} (GeV)");
	TH1F* hPztot= new TH1F("hPztot","",100,-5000.,5000);
	hPztot->GetXaxis()->SetTitle("P_{z}^{tot} (GeV)");

	
	TH1F *hPtqqll = new TH1F("hPtqqll","",50.,0.,500.);
	hPtqqll->GetXaxis()->SetTitle("p_{T} of qqll system (GeV)");
	TH1F *hPhiqqll = new TH1F("hPhiqqll","",32,-3.2,3.2);
	hPhiqqll->GetXaxis()->SetTitle("-#phi of qqll system");
	TH1F *hEtaqqll = new TH1F("hEtaqqll","",160,0,8);
	hEtaqqll->GetXaxis()->SetTitle("#eta of qqll system");

	TH1F *hnPVs = new TH1F("hPVs","",50,0,50);
	hnPVs->GetXaxis()->SetTitle("nPVs");
	
    TH1F *hdeltaM=new TH1F("hdeltaM","",20,0,20);
    hdeltaM->GetXaxis()->SetTitle("hdeltaM");
    TH1F *hdeltaMRel=new TH1F("hdeltaMRel","",20,0,0.2);
    hdeltaMRel->GetXaxis()->SetTitle("hdeltaMRel");
    TH1F *hnormalizedDistance_from_mH=new TH1F("hdnormalizedDistance_from_mH","",40,0,3);
    hnormalizedDistance_from_mH->GetXaxis()->SetTitle("hnormalizedDistance_from_mH");
    
	TH1F* hZll_mass = new TH1F("hZll_mass","",50,100.,150.);
	hZll_mass->GetXaxis()->SetTitle("m(ll) (GeV)");
	TH1F* hZll_pt = new TH1F("hZll_pt","",40,0.,400.);
	hZll_pt->GetXaxis()->SetTitle("p_{T}(ll) (GeV)");
	TH1F* hZll_eta = new TH1F("hZll_eta","",20,-5,5.);
	hZll_eta->GetXaxis()->SetTitle("#eta(ll)");
	TH1F* hZll_phi = new TH1F("hZll_phi","",32,-3.2,3.2);
	hZll_phi->GetXaxis()->SetTitle("#phi(ll)");

	TH1F *hJet1q_pt = new TH1F("hJet1q_pt","",60,0,500.);
	hJet1q_pt->GetXaxis()->SetTitle("p_{T} 1^{st} q-jet");
	TH1F *hJet1q_eta = new TH1F("hJet1q_eta","",20,-5,5);
	hJet1q_eta->GetXaxis()->SetTitle("#eta 1^{st} q-jet");
	TH1F *hJet1q_ptd = new TH1F("hJet1q_ptd","",100,0,1);
	hJet1q_ptd->GetXaxis()->SetTitle("ptd 1^{st} q-jet");
	TH1F *hJet1q_axis2= new TH1F("hJet1q_axis2","",80,0.,0.16);
	hJet1q_axis2->GetXaxis()->SetTitle("#sigma_{2} 1^{st} q-jet");
	TH1F *hJet1q_mult= new TH1F("hJet1q_mult","",30,0,30.);
	hJet1q_mult->GetXaxis()->SetTitle("N 1^{st} q-jet");
	TH1F *hJet1q_leadTrackPt= new TH1F("hJet1q_leadTrackPt","",20,0,100.);
	hJet1q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 1^{st} q-jet");
	TH1F *hJet1q_leadTrackEta= new TH1F("hJet1q_leadTrackEta","",20,-5,5.);
	hJet1q_leadTrackEta->GetXaxis()->SetTitle("leading track #eta 1^{st} q-jet");
	TH1F *hJets12_pt = new TH1F("hJets12_pt","",60,0,600.);
	hJets12_pt->GetXaxis()->SetTitle("|p_{T}_{1q} + p_{T}_{2q}| (GeV)");
	TH1F *hJets12_pt_log = new TH1F("hJets12_pt_log","",60,3,9);
	hJets12_pt_log->GetXaxis()->SetTitle("ln|p_{T}_{1q} + p_{T}_{2q}| (GeV)");
	TH1F *hJet1q_phi = new TH1F("hJet1q_phi","",32,-3.2,3.2);
	hJet1q_phi->GetXaxis()->SetTitle("#phi 1^{st} q-jet");
	TH1F *hJet2q_phi = new TH1F("hJet2q_phi","",32,-3.2,3.2);
	hJet2q_phi->GetXaxis()->SetTitle("#phi 2^{nd} q-jet");


	TH1F *hJet2q_pt = new TH1F("hJet2q_pt","",60,0,300);
	hJet2q_pt->GetXaxis()->SetTitle("p_{T} 2^{nd} q-jet");
	TH1F *hJet2q_eta = new TH1F("hJet2q_eta","",20,-5,5);
	hJet2q_eta->GetXaxis()->SetTitle("#eta 2^{nd} q-jet");
	TH1F *hJet2q_ptd = new TH1F("hJet2q_ptd","",100,0,1);
	hJet2q_ptd->GetXaxis()->SetTitle("ptd 2^{nd} q-jet");
	TH1F *hJet2q_axis2= new TH1F("hJet2q_axis2","",80,0.,0.16);
	hJet2q_axis2->GetXaxis()->SetTitle("#sigma_{2} 2^{nd} q-jet");

	TH1F *hist_bins = new TH1F("bins","",80,0.,0.16);
	
	TH1F *hJet2q_pt_log = new TH1F("hJet2q_pt_log","",60,2,7);
	hJet2q_pt_log->GetXaxis()->SetTitle("log p_{T} 2^{nd} q-jet");
	TH1F *hJet1q_pt_log = new TH1F("hJet1q_pt_log","",60,2,7);
	hJet1q_pt_log->GetXaxis()->SetTitle("log p_{T} 1^{st} q-jet");
	


	TH1F *hJet2q_mult= new TH1F("hJet2q_mult","",30,0,30.);
	hJet2q_mult->GetXaxis()->SetTitle("N 2^{nd} q-jet");
	TH1F *hJet2q_leadTrackPt= new TH1F("hJet2q_leadTrackPt","",20,0,100.);
	hJet2q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 2^{nd} q-jet");
	
	TH1F *hJet3_pt = new TH1F("hJet3_pt","",18,20,200);
	hJet3_pt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
        TH1F *hJet3_pt_log = new TH1F("hJet3_pt_log","",60,2,7);
	hJet3_pt_log->GetXaxis()->SetTitle("log(p_{T} 3^{rd} jet)");
        
	TH1F *hJet3_pt_new = new TH1F("hJet3_pt_new","",13,0,195);
	hJet3_pt_new->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");
	TH1F *hJet3_eta = new TH1F("hJet3_eta","",20,-5,5);
	hJet3_eta->GetXaxis()->SetTitle("#eta 3^{rd} jet");
	TH1F *hJet3_eta_bdt = new TH1F("hJet3_eta_bdt","",20,-5,5);
	hJet3_eta_bdt->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.92");
	TH1F *hJet3_eta_bdt2 = new TH1F("hJet3_eta_bdt2","",20,-5,5);
	hJet3_eta_bdt2->GetXaxis()->SetTitle("#eta 3^{rd} jet, BDT > 0.84");
        
        TH1F *hJet3_etaRatio = new TH1F("hJet3_etaRatio","",20,-5,5);
	hJet3_etaRatio->GetXaxis()->SetTitle("#etaRatio 3^{rd} jet");
	
	TH1F *hsoftleadTrackPt= new TH1F("hsoftleadTrackPt","",40,0,200.);
	hsoftleadTrackPt->GetXaxis()->SetTitle("leading track p_{T}");
	TH1F *hsoftleadTrackEta= new TH1F("hsoftleadTrackEta","",20,-5,5.);
	hsoftleadTrackEta->GetXaxis()->SetTitle("leading s track #eta ");
	
	TH1F *hAdJetHT = new TH1F("hAdJetHT","",62,0,500);
	hAdJetHT->GetXaxis()->SetTitle("additional jets H_{T} (GeV)");


	TH1F *hmet = new TH1F("hmet","",40,0.,250.);
	hmet->GetXaxis()->SetTitle("MET p_{T} (GeV)");
	TH1F *hrho = new TH1F("hrho","",60,0.,30.);
	hrho->GetXaxis()->SetTitle("rho");
	TH1F *hHT = new TH1F("hHT","",50,0.,1000.);
	hHT->GetXaxis()->SetTitle("lhe H_{T} (GeV)" );
	TH1F *hlheHT_log = new TH1F("hlheHT_log","",100,0.,10.);
	hlheHT_log->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );
	TH1F *hlheNj = new TH1F("hlheNj","",6,0.,6);
	hlheNj->GetXaxis()->SetTitle("lhe N jets" );

	TH1F *hDeltaRelQQ = new TH1F("hDeltaRelQQ","",25.,0.,1.);
	hDeltaRelQQ->GetXaxis()->SetTitle("#Delta_{rel}(qq)");
	TH1F *hRpt = new TH1F("hRpt","",25.,0.,1.);
	hRpt->GetXaxis()->SetTitle("R(p_{T})");
        TH1F *hRptAtanh = new TH1F("hRptAtanh","",35.,-3.5,3.5);
	hRptAtanh->GetXaxis()->SetTitle("tanh^{-1}(R(p_{T}))");
        
	TH1F *hEtaQQSum = new TH1F("hEtaQQSum","",90,0.,9.);
	hEtaQQSum->GetXaxis()->SetTitle("|#eta_{q_{1}}| + |#eta_{q_{2}}| ");
	TH1F *hPhiZQ1 = new TH1F("hPhiZQ1","",32,0.,3.2);
	hPhiZQ1->GetXaxis()->SetTitle("|#Delta#phi(ll,q_{1})|");
	TH1F* hZll_y = new TH1F("hZll_y","",20,-4,4.);
	hZll_y->GetXaxis()->SetTitle("y(ll)");
	TH1F* hZll_ystar = new TH1F("hZll_ystar","",20,-6,6.);
	hZll_ystar->GetXaxis()->SetTitle("y*(ll)");
	TH1F* hZll_zstar = new TH1F("hZll_zstar","",20,-8,3.);
	hZll_zstar->GetXaxis()->SetTitle("log(z*(ll))");
	TH1F* hlheV_pt = new TH1F("hlheV_pt","",40,0.,400.);
	hlheV_pt->GetXaxis()->SetTitle("lheV_pt (GeV)");


	TH1F *hJet3_pt_bdt = new TH1F("hJet3_pt_bdt","",13,0,195);
	hJet3_pt_bdt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.92 (GeV)");
	TH1F *hAdJetHT_bdt = new TH1F("hAdJetHT_bdt","",30,0,450);
	hAdJetHT_bdt->GetXaxis()->SetTitle("additional jets HT, BDT>0.92 (GeV)");
	TH1F *hNAdJets_bdt = new TH1F("hNAdJets_bdt","",10,0,10);
	hNAdJets_bdt->GetXaxis()->SetTitle("N of jets, BDT>0.92");
	TH1F *hNAdJets = new TH1F("hNAdJets","",10,0,10);
	hNAdJets->GetXaxis()->SetTitle("N of jets");

	TH1F *hJet3_pt_bdt2 = new TH1F("hJet3_pt_bdt2","",13,0,195);
	hJet3_pt_bdt2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, BDT>0.84 (GeV)");
	TH1F *hAdJetHT_bdt2 = new TH1F("hAdJetHT_bdt2","",30,0,450);
	hAdJetHT_bdt2->GetXaxis()->SetTitle("additional jets HT, BDT>0.84 (GeV)");
	TH1F *hNAdJets_bdt2 = new TH1F("hNAdJets_bdt2","",10,0,10);
	hNAdJets_bdt2->GetXaxis()->SetTitle("N of jets, BDT>0.84");
// 	TH1F *hJet3_pt_mjj1 = new TH1F("hJet3_pt_mjj1","",13,0,195);
// 	hJet3_pt_mjj1->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 1500 (GeV)");
// 	TH1F *hAdJetHT_mjj1 = new TH1F("hAdJetHT_mjj1","",30,0,450);
// 	hAdJetHT_mjj1->GetXaxis()->SetTitle("additional jets HT, m(qq) > 1500 (GeV)");
// 	TH1F *hNAdJets_mjj1 = new TH1F("hNAdJets_mjj1","",10,0,10);
// 	hNAdJets_mjj1->GetXaxis()->SetTitle("N of jets, m(qq) > 1500");
// 	TH1F *hJet3_pt_mjj2 = new TH1F("hJet3_pt_mjj2","",13,0,195);
// 	hJet3_pt_mjj2->GetXaxis()->SetTitle("p_{T} 3^{rd} jet, m(qq) > 2500 (GeV)");
// 	TH1F *hAdJetHT_mjj2 = new TH1F("hAdJetHT_mjj2","",30,0,450);
// 	hAdJetHT_mjj2->GetXaxis()->SetTitle("additional jets HT, m(qq) > 2500 (GeV)");
// 	TH1F *hNAdJets_mjj2 = new TH1F("hNAdJets_mjj2","",10,0,10);
// 	hNAdJets_mjj2->GetXaxis()->SetTitle("N of jets, m(qq) > 2500 ");


	TH1F *hJet1q_eta_bdt = new TH1F("hJet1q_eta_bdt","",20,-5,5);
	hJet1q_eta_bdt->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.92");
	TH1F *hJet2q_eta_bdt = new TH1F("hJet2q_eta_bdt","",20,-5,5);
	hJet2q_eta_bdt->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.92");
	
	TH1F *hJet1q_eta_bdt2 = new TH1F("hJet1q_eta_bdt2","",20,-5,5);
	hJet1q_eta_bdt2->GetXaxis()->SetTitle("#eta 1^{st} q-jet, BDT > 0.84");
	TH1F *hJet2q_eta_bdt2 = new TH1F("hJet2q_eta_bdt2","",20,-5,5);
	hJet2q_eta_bdt2->GetXaxis()->SetTitle("#eta 2^{nd} q-jet, BDT > 0.84");
	
	TH1F *hpdgId= new TH1F("hpdgId","",40,-20,20);
	hpdgId->GetXaxis()->SetTitle("type of particle");
    
    TH1F *hweights= new TH1F("hweights","",300,-1.5,1.5);
    hweights->GetXaxis()->SetTitle("sign(weight)");
	 
    TH1F *hweights_weighted= new TH1F("hweights_weighted","",300,-1.5,1.5);
    hweights_weighted->GetXaxis()->SetTitle("sign(weight)weighted");
    
	TH1F *hveto_jet3pt_nom = new TH1F("hveto_jet3pt_nom","",9,0,270);
	TH1F *hveto_jet3pt_denom = new TH1F("hveto_jet3pt_denom","",9,0,270);
	TH1F *hveto_ht_nom = new TH1F("hveto_ht_nom","",14,0,420);
	TH1F *hveto_ht_denom = new TH1F("hveto_ht_denom","",14,0,420);
	TH1F *hveto_softht_nom = new TH1F("hveto_softht_nom","",8,0,320);
	TH1F *hveto_softht_denom = new TH1F("hveto_softht_denom","",8,0,320);
	TH1F *hveto_softpt_nom = new TH1F("hveto_softpt_nom","",6,0,180);
	TH1F *hveto_softpt_denom = new TH1F("hveto_softpt_denom","",6,0,180);


        
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
        
        const int binNUmberBDT = 10;
        float BDT_bin[binNUmberBDT] = {0, 0.17175,        0.32175,        0.47175,        0.62175,        0.77175,        0.92175,      1.07175,         1.22175,        1.5};
        
        
        
     TH1F *hBDT_VBF_atanh_findBinning = new TH1F("hBDT_VBF_atanh_findBinning","",10000, 0, 1.5);
     hBDT_VBF_atanh_findBinning->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
        
    TH1F *hBDT_VBF_atanh = new TH1F("hBDT_VBF_atanh","",binNUmberBDT-1, BDT_bin);
    hBDT_VBF_atanh->GetXaxis()->SetTitle("tanh^{-1}(( BDT output + 1.)/2.) ");
    

    
    int bin_new_NumberBDT = 50;
    float upperLimitBDT = 1.;//bin_new_NumberBDT*0.2;
    TH1F *hBDT_VBF = new TH1F("hBDT_VBF","",bin_new_NumberBDT,0.,upperLimitBDT);
    hBDT_VBF->GetXaxis()->SetTitle(" BDT output  ");
    
    TH1F *hThetaPlanes = new TH1F("hThetaPlanes","",60,-0.1,1.1);
    hThetaPlanes->GetXaxis()->SetTitle("cos(#theta_{(#mu#mu)(jj)})"); 
    TH1F *hThetaPlanesAtanh = new TH1F("hThetaPlanesAtanh","",40,-6.,6.);
    hThetaPlanesAtanh->GetXaxis()->SetTitle("tanh^{-1}(cos(#theta_{(#mu#mu)(jj)}))"); 
    
    TH1F *hThetaStarJet = new TH1F("hThetaStarJet","",60,-0.1,1.1);
    hThetaStarJet->GetXaxis()->SetTitle("cos(#theta_{jj}*)"); 
    TH1F *hThetaStarJetAtanh = new TH1F("hThetaStarJetAtanh","",40,-6.,6.);
    hThetaStarJetAtanh->GetXaxis()->SetTitle("tanh^{-1}(cos(#theta_{jj}*))"); 
    
    TH1F *hThetaStar = new TH1F("hThetaStar","",100,-1.1,1.1);
    hThetaStar->GetXaxis()->SetTitle("cos(#theta*)");
    
    TH1F *hThetaStarAbs = new TH1F("hThetaStarAbs","",50,0,1.1);
    hThetaStarAbs->GetXaxis()->SetTitle("|cos(#theta*)|");


    TH1F *hBDiscriminator_CSV = new TH1F("hBDiscriminator_CSV","",60,-0.3,1.);
    hBDiscriminator_CSV->GetXaxis()->SetTitle("max CSV");
    
    
    
    TH1F *hMaxJetBTagCSV = new TH1F("hMaxJetBTagCSV","",60,-0.3,1.);
    hMaxJetBTagCSV->GetXaxis()->SetTitle("max CSV");
    TH1F *hMaxSecondJetBTagCSV = new TH1F("hMaxSecondJetBTagCSV","",60,-0.3,1.);
    hMaxSecondJetBTagCSV->GetXaxis()->SetTitle("second max CSV");

    TH1F *hMaxJetBTagCMVA = new TH1F("hMaxJetBTagCMVA","",50,-1,1.);
    hMaxJetBTagCMVA->GetXaxis()->SetTitle("tanh^{-1}( max CMVA )");
    TH1F *hMaxSecondJetBTagCMVA = new TH1F("hMaxSecondJetBTagCMVA","",50,-1.,1.);
    hMaxSecondJetBTagCMVA->GetXaxis()->SetTitle("second max CMVA ");
    
    
    TH1F *hBDiscriminator_CMVA = new TH1F("hBDiscriminator_CMVA","",50,-1,1.);
    hBDiscriminator_CMVA->GetXaxis()->SetTitle("tanh^{-1}( max CMVA )");
    
    
    TH1F * hminAbsEta = new TH1F("hminAbsEta","",40,0,4.);
    hminAbsEta->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    TH1F * hminAbsEta_cut06 = new TH1F("hminAbsEta_cut06","",40,0,4.);
    hminAbsEta_cut06->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1F * hminAbsEta_cut08 = new TH1F("hminAbsEta_cut08","",40,0,4.);
    hminAbsEta_cut08->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1F * hminAbsEta_cut10 = new TH1F("hminAbsEta_cut10","",40,0,4.);
    hminAbsEta_cut10->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    TH1F * hminAbsEta_cut12 = new TH1F("hminAbsEta_cut12","",40,0,4.);
    hminAbsEta_cut12->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    
    TH1F * hminAbsEta_GEN = new TH1F("hminAbsEta_GEN","",40,0,4.);
    hminAbsEta_GEN->GetXaxis()->SetTitle("min(|#eta(j1)|, |#eta(j2)|)");
    
    
    ///Agnese
    
    TH1F *hInvariant_Mass=new TH1F("hInvariant_Mass","",60,0,6000);
    hInvariant_Mass->GetXaxis()->SetTitle("M(jj#mu#mu) [GeV]");
    TH1F *hInvariant_MassLog=new TH1F("hInvariant_MassLog","",50,5,10);
    hInvariant_MassLog->GetXaxis()->SetTitle("log(M(jj#mu#mu)) [GeV]");
    
    TH1F *hTotalEnergy = new TH1F("hTotalEnergy","",50,0,6000);
    hTotalEnergy->GetXaxis()->SetTitle("E(jj#mu#mu) [GeV]");
    TH1F *hTotalEnergylog = new TH1F("hTotalEnergylog","",50,4,10);
    hTotalEnergylog->GetXaxis()->SetTitle("E(jj#mu#mu) [GeV]");
    
    TH1F *hPz= new TH1F("hPz","",50,-5000,5000);
    hPz->GetXaxis()->SetTitle("Pz(jj#mu#mu) [GeV]");
    TH1F *hPzAbs= new TH1F("hPzAbs","",50,0,5000);
    hPzAbs->GetXaxis()->SetTitle("|Pz(jj#mu#mu)| [GeV]");
    TH1F *hPzAbsLog= new TH1F("hPzAbsLog","",50,0,10);
    hPzAbsLog->GetXaxis()->SetTitle("log(|Pz(jj#mu#mu)|)");
   
    
    
    TH1F *hEnergy_fraction_Parton1=new TH1F("hEnergy_fraction_Parton1","",100,0,1);
    hEnergy_fraction_Parton1->GetXaxis()->SetTitle("X(parton1)");
    TH1F *hEnergy_fraction_Parton2=new TH1F("hEnergy_fraction_Parton2","",100,0,1);
    hEnergy_fraction_Parton2->GetXaxis()->SetTitle("X(parton2)");
    TH1F *hEnergy_fraction_Parton1_log=new TH1F("hEnergy_fraction_Parton1_log","",100,-6,0);
    hEnergy_fraction_Parton1_log->GetXaxis()->SetTitle("log(X(parton1))");
    TH1F *hEnergy_fraction_Parton2_log=new TH1F("hEnergy_fraction_Parton2_log","",100,-6,0);
    hEnergy_fraction_Parton2_log->GetXaxis()->SetTitle("log(X(parton2))");
    
    TH1F *hVirtual_Wmass1=new TH1F("hVirtual_Wmass1","",100,-600,100);
    hVirtual_Wmass1->GetXaxis()->SetTitle("M(virtual_mass1) [GeV]");
    TH1F *hVirtual_Wmass2=new TH1F("hVirtual_Wmass2","",100,-600,100);
    hVirtual_Wmass2->GetXaxis()->SetTitle("M(virtual_mass2) [GeV]");
    
    TH1F *hVirtual_Wmass1_log=new TH1F("hVirtual_Wmass1_log","",100,2,8);
    hVirtual_Wmass1_log->GetXaxis()->SetTitle("M(virtual_mass1_log) [GeV]");
    TH1F *hVirtual_Wmass2_log=new TH1F("hVirtual_Wmass2_log","",100,2,8);
    hVirtual_Wmass2_log->GetXaxis()->SetTitle("M(virtual_mass2_log) [GeV]");   
    TH1F *hWWmass=new TH1F("hWWmass","",100,-100,300);
    hWWmass->GetXaxis()->SetTitle("M(W+W) [GeV]");
    TH1F *hDiffmass=new TH1F("hDiffmass","",100,-250,300);
    hDiffmass->GetXaxis()->SetTitle("M(W+W)-M_ll [GeV]");
    
    TH1F *hVirtual_Pt1=new TH1F("hVirtual_Pt1","",80,0,400);
    hVirtual_Pt1->GetXaxis()->SetTitle("Pt(W1) [GeV]");
    TH1F *hVirtual_Pt2=new TH1F("hVirtual_Pt2","",80,0,400);
    hVirtual_Pt2->GetXaxis()->SetTitle("Pt(W2) [GeV]");
    
    TH1F *hVirtual_Pt1_log=new TH1F("hVirtual_Pt1_log","",80,2,8);
    hVirtual_Pt1_log->GetXaxis()->SetTitle("log(Pt(W1))");
    TH1F *hVirtual_Pt2_log=new TH1F("hVirtual_Pt2_log","",80,2, 8);
    hVirtual_Pt2_log->GetXaxis()->SetTitle("log(Pt(W2))");
    
    
    TH1F *hVirtual_eta1=new TH1F("hVirtual_eta1","",50,-5,5);
    hVirtual_eta1->GetXaxis()->SetTitle("#eta W1");
    TH1F *hVirtual_eta2=new TH1F("hVirtual_eta2","",50,-5,5);
    hVirtual_eta2->GetXaxis()->SetTitle("#eta W2");
    
    TH1F *hVirtual_phi1=new TH1F("hVirtual_phi1","",50,-5,5);
    hVirtual_phi1->GetXaxis()->SetTitle("#phi W1");
    TH1F *hVirtual_phi2=new TH1F("hVirtual_phi2","",50,-5,5);
    hVirtual_phi2->GetXaxis()->SetTitle("#phi W2");
    
    TH1F *hParton_M1=new TH1F("hParton_M1","",50,-5,5);
    hParton_M1->GetXaxis()->SetTitle("M(parton1) [GeV]");
    TH1F *hParton_M2=new TH1F("hParton_M2","",50,-5,5);
    hParton_M2->GetXaxis()->SetTitle("M(parton2) [GeV]");
    
    TH1F *hTheta_HiggsJ1=new TH1F ("hTheta_HiggsJ1","",50,-1.1,1.1);
    hTheta_HiggsJ1->GetXaxis()->SetTitle("cos(#theta)(h,q1)");
    
    TH1F *hTheta_HiggsJ2=new TH1F ("hTheta_HiggsJ2","",50,-1.1,1.1);
    hTheta_HiggsJ2->GetXaxis()->SetTitle("cos(#theta)(h,q2)");
    
    TH1F *hthetastar_W1=new TH1F("hthetastar_W1","",50,-1.1,1.1);
    hthetastar_W1->GetXaxis()->SetTitle("cos(#theta*)(p1,q1)");
    TH1F *hthetastar_W2=new TH1F("hthetastar_W2","",50,-1.1,1.1);
    hthetastar_W2->GetXaxis()->SetTitle("cos(#theta*)(p2,q2)");
    
    TH1F *hthetastar_W2toHW1= new TH1F("hthetastar_W2toHW1","",50,-1.1,1.1);
    hthetastar_W2toHW1->GetXaxis()->SetTitle("cos(#theta*)(h,w1)");
    TH1F *hthetastar_W1toHW2=new TH1F("hthetastar_W1toHW2","",50,-1.1,1.1);
    hthetastar_W1toHW2->GetXaxis()->SetTitle("cos(#theta*)(h,w2)");
    TH1F *hthetastar_HtoWW=new TH1F("hthetastar_HtoWW","",50,-1.1,1.1);
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
    
    TH1F *hgen_mass=new TH1F ("hgen_mass","",100,80,180);  
    hgen_mass->GetXaxis()->SetTitle("GenLeptons Mass [GeV]");
    
    TH1F *hgenJetMass=new TH1F ("hgenJetMass","",100,0,3000);  
    hgenJetMass->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1F *hgenJetMass_cut06=new TH1F ("hgenJetMass_cut06","",100,0,3000);  
    hgenJetMass_cut06->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1F *hgenJetMass_cut08=new TH1F ("hgenJetMass_cut08","",100,0,3000);  
    hgenJetMass_cut08->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1F *hgenJetMass_cut10=new TH1F ("hgenJetMass_cut10","",100,0,3000);  
    hgenJetMass_cut10->GetXaxis()->SetTitle("GenJet Mass [GeV]");
    TH1F *hgenJetMass_cut12=new TH1F ("hgenJetMass_cut12","",100,0,3000);  
    hgenJetMass_cut12->GetXaxis()->SetTitle("GenJet Mass [GeV]");

    
    TH1F *hgenJetMass_matched=new TH1F ("hgenJetMass_matched","",100,0,500);  
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
        hMll_deltaM->GetXaxis()->SetTitle("m(ll) (GeV)");
        hMll_deltaM->GetYaxis()->SetTitle("#DeltaM (GeV)");

        std::vector<TH2F*> histo2D_vector;
        std::vector<std::string> variablesName_in_2D_plot = {"Zll_mass", "deltaM", "deltaR1", "deltaR2"/*, "Xparton1Log", "Xparton2Log"*/, "RpT", "zStar", "BDToutput", "genJetMassLeading", "hgenJetMassMatched", "Mqq" };

//         std::vector<float> limitDown  = {115,   0.,     -6,     -6,     0.,     -8,     0,      0,      0};
//         std::vector<float> limitUp    = {135,   0.1,    0.,     0.,     0.5,    3.,     1.6,    1500,    1500};
//         std::vector<int> binNumber  =   {12,    5,      6,      6,      10,     5,      16,     100,     100};
        
        std::vector<float> limitDown  = {115,   0.,     0,     0,     0.,     -8,     0,      0,      0,      0};
        std::vector<float> limitUp    = {135,   0.1,    6,     6,     0.5,    3.,     1.6,    1500,    1500,    1500};
        std::vector<int> binNumber  =   {12,    5,      20,      20,      10,     5,      16,     100,     100,     100};
        
        
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
//         TH1F* histArray[numArray] = { hMqq, hEtaQQ,hHTsoft,hSoft_n2,hSoft_n5,hSoft_n10,hHTsoftEWK,hSoft_n2EWK,hSoft_n5EWK,hSoft_n10EWK,hHTsoftEWK_bdt,hSoft_n2EWK_bdt,hSoft_n5EWK_bdt, hSoft_n10EWK_bdt,hnPVs, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hmet,   hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt,hV_mass, hqgl, hqgl2, hZll_mass, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRptHard, hEtaQQSum, hPhiZQ1, hZll_y, hZll_ystar, hZll_zstar, hMqq_log, hlheV_pt, hJet3_pt, hlheHT_log, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh,hbdt_atanh2 , hlepton1_iso03, hlepton2_iso03, hveto_jet3pt_nom, hveto_jet3pt_denom, hveto_ht_nom, hveto_ht_denom, hveto_softht_nom, hveto_softht_denom, hveto_softpt_nom, hveto_softpt_denom, hJet2q_phi, hJet1q_pffhi, hNAdJets, hNAdJets_bdt, hJet3_pt_bdt, hAdJetHT_bdt, hNAdJets_bdt2, hJet3_pt_bdt2, hAdJetHT_bdt2,hNAdJets_mjj1, hJet3_pt_mjj1, hAdJetHT_mjj1,hNAdJets_mjj2, hJet3_pt_mjj2, hAdJetHT_mjj2, hHTsoftEWK_bdt2,hSoft_n2EWK_bdt2,hSoft_n5EWK_bdt2, hSoft_n10EWK_bdt2,hHTsoftEWK_mjj1, hSoft_n2EWK_mjj1,hSoft_n5EWK_mjj1,hSoft_n10EWK_mjj1, hHTsoftEWK_mjj2,hSoft_n2EWK_mjj2,hSoft_n5EWK_mjj2,hSoft_n10EWK_mjj2 ,hJet1q_eta_bdt, hJet1q_eta_bdt2, hJet2q_eta_bdt, hJet2q_eta_bdt2, hsoftleadTrackPt, hsoftleadTrackEta, hAdJetHT, hJet3_eta , hJet3_pt_new , hJet3_eta_bdt, hJet3_eta_bdt2, hThetaStar, hMaxJetBTag};

        
        const int numArray= 185;//64+8 
        TH1F* histArray[numArray] = { hMqq,hMqq_cut06, hMqq_cut08, hMqq_cut10, hMqq_cut12, hgenJetMass_cut06, hgenJetMass_cut08, hgenJetMass_cut10, hgenJetMass_cut12, hEtaQQ,hSoft_n2,hSoft_n5,hSoft_n10,hHTsoftEWK,hSoft_n2EWK,hSoft_n5EWK,hSoft_n10EWK,hHTsoftEWK_bdt,hSoft_n2EWK_bdt,hSoft_n5EWK_bdt, hSoft_n10EWK_bdt,hnPVs, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hVtype, hVtypeSim, hmet, hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt, hqgl, hqgl2,hqglAtanh, hqgl2Atanh, hZll_mass, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRpt, hRptAtanh, hEtaQQSum, hPhiZQ1, hZll_y, hZll_ystar, hZll_zstar, hMqq_log, hlheV_pt, hlheNpNLO, hJet3_pt, hJet3_pt_log, hlheHT_log, hlheNj, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh,hbdt_atanh2 , hlepton1_iso03, hlepton2_iso03, hveto_jet3pt_nom, hveto_jet3pt_denom, hveto_ht_nom, hveto_ht_denom, hveto_softht_nom, hveto_softht_denom, hveto_softpt_nom, hveto_softpt_denom, hJet2q_phi, hJet1q_phi, hNAdJets, hNAdJets_bdt, hJet3_pt_bdt, hAdJetHT_bdt, hNAdJets_bdt2, hJet3_pt_bdt2, hAdJetHT_bdt2, hHTsoftEWK_bdt2,hSoft_n2EWK_bdt2,hSoft_n5EWK_bdt2, hPtSoftJets, hJet1q_eta_bdt, hJet1q_eta_bdt2, hJet2q_eta_bdt, hJet2q_eta_bdt2, hsoftleadTrackPt, hsoftleadTrackEta, hAdJetHT, hJet3_eta, hJet3_etaRatio, hJet3_pt_new , hJet3_eta_bdt, hJet3_eta_bdt2, hBDT_VBF, hBDT_VBF_atanh, hBDT_VBF_atanh_findBinning, hThetaStarJet, hThetaPlanes, hThetaStarJetAtanh, hThetaPlanesAtanh, hThetaStar, hThetaStarAbs, hMaxJetBTagCSV,hBDiscriminator_CSV, hBDiscriminator_CMVA, hweights_weighted,hweights,hdeltaMRel,hdeltaM, hnormalizedDistance_from_mH, hminAbsEta, hminAbsEta_GEN, hminAbsEta_cut06, hminAbsEta_cut08, hminAbsEta_cut10, hminAbsEta_cut12,  hHiggsSister1_Leading_eta,hHiggsSister2_Subleading_eta,hHiggsSister1_Leading_phi,hHiggsSister2_Subleading_phi,hHiggsSister1_Leading_R,hHiggsSister2_Subleading_R, hdeltaR1, hdeltaR2 ,  hMaxJetBTagCMVA,hTotalEnergy,hTotalEnergylog,hWWmass,hDiffmass,hpdgId,hgen_mass,hgenJetMass, hgenJetMass_matched, hEnergy_fraction_Parton2_log,hEnergy_fraction_Parton1_log, hmumujj_pt, hmumujj_ptLog, hEnergy_fraction_Parton1,hPz,hPzAbs, hPzAbsLog, hInvariant_MassLog, hInvariant_Mass, hthetastar_W2toHW1, hthetastar_W1toHW2,hthetastar_HtoWW, hEnergy_fraction_Parton2, hVirtual_Wmass1,hVirtual_Wmass2,hVirtual_Wmass1_log,hVirtual_Wmass2_log,hVirtual_Pt1,hVirtual_Pt2, hVirtual_Pt1_log, hVirtual_Pt2_log, hVirtual_eta1,hVirtual_eta2,hTheta_HiggsJ1,hTheta_HiggsJ2,hthetastar_W1,hthetastar_W2, hVirtual_phi1, hVirtual_phi2, hParton_M1,hParton_M2,hMaxSecondJetBTagCSV, hMaxSecondJetBTagCMVA, hSelectionCuts};
     
        
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


//	cout<<nentries<<endl;


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
	TH1F *hbdtUp[20];
	TH1F *hbdt_atanhUp[20];
	TH1F *hbdt_atanhUp_pdf[20];
	TH1F *hbdtDown[20];
	TH1F *hbdt_atanhDown[20];
	TH1F *hbdt_atanhDown_pdf[20];
    TString file_tag_hBDT_name = file_tag;
    file_tag_hBDT_name = Correct_fileTag (file_tag_hBDT_name);
	for (int i=0;i<Nsyst;i++){
		std::string histTitleUp;
		std::string histTitleDown;

		histTitleUp = ((std::string) hBDT_VBF->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] + "_Up";
		hbdtUp[i] = (TH1F*) hBDT_VBF->Clone(histTitleUp.c_str());
		histTitleUp = ((std::string) hBDT_VBF_atanh->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] +"_Up";
		hbdt_atanhUp[i] = (TH1F*) hBDT_VBF_atanh->Clone(histTitleUp.c_str());


//        hbdtUp[i]->Sumw2();
//        hbdt_atanhUp[i]->Sumw2();

		if (i!=0) {
			histTitleDown = ((std::string) hBDT_VBF->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] + "_Down";
			hbdtDown[i] = (TH1F*) hBDT_VBF->Clone(histTitleDown.c_str());
			histTitleDown = ((std::string) hBDT_VBF_atanh->GetName()) +"_"+ region +"_"+ file_tag_hBDT_name + "_" + uncertainty_name[i] +"_Down";
			hbdt_atanhDown[i] = (TH1F*) hBDT_VBF_atanh->Clone(histTitleDown.c_str());
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



        TFile fileMVA("mvaTree/main_tmva_tree_"+file_tag+"_v25"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+".root","recreate");

                
        TMVAstruct TMVA;
        TTree *treeMVA = new TTree("TMVAtree","TMVAtree");

        
        treeMVA->Branch("ll_mass",&TMVA.ll_mass,"ll_mass/F");
        treeMVA->Branch("diffMassWWH",&TMVA.diffMassWWH,"diffMassWWH/F");
	treeMVA->Branch("ll_pt",&TMVA.ll_pt,"ll_pt/F");
	treeMVA->Branch("ll_eta",&TMVA.ll_eta,"ll_eta/F");
	treeMVA->Branch("ll_ystar",&TMVA.ll_ystar,"ll_ystar/F");
        treeMVA->Branch("ll_zstar",&TMVA.ll_zstar,"ll_zstar/F");   
        
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
	treeMVA->Branch("evt",&evt,"evt/I");
	treeMVA->Branch("DeltaEtaQQ",&TMVA.DeltaEtaQQ,"DeltaEtaQQ/F");
	treeMVA->Branch("Jet2q_pt",&TMVA.Jet2q_pt,"Jet2q_pt/F");
	treeMVA->Branch("Jet1q_pt",&TMVA.Jet1q_pt,"Jet1q_pt/F");
        
        treeMVA->Branch("Jet1q_ptLog",&TMVA.Jet1q_ptLog,"Jet1q_ptLog/F");
        treeMVA->Branch("Jet2q_ptLog",&TMVA.Jet2q_ptLog,"Jet2q_ptLog/F");
        treeMVA->Branch("Jet3q_ptLog",&TMVA.Jet3q_ptLog,"Jet3q_ptLog/F");
        
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
        
        
        treeMVA->Branch("BDToutput",&TMVA.BDToutput,"BDToutput/F");
        
        treeMVA->Branch("weightMVA",  &TMVA.weightMVA,"weightMVA/F"); 
        treeMVA->Branch("genweight",  &TMVA.genweight,"genweight/F"); 
  

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


//         if  ((file_tag.CompareTo("TT")==0)) MVAcountMAX = 3940;  
//         if  ((file_tag.CompareTo("DYJetstoLL_madgraph")==0)) MVAcountMAX = 100000;
//         if  ((file_tag.CompareTo("DY0JetsToLL_M")==0)) MVAcountMAX = 1500;    //   713 events
//         if  ((file_tag.CompareTo("DY1JetsToLL_M")==0)) MVAcountMAX = 5700;
//         if  ((file_tag.CompareTo("DY2JetsToLL_M")==0)) MVAcountMAX = 4600;
//         if  ((file_tag.CompareTo("DY3JetsToLL_M")==0)) MVAcountMAX = 2550;
//         if  ((file_tag.CompareTo("DY4JetsToLL_M")==0)) MVAcountMAX = 1970;
//         if  ((file_tag.CompareTo("VBF_HToMuMu")==0)) MVAcountMAX = 30000;  



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
        reader->AddVariable("ll_mass",&TMVA.ll_mass);   
        reader->AddVariable("MqqLog",&TMVA.MqqLog);
        reader->AddVariable("Rpt",&TMVA.Rpt);
//         reader->AddVariable("W_mass_virtual1Log",&TMVA.W_mass_virtual1Log); 
        reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1Log); 
        reader->AddVariable("ll_zstar",&TMVA.ll_zstar);
        reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);   // this is float
        reader->AddVariable("ll_pt",&TMVA.ll_pt);
//         reader->AddVariable("W_mass_virtual2Log",&TMVA.W_mass_virtual2Log);
        reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2Log);
        
        
//         {"ll_mass","MqqLog", "Rpt","W_mass_virtual1Log","ll_zstar","softActivityEWK_njets5", "ll_pt","W_mass_virtual2Log"}
        
        
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


    
    
/*    reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
    reader->AddVariable("thetastarHtoWW",&TMVA.thetastarHtoWW);    
    reader->AddVariable("W_Pt_virtual1",&TMVA.W_Pt_virtual1);*/ 
   //reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1); 
   // reader->AddVariable("W_mass_virtual2",&TMVA.W_mass_virtual2); 
   // 
/*        
      reader->AddVariable("ll_mass",&TMVA.ll_mass);
       reader->AddVariable("Mqq",&TMVA.Mqq);
//    reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
//    reader->AddVariable("q1_eta",&TMVA.q1_eta);
 //   reader->AddVariable("ll_pt",&TMVA.ll_pt);
 
   // reader->AddVariable("met_pt",&TMVA.met_pt);
   
//    reader->AddVariable("qq_pt",&TMVA.qq_pt);
      reader->AddVariable("RptHard",&TMVA.RptHard);
      reader->AddVariable("thetastarW1",&TMVA.thetastarW1);
       //reader->AddVariable("theta1",&TMVA.theta1);
      reader->AddVariable("W_mass_virtual1",&TMVA.W_mass_virtual1);
     reader->AddVariable("W_Pt_virtual1",&TMVA.W_Pt_virtual1);
   
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
    reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
     // reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
    //reader->AddVariable("W_eta_virtual1",&TMVA.W_eta_virtual1);
    //reader->AddVariable("E_parton1",&TMVA.E_parton1);
     //reader->AddVariable("ll_pt",&TMVA.ll_pt);
    //reader->AddVariable("energytot",&TMVA.energytot);
   //reader->AddVariable("deltaMRel",&TMVA.deltaMRel);*/
   
//    reader->BookMVA("BDTG", "/gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkimOutputMVA/Classification_BDTG.weights.xml");
//     reader->BookMVA("BDTG", "BDTClassification/Classification_BDTG.weights.xml");
//     reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/training/CMSSW_8_0_28/src/training/TMVA-v4.2.0/test/weights/TMVAClassification_BDTG_nomoremuaxis2jet2q_v25/Classification_BDTG.weights.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/Classification_BDTG.weights_filerootJune_ll_mass_Mqq_DeltaEtaQQ_ll_eta_RptHard_EWKHTsoft_ll_pt_Jet2q_pt.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/trainingForRecoveringJune/Classification_BDTG.weights_JuneOption_ll_mass_Mqq_RptHard_DeltaEtaQQ_llPt_llEta_Jet2qPt_EWKHTsoft.xml");
//         reader->BookMVA("BDTG", "/afs/cern.ch/user/a/abonavit/private/tesi/CMSSW_8_0_28/src/code/BDTClassification/januaryTraining/Classification_BDTG.weights_80000Event_200Tree_4Deep_mll_Mqq_RptHard_llZstar_softN5_llPt_Wmass2_Wpt1.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/marchTraining/Classification_BDTG.weights_80000Event_100Tree_2Deep_mll_Mqq_RptHard_llZstar_softN5_Wmass1.xml");
        reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt_Wmass2.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_llZstar_softN5_llPt_qgl1qAtanh.xml");
//         reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/aprilTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_llZstar_softN5_llPt_qgl1qAtanh_Wmass1.xml");

//             reader->BookMVA("BDTG", "/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/BDTxml/mayTraining/Classification_BDTG.weights_35000Event_100Tree_2Deep_mll_MqqLog_Rpt_Wmass1_llZstar_softN5_llPt_Wmass2.xml");

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




// nentries = 2000;

 for (int entry=0; entry<nentries;++entry){
 //  for (int entry=0; entry<2000;++entry){

              
//      if ( (entry != 4787) && (entry != 8308) && (entry != 22761) && (entry != 32720) && (entry != 46473) && (entry != 55252) && (entry != 56360) && (entry != 57195) && (entry != 80033) && (entry != 87393) && (entry != 148917) && (entry != 154148) && (entry != 161992) && (entry != 164665) && (entry != 194110) && (entry != 229504) && (entry != 237120) && (entry != 249155) && (entry != 252616) && (entry != 258547) && (entry != 265599) && (entry != 273465) && (entry != 280313) && (entry != 302077) && (entry != 303664) && (entry != 308495) && (entry != 310651) && (entry != 310748) && (entry != 324889) && (entry != 328329) && (entry != 329909) && (entry != 355420) && (entry != 366617) && (entry !=   390667) && (entry !=   392749) && (entry !=   394375) && (entry !=   401849) && (entry !=   404989) && (entry !=   408357) && (entry !=   425786) && (entry !=   427628) && (entry !=   436848) && (entry !=   437986) && (entry !=   443537) && (entry !=   467547) && (entry !=   477715) && (entry !=   481352) && (entry !=   484553) && (entry !=   511362) && (entry !=   528665) && (entry !=   536729) && (entry !=   539025) && (entry !=   560850) && (entry !=   561803) && (entry !=   562517) && (entry !=   575788) && (entry !=   593305) && (entry !=   593800) && (entry !=   602564) && (entry !=   602824) && (entry !=   603028) && (entry !=   626380) && (entry !=   638724) && (entry !=   649698) && (entry !=   655687) && (entry !=   663670) && (entry !=   681660) && (entry !=   682601) && (entry !=   699267) && (entry !=   702234) && (entry !=   703302) && (entry !=   716756) && (entry !=   724288) && (entry !=   727431) && (entry !=   765586) && (entry !=   770406) && (entry !=   783386) && (entry !=   789267) && (entry !=   790730) && (entry !=   797479) && (entry !=   805647) && (entry !=   809129) && (entry !=   823961) && (entry !=   832809) && (entry !=   834304) && (entry !=   849385) && (entry !=   852840) && (entry !=   852944) && (entry !=   858622) && (entry !=   890201) && (entry !=   899894) && (entry !=   928382) && (entry !=   932415) && (entry !=   947356) && (entry !=   949146) && (entry !=   969778) && (entry !=   972886) && (entry !=   978168) && (entry !=   988717) && (entry !=   991104) && (entry !=   994428) && (entry !=   1004135) && (entry !=   1006420) && (entry !=   1019709) && (entry !=   1035251) && (entry !=   1050713) && (entry !=   1060715) && (entry !=   1061896) && (entry !=   1073326) && (entry !=   1073567) && (entry !=   1090780) && (entry !=   1091600) && (entry !=   1098146) && (entry !=   1122343) && (entry !=   1134366) && (entry !=   1139646) && (entry !=   1146350) && (entry !=   1154885) && (entry !=   1155204) && (entry !=   1157237) && (entry !=   1159402) && (entry !=   1168739) && (entry !=   1186210) && (entry !=   1223647) && (entry !=   1232887) && (entry !=   1244130) && (entry !=   1267114) && (entry !=   1299130) && (entry !=   1320667) && (entry !=   1327994) && (entry !=   1342187) && (entry !=   1347575) && (entry !=   1398908) && (entry !=   1406095) && (entry !=   1407029) && (entry !=   1409048) && (entry !=   1426614) && (entry !=   1431689) && (entry !=   1461306) && (entry !=   1466498)) continue;
//      
//          
//          
//          
// 
// std::cout << "entry       " << entry  << std::endl;


        if (entry%10000 == 0) std::cout << "Writing " << entry << "th event" << std::endl;
        tree_initial->GetEntry(entry);
        ++count1;

// if (VtypeSim  == -1) cout << "v_typeSim 0 " << v_typeSim << "  " << VtypeSim << endl;
//   cout << "entry   " << entry << endl;

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
	TMVA.Jet1q_leadTrackPt=0;
	TMVA.ll_ystar=0;
	TMVA.ll_zstar=0;
	TMVA.Rpt=0;
	TMVA.ll_pt=0;
        TMVA.ll_eta=0;
	TMVA.qgl_1q=0;
	TMVA.qgl_2q=0;
        TMVA.qgl_1qAtanh=0;
	TMVA.qgl_2qAtanh==0;
        
        
        TMVA.cosThetaPlaneAtanh =0;
        TMVA.absCosThetaStarJetAtanh =0;
        TMVA.X_parton1Log=0;
        TMVA.X_parton2Log=0;
        TMVA.W_mass_virtual1Log=0;
        TMVA.W_mass_virtual2Log=0;
        TMVA.W_Pt_virtual1Log=0;
        TMVA.W_Pt_virtual2Log=0;
                    
                    
    
        TMVA.softActivityEWK_njets2=0;
        TMVA.softActivityEWK_njets5=0;
	TMVA.softActivityEWK_njets10=0;
        TMVA.randomVariable=0;
        TMVA.xSection=0;
        TMVA.weightMVA=0;
        TMVA.genweight=0;
        TMVA.deltaM=0;
        
        TMVA.BDToutput=0;
        
        float maxBTagCSV = -0.2;
        float maxBTagCMVA = -1.2;
        float maxSecondBTagCSV = -0.2;
        float maxSecondBTagCMVA = -1.2;
        lheNj = (int) lheNj_f;
        

        for (int i=0;i<nJets;i++){
            if (whichJESWeight==1) if(Jet.JEC_corr_up[i] > 0 && Jet.JEC_corr[i] > 0)    {Jet.pt[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i];    Jet.mass[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i];}	
            if (whichJESWeight==2) if(Jet.JEC_corr_down[i] > 0 && Jet.JEC_corr[i] > 0)  {Jet.pt[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];  Jet.mass[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];}	
            if (whichJERWeight==1) if(Jet.JER_corr_up[i] > 0 && Jet.JER_corr[i] > 0)    {Jet.pt[i]*= Jet.JER_corr_up[i]/Jet.JER_corr[i];    Jet.mass[i]*= Jet.JER_corr_up[i]/Jet.JER_corr[i];}	
            if (whichJERWeight==2) if(Jet.JER_corr_down[i] > 0 && Jet.JER_corr[i] > 0)  {Jet.pt[i]*= Jet.JER_corr_down[i]/Jet.JER_corr[i];  Jet.mass[i]*= Jet.JER_corr_down[i]/Jet.JER_corr[i];}	

        }
        


        if (JSON!=1) {
            continue;
        }

	//	if ((file_tag.CompareTo("EWK_LLJJ")==0) &&(evt%2!=0)) continue;

// 		if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
	//	if (region.CompareTo("el")==0) if (!(v_type==1)) continue;



	
	

	
	
		if (data==1) PU=1.;
		else PU=puweight;
		genweight0 = genweight/TMath::Abs(genweight);
		genweight=genweight/TMath::Abs(genweight)*PU;   
// 		cout << "PU  " << PU << " \t   genweight     " << genweight << endl; 

//		genweight/=events_generated/xsec[file_tag]; 
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_top")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==1)  genweight*=LHE_weights_scale_wgt[4];
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)&& (file_tag.CompareTo("ST_tW_antitop")!=0)) if (whichQCDScaleWeight==2)  genweight*=LHE_weights_scale_wgt[5];

		
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
        /////////////////////////////////////////////////////////////////////////// FILTER GENJET INVARIANT MASS ////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<TLorentzVector> genLeptons;
        
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
        
        
           
        std::vector<TLorentzVector> genJet;
        
        for (int n=0; n<nGenJet; n++) {
            TLorentzVector genJet_tmp;
            genJet_tmp.SetPtEtaPhiM(GenJet_pt[n],GenJet_eta[n], GenJet_phi[n], GenJet_mass[n]);
            genJet.push_back(genJet_tmp);
        }
        
        
//         float genJetMass = VBFFilter(genLeptons, genJet);
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
                            
                            
                
        for (int i=0;i<nJets;i++){

            TLorentzVector jet0;
            if (!((Jet.id[i]>2)&&(Jet.puId[i]>0))) continue;
            
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

//             if (maxBTagCMVA < -1) maxBTagCMVA = Jet.btagCMVA[i];
//             if (maxBTagCMVA > -1 && maxSecondBTagCMVA < -1) maxSecondBTagCMVA = Jet.btagCMVA[i];
//             
//             if (maxBTagCSV < 0) maxBTagCSV = Jet.btagCSV[i];
//             if (maxBTagCSV > 0 && maxSecondBTagCSV < 0) maxSecondBTagCSV = Jet.btagCSV[i];
                
            bool condition = false;
            if (good_jets < 2) condition = true;
            else if (jet0.Eta() < max(jets_pv[0].Eta(), jets_pv[1].Eta()) && jet0.Eta() > min(jets_pv[0].Eta(), jets_pv[1].Eta())) condition = true;  

            if (condition) {
                jets_pv.push_back(jet0);
                jets_indices.push_back(i);
                good_jets++;
            }
        }
 
 
        

        
        if (good_jets<2) continue;
        else  ++count2;    

        

//         if (maxSecondBTagCMVA > maxBTagCMVA) {float tmpBTag = maxSecondBTagCMVA; maxSecondBTagCMVA = maxBTagCMVA; maxBTagCMVA = tmpBTag;}
//         if (maxSecondBTagCSV > maxBTagCSV) {float tmpBTag = maxSecondBTagCSV; maxSecondBTagCSV = maxBTagCSV; maxBTagCSV = tmpBTag;}

            
//         if(Jet.qgl[jets_indices[0]]<-0.5) std::cout << good_jets << " \t " << entry << " \t jets_indices[0]:  " << jets_indices[0] << " \t " << Jet.qgl[jets_indices[0]] << std::endl;
//         if(Jet.qgl[jets_indices[1]]<-0.5) std::cout << good_jets << " \t " << entry << " \t jets_indices[1]:  " << jets_indices[1] << " \t " << Jet.qgl[jets_indices[1]] << std::endl;

        if(Jet.qgl[jets_indices[0]]<-0.5) Jet.qgl[jets_indices[0]] = 0.0000001;
        if(Jet.qgl[jets_indices[1]]<-0.5) Jet.qgl[jets_indices[1]] = 0.0000001;
        
        
        
        if (maxBTagCMVA>-0.853213 && maxBTagCMVA <-0.853211) maxBTagCMVA = -0.999;
        if (maxSecondBTagCMVA>-0.853213 && maxSecondBTagCMVA <-0.853211) maxSecondBTagCMVA = -0.999;
        Int_t EWKnsoft2_toSubtract = 0;
        Int_t EWKnsoft5_toSubtract = 0;
        Int_t EWKnsoft10_toSubtract = 0;
        float EWKHTsoft_toSubtract = 0.;

//      SOFT EWK 3rd JET IN THE RAPIDITY GAP
         for (int i=0;i<Jet.nsoftActivityEWKJets;i++){
             TLorentzVector jet0;
             bool found_softLeadingJet = false;
             jet0.SetPtEtaPhiM(Jet.EWKsoft_pt[i],Jet.EWKsoft_eta[i],Jet.EWKsoft_phi[i],Jet.EWKsoft_mass[i]);
             if (jet0.Eta() > max(jets_pv[0].Eta(), jets_pv[1].Eta()) || jet0.Eta() < min(jets_pv[0].Eta(), jets_pv[1].Eta())) {
                 if (Jet.EWKsoft_pt[i] > 1)  EWKHTsoft_toSubtract += Jet.EWKsoft_pt[i];
                 if (Jet.EWKsoft_pt[i] > 2)  EWKnsoft2_toSubtract++;
                 if (Jet.EWKsoft_pt[i] > 5)  EWKnsoft5_toSubtract++;
                 if (Jet.EWKsoft_pt[i] > 10)  EWKnsoft10_toSubtract++;
             }
             else 
                 if (!found_softLeadingJet) {
                     Jet.softLeadingJet_pt = Jet.EWKsoft_pt[i];
                     found_softLeadingJet = true;
                 }
         }

        Qjet1 = jets_pv[0];
        Qjet2 = jets_pv[1];
		if (good_jets>=3) {
			jet3_pt=jets_pv[2].Pt();
			jet3_eta=jets_pv[2].Eta();
		}
		qq=Qjet1+Qjet2;
		Float_t Mqq = qq.M();
		Float_t qq_pt = qq.Pt();
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));

//////////////////leptons////////////////
		TLorentzVector lepton1;
		TLorentzVector lepton2;
		TLorentzVector Zll;

		
                
                



                
		int idx_1stLepton = 0;
		int idx_2ndLepton = 0;
		int count_l=0;


                if (region.CompareTo("mu")==0) {
		    count_l=0;
		    idx_1stLepton = 0;
		    idx_2ndLepton = 0;
		    for (int i=0; i<nselLeptons;i++ ){
			    if (!((selLeptons_looseIdPOG[i]>0) && (selLeptons_relIso04[i]<0.25) && (TMath::Abs(selLeptons_pdgId[i])==13 )) ) continue;
                            if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
			    if (count_l==1) {
				    idx_2ndLepton=i;
				    lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				    count_l++;
				    break;
			    }
			    if (count_l==0) {
				    idx_1stLepton=i;
				    lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				    count_l++;
			    }
		    }
		}


        
        
//        TLorentzVector Hll = lepton1+lepton2;
//        if  ( Hll.M()<100 )  continue;
//        else  ++count14;
//        std::cout << "Hll.M()   " << Hll.M()<< std::endl;


/*
		lepton1.SetPtEtaPhiM(vLeptons_pt[0], vLeptons_eta[0], vLeptons_phi[0], vLeptons_mass[0]);	
		int idx_2ndLepton = 0;
		for (int i=1; i<nvLeptons;i++ ){
			if (vLeptons_charge[0]*vLeptons_charge[i] < 0) {
				idx_2ndLepton=i;
				break;
			}
		}
		lepton2.SetPtEtaPhiM(vLeptons_pt[idx_2ndLepton], vLeptons_eta[idx_2ndLepton], vLeptons_phi[idx_2ndLepton], vLeptons_mass[idx_2ndLepton]);

		float qter1 = 1.0;
		float qter2 = 1.0;
		float mu_correction1 = 1.0;
		float mu_correction2 = 1.0;
		if (data!=1) 	if (region.CompareTo("mu")==0) {
			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
			}
		if (data==1) 	if (region.CompareTo("mu")==0){
			rmcor->momcor_data(lepton1, vLeptons_charge[0],  0, qter1);
			rmcor->momcor_data(lepton2, vLeptons_charge[idx_2ndLepton], 0, qter2);
			}
*/
///////////////muon corrections 2016 calibration////////////////
// RoccoR  *rc = new RoccoR("/scratch/mandorli/Hmumu/CMSSW_8_0_25/src/giulioMandorli/2016/rcdata.2016.v3/");



		if (region.CompareTo("mu")==0) {
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
			//	genweight*= eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}

		}


                
		string file_tag_str = file_tag.Data();  //ALREADY DECLARED: CHECK!
// 		if  (file_tag_str.find/*("DYJetstoLL"*/)!=std::string::npos)  genweight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 


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


                
        qgl_weight=1.;
        if (!( (data==1)|| (Jet.partonFlavour[jets_indices[1]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[1]])>=2) || (Jet.qgl[jets_indices[1]] < 0)) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) < 4 ) qgl_weight=func_qgl_q->Eval(Jet.qgl[jets_indices[1]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) ==21 ) qgl_weight=func_qgl_g->Eval(Jet.qgl[jets_indices[1]]);
		}
		if (qgl_weight!=1.) apply_qgl+=1;
		if (data!=1) genweight*=qgl_norm[file_tag];

		
		genweight*=qgl_weight;
	
// 	        cout  << "4 \t\t " << genweight << "  *  " << xsec[file_tag] << "/" << events_generated << endl;

		
//         std::cout << entry << "\t lepton1  " << lepton1.Pt() << " \t  " << lepton1.Phi() << " \t  " << lepton1.Eta() << " \t  " << lepton1.M() << " \t  "  << std::endl;
//         std::cout << entry << "\t lepton2  " << lepton2.Pt() << " \t  " << lepton2.Phi() << " \t  " << lepton2.Eta() << " \t  " << lepton2.M() << " \t  "  << std::endl;

		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();
		
//         std::cout << "Zll_mass  " << Zll_mass << std::endl;


//         if (genweight<0.)  continue;
            
//             if(genweight > 0.) genweight = 1.;
//             else genweight = -1.;
                
        genweight/=events_generated/xsec[file_tag]; 
                
                
        hSelectionCuts->Fill(0., genweight);
/*        if ((file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0) || (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0) || (file_tag.CompareTo("DY0JetsToLL_M")==0) || (file_tag.CompareTo("DY1JetsToLL_M")==0) || (file_tag.CompareTo("DY2JetsToLL_M")==0) || (file_tag.CompareTo("DY3JetsToLL_M")==0) || (file_tag.CompareTo("DY4JetsToLL_M")==0) || (file_tag.CompareTo("DYJetstoLL_amc_M-50")==0)) 
        if((gen_mass < 105) || (gen_mass > 160)) continue; */ //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        hSelectionCuts->Fill(1, genweight);
        

  
                        
// 	if (Zll_mass < 105 ) continue;
//         if (Zll_mass < 110 ) continue;
        if (Zll_mass < 115 ) continue;      //----------------------------------------------------------------------------------------------------
        else  ++count12;
        hSelectionCuts->Fill(2, genweight);
                		
//         if (Zll_mass > 145 ) continue;
        if (Zll_mass > 135 ) continue;      //----------------------------------------------------------------------------------------------------
        hSelectionCuts->Fill(3, genweight);
        
        
        
        
        
//         if (Zll_mass > 118 && Zll_mass < 128 ) continue;
        
        
        
                
        if  (region.CompareTo("mu")==0) if (!((HLT_IsoMu24==1) || (HLT_IsoTkMu24==1)  )) continue; 
        else  ++count4;
        hSelectionCuts->Fill(4, genweight);  //Trigger only      
                
                
        if (good_jets<2) continue;
        else  ++count2;    
        hSelectionCuts->Fill(5, genweight);
        
        if (qqDeltaEta < 2.5) continue;
        hSelectionCuts->Fill(6, genweight);
        
        if (count_l<2)  continue;
        else  ++count3;
        hSelectionCuts->Fill(7, genweight);
        if ((selLeptons_charge[idx_1stLepton]*selLeptons_charge[idx_2ndLepton]) >0) continue;                                   // redundant
        hSelectionCuts->Fill(8, genweight);
        

          
		cut_flow[1]+=genweight;
		if (Qjet1.Pt() < 35) continue;
                else  ++count5;
                hSelectionCuts->Fill(9, genweight);
                
		cut_flow[2]+=genweight;
 		if (Qjet2.Pt() < 25) continue;
// 		if (Qjet2.Pt() < 22) continue;
                else  ++count6;
                hSelectionCuts->Fill(10, genweight);
                
		cut_flow[3]+=genweight;
// 		if (Mqq<200) continue;
		if (Mqq<250) continue;      //----------------------------------------------------------------------------------------------------
                else  ++count7;
                hSelectionCuts->Fill(11, genweight);
                
		cut_flow[4]+=genweight;

	 	if (region.CompareTo("mu")==0) {
			if (lepton1.Pt()<30) continue;
                        else  ++count8;
                        hSelectionCuts->Fill(12, genweight);
			if (lepton2.Pt()<10) continue;
                        else  ++count9;
                        hSelectionCuts->Fill(13, genweight);
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
	 	if (region.CompareTo("mu")==0) {
			if (TMath::Abs(lepton1.Eta())>2.4) continue;
                        else  ++count10;
                        hSelectionCuts->Fill(14, genweight);
			if (TMath::Abs(lepton2.Eta())>2.4) continue;
                        else  ++count11;
                        hSelectionCuts->Fill(15, genweight);
		}	


                if (maxBTagCMVA > 0.95) continue;
                if (maxBTagCMVA > 0.8) continue;      //----------------------------------------------------------------------------------------------------
                else  ++count13;
                hSelectionCuts->Fill(16, genweight);


		presel+=genweight;

		presel_vtype[(int)(v_type+1)]+=genweight;

		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		float Zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t Rpt = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 


//ADDITIONAL CUTS
		if (Zll_zstar > 2.5 ) continue;
                hSelectionCuts->Fill(17, genweight);
                if (Rpt>0.4) continue;      //----------------------------------------------------------------------------------------------------

                

//                 float genJetMass = VBFFilter(genLeptons, genJet);
            float genJetMass_leading = 0;
//             std::vector<TLorentzVector> genJetsWithoutLeptonsP4= VBFFilter(genLeptons, genJet);
            std::vector<TLorentzVector> recoLepton;
            recoLepton.push_back(lepton1);
            recoLepton.push_back(lepton2);
            std::vector<TLorentzVector> genJetsWithoutLeptonsP4= VBFFilter(recoLepton, genJet);
            if (genJetsWithoutLeptonsP4.size() > 1) genJetMass_leading = (genJetsWithoutLeptonsP4[0] + genJetsWithoutLeptonsP4[1]).M();
                 
            float minAbsEta_GEN = 10;
            if (genJetsWithoutLeptonsP4.size() > 1) minAbsEta_GEN = min(abs(genJetsWithoutLeptonsP4[0].Eta()), abs(genJetsWithoutLeptonsP4[1].Eta()));
            
            float minAbsEta = 10;
            minAbsEta = min(abs(Qjet1.Eta()), abs(Qjet2.Eta()));
            
            
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
            
            hHiggsSister1_Leading_R->Fill(HiggsSister1_Leading_R, genweight);
            hHiggsSister2_Subleading_R->Fill(HiggsSister2_Subleading_R, genweight);
            
            
            
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
				hlheNj->Fill(lheNj ,genweight);
		        hgen_mass->Fill(gen_mass,genweight);
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
/*/
           if(file_tag.CompareTo("DYJetstoLL_amc_0J")==0){
            if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
                if((M_ll_cutbin==1)&&(lheNj==0))
           genweight=genweight*(0.042/4.779);
           }
           }
           if(file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0){
            if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
                if((M_ll_cutbin==1)&&(lheNj==0))
           genweight=genweight*(0.042/4.779);
           else genweight=0;
           }
           else genweight=0;
           }

/*/

//---------------------STITCHING------------------------------------------------------------

            


           float ratio = 1.;/*/
            if(file_tag.CompareTo("DYJetstoLL_amc_Filter105")==0){
                if(abs(GenLep_pdgId[0])==11 && abs(GenLep_pdgId[1])==11 ){
                        if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==2))  ratio = 1.;
                    }
                if(abs(GenLep_pdgId[0])==13 && abs(GenLep_pdgId[1])==13 ){
                        if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 0.00290597/0.00290597;
                        if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 0.0172459/0.0171639;
                        if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 0.0196249/0.0192204;
                        if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 0.0429222/4.77991;
                        if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 0.189544/4.31358;
                        if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 0.211394/1.93457;
                        if((M_ll_cutbin==2)&&(lheNj==0))  ratio =  0.000469661/0.000469661;
                        if((M_ll_cutbin==2)&&(lheNj==1))  ratio =  0.00111306/0.00111306;
                        if((M_ll_cutbin==2)&&(lheNj==2))  ratio =  0.000978254/0.000918971;
                    }
                if(abs(GenLep_pdgId[0])==15 && abs(GenLep_pdgId[1])==15 ){
                        if((M_ll_cutbin==0)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==0)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==0)&&(lheNj==2))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==1)&&(lheNj==2))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==0))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==1))  ratio = 1.;
                        if((M_ll_cutbin==2)&&(lheNj==2))  ratio = 1.;
                    }
            }/*/

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
                        TMVA.jets12=jets_pv[0].Pt() + jets_pv[1].Pt();
                        
                        TMVA.deltaMRel=deltaMRel;
                    TMVA.deltaM=deltaM;
                    TMVA.normalizedDistance_from_mH=normalizedDistance_from_mH;
                    TMVA.ll_mass = Zll_mass;
                    TMVA.diffMassWWH=WWmass-Zll_mass;
//                     TMVA.ll_pt = log(Zll_pt);
                    TMVA.q1_eta = jets_pv[0].Eta();
                    TMVA.met_pt = met_pt;
                    TMVA.EWKHTsoft = Jet.EWKHTsoft - EWKHTsoft_toSubtract;
                    TMVA.btagCMVA = maxBTagCMVA;
                    TMVA.btagCSV= maxBTagCSV;
                    TMVA.btagCMVA_second = maxSecondBTagCMVA;
                    TMVA.btagCSV_second= maxSecondBTagCSV;
                    TMVA.softLeadingJet_pt = Jet.softLeadingJet_pt;

                    TMVA.cosThetaStar=cosThetaStar;
                    TMVA.cosThetaStarAbs=abs(cosThetaStar);
                    TMVA.cosThetaPlane = cosThetaPlane;
                    TMVA.absCosThetaStarJet = absCosThetaStarJet;

                    
                    TMVA.Mqq = Mqq;
                    TMVA.MqqLog = log(Mqq);
                    
                    TMVA.Jet3_eta=jet3_eta;
                    
                    TMVA.Jet1q_pt = Qjet1.Pt();
                    TMVA.Jet2q_pt = Qjet2.Pt();
                    TMVA.Jet3q_pt = jet3_pt;
                    TMVA.Jet3q_ptLog = log(jet3_pt);
                    TMVA.Jet2q_ptLog = log(Qjet2.Pt());
                    TMVA.Jet1q_ptLog = log(Qjet1.Pt());
                    //TMVA.PhiZQ1 =PhiZQ1;

                    TMVA.ll_zstar =  (Zll_zstar > 0.000001) ? log(Zll_zstar) : -13.8155105579642736;
                    TMVA.ll_ystar =  Zll_ystar;
    // 		if (TMath::Exp((-1)*Jet_axis2[jets_indices[0]]) < 0.2)  TMVA.axis2_jet1 = TMath::Exp((-1)*Jet_axis2[jets_indices[0]]);
    // 		if (TMath::Exp((-1)*Jet_axis2[jets_indices[1]]) < 0.2) TMVA.axis2_jet2 = TMath::Exp((-1)*Jet_axis2[jets_indices[1]]);
    // 		TMVA.Jet1q_leadTrackPt = Jet_leadTrackPt[jets_indices[0]];
                    TMVA.Rpt = Rpt;
                    TMVA.DeltaEtaQQ = qqDeltaEta;
    // 		TMVA.Jet2q_leadTrackPt = Jet_leadTrackPt[jets_indices[1]];
                    TMVA.qq_pt = log(qq_pt);
                    TMVA.ll_pt = Zll_pt;
                    TMVA.ll_eta = Zll_eta;
                    
                    TMVA.Jet3_pt = jet3_pt;
                    TMVA.qgl_1q = Jet.qgl[jets_indices[0]];
                    TMVA.qgl_2q = Jet.qgl[jets_indices[1]];
                    
                    if(Jet.qgl[jets_indices[0]]<-0.5) std::cout << entry << " \t jets_indices[0]:  " << jets_indices[0] << " \t " << Jet.qgl[jets_indices[0]] << std::endl;
                    if(Jet.qgl[jets_indices[1]]<-0.5) std::cout << entry << " \t jets_indices[1]:  " << jets_indices[1] << " \t " << Jet.qgl[jets_indices[1]] << std::endl;
                    
                    
                    TMVA.softActivityEWK_njets2 = Jet.EWKnsoft2 - EWKnsoft2_toSubtract;
                    TMVA.softActivityEWK_njets5 = Jet.EWKnsoft5 - EWKnsoft5_toSubtract;
                    TMVA.softActivityEWK_njets10 = Jet.EWKnsoft10 - EWKnsoft10_toSubtract;
//                    TMVA.weightMVA = genweight; 
                   TMVA.genweight = genweight; 
//                    TMVA.weightMVA = weightMVA;
                   TMVA.weightMVA = 1.;

                   
//                     TMVA.randomVariable = gRandom->Rndm();
                    TMVA.xSection = xsec[file_tag];
                    MVAcount +=1;
//                     if ( MVAtree_to_fill && MVAcount < MVAcountMAX) treeMVA->Fill();   

                    
                    if  ((file_tag.CompareTo("DY0JetsToLL_M")==0) && lheNj == 0) {
                        TMVA.weightMVA = 2.1;
//                        std::cout << "weightMVA: " << TMVA.weightMVA << std::endl;
                    }
                    
                    temp_softActivityEWK_njets5 = (float) TMVA.softActivityEWK_njets5;
                    if (plotOutput)  {
//                     if (0)  {
//                        std::cout << "plotOutput  "  << std::endl;
                        BDT_VBF = reader->EvaluateMVA("BDTG");
//                         BDT_VBF = reader->EvaluateMVA("MLP");
                     
//                         if(atanh((BDT_VBF+1.)/2.)<0.9) continue;
                        
                        TMVA.BDToutput = BDT_VBF;
                        
                        hBDT_VBF->Fill(BDT_VBF-0.01,genweight);
                        hBDT_VBF_atanh->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        hBDT_VBF_atanh_findBinning->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        
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
                        
                        
                    }
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
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
                    
                    if ( MVAtree_to_fill && MVAcount < MVAcountMAX) treeMVA->Fill();   
                    
                    
                    
                    
                    
                    

//                     std::cout << entry << "\t1:  " << reader->EvaluateMVA("BDTG") << " \t " << TMVA.ll_mass << " \t " << TMVA.Mqq << " \t " << TMVA.RptHard << " \t " << TMVA.DeltaEtaQQ << " \t " << TMVA.ll_pt << " \t " << TMVA.ll_eta << " \t " << TMVA.Jet2q_pt << " \t " << TMVA.EWKHTsoft << std::endl;
//                     
//                     TMVA.ll_mass = 126.13;
//                     TMVA.Mqq = 1879.31;
//                     TMVA.RptHard = 0.0314452;
//                     TMVA.DeltaEtaQQ = 6.91397;
//                     
//                     TMVA.ll_pt = 126.873;
//                     TMVA.ll_eta = 0.124247;
// 
//                     TMVA.qgl_2q = 49.9826;
//                     TMVA.EWKHTsoft=3.25488;
//                     
//                     
//                     
//                     std::cout << entry << "\t2:  " << atanh((reader->EvaluateMVA("BDTG")+1.)/2.) << " \t " << TMVA.ll_mass << " \t " << TMVA.Mqq << " \t " << TMVA.RptHard << "  \t " << TMVA.DeltaEtaQQ << " \t " << TMVA.ll_pt << " \t " << TMVA.ll_eta << " \t " << TMVA.Jet2q_pt << " \t " << TMVA.EWKHTsoft << std::endl;
                    
                    
//                     lwt::ValueMap nnout;
//                     float vars[18] = { TMVA.Mqq, TMVA.DeltaEtaQQ, TMVA.q1_eta, TMVA.ll_pt, TMVA.ll_mass, TMVA.met_pt, TMVA.EWKHTsoft, TMVA.qq_pt, TMVA.RptHard, (float) TMVA.softActivityEWK_njets5, TMVA.btagCMVA, TMVA.btagCSV, TMVA.qgl_1q, TMVA.qgl_2q, TMVA.cosThetaStar, TMVA.cosThetaStarAbs};
                    float vars[18] = { TMVA.Mqq, TMVA.DeltaEtaQQ, TMVA.q1_eta, TMVA.ll_pt, TMVA.ll_mass, TMVA.met_pt, TMVA.EWKHTsoft, TMVA.qq_pt, TMVA.Rpt, (float) TMVA.softActivityEWK_njets5, TMVA.btagCSV, TMVA.qgl_1q, TMVA.qgl_2q, TMVA.ll_zstar, TMVA.ll_ystar, TMVA.Jet1q_pt, TMVA.Jet2q_pt, TMVA.ll_eta};
                    
                    
                    

                    
//                     int count = 0;
//                     inputs.clear();
//                     for(auto& input : config_TrackPair.inputs) {	
//                         std::cout << input.name <<"  "<<  vars[count] <<std::endl;                
//                         inputs[input.name] = vars[count]; 
//                         count++;
//                     }
//         float vars[49] = {0};

                
                   if (plotOutput_NN)  {

                       
                        BDT_VBF = nnResult.result(vars) - 0.0000001;
//                         BDT_VBF = BDT_VBF + 4.*reader->EvaluateMVA("BDTG");
//                         std::cout << "BDT_VBF:   " << BDT_VBF << std::endl;

                        
                        float  atanhBDT = atanh(BDT_VBF) + 4*atanh((reader->EvaluateMVA("BDTG")+1.)/2.);
//                         hBDT_VBF->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        hBDT_VBF->Fill(atanhBDT,genweight);
//                         hBDT_VBF_atanh->Fill(atanh((BDT_VBF+1.)/2.),genweight);
                        hBDT_VBF_atanh->Fill(atanhBDT,genweight);
                        hBDT_VBF_atanh_findBinning->Fill(atanhBDT,genweight);
                    }
        
        
        
                }
                
                
                
                

            
            
            
            
/////////////////////////////////////////////// BTAG FILIPPO ////////////////////////////////////////////////////////////////////////////////
            hBDiscriminator_CSV->Fill(Jet.btagCSV[jets_indices[0]], genweight);
            hBDiscriminator_CSV->Fill(Jet.btagCSV[jets_indices[1]], genweight);
            hBDiscriminator_CMVA->Fill(Jet.btagCMVA[jets_indices[0]], genweight);
            hBDiscriminator_CMVA->Fill(Jet.btagCMVA[jets_indices[1]], genweight);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       
            
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
           hVirtual_eta1->Fill(W_eta_virtual1,genweight);
           hVirtual_eta2->Fill(W_eta_virtual2,genweight);
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
		   	hZll_mass->Fill(Zll_mass,genweight);
		   	hZll_pt->Fill(Zll_pt,genweight);
		   	hZll_phi->Fill(Zll_phi,genweight);
		   	hZll_eta->Fill(Zll_eta,genweight);
// 			   hHTsoft->Fill(Jet.HTsoft,genweight);
// 			   hSoft_n2->Fill(Jet.nsoft2, genweight);
// 			   hSoft_n5->Fill(Jet.nsoft5, genweight);
// 			   hSoft_n10->Fill(Jet.nsoft10, genweight);
			   hHTsoftEWK->Fill(Jet.EWKHTsoft - EWKHTsoft_toSubtract,genweight);
			   hSoft_n2EWK->Fill(Jet.EWKnsoft2 - EWKnsoft2_toSubtract, genweight);
			   hSoft_n5EWK->Fill(Jet.EWKnsoft5 - EWKnsoft5_toSubtract, genweight);
			   hSoft_n10EWK->Fill(Jet.EWKnsoft10 - EWKnsoft10_toSubtract, genweight);
                hPtSoftJets->Fill(Jet.softLeadingJet_pt, genweight);
				hnPVs->Fill(nPVs,genweight);
				hqgl->Fill(Jet.qgl[jets_indices[0]],genweight);
				hqgl2->Fill(Jet.qgl[jets_indices[1]],genweight);
                                hqglAtanh->Fill(atanh(2.*Jet.qgl[jets_indices[0]]-1.),genweight);
                                hqgl2Atanh->Fill(atanh(2.*Jet.qgl[jets_indices[1]]-1.),genweight);
                        
				hJet1q_pt->Fill(jets_pv[0].Pt(),genweight);
				hJet1q_eta->Fill(jets_pv[0].Eta(),genweight);
				hJet1q_ptd->Fill(Jet.ptd[jets_indices[0]],genweight);
				hJet1q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[0]]),genweight);
				hJet1q_mult->Fill(Jet.mult[jets_indices[0]],genweight);
				hJet1q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[0]],genweight);
				hJet1q_phi->Fill(jets_pv[0].Phi(),genweight);
				hJet2q_pt->Fill(jets_pv[1].Pt(),genweight);
				hJet2q_eta->Fill(jets_pv[1].Eta(),genweight);
				hJet2q_ptd->Fill(Jet.ptd[jets_indices[1]],genweight);
				hJet2q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[1]]),genweight);
				hJet2q_mult->Fill(Jet.mult[jets_indices[1]],genweight);
				hJet2q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[1]],genweight);
				hJet2q_phi->Fill(jets_pv[1].Phi(),genweight);
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
				hZll_y->Fill(Zll.Rapidity(),genweight);
				hZll_ystar->Fill(Zll_ystar   ,genweight);
				hZll_zstar->Fill(log(Zll_zstar),genweight);
				hlheV_pt->Fill(lheV_pt  ,genweight);
				hJet3_pt->Fill(jet3_pt ,genweight);	
				if (good_jets>=3) hJet3_pt_log->Fill(log(jets_pv[2].Pt()) ,genweight);	
				if (good_jets>=3) hJet3_eta->Fill(jet3_eta ,genweight);	
                                if (good_jets>=3) hJet3_etaRatio->Fill((2.*jet3_eta + jets_pv[0].Eta() + jets_pv[1].Eta()/abs(jets_pv   [0].Eta() - jets_pv[1].Eta())) ,genweight);	
				if (good_jets>=3) hJet3_pt_new->Fill(jets_pv[2].Pt(),genweight);
				if (good_jets==2) hJet3_pt_new->Fill(10.,genweight);
				float AdJetHT = 0;
				if (good_jets>=3)
					for (int i=2;i<good_jets;i++)
						if (jets_pv[i].Pt() > 15 ) AdJetHT+=jets_pv[i].Pt();
				if (good_jets==2) hAdJetHT->Fill(10.,genweight);
				if (good_jets>=3) hAdJetHT->Fill(AdJetHT,genweight);
				hsoftleadTrackPt->Fill(Jet.EWKsoft_pt[0],genweight);	
				hsoftleadTrackEta->Fill(Jet.EWKsoft_eta[0],genweight);	
			
				hJets12_pt->Fill((jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJets12_pt_log->Fill(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJet1q_pt_log->Fill(TMath::Log(jets_pv[0].Pt()),genweight);
				hJet2q_pt_log->Fill(TMath::Log(jets_pv[1].Pt()),genweight);
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




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////// 2D HISTOS ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<float> variables_in_2D_plot = {Zll_mass, deltaM/*, log(X_parton1), log(X_parton2)*/, deltaR1, deltaR2, Rpt, log(Zll_zstar), atanh((TMVA.BDToutput+1)/2), genJetMass_leading, genJetMass_matched, Mqq};

        
//         for (int i = 0; i < variables_in_2D_plot.size(); i++) variablesName_in_2D_plot.push_back(getName(variables_in_2D_plot[i]));
        
        int idx = 0;
        for (int i = 0; i < variables_in_2D_plot.size(); i++) {
            for (int j = i+1; j < variables_in_2D_plot.size(); j++) {
                histo2D_vector[idx]->Fill( variables_in_2D_plot[i], variables_in_2D_plot[j], genweight);;
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






        }
// ------------------------------------------------------------------------------------END EventLoop----------------------------

//         hSelectionCuts->Scale(xsec[file_tag]/events_generated);



	treeMVA->Write();
	fileMVA.Close();


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
		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".root","recreate");
// 		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+"_noQGLcorrection.root","recreate");



            for (int i=0;i<numArray;++i){

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

	 ofstream out(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".txt");
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

return 0;
    
}
