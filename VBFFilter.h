#ifndef VBFFilter_h
#define VBFFilter_h

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


const float deltaRNoLep = 0.3;








// float VBFFilter(std::vector<TLorentzVector> genLeptons, std::vector<TLorentzVector> genJet) {
std::vector<TLorentzVector> VBFFilter(std::vector<TLorentzVector> genLeptons, std::vector<TLorentzVector> genJet) {
    
    

//     std::cout << "lep and jet: " << genLeptons.size() << " \t " << genJet.size() << std::endl;
      
      // Getting p4 of jet with no lepton
      std::vector<TLorentzVector> genJetsWithoutLeptonsP4;
      unsigned int jetIdx = 0;
      
      
      

//       while(genJetsWithoutLeptonsP4.size()<2 && jetIdx < genJet.size()) {
      while(jetIdx < genJet.size()) {
          bool jetWhitoutLep = true;
          const TLorentzVector & p4J= genJet[jetIdx];
          for(unsigned int i = 0; i < genLeptons.size() && jetWhitoutLep; ++i) {
              if(genLeptons[i].DeltaR(p4J) < deltaRNoLep)
              
                  jetWhitoutLep = false;
          }
          
          if (jetWhitoutLep)  genJetsWithoutLeptonsP4.push_back(p4J);
          ++jetIdx;
      }
      
      
      
      // Checking the invariant mass of the leading jets
//       if (genJetsWithoutLeptonsP4.size() < 2) return 0.;
      
//       float invMassLeadingJet = (genJetsWithoutLeptonsP4[0] + genJetsWithoutLeptonsP4[1]).M();
//       std::cout << "invMassLeadingJet: " << invMassLeadingJet << std::endl;
//       return invMassLeadingJet;
      return genJetsWithoutLeptonsP4;


    
      
      
      
      
      
      
      
      
      
      ///////////////////////////////////////////////////////////////// ALTRO FILTRO //////////////////////////////////////////////////////////////////////////////////
      
//     float maxDiJetMass = 0;
//     
//     for(unsigned a=0; a<genJet.size(); a++){
//         for(unsigned b=a+1; b<genJet.size(); b++){    
//         
//         TLorentzVector pA = genJet[a];
//         TLorentzVector pB = genJet[b];
//         
//         // Getting the dijet vector
//         TLorentzVector diJet = pA + pB;
//         
// 
//         
//         // Testing dijet mass
//         double invMass = diJet.M();
//         if(maxDiJetMass<=invMass) maxDiJetMass = invMass;
//         
//         }
//     }
//     
//     return maxDiJetMass;
    
    //////////////////////////////////////////////////////////////////////////77777///////////////////////////////////////////////////////////////////////////////////



}








#endif
