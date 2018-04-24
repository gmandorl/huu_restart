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







void SetStyle (TH1 * histo, float x) {

    
    float axisTitleSize=0.05*x;
    histo->GetYaxis()->SetTitleSize(axisTitleSize);
    histo->GetXaxis()->SetTitleSize(axisTitleSize);

    
    float axisLabelSize=0.05*x;
    histo->GetYaxis()->SetLabelSize(axisLabelSize);
    histo->GetXaxis()->SetLabelSize(axisLabelSize);

}




float computeSignificance(TH1F * hSignal, TH1F * hBackground) {
    
    float totalSensitivity = 0;

    for (int n = 1; n <= hSignal->GetXaxis()->GetNbins(); ++n) {
            if(hBackground->GetBinContent(n) > 0.000001) totalSensitivity += hSignal->GetBinContent(n)*hSignal->GetBinContent(n)/hBackground->GetBinContent(n);
        }
    
    return sqrt(totalSensitivity);
}




void macroSignificanceDY_LO_NLO () {
    std::cout << "I can see here" << std::endl;

    float lumi = 35900.;


    TFile* fileDYNLO= new TFile ("DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScalenom_JESnom_v25_reskim.root");
    TFile* fileDY_LO= new TFile ("DYJetsToLL_M-105To160-madgraphMLM_mu_QCDScalenom_JESnom_v25_reskim.root");
    TFile* fileHmumu= new TFile ("VBF_HToMuMu_mu_QCDScalenom_JESnom_v25_reskim.root");
    

    TH1F * hHmumu = (TH1F *) fileHmumu->Get("hBDT_VBF_atanh");
    TH1F * hDY_LO = (TH1F *) fileDY_LO->Get("hBDT_VBF_atanh");
    TH1F * hDYNLO = (TH1F *) fileDYNLO->Get("hBDT_VBF_atanh");

    hHmumu->SetLineColor(2);
    hDY_LO->SetLineColor(4);
    hDYNLO->SetLineColor(8);

//     hHmumu->Scale(lumi);
//     hDY_LO->Scale(lumi);
//     hDYNLO->Scale(lumi);


    hHmumu->SetLineWidth(2);
    hDY_LO->SetLineWidth(2);
    hDYNLO->SetLineWidth(2);


    float DY_normalization = 6071.;
    hHmumu->Scale(8.4/hHmumu->Integral());
    hDY_LO->Scale(DY_normalization/hDY_LO->Integral());
    hDYNLO->Scale(DY_normalization/hDYNLO->Integral());

    hDYNLO->GetYaxis()->SetRangeUser(0.01, 1000);

    
    float significance_LO   = computeSignificance(hHmumu, hDY_LO);
    float significance_NLO  = computeSignificance(hHmumu, hDYNLO);

    float Xbottom = 0.5;
    float Ybottom = 0.75;
    float Xtop = 0.9;
    float Ytop = 0.9;
    TLegend *myLegend=new TLegend(Xbottom, Ybottom, Xtop, Ytop, "");
    myLegend->SetBorderSize(0);
    myLegend->SetFillColor(0);
    myLegend->SetFillStyle(0);
    myLegend->SetTextFont(42);
    myLegend->SetTextSize(0.03);
    myLegend->SetBorderSize(0);                  // without border
    myLegend->SetFillColorAlpha(1,0);           // trasparent

    myLegend->AddEntry(hHmumu,"Signal","l");
    myLegend->AddEntry(hDY_LO,("DY LO, significance: "+to_string(significance_LO)).c_str(),"l");
    myLegend->AddEntry(hDYNLO,("DY NLO, significance: "+to_string(significance_NLO)).c_str(),"l");

    
    TCanvas * canv = new TCanvas("canv", "", 800, 600);
    canv->cd();
    gStyle->SetOptStat(0000); 
    gPad->SetLogy();
    
    
    
    
    hDYNLO->Draw();
    hHmumu->Draw("same");
    hDY_LO->Draw("same");
    
    myLegend->Draw();
    
    canv->Print("../figure_LO_NLO/LO_NLO_comparison/BDT_comparison.png");


    
    
    
    std::vector<std::string> histonames;
    
    histonames.push_back("hMqq");
    histonames.push_back("hRpt");
     
    histonames.push_back("hZll_mass");
    histonames.push_back("hZll_pt");
    
    histonames.push_back("hlepton1_pt");
    histonames.push_back("hlepton2_pt");
    
    histonames.push_back("hJet1q_pt");
    histonames.push_back("hJet2q_pt");
    
    histonames.push_back("hPzAbsLog");
    histonames.push_back("hVirtual_Pt1_log");
    histonames.push_back("hVirtual_Pt2_log");
    histonames.push_back("hVirtual_Wmass1_log");
    histonames.push_back("hVirtual_Wmass2_log");
    
    histonames.push_back("hqgl");
    histonames.push_back("hqgl2");
    histonames.push_back("hqglAtanh");
    histonames.push_back("hqgl2Atanh");
    
    histonames.push_back("hqq_pt");
    histonames.push_back("hmet");
    
    
    
    histonames.push_back("hJet3_pt");
    histonames.push_back("hlheNj");
    histonames.push_back("hJets12_pt");
    
    histonames.push_back("hThetaStarJetAtanh");
    histonames.push_back("hThetaPlanesAtanh");
    histonames.push_back("hdeltaMRel");
    
    histonames.push_back("hMaxJetBTagCSV");
    histonames.push_back("hMaxJetBTagCMVA");
    
    histonames.push_back("hWWmass");
    histonames.push_back("hgen_mass");
    
    histonames.push_back("hEnergy_fraction_Parton1_log");
    histonames.push_back("hEnergy_fraction_Parton2_log");
    
    
    histonames.push_back("hmumujj_ptLog");
    histonames.push_back("hInvariant_MassLog");
    histonames.push_back("hPzAbs");
    
    
    


    for (int n = 0; n < histonames.size();n++) {
        
        TCanvas * canv_ = new TCanvas("canv_", "", 800, 600);
        canv_->cd();
        
        TPad * pad1 = new TPad("pad1", "",0, 0.3, 1, 1.0);
        pad1->SetFillStyle(4000);
        pad1->SetTickx(1);
        pad1->SetTicky(1);
        pad1->SetBottomMargin(0.0);

        pad1->SetTopMargin(0.12);
        pad1->Draw();
        pad1->cd();
    
        gStyle->SetOptStat(0000); 
    
    
        TH1F * hDY_LO = (TH1F *) fileDY_LO->Get(histonames[n].c_str());
        TH1F * hDYNLO = (TH1F *) fileDYNLO->Get(histonames[n].c_str());
        
        hDY_LO->Scale(DY_normalization/hDY_LO->Integral());
        hDYNLO->Scale(DY_normalization/hDYNLO->Integral());
    
    
        hDY_LO->SetLineColor(4);
        hDYNLO->SetLineColor(8);


        hDY_LO->SetLineWidth(2);
        hDYNLO->SetLineWidth(2);
        
        
        
        
        
        float Xbottom = 0.5;
        float Ybottom = 0.75;
        float Xtop = 0.9;
        float Ytop = 0.9;
        TLegend *myLegend=new TLegend(Xbottom, Ybottom, Xtop, Ytop, "");
        myLegend->SetBorderSize(0);
        myLegend->SetFillColor(0);
        myLegend->SetFillStyle(0);
        myLegend->SetTextFont(42);
        myLegend->SetTextSize(0.03);
        myLegend->SetBorderSize(0);                  // without border
        myLegend->SetFillColorAlpha(1,0);           // trasparent

        myLegend->AddEntry(hDY_LO,"DY LO","l");
        myLegend->AddEntry(hDYNLO,"DY NLO","l");
    
    
    
        hDYNLO->GetYaxis()->SetRangeUser(0.1, hDYNLO->GetMaximum()*1.2);
        
        hDYNLO->Draw();
        hDY_LO->Draw("same");
        myLegend->Draw();
    
        
        
        canv_->cd();
        TPad * pad2 = new TPad("pad2", "",0, 0.01, 1, 0.3);
        pad2->SetFillStyle(4000);
        pad2->SetTickx(1);
        pad2->SetTicky(1);
        pad2->SetBottomMargin(0.25);

        pad2->SetTopMargin(0.);
        pad2->Draw();
        pad2->cd();
    
        gStyle->SetOptStat(0000);
        gPad->SetGrid();  
    
        TH1F * hRatio = (TH1F*) hDYNLO->Clone(histonames[n].c_str());
        hRatio->SetLineColor(1);
        
        hRatio->SetTitle("");
        hRatio->GetYaxis()->SetTitle("NLO/LO");
        hRatio->GetYaxis()->SetTitleSize(0.05);
        hRatio->GetYaxis()->SetTitleOffset(1.);
        
        hRatio->GetYaxis()->SetRangeUser(0.6, 1.4);
        hRatio->Divide(hDY_LO);
        hRatio->GetYaxis()->SetTitleOffset(0.4);
        SetStyle(hRatio, 2.5);
        hRatio->Draw(); 
    
    
    
        canv_->Print(("../figure_LO_NLO/LO_NLO_comparison/"+histonames[n]+".png").c_str());
        
    }
    
    
}



        
        
        
        
