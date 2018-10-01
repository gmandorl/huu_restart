#!/bin/bash


# QCDcorrection=nom
# QCDcorrection=up
# QCDcorrection=down


# JEScorrection=nom
# JEScorrection=up
# JEScorrection=down



QCDcorrectionARRAY=(nom up down nom nom nom nom nom nom)
JEScorrectionARRAY=(nom nom nom up down nom nom nom nom)
JERcorrectionARRAY=(nom nom nom nom nom up down nom nom)
PUcorrectionARRAY=(nom nom nom nom nom nom nom up down)













for i in $(seq 0 1 8); do QCDcorrection=${QCDcorrectionARRAY[$i]}; JEScorrection=${JEScorrectionARRAY[$i]}; JERcorrection=${JERcorrectionARRAY[$i]}; PUcorrection=${PUcorrectionARRAY[$i]};


echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection "     " JERcorrection $JERcorrection "     " PUcorrection $PUcorrection


# hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetsToLL_Pt-50To100_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_Pt-100To250_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetsToLL_Pt-250To400_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_Pt-400To650_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root



# hadd -f SingleMuon_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB1_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonB2_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonC_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonD_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonE_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonF_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonG_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonH2_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonH3_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root


# mv TT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTinclusive_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root



cp -f DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root

hadd -f  WJetsToLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root




hadd -f TT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTToSemilepton_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root



hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root   ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root

# hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root   ST_tW_top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root


done






# hadd -f SingleMuon_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB1_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB2_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonC_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonD_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonE_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonF_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonG_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonH2_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonH3_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root



