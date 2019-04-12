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









echo "this is 1  " $1
if [ $1 = Z ]; then
    hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root 
else
    hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root 
fi


hadd -f  Other_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WW_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WZ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  ZZ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WWTo2L2Nu_DoubleScattering_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WLLJJ_WToLNu_EWK_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WWJJToLNuLNu_EWK_noTop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root &


cp EWK_LLJJ_herwig_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  &

hadd -f Top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTToSemiLeptonic_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root &



# DA ELIMINARE





# hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetsToLL_Pt-50To100_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_Pt-100To250_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetsToLL_Pt-250To400_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_Pt-400To650_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root



# hadd -f SingleMuon_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB1_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonB2_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonC_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonD_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonE_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonF_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonG_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonH2_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root SingleMuonH3_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root


# mv TT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTinclusive_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root


# echo "this is 1  " $1
# if [ $1 = Z ]; then
#     hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root 
#     cp DYJetsToTauTau_ForcedMuDecay_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoTauTau_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root
# else
#     # cp -f DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root
#     hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root
#     # cp -f DYJetsToLL_M-105To160-madgraphMLM_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root
#     hadd -f DYJetstoTauTau_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root 
# fi



# hadd -f  WJetsToLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root &
# hadd -f  WJetsToLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root &



# hadd -f  EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EWK_LLJJ_herwig_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EWK_LLJJ_INT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root

# hadd -f  EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EKW_LLJJ_pythia8_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EWK_LLJJ_INT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  &

# cp EWK_LLJJ_herwig_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  &

# hadd -f  WWJJToLNuLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WpWpJJ_EWK_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WWJJToLNuLNu_EWK_noTop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  &
# hadd -f  WWJJToLNuLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WpWpJJ_EWK_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root WWJJToLNuLNu_EWK_noTop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root  &




# cp -f DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root DYJetstoTauTau_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root   


# hadd -f TT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTToSemilepton_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root TTTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root




# hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root   ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root &

# hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root   ST_tW_top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2016.root


done






# hadd -f SingleMuon_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB1_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonB2_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonC_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonD_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonE_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonF_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonG_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonH2_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root SingleMuonH3_mu_QCDScalenom_JESnom_JERnom_PUnom_nano_2016.root



