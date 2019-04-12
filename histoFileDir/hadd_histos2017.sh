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



hadd -f Top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root TTToHadronic_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root   TTToSemilepton_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root TTTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root ST_tW_top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root


cp EWK_LLJJ_herwig_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root   &




echo "this is 1  " $1
if [ $1 = Z ]; then
    hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  DYJetsToTauTau_ForcedMuDecay_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
else
    hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  DYJetsToTauTau_ForcedMuDecay_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
fi



# hadd -f  WJetsToLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root &


hadd -f  Other_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WWTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WZ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  ZZ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WWTo2L2Nu_DoubleScattering_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WLLJJ_WToLNu_EWK_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WWJJToLNuLNu_EWK_noTop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root &



done


# DA ELIMINARE

# hadd -f TT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root TTToHadronic_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root   TTToSemilepton_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root TTTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root

# hadd -f  EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root EKW_LLJJ_pythia8_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root EWK_LLJJ_INT_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  &
# cp EWK_LLJJ_herwig_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  EWK_LLJJ_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root   &



# # # hadd -f DYJetstoTauTau_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root 
# cp DYJetsToTauTau_ForcedMuDecay_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  DYJetstoTauTau_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root




# hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
# hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root

# cp WW_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WWInclusive_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root

# 
# echo "this is 1  " $1
# if [ $1 = Z ]; then
#     hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
# else
#     hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
# fi
# 
# 
# hadd -f  WJetsToLNu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WToLNu_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root &
# 
# 
# hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root   ST_tW_top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
# 
# 
# cp WWTo2L2Nu_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root WW_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
# 
# done





# cp -f SingleMuon_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  SingleMuon_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_v25_reskim.root


# hadd -f SingleMuon_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root SingleMuonB_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root SingleMuonC_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root SingleMuonD_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root SingleMuonE_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root  SingleMuonF_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_JER${JERcorrection}_PU${PUcorrection}_nano_2017.root
 




