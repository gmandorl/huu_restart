#!/bin/bash

# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  

QCDcorrection=$1
JEScorrection=$2
JERcorrection=$3
PUcorrection=$4

executable=./$5
v=$6
label=$7

echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection "     " JERcorrection $JERcorrection "     " PUcorrection $PUcorrection "     " $executable "     " $v "     " $label ;




#############################################################################################################################################################################
##################################################################     FILESKIM    ##########################################################################################
#############################################################################################################################################################################




 
 
 
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-105To160-amcatnloFXFX_v25_reskim.root DYJetsToLL_M-105To160-amcatnloFXFX mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection v25 reskim histoFileDir &


# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-105To160-madgraphMLM_v25_reskim.root DYJetsToLL_M-105To160-madgraphMLM mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;



# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-50_VBFFilter-amcatnloFXFX_v25_reskim.root DYJetsToLL_M-50_VBFFilter-amcatnloFXFX mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-50_VBFFilter-madgraphMLM_v25_reskim.root DYJetsToLL_M-50_VBFFilter-madgraphMLM mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_v25_reskim.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-105To160_VBFFilter-madgraphMLM_v25_reskim.root DYJetsToLL_M-105To160_VBFFilter-madgraphMLM mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;




# echo "DYJetstoLL_amc_2J"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_0J_v25_reskim.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_1J_v25_reskim.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_2J_v25_reskim.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# echo "End DYJetstoLL_amc_2J"
 


#  echo "Data"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim2/SingleMuon_reminiaod_v25.root SingleMuon mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;   #  <------------- there is fileSkim2
#  echo "End Data"
# 
# 
#   echo "Single Top"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_tW_top_v25_reskim.root ST_tW_top mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_tW_antitop_v25_reskim.root ST_tW_antitop mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_s_v25_reskim.root ST_s-channel mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_t_top_v25_reskim.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_t_antitop_v25_reskim.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
#   echo "End Single Top"
#   
#   
# # # # #  ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/Hmumu_nTuples/TTToSemilepton_v25_reskim.root TTToSemilepton mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# # # # #  ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/Hmumu_nTuples/TTTo2L2Nu_v25_reskim.root TTTo2L2Nu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# 
# 
#  
#  
#  
#   echo "TT and Diboson"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/TT_v25_reskim.root TT mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WW_v25_reskim.root WW mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WZ_v25_reskim.root WZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ZZ_v25_reskim.root ZZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WJetsToLnu_madgraph_v25_reskim.root WJetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
#  echo "End TT and Diboson"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WToLNu_0J_v25_reskim.root WToLNu_0J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
#  ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WToLNu_1J_v25_reskim.root WToLNu_1J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
#  ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WToLNu_2J_v25_reskim.root WToLNu_2J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# 
# 
#  
#  echo "Signal"
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/VBF_HToMuMu_v25_reskim.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir &
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/GluGlu_HToMuMu_v25_reskim.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection v25 reskim histoFileDir;
# # $executable /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/GluGlu_HToMuMu_v25_reskim.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection $v reskim histoFileDir;
#  echo "End Signal"
 
 
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################








# ----------------------------------------------------------------- nanoAOD 2017 -----------------------------------------------------------------------------------------------------



# echo "DYJetstoLL"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DYJetsToTauTau_ForcedMuDecay_nano2017.root DYJetsToTauTau_ForcedMuDecay mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DYJetsToLL_M-50_nano2017.root DYJetstoLL_amc_M-50 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# 


# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY1JetsToLL_M-50_LHEZpT_50-150_nano2017.root  DY1JetsToLL_M-50_LHEZpT_50-150 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY1JetsToLL_M-50_LHEZpT_150-250_nano2017.root  DY1JetsToLL_M-50_LHEZpT_150-250 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY1JetsToLL_M-50_LHEZpT_250-400_nano2017.root  DY1JetsToLL_M-50_LHEZpT_250-400 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY1JetsToLL_M-50_LHEZpT_400-inf_nano2017.root  DY1JetsToLL_M-50_LHEZpT_400-inf mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# 
# 
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY2JetsToLL_M-50_LHEZpT_50-150_nano2017.root DY2JetsToLL_M-50_LHEZpT_50-150 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY2JetsToLL_M-50_LHEZpT_150-250_nano2017.root DY2JetsToLL_M-50_LHEZpT_150-250 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY2JetsToLL_M-50_LHEZpT_250-400_nano2017.root DY2JetsToLL_M-50_LHEZpT_250-400 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/DY2JetsToLL_M-50_LHEZpT_400-inf_nano2017.root DY2JetsToLL_M-50_LHEZpT_400-inf mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;

# echo "End DYJetstoLL"



# if [ "$QCDcorrection" = "nom" ]; then
# if [ "$JEScorrection" = "nom" ]; then
# if [ "$JERcorrection" = "nom" ]; then
# if [ "$PUcorrection"  = "nom" ]; then
#     echo "DATA:     " QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection "     " JERcorrection $JERcorrection "     " PUcorrection $PUcorrection  ;
#     echo "Data"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_nano2017.root SingleMuon mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;  
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_Run2017B_nano2017.root SingleMuonB mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;  
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_Run2017C_nano2017.root SingleMuonC mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;   
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_Run2017D_nano2017.root SingleMuonD mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;   
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_Run2017E_nano2017.root SingleMuonE mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;   
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_Run2017F_nano2017.root SingleMuonF mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;  
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/SingleMuon_ultimo.root SingleMuonG mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;   
# echo "End Data"
# fi
# fi
# fi
# fi



# echo "Single Top"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ST_tW_top_nano2017.root ST_tW_top mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ST_tW_antitop_nano2017.root ST_tW_antitop mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ST_s_nano2017.root ST_s-channel mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ST_t_top_nano2017.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ST_t_antitop_nano2017.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# echo "End Single Top"
#   
# 
#   
# 
# echo "TT and Diboson"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/TTToHadronic_nano2017.root TTToHadronic mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/TTToSemiLeptonic_nano2017.root TTToSemilepton mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/TTTo2L2Nu_nano2017.root TTTo2L2Nu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# 
# 
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WW_nano2017.root WW mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WZ_nano2017.root WZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/ZZ_nano2017.root ZZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# echo "End TT and Diboson"
# 
# 
# echo "W+Jet"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLnu_madgraph_nano2017.root WJetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-100To200_nano2017.root WJetsToLNu_HT-100To200 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-200To400_nano2017.root WJetsToLNu_HT-200To400 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-400To600_nano2017.root WJetsToLNu_HT-400To600 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-600To800_nano2017.root WJetsToLNu_HT-600To800 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-800To1200_nano2017.root WJetsToLNu_HT-800To1200 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-1200To2500_nano2017.root WJetsToLNu_HT-1200To2500 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/WJetsToLNu_HT-2500ToInf_nano2017.root WJetsToLNu_HT-2500ToInf mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# echo "End W+Jet"
# 
#  
# echo "Signal"
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/VBF_HToMuMu_nano2017.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2017/GluGlu_HToMuMu_nano2017.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# echo "End Signal"
 
 
 
 
 
 
 
 
 
 
 
 
 # ----------------------------------------------------------------- nanoAOD 2016 -----------------------------------------------------------------------------------------------------

 
#  mu  0 0 nom 0 nom nom nom nano 2016 histoFileDir
 

 


# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-100to200_nano2016.root DYJetsToLL_M-100to200 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-200to400_nano2016.root DYJetsToLL_M-200to400 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-400to500_nano2016.root DYJetsToLL_M-400to500 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-500to700_nano2016.root DYJetsToLL_M-500to700 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-700to800_nano2016.root DYJetsToLL_M-700to800 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # 
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_Pt-50To100_nano2016.root DYJetsToLL_Pt-50To100 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_Pt-100To250_nano2016.root DYJetsToLL_Pt-100To250 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_Pt-250To400_nano2016.root DYJetsToLL_Pt-250To400 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_Pt-400To650_nano2016.root DYJetsToLL_Pt-400To650 mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_Pt-650ToInf_nano2016.root DYJetsToLL_Pt-650ToInf mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # # # # 


# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX2_nano2016.root DYJetsToLL_M-105To160-amcatnloFXFX1 mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX1_nano2016.root DYJetsToLL_M-105To160-amcatnloFXFX2 mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir &
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX_nano2016.root DYJetsToLL_M-105To160-amcatnloFXFX mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir &


./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-madgraphMLM_nano2016.root DYJetsToLL_M-105To160-madgraphMLM mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir  &
 ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160_VBFFilter-madgraphMLM_nano2016.root DYJetsToLL_M-105To160_VBFFilter-madgraphMLM mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir 
 
 
 
echo "DYJetstoLL_amc_2J"
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetstoLL_amc_0J_v25_nano2016.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetstoLL_amc_1J_v25_nano2016.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetstoLL_amc_2J_v25_nano2016.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End DYJetstoLL_amc_2J"



# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/EWK_LLJJ_MLL-50_MJJ-120-madgraph-herwigpp_nano2016.root EWK_LLJJ_herwig mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/EKW_LLJJ_pythia8_nano2016.root EKW_LLJJ_pythia8 mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/EWK_LLJJ_INT_nano2016.root EWK_LLJJ_INT mu  0 0 $QCDcorrection 0 $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir




./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/TTTo2L2Nu_nano2016.root TTTo2L2Nu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
 
 echo "Data"
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/SingleMuon_nano2016.root SingleMuon mu  1 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir
echo "End Data"


echo "Single Top"
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ST_tW_top_nano2016.root ST_tW_top mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ST_tW_antitop_nano2016.root ST_tW_antitop mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ST_s_nano2016.root ST_s-channel mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ST_t_top_nano2016.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ST_t_antitop_nano2016.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End Single Top"
  
echo "WToLnu"
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WToLNu_0J_nano2016.root WToLNu_0J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WToLNu_1J_nano2016.root WToLNu_1J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WToLNu_2J_nano2016.root WToLNu_2J mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

# # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/W4JetsToLNu_nano2016.root W4JetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/W3JetsToLNu_nano2016.root W3JetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/W2JetsToLNu_nano2016.root W2JetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/W1JetsToLNu_nano2016.root W1JetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End WToLnu"


echo "TT and Diboson"

./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/TT_nano2016.root TT mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/TTToSemilepton_nano2016.root TTToSemilepton mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;


./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WpWpJJ_EWK_nano2016.root WpWpJJ_EWK mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WWTo2L2Nu_DoubleScattering_nano2016.root WWTo2L2Nu_DoubleScattering mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WLLJJ_WToLNu_EWK_nano2016.root WLLJJ_WToLNu_EWK mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WWJJToLNuLNu_EWK_noTop_nano2016.root WWJJToLNuLNu_EWK_noTop mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;




./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WW_nano2016.root WW mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WZ_nano2016.root WZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/ZZ_nano2016.root ZZ mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/WJetsToLnu_madgraph_nano2016.root WJetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End TT and Diboson"

 
echo "Signal"
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_nano2016.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# # # ./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/VBF_HToMuMu_P02_nano2016.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
./plot /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/GluGlu_HToMuMu_nano2016.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End Signal"


 
 
 
 
 
