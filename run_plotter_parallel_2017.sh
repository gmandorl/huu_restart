#!/bin/bash

# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  

QCDcorrection=$1
JEScorrection=$2
JERcorrection=$3
PUcorrection=$4


Z=$5

fileSkim2016=fileSkim2016_new
fileSkim2017=fileSkim2017
fileSkim2018=fileSkim2018

if [ $Z = 1 ]; then
    fileSkim2016=fileSkim2016_Z
    fileSkim2017=fileSkim2017_Z
    fileSkim2018=fileSkim2018_Z
fi



echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection "     " JERcorrection $JERcorrection "     " PUcorrection $PUcorrection "     " $executable "     " $v "     " $label "    " $era=$8;







 
#  mu  0 0 nom 0 nom nom nom nano 2017 histoFileDir
#  mu  0 0 nom 1 nom nom nom nano 2017 histoFileDir
 



 ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX_nano2016.root DYJetsToLL_M-105To160-amcatnloFXFX mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir &
# 
# 
# # ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_M-105To160-madgraphMLM_nano2017.root DYJetsToLL_M-105To160-madgraphMLM mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir &
# # ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_M-105To160_VBFFilter-madgraphMLM_nano2016.root DYJetsToLL_M-105To160_VBFFilter-madgraphMLM mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir 
#  ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2018/DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX_nano2018.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir &
 

# echo "DYJetstoLL_amc_2J"
# ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/DYJetsToLL_0J_nano2017.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir &
# ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/DYJetsToLL_1J_nano2017.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/DYJetsToLL_2J_nano2017.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir &
# echo "End DYJetstoLL_amc_2J"

./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/DYJetsToTauTau_ForcedMuDecay_nano2017.root DYJetsToTauTau_ForcedMuDecay mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;



./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/EWK_LLJJ_MLL-50_MJJ-120-madgraph-herwigpp_nano2016.root EWK_LLJJ_herwig mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/EWK_LLJJ_pythia8_nano2016.root EWK_LLJJ_pythia8 mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/EWK_LLJJ_INT_nano2018.root EWK_LLJJ_INT mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2017 histoFileDir



 ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/TTTo2L2Nu_nano2017.root TTTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir 


 
 echo "Data"
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/SingleMuon_nano2017.root SingleMuon mu  1 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir
echo "End Data"


echo "Single Top"
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ST_tW_top_nano2017.root ST_tW_top mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ST_tW_antitop_nano2017.root ST_tW_antitop mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_s_nano2016.root ST_s-channel mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ST_t_top_nano2017.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ST_t_antitop_nano2017.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
echo "End Single Top"
  
echo "WToLnu"
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WToLNu_0J_nano2017.root WToLNu_0J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WToLNu_1J_nano2016.root WToLNu_1J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WToLNu_2J_nano2017.root WToLNu_2J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
echo "End WToLnu"   


echo "TT and Diboson"

./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/TT_nano2017.root TT mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/TTToSemiLeptonic_nano2017.root TTToSemiLeptonic mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/TTToHadronic_nano2017.root TTToHadronic mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;


./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWTo2L2Nu_DoubleScattering_nano2017.root WWTo2L2Nu_DoubleScattering mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WLLJJ_WToLNu_EWK_nano2017.root WLLJJ_WToLNu_EWK mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWJJToLNuLNu_EWK_nano2017.root WWJJToLNuLNu_EWK mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWJJToLNuLNu_EWK_noTop_nano2017.root WWJJToLNuLNu_EWK_noTop mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;




./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWTo1L1Nu2Q_nano2017.root WWTo1L1Nu2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWTo4Q_nano2017.root WWTo4Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WWTo2L2Nu_nano2017.root WWTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;

./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZTo3LNu_nano2017.root WZTo3LNu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZTo2L2Q_nano2017.root WZTo2L2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZTo1L1Nu2Q_nano2017.root WZTo1L1Nu2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZTo1L3Nu_nano2017.root WZTo1L3Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
# ./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZJToLLLNu_nano2017.root WZJToLLLNu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;




./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo2L2Q_nano2017.root ZZTo2L2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo4L_nano2017.root ZZTo4L mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo2L2Nu_nano2017.root ZZTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;



./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WW_nano2016.root WW mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WZ_nano2017.root WZ mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZ_nano2017.root ZZ mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
echo "End TT and Diboson"

 
echo "Signal"
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/VBF_HToMuMu_nano2017.root VBF_HToMuMu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir 
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/GluGluHToMuMu_nano2017.root GluGlu_HToMuMu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2017 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ttHToMuMu_nano2017.root ttHToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WplusH_HToMuMu_nano2017.root WplusH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/WminusH_HToMuMu_nano2017.root WminusH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZH_HToMuMu_nano2017.root ZH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/GluGlu_HToMuMu_amc_nano2017.root GluGlu_HToMuMu_amc mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot17 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/VBF_HToMuMu_amc_nano2017.root VBF_HToMuMu_amc mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir 

echo "End Signal"


 
 
 
 
 
