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
#     fileSkim2016=provaJEC
    fileSkim2017=fileSkim2017_Z
    fileSkim2018=fileSkim2018_Z
fi


echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection "     " JERcorrection $JERcorrection "     " PUcorrection $PUcorrection "     " $executable "     " $v "     " $label "    " $era=$8;







 
#  mu  0 0 nom 0 nom nom nom nano 2016 histoFileDir
#  mu  0 0 nom 1 nom nom nom nano 2016 histoFileDir



# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX_nano2016.root DYJetsToLL_M-105To160-amcatnloFXFX mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir &


# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-madgraphMLM_nano2016.root DYJetsToLL_M-105To160-madgraphMLM mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir  &
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160-madgraphMLM_nano2016.root DYJetsToLL_M-105To160-madgraphMLM mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/fileSkim2016/DYJetsToLL_M-105To160_VBFFilter-madgraphMLM_nano2016.root DYJetsToLL_M-105To160_VBFFilter-madgraphMLM mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir 
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_M-105To160-amcatnloFXFX_VBFFilter_nano2016.root DYJetsToLL_M-105To160_VBFFilter-amcatnloFXFX mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir 
 
 
# echo "DYJetstoLL_amc_2J"
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_2J_nano2016.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_1J_nano2016.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/DYJetsToLL_0J_nano2016.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
# echo "End DYJetstoLL_amc_2J"



./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/EWK_LLJJ_MLL-50_MJJ-120-madgraph-herwigpp_nano2016.root EWK_LLJJ_herwig mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/EWK_LLJJ_pythia8_nano2016.root EWK_LLJJ_pythia8 mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/EWK_LLJJ_INT_nano2016.root EWK_LLJJ_INT mu  0 0 $QCDcorrection $Z $JEScorrection $JERcorrection  $PUcorrection nano 2016 histoFileDir




./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/TTTo2L2Nu_nano2016.root TTTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &
 
 echo "Data"
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/SingleMuon_nano2016.root SingleMuon mu  1 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir
echo "End Data"


echo "Single Top"
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_tW_top_nano2016.root ST_tW_top mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_tW_antitop_nano2016.root ST_tW_antitop mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_s_nano2016.root ST_s-channel mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_t_top_nano2016.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ST_t_antitop_nano2016.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End Single Top"
  
echo "WToLnu"
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WToLNu_0J_nano2016.root WToLNu_0J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WToLNu_1J_nano2016.root WToLNu_1J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WToLNu_2J_nano2016.root WToLNu_2J mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

echo "End WToLnu"


echo "TT and Diboson"

./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/TT_nano2016.root TT mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/TTToSemiLeptonic_nano2016.root TTToSemiLeptonic mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/TTJets_DiLept_nano2016.root TTJets_DiLept mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir &


./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WWTo2L2Nu_DoubleScattering_nano2016.root WWTo2L2Nu_DoubleScattering mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WLLJJ_WToLNu_EWK_nano2016.root WLLJJ_WToLNu_EWK mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WWJJToLNuLNu_EWK_noTop_nano2016.root WWJJToLNuLNu_EWK_noTop mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WWJJToLNuLNu_EWK_nano2016.root WWJJToLNuLNu_EWK mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;


./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WWTo1L1Nu2Q_nano2016.root WWTo1L1Nu2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WWTo2L2Nu_nano2016.root WWTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZTo3LNu_nano2016.root WZTo3LNu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZTo2L2Q_nano2016.root WZTo2L2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZTo1L1Nu2Q_nano2016.root WZTo1L1Nu2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZTo1L3Nu_nano2016.root WZTo1L3Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZJToLLLNu_nano2016.root WZJToLLLNu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo2L2Q_nano2017.root ZZTo2L2Q mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
# # # ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo4L_nano2017.root ZZTo4L mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;    non esiste
# # # ./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2017/ZZTo2L2Nu_nano2017.root ZZTo2L2Nu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;     non esiste


./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WW_nano2016.root WW mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WZ_nano2016.root WZ mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ZZ_nano2016.root ZZ mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
echo "End TT and Diboson"

 
echo "Signal"
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/VBF_HToMuMu_nano2016.root VBF_HToMuMu mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir 
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/GluGlu_HToMuMu_nano2016.root GluGlu_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ttHToMuMu_nano2016.root ttHToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WplusH_HToMuMu_nano2016.root WplusH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/WminusH_HToMuMu_nano2016.root WminusH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/ZH_HToMuMu_nano2016.root ZH_HToMuMu mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;

./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/GluGlu_HToMuMu_amc_nano2016.root GluGlu_HToMuMu_amc mu  0 0 $QCDcorrection $Z     $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir;
./plot16 /gpfs/ddn/cms/user/mandorli/Hmumu/CMSSW_9_4_6/src/Skim0/$fileSkim2016/VBF_HToMuMu_amc_nano2016.root VBF_HToMuMu_amc mu  0 0 $QCDcorrection $Z $JEScorrection  $JERcorrection  $PUcorrection nano 2016 histoFileDir 

echo "End Signal"


 
 
 
 
 
