import ROOT 
import numpy as np
import root_numpy

#remake the trainign fiel to get the normalizations, 
import root_numpy
from keras.layers import Dense, Dropout, Flatten, Convolution2D, merge, Convolution1D,concatenate
from keras.models import Model, load_model
from keras.layers import Input
#model=load_model("WEIGHTS_at_epoch99.h5")
#model=load_model("WEIGHTS_at_epoch99.h5")

#f=ROOT.TFile("tmva_tree_DYJetsToLL_M-105To160-madgraphMLM_allEvents.root")
#f=ROOT.TFile("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/mvaTree/main_tmva_tree_DYJetsToLL_M-105To160-madgraphMLM_v25mu_QCDScaleup_JESnom.root")
#f=ROOT.TFile("/afs/cern.ch/work/g/gimandor/public/sampleAndreaPerNN/training_29gennaio/tmva_tree_DYJetsToLL_M-105To160-madgraphMLM_allEvents.root")
#f=ROOT.TFile("rootToComputeAverages/main_tmva_tree_DYJetsToLL_M-105To160-madgraphMLM_v25mu_QCDScalenom_JESnom.root")
f=ROOT.TFile("rootToComputeAverages/main_tmva_tree_DYJetsToLL_M-105To160-madgraphMLM_v25mu_QCDScalenom_JESnom_JERnom_PUnom.root")
tree=f.Get("TMVAtree")
t=root_numpy.tree2array(tree)
t2=root_numpy.rec2array(t)

#fs=ROOT.TFile("tmva_tree_VBF_HToMuMu_allEvents.root")
#fs=ROOT.TFile("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/mvaTree/main_tmva_tree_VBF_HToMuMu_v25mu_QCDScalenom_JESnom.root")
#fs=ROOT.TFile("/afs/cern.ch/work/g/gimandor/public/sampleAndreaPerNN/training_29gennaio/tmva_tree_VBF_HToMuMu_allEvents.root")
#fs=ROOT.TFile("rootToComputeAverages/main_tmva_tree_VBF_HToMuMu_v25mu_QCDScalenom_JESnom.root")
fs=ROOT.TFile("rootToComputeAverages/main_tmva_tree_VBF_HToMuMu_v25mu_QCDScalenom_JESnom_JERnom_PUnom.root")
trees=fs.Get("TMVAtree")
ts=root_numpy.tree2array(trees)



print ts
print ts[0]

ts2=root_numpy.rec2array(ts)

print "Lenghts:"
print len(t2)
print len(ts2)
data=np.concatenate((t2,ts2))
print data[0]

#print "Prima ", data[0,80]
#data[:,80] = data[:,0]
#print "Dopo  ", data[0,80]

#print "Prima ", data[0,81]
#data[:,81] = data[:,0]
#print "Dopo   ", data[0,81]


col_mean = np.nanmean(data,axis=0)  
col_std = np.nanstd(data,axis=0)
col_mean[-1]=0
col_std[-1]=1
print col_mean, col_std




#[["ll_mass", "ll_pt", "ll_zstar","Mqq","qgl_1q","qgl_2q","Rpt","softActivityEWK_njets5","W_mass_virtual1","theta1","theta2","deltaM","btagCMVA", "genweight"]]

##########################   THE FIRST TIME THOSE WERE THE VARIABLES NUMBER   #########################################
###col_mean[0] = 124.313512     # ll_mass
###col_mean[2] = 107.231412    # ll_pt
###col_mean[5] = -1.51578777     # ll_zstar
###col_mean[18] = 928.204296   # Mqq
###col_mean[28] = 0.605277257   # qgl_1q
###col_mean[29] = 0.563018216   # qgl_2q
###col_mean[31] = 0.116890493    # RptHard
###col_mean[35] = 0.813349207   # softActivityEWK_njets5
###col_mean[43] = -95.3516522   # W_mass_virtual1
###col_mean[54] = -0.0837856276   # theta1
###col_mean[55] = -0.00891449648    # theta2
###col_mean[65] = 3.46269105   # deltaM
####col_mean[] = 0.767483536    # btagCMVA


###col_std[0] = 3.25459306
###col_std[2] = 71.2830539
###col_std[5] = 1.16217068
###col_std[18] = 698.935084   # Mqq
###col_std[28] = 0.334484925   # qgl_1q
###col_std[29] = 0.311992516   # qgl_2q
###col_std[31] = 0.0860204832    # RptHard
###col_std[35] = 1.18478729   # softActivityEWK_njets5
###col_std[43] = 67.5337482   # W_mass_virtual1
###col_std[54] = 0.742478121   # theta1
###col_std[55] = 0.747562268    # theta2
###col_std[65] = 1.66275975   # deltaM
####col_std[] = 0.382634284   # btagCMVA



#54 thetastarW1
#59 theta1


col_mean[0] = 124.313512     # ll_mass
col_mean[2] = 107.231412    # ll_pt
col_mean[5] = -1.51578777     # ll_zstar
col_mean[18] = 928.204296   # Mqq
col_mean[28] = 0.605277257   # qgl_1q
col_mean[29] = 0.563018216   # qgl_2q
col_mean[31] = 0.116890493    # RptHard
col_mean[35] = 0.813349207   # softActivityEWK_njets5
col_mean[43] = -95.3516522   # W_mass_virtual1
col_mean[59] = -0.0837856276   # theta1
col_mean[60] = -0.00891449648    # theta2
col_mean[66] = 3.46269105   # deltaM
#col_mean[] = 0.767483536    # btagCMVA


col_std[0] = 3.25459306
col_std[2] = 71.2830539
col_std[5] = 1.16217068
col_std[18] = 698.935084   # Mqq
col_std[28] = 0.334484925   # qgl_1q
col_std[29] = 0.311992516   # qgl_2q
col_std[31] = 0.0860204832    # RptHard
col_std[35] = 1.18478729   # softActivityEWK_njets5
col_std[43] = 67.5337482   # W_mass_virtual1
col_std[59] = 0.742478121   # theta1
col_std[60] = 0.747562268    # theta2
col_std[66] = 1.66275975   # deltaM
#col_std[] = 0.382634284   # btagCMVA




#for n in range(col_mean.size-1) :
    #print n, " \t " , col_mean[n] , " \t " , col_std[n]




ST_list = ["ST_t-channel_antitop_4f_inclusiveDecays","ST_t-channel_top_4f_inclusiveDecays","ST_tW_antitop","ST_tW_top"]
#ST_list = ["ST_t-channel_top_4f_inclusiveDecays","ST_tW_antitop","ST_tW_top"]
#ST_list = ["ST_t-channel_top_4f_inclusiveDecays","ST_tW_antitop"]
DY_list = ["DYJetsToLL_M-105To160-amcatnloFXFX","DYJetsToLL_M-105To160-madgraphMLM"]
fileToRewrite = ["SingleMuon","WToLNu_2J","W3JetsToLNu", "W4JetsToLNu", "WW","ZZ","WZ","TT", "TTToSemilepton", "TTTo2L2Nu","VBF_HToMuMu", "GluGlu_HToMuMu"] + ST_list + DY_list
#fileToRewrite = ["SingleMuon","WToLNu_2J","WW","TT","VBF_HToMuMu", "GluGlu_HToMuMu"] + ST_list + DY_list  # only if Mqq > 600
#fileToRewrite = ["SingleMuon","WToLNu_2J"] 
QCD_varation = ["nom", "up", "down", "nom", "nom", "nom", "nom", "nom", "nom"]
JES_varation = ["nom", "nom", "nom", "up", "down", "nom", "nom", "nom", "nom"]
JER_varation = ["nom", "nom", "nom", "nom", "nom", "up", "down", "nom", "nom"]
PU_varation  = ["nom", "nom", "nom", "nom", "nom", "nom", "nom", "up", "down"]


#QCD_varation = ["nom"]
#JES_varation = ["nom"]



#ind=np.r_[2:3,5:6,18:19,28:29,29:30,31:32,35:36,43:44,59:60,60:61]
#ind2=np.r_[0:1,65:66]


ind=np.r_[2:3,5:6,18:19,28:29,29:30,31:32,35:36,43:44,59:60,60:61]
ind2=np.r_[0:1,66:67]

dropoutRate=0.1

dim=len(ind)
print dim,ind,ind2
dim2=len(ind2)



Inputs=Input(shape=(dim,))
newInput=Input(shape=(dim2,))



x=  Dense(50, activation='relu',init='lecun_uniform',trainable=True)(Inputs)
x = Dropout(dropoutRate)(x)
x=  Dense(50, activation='relu',init='lecun_uniform',trainable=True)(x)
x = Dropout(dropoutRate)(x)
x=  Dense(30, activation='relu',init='lecun_uniform',trainable=True)(x)
x = Dropout(dropoutRate)(x)
x=  Dense(4, activation='relu',init='lecun_uniform',trainable=True)(x)
x = Dense(1, activation='sigmoid',init='lecun_uniform', trainable=False)(x)

x= concatenate([newInput,x])
x=  Dense(50, activation='relu',init='lecun_uniform')(x)
x = Dropout(dropoutRate)(x)
x=  Dense(30, activation='relu',init='lecun_uniform')(x)
x = Dropout(dropoutRate)(x)
predictions = Dense(1, activation='sigmoid',init='lecun_uniform')(x)
model = Model(input=[Inputs,newInput], output=predictions)

model.load_weights("WEIGHTS_at_epoch99.h5")


for fileName in fileToRewrite :
    for n in range(9) :
    #for qcd in QCD_varation :
        #for jes in JES_varation :
                qcd = QCD_varation[n]
                jes = JES_varation[n]
                jer = JER_varation[n]
                pu  = PU_varation[n]
                if fileName == "SingleMuon" and not (qcd == "nom" and jes == "nom" and jer == "nom" and pu == "nom") : continue
                print "file   \t ", fileName , qcd , jes, jer , pu

                #open the actual tree
                f_tocompute=ROOT.TFile("/scratch/mandorli/Hmumu/restartFromFilippo/CMSSW_8_0_28/src/code/mvaTree/main_tmva_tree_"+fileName+"_2016"+"_mu_QCDScale"+qcd+"_JES"+jes+"_JER"+jer+"_PU"+pu+".root")
                t_tocompute=f_tocompute.Get("TMVAtree")
                t_tocompute=root_numpy.tree2array(t_tocompute)
                #t_tocompute=t_tocompute[["ll_mass", "ll_pt", "ll_zstar","Mqq","qgl_1q","qgl_2q","Rpt","softActivityEWK_njets5","W_mass_virtual1","theta1","theta2","deltaM"]]
                
                
                
                array_tocompute=root_numpy.rec2array(t_tocompute)
                array_tocompute[:,43]=-array_tocompute[:,43]
                array_tocompute[:,65]=(2.**0.5)*array_tocompute[:,65]
                array_tocompute[:]=array_tocompute[:]-col_mean
                array_tocompute[:]=array_tocompute[:]/col_std

                #variables used in training
                #ind=np.r_[2:3,5:6,18:19,24:25,25:26,27:28,30:31,38:39,54:55,55:56]
                #ind2=np.r_[0:1,63:64]
                
                
                #ind=np.r_[1:2,2:3,3:4,4:5,5:6,6:7,7:8,8:9,9:10,10:11]
                #ind2=np.r_[0:1,11:12]
                
                
                predictions=model.predict([array_tocompute[:,ind],array_tocompute[:,ind2]])

                #print fileName , qcd , jes
                print "shape array is   ", array_tocompute.shape
                print "shape is   ", predictions.shape

                ###continue here
                import numpy.lib.recfunctions as rfn
                fnew=ROOT.TFile("ROOT/main_tmva_tree_"+fileName+"_v25mu_QCDScale"+qcd+"_JES"+jes+"_JER"+jer+"_PU"+pu+".root","RECREATE")
                fnew.cd()
                new = rfn.append_fields(t_tocompute, ['NNout'], [predictions[:,0]])
                final_tree=root_numpy.array2tree(new)
                final_tree.Write()
                fnew.Close()


                f_tocompute.Close()
                fs.Close()
                f.Close()






