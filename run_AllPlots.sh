#!/bin/bash

executable=plot
v=v25
lab=reskim
# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA

era=$1
Z=$2
boolZ=0

if [ $Z = "Z" ]; then
        boolZ=1
fi


# g++ ./plotter_vbfzll.C -g -o $executable `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/include  -lboost_thread -llwtnn -lboost_system -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2-ikhhed/include/eigen3 


# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/include  -lboost_thread -llwtnn -lboost_system -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2-ikhhed/include/eigen3 



# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  


# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/include  -lboost_thread -llwtnn -lboost_system -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2-ikhhed/include/eigen3



QCDcorrectionARRAY=(nom up down nom nom nom nom nom nom)
JEScorrectionARRAY=(nom nom nom up down nom nom nom nom)
JERcorrectionARRAY=(nom nom nom nom nom up down nom nom)
PUcorrectionARRAY=( nom nom nom nom nom nom nom up down)


if [ $era = 2016 ]; then
                g++ ./plotter_vbfzll.C -g -o plot16 `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA
fi

if [ $era = 2017 ]; then
                g++ ./plotter_vbfzll.C -g -o plot17 `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA
fi

if [ $era = 2018 ]; then
                g++ ./plotter_vbfzll.C -g -o plot18 `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA
fi



for i in $(seq 0 1 8 ); do QCDcorrection=${QCDcorrectionARRAY[$i]}; JEScorrection=${JEScorrectionARRAY[$i]}; JERcorrection=${JERcorrectionARRAY[$i]}; PUcorrection=${PUcorrectionARRAY[$i]};

# echo run_plotter_parallel.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $executable $v $lab $era;
# source run_plotter_parallel.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $executable $v $lab &
if [ $era = 2016 ]; then
               echo run_plotter_parallel_2016.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $executable $v $lab 2016
               source run_plotter_parallel_2016.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $boolZ &
fi

if [ $era = 2017 ]; then
               echo run_plotter_parallel_2017.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $executable $v $lab 2017
               source run_plotter_parallel_2017.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $boolZ &
fi

if [ $era = 2018 ]; then
               echo run_plotter_parallel_2018.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $executable $v $lab 2018
               source run_plotter_parallel_2018.sh $QCDcorrection $JEScorrection $JERcorrection $PUcorrection $boolZ &
fi

done







