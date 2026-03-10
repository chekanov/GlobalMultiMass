#!/bin/bash
# S.Chekanov (ANL)

echo "Set ROOT enviroment for Dijet+Lepton program"

HH=`hostname -A` 

echo "HOST=$HH"

export STORAGE="/lcrc/group/ATLAS/atlasfs/local/chekanov/AnomalyDetect23/data/"
export STORAGE_MC="/lcrc/group/ATLAS/atlasfs/local/chekanov/AnomalyDetect23/mc/"

if [[ $HH =~ .*atlaslogin.* ]]; then
  echo "HEP ANL ROOT setup"
  # source /users/admin/share/sl7/setup.sh 
  source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos9-gcc11-opt/setup.sh  
  export STORAGE="/data/atlasfs02/a/users/chekanov/Dijet2021lepton/data_wjets/"
fi

if [[ $HH =~ .*atlas24* ]]; then
  echo "HEP ANL ROOT setup"
  # source /users/admin/share/sl7/setup.sh
  source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos9-gcc11-opt/setup.sh
  export STORAGE="/data/atlasfs02/a/users/chekanov/Dijet2021lepton/data_wjets/"
fi


if [[ $HH =~ .*lcrc.anl.* ]]; then
  echo "LCRC ANL ROOT setup"
  # source /soft/hep/hep_setup.sh
  # source $HOME/soft/setup.sh 
  source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos8-gcc11-opt/setup.sh
  # source ./setup_cpp.sh
  #export STORAGE="/lcrc/group/ATLAS/atlasfs/local/chekanov/Dijet2021lepton/data_wjets/"
  # use this for systematic 94 (data) for highPt muon working point
  export STORAGE="/lcrc/group/ATLAS/atlasfs/local/chekanov/AnomalyDetect23/data/"
  export STORAGE_MC="/lcrc/group/ATLAS/atlasfs/local/chekanov/AnomalyDetect23/mc/"
  source /lcrc/group/ATLAS/users/software/el8/hep_texlive.sh
fi


FANN=./lib/fann/src/
export LD_LIBRARY_PATH=$FANN:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./lib/src/:$LD_LIBRARY_PATH

echo "STORAGE=$STORAGE"
