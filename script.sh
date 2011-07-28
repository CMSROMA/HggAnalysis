#!/bin/bash

. $HOME/.bash_profile

#export VO_CMS_SW_DIR=/swshare/cms
#. $VO_CMS_SW_DIR/cmsset_default.sh

cd /shome/meridian/software/CMSSW423/src/Analysis/Higgs

#echo "SCRAM_ARCH = $SCRAM_ARCH"
eval `scramv1 runtime -sh`

#env

echo dir is $CMSSW_BASE file is $1 $2 $3 $4 $5 $6

${CMSSW_BASE}/src/Analysis/Higgs/tmp/redntpApp $1 $2 $3 $4 $5 $6
