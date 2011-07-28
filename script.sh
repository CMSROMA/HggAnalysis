#!/bin/bash

domain=`dnsdomainname`
redirector=pccmsrm23.cern.ch:1094//u2/xrootd

if [ "$domain" == "psi.ch" ]; then
    . $HOME/.bash_profile
fi

#export VO_CMS_SW_DIR=/swshare/cms
#. $VO_CMS_SW_DIR/cmsset_default.sh
cd $1
#echo "SCRAM_ARCH = $SCRAM_ARCH"
eval `scramv1 runtime -sh`

if [ "$domain" == "cern.ch" ]; then
    castordir=`dirname $3`
    filename=`basename $3`
    rfmkdir ${castordir}
    XROOTLIB=`scram tool info xrootd | grep LIBDIR | awk -F "=" '{print $2}'`
    xrootdir=`echo $castordir | awk -F '/' '{for (i=NF-3; i<=NF; i++) { printf "%s/",$i};}'`
    export LD_PRELOAD=${XROOTLIB}/libXrdPosixPreload.so 
    echo "Creating dir  xroot://${redirector}/${xrootdir}" 
    mkdir xroot://${redirector}/${xrootdir}
#    env --unset=LD_PRELOAD
    echo "Output ${filename} will be copied in directory ${castordir} and in xroot://${redirector}/${xrootdir}"
else
    filename=$2
fi

#env
if [ "$domain" == "cern.ch" ]; then
    cd -
fi
echo dir is $CMSSW_BASE file is $1 $2 $3 $4 $5 $6 $7

${CMSSW_BASE}/src/Analysis/Higgs/tmp/redntpApp $2 ${filename} $4 $5 $6 $7

if [ "$domain" == "cern.ch" ]; then
    rfcp ${filename} ${castordir}/${filename}
    export LD_PRELOAD=${XROOTLIB}/libXrdPosixPreload.so 
    cp ${filename} xroot://${redirector}/${xrootdir}/${filename}
fi