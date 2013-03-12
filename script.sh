#!/bin/bash

domain=`dnsdomainname`


if [ "$domain" == "psi.ch" ]; then
    . $HOME/.bash_profile
fi

#export VO_CMS_SW_DIR=/swshare/cms
#. $VO_CMS_SW_DIR/cmsset_default.sh
cd $1
#echo "SCRAM_ARCH = $SCRAM_ARCH"
eval `scramv1 runtime -sh`

if [ "$domain" == "cern.ch" ]; then
    redirector=pccmsrm27.cern.ch
    castordir=`dirname $3`
    filename=`basename $3`
    #rfmkdir ${castordir}
 #   XROOTLIB=`scram tool info xrootd | grep LIBDIR | awk -F "=" '{print $2}'`
    xrootdir=/cms/local/`echo $castordir | awk -F '/' '{for (i=NF-3; i<=NF; i++) { printf "%s/",$i};}'`
 #   export LD_PRELOAD=${XROOTLIB}/libXrdPosixPreload.so 
    echo "Creating dir root://${redirector}//${xrootdir}" 
    xrd mkdir ${redirector} ${xrootdir}
#    unset LD_PRELOAD
    echo "Output ${filename} will be copied in root://${redirector}/${xrootdir}"
else
    filename=$2
fi

#env
if [ "$domain" == "cern.ch" ]; then
    cd -
fi
echo dir is $CMSSW_BASE file is $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
echo ${CMSSW_BASE}/src/Analysis/Higgs/tmp/redntpApp $2 ${filename} $4 ${11} ${12} ${13} $5 $6 $7 $8 $9 ${10}

${CMSSW_BASE}/src/Analysis/Higgs/tmp/redntpApp $2 ${filename} $4 ${11} ${12} ${13} $5 $6 $7 $8 $9 ${10}
exit_stat=$?

if [ "$domain" == "cern.ch" ]; then
    if [ ${exit_stat} != 0 ]; then
	echo `date` ${xrootdir}/${filename} ${exit_stat} >> $1/log/runerror.jobs
	exit ${exit_stat}
    else
	echo `date` ${xrootdir}/${filename} >> $1/log/runsuccess.jobs
    fi
fi

if [ "$domain" == "cern.ch" ]; then
#    rfcp ${filename} ${castordir}/${filename}
#    exit_stat=$?
#    export LD_PRELOAD=${XROOTLIB}/libXrdPosixPreload.so 
    xrdcp ${filename} root://${redirector}//${xrootdir}/${filename}
    exit_stat1=$?
#    if [ ${exit_stat} != 0 ]; then
#	echo `date` ${xrootdir}/${filename} ${exit_stat} >> $1/log/castorcopyerror.jobs
    if [ ${exit_stat1} != 0 ]; then
	echo `date` ${xrootdir}/${filename} ${exit_stat1} >> $1/log/xrootcopyerror.jobs
    else
	echo `date` ${xrootdir}/${filename} >> $1/log/copysuccess.jobs
    fi
fi
