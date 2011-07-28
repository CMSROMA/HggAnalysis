#!/bin/csh
# $Id: runPrepareList.csh,v 1.3 2011/07/07 08:54:48 meridian Exp $

if( $#argv<3  ) then
  echo "usage:  runPrepareList.csh  <list dir>  <directory> <location>   [run if 1]"
  exit 0
endif

set prepareListCommand = preparelist_eth.csh

set run = 0

if( $#argv>2 ) then
  set run = $4
endif

set listdir = $1

set srmdir = "$2"

set location = $3

if ( "$location" == "cern" ) then
    set lsCommand="rfdir"
else if ( "$location" == "xrootd" ) then
    set lsCommand="find"
else if ( "$location" == "eth" ) then
    set lsCommand="lcg-ls"
endif 

rm -rf ${listdir}/allFiles.txt
mkdir -p ${listdir}
touch ${listdir}/allFiles.txt
if ($location != "xrootd") then 
    ${lsCommand} "${srmdir}" | awk -F '/' '{print $NF}' | xargs -I {} ${lsCommand} "${srmdir}/{}" >> ${listdir}/allFiles.txt
else
    ${lsCommand} "${srmdir}" -type f >> ${listdir}/allFiles.txt
endif

cd ${listdir}/

if ($location != "xrootd") then 
    ${lsCommand} "${srmdir}" | awk -F '/' '{print $NF}' | xargs -I {} ../${prepareListCommand} allFiles.txt {} ${location} ${run} >! makeLists.log
else    
    ${lsCommand} "${srmdir}" -type d | awk -F '/' '{print $NF}' | xargs -I {} ../${prepareListCommand} allFiles.txt {}  ${location} ${run} >! makeLists.log
endif

if( $run != 1 ) then
  rm -rf ${listdir}
endif
