#!/bin/tcsh

set data_json = "/shome/meridian/software/CMSSW423/src/Analysis/Higgs/jsonFiles/2010_Dec22_Golden_May10_withMissingEBMinus17_PromptCertified_UpTo_165525_DCSOnly_Upto165620.json"
set puweight_41x = "/shome/meridian/software/CMSSW423/src/Analysis/Higgs/mc_41x_PUweight.root"
set puweight_42x = "/shome/meridian/software/CMSSW423/src/Analysis/Higgs/mc_42x_PUweight.root"
set ptweightfile_template = "/shome/meridian/software/CMSSW423/src/Analysis/Higgs/kfactors/Kfactors_MASSVALUE_AllScales.root"

set location = "eth"
set version = "v1"
set run = 0

if($#argv == 0 || $#argv < 3 || $#argv > 6 ) then
  echo "usage:  runAllAnalysis2011.csh  <location> <version> <run if 1> <jsonfile> <pureweight> <ptweight>"
  echo "        locations: cern roma eth"
  echo "        version: version string for redntp"
  echo "        run: default=0  set to 1 to execute"
  echo "        jsonfile: optional json to select good RUN_LS used for data (full path is required)"
  echo "        pu weight for MC: default=-1  set to 1 to store them"
  echo "        pt weight for MC (HiggsPt): default=-1  set to 1 to store them"
  exit -1
endif

set location = $1
echo "location : $location "

set version = $2
echo "version : $version "

set run = $3
echo "run : $run "

if ($#argv > 3) then
  set data_json = $4
  echo "json : ${data_json}"
endif 

set puweight = -1
if ($#argv > 4) then
  set puweight = $5
  echo "pu weight: ${puweight}"
endif 

set ptweight = -1
if ($#argv > 5) then
  set ptweight = $6
  echo "pt weight: ${ptweight}"
endif 

foreach class ( 41xv10  41xv10_data 42xv2 42xv1_data ) 
#foreach class ( 42xv2 ) 
    foreach preseltype ( preselectionCS cicloose ) 
#    foreach preseltype ( preselectionCS ) 
	if ( "`echo ${class} | grep data`XXX" != "XXX" ) then
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${version} ${preseltype} ${location} ${run} $data_json -1 -1"
	else 
	    if ( $puweight !=  -1 ) then
		if ( "`echo ${class} | grep 41x`XXX" != "XXX" ) then
		    set puweightFile = ${puweight_41x}
		else if ( "`echo ${class} | grep 42x`XXX" != "XXX" ) then
			set puweightFile = ${puweight_42x}
		endif
	    else
		set puweightFile = -1
	    endif

	    if ( $ptweight !=  -1 ) then
		set ptweightFile = ${ptweightfile_template}
	    else
		set ptweightFile = -1
	    endif
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${version} ${preseltype} ${location} ${run} -1 ${puweightFile} ${ptweightFile}"
	endif
	echo ${command}
	if ( $run == 1 ) then
	   ${command}
	endif
    end
end
