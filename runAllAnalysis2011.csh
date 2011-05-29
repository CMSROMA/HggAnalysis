#!/bin/tcsh

set data_json = "/shome/meridian/software/CMSSW423/src/Analysis/Higgs/jsonFiles/2010_Dec22_Golden_May10_withMissingEBMinus17_PromptCertified_UpTo_165525_DCSOnly_Upto165620.json"
set location = "eth"
set version = "v1"
set run = 0

if($#argv == 0 || $#argv < 3 || $#argv > 4 ) then
  echo "usage:  runAllAnalysis2011.csh  <location> <version> <run if 1> <jsonfile>"
  echo "        locations: cern roma eth"
  echo "        version: version string for redntp"
  echo "        run: default=0  set to 1 to execute"
  echo "        jsonfile: optional json to select good RUN_LS used for data (full path is required)"
  exit 0
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

foreach class ( 41xv10  41xv10_data 42xv1 42xv1_data ) 
    foreach preseltype ( preselection cicloose ) 
	if ( "`echo ${class} | grep data`XXX" != "XXX" ) then
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${version} ${preseltype} ${location} ${run} $data_json"
	else
	    set command="./makeRedNtp.csh list.${class}/ redntp.${class}.${preseltype}.${version} ${preseltype} ${location} ${run}"
	endif
	echo ${command}
	if ( $run == 1 ) then
	   ${command}
	endif
    end
end
