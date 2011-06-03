#!/bin/tcsh

set run = 0
set version = "v1"

if($#argv == 0 || $#argv < 2) then
  echo "usage:  mergeAll2011.csh <version> <run if 1>"
  echo "        version: version string for redntp"
  echo "        run: default=0  set to 1 to execute"
  exit 0
endif

set version = $1
echo "version : $version "

set run = $2
echo "run : $run "

foreach class ( 41xv10  41xv10_data 42xv1 42xv1_data ) 
#foreach class ( 42xv1_data ) 
    foreach preseltype ( preselection cicloose ) 
#    foreach preseltype ( cicsuper ) 
	set command="./mergeRedNtp.csh redntp.${class}.${preseltype}.${version} ${run}"
	echo ${command}
	if ( $run == 1 ) then
	   ${command}
	endif
    end
end
