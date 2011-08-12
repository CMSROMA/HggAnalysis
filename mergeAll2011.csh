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

foreach class ( 42xv3 42xv3_data ) 
#foreach class ( 42xv1_data ) 
    foreach preseltype ( preselectionCS ) 
#    foreach preseltype ( cicsuper ) 
	set command="./mergeRedNtp.csh redntp.${class}.${preseltype}.${version} ${run}"
	echo ${command}
	${command}
    end
end
