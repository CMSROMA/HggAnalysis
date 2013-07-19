#!/bin/tcsh

set run = 0

set indir = ""
set version = "v1"

if($#argv == 0 || $#argv < 3 || $#argv > 4) then
  echo "usage:  mergeAll2011.csh <version> <indir>  <run if 1> <energycorrections>"
  echo "        version: version string for redntp"
  echo "        indir: input directory (also on xrootd)"
  echo "        run: default=0  set to 1 to execute"
  exit 0
endif

set version = $1
echo "version : $version "

set indir = $2
echo "indir : $indir "

set run = $3
echo "run : $run "

set energyCorrection = -1
if ($#argv > 3) then
  set energyCorrection = $4
  echo "energyCorrection: ${energyCorrection}"
endif 

foreach class ( 42xv4 42xv4_data_new  ) 
#foreach class ( 42xv1_data ) 
    foreach preseltype ( preselectionCS cicloose) 
#    foreach preseltype ( cicloosenoeleveto ) 
#    foreach preseltype ( cicsuper ) 
	set command="./mergeRedNtp.csh ${indir}/redntp.${class}.${preseltype}.${energyCorrection}.${version} ${run}"
	echo ${command}
	${command}
    end
end
