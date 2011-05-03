#!/bin/tcsh

# change if needed
set castordir = /castor/cern.ch/user/d/delre/reduced/

set preselections = ( looseeg  tighteg  hggtighteg looseegpu  tightegpu  hggtightegpu isem superloose loose medium mcass )


if($#argv == 0 || $#argv < 5) then
  echo "usage:  makeRedNtp.csh  <listdir>  <outdir>  <pre-selection>  <location>   <run if 1>"
  echo "        listdir: valid directory containing list files"
  echo "        outdir: will be created in current directory in rome or on castor at cern"
  echo "                check |castordir| at the beginning of script"
  echo "        preselection:  $preselections"
  echo "        locations: cern roma"
  echo "        run: default=0  set to 1 to execute"
  exit 0
endif

# submit the job only if the 2nd argument is 1
set listdir  = list.V27.41x 
if ($#argv > 0) then
  set listdir = $1
  echo "listdir : $listdir "
  if(! -d $listdir ) then
    echo "<$listdir> does not exist... check again"
    exit -1
  endif
endif 

set outdir = ""
if ($#argv > 1) then
  set outdir = $2
  echo "outdir : $outdir "
endif 

set selection = ""
if ($#argv > 2) then
  set selection = $3
  set found = 0
  foreach i  ( $preselections )
    if( $selection == $i ) set found = 1
  end
  if($found == 0) then
    echo "bad preselection <$selection>. Choose from $preselections"
    exit -1
  end
  echo "selection : $selection "
endif 

# location: roma or cern
set location = ""
if ($#argv > 3) then
  set location = $4
  echo "location : $location "
  if( $location != roma && $location != cern ) then
    echo "bad location. options: roma or cern"
    exit -1
endif 

set run = 0
if ($#argv > 4) then
  set run = $5
  echo "run : $run "
endif 

echo "------   ready ro run at $location ------------------"

# logfiles always stored locally
set logdir = "./log/$outdir"
if($run == 1) mkdir -p $logdir

# choose queue, location based on location
if ($location == "cern" ) then
  set queue = 8nh
  set outdir = $castordir/$outdir
  set prefix = "rfio:"
  if($run == 1) rfmkdir $outdir
  echo "$outdir created on castor"
else if ($location == "roma" ) then
  set queue = "cmsshort"
  set outdir = ./$outdir
  set prefix = ""
  if($run == 1) mkdir -p $outdir
endif 

echo "queue : $queue "
echo "output files: $outdir"
echo "log files: $logdir"
echo "pre-selection: $selection"

# app to run
set app = ./tmp/redntpApp

if(! -e $app ) then
  echo "missing executable $app"
  exit 0
endif 

#set samples = ( $datasamples )

set samples  =  `/bin/ls -1 ${listdir} | awk 'BEGIN{FS="."}{print $1}'` 


foreach sample ( $samples )
   set rootfile = "${prefix}${outdir}/redntp_${sample}.root"
   set jobname = "${sample}"
   set logfile = "${logdir}/${sample}.txt"
   set listfile = "${listdir}/${sample}.list"
   if(! -e $listfile ) then
      echo "skipping non-existent file $listfile"
      continue
   endif
   set command = "bsub -q ${queue} -o $logfile -J ${jobname} cd ${PWD}; ${app} ${listdir}/${sample}.list ${rootfile} ${selection}"
   echo "---------------------------"
   echo "job name: ${jobname}"
   #echo " command: ${command}"
   if($run == 1) then
     ${command}
     #sleep 2 
   endif
end
