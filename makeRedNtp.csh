#!/bin/tcsh
# $Id: makeRedNtp.csh,v 1.23 2011/06/03 15:12:44 meridian Exp $

# change if needed
set castordir = /castor/cern.ch/user/d/delre/reduced/

set preselections = ( looseeg  tighteg  hggtighteg looseegpu  tightegpu  hggtightegpu isem superloose loose medium cicloose cicmedium cictight cicsuper cichyper mcass preselection preselectionCS)


if($#argv == 0 || $#argv < 6 || $#argv > 7 ) then
  echo "usage:  makeRedNtp.csh  <inlist>  <outdir>  <pre-selection>  <location>  <run if 1> <jsonfile> <puweight>"
  echo "        inlist: valid directory containing list files OR valid list file"
  echo "        outdir: will be created in current directory in rome or on castor at cern"
  echo "                check |castordir| at the beginning of script"
  echo "        preselection:  $preselections"
  echo "        locations: cern roma eth"
  echo "        run: default=0  set to 1 to execute"
  echo "        jsonfile: optional json to select good RUN_LS"
  exit 0
endif

# submit the job only if the 2nd argument is 1
set listdir  = list.V27.41x 
if ($#argv > 0) then
  set listdir = $1
  #echo "listdir : $listdir "
  if(-f $listdir) then
       echo "<$listdir> is a single file"
  else if( -d $listdir ) then
       echo "<$listdir> is a directory"
    #echo "<$listdir> does not exist... check again"
   else 
    exit -1
  endif
endif 

set outdir = ""
if ($#argv > 1) then
  set outdir = $2
  echo "outdir : $outdir "
endif 

setenv selection  ""
if ($#argv > 2) then
  setenv selection  $3
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

# location: roma or cern or eth
set location = ""
if ($#argv > 3) then
  set location = $4
  echo "location : $location "
  if( $location != roma && $location != cern && $location != eth) then
    echo "bad location. options: roma or cern or eth"
    exit -1
endif 

set run = 0
if ($#argv > 4) then
  set run = $5
  echo "run : $run "
endif 

set json = 0
if ($#argv > 5) then
  set json = $6
  echo "json : $json "
endif 

set puweight = 0
if ($#argv > 6) then
  set puweight = $7
  echo "puweight : $puweight "
endif 

echo "------   ready ro run at $location ------------------"

# logfiles always stored locally
set logdir = "$PWD/log/$outdir"
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
else if ($location == "eth" ) then
  set queue = "short.q"
  set outdir = $outdir
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


## if listdir is one file only
if(-f $listdir) then
   set sample = `echo $listdir  | awk  'BEGIN{FS="/"}{print $NF}' | awk 'BEGIN{FS="."}{print $1}'`
   set listdir = `echo $listdir  | awk  'BEGIN{FS="/"}{print $1}'`
   setenv listfile $listdir
   setenv rootfile "${prefix}${outdir}/redntp_${sample}.root"
   set jobname = "${sample}"
   set logfile = "${logdir}/${sample}.txt"
   set logerrfile = "${logdir}/${sample}_err.txt"
   set listfile = "${listdir}/${sample}.list"
   if ($location == "cern" || $location == "roma") then  
     set command = "bsub -q ${queue} -o $logfile -J ${jobname} cd ${PWD}; ${app} ${listdir}/${sample}.list ${rootfile} ${selection} ${json} ${puweight}"
   else if ($location == "eth" ) then
     set command = "qsub -q ${queue} -o $logfile -e $logerrfile script.sh ${listfile} ${rootfile} ${selection} ${json} ${puweight}"
   endif  

   echo "---------------------------"
   echo "job name: ${jobname}"
   echo " command: ${command}"
   if($run == 1) then
     ${command}
     #sleep 2 
   endif

  ## if input is a directory
else if(-d $listdir) then
  #set samples  =  `/bin/ls -1 ${listdir} | awk 'BEGIN{FS="."}{print $1}'` 
  foreach i ( `/bin/ls -1 ${listdir} | grep ".list" | awk 'BEGIN{FS="."}{print $1}' | xargs  -I sample echo sample ` )
   set  sample = $i
   echo "sample : $sample"

   if( "`echo $sample | egrep -e "^G_Pt_"`XXX" != "XXX" ) then
     echo "skipping $sample"
     continue
   endif


   setenv rootfile "${prefix}${outdir}/redntp_${sample}.root"
   set jobname = "${sample}"
   set logfile = "${logdir}/${sample}.txt"
   set logerrfile = "${logdir}/${sample}_err.txt"
   setenv listfile  "${listdir}/${sample}.list"
   #if(! -e $listfile ) then
   #   echo "skipping non-existent file $listfile"
   #   continue
   #endif

   if ($location == "cern" || $location == "roma") then  
     set command = "bsub -q ${queue} -o $logfile -J ${jobname} cd ${PWD}; ${app} ${listdir}/${sample}.list ${rootfile} ${selection} ${json}"
   else if ($location == "eth" ) then
     set command = "qsub -q ${queue} -o $logfile -e $logerrfile script.sh ${listfile} ${rootfile} ${selection} ${json}"
   endif  

   echo "---------------------------"
   echo "job name: ${jobname}"
   echo " command: ${command}"
   if($run == 1) then
     ${command}
     #sleep 2 
   endif


  end

endif
