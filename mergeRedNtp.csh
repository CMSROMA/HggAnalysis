#!/bin/csh
## $Id: $

if($#argv == 0 || $#argv < 2) then
  echo "usage:  mergeRedNtp.csh  <indir>  <run if 1>"
  echo "       output will be in  indir/merged"
  exit 0
endif

set indir = $1
echo "merging files in $indir"

set run = 0
if($#argv > 1)  set run = $2

# output dir
set outdir = $indir/merged
echo "creating $outdir for output"
if($run == 1) mkdir  $outdir

## list of samples to merge
set samples  = ( `/bin/ls -1 $indir | awk 'BEGIN{FS="."}{print  $1}' | awk 'BEGIN{FS="_[0-9][0-9]$"}{print  $1}' | grep -v merged | uniq ` )

foreach i ( $samples )
  echo "-----------  sample: $i"
  set outname = ${i}.root
  set outfile = $outdir/$outname
  
  # remove old merged file if exists
  if(-f $outfile)  then
      echo "removing old merged file $outfile"
      if($run == 1) rm $outfile
  endif

  # hadd command
  set comm = "hadd  $outfile $indir/${i}_??.root"
  if($run == 1) then
     $comm
  else 
    echo "$comm"
  endif
  echo "created $outname in $outdir"
end

echo "list of all merged files in $outdir"
/bin/ls -1 $outdir


