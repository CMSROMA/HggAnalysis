#!/bin/csh

if( $#argv<2  ) then
  echo "usage:  prepareList.csh  <valid directory>   <listname>   [run if 1]"
  exit 0
endif

set run = 0
if( $#argv>2 ) then
  set run = $3
endif

set indir = $1
if( ! -d $indir ) then
  echo "invalid directory <$indir>"
  exit -1
endif

set listname = $2

# num of files per list file
set filexlist  = 50

set prepend="dcap://cmsrm-se01.roma1.infn.it"
 
set files = ( `/bin/ls $indir` )

echo "# of root files in directory: $#files"
echo "# of files per list: $filexlist"

set tmpfile = /tmp/tmpfile
rm -Rf $tmpfile
touch $tmpfile

foreach i ( $files )
  echo "$prepend/$indir/$i" >> $tmpfile
end

set suffixlen = 2

split -l $filexlist -d -a $suffixlen  $tmpfile  ${listname}_

foreach i ( ${listname}_?? )
  mv $i ${i}.list
  echo "new list:   ${i}.list"
end

if( $run != 1 ) then
  rm -Rf ${listname}_???.list
endif

rm -Rf $tmpfile
