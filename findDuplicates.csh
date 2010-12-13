#!/bin/csh
#

if( $#argv == 0 ) then
  echo "usage: ./findDuplicates.csh  <directory containing crab output files>"
  exit 0
endif

set indir = $1

set files = ( `/bin/ls $indir` )

set tmpfile = /tmp/allfiles.crab
rm -Rf $tmpfile

foreach i ( $files )
  echo $i | awk 'BEGIN{FS="_"}{ file = $0; if($2<10) j="000"$2;else if($2<100) j = "00"$2;else if($2<1000) j = "0"$2;else j=$2;  print j":"file }' >> $tmpfile
end

sort $tmpfile | uniq -w 4 | awk '{FS=":"}{print $2}' > /tmp/uniq.files
sort $tmpfile | uniq -d -w 4 | awk '{FS=":"}{print $2}' > /tmp/duplicates.files


wc -l  $tmpfile /tmp/{uniq.files,duplicates.files}

