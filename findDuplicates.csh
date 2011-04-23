#!/bin/csh
#

if( $#argv == 0 ) then
  echo "usage: ./findDuplicates.csh  <directory containing crab output files> <name of the output with unique files>"
  exit 0
endif

set indir = $1 

#set files = ( `/bin/ls $indir` )
set files = ( `nsls $indir` )

set tmpfile = /tmp/${USER}/allfiles.crab
rm -Rf $tmpfile
touch  $tmpfile
rm -Rf /tmp/uniq.files /tmp/uniq.files.tmp /tmp/duplicates.files
touch  /tmp/uniq.files

foreach i ( $files )
  echo $i | awk 'BEGIN{FS="_"}{ file = $0; if($2<10) j="000"$2;else if($2<100) j = "00"$2;else if($2<1000) j = "0"$2;else j=$2;  print j":"file }' >> $tmpfile
end


sort $tmpfile | uniq -w 4 | awk -F: '{print $2}' > /tmp/uniq.files.tmp
sort $tmpfile | uniq -d -w 4 | awk -F: '{print $2}' > /tmp/duplicates.files

foreach i ( `cat /tmp/uniq.files.tmp`  ) 
  echo $indir/$i >> /tmp/uniq.files
end

wc -l  $tmpfile /tmp/{uniq.files,duplicates.files}
cp /tmp/uniq.files $2
