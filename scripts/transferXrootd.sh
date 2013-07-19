#!/bin/bash

usage(){
    echo "Usage: `basename $0` sample.list output_dir parallel_streams" > /dev/stderr
}

case $# in
    3)
	;;
    *)
	echo -n "Error: "
	usage
	exit 1
	;;
esac

function rfcopy { 
    (rfcp $1/$3/$file $2/$3/$file)&
    pid=$!
    wait $pid
    status=$?
    if [ $status -ne 0 ]; then
	echo "$3/$file $status" >> $3_error.log
    else
	echo "$3/$file" >> $3_success.log
    fi
}

for dir in `rfdir $1 | awk '{print $9}'`; do
    echo $dir
    mkdir -p $2/$dir
    rfdir $1/$dir | grep root | awk '{print $9}' > /tmp/filelist.$dir.txt
    IFSOLD=$IFS
    IFS=$'\n'
    for file in `cat /tmp/filelist.$dir.txt`; do
	echo "Copying $dir/$file"
	while [ `jobs | wc -l` -ge $3 ];
	  do
#	  echo "Sleeping 5"
	  sleep 5 
	done 
	if [ ! -r "$2/$dir/$file" ]; then 
	    (rfcopy $1 $2 $dir $file)&
	fi
    done
    IFS=$IFSOLD
done
