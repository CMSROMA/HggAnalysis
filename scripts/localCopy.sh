#!/bin/bash
#PREFIX=Higgs

usage(){
    echo "Usage: `basename $0` ./filelist/sample.list" > /dev/stderr
}

case $# in
    0)
	echo -n "Error: "
	usage
	exit 1
	;;
    *)
	;;
esac

for i in $@; do
    echo -n "Check file $i: "
    if [ ! -r "$i" ]; then 
	echo "File $i not found or not readable" 
	exit 1
    fi
    echo
done

for ciao in $@; do
    filelist=(`cat $ciao`)
    sample=`basename $ciao .list`
    
# creo la cartella
    if [ ! -d "/xrootdfs/u2/xrootd/" ]; then 
	dir=/tmp/${USER}/$sample
    else
	dir=/xrootdfs/u2/xrootd/meridian/Higgs/reduced_bck/$sample
    fi

    if [ ! -d "$dir" ]; then
	mkdir $dir -p
    fi
    
    echo "Copy $ciao"
    for i in ${filelist[@]}; do
	echo "- copy: $i"
	while [ `jobs | wc -l` -gt 5 ];
	  do
#	  echo "Sleeping 5"
	  sleep 5 
	done 

	file=`basename $i`
	if [ ! -r "$dir/$file" ]; then 
	    (rfcp $i $dir)&
	fi
    done
done

wait
