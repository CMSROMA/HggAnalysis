#!/bin/tcsh


foreach dir (`rfdir $1 | awk '{print $9}'`)
    echo $dir
    rfdir $1/$dir | grep root | awk '{print $9}' | xargs -I file stager_get -M $1/$dir/file
end
