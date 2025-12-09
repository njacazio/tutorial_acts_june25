#!/bin/bash

outfile=geom/iris_segm_june25/material-map.json
if [ ! -f $outfile ]; then
    wget https://cernbox.cern.ch/remote.php/dav/public-files/t8ZTuM5vZuehFkY/material-map.json -O $outfile
fi
