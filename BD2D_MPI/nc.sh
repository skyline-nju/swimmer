#!/bin/bash
function check() {
    if [ ! -d $1 ]; then
        echo "Error, $1 not exists"
    else
        echo "$1"
    fi
}

echo "user: ${USER}"
if [ ${USER} == "nsyw449_YK" ]; then
    source /lsfhome/env/intel-2016.sh
    source /lsfhome/env/mvapich2-2.0-intel.sh
    export NCDIR=/home-yw/users/nsyw449_YK/dy/Program/netCDF-4.6.1-hdf5-1.8.20
    export H5DIR=/home-yw/users/nsyw449_YK/dy/Program/hdf5-1.8.20-mvapich2-2.0-intel
    export CURLDIR=/home-yw/users/nsyw449_YK/dy/Program/anaconda3/lib
    export CMDLINE_DIR=/home-yw/users/nsyw449_YK/dy/Program/include
elif [ ${USER} == "dy" ]; then
    export NCDIR=/cluster/tool/netcdf-4.6.1-mpicc-3.0.4
    export H5DIR=/cluster/tool/hdf5-1.8.20-mpicc-3.0.4
    export CURLDIR=/home/dy/local/anaconda3/lib
    export CMLINE_DIR=/home/dy/local/include
else
    echo "export nothing"
fi

check $NCDIR
check $H5DIR
check $CURLDIR
check $CMDLINE_DIR
