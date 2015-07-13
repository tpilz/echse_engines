#!/bin/bash

# Change to working directory
wd="/home/dkneis/cpp/dklib/hydro_snow/r_wrapper" 
cd $wd

# Remove old objects
old=$wd"/*.o"
ls $old > /dev/null 2>&1
if [ "$?" == "0" ]; then rm -I $old; fi
old=$wd"/*.so"
ls $old > /dev/null 2>&1
if [ "$?" == "0" ]; then rm -I $old; fi

# Create shared library (*.so)
R CMD SHLIB -lm snow_wrapper.cpp
