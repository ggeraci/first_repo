#!/bin/bash

cd data
./clear.sh
cd ../data_IO
./clear.sh
cd ../stat
rm *.dat
rm *.bin
cd ../stat_IO
rm *.dat
rm *.bin
cd ..