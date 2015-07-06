#! /bin/bash

mpiexec -np 64 ./box params1.in
mkdir NG1
mv *.dat NG1
mv *.bin NG1
mv info.txt NG1
mv NG1/U.bin ./
mv NG1/V.bin ./
mv NG1/W.bin ./

mpiexec -np 64 ./box params2.in
mkdir NG2
mv *.dat NG2
mv *.bin NG2
mv info.txt NG2
mv NG2/U.bin ./
mv NG2/V.bin ./
mv NG2/W.bin ./

mpiexec -np 64 ./box params3.in
mkdir NG3
mv *.dat NG3
mv *.bin NG3
mv info.txt NG3
mv NG3/U.bin ./
mv NG3/V.bin ./
mv NG3/W.bin ./

mpiexec -np 64 ./box params4.in
mkdir NG1
mv *.dat NG4
mv *.bin NG4
mv info.txt NG4
mv NG4/U.bin ./
mv NG4/V.bin ./
mv NG4/W.bin ./

mpiexec -np 64 ./box params.in
mkdir Cooling_long
mv *.dat Cooling_long
mv *.bin Cooling_long
mv info.txt Cooling_long
mv Cooling_long/U.bin ./
mv Cooling_long/V.bin ./
mv Cooling_long/W.bin ./
