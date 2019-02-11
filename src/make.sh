#!/bin/sh
#$ -S /bin/sh
echo Start Time : 
date
echo g++	-g -O2 RRSelection.cpp		-lz  -L/usr/lib/ -L./include/zlib/	-o	../bin/RRSelection
g++	-g -O2 RRSelection.cpp 		-lz -L/usr/lib64/   -L/usr/lib/	-L./include/zlib/	-o	../bin/RRSelection
echo done see the [ ../bin/RRSelection ]
echo End Time : 
date
