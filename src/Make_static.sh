#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2018-12-27
echo Start Time : 
date
g++	-g	-O2	*.cpp	-lz	-L/usr/lib64/	-L/usr/lib/	-L./include/zlib/	-o	../bin/RRSelection_static	-static	-lz	-I	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/include/		-L	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/lib	
echo End Time : 
date
