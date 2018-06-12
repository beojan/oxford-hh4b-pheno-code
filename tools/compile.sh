#!/bin/sh
source ../setupRoot.sh
g++ -O2 `root-config --cflags` -o csv2root csv2root.cpp -lstdc++fs `root-config --libs`

