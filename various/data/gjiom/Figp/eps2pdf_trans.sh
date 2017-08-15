#!/bin/sh

epslist=`more tt`

for i in ${epslist[@]}
do

convert $i.eps $i.png
convert $i.png $i.pdf
done


