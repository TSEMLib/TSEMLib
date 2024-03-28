#! /bin/bash


for i in {0..7..1}
do
    mesh2opendx 3 result/mesh$i.mesh result/mesh$i.dx
done
