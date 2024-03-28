#! /bin/bash

mkdir result
./run.sh.sample cube.mesh 4 7 2>/dev/null
./mesh2dx.sh 2>/dev/null

