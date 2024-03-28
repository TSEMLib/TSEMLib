#! /bin/bash

for polynomial_order in {4..10}
do
    ./run.sh.sample cube.mesh $polynomial_order 0 2>/dev/null
done
