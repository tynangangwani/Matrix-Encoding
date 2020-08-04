#!/bin/sh
mpirun -np 17 python3 MatrixMultErrorDetection17.py
mpirun -np 33 python3 MatrixMultErrorDetection33.py
mpirun -np 65 python3 MatrixMultErrorDetection65.py
mpirun -np 129 python3 MatrixMultErrorDetection129.py
