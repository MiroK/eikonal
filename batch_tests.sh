#!/bin/bash

ROOT=/home/miro3/Documents/Programming/Cpp/Eikonal/test_results

cd $ROOT/lin_geometric;
./run.sh;

cd $ROOT/lin_brent/digits_double;
./run.sh;

cd $ROOT/lin_brent/digits_double_2;
./run.sh;

cd $ROOT/lin_brent/digits_double_3;
./run.sh;

cd $ROOT/lin_newton/digits_double;
./run.sh;

cd $ROOT/lin_newton/digits_double_2;
./run.sh;

cd $ROOT/lin_newton/digits_double_3;
./run.sh;

cd $ROOT/hermite/corners;
./run.sh;

cd $ROOT/hermite/surface;
./run.sh;

cd $ROOT/hermite/distance;
./run.sh;


