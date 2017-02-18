#!/bin/bash
Path=`pwd`

cd numerics/LUDe*
ln -s ../../../../src/numerics/LUDe*/* .
cd $Path

cd tensors/tensor
ln -s ../../../../src/tensors/tensor/* .
cd $Path

cd tensors/matrix
ln -s ../../../../src/tensors/matrix/* .
cd $Path

cd tensors/vector
ln -s ../../../../src/tensors/vector/* .
cd $Path

cd typedef
ln -s ../../../src/typedef/* .
cd $Path
