
sudo apt install libeigen3-dev libpcl1.7

cd ppe/cpp/
mkdir build
cd build/
cmake ..
make -j8

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src
#cd ../mex/
#ln -s ../build/src/*.so ./

cd ../mex/
matlab
mex -I../include -L../build/src -linc_test pcextrplnc.c
a = pcextrplnc([1, 2, 4])
