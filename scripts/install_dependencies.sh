#/bin/bash

# Install boost
sudo apt-get update
sudo apt-get install libboost1.71-all-dev

# Install LCM
mkdir external
cd external
wget https://github.com/lcm-proj/lcm/releases/download/v1.4.0/lcm-1.4.0.zip
unzip lcm-1.4.0.zip
rm -f lcm-1.4.0.zip
cd lcm-1.4.0
mkdir build
cd build
cmake ..
make
sudo make install
cd ../..

# Install Eigen
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
unzip eigen-3.4.0.zip
rm -f eigen-3.4.0.zip
cd eigen-3.4.0
mkdir build
cd build && cmake ..
make
sudo make install
cd ../../..

# OpenMP
sudo apt-get install --reinstall openmpi-bin libopenmpi-dev


