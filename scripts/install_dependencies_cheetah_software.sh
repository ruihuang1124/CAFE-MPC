#/bin/bash

# Dependencies of Cheetah Software
sudo apt install mesa-common-dev freeglut3-dev libblas-dev liblapack-dev gfortran  build-essential libglib2.0-dev

# Install qt5
# By default, Ubuntu 20.04.3 runs version 5.12.8 of qt5
sudo apt update
sudo apt install qt5-default libqt5gamepad5-dev