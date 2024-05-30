## **Dependency**
This implementation uses [Eigen](https://gitlab.com/libeigen/eigen) for linear algegra, [Pinocchio](git@github.com:stack-of-tasks/pinocchio.git) for analytical derivatives, [LCM](https://github.com/lcm-proj/lcm/releases) for communications to low-level controllers, and [boost](https://www.boost.org/users/history/) for reading configuration files. A customized **Hybrid-Systems DDP (HS-DDP)** solver is employed to solve the nonlinear trajectory optimization problem. A C++ implementation of **HS-DDP** is included in this repo and can be found [here](https://github.com/heli-sudoo/HKD-MPC/tree/ICRA22%2BIROS23/MPC_Controller/HSDDPSolver).

### **Installation instructions for the third-party libraries**
- [Pinocchio 2.6.10](git@github.com:stack-of-tasks/pinocchio.git): Folow the instructions [here](https://stack-of-tasks.github.io/pinocchio/download.html) to install Pinocchio via robotpkg. As a reminder, do NOT forget configuring the environment variables, so that CMake can find Pinocchio while building.
- [LCM1.4.0](https://github.com/lcm-proj/lcm/releases/tag/v1.4.0):
    
    Download lcm v1.4.0 and unzip to lcm
    ```
    cd lcm
    mkdir build
    cd build && cmake ..
    make
    sudo make install
    ```
- [Boost1.71](https://www.boost.org/users/history/)
    ```
    sudo apt-get update
    sudo apt-get install libboost1.71-all-dev
    ```
- [Eigen3.4](https://gitlab.com/libeigen/eigen/-/releases)

    Download the source code of [Eigen3.4](https://gitlab.com/libeigen/eigen/-/releases) and unzip it to Eigen3.4
    ```
    mkdir build
    cd build && cmake ..
    make
    sudo make install
    ```
- Note: By default, Pinocchio2 would install Eigen3.3. If you install Eigen3.4 first, Pinocchio2 will overwrite it. Therefore, it is suggested install Eigen3.4 after the installation of Pinocchio.



## **Build**

Once Eigen and LCM are successfully installed, generate necessary lcm types

```bash
cd scripts
./make_types.sh
```

To build the MPC controller

```bash
mkdir build && cd build
cmake ..
make -j4
```

To run the CAFE-MPC controller
```bash
cd build
MHPC/mhpc_run
```

## **Configure MPC**
To be done