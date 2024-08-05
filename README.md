## **Dependency**
This implementation uses [Eigen](https://gitlab.com/libeigen/eigen) for linear algegra, [Pinocchio](git@github.com:stack-of-tasks/pinocchio.git) for analytical derivatives, [LCM](https://github.com/lcm-proj/lcm/releases) for communications to low-level controllers, and [boost](https://www.boost.org/users/history/) for reading configuration files. A customized **Hybrid-Systems DDP (HS-DDP)** solver is employed to solve the nonlinear trajectory optimization problem. A C++ implementation of **HS-DDP** is included in this repo and can be found [here](https://github.com/heli-sudoo/HKD-MPC/tree/ICRA22%2BIROS23/MPC_Controller/HSDDPSolver).

### **Installation instructions for the third-party libraries**
- [Pinocchio 2.6.10](git@github.com:stack-of-tasks/pinocchio.git): Folow the instructions [here](https://stack-of-tasks.github.io/pinocchio/download.html) to install Pinocchio via robotpkg. As a reminder, do NOT forget configuring the environment variables, so that CMake can find Pinocchio while building.
- [LCM1.4.0](https://github.com/lcm-proj/lcm/releases/tag/v1.4.0), [Boost1.71](https://www.boost.org/users/history/), and [Eigen3.4](https://gitlab.com/libeigen/eigen/-/releases) could be installed by running the scripts

    ```
        ./install_dependencies.sh
    ```


- Note: By default, Pinocchio2 would install Eigen3.3. If you install Eigen3.4 first, Pinocchio2 will overwrite it. Therefore, it is suggested install Eigen3.4 after the installation of Pinocchio.



## **Build CAFE-MPC**

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

So far, you have built the CAFE-MPC controller. It has to function with a simulator and a whole-body controller.

## **Simulation and Value-Based Whole-Body Controller**
We use Cheetah-Software for dynamics simulation. The Value-Based Whole-body controller is implemented in Cheetah-Software as well. To buld Cheetah-Software, you need to install dependencies with the following script

```
./install_dependencies_cheetah_software.sh
```
Clone the Cheetah-Software repo:
```
cd ..
git clone -b -b cafempc_low_level https://github.com/heli-sudoo/Cheetah-Software-He.git Cheetah-Software
```
Build Cheetah-Software:
```
cd Cheetah-Software/scripts
./make_types.sh
cd .. 
mkdir build && cd build
cmake ..
make -j4
```

## **Run Simulator, VWBC, CAFE-MPC**
Open three terminals, one for simulation, one for low-level VWBC, and one for CAFE-MPC. 

In the **first terminal** 
```
cd Cheetah-Software/build
sim/sim
```
This opens two windows, one for simulation, the other one is a control panel.

In the **second terminal**,
```
user/MHPC_LLController/mhpc_llctrl m s
```
You will see the robot moves its legs to a zero configuration and stand up.

In the **third terminal**,
```bash
cd CAFE-MPC/build
MHPC/mhpc_run
```
In the control panel, switch the `contrl_mode` to 2. You will see the robot starts performing locomotion.

