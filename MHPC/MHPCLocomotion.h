#ifndef MHPCLOCOMOTION_H
#define MHPCLOCOMOTION_H

#include <thread>
#include <mutex>
#include <lcm/lcm-cpp.hpp>

#include "utilities.h"
#include "solver_info_lcmt.hpp"

#include "MHPC_Command_lcmt.hpp"
#include "MHPC_Data_lcmt.hpp"
#include "MHPCProblem.h"
#include "MultiPhaseDDP.h"
#include "MHPCUtils.h"

template<typename T>
class MHPCLocomotion
{
public:
    MHPCLocomotion() : mpc_lcm(getLcmUrl(255)),  
                       solve_time(0)
    {                          
        // Check LCM initialization
        if (!mpc_lcm.good())
        {
            printf(RED);
            printf("Failed to inialize lcm for MHPC\n");
            return;
        }       

        is_first_mpc = true;

        // Subscribe to MHPC channel
        mpc_lcm.subscribe("MHPC_DATA", &MHPCLocomotion::mpcdata_lcm_handler, this);               
    }
    void mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                             const MHPC_Data_lcmt *msg);
    void publish_mpc_cmd();
    void initialize();
    void update();    
    void run(){
        while (mpc_lcm.handle()==0){}
    }

public:
    // MPC
    MHPCProblem<T> opt_problem;
    MHPCProblemData<T> opt_problem_data;
    MHPCConfig mpc_config;
    HSDDP_OPTION ddp_setting;    

    T dt_mpc;
    T mpc_time;
    T mpc_time_prev;
    int mpc_iter;
    bool is_first_mpc;
    
    DVec<T> xinit;
    VecM<T, WBM::xs> x_init_wb;    

    VecM<T, WBM::nu> qJ, qJd;
    VecM<T, WBM::nu> torque;
    VecM<T, 3> pos, eul, vWorld, eulrate; // Euler angle in local cooridinates in yaw-pitch-roll order

    // LCM message
    MHPC_Data_lcmt mpc_data;
    MHPC_Command_lcmt mpc_cmd;
    
    lcm::LCM mpc_lcm;

    // mutex lock
    std::mutex mpc_mutex;    

    // solve time
    float solve_time;

    MHPCVisualization mhpc_viz;    
};



#endif