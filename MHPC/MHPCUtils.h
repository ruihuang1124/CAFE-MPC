#ifndef MHPC_UTILS_H
#define MHPC_UTILS_H

#include <lcm/lcm-cpp.hpp>

#include "utilities.h"
#include "MHPC-Trajopt/MHPCProblem.h"
#include "lcmtypes/cpp/wbTraj_lcmt.hpp"

class MHPCVisualization
{
public:
    MHPCVisualization() : viz_lcm(getLcmUrl(255))
    {
        if (!viz_lcm.good())
        {
            printf("Failed to initilize lcm for MHPC visualization \n");
        }
        
    }

    void publishWBTrajectory(const MHPCProblemData<double>* problem_data);

public:
    wbTraj_lcmt wbtraj_lcmt;
    lcm::LCM viz_lcm;
};



#endif