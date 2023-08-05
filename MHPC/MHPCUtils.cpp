#include "MHPCUtils.h"

void MHPCVisualization::publishWBTrajectory(const MHPCProblemData<double>* problem_data)
{
    const auto& wb_trajs = problem_data->wb_trajs;
    
    // Clear wbtraj_lcmt in case it is not empty
    if  (wbtraj_lcmt.pos.size() > 0)
    {
        wbtraj_lcmt.pos.clear();
        wbtraj_lcmt.eul.clear();    
        wbtraj_lcmt.qJ.clear();
        wbtraj_lcmt.vWorld.clear();
        wbtraj_lcmt.eulrate.clear();
        wbtraj_lcmt.qJd.clear();
        wbtraj_lcmt.torque.clear();
        wbtraj_lcmt.sz = 0;
    }
    
    std::vector<double> pos_vec(3), eul_vec(3), vWorld_vec(3), eulrate_vec(3);
    std::vector<double> qJ_vec(12), qJd_vec(12), torque_vec(12);

    for (const auto& traj: wb_trajs)
    {
        for (size_t k = 0; k < traj->size() - 1; k++)
        {
            const auto& xk = traj->Xbar[k];
            const auto& uk = traj->Ubar[k];

            std::copy(xk.data(), xk.data() + 3, pos_vec.data());
            std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());
            std::copy(xk.data() + 6, xk.data() + 18, qJ_vec.data());
            std::copy(xk.data() + 18, xk.data() + 21, vWorld_vec.data());
            std::copy(xk.data() + 21, xk.data() + 24, eulrate_vec.data());
            std::copy(xk.data() + 24, xk.data() + 36, qJd_vec.data());
            std::copy(uk.data(), uk.data() + 12, torque_vec.data());

            wbtraj_lcmt.pos.push_back(pos_vec);
            wbtraj_lcmt.eul.push_back(eul_vec);
            wbtraj_lcmt.qJ.push_back(qJ_vec);            
            wbtraj_lcmt.vWorld.push_back(vWorld_vec);
            wbtraj_lcmt.eulrate.push_back(eulrate_vec);  
            wbtraj_lcmt.qJd.push_back(qJd_vec);                      
            wbtraj_lcmt.torque.push_back(torque_vec);
            wbtraj_lcmt.sz ++;
        }
        
    }        
    viz_lcm.publish("visualize_wb_traj", &wbtraj_lcmt);
    printf("Published a visualization_wb_traj lcm message \n");
}