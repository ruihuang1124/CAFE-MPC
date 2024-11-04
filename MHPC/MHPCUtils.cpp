#include "MHPCUtils.h"

void MHPCVisualization::publishWBTrajectory(const MHPCProblemData<double>* pdata,
                                            const MHPCConfig& config)
{
    std::vector<double> pos_vec(3), eul_vec(3), vWorld_vec(3), eulrate_vec(3);
    std::vector<double> qJ_vec(12), qJd_vec(12), pFoot_vec(12), torque_vec(12);
    std::vector<double> hg_vec(3), dhg_vec(3);
    std::vector<int> phase_contact(4);
    wbtraj_lcmt.sz = 0;
    wbtraj_lcmt.wb_sz = 0;
    wbtraj_lcmt.pos.clear();
    wbtraj_lcmt.vWorld.clear();
    wbtraj_lcmt.eul.clear();
    wbtraj_lcmt.eulrate.clear();
    wbtraj_lcmt.qJ.clear();
    wbtraj_lcmt.qJd.clear(); 
    wbtraj_lcmt.torque.clear();
    wbtraj_lcmt.contact.clear();
    wbtraj_lcmt.defect.clear();
    wbtraj_lcmt.hg.clear(); 
    wbtraj_lcmt.dhg.clear();
    wbtraj_lcmt.time.clear();
    double t = 0;
    for (int i(0); i < pdata->wb_phase_horizons.size(); i++)
    {
        const int h = pdata->wb_phase_horizons[i];
        const auto& tau = pdata->wb_trajs[i];
        std::copy(pdata->wb_phase_contacts[i].data(), pdata->wb_phase_contacts[i].data() + pdata->wb_phase_contacts[i].size(), phase_contact.data());
        for (int k = 0; k < h; k++)
        {
            wbtraj_lcmt.sz ++;
            wbtraj_lcmt.wb_sz ++;

            const auto& xk = tau->Xsim[k];
            const auto& uk = tau->U[k];
            const auto& defect_k = tau->Defect[k].lpNorm<Eigen::Infinity>();                       

            std::copy(xk.data(), xk.data() + 3, pos_vec.data());
            std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());
            std::copy(xk.data() + 6, xk.data() + 18, qJ_vec.data());
            std::copy(xk.data() + 18, xk.data() + 21, vWorld_vec.data());
            std::copy(xk.data() + 21, xk.data() + 24, eulrate_vec.data());
            std::copy(xk.data() + 24, xk.data() + 36, qJd_vec.data());
            std::copy(uk.data(), uk.data() + 12, torque_vec.data());

            wbtraj_lcmt.pos.push_back(pos_vec);
            wbtraj_lcmt.vWorld.push_back(vWorld_vec);
            wbtraj_lcmt.eul.push_back(eul_vec);
            wbtraj_lcmt.eulrate.push_back(eulrate_vec);
            wbtraj_lcmt.qJ.push_back(qJ_vec);
            wbtraj_lcmt.qJd.push_back(qJd_vec); 
            wbtraj_lcmt.torque.push_back(torque_vec);
            wbtraj_lcmt.contact.push_back(phase_contact);
            wbtraj_lcmt.defect.push_back(defect_k);
            wbtraj_lcmt.hg.push_back(hg_vec);
            wbtraj_lcmt.dhg.push_back(dhg_vec);
            wbtraj_lcmt.time.push_back(t);
            t+=config.dt_wb;
            
        }        
    }

    const int h = pdata->srb_phase_horizon;
    const auto& tau = pdata->srb_traj;

    std::cout << "srb phase horizon = " << h << "\n";
    hg_vec[0] = 0;hg_vec[2] = 0; hg_vec[2] = 0;    
    dhg_vec[0] = 0;dhg_vec[2] = 0; dhg_vec[2] = 0;    
    for (int k = 0; k < h; k++)
    {
        wbtraj_lcmt.sz ++;
        const auto& xk = tau->Xsim[k];
        const auto& uk = tau->U[k];
        const auto& defect_k = tau->Defect[k].lpNorm<Eigen::Infinity>();

        std::copy(xk.data(), xk.data() + 3, pos_vec.data());
        std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());        
        std::copy(xk.data() + 6, xk.data() + 9, vWorld_vec.data());
        std::copy(xk.data() + 9, xk.data() + 12, eulrate_vec.data());

        wbtraj_lcmt.pos.push_back(pos_vec);
        wbtraj_lcmt.vWorld.push_back(vWorld_vec);
        wbtraj_lcmt.eul.push_back(eul_vec);
        wbtraj_lcmt.eulrate.push_back(eulrate_vec);
        wbtraj_lcmt.qJ.push_back(qJ_vec);
        wbtraj_lcmt.qJd.push_back(qJd_vec);
        wbtraj_lcmt.torque.push_back(torque_vec);
        wbtraj_lcmt.contact.push_back(phase_contact);
        wbtraj_lcmt.defect.push_back(defect_k); 
        wbtraj_lcmt.hg.push_back(hg_vec);       
        wbtraj_lcmt.dhg.push_back(dhg_vec);
        wbtraj_lcmt.time.push_back(t);
        t += config.dt_srb;
    }      

    viz_lcm.publish("visualize_wb_traj", &wbtraj_lcmt);
    printf("Published a visualization_wb_traj lcm message \n");
}