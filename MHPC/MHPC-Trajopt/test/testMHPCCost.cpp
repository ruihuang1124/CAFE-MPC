#include "MHPCCost.h"
#include "MHPCReference.h"

int main()
{     
    SRBReference srb_ref;

    std::string reference_file_path = "../Reference/Data/quad_reference.csv";       
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    quad_ref->load_top_level_data(reference_file_path, true);
    quad_ref->initialize(1.0);
    srb_ref.set_quadruped_reference(quad_ref);

    // contact    
    SRBTrackingCost<double> srb_cost;
    srb_cost.set_reference(&srb_ref);

    VecM<double, SRBM::xs> x;
    VecM<double, SRBM::us> u;
    VecM<double, SRBM::ys> y;
    RCostData<double,SRBM::xs,SRBM::us,SRBM::ys> rcostData;
    TCostData<double,SRBM::xs> tcostData;
    x.setOnes();
    u.setOnes();
    srb_cost.running_cost(rcostData, x, u, y, 0.01, 0.0);
    srb_cost.running_cost_par(rcostData, x, u, y, 0.01, 0.0);
    srb_cost.terminal_cost(tcostData, x, 0.0);
    srb_cost.terminal_cost_par(tcostData, x, 0.0);

    std::cout << "l = " << rcostData.l << "\n";
    std::cout << "lx = \n" << rcostData.lx.transpose() << "\n";
    std::cout << "lu = \n" << rcostData.lu.transpose() << "\n";
    std::cout << "lxx = \n" << rcostData.lxx << "\n";
    std::cout << "luu = \n" << rcostData.luu << "\n";
    std::cout << "Phi = \n" << tcostData.Phi << "\n";
    std::cout << "Phi_x = \n" << tcostData.Phix.transpose() << "\n";
    std::cout << "Phi_xx = \n" << tcostData.Phixx << "\n";
}