#include"WBM.h"
#include"SRBM.h"
#include"PinocchioInteface.h"
#include"QuadReference.h"
#include"MHPCFootStep.h"

int main()
{
    /* Create quadruped reference */
    std::string reference_file_path = "../Reference/Data/quad_reference.csv";       
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    quad_ref->load_top_level_data(reference_file_path, true);
    quad_ref->initialize(0.1);

    /* Build a pinocchio model and a wbm */
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";
    pinocchio::Model pin_model;
    buildPinModelFromURDF(urdf_filename, pin_model);
    WBM::Model<double> wb_model(pin_model);

    /* Create the footstepplanner */
    MHPCFootStep<double> footStepPlanner(&wb_model, quad_ref.get());

    /* Build and initialize the srb model */
    SRBM::Model<double> srb;
    srb.set_footStepPlanner(&footStepPlanner);


    VecM<double, SRBM::xs> x, xdot;
    VecM<double, SRBM::us> u;
    VecM<double, SRBM::ys> y;
    MatMN<double, SRBM::xs, SRBM::xs> A;
    MatMN<double, SRBM::xs, SRBM::us> B;
    MatMN<double, SRBM::ys, SRBM::xs> C;
    MatMN<double, SRBM::ys, SRBM::us> D;
    x.setZero();
    u.setOnes();
    x[2] = 0.28;
    
    srb.dynamics(xdot,y,x,u,0.5,0.01);
    srb.dynamics_partial(A,B,C,D,x,u,0.5,0.01);

    VecM<double, 12> EE_pos;
    footStepPlanner.getFootPositions(EE_pos);

    std::cout << "xnext = \n" << xdot.transpose() << "\n";
    std::cout << "Foot position = \n" << EE_pos.transpose() << "\n";
    std::cout << "A=\n" << A << "\n";
    std::cout << "B=\n" << B << "\n";
}