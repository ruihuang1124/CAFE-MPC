#include "MHPCReset.h"
#include "WBM.h"
#include "PinocchioInteface.h"
#include "MHPCFootStep.h"
#include "QuadReference.h"

int main()
{
    // Build a whole-body model
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";
    pinocchio::Model pin_model;
    buildPinModelFromURDF(urdf_filename, pin_model);
    WBM::Model<double> wb_model(pin_model);
    
    // Build a quadruped reference model
    std::string reference_file_path = "../Reference/Data/quad_reference.csv";       
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    quad_ref->load_top_level_data(reference_file_path);
    quad_ref->initialize(0.5);

    // Create a foot step planner
    MHPCFootStep<double> footStepPlanner(&wb_model, quad_ref.get());
    MHPCReset<double> reset(&wb_model, &footStepPlanner);

    // Execute reset
    DVec<double> x_wb(WBM::xs);
    DVec<double> x_srb(SRBM::xs);
    VecM<double, 3> pos, eul, vel, eulrate;
    VecM<double, 12> qJ, qJd;
    pos.setOnes();
    eul.setOnes();
    vel.setConstant(0.5);
    eulrate.setConstant(0.1);
    qJ.setZero();
    qJd.setZero();    
    x_wb << pos, eul, qJ, vel, eulrate, qJd;
    VecM<int, 4> contact;
    VecM<int, 4> contact_next;
    ModelType mtype_cur = ModelType::WB;
    ModelType mtype_next = ModelType::SRB;
    reset.reset_map(x_srb, x_wb, contact, contact_next, mtype_cur, mtype_next);
    std::cout << "x_wb = \n" << x_wb.transpose() << "\n";
    std::cout << "x_srb = \n" << x_srb.transpose() << "\n";
}