#include "PinocchioInteface.h"
#include "MHPCConstraint.h"
#include "visualize_quadState_lcmt.hpp"
#include <lcm/lcm-cpp.hpp>

int main()
{
    WBM::Model<double>::CtactStatusType contact;
    contact << 1, 0, 0, 0;

    /* Build a whole-body model */
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";
    pinocchio::Model pin_model;
    buildPinModelFromURDF(urdf_filename, pin_model);
    WBM::Model<double> wb_model(pin_model);
    WBM::Model<double>::CtactStatusType tdStatus;

    /* Create constraints */
    tdStatus << 1, 0, 0, 0;
    MHPCConstraints::WBTouchDown<double> tdConstraint(tdStatus, &wb_model);
    MHPCConstraints::WBGRF<double> grfConstraint(contact);
    MHPCConstraints::TorqueLimit<double> torqueLimit;
    MHPCConstraints::JointLimit<double> jointLimit;

    /* Test kinematics */
    VecM<double, WBM::xs> x;
    VecM<double, 3> pos, eul, vel, eulrate;
    VecM<double, 12> qJ, qJd;
    VecM<double, 12> pFoot;

    pos.setZero();
    eul.setZero();
    vel.setZero();
    eulrate.setZero();
    pos[2] = 0.404;
    qJ.setZero();
    qJ.head<3>() << 3.14 / 2, 0, 3.14 / 2;
    x << pos, eul, qJ, vel, eulrate, qJd;

    wb_model.get_footPositions(pFoot, x);

    /* Publish robot state via LCM */
    lcm::LCM visualization_lcm;
    visualize_quadState_lcmt quad_state_lcm_data;

    // Check LCM initialization
    if (!visualization_lcm.good())
    {        
        printf("Failed to inialize mpc_lcm for hkd command\n");        
    }

    for (size_t i = 0; i < 3; i++)
    {
        quad_state_lcm_data.pos[i] = pos[i];
        quad_state_lcm_data.eul[i] = eul[i];
        quad_state_lcm_data.vWorld[i] = vel[i];
        quad_state_lcm_data.eulrate[i] = eulrate[i];
    }

    for (size_t i = 0; i < 12; i++)
    {
        quad_state_lcm_data.qJ[i] = qJ[i];
        quad_state_lcm_data.pFoot[i] = pFoot[i];
    }
    visualization_lcm.publish("visualize_mc", &quad_state_lcm_data);
}