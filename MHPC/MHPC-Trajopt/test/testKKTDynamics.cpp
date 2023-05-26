
#include"PinocchioInteface.h"
#include"WBM.h"
#include <iomanip>

int main()
{
    // Build a whole-body model
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";
    pinocchio::Model pin_model;
    buildPinModelFromURDF(urdf_filename, pin_model);
    WBM::Model<double> wb_model(pin_model);
    wb_model.printModelInfo();

    // Set configuration, velocity, torque
    VecM<double, WBM::xs> xnext, x;
    VecM<double, WBM::ys> GRF;
    VecM<double, WBM::nq> q; 
    VecM<double, WBM::nv> qd;
    VecM<double, WBM::nu> u;
    srand((unsigned int) time(0));
    q.setRandom();
    qd.setRandom();    
    x << q, qd;
    u.setRandom();

    WBM::Model<double>::CtactStatusType contact, nextContact;
    contact << 1,1,0,1;
    nextContact << 1,1,1,1;
    float t(0), dt(0.005);
    
    // Test dynamics partial
    WBM::Model<double>::StateMapType A;
    WBM::Model<double>::ContrlMapType B;
    WBM::Model<double>::OutputMapType C;
    WBM::Model<double>::DirectMapType D;
    wb_model.dynamics_partial(A, B, C, D, x, u, t, contact, dt);    

    // Finite difference method for df_dx, dGRF_dx
    double eps = 1e-8;
    WBM::Model<double>::StateType x_eps, x_FD, xnext_post;
    WBM::Model<double>::OutputType GRF_post;
    WBM::Model<double>::StateMapType A_FD;
    WBM::Model<double>::OutputMapType C_FD;
    wb_model.dynamics(xnext, GRF, x, u, t, contact, dt);
    for (size_t i = 0; i < WBM::xs; i++)
    {
        x_eps.setZero();
        x_eps[i] = eps;
        x_FD = x + x_eps;
        wb_model.dynamics(xnext_post, GRF_post, x_FD, u, t, contact, dt);
        A_FD.col(i) = (xnext_post - xnext)/eps;
        C_FD.col(i) = (GRF_post - GRF)/eps;
    }        
    std::cout << "Are A and A_FD close ? " << A.isApprox(A_FD, 1e-5) << "\n";    
    std::cout << "Are C and C_FD close ? " << C.isApprox(C_FD, 1e-5) << "\n";

    // Finite difference for df_du
    WBM::Model<double>::ContrlType u_eps, u_FD;
    WBM::Model<double>::ContrlMapType B_FD;
    for (size_t i = 0; i < WBM::us; i++)
    {
        u_eps.setZero();
        u_eps[i] = eps;
        u_FD = u + u_eps;
        wb_model.dynamics(xnext_post, GRF_post, x, u_FD, t, contact, dt);
        B_FD.col(i) = (xnext_post - xnext)/eps;
    }
    std::cout << "Are B and B_FD close ? " << B.isApprox(B_FD, 1e-5) << "\n";        

    // Test Impact partial
    DMat<double> dP_dx;
    DMat<double> dImp_dx;
    DVec<double> x_postImp, x_preImp, yImp;
    x_preImp = x;    
    wb_model.impact_partial(dP_dx, dImp_dx, x_preImp, contact, nextContact);

    // Finite difference for dP_dx
    DVec<double> x_postImp_post, yImp_post;
    WBM::Model<double>::StateMapType dP_dx_FD;
    WBM::Model<double>::OutputMapType dImp_dx_FD;
    wb_model.impact(x_postImp, yImp, x_preImp, contact, nextContact);
    for (size_t i = 0; i < WBM::xs; i++)
    {
        x_eps.setZero();
        x_eps[i] = eps;
        x_FD = x_preImp + x_eps;
        wb_model.impact(x_postImp_post, yImp_post, x_FD, contact, nextContact);
        dP_dx_FD.col(i) = (x_postImp_post - x_postImp)/eps;
        dImp_dx_FD.col(i) = (yImp_post - yImp)/eps;
    }
    std::cout << "Are dP_dx and dP_dx_FD close ? " << dP_dx.isApprox(dP_dx_FD, sqrt(eps)) << "\n";    
    std::cout << "Are dImp_dx and dImp_dx_FD close ? " << dImp_dx.isApprox(dImp_dx_FD, sqrt(eps)) << "\n";

    /* Test dynamics */
    q.setOnes();
    qd.setOnes();    
    u.setOnes();
    x << q, qd;    
    contact << 1,1,1,1;
    WBM::Model<double>::StateType xdot;
    wb_model.dynamics_continuousTime(xdot, GRF, x, u, contact);
    std::cout << "qd = \n" << std::setprecision(4) << xdot.head(18).transpose() << "\n";
    std::cout << "qdd = \n" << std::setprecision(4) << xdot.tail(18).transpose() << "\n";
    std::cout << "GRF = \n" << std::setprecision(4) << GRF.transpose() << "\n";

    /* Test kinematics */
    VecM<double, 12> EE_pos;
    wb_model.get_footPositions(EE_pos, x);
    std::cout << "EE_pos = \n" << std::setprecision(5) << EE_pos.transpose() << "\n";
}