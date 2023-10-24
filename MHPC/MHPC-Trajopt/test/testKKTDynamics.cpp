
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
    double eps = 1e-6;
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
    std::cout << "Are A and A_FD close ? " << A.isApprox(A_FD, 1e-4) << "\n";    
    std::cout << "Are C and C_FD close ? " << C.isApprox(C_FD, 1e-4) << "\n";

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
    std::cout << "Are B and B_FD close ? " << B.isApprox(B_FD, 1e-4) << "\n";        

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
    std::cout << "Test full contact dynamics \n";
    contact.setOnes();
    q.setOnes();
    qd.setOnes();    
    u.setZero();
    x << q, qd;        
    WBM::Model<double>::StateType xdot;
    wb_model.dynamics_continuousTime(xdot, GRF, x, u, contact);
    VecM<double, 18> qdd_ref;
    VecM<double, 12> GRF_ref;    
    qdd_ref.head<6>() << -6.3095,   -4.2604,  -14.1384,   22.9058,   17.9408,   39.8478;
    qdd_ref.tail<12>() <<  90.4579,  -65.9947,   66.0292,  138.4558,   22.0186,    7.8347,
                          -434.0716,    5.2086,   99.9593, -386.0737,  -60.6831,  -18.4982;                       
    GRF_ref << 5.6718,    3.3482,    4.9412,    9.5903,    5.4040,    6.7458 ,
               -40.7861,  -21.8598,  -24.2717,  -37.9467,  -20.9369,  -24.4426;    
    std::cout << "(qdd - qdd_ref).norm = " << (xdot.tail(18) - qdd_ref).norm() << "\n";
    std::cout << "(GRF - GRF_ref).norm = " << (GRF - GRF_ref).norm() << "\n";

    std::cout << "Test free fall dynamics \n";
    contact.setZero();
    x << q, qd;  
    wb_model.dynamics_continuousTime(xdot, GRF, x, u, contact);
    qdd_ref.head<6>() << 0.0167,    0.0347,   -9.8007,    3.0514,   -0.7017,    4.1830;
    qdd_ref.tail<12>() << 0.8268,    0.6603,   -5.2010,   -0.2201,    1.3537,   -4.9655,
                          1.0438,   -1.0924,   -0.5153,    0.1417,   -0.2949,   -0.2164;
    std::cout << "(qdd - qdd_ref).norm = " << (xdot.tail(18) - qdd_ref).norm() << "\n";    

    /* Test kinematics */
    VecM<double, 3> EE_pos[4];
    wb_model.get_footPositions(EE_pos, x);    
}