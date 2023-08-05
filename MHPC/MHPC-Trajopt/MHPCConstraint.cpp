#include "MHPCConstraint.h"
#include "casadiGen_MHPC.h"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
/*
    Whole-Body GRF Constraint A*y - b >= 0 where y is the GRF of dim 12x1
*/
template<typename T>
MHPCConstraints::WBGRF<T>::WBGRF(const VecM<int, 4> &ctact)
    :PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>("GRF"),
     mu_fric(0.6)
{
    ctact_status = ctact;
    const int& n_constrained_foot = ctact.cwiseEqual(1).cast<int>().sum();
    this->update_constraint_size(5 * n_constrained_foot);
   
    const auto &mu = mu_fric;
    A.setZero(this->size, WBM::ys);
    b.setZero(this->size);

    MatMN<T, 5, 3> A_leg;
    A_leg << 0, 0, 1,
        -1, 0, mu,
        1, 0, mu,
        0, -1, mu,
        0, 1, mu;

    int i(0);
    for (size_t leg = 0; leg < 4; ++leg)
    {
        if (ctact_status[leg] > 0)
        {
            A.block(5 * i, 3 * leg, 5, 3) = A_leg;
            ++i;
        }               
    }
}    

template <typename T>
void MHPCConstraints::WBGRF<T>::compute_violation(const State &x, 
    const Contrl &u, const Output &y, int k)
{    
    (void)(x);
    (void)(u);
    if (this->data.empty())
    {
        printf("Constraint data is empty \n");
        printf("Creating constaint data \n");
        this->create_data();
    }    
    for (size_t i = 0; i < this->data[k].size(); i++)
    {
        this->data[k][i].g = A.row(i) * y - b[i];        
    }
    this->update_max_violation(k);
}

template <typename T>
void MHPCConstraints::WBGRF<T>::compute_partial(const State &x, 
const Contrl &u, const Output &y, int k)
{
    (void)(x);
    (void)(y);
    (void)(u);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gy = A.row(i).transpose();
    }
}


/*
    Whole-Body Torque Limit Constraint C*u - b >= 0 (upper and louwer toruqe bound)
*/
template <typename T>
MHPCConstraints::TorqueLimit<T>::TorqueLimit(): torque_limit(17.0)                                    
{
    this->update_constraint_size(2 * WBM::us);  // Each actuated joint has an upper and a lower limit

    C.setZero();
    b.setZero();
    C.template topRows<WBM::us>() =  -MatMN<T, WBM::us, WBM::us>::Identity(); 
    C.template bottomRows<WBM::us>().setIdentity();
    b.setConstant(-torque_limit);
}

template <typename T>
void MHPCConstraints::TorqueLimit<T>::compute_violation(const State &x, const Contrl &u, const Output &y, int k)
{
    (void) (x);
    (void) (y);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {
        this->data[k][i].g = C.row(i) * u - b[i];        
    }
    this->update_max_violation(k);
}

template <typename T>
void MHPCConstraints::TorqueLimit<T>::compute_partial(const State &x, const Contrl &u, const Output &y, int k)
{
    (void)(x);
    (void)(y);
    (void)(u);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gu = C.row(i).transpose();
    }
}


/*
    Whole-Body Joint Limit Constraint C*qJ - b > = 0 (upper and lower joint limits)
*/
template <typename T>
MHPCConstraints::JointLimit<T>::JointLimit()
{
    this->update_constraint_size(2 * WBM::nu);  // Each joint has an upper and a lower limit

    C.setZero();
    b.setZero();
    C.template topRows<WBM::nu>() = -MatMN<T, WBM::nu, WBM::nu>::Identity();
    C.template bottomRows<WBM::nu>().setIdentity();
    
    T PI = 3.1415926;
    // To Do: check the the uppper and lower bound
    b << PI / 4, -0.1, 1.15 * PI, -0.1,
         PI, PI - 0.2, .1, PI - 0.2;
}

template <typename T>
void MHPCConstraints::JointLimit<T>::compute_violation(const State &x, const Contrl &u, const Output &y, int k)
{
    (void) (u);
    (void) (y);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {
        const auto& qJ= x.template segment<12>(6);
        this->data[k][i].g = C.row(i) * qJ - b[i];        
    }
    this->update_max_violation(k);
}

template <typename T>
void MHPCConstraints::JointLimit<T>::compute_partial(const State &x, const Contrl &u, const Output &y, int k)
{
    (void)(x);
    (void)(y);
    (void)(u);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gx.template segment<12>(6) = C.row(i).transpose();
    }
}


/*
    Whole-Body Touchdown Constraint
*/
template <typename T>
MHPCConstraints::WBTouchDown<T>::WBTouchDown(const CtactStatusType &TDStatus_in, WBM::Model<T>* wbm_ptr_in)
    :TerminalConstraintBase<T, WBM::xs>("TouchDwon")
{
    ground_height = 0;
    TDStatus = TDStatus_in;
    wbm_ptr = wbm_ptr_in;

    int n_contacts = TDStatus.cwiseEqual(1).template cast<int>().sum();

    this->update_constraint_size(n_contacts);

    TD_EE_heights.setZero(n_contacts);
    TD_EE_Jacobians.setZero(3*n_contacts, WBM::nq);
}

template <typename T>
void MHPCConstraints::WBTouchDown<T>::compute_violation(const State &x)
{
    if (this->size != this->data.size())
    {
        printf("constraint size and constraint data size do not match \n");
        return;
    }

    // If constraint data is empty, create data
    if (this->size == 0)
    {
        printf("constraint size is zero \n");
        this->create_data();
    }

    wbm_ptr->get_contactFootHeights(TD_EE_heights, x, TDStatus);

    for (size_t i = 0; i < this->size; i++)
    {
        this->data[i].h = TD_EE_heights[i] - ground_height;
    }
            
    this->update_max_violation();
}

template <typename T>
void MHPCConstraints::WBTouchDown<T>::compute_partial(const State &x)
{
    wbm_ptr->get_contactFootJacobians(TD_EE_Jacobians, x, TDStatus);
    
    for (size_t i = 0; i < this->size; i++)
    {                
        this->data[i].hx.template head<WBM::nq>() = TD_EE_Jacobians.row(3*i+2);
    }
}

template <typename T>
MHPCConstraints::SRBGRF<T>::SRBGRF(const VecM<int, 4> &ctact)
    : PathConstraintBase<T, SRBM::xs, SRBM::us, SRBM::ys>("GRF"),
      mu_fric(0.7)
{
    ctact_status = ctact;
    int n_constrained_foot = ctact.cwiseEqual(1).cast<int>().sum();
    this->update_constraint_size(5 * n_constrained_foot);
   
    auto &mu = mu_fric;
    A.setZero(this->size, 24);
    b.setZero(this->size);
    MatMN<T, 5, 3> A_leg;
    A_leg << 0, 0, 1,
        -1, 0, mu,
        1, 0, mu,
        0, -1, mu,
        0, 1, mu;

    int i(0);
    for (size_t leg = 0; leg < 4; ++leg)
    {
        if (ctact_status[leg] > 0)
        {
            A.block(5 * i, 3 * leg, 5, 3) = A_leg;
            ++i;
        }               
    }
}

template <typename T>
void MHPCConstraints::SRBGRF<T>::compute_violation(const State& x, const Contrl& u, const Output& y, int k)
{
     if (this->data.empty())
    {
        printf("Constraint data is empty \n");
        printf("Creating constaint data \n");
        this->create_data();
    }

    (void)(x);
    (void)(y);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {
        this->data[k][i].g = A.row(i) * u;        
    }
    this->update_max_violation(k);
}

template <typename T>
void MHPCConstraints::SRBGRF<T>::compute_partial(const State& x, const Contrl& u, const Output& y, int k)
{
    (void)(x);
    (void)(y);
    (void)(u);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gu = A.row(i).transpose();
    }
}

template class MHPCConstraints::WBGRF<double>;
template class MHPCConstraints::TorqueLimit<double>;
template class MHPCConstraints::JointLimit<double>;
template class MHPCConstraints::WBTouchDown<double>;
template class MHPCConstraints::SRBGRF<double>;
