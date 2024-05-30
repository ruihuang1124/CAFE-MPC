#include "SinglePhaseInterface.h"
#include "HSDDP_Utils.h"
#include <cassert>

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::update_weighting_matrix(std::vector<T>& weights)
{
    // weights = [q, r, qf]
    assert((weights.size()==xs_+xs_+us_));

    Eigen::Map<DVec<T>> qw(weights.data(), xs_);
    Eigen::Map<DVec<T>> rw(weights.data()+xs_, us_);
    Eigen::Map<DVec<T>> qfw(weights.data()+xs_+us_, xs_);
    
    Q = qw.asDiagonal();
    R = rw.asDiagonal();
    Qf = qfw.asDiagonal();
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::running_cost(RCost& rcost, const State& x, const Contrl& u, 
                                                   const Output& y, T dt, float t)
{    
    (void) (t);
    const auto& dx = x ;
    const auto& du = u ;
    const auto& dy = y ;
    rcost.l = 0.5*dx.transpose() * Q * dx;
    rcost.l += 0.5*du.transpose() * R * du;
    rcost.l += 0.5*dy.transpose() * S * dy;
    rcost.l *= dt;
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::running_cost_par(RCost& rcost, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, float t)
{
    (void) (t);
    const auto& dx = x ;
    const auto& du = u ;
    const auto& dy = y ;    
    rcost.lx = dt * Q * dx;
    rcost.lu = dt * R * du;
    rcost.ly = dt * S * dy;
    rcost.lxx = dt * Q;
    rcost.luu = dt * R;
    rcost.lux.setZero();
    rcost.lyy = dt * S;    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::terminal_cost(TCost&tcost, const State& x, float tend)
{
    (void) (tend);
    const auto& dx = x ;
    tcost.Phi = dx.transpose() * Qf * dx;
    tcost.Phi *= 0.5;
}                                                   

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::terminal_cost_par(TCost&tcost, const State& x, float tend)
{  
    (void) (tend);
    const auto& dx = x ;
    tcost.Phix = Qf * dx;
    tcost.Phixx = Qf;
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticTrackingCost<T, xs_, us_, ys_>::running_cost(RCost& rcost, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, float t)
{        
#ifdef DEBUG_MODE    
    assert((reference != nullptr));    
#endif
    if (nullptr != reference)
    {
        reference->get_reference_at_t(xr_t, ur_t, yr_t, t);
    } 

    const auto& dx = x - xr_t.template cast<T>();
    const auto& du = u - ur_t.template cast<T>();
    const auto& dy = y - yr_t.template cast<T>();
    QuadraticCost<T, xs_, us_, ys_>::running_cost(rcost, dx, du, dy, dt, t);
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticTrackingCost<T, xs_, us_, ys_>::running_cost_par(RCost& rcost, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, float t)
{        
#ifdef DEBUG_MODE    
    assert((reference != nullptr));    
#endif
    if (nullptr != reference)
    {
        reference->get_reference_at_t(xr_t, ur_t, yr_t, t);
    } 

    const auto& dx = x - xr_t.template cast<T>();
    const auto& du = u - ur_t.template cast<T>();
    const auto& dy = y - yr_t.template cast<T>();
    QuadraticCost<T, xs_, us_, ys_>::running_cost_par(rcost, dx, du, dy, dt, t);
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticTrackingCost<T, xs_, us_, ys_>::terminal_cost(TCost&tcost, const State& x, float tend)
{
#ifdef DEBUG_MODE    
    assert((reference != nullptr));    
#endif
    if (nullptr != reference)
    {
        reference->get_reference_at_t(xr_t, ur_t, yr_t, tend);
    }    

    const auto& dx = x - xr_t.template cast<T>();
    QuadraticCost<T, xs_, us_, ys_>::terminal_cost(tcost, dx, tend);
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticTrackingCost<T, xs_, us_, ys_>::terminal_cost_par(TCost&tcost, const State& x, float tend)
{
#ifdef DEBUG_MODE        
    assert((reference != nullptr));    
#endif
    if (nullptr != reference)
    {
        reference->get_reference_at_t(xr_t, ur_t, yr_t, tend);
    } 

    const auto& dx = x - xr_t.template cast<T>();    
    QuadraticCost<T, xs_, us_, ys_>::terminal_cost_par(tcost, dx, tend);
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::running_cost(RCost& rcostdata, const State& x, const Contrl& u, 
                                                   const Output& y, T dt, float t)
{
    rcostdata.Zeros();
    for (auto cost_ptr:cost_ptrs)
    {
        rcostdata_temp.Zeros();
        cost_ptr->running_cost(rcostdata_temp, x, u, y, dt, t);     
        rcostdata.add(rcostdata_temp);
    }    
}                                                

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::running_cost_par(RCost& rcostdata, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, float t)
{
    for (auto cost_ptr:cost_ptrs)
    {
        rcostdata_temp.Zeros();
        cost_ptr->running_cost_par(rcostdata_temp, x, u, y, dt, t);     
        rcostdata.add(rcostdata_temp);
    }    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::terminal_cost(TCost&tcostdata, const State& x, float tend)
{
    tcostdata.Zeros();
    for (auto cost_ptr:cost_ptrs)
    {
        tcostdata_temp.Zeros();
        cost_ptr->terminal_cost(tcostdata_temp, x, tend);
        tcostdata.add(tcostdata_temp);
    }    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::terminal_cost_par(TCost&tcostdata, const State& x, float tend)
{
    for (auto cost_ptr:cost_ptrs)
    {
        tcostdata_temp.Zeros();
        cost_ptr->terminal_cost_par(tcostdata_temp, x, tend);
        tcostdata.add(tcostdata_temp);
    }    
}

// Explicit specilization for QuadraticCost when implemend separately in cpp
template class QuadraticCost<double,24,24,0>;
template class QuadraticTrackingCost<double,24,24,0>;
template class CostContainer<double,24,24,0>;

template class QuadraticCost<double,36,12,12>;
template class QuadraticTrackingCost<double,36,12,12>;
template class CostContainer<double,36,12,12>;

template class QuadraticCost<double,12,12,0>;
template class QuadraticTrackingCost<double,12,12,0>;
template class CostContainer<double,12,12,0>;
