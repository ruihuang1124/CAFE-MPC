#ifndef MHPCREFERENCE_H
#define MHPCREFERENCE_H

#include <deque>
#include <memory>
#include <algorithm>

#include "QuadReference.h"
#include "WBM.h"
#include "SRBM.h"

#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"
#include "SinglePhaseInterface.h"


        
class WBReference : public SinglePhaseReferenceAbstract<WBM::xs, WBM::us, WBM::ys>
{
public:
    WBReference() {}

    void get_reference_at_t(VecM<double, WBM::xs>& xt,  float t) override;

    void get_reference_at_t(VecM<double, WBM::xs>& xt, VecM<double, WBM::us>& ut, VecM<double, WBM::ys>& yt, float t) override;

    void set_quadruped_reference(shared_ptr<QuadReference> quad_ref) {quad_ref_ = quad_ref;}


private:
    shared_ptr<QuadReference> quad_ref_ = nullptr;  

    QuadAugmentedState* quad_state_t_ptr = nullptr;  

};

class SRBReference : public SinglePhaseReferenceAbstract<SRBM::xs, SRBM::us, SRBM::ys>
{
public:
    SRBReference() {}

    void get_reference_at_t(VecM<double, SRBM::xs>& xt,  float t) override;

    void get_reference_at_t(VecM<double, SRBM::xs>& xt, VecM<double, SRBM::us>& ut, VecM<double, SRBM::ys>& yt, float t) override;

    void set_quadruped_reference(shared_ptr<QuadReference> quad_ref) {quad_ref_ = quad_ref;}


private:
    shared_ptr<QuadReference> quad_ref_ = nullptr;  

    QuadAugmentedState* quad_state_t_ptr = nullptr;  

};



#endif //MHPCREFERENCE_H