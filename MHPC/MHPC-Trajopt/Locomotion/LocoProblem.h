#ifndef LOCO_PROBLEM_H
#define LOCO_PROBLEM_H

#include "MHPCProblem.h"

template <typename T>
class LocoProblem : public MHPCProblem<T>
{
public:
    using WBPhase_T = typename MHPCProblem<T>::WBPhase_T;
    using WBState = typename MHPCProblem<T>:: WBState;
    using WBContrl = typename MHPCProblem<T>:: WBContrl;
    using WBOutput = typename MHPCProblem<T>:: WBOutput;
    using WBStateMap = typename MHPCProblem<T>:: WBStateMap;
    using WBContrlMap = typename MHPCProblem<T>:: WBContrlMap;
    using WBOutputMap = typename MHPCProblem<T>:: WBOutputMap;
    using WBDirectMap = typename MHPCProblem<T>:: WBDirectMap;    

public:    
    LocoProblem() {}
    
    virtual void initialize_parameters() override;
    
    virtual void create_problem_one_phase(shared_ptr<WBPhase_T>, int phase_idx) override;

};

#endif