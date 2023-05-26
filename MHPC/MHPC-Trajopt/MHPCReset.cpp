 #include "MHPCReset.h"

 template <typename T>
 void MHPCReset<T>::reset_map(DVec<T>& xnext, 
                              DVec<T>& x,
                              const VecM<int, 4>& contact_cur,
                              const VecM<int, 4>& contact_next,
                              const ModelType& mtype_cur,
                              const ModelType& mtype_next)
 {
     const auto& contact_change = contact_next - contact_cur;
     if (contact_change.cwiseEqual(1).any())
     {
         wbm_ptr->impact(xnext, x, contact_cur, contact_next);
     }else{
         xnext = x;
     }
     
     // If there is a model change at transition, apply the low-rank state projection
     // Also, updates the foot step location for contact foot using WB kinematics
     if (mtype_cur == ModelType::WB &&
         mtype_next == ModelType::SRB)
     {
         xnext = StateProjection * xnext;
         footStepPlanner_ptr->updateFootPosAtTransition(x);
     }  
     
 }

 template <typename T>
 void MHPCReset<T>::reset_map_partial(DMat<T>& df_dx, 
                           DVec<T>& x,
                           const VecM<int, 4>& contact_cur,
                           const VecM<int, 4>& contact_next,
                           const ModelType& mtype_cur,
                           const ModelType& mtype_next)
{
    // Check mode change
    const auto& contact_change = contact_next - contact_cur;
     if (contact_change.cwiseEqual(1).any())
     {
         wbm_ptr->impact_partial(df_dx, x, contact_cur, contact_next);
     }else{
         df_dx = I_wb;
     }
     
     // Check model change
     if (mtype_cur == ModelType::WB &&
         mtype_next == ModelType::SRB)
     {
         df_dx = StateProjection * df_dx;
     }  
}          

template class MHPCReset<double>;