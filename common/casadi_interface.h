#pragma once
#ifndef CASADI_INTERFACE
#define CASADI_INTERFACE

#define int_T long long int

#include <vector>

/*
  @brief: Get the numerical evaluation of a CasadiGen function and the output sparsity pattern
  @params: 
          arg: T ptr to an array of pointers whose element points to an input variable
          res: T ptr to an array of pointers whose element points to an output variable
          max_sz_res: maximum size of output variables
*/
template<typename T>
void casadi_interface(std::vector<T *> ARG, std::vector<T *> RES, int max_sz_res,
                      int f(const T **, T **, int_T *, T *, int),
                      const int_T *f_sparse_out(int_T),
                      int f_work(int_T *, int_T *, int_T *, int_T *));

//  template <typename T>
//   void casadi_interface(
//       std::vector<DVec<typename T::Scalar>> &arg, std::vector<Eigen::MatrixBase<T> *> &res,
//       int f(const typename T::Scalar **, typename T::Scalar **, int_T *, typename T::Scalar *, int),
//       const int_T *f_sparse_out(int_T),
//       int f_work(int_T *, int_T *, int_T *, int_T *))
//   {
//     std::vector<typename T::Scalar *> ARG;
//     for (size_t idx_arg = 0; idx_arg < arg.size(); idx_arg++)
//     {
//       ARG.push_back(arg[idx_arg].data());
//     }

//     std::vector<typename T::Scalar *> RES;
//     int max_res_size = 0;
//     for (size_t idx_res = 0; idx_res < res.size(); idx_res++)
//     {
//       Eigen::MatrixBase<T> *res_ptr = res[idx_res];
//       res_ptr->setZero();

//       RES.push_back(res_ptr->derived().data());

//       if (res_ptr->size() > max_res_size)
//       {
//         max_res_size = res_ptr->size();
//       }
//     }

//     casadi_interface(ARG, RES, max_res_size, f, f_sparse_out, f_work);
//   }

//   template <typename T>
//   void casadi_interface(
//       std::vector<DVec<typename T::Scalar>> &arg, Eigen::MatrixBase<T> &res,
//       int f(const typename T::Scalar **, typename T::Scalar **, int_T *, typename T::Scalar *, int),
//       const int_T *f_sparse_out(int_T),
//       int f_work(int_T *, int_T *, int_T *, int_T *))
//   {
//     std::vector<Eigen::MatrixBase<T> *> res_vec = {&res};
//     casadi_interface(arg, res_vec, f, f_sparse_out, f_work);
//   }
#endif //CASADI_INTERFACE