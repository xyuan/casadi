/// \cond INTERNAL

extern "C" {

  extern void sqic(
   const int32_t  *m, // Number of constraints + 1 (for the objective)
   const int32_t * n, // Number of decision variables
   const int32_t * nnzA, // Number of nonzeros in objective-augmented linear constraint matrix  A
   const int32_t  *indA, // colind of Compressed Column Storage A , length: nnzA
   const int32_t  *locA, // row of  Compressed Column Storage A, length n + 1
   const double *valA, // Values of A
   const double* bl, // Lower bounds to decision variables + objective
   const double* bu, // Upper bounds to decision variables + objective
   const int32_t  *hEtype, // ?
   const int32_t  *hs, // ?
   double *x,  // Decision variables + evaluated linear constraints ((initial + optimal), length n+m
   double *pi, // ?
   double *rc, // Multipliers (initial + optimal), length n+m
   const int32_t * nnzH, // Number of nonzeros in full hessian H
   const int32_t * indH, // colind of Compressed Column Storage H , length: nnzH
   const int32_t * locH, // row of  Compressed Column Storage H, length n + 1
   double* valH
   );

  extern void sqicSolve(
   double* Obj // Output: hessian part of the resulting objective
  );

  extern void sqicSolveStabilized(
   double* Obj, // Output: hessian part of the resulting objective
   double *mu,
   int32_t  *lenpi,
   double* piE
  );

  extern void sqicDestroy();
}
/// \endcond
