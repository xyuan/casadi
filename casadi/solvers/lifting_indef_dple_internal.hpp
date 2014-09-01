/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_LIFTING_INDEF_DPLE_INTERNAL_HPP
#define CASADI_LIFTING_INDEF_DPLE_INTERNAL_HPP

#include "../core/function/dple_internal.hpp"
#include "../core/function/dle_solver.hpp"
#include <casadi/solvers/casadi_dplesolver_lifting_export.h>

/** \defgroup plugin_DpleSolver_lifting
 Solving the Discrete Periodic Lyapunov Equations by
 lifting the entire period to a single Discrete Lyapunov Equation

*/
/** \pluginsection{DpleSolver,lifting} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{DpleSolver,lifting}
  
   @copydoc DPLE_doc
   @copydoc plugin_DpleSolver_lifting
   
       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_LIFTING_EXPORT LiftingIndefDpleInternal : public DpleInternal {
  public:
    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    LiftingIndefDpleInternal(const DpleStructure & st);

    /** \brief  Destructor */
    virtual ~LiftingIndefDpleInternal();

    /** \brief  Clone */
    virtual LiftingIndefDpleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual LiftingIndefDpleInternal* create(const DpleStructure & st) const {
        return new LiftingIndefDpleInternal(st);}

    /** \brief  Create a new DPLE Solver */
    static DpleInternal* creator(const DpleStructure & st)
    { return new LiftingIndefDpleInternal(st);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
    */
    virtual Function getDerivative(int nfwd, int nadj);

    /// A documentation string
    static const std::string meta_doc;

  private:
    /// Main implementation as MXFunction
    Function f_;
    
    /// DleSolver solving the condensed Lyapunov form
    DleSolver dlesolver_;
    
    /// State space dimension
    int n_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_LIFTING_INDEF_DPLE_INTERNAL_HPP
