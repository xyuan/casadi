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

#include "lifting_indef_dple_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DPLESOLVER_LIFTING_EXPORT
  casadi_register_dplesolver_lifting(DpleInternal::Plugin* plugin) {
    plugin->creator = LiftingIndefDpleInternal::creator;
    plugin->name = "lifting";
    plugin->doc = LiftingIndefDpleInternal::meta_doc.c_str();
    plugin->version = 20;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_LIFTING_EXPORT casadi_load_dplesolver_lifting() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_lifting);
  }

  LiftingIndefDpleInternal::LiftingIndefDpleInternal(
      const DpleStructure & st) : DpleInternal(st) {

    // set default options
    setOption("name", "unnamed_lifting_indef_dple_solver"); // name of the function

    addOption("dle_solver",            OT_STRING, GenericType(),
              "User-defined Dle solver class. Needed for sensitivities.");
    addOption("dle_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the Dle solver.");

  }

  LiftingIndefDpleInternal::~LiftingIndefDpleInternal() {

  }

  void LiftingIndefDpleInternal::init() {

    DpleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
      "const_dim option set to False: Solver only handles the True case.");
    std::cout << "bar" << std::endl;
    std::cout << "test" << range(1,A_.size()) << std::endl;
    
    Sparsity A;
    if (K_==1) {
      A = A_[0];
    } else {
      std::cout << vector_slice(A_,range(1,A_.size())) << std::endl;
      Sparsity AL = blkdiag(vector_slice(A_,range(1,A_.size())));
      
      Sparsity AL2 = horzcat(AL,Sparsity::sparse(AL.size1(),A_[0].size2()));
      Sparsity AT = horzcat(Sparsity::sparse(A_[0].size1(),AL.size2()),A_[0]);

      A = vertcat(AT,AL2);
    }
    n_ = A_[0].size1();

    A.spy();
    
    // Create an dlesolver instance
    std::string dlesolver_name = getOption("dle_solver");
    dlesolver_ = DleSolver(dlesolver_name, dleStruct("a",A,"v",blkdiag(V_)));

    if (hasSetOption("dle_solver_options")) {
      dlesolver_.setOption(getOption("dle_solver_options"));
    }

    // Initialize the NLP solver
    dlesolver_.init();
    
    dlesolver_.output().sparsity().spy();

  }



  void LiftingIndefDpleInternal::evaluate() {
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), dlesolver_.input(i).begin());
    }
    dlesolver_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(dlesolver_.output(i).begin(), dlesolver_.output(i).end(), output(i).begin());
    }
  }

  Function LiftingIndefDpleInternal::getDerivative(int nfwd, int nadj) {
    return dlesolver_.derivative(nfwd, nadj);
  }


  void LiftingIndefDpleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DpleInternal::deepCopyMembers(already_copied);
  }

  LiftingIndefDpleInternal* LiftingIndefDpleInternal::clone() const {
    // Return a deep copy
    LiftingIndefDpleInternal* node = new LiftingIndefDpleInternal(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


