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

#include "homotopy_nlp_solver.hpp"
#include "homotopy_nlp_internal.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"

using namespace std;
namespace casadi{

  HomotopyNLPSolver::HomotopyNLPSolver(){
  }

  HomotopyNLPInternal* HomotopyNLPSolver::operator->(){
    return static_cast<HomotopyNLPInternal*>(Function::operator->());
  }

  const HomotopyNLPInternal* HomotopyNLPSolver::operator->() const{
    return static_cast<const HomotopyNLPInternal*>(Function::operator->());
  }

  bool HomotopyNLPSolver::checkNode() const{
    return dynamic_cast<const HomotopyNLPInternal*>(get())!=0;
  }


} // namespace casadi

