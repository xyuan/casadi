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


      #include "worhp_interface.hpp"
      #include <string>

      const std::string casadi::WorhpInterface::meta_doc=
      "\n"
"WORHP interface\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| print_time      | OT_BOOLEAN      | true            | Print           |\n"
"|                 |                 |                 | information     |\n"
"|                 |                 |                 | about execution |\n"
"|                 |                 |                 | time            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipBarrier    | OT_REAL         | worhp_p_.qp.ipB | IP barrier      |\n"
"|                 |                 | arrier          | parameter.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipComTol     | OT_REAL         | worhp_p_.qp.ipC | IP              |\n"
"|                 |                 | omTol           | complementarity |\n"
"|                 |                 |                 | tolerance.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipFracBound  | OT_REAL         | worhp_p_.qp.ipF | IP fraction-to- |\n"
"|                 |                 | racBound        | the-boundary    |\n"
"|                 |                 |                 | parameter.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipLsMethod   | OT_STRING       | GenericType()   | Select the      |\n"
"|                 |                 |                 | direct linear   |\n"
"|                 |                 |                 | solver used by  |\n"
"|                 |                 |                 | the IP method.  |\n"
"|                 |                 |                 | (LAPACK::0|MA57 |\n"
"|                 |                 |                 | : only          |\n"
"|                 |                 |                 | available if    |\n"
"|                 |                 |                 | provided by the |\n"
"|                 |                 |                 | user:1|SuperLU: |\n"
"|                 |                 |                 | :2|PARDISO:     |\n"
"|                 |                 |                 | only available  |\n"
"|                 |                 |                 | if provided by  |\n"
"|                 |                 |                 | the user,       |\n"
"|                 |                 |                 | subject to      |\n"
"|                 |                 |                 | license availab |\n"
"|                 |                 |                 | ility:3|MUMPS:  |\n"
"|                 |                 |                 | currently Linux |\n"
"|                 |                 |                 | platforms       |\n"
"|                 |                 |                 | only:5|WSMP:    |\n"
"|                 |                 |                 | subject to      |\n"
"|                 |                 |                 | license availab |\n"
"|                 |                 |                 | ility:6|MA86:   |\n"
"|                 |                 |                 | experimental,   |\n"
"|                 |                 |                 | only available  |\n"
"|                 |                 |                 | if provided by  |\n"
"|                 |                 |                 | the user:7|MA97 |\n"
"|                 |                 |                 | :experimental,  |\n"
"|                 |                 |                 | only available  |\n"
"|                 |                 |                 | if provided by  |\n"
"|                 |                 |                 | the user:8)     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipMinAlpha   | OT_REAL         | worhp_p_.qp.ipM | IP line search  |\n"
"|                 |                 | inAlpha         | minimum step    |\n"
"|                 |                 |                 | size.           |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipRelaxDiv   | OT_REAL         | worhp_p_.qp.ipR | The relaxation  |\n"
"|                 |                 | elaxDiv         | term is divided |\n"
"|                 |                 |                 | by this value   |\n"
"|                 |                 |                 | if successful.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipRelaxMax   | OT_REAL         | worhp_p_.qp.ipR | Maximum         |\n"
"|                 |                 | elaxMax         | relaxation      |\n"
"|                 |                 |                 | value.          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipRelaxMin   | OT_REAL         | worhp_p_.qp.ipR | Mimimum         |\n"
"|                 |                 | elaxMin         | relaxation      |\n"
"|                 |                 |                 | value.          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipRelaxMult  | OT_REAL         | worhp_p_.qp.ipR | The relaxation  |\n"
"|                 |                 | elaxMult        | term is         |\n"
"|                 |                 |                 | multiplied by   |\n"
"|                 |                 |                 | this value if   |\n"
"|                 |                 |                 | unsuccessful.   |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipResTol     | OT_REAL         | worhp_p_.qp.ipR | IP residuals    |\n"
"|                 |                 | esTol           | tolerance.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_ipTryRelax   | OT_BOOLEAN      | worhp_p_.qp.ipT | Enable          |\n"
"|                 |                 | ryRelax         | relaxation      |\n"
"|                 |                 |                 | strategy when   |\n"
"|                 |                 |                 | encountering an |\n"
"|                 |                 |                 | error.          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsItMaxIter  | OT_INTEGER      | worhp_p_.qp.lsI | Maximum number  |\n"
"|                 |                 | tMaxIter        | of iterations   |\n"
"|                 |                 |                 | of the          |\n"
"|                 |                 |                 | iterative       |\n"
"|                 |                 |                 | linear solvers. |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsItMethod   | OT_STRING       | GenericType()   | Select the      |\n"
"|                 |                 |                 | iterative       |\n"
"|                 |                 |                 | linear solver.  |\n"
"|                 |                 |                 | (none:Deactivat |\n"
"|                 |                 |                 | e; use a direct |\n"
"|                 |                 |                 | linear solver.: |\n"
"|                 |                 |                 | 0|CGNR::1|CGNE: |\n"
"|                 |                 |                 | :2|CGS::3|BiCGS |\n"
"|                 |                 |                 | TAB::4)         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsItPrecondM | OT_STRING       | GenericType()   | Select          |\n"
"| ethod           |                 |                 | preconditioner  |\n"
"|                 |                 |                 | for the         |\n"
"|                 |                 |                 | iterative       |\n"
"|                 |                 |                 | linear solver.  |\n"
"|                 |                 |                 | (none:No precon |\n"
"|                 |                 |                 | ditioner.:0|sta |\n"
"|                 |                 |                 | tic:Static      |\n"
"|                 |                 |                 | preconditioner  |\n"
"|                 |                 |                 | (KKT-matrix     |\n"
"|                 |                 |                 | with constant   |\n"
"|                 |                 |                 | lower-right blo |\n"
"|                 |                 |                 | ck).:1|full:Ful |\n"
"|                 |                 |                 | l KKT-          |\n"
"|                 |                 |                 | matrix.:2)      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsRefineMaxI | OT_INTEGER      | worhp_p_.qp.lsR | Maximum number  |\n"
"| ter             |                 | efineMaxIter    | of iterative    |\n"
"|                 |                 |                 | refinement      |\n"
"|                 |                 |                 | steps of the    |\n"
"|                 |                 |                 | direct linear   |\n"
"|                 |                 |                 | solvers.        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsScale      | OT_BOOLEAN      | worhp_p_.qp.lsS | Enables scaling |\n"
"|                 |                 | cale            | on linear       |\n"
"|                 |                 |                 | solver level.   |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsTol        | OT_REAL         | worhp_p_.qp.lsT | Tolerance for   |\n"
"|                 |                 | ol              | the linear      |\n"
"|                 |                 |                 | solver.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_lsTrySimple  | OT_BOOLEAN      | worhp_p_.qp.lsT | Some matrices   |\n"
"|                 |                 | rySimple        | can be solved   |\n"
"|                 |                 |                 | without calling |\n"
"|                 |                 |                 | a linear        |\n"
"|                 |                 |                 | equation solver |\n"
"|                 |                 |                 | .Currently only |\n"
"|                 |                 |                 | diagonal        |\n"
"|                 |                 |                 | matrices are    |\n"
"|                 |                 |                 | supported.Non-  |\n"
"|                 |                 |                 | diagonal        |\n"
"|                 |                 |                 | matrices will   |\n"
"|                 |                 |                 | besolved with   |\n"
"|                 |                 |                 | the chosen      |\n"
"|                 |                 |                 | linear equation |\n"
"|                 |                 |                 | solver.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_maxIter      | OT_INTEGER      | worhp_p_.qp.max | Imposes an      |\n"
"|                 |                 | Iter            | upper limit on  |\n"
"|                 |                 |                 | the number of   |\n"
"|                 |                 |                 | minor solver    |\n"
"|                 |                 |                 | iterations,     |\n"
"|                 |                 |                 | i.e. for the    |\n"
"|                 |                 |                 | quadratic       |\n"
"|                 |                 |                 | subproblem      |\n"
"|                 |                 |                 | solver.If the   |\n"
"|                 |                 |                 | limit is        |\n"
"|                 |                 |                 | reached before  |\n"
"|                 |                 |                 | convergence,    |\n"
"|                 |                 |                 | WORHP will      |\n"
"|                 |                 |                 | activate QP     |\n"
"|                 |                 |                 | recovery        |\n"
"|                 |                 |                 | strategies to   |\n"
"|                 |                 |                 | prevent a       |\n"
"|                 |                 |                 | solver          |\n"
"|                 |                 |                 | breakdown.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_method       | OT_STRING       | GenericType()   | Select the      |\n"
"|                 |                 |                 | solution method |\n"
"|                 |                 |                 | used by the QP  |\n"
"|                 |                 |                 | solver. (ip     |\n"
"|                 |                 |                 | :Interior-Point |\n"
"|                 |                 |                 | method.:1|nsn   |\n"
"|                 |                 |                 | :Nonsmooth-     |\n"
"|                 |                 |                 | Newton method.: |\n"
"|                 |                 |                 | 2|automatic:    |\n"
"|                 |                 |                 | Prefer IP and   |\n"
"|                 |                 |                 | fall back to    |\n"
"|                 |                 |                 | NSN on          |\n"
"|                 |                 |                 | error.:12)      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnBeta      | OT_REAL         | worhp_p_.qp.nsn | NSN stepsize    |\n"
"|                 |                 | Beta            | decrease        |\n"
"|                 |                 |                 | factor.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnGradStep  | OT_BOOLEAN      | worhp_p_.qp.nsn | Enable gradient |\n"
"|                 |                 | GradStep        | steps in the    |\n"
"|                 |                 |                 | NSN method.     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnKKT       | OT_REAL         | worhp_p_.qp.nsn | NSN KKT         |\n"
"|                 |                 | KKT             | tolerance.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnLsMethod  | OT_STRING       | GenericType()   | Select the      |\n"
"|                 |                 |                 | direct linear   |\n"
"|                 |                 |                 | solver used by  |\n"
"|                 |                 |                 | the NSN method. |\n"
"|                 |                 |                 | (SuperLU::2|MA4 |\n"
"|                 |                 |                 | 8: only         |\n"
"|                 |                 |                 | available if    |\n"
"|                 |                 |                 | provided by the |\n"
"|                 |                 |                 | user:4)         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnMinAlpha  | OT_REAL         | worhp_p_.qp.nsn | NSN line search |\n"
"|                 |                 | MinAlpha        | minimum step    |\n"
"|                 |                 |                 | size.           |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_nsnSigma     | OT_REAL         | worhp_p_.qp.nsn | NSN line search |\n"
"|                 |                 | Sigma           | slope           |\n"
"|                 |                 |                 | parameter.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_printLevel   | OT_STRING       | GenericType()   | Controls the    |\n"
"|                 |                 |                 | amount of QP    |\n"
"|                 |                 |                 | solver output.  |\n"
"|                 |                 |                 | (none:No output |\n"
"|                 |                 |                 | .:0|warn:Print  |\n"
"|                 |                 |                 | warnings and er |\n"
"|                 |                 |                 | rors.:1|iterati |\n"
"|                 |                 |                 | ons:Print       |\n"
"|                 |                 |                 | iterations.:2)  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_scaleIntern  | OT_BOOLEAN      | worhp_p_.qp.sca | Enable scaling  |\n"
"|                 |                 | leIntern        | on QP level.    |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_strict       | OT_BOOLEAN      | worhp_p_.qp.str | Use strict      |\n"
"|                 |                 | ict             | termination     |\n"
"|                 |                 |                 | criteria in IP  |\n"
"|                 |                 |                 | method.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
">List of available monitors\n"
"\n"
"+-------------+\n"
"|     Id      |\n"
"+=============+\n"
"| eval_f      |\n"
"+-------------+\n"
"| eval_g      |\n"
"+-------------+\n"
"| eval_grad_f |\n"
"+-------------+\n"
"| eval_h      |\n"
"+-------------+\n"
"| eval_jac_g  |\n"
"+-------------+\n"
"\n"
"\n"
">List of available stats\n"
"\n"
"+--------------------+\n"
"|         Id         |\n"
"+====================+\n"
"| iter_count         |\n"
"+--------------------+\n"
"| iteration          |\n"
"+--------------------+\n"
"| iterations         |\n"
"+--------------------+\n"
"| n_eval_f           |\n"
"+--------------------+\n"
"| n_eval_g           |\n"
"+--------------------+\n"
"| n_eval_grad_f      |\n"
"+--------------------+\n"
"| n_eval_h           |\n"
"+--------------------+\n"
"| n_eval_jac_g       |\n"
"+--------------------+\n"
"| return_code        |\n"
"+--------------------+\n"
"| return_status      |\n"
"+--------------------+\n"
"| t_callback_fun     |\n"
"+--------------------+\n"
"| t_callback_prepare |\n"
"+--------------------+\n"
"| t_eval_f           |\n"
"+--------------------+\n"
"| t_eval_g           |\n"
"+--------------------+\n"
"| t_eval_grad_f      |\n"
"+--------------------+\n"
"| t_eval_h           |\n"
"+--------------------+\n"
"| t_eval_jac_g       |\n"
"+--------------------+\n"
"| t_mainloop         |\n"
"+--------------------+\n"
"\n"
"\n"
"\n"
"\n"
;