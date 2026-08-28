(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  functions_common.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Defines reusable physics and mathematical helper functions.

  Content:
    - Analytical formulae (e.g. amplitude relations)
    - Solvers for hadronic parameters
    - Phase-space functions
    - Auxiliary trigonometric utilities

  Design Principle:
    - Functions should be pure and independent of global state.
    - No numerical scanning or plotting in this file.
    - This file defines the theoretical building blocks.

===========================================================
*)


(*-------- Two-body phase-space factor --------*)
Phi[x_,y_]:=Sqrt[(1-(x+y)^2) (1-(x-y)^2)];

(*-------- Solve d & theta --------*)
Solvedtheta[AcpPM_, ScpPM_, gamma_, phid_] :=
      Module[
         {d, theta, eq1, eq2, sol},
         eq1 = AcpPM == (2 d Sin[theta] Sin[gamma])/(1 - 2 d Cos[theta] Cos[gamma] + d^2); 
         eq2 = ScpPM == -(d^2 Sin[phid] - 2 d Cos[theta] Sin[phid + gamma] + Sin[phid + 2 gamma])/(1 - 2 d Cos[theta] Cos[gamma] + d^2);
         sol = FindRoot[{eq1, eq2}, {{d, 0.6}, {theta, 150 Degree}}];
         {d /. sol, theta /. sol}
         ];

(*-------- Solve x & Delta --------*)
rpi[d_, theta_, gamma_] := 1 - 2 d Cos[theta] Cos[gamma] + d^2;
xPM[Delta_, d_, theta_, gamma_, Rpm_] := -Cos[Delta] + Sqrt[rpi[d, theta, gamma] Rpm - Sin[Delta]^2];
x00[Delta_, d_, theta_, gamma_, R00_] := -d Cos[gamma] Cos[Delta - theta] + Sqrt[rpi[d, theta, gamma] R00 - (1 - Cos[gamma]^2 Cos[Delta - theta]^2) d^2];

SolvexDelta[d_, theta_, gamma_, Rpm_, R00_] :=
      Module[
         {x, Delta, eq, xsol, Deltasol},
         eq = 0 == xPM[Delta, d, theta, gamma, Rpm] - x00[Delta, d, theta, gamma, R00];
         Deltasol = Delta /. FindRoot[eq, {Delta, -50 Degree}];
         xsol = xPM[Deltasol, d, theta, gamma, Rpm];
         {xsol, Deltasol}
         ];

(*-------- Determine x & Delta using Acp00 --------*)
ClearAll[xAcpBranches];
xAcpBranches[Delta_, d_, theta_, gamma_, Acp00_] :=
 Module[{T, x1, x2},
  T = Cos[theta - Delta] Cos[gamma] + Sin[theta - Delta] Sin[gamma]/Acp00;
  x1 = d (-T + Sqrt[T^2 - 1]);
  x2 = d (-T - Sqrt[T^2 - 1]);
  {x1, x2}
 ];

(*-------- S_CP^(pi0 pi0) observable --------*)
ClearAll[SCPpi0pi0];
SCPpi0pi0[x_, Delta_, d_, theta_, gamma_, phid_] :=
 -(
    d^2 Sin[phid]
    + 2 d x Cos[theta - Delta] Sin[phid + gamma]
    + x^2 Sin[phid + 2 gamma]
   )/(
    d^2
    + 2 d x Cos[theta - Delta] Cos[gamma]
    + x^2
   );
