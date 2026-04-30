(* ::Package:: *)

(*
============================================================
  parameters.wl
  Master Thesis Project
------------------------------------------------------------

  Purpose:
    Contains numerical theory and hadronic input parameters
    used throughout the project.

  Content:
    - B -> pipi hadronic parameters
    - B -> piK hadronic parameters
    - Electroweak penguin parameters

  Notes:
    - This file contains only numerical inputs.
    - No physics calculations should be implemented here.
    - Errors are stored with "Err" suffix.
    - Angles are stored in radians using Mathematica's Degree multiplier.

============================================================
*)

(*---------- B -> pipi parameters ----------*)

x0 = 1.0326101304342625; xErr = 0.061259490916435864;
Delta0 = -60.428605716100854 Degree; DeltaErr = 9.349804127617578 Degree;
d0 = 0.553546153238287; dErr = 0.07167491286103594;
theta0 = 148.57307378289443 Degree; thetaErr = 4.09561037569086 Degree;


(*---------- B -> piK parameters ----------*)

r0 = 0.09634164769843417; rErr = 0.02362056559040557;
delta0 = 31.426926217105564 Degree; deltaErr = 20.415044069251152 Degree;
rc0 = 0.1692294114576668; rcErr = 0.026675616673925817;
deltac0 = 0.6773308835595112 Degree; deltacErr = 20.634575444652484 Degree;


(*---------- Electroweak penguin parameters ----------*)

q0 = 0.68; qErr = 0.2;
phi0 = 0; phiErr = 0;
omega0 = 0; omegaErr = 0;

