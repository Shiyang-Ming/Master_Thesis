(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  inputs_exp.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Contains experimental input parameters and numerical
    constants used throughout the project.

  Content:
    - Particle masses and lifetimes
    - Branching fractions (CP-averaged)
    - CP violation observables
    - CKM parameters

  Notes:
    - This file contains only numerical inputs.
    - No physics calculations should be implemented here.
    - Errors are stored with "Err" suffix.

===========================================================
*)


(*-------- Partical Mass (MeV) --------*)
MBd0 = 5279.72; MBdErr=0.08; 
MBp0 = 5279.41; BBpErr=0.07;
Mpip0 = 139.57039; MpipErr=0.00018;
Mpi00 = 134.9768; Mpi0Err = 0.0005;
MK00 = 497.611; MK0Err = 0.013;
MKp0 = 493.677; MKpErr = 0.013;

(*-------- Particle Life Time (s) --------*)
TBp0 = 1.638*10^-12; TBpErr = 0.004*10^-12;
TBd0 = 1.517*10^-12; TBdErr = 0.004*10^-12;

(*-------- CP-averaged Branch Fraction --------*)
BrB0pi0pi00 = 1.46*10^-6; 
BrB0pi0pi0Err = 0.19*10^-6; (*B0 p19*)

BrB0pippim0 = 5.37*10^-6; 
BrB0pippimErr = 0.20*10^-6; (*B0 p19*)

BrBppippi00 = 5.31*10^-6; 
BrBppippi0Err = 0.26*10^-6; (*Bp p18*)

BrBppipK00= 2.39*10^-5;
BrBppipK0Err = 0.06*10^-5; (*Bp p14*)

BrB0pi0K00 = 1.01*10^-5;
BrB0pi0K0Err = 0.04*10^-5; (*B0 p15*)

BrB0pimKp0 = 2.00*10^-5;
BrB0pimKpErr = 0.04*10^-5; (*B0 p15*)

BrBppi0Kp0 = 1.32*10^-5;
BrBppi0KpErr = 0.04*10^-5; (*Bp p15*)

(*-------- CKM Matrix Parameters --------*)
lambdaCKM0 = 0.22501; lambdaCKMErr = 0.00068; (*rpp2024-rev-ckm-matrix*)
gamma0 = 65.6 Degree; gammaErr = 3.0 Degree; (*HFLAV Summer 2025*)
Vud0 = 0.97367; VudErr = 0.00032;  
Vus0 = 0.22431; VusErr = 0.00085; (*rpp2024-rev-ckm-matrix*)

(*-------- CPV Parameters --------*)
Acppippim0 = 0.314; AcppippimErr = 0.030; (*B0 p227 without minus sign as a convention*)
Scppippim0 = -0.670; ScppippimErr = 0.030; (*B0 p228*)
Acppi0pi00 = 0.23; Acppi0pi0Err = 0.18; (*B0 p228 without minus sign as a convention*)
Scppi0pi00 = 0.61; Scppi0pi0Err = 0.78; (*Belle II*)
Acppi0Ks0 = 0; Acppi0KsErr = 0.08; (*B0 p216*)
Scppi0Ks0 = 0.64; Scppi0KsErr = 0.13; (*B0 p217*)
AcppimKp0 = -0.0831; AcppimKpErr = 0.0031; (*B0 p201*)
AcppipKs0 = -0.003;AcppipKsErr = 0.015; (*Bp p194*)
Acppi0Kp0 = 0.027;Acppi0KpErr = 0.012; (*Bp p194*)
Acppi0KsSM0 = -0.106567; Acppi0KsSMErr = 0.0315165;

(*-------- Hadronic Parameters --------*)
rhoc0 = 0.03; rhocErr = 0.01;
thetac0 = 2.6 Degree; thetacErr = 4.6 Degree;

(*-------- Phid --------*)
phid0 = 45.7 Degree; phidErr = 1 Degree; (*HFLAV 2*beta*)

(*-------- SU(3) Breaking Effect -------*)
rSU3Err = 0.02; thetaSU3Err = 20; (*estimation*)
RTC0 = 1.2; RTCErr = 0.2; (*consider the 20% non-factorizable effects, waiting for lattice results(0308297)*)

