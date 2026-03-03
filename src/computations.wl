(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  computations.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Implements high-level physics computations and derived
    quantities used in numerical analyses.

  Content:
    - Derived observables
    - Core computational routines
    - Structured outputs

  Design Principle:
    - Uses functions defined in functions.wl
    - Performs physics calculations but no plotting
    - Returns structured data for further processing

===========================================================
*)


(*===============================================*)
(*------------- Compute R00 and Rpm -------------*)
(*===============================================*)
ClearAll[ComputeR00];
ComputeR00[] :=
 Module[{phiratio, R00, dRdB00, dRdBpm, R00Err},

  phiratio = Phi[Mpip0/MBd0, Mpip0/MBd0]/Phi[Mpi00/MBd0, Mpi00/MBd0];
  R00 = 2*phiratio*(BrB0pi0pi00/BrB0pippim0);

  dRdB00 = 2*phiratio/BrB0pippim0;
  dRdBpm = -2*phiratio*BrB0pi0pi00/(BrB0pippim0^2);

  R00Err = Sqrt[(dRdB00*BrB0pi0pi0Err)^2 + (dRdBpm*BrB0pippimErr)^2];

  <|"Value" -> R00, "Error" -> R00Err|>
 ];

ClearAll[ComputeRpm];
ComputeRpm[] :=
 Module[{phiratio, Rpm, dRdBp, dRdB0, RpmErr},

  phiratio = Phi[Mpip0/MBd0, Mpip0/MBd0]/Phi[Mpi00/MBp0, Mpip0/MBp0];
  Rpm = 2*MBd0/MBp0*phiratio*(BrBppippi00/BrB0pippim0)*TBd0/TBp0;

  dRdBp = 2*MBd0/MBp0*phiratio/BrB0pippim0*TBd0/TBp0;
  dRdB0 = -2*MBd0/MBp0*phiratio*BrBppippi00/(BrB0pippim0^2)*TBd0/TBp0;

  RpmErr = Sqrt[(dRdB0*BrB0pippimErr)^2 + (dRdBp*BrBppippi0Err)^2];

  <|"Value" -> Rpm, "Error" -> RpmErr|>
 ];

ClearAll[ComputeR00Rpm];
ComputeR00Rpm[] :=
 <|
  "R00" -> ComputeR00[],
  "Rpm" -> ComputeRpm[]
 |>;
 
(*===============================================*)
(*------------- Compute x and Delta -------------*)
(*===============================================*)
ClearAll[ComputeXDeltaCenter];
Options[ComputeXDeltaCenter] = {"DeltaStepDeg" -> 1};

ComputeXDeltaCenter[inputs_Association, OptionsPattern[]] :=
 Module[
  {
   A, S, gamma, phid, Rpm, R00, Acp00,
   d0, theta0, x0, Delta0,
   deltas, centerRpm, centerR00, centerAcp1, centerAcp2,
   step
  },

  (* required inputs *)
  A = inputs["A"];
  S = inputs["S"];
  gamma = inputs["gamma"];
  phid = inputs["phid"];
  Rpm = inputs["Rpm"];
  R00 = inputs["R00"];

  step = OptionValue["DeltaStepDeg"];
  deltas = Range[0, 360, step] Degree;

  {d0, theta0} = Solvedtheta[A, S, gamma, phid];
  {x0, Delta0} = SolvexDelta[d0, theta0, gamma, Rpm, R00];

  centerRpm = Table[{\[CapitalDelta]/Degree, xPM[\[CapitalDelta], d0, theta0, gamma, Rpm]}, {\[CapitalDelta], deltas}];
  centerR00 = Table[{\[CapitalDelta]/Degree, x00[\[CapitalDelta], d0, theta0, gamma, R00]}, {\[CapitalDelta], deltas}];

  (* optional blue curves if Acp00 exists *)
  centerAcp1 = Missing["NotProvided"];
  centerAcp2 = Missing["NotProvided"];
  If[KeyExistsQ[inputs, "Acp00"],
   Acp00 = inputs["Acp00"];
   centerAcp1 = Table[{\[CapitalDelta]/Degree, xAcpBranches[\[CapitalDelta], d0, theta0, gamma, Acp00][[1]]}, {\[CapitalDelta], deltas}];
   centerAcp2 = Table[{\[CapitalDelta]/Degree, xAcpBranches[\[CapitalDelta], d0, theta0, gamma, Acp00][[2]]}, {\[CapitalDelta], deltas}];
  ];

  <|
   "Center" -> <|"x" -> x0, "Delta" -> Delta0, "d" -> d0, "theta" -> theta0|>,
   "GridDeg" -> (deltas/Degree),
   "Curves" -> <|
     "Rpm" -> centerRpm,
     "R00" -> centerR00,
     "Acp1" -> centerAcp1,
     "Acp2" -> centerAcp2
   |>
  |>
 ];
