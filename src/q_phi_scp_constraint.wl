(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  q_phi_scp_constraint.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Implements reusable helpers for the diagnostic
    S_CP(B_d -> pi0 K_S) constraint in the q-\[Phi] plane.

  Content:
    - R observable computation
    - \[Phi]00 branch generation
    - Local atildeC / atildeS construction
    - Purple-constraint equation assembly and solving
    - Diagnostic branch scanning utilities

  Design Principle:
    - Keeps the purple-constraint logic separate from the
      existing blue/green q-\[Phi] implementation
    - Uses centralized project inputs
    - Avoids self-referential q = f(q) definitions by
      solving symbolically in qSym

===========================================================
*)


(*=========================================================*)
(*--------------------- R observable ----------------------*)
(*=========================================================*)
ClearAll[ComputeRPiK];
ComputeRPiK[brB0pimKp_, brBpPiPlusK0_, tauBp_, tauBd_] :=
  (brB0pimKp/brBpPiPlusK0) (tauBp/tauBd);

ClearAll[ComputeRPiKFromInputs];
ComputeRPiKFromInputs[] :=
  ComputeRPiK[BrB0pimKp0, BrBppipK00, TBp0, TBd0];

ClearAll[ComputeRPiKError];
ComputeRPiKError[
  brB0pimKp_, brB0pimKpErr_,
  brBpPiPlusK0_, brBpPiPlusK0Err_,
  tauBp_, tauBpErr_,
  tauBd_, tauBdErr_
] :=
 Module[{R0, dRdb0, dRdbp, dRdtp, dRdtd},
  R0 = ComputeRPiK[brB0pimKp, brBpPiPlusK0, tauBp, tauBd];
  dRdb0 = R0/brB0pimKp;
  dRdbp = -R0/brBpPiPlusK0;
  dRdtp = R0/tauBp;
  dRdtd = -R0/tauBd;
  Sqrt[
    (dRdb0 brB0pimKpErr)^2 +
    (dRdbp brBpPiPlusK0Err)^2 +
    (dRdtp tauBpErr)^2 +
    (dRdtd tauBdErr)^2
  ]
 ];

ClearAll[ComputeRPiKFromInputsWithError];
ComputeRPiKFromInputsWithError[] :=
 <|
  "Value" -> ComputeRPiKFromInputs[],
  "Error" -> ComputeRPiKError[
    BrB0pimKp0, BrB0pimKpErr,
    BrBppipK00, BrBppipK0Err,
    TBp0, TBpErr,
    TBd0, TBdErr
  ]
 |>;


(*=========================================================*)
(*-------------------- Phi00 branches ---------------------*)
(*=========================================================*)
ClearAll[NormalizeAngleToPi];
NormalizeAngleToPi[ang_] := Mod[ang + Pi, 2 Pi] - Pi;

ClearAll[NormalizeAngleToTanBranch];
NormalizeAngleToTanBranch[ang_] := Mod[ang + Pi/2, Pi] - Pi/2;

ClearAll[UniqueAnglesByTolerance];
UniqueAnglesByTolerance[angles_List, tol_ : 10^-10] :=
 Module[{sorted, unique = {}},
  sorted = Sort[angles];
  Do[
    If[unique === {} || Abs[ang - Last[unique]] > tol,
      AppendTo[unique, ang]
    ],
    {ang, sorted}
  ];
  unique
 ];

ClearAll[ComputePhi00Branches];
ComputePhi00Branches[scp_, acp_, phid_] :=
 Module[{root, alpha, raw},
  root = 1 - acp^2;
  If[root <= 0, Return[{}]];

  alpha = scp/Sqrt[root];
  If[!NumericQ[alpha] || Abs[alpha] > 1, Return[{}]];

  raw = {
    phid - ArcSin[alpha],
    phid - (Pi - ArcSin[alpha])
  };

  UniqueAnglesByTolerance[NormalizeAngleToTanBranch /@ raw]
 ];


(*=========================================================*)
(*--------------- Purple-constraint helpers ---------------*)
(*=========================================================*)
ClearAll[ComputeAtildeCSym];
ComputeAtildeCSym[Robs_, r_, delta_, gamma_, rhoc_, thetac_, rc_, qSym_, phiScan_] :=
  (Robs - 1 + 2 r Cos[delta] Cos[gamma] + 2 rhoc Cos[thetac] Cos[gamma])/
  (2 rc qSym Cos[phiScan]);

ClearAll[ComputeAtildeSSym];
ComputeAtildeSSym[AcpPimKp_, r_, delta_, gamma_, rc_, qSym_, phiScan_] :=
  (AcpPimKp + 2 r Sin[delta] Sin[gamma])/
  (((4/3) rc qSym Sin[phiScan]));

ClearAll[BuildPurpleConstraintPolynomial];
BuildPurpleConstraintPolynomial[qSym_, phiScan_, phi00_, pars_Association, obs_Association] :=
 Module[
  {
   r, delta, rc, deltac, gamma, rhoc, thetac, Robs, AcpPimKp,
   atildeC, atildeS, chatPlus, quadAc, quadBc, quadDc, expr
  },

  r = pars["r"];
  delta = pars["delta"];
  rc = pars["rc"];
  deltac = pars["deltac"];
  gamma = pars["gamma"];
  rhoc = pars["rhoc"];
  thetac = pars["thetac"];

  Robs = obs["R"];
  AcpPimKp = obs["AcpPimKp"];

  atildeC = ComputeAtildeCSym[Robs, r, delta, gamma, rhoc, thetac, rc, qSym, phiScan];
  atildeS = ComputeAtildeSSym[AcpPimKp, r, delta, gamma, rc, qSym, phiScan];

  chatPlus = atildeC qSym Cos[deltac] + atildeS qSym Sin[deltac];

  quadAc = rc^2 (-Tan[phi00] Cos[2 phiScan] - Sin[2 phiScan]);

  quadBc =
    2 rc Cos[deltac] (Tan[phi00] Cos[phiScan] + Sin[phiScan])
    - (4/3) chatPlus quadAc
    - (2 rc^2 - 2 rc r Cos[deltac - delta])
      (-Tan[phi00] Cos[gamma + phiScan] - Sin[gamma + phiScan]);

  quadDc =
    -Tan[phi00]
    - (2 rc Cos[deltac] - 2 r Cos[delta])
      (Tan[phi00] Cos[gamma] + Sin[gamma])
    + (rc^2 + r^2 - 2 rc r Cos[deltac - delta])
      (-Tan[phi00] Cos[2 gamma] - Sin[2 gamma])
    + (4/3) atildeC qSym rc
      (-Tan[phi00] Cos[phiScan] - Sin[phiScan])
    + (4/9) qSym^2 (atildeS^2 + atildeC^2) quadAc
    + (4/3) (-Tan[phi00] Cos[gamma + phiScan] - Sin[gamma + phiScan])
      (rc^2 chatPlus - rc r (atildeC Cos[delta] + atildeS Sin[delta]));

  expr = Expand[qSym^2 quadAc + qSym quadBc + quadDc];
  Expand @ Numerator @ Together[expr]
 ];

ClearAll[BuildPurpleConstraintEquation];
BuildPurpleConstraintEquation[qSym_, phiScan_, phi00_, pars_Association, obs_Association] :=
 Module[
  {
   r, delta, rc, deltac, gamma, rhoc, thetac, Robs, AcpPimKp,
   atildeC, atildeS, chatPlus, quadAc, quadBc, quadDc
  },

  r = pars["r"];
  delta = pars["delta"];
  rc = pars["rc"];
  deltac = pars["deltac"];
  gamma = pars["gamma"];
  rhoc = pars["rhoc"];
  thetac = pars["thetac"];

  Robs = obs["R"];
  AcpPimKp = obs["AcpPimKp"];

  atildeC = ComputeAtildeCSym[Robs, r, delta, gamma, rhoc, thetac, rc, qSym, phiScan];
  atildeS = ComputeAtildeSSym[AcpPimKp, r, delta, gamma, rc, qSym, phiScan];

  chatPlus = atildeC qSym Cos[deltac] + atildeS qSym Sin[deltac];

  quadAc = rc^2 (-Tan[phi00] Cos[2 phiScan] - Sin[2 phiScan]);

  quadBc =
    2 rc Cos[deltac] (Tan[phi00] Cos[phiScan] + Sin[phiScan])
    - (4/3) chatPlus quadAc
    - (2 rc^2 - 2 rc r Cos[deltac - delta])
      (-Tan[phi00] Cos[gamma + phiScan] - Sin[gamma + phiScan]);

  quadDc =
    -Tan[phi00]
    - (2 rc Cos[deltac] - 2 r Cos[delta])
      (Tan[phi00] Cos[gamma] + Sin[gamma])
    + (rc^2 + r^2 - 2 rc r Cos[deltac - delta])
      (-Tan[phi00] Cos[2 gamma] - Sin[2 gamma])
    + (4/3) atildeC qSym rc
      (-Tan[phi00] Cos[phiScan] - Sin[phiScan])
    + (4/9) qSym^2 (atildeS^2 + atildeC^2) quadAc
    + (4/3) (-Tan[phi00] Cos[gamma + phiScan] - Sin[gamma + phiScan])
      (rc^2 chatPlus - rc r (atildeC Cos[delta] + atildeS Sin[delta]));

  Expand[qSym^2 quadAc + qSym quadBc + quadDc]
 ];

ClearAll[SolvePurpleConstraintAtPhi];
Options[SolvePurpleConstraintAtPhi] = {
  "QRange" -> {0, 3.05},
  "QMinCut" -> 0.05,
  "DenominatorEps" -> 10^-8,
  "ImagTolerance" -> 10^-8,
  "ResidualTolerance" -> 10^-7
};

SolvePurpleConstraintAtPhi[
  phi00_, phiScan_, pars_Association, obs_Association,
  OptionsPattern[]
] :=
 Module[
  {
   qSym, qRange, qMinCut, eps, imagTol, residualTol,
   poly, eqn, rawRoots, candidateRoots, nearZeroRejected, validRoots, residual
  },

  qRange = OptionValue["QRange"];
  qMinCut = OptionValue["QMinCut"];
  eps = OptionValue["DenominatorEps"];
  imagTol = OptionValue["ImagTolerance"];
  residualTol = OptionValue["ResidualTolerance"];

  If[Abs[Sin[phiScan]] < eps || Abs[Cos[phiScan]] < eps,
    Return[
      <|
        "RawRoots" -> {},
        "NearZeroRejectedCount" -> 0,
        "ValidRoots" -> {}
      |>
    ]
  ];

  poly = BuildPurpleConstraintPolynomial[qSym, phiScan, phi00, pars, obs];
  eqn = BuildPurpleConstraintEquation[qSym, phiScan, phi00, pars, obs];

  rawRoots = qSym /. Quiet @ NSolve[poly == 0, qSym];
  rawRoots = Select[rawRoots, NumericQ];
  rawRoots = Select[rawRoots, Abs[Im[#]] < imagTol &];
  rawRoots = Re /@ rawRoots;
  rawRoots = UniqueAnglesByTolerance[rawRoots, 10^-8];

  candidateRoots = Select[rawRoots, qRange[[1]] < # < qRange[[2]] &];
  nearZeroRejected = Count[candidateRoots, q_ /; q <= qMinCut];
  candidateRoots = Select[candidateRoots, # > qMinCut &];

  validRoots = Select[
    candidateRoots,
    Function[qVal,
      If[Abs[qVal Sin[phiScan]] < eps || Abs[qVal Cos[phiScan]] < eps,
        False,
        residual = Quiet @ Check[N[eqn /. qSym -> qVal, 30], Indeterminate];
        NumericQ[residual] && Abs[residual] < residualTol
      ]
    ]
  ];

  <|
    "RawRoots" -> rawRoots,
    "NearZeroRejectedCount" -> nearZeroRejected,
    "ValidRoots" -> validRoots
  |>
 ];

ClearAll[ClassifyPurpleRoot];
Options[ClassifyPurpleRoot] = {
  "QRange" -> {0, 3.05},
  "QMinCut" -> 0.05,
  "DenominatorEps" -> 10^-8,
  "ImagTolerance" -> 10^-8,
  "ResidualTolerance" -> 10^-7
};

ClassifyPurpleRoot[
  qRoot_, phiScan_, phi00_, pars_Association, obs_Association,
  OptionsPattern[]
] :=
 Module[
  {
   qMinCut, eps, imagTol, residualTol, qSym, eqn, qNum, residual,
   qRange, rootType, inPlotRangeQ, denomSingularQ, residualPassQ, acceptedQ,
   rejectionReason
  },

  qMinCut = OptionValue["QMinCut"];
  qRange = OptionValue["QRange"];
  eps = OptionValue["DenominatorEps"];
  imagTol = OptionValue["ImagTolerance"];
  residualTol = OptionValue["ResidualTolerance"];

  eqn = BuildPurpleConstraintEquation[qSym, phiScan, phi00, pars, obs];

  If[!NumericQ[qRoot],
    Return[
      <|
        "Root" -> qRoot,
        "Residual" -> Missing["NotNumeric"],
        "Class" -> "non-numeric",
        "DenominatorSingular" -> Missing["NotApplicable"],
        "Accepted" -> False,
        "RejectionReason" -> "non-numeric root"
      |>
    ]
  ];

  If[Abs[Im[qRoot]] >= imagTol,
    Return[
      <|
        "Root" -> qRoot,
        "Residual" -> Missing["ComplexRoot"],
        "Class" -> "complex",
        "DenominatorSingular" -> Missing["NotApplicable"],
        "Accepted" -> False,
        "RejectionReason" -> "complex root"
      |>
    ]
  ];

  qNum = Re[qRoot];
  residual = Quiet @ Check[N[eqn /. qSym -> qNum, 30], Indeterminate];
  inPlotRangeQ = TrueQ[qRange[[1]] <= qNum && qNum <= qRange[[2]]];
  denomSingularQ = TrueQ[Abs[qNum Sin[phiScan]] < eps || Abs[qNum Cos[phiScan]] < eps];
  residualPassQ = TrueQ[NumericQ[residual] && Abs[residual] < residualTol];

  rootType = Which[
    qNum < 0, "negative",
    qNum == 0, "zero",
    0 < qNum < qMinCut, "near-zero",
    qMinCut <= qNum < 0.3, "small-positive",
    qNum >= 0.3, "normal-positive",
    True, "unclassified"
  ];

  acceptedQ = TrueQ[
    qNum > 0 &&
    qNum > qMinCut &&
    inPlotRangeQ &&
    !denomSingularQ &&
    residualPassQ
  ];

  rejectionReason = Which[
    TrueQ[acceptedQ], "accepted",
    qNum == 0, "exactly zero root",
    0 < qNum < qMinCut, "below qMinCut",
    TrueQ[!inPlotRangeQ], "outside plot range",
    TrueQ[denomSingularQ], "denominator singular",
    TrueQ[!residualPassQ] && !NumericQ[residual], "residual not numeric",
    TrueQ[!residualPassQ], "residual too large",
    qNum < 0, "negative root",
    True, "filtered"
  ];

  <|
    "Root" -> qNum,
    "Residual" -> residual,
    "Class" -> rootType,
    "InPlotRange" -> inPlotRangeQ,
    "DenominatorSingular" -> denomSingularQ,
    "ResidualPassed" -> residualPassQ,
    "Accepted" -> acceptedQ,
    "RejectionReason" -> rejectionReason
  |>
 ];

ClearAll[DiagnosePurpleRootsAtPhi];
Options[DiagnosePurpleRootsAtPhi] = Options[ClassifyPurpleRoot];

DiagnosePurpleRootsAtPhi[
  phi00_, phiScan_, pars_Association, obs_Association,
  OptionsPattern[]
] :=
 Module[
  {
   qSym, poly, rawRoots, classified, qMinCut, eps
  },

  qMinCut = OptionValue["QMinCut"];
  eps = OptionValue["DenominatorEps"];
  poly = BuildPurpleConstraintPolynomial[qSym, phiScan, phi00, pars, obs];
  rawRoots = qSym /. Quiet @ NSolve[poly == 0, qSym];

  classified = ClassifyPurpleRoot[
      #, phiScan, phi00, pars, obs,
        "QMinCut" -> qMinCut,
        "QRange" -> OptionValue["QRange"],
        "DenominatorEps" -> eps,
        "ImagTolerance" -> OptionValue["ImagTolerance"],
        "ResidualTolerance" -> OptionValue["ResidualTolerance"]
    ] & /@ rawRoots;

  <|
    "RawRoots" -> rawRoots,
    "Diagnostics" -> classified,
    "AcceptedRoots" -> Cases[classified, a_ /; TrueQ[a["Accepted"]] :> a["Root"]]
  |>
 ];


(*=========================================================*)
(*---------------- Diagnostic branch scan -----------------*)
(*=========================================================*)
ClearAll[ScanPurpleConstraintBranch];
Options[ScanPurpleConstraintBranch] = {
  "PhiScanRangeDeg" -> {-85, 85},
  "PhiStepDeg" -> 2,
  "QRange" -> {0, 3.05},
  "PhiPlotRangeDeg" -> {-75, 78},
  "DenominatorEps" -> 10^-8,
  "QMinCut" -> 0.05,
  "ResidualTolerance" -> 10^-7
};

ScanPurpleConstraintBranch[
  phi00_, pars_Association, obs_Association,
  OptionsPattern[]
] :=
 Module[
  {
   phiScanRangeDeg, phiStepDeg, qRange, phiPlotRangeDeg, phiGridDeg,
   pts, phiDeg, solveRes, qVals, targetPts,
   rawRootCount, nearZeroRejectedCount, validRootCount
  },

  phiScanRangeDeg = OptionValue["PhiScanRangeDeg"];
  phiStepDeg = OptionValue["PhiStepDeg"];
  qRange = OptionValue["QRange"];
  phiPlotRangeDeg = OptionValue["PhiPlotRangeDeg"];

  rawRootCount = 0;
  nearZeroRejectedCount = 0;
  validRootCount = 0;

  phiGridDeg = Range[phiScanRangeDeg[[1]], phiScanRangeDeg[[2]], phiStepDeg];

  pts = Flatten[
    Table[
      solveRes = SolvePurpleConstraintAtPhi[
        phi00, phiDeg Degree, pars, obs,
        "QRange" -> qRange,
        "QMinCut" -> OptionValue["QMinCut"],
        "DenominatorEps" -> OptionValue["DenominatorEps"],
        "ResidualTolerance" -> OptionValue["ResidualTolerance"]
      ];
      rawRootCount += Length[solveRes["RawRoots"]];
      nearZeroRejectedCount += solveRes["NearZeroRejectedCount"];
      qVals = solveRes["ValidRoots"];
      validRootCount += Length[qVals];
      If[qVals === {},
        Nothing,
        ({phiDeg, #} &) /@ qVals
      ],
      {phiDeg, phiGridDeg}
    ],
    1
  ];

  targetPts = Select[
    pts,
    phiPlotRangeDeg[[1]] <= #[[1]] <= phiPlotRangeDeg[[2]] &&
    qRange[[1]] <= #[[2]] <= qRange[[2]] &
  ];

  <|
    "Phi00" -> phi00,
    "Points" -> pts,
    "RawRootCount" -> rawRootCount,
    "NearZeroRejectedCount" -> nearZeroRejectedCount,
    "ValidRootCount" -> validRootCount,
    "PointCount" -> Length[pts],
    "PhiRangeDeg" -> If[pts === {}, Missing["NoPoints"], MinMax[pts[[All, 1]]]],
    "QRange" -> If[pts === {}, Missing["NoPoints"], MinMax[pts[[All, 2]]]],
    "InTargetWindow" -> (targetPts =!= {})
  |>
 ];
