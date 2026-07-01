(* ::Package:: *)

(*============== Debug script for purple S_CP(pi0 K_S) MC pipeline ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "uncertainty_MC.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "q_phi_scp_constraint.wl"}]];

  (* -------- settings -------- *)
  phiPlotRangeDeg = {-75, 78};
  qPlotRange = {0, 3.05};
  phiScanRangeDeg = {-85, 85};
  phiStepDeg = 1.0;
  qMinCut = 0.05;
  residualTolerance = 10^-7;
  denominatorEps = 10^-8;
  fastMode = True;
  nSamples = If[fastMode, 120, 400];
  minAcceptedPerPhiBin = 5;

  (* -------- local helpers -------- *)
  ClearAll[PhiKey];
  PhiKey[phiDeg_] := Round[1000 phiDeg];

  ClearAll[BuildPurpleSampleState];
  BuildPurpleSampleState[s_Association] :=
   Module[{rObsLocal},
    rObsLocal = ComputeRPiK[
      s["BrB0pimKp"],
      s["BrBppipK0"],
      s["TBp"],
      s["TBd"]
    ];

    <|
      "Pars" -> <|
        "r" -> s["r"],
        "delta" -> s["delta"],
        "rc" -> s["rc"],
        "deltac" -> s["deltac"],
        "gamma" -> s["gamma"],
        "rhoc" -> s["rhoc"],
        "thetac" -> s["thetac"]
      |>,
      "Obs" -> <|
        "R" -> rObsLocal,
        "AcpPimKp" -> s["AcpPimKp"],
        "ScpPi0Ks" -> s["ScpPi0Ks"],
        "AcpPi0Ks" -> s["AcpPi0Ks"],
        "phid" -> s["phid"]
      |>
    |>
   ];

  ClearAll[ComputePurpleDensePointsDetailed];
  ComputePurpleDensePointsDetailed[parsLocal_Association, obsLocal_Association] :=
   Module[{phi00Local, branchRes, points},
    phi00Local = ComputePhi00Branches[
      obsLocal["ScpPi0Ks"],
      obsLocal["AcpPi0Ks"],
      obsLocal["phid"]
    ];
    If[phi00Local === {}, Return[$Failed]];

    branchRes = Table[
      ScanPurpleConstraintBranch[
        phi00,
        parsLocal,
        obsLocal,
        "PhiScanRangeDeg" -> phiScanRangeDeg,
        "PhiStepDeg" -> phiStepDeg,
        "QRange" -> qPlotRange,
        "PhiPlotRangeDeg" -> phiPlotRangeDeg,
        "QMinCut" -> qMinCut,
        "DenominatorEps" -> denominatorEps,
        "ResidualTolerance" -> residualTolerance
      ],
      {phi00, phi00Local}
    ];

    points = Flatten[branchRes[[All, "Points"]], 1];

    <|
      "Phi00Branches" -> phi00Local,
      "BranchResults" -> branchRes,
      "Points" -> points,
      "AcceptedRootCount" -> Total[branchRes[[All, "ValidRootCount"]]]
    |>
   ];

  ClearAll[ComputePurpleCurveForParameters];
  ComputePurpleCurveForParameters[paramAssoc_Association, phiGridDeg_List] :=
   Module[{sampleState, sampleRes},
    sampleState = BuildPurpleSampleState[paramAssoc];
    sampleRes = ComputePurpleDensePointsDetailed[sampleState["Pars"], sampleState["Obs"]];
    If[sampleRes === $Failed, Return[$Failed]];
    sampleRes["Points"]
   ];

  ClearAll[CurveSummary];
  CurveSummary[curve_] :=
   If[curve === $Failed || curve === {},
    <|
      "PointCount" -> 0,
      "First5" -> {},
      "PhiRange" -> Missing["NoPoints"],
      "QRange" -> Missing["NoPoints"]
    |>,
    <|
      "PointCount" -> Length[curve],
      "First5" -> Take[curve, UpTo[5]],
      "PhiRange" -> MinMax[curve[[All, 1]]],
      "QRange" -> MinMax[curve[[All, 2]]]
    |>
   ];

  ClearAll[FormatAssoc];
  FormatAssoc[a_Association] := ToString[Normal[a], InputForm];

  (* -------- central input -------- *)
  centralParamAssoc = <|
    "ScpPi0Ks" -> Scppi0Ks0,
    "AcpPi0Ks" -> Acppi0Ks0,
    "phid" -> phid0,
    "gamma" -> gamma0,
    "r" -> r0,
    "delta" -> delta0,
    "rc" -> rc0,
    "deltac" -> deltac0,
    "AcpPimKp" -> AcppimKp0,
    "BrB0pimKp" -> BrB0pimKp0,
    "BrBppipK0" -> BrBppipK00,
    "TBp" -> TBp0,
    "TBd" -> TBd0,
    "rhoc" -> rhoc0,
    "thetac" -> thetac0
  |>;

  (* -------- MC parameter specs -------- *)
  mcParamSpecs = {
    {"ScpPi0Ks", Scppi0Ks0, Scppi0KsErr, "Any"},
    {"AcpPi0Ks", Acppi0Ks0, Acppi0KsErr, "Any"},
    {"phid", phid0, phidErr, "Any"},
    {"gamma", gamma0, gammaErr, "Any"},
    {"r", r0, rErr, "Positive"},
    {"delta", delta0, deltaErr, "Any"},
    {"rc", rc0, rcErr, "Positive"},
    {"deltac", deltac0, deltacErr, "Any"},
    {"AcpPimKp", AcppimKp0, AcppimKpErr, "Any"},
    {"BrB0pimKp", BrB0pimKp0, BrB0pimKpErr, "Positive"},
    {"BrBppipK0", BrBppipK00, BrBppipK0Err, "Positive"},
    {"TBp", TBp0, TBpErr, "Positive"},
    {"TBd", TBd0, TBdErr, "Positive"},
    {"rhoc", rhoc0, rhocErr, "Positive"},
    {"thetac", thetac0, thetacErr, "Any"}
  };

  sampler = Function[{}, SampleParametersMC[mcParamSpecs]];
  compute = Function[paramAssoc, ComputePurpleCurveForParameters[paramAssoc, mcPhiGridDeg]];

  mcPhiGridDeg = Range[phiScanRangeDeg[[1]], phiScanRangeDeg[[2]], phiStepDeg];
  mcPhiGridPlotDeg = Select[mcPhiGridDeg, phiPlotRangeDeg[[1]] <= # <= phiPlotRangeDeg[[2]] &];

  debugLines = {};
  AppendTo[debugLines, "Purple MC pipeline debug"];
  AppendTo[debugLines, "fastMode = " <> ToString[fastMode]];
  AppendTo[debugLines, "nSamples = " <> ToString[nSamples]];
  AppendTo[debugLines, "phiScanRangeDeg = " <> ToString[phiScanRangeDeg, InputForm]];
  AppendTo[debugLines, "phiStepDeg = " <> ToString[phiStepDeg, InputForm]];
  AppendTo[debugLines, "phiPlotRangeDeg = " <> ToString[phiPlotRangeDeg, InputForm]];
  AppendTo[debugLines, "qPlotRange = " <> ToString[qPlotRange, InputForm]];
  AppendTo[debugLines, "qMinCut = " <> ToString[qMinCut, InputForm]];

  (* -------- 1. MC sample generation -------- *)
  sampledParamSets = Table[sampler[], {nSamples}];
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "1. MC sample generation"];
  AppendTo[debugLines, "samplesRequested = " <> ToString[nSamples]];
  AppendTo[debugLines, "sampleParameterSetsGenerated = " <> ToString[Length[sampledParamSets]]];
  Do[
    AppendTo[debugLines, "sampleAssoc[" <> ToString[i] <> "] = " <> FormatAssoc[sampledParamSets[[i]]]],
    {i, Min[3, Length[sampledParamSets]]}
  ];

  (* -------- 2. Parameter passing into purple solver -------- *)
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "2. Parameter passing into purple solver"];
  Do[
    state = BuildPurpleSampleState[sampledParamSets[[i]]];
    AppendTo[debugLines, "sampleState[" <> ToString[i] <> "] = <|" <>
      "Scppi0Ks -> " <> ToString[state["Obs"]["ScpPi0Ks"], InputForm] <> ", " <>
      "Acppi0Ks -> " <> ToString[state["Obs"]["AcpPi0Ks"], InputForm] <> ", " <>
      "phid -> " <> ToString[state["Obs"]["phid"], InputForm] <> ", " <>
      "gamma -> " <> ToString[state["Pars"]["gamma"], InputForm] <> ", " <>
      "r -> " <> ToString[state["Pars"]["r"], InputForm] <> ", " <>
      "delta -> " <> ToString[state["Pars"]["delta"], InputForm] <> ", " <>
      "rc -> " <> ToString[state["Pars"]["rc"], InputForm] <> ", " <>
      "deltac -> " <> ToString[state["Pars"]["deltac"], InputForm] <> ", " <>
      "AcppimKp -> " <> ToString[state["Obs"]["AcpPimKp"], InputForm] <> ", " <>
      "R -> " <> ToString[state["Obs"]["R"], InputForm] <> ", " <>
      "rhoc -> " <> ToString[state["Pars"]["rhoc"], InputForm] <> ", " <>
      "thetac -> " <> ToString[state["Pars"]["thetac"], InputForm] <> "|>"],
    {i, Min[3, Length[sampledParamSets]]}
  ];

  (* -------- 3. Single-sample purple curve generation -------- *)
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "3. Single-sample purple curve generation"];
  centralState = BuildPurpleSampleState[centralParamAssoc];
  directCentral = ComputePurpleDensePointsDetailed[centralState["Pars"], centralState["Obs"]];
  wrapperCentralCurve = ComputePurpleCurveForParameters[centralParamAssoc, mcPhiGridDeg];
  centralSummary = CurveSummary[wrapperCentralCurve];
  AppendTo[debugLines, "centralCurvePointCount = " <> ToString[centralSummary["PointCount"]]];
  AppendTo[debugLines, "centralCurveFirst5 = " <> ToString[centralSummary["First5"], InputForm]];
  AppendTo[debugLines, "centralCurvePhiRange = " <> ToString[centralSummary["PhiRange"], InputForm]];
  AppendTo[debugLines, "centralCurveQRange = " <> ToString[centralSummary["QRange"], InputForm]];
  AppendTo[debugLines, "centralPhi00BranchesTried = " <> ToString[If[directCentral === $Failed, 0, Length[directCentral["Phi00Branches"]]]]];
  AppendTo[debugLines, "centralAcceptedRoots = " <> ToString[If[directCentral === $Failed, 0, directCentral["AcceptedRootCount"]]]];

  firstSampleCurves = Table[
    Module[{sampleState, sampleDetailed, sampleCurve, sampleSummary},
      sampleState = BuildPurpleSampleState[sampledParamSets[[i]]];
      sampleDetailed = ComputePurpleDensePointsDetailed[sampleState["Pars"], sampleState["Obs"]];
      sampleCurve = ComputePurpleCurveForParameters[sampledParamSets[[i]], mcPhiGridDeg];
      sampleSummary = CurveSummary[sampleCurve];
      <|
        "SampleIndex" -> i,
        "PointCount" -> sampleSummary["PointCount"],
        "First5" -> sampleSummary["First5"],
        "PhiRange" -> sampleSummary["PhiRange"],
        "QRange" -> sampleSummary["QRange"],
        "Phi00BranchesTried" -> If[sampleDetailed === $Failed, 0, Length[sampleDetailed["Phi00Branches"]]],
        "AcceptedRoots" -> If[sampleDetailed === $Failed, 0, sampleDetailed["AcceptedRootCount"]]
      |>
    ],
    {i, Min[3, Length[sampledParamSets]]}
  ];
  Do[
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].pointCount = " <> ToString[row["PointCount"]]];
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].first5 = " <> ToString[row["First5"], InputForm]];
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].phiRange = " <> ToString[row["PhiRange"], InputForm]];
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].qRange = " <> ToString[row["QRange"], InputForm]];
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].phi00BranchesTried = " <> ToString[row["Phi00BranchesTried"]]];
    AppendTo[debugLines, "sampleCurve[" <> ToString[row["SampleIndex"]] <> "].acceptedRoots = " <> ToString[row["AcceptedRoots"]]];
    ,
    {row, firstSampleCurves}
  ];

  (* -------- 4. Compare central solver vs MC wrapper -------- *)
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "4. Compare central solver vs MC wrapper"];
  AppendTo[debugLines, "Both central and MC wrapper call ComputePurpleDensePointsDetailed -> ScanPurpleConstraintBranch with the same q/denominator/residual cuts."];
  AppendTo[debugLines, "centralWrapperMatchesDirectPoints = " <> ToString[wrapperCentralCurve === If[directCentral === $Failed, $Failed, directCentral["Points"]]]];

  (* -------- 5/6. Per-phi accumulation -------- *)
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "5/6. Per-phi bin accumulation"];
  sampleCurvesRaw = compute /@ sampledParamSets;
  validSampleCurves = Select[sampleCurvesRaw, ListQ];
  nonEmptySampleCurves = Select[validSampleCurves, Length[#] > 0 &];
  emptyListSampleCurveCount = Count[validSampleCurves, curve_ /; curve === {}];
  rejectedSampleCurves = Length[sampleCurvesRaw] - Length[validSampleCurves];
  AppendTo[debugLines, "sampleCurvesRequested = " <> ToString[Length[sampledParamSets]]];
  AppendTo[debugLines, "sampleCurvesValidListQ = " <> ToString[Length[validSampleCurves]]];
  AppendTo[debugLines, "sampleCurvesNonEmpty = " <> ToString[Length[nonEmptySampleCurves]]];
  AppendTo[debugLines, "sampleCurvesEmptyList = " <> ToString[emptyListSampleCurveCount]];
  AppendTo[debugLines, "sampleCurvesRejected = " <> ToString[rejectedSampleCurves]];

  qByPhi = AssociationMap[{} &, PhiKey /@ mcPhiGridPlotDeg];
  Do[
    grouped = GroupBy[
      nonEmptySampleCurves[[i]],
      PhiKey[#[[1]]] &,
      (Last /@ #) &
    ];
    Do[
      If[KeyExistsQ[grouped, key],
        qByPhi[key] = Join[qByPhi[key], grouped[key]]
      ],
      {key, Keys[qByPhi]}
    ],
    {i, Length[nonEmptySampleCurves]}
  ];

  phiBinCounts = Table[
    <|
      "PhiDeg" -> mcPhiGridPlotDeg[[i]],
      "Count" -> Length[qByPhi[PhiKey[mcPhiGridPlotDeg[[i]]]]]
    |>,
    {i, Length[mcPhiGridPlotDeg]}
  ];
  nonEmptyPhiBins = Select[phiBinCounts, #["Count"] > 0 &];
  binsAtLeast5 = Count[phiBinCounts, row_ /; row["Count"] >= 5];
  binsAtLeast10 = Count[phiBinCounts, row_ /; row["Count"] >= 10];

  AppendTo[debugLines, "phiBinsTotal = " <> ToString[Length[mcPhiGridPlotDeg]]];
  AppendTo[debugLines, "first20NonEmptyPhiBins = " <> ToString[Take[nonEmptyPhiBins, UpTo[20]], InputForm]];
  AppendTo[debugLines, "phiBinsAtLeast5 = " <> ToString[binsAtLeast5]];
  AppendTo[debugLines, "phiBinsAtLeast10 = " <> ToString[binsAtLeast10]];

  (* -------- 7. Envelope table construction -------- *)
  AppendTo[debugLines, ""];
  AppendTo[debugLines, "7. Envelope table construction"];
  envelopeRows = Table[
    Module[{phiDeg, vals, qLow16, qMedian, qHigh84},
      phiDeg = mcPhiGridPlotDeg[[i]];
      vals = qByPhi[PhiKey[phiDeg]];
      If[vals =!= {} && !VectorQ[vals, NumericQ],
        AppendTo[debugLines, "ERROR bad qValues at phiDeg = " <> ToString[phiDeg, InputForm] <> ": " <> ToString[vals, InputForm]];
        Abort[]
      ];
      If[Length[vals] >= minAcceptedPerPhiBin,
        qLow16 = N @ Quantile[vals, 0.16];
        qMedian = N @ Quantile[vals, 0.50];
        qHigh84 = N @ Quantile[vals, 0.84];
        <|
          "PhiDeg" -> phiDeg,
          "QLow16" -> qLow16,
          "QMedian" -> qMedian,
          "QHigh84" -> qHigh84,
          "NAcceptedAtPhi" -> Length[vals]
        |>,
        Nothing
      ]
    ],
    {i, Length[mcPhiGridPlotDeg]}
  ];
  envelopeRows = DeleteCases[envelopeRows, Nothing];
  If[!And @@ (NumericQ /@ envelopeRows[[All, "QLow16"]]) ||
     !And @@ (NumericQ /@ envelopeRows[[All, "QMedian"]]) ||
     !And @@ (NumericQ /@ envelopeRows[[All, "QHigh84"]]),
    AppendTo[debugLines, "ERROR non-numeric envelope rows = " <> ToString[
      Select[
        envelopeRows,
        !(NumericQ[#["QLow16"]] && NumericQ[#["QMedian"]] && NumericQ[#["QHigh84"]]) &
      ],
      InputForm
    ]];
    Abort[]
  ];

  AppendTo[debugLines, "envelopeRowCount = " <> ToString[Length[envelopeRows]]];
  AppendTo[debugLines, "first10EnvelopeRows = " <> ToString[Take[envelopeRows, UpTo[10]], InputForm]];
  If[envelopeRows =!= {},
    AppendTo[debugLines, "qLow16MinMax = " <> ToString[MinMax[envelopeRows[[All, "QLow16"]]], InputForm]];
    AppendTo[debugLines, "qHigh84MinMax = " <> ToString[MinMax[envelopeRows[[All, "QHigh84"]]], InputForm]];
    AppendTo[debugLines, "maxBandWidth = " <> ToString[Max[envelopeRows[[All, "QHigh84"]] - envelopeRows[[All, "QLow16"]]], InputForm]],
    AppendTo[debugLines, "qLow16MinMax = Missing[\"NoEnvelope\"]"];
    AppendTo[debugLines, "qHigh84MinMax = Missing[\"NoEnvelope\"]"];
    AppendTo[debugLines, "maxBandWidth = Missing[\"NoEnvelope\"]"]
  ];

  (* -------- export -------- *)
  outRoot = FileNameJoin[{projectDir, "results", "q_phi_scp"}];
  If[!DirectoryQ[outRoot],
    CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
  ];

  debugPath = FileNameJoin[{outRoot, "purple_mc_pipeline_debug.txt"}];
  Export[debugPath, StringRiffle[debugLines, "\n"], "String"];
  Print["[OK] Saved: ", debugPath];

  Exit[0];

,
  Print["[ERROR] Purple MC pipeline debug script failed."];
  Exit[1];
];
