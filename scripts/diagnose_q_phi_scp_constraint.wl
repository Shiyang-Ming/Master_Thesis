(* ::Package:: *)

(*============== Diagnostic script for S_CP(pi0 K_S) q-\[Phi] root structure ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "uncertainty_MC.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "q_phi_scp_constraint.wl"}]];

  (* -------- diagnostic settings -------- *)
  phiPlotRangeDeg = {-75, 78};
  qPlotRange = {0, 3.05};
  phiScanDegList = {5, 10, 15, 20, 30, 40, 60, 75};
  phiScanRangeDeg = {-85, 85};
  phiStepDeg = 1.0;
  qMinCut = 0.05;
  residualTolerance = 10^-7;
  denominatorEps = 10^-8;
  fastMode = False;
  minAcceptedPerPhiBin = 5;
  nSamples = If[fastMode, 150, 1500];
  dqMax = Infinity;
  useDQMaxCut = NumericQ[dqMax] && dqMax < Infinity;

  If[fastMode,
    phiStepDeg = 1.0;
    minAcceptedPerPhiBin = 5;
  ];

  (* -------- local helpers -------- *)
  ClearAll[PhiKey];
  PhiKey[phiDeg_] := Round[1000 phiDeg];

  ClearAll[BuildPurpleSampleState];
  BuildPurpleSampleState[s_Association] :=
   Module[{rObsLocal},
    rObsLocal = <|
      "Value" -> ComputeRPiK[
        s["BrB0pimKp"],
        s["BrBppipK0"],
        s["TBp"],
        s["TBd"]
      ]
    |>;

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
        "R" -> rObsLocal["Value"],
        "AcpPimKp" -> s["AcpPimKp"],
        "ScpPi0Ks" -> s["ScpPi0Ks"],
        "AcpPi0Ks" -> s["AcpPi0Ks"],
        "phid" -> s["phid"]
      |>
    |>
   ];

  ClearAll[ComputePurpleDensePoints];
  ComputePurpleDensePoints[parsLocal_Association, obsLocal_Association] :=
   Module[{phi00Local, branchRes},
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

    <|
      "Phi00Branches" -> phi00Local,
      "BranchResults" -> branchRes,
      "Points" -> Flatten[branchRes[[All, "Points"]], 1]
    |>
   ];

  ClearAll[ComputePurpleCurveForParameters];
  ComputePurpleCurveForParameters[paramAssoc_Association, phiGridDeg_List] :=
   Module[{sampleState, sampleRes},
    sampleState = BuildPurpleSampleState[paramAssoc];
    sampleRes = ComputePurpleDensePoints[sampleState["Pars"], sampleState["Obs"]];
    If[sampleRes === $Failed, Return[$Failed]];
    sampleRes["Points"]
   ];

  ClearAll[ComputeAcceptedRootsByPhi];
  ComputeAcceptedRootsByPhi[parsLocal_Association, obsLocal_Association, phiGridDeg_List] :=
   Module[{phi00Local, qByPhi, totalAcceptedRoots = 0, phiDeg, phi00, solveRes, qVals, key},
    phi00Local = ComputePhi00Branches[
      obsLocal["ScpPi0Ks"],
      obsLocal["AcpPi0Ks"],
      obsLocal["phid"]
    ];
    If[phi00Local === {}, Return[$Failed]];

    qByPhi = AssociationMap[{} &, PhiKey /@ phiGridDeg];
    Do[
      Do[
        solveRes = SolvePurpleConstraintAtPhi[
          phi00,
          phiDeg Degree,
          parsLocal,
          obsLocal,
          "QRange" -> qPlotRange,
          "QMinCut" -> qMinCut,
          "DenominatorEps" -> denominatorEps,
          "ResidualTolerance" -> residualTolerance
        ];
        qVals = solveRes["ValidRoots"];
        totalAcceptedRoots += Length[qVals];
        key = PhiKey[phiDeg];
        If[qVals =!= {} && KeyExistsQ[qByPhi, key],
          qByPhi[key] = Join[qByPhi[key], qVals]
        ];
        ,
        {phiDeg, phiGridDeg}
      ],
      {phi00, phi00Local}
    ];

    <|
      "Phi00Branches" -> phi00Local,
      "QByPhi" -> qByPhi,
      "AcceptedRootCount" -> totalAcceptedRoots
    |>
   ];

  ClearAll[BuildCentralTargetBranch];
  BuildCentralTargetBranch[qByPhi_Association, phiGridDeg_List] :=
   Module[
    {
      phiCandidates, positiveCandidates, seedIndex, seedPhiDeg, seedQ, branchJumpMax,
      selectedAssoc = <||>, phiDeg, qVals, nearestQ, orderedKeys, backwardPhis, forwardPhis,
      prevQ
    },
    phiCandidates = Table[
      phiDeg = phiGridDeg[[i]];
      <|
        "PhiDeg" -> phiDeg,
        "QValues" -> Sort @ Select[qByPhi[PhiKey[phiDeg]], NumericQ]
      |>,
      {i, Length[phiGridDeg]}
    ];
    positiveCandidates = Select[
      phiCandidates,
      #["PhiDeg"] > 0 && #["QValues"] =!= {} &
    ];
    If[positiveCandidates === {}, Return[{}]];

    seedIndex = FirstPosition[
      positiveCandidates[[All, "QValues"]],
      qs_ /; Max[qs] >= 0.8 qPlotRange[[2]],
      Missing["NotFound"]
    ];
    If[seedIndex === Missing["NotFound"],
      seedPhiDeg = positiveCandidates[[First @ Ordering[Max /@ positiveCandidates[[All, "QValues"]], -1], "PhiDeg"]];
      seedQ = Max @ positiveCandidates[[First @ Ordering[Max /@ positiveCandidates[[All, "QValues"]], -1], "QValues"]],
      seedPhiDeg = positiveCandidates[[seedIndex[[1]], "PhiDeg"]];
      seedQ = Max @ positiveCandidates[[seedIndex[[1]], "QValues"]]
    ];

    branchJumpMax = 0.75;
    selectedAssoc[PhiKey[seedPhiDeg]] = seedQ;

    forwardPhis = Select[phiGridDeg, # > seedPhiDeg &];
    prevQ = seedQ;
    Do[
      qVals = Sort @ Select[qByPhi[PhiKey[phiDeg]], NumericQ];
      If[qVals === {}, Continue[]];
      nearestQ = First @ SortBy[qVals, Abs[# - prevQ] &];
      If[Abs[nearestQ - prevQ] <= branchJumpMax,
        selectedAssoc[PhiKey[phiDeg]] = nearestQ;
        prevQ = nearestQ
      ];
      ,
      {phiDeg, forwardPhis}
    ];

    backwardPhis = Reverse @ Select[phiGridDeg, 0 < # < seedPhiDeg &];
    prevQ = seedQ;
    Do[
      qVals = Sort @ Select[qByPhi[PhiKey[phiDeg]], NumericQ];
      If[qVals === {}, Continue[]];
      nearestQ = First @ SortBy[qVals, Abs[# - prevQ] &];
      If[Abs[nearestQ - prevQ] <= branchJumpMax,
        selectedAssoc[PhiKey[phiDeg]] = nearestQ;
        prevQ = nearestQ
      ];
      ,
      {phiDeg, backwardPhis}
    ];

    orderedKeys = Sort @ Keys[selectedAssoc];
    ({#/1000., selectedAssoc[#]} &) /@ orderedKeys
   ];

  (* -------- inputs -------- *)
  pars = <|
    "r" -> r0,
    "delta" -> delta0,
    "rc" -> rc0,
    "deltac" -> deltac0,
    "gamma" -> gamma0,
    "rhoc" -> rhoc0,
    "thetac" -> thetac0
  |>;

  rObs = ComputeRPiKFromInputsWithError[];

  obs = <|
    "R" -> rObs["Value"],
    "AcpPimKp" -> AcppimKp0,
    "ScpPi0Ks" -> Scppi0Ks0,
    "AcpPi0Ks" -> Acppi0Ks0,
    "phid" -> phid0
  |>;

  Print["[INFO] Root-structure diagnostics"];
  Print["[INFO] phiScanDegList = ", phiScanDegList];
  Print["[INFO] phiScanRangeDeg = ", phiScanRangeDeg];
  Print["[INFO] phiStepDeg = ", phiStepDeg];
  Print["[INFO] phiPlotRangeDeg = ", phiPlotRangeDeg];
  Print["[INFO] qPlotRange = ", qPlotRange];
  Print["[INFO] qMinCut = ", qMinCut];
  Print["[INFO] residualTolerance = ", residualTolerance];
  Print["[INFO] fastMode = ", fastMode];
  Print["[INFO] MC samples = ", nSamples];
  Print["[INFO] R = ", N[rObs["Value"], 8], " +/- ", N[rObs["Error"], 4]];

  (* -------- phi00 branches -------- *)
  phi00Branches = ComputePhi00Branches[
    obs["ScpPi0Ks"],
    obs["AcpPi0Ks"],
    obs["phid"]
  ];

  If[phi00Branches === {},
    Print["[ERROR] No valid phi00 branches were generated."];
    Abort[];
  ];

  Print["[INFO] Number of phi00 branches = ", Length[phi00Branches]];

  diagnosticRows = Flatten[
    Table[
      Module[
        {
          phi00, phiDeg, diagRes, classified, rawRoots, nearZeroRoots,
          normalPositiveRootsBeforePlotRange, plotRangeAcceptedRoots,
          denomRejected, residualRejected, plotRangeRejected
        },

        phi00 = phi00Branches[[branchIndex]];
        phiDeg = phiScanDegList[[j]];

        diagRes = DiagnosePurpleRootsAtPhi[
          phi00,
          phiDeg Degree,
          pars,
          obs,
          "QRange" -> qPlotRange,
          "QMinCut" -> qMinCut,
          "DenominatorEps" -> denominatorEps,
          "ResidualTolerance" -> residualTolerance
        ];

        classified = diagRes["Diagnostics"];
        rawRoots = classified[[All, "Root"]];
        nearZeroRoots = Cases[
          classified,
          a_ /; a["Class"] === "near-zero" || a["RejectionReason"] === "below qMinCut" :> a["Root"]
        ];
        normalPositiveRootsBeforePlotRange = Cases[
          classified,
          a_ /; (a["Class"] === "small-positive" || a["Class"] === "normal-positive") &&
                TrueQ[a["ResidualPassed"]] &&
                !TrueQ[a["DenominatorSingular"]] :> a["Root"]
        ];
        plotRangeAcceptedRoots = Cases[classified, a_ /; TrueQ[a["Accepted"]] :> a["Root"]];
        denomRejected = Cases[classified, a_ /; a["RejectionReason"] === "denominator singular" :> a["Root"]];
        residualRejected = Cases[classified, a_ /; a["RejectionReason"] === "residual too large" :> a["Root"]];
        plotRangeRejected = Cases[classified, a_ /; a["RejectionReason"] === "outside plot range" :> a["Root"]];

        Table[
          <|
            "phi00BranchIndex" -> branchIndex,
            "phi00Deg" -> N[phi00/Degree, 12],
            "phiScanDeg" -> phiDeg,
            "rawRootIndex" -> k,
            "rawRoot" -> classified[[k, "Root"]],
            "rawRootClass" -> classified[[k, "Class"]],
            "residual" -> classified[[k, "Residual"]],
            "inPlotRange" -> classified[[k, "InPlotRange"]],
            "denominatorSingular" -> classified[[k, "DenominatorSingular"]],
            "residualPassed" -> classified[[k, "ResidualPassed"]],
            "accepted" -> classified[[k, "Accepted"]],
            "rejectionReason" -> classified[[k, "RejectionReason"]],
            "rawRoots" -> ToString[rawRoots, InputForm],
            "nearZeroRoots" -> ToString[nearZeroRoots, InputForm],
            "normalPositiveRootsBeforePlotRange" -> ToString[normalPositiveRootsBeforePlotRange, InputForm],
            "plotRangeAcceptedRoots" -> ToString[plotRangeAcceptedRoots, InputForm],
            "rootsRejectedByResidual" -> ToString[residualRejected, InputForm],
            "rootsRejectedByDenominator" -> ToString[denomRejected, InputForm],
            "rootsRejectedByPlotRange" -> ToString[plotRangeRejected, InputForm]
          |>,
          {k, Length[classified]}
        ]
      ],
      {branchIndex, Length[phi00Branches]}, {j, Length[phiScanDegList]}
    ],
    2
  ];

  summaryRows = Table[
    Module[
      {
        branchIndex, phiDeg, rows, nearZeroRoots,
        normalPositiveRootsBeforePlotRange, plotRangeAcceptedRoots,
        rootsRejectedByResidual, rootsRejectedByDenominator, rootsRejectedByPlotRange
      },
      branchIndex = bi;
      phiDeg = pd;
      rows = Select[diagnosticRows, #["phi00BranchIndex"] == branchIndex && #["phiScanDeg"] == phiDeg &];
      If[rows === {},
        Return[
          <|
            "phi00BranchIndex" -> branchIndex,
            "phi00Deg" -> N[phi00Branches[[branchIndex]]/Degree, 12],
            "phiScanDeg" -> phiDeg,
            "nRawRoots" -> 0,
            "nAcceptedRoots" -> 0,
            "nearZeroRoots" -> "{}",
            "normalPositiveRootsBeforePlotRange" -> "{}",
            "plotRangeAcceptedRoots" -> "{}",
            "rootsRejectedByResidual" -> "{}",
            "rootsRejectedByDenominator" -> "{}",
            "rootsRejectedByPlotRange" -> "{}"
          |>
        ]
      ];
      nearZeroRoots = DeleteDuplicates @ Cases[rows, row_ /; row["rawRootClass"] === "near-zero" :> row["rawRoot"]];
      normalPositiveRootsBeforePlotRange = DeleteDuplicates @ Cases[
        rows,
        row_ /; (row["rawRootClass"] === "small-positive" || row["rawRootClass"] === "normal-positive") &&
                TrueQ[row["residualPassed"]] &&
                !TrueQ[row["denominatorSingular"]] :> row["rawRoot"]
      ];
      plotRangeAcceptedRoots = DeleteDuplicates @ Cases[rows, row_ /; TrueQ[row["accepted"]] :> row["rawRoot"]];
      rootsRejectedByResidual = DeleteDuplicates @ Cases[rows, row_ /; row["rejectionReason"] === "residual too large" :> row["rawRoot"]];
      rootsRejectedByDenominator = DeleteDuplicates @ Cases[rows, row_ /; row["rejectionReason"] === "denominator singular" :> row["rawRoot"]];
      rootsRejectedByPlotRange = DeleteDuplicates @ Cases[rows, row_ /; row["rejectionReason"] === "outside plot range" :> row["rawRoot"]];
      <|
        "phi00BranchIndex" -> branchIndex,
        "phi00Deg" -> rows[[1, "phi00Deg"]],
        "phiScanDeg" -> phiDeg,
        "nRawRoots" -> Length[rows],
        "nAcceptedRoots" -> Length[plotRangeAcceptedRoots],
        "nearZeroRoots" -> ToString[nearZeroRoots, InputForm],
        "normalPositiveRootsBeforePlotRange" -> ToString[normalPositiveRootsBeforePlotRange, InputForm],
        "plotRangeAcceptedRoots" -> ToString[plotRangeAcceptedRoots, InputForm],
        "rootsRejectedByResidual" -> ToString[rootsRejectedByResidual, InputForm],
        "rootsRejectedByDenominator" -> ToString[rootsRejectedByDenominator, InputForm],
        "rootsRejectedByPlotRange" -> ToString[rootsRejectedByPlotRange, InputForm]
      |>
    ],
    {bi, Length[phi00Branches]}, {pd, phiScanDegList}
  ];
  summaryRows = Flatten[summaryRows, 1];

  Do[
    branchRows = Select[summaryRows, #["phi00BranchIndex"] == i &];
    acceptedAll = DeleteDuplicates @ Flatten[ToExpression /@ branchRows[[All, "plotRangeAcceptedRoots"]]];
    Print["[BRANCH ", i, "] phi00 = ", branchRows[[1, "phi00Deg"]], " deg"];
    Print["[BRANCH ", i, "] accepted-root counts per phi = ", branchRows[[All, "nAcceptedRoots"]]];
    Print["[BRANCH ", i, "] plot-range accepted roots seen = ", acceptedAll];
    ,
    {i, Length[phi00Branches]}
  ];

  debugRows = Select[diagnosticRows, #["phi00BranchIndex"] == 1 && #["phiScanDeg"] == 15 &];
  If[debugRows =!= {},
    Print["[DEBUG] Branch 1, phi = 15 deg"];
    Do[
      Print[
        "[DEBUG] q = ", row["rawRoot"],
        ", q>qMinCut = ", If[NumericQ[row["rawRoot"]], row["rawRoot"] > qMinCut, False],
        ", inPlotRangeQ = ", row["inPlotRange"],
        ", denominatorSafe = ", !TrueQ[row["denominatorSingular"]],
        ", residualPassed = ", row["residualPassed"],
        ", finalAccepted = ", row["accepted"]
      ],
      {row, debugRows}
    ]
  ];

  denseBranchResults = Table[
    ScanPurpleConstraintBranch[
      phi00,
      pars,
      obs,
      "PhiScanRangeDeg" -> phiScanRangeDeg,
      "PhiStepDeg" -> phiStepDeg,
      "QRange" -> qPlotRange,
      "PhiPlotRangeDeg" -> phiPlotRangeDeg,
      "QMinCut" -> qMinCut,
      "DenominatorEps" -> denominatorEps,
      "ResidualTolerance" -> residualTolerance
    ],
    {phi00, phi00Branches}
  ];

  Do[
    denseRes = denseBranchResults[[i]];
    Print["[DENSE ", i, "] phi00 = ", N[denseRes["Phi00"]/Degree, 8], " deg"];
    Print["[DENSE ", i, "] valid points = ", denseRes["PointCount"]];
    Print["[DENSE ", i, "] phi range (deg) = ", denseRes["PhiRangeDeg"]];
    Print["[DENSE ", i, "] q range = ", denseRes["QRange"]];
    ,
    {i, Length[denseBranchResults]}
  ];

  mcPhiGridDeg = Range[phiScanRangeDeg[[1]], phiScanRangeDeg[[2]], phiStepDeg];
  mcPhiGridPlotDegFull = Select[mcPhiGridDeg, phiPlotRangeDeg[[1]] <= # <= phiPlotRangeDeg[[2]] &];
  centralAcceptedRootsByPhiRes = ComputeAcceptedRootsByPhi[pars, obs, mcPhiGridDeg];
  If[centralAcceptedRootsByPhiRes === $Failed,
    Print["[ERROR] Failed to build central accepted roots by phi."];
    Abort[]
  ];
  centralTargetBranch = BuildCentralTargetBranch[
    centralAcceptedRootsByPhiRes["QByPhi"],
    mcPhiGridPlotDegFull
  ];
  If[centralTargetBranch === {},
    Print["[ERROR] Failed to identify central target branch."];
    Abort[]
  ];
  centralTargetQByKey = Association @ Map[PhiKey[#[[1]]] -> #[[2]] &, centralTargetBranch];
  centralTargetPhiGridDeg = centralTargetBranch[[All, 1]];
  acceptedPlotPoints = centralTargetBranch;

  Print["[TARGET] central target branch phi range (deg) = ", MinMax[centralTargetBranch[[All, 1]]]];
  Print["[TARGET] central target branch q range = ", MinMax[centralTargetBranch[[All, 2]]]];
  Print["[TARGET] central target branch point count = ", Length[centralTargetBranch]];

  (* -------- MC purple envelope -------- *)
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

  mcPhiGridPlotDeg = centralTargetPhiGridDeg;
  mcQByPhi = AssociationMap[{} &, PhiKey /@ mcPhiGridPlotDeg];
  mcSampler = Function[{}, SampleParametersMC[mcParamSpecs]];

  Print["[INFO] Building purple MC target-branch envelope..."];

  sampledParamSets = Table[mcSampler[], {nSamples}];
  mcSampleAcceptedCount = 0;
  mcSampleRejectedCount = 0;
  mcRejectedNonTargetRootsCount = 0;
  mcRejectedOutsidePhiSupportCount = 0;
  mcRejectedOutsideQRangeCount = 0;
  mcRejectedNotClosestCount = 0;
  mcRejectedDistanceCutCount = 0;
  mcTargetSampleCurves = {};

  Do[
    Module[{sampleState, sampleRootsRes, samplePts = {}, key, phiDeg, candidates, bestQ, centralQ},
      sampleState = BuildPurpleSampleState[sampledParamSets[[i]]];
      sampleRootsRes = ComputeAcceptedRootsByPhi[sampleState["Pars"], sampleState["Obs"], mcPhiGridDeg];
      If[sampleRootsRes === $Failed,
        mcSampleRejectedCount++;
        Continue[]
      ];

      Do[
        key = PhiKey[phiDeg];
        candidates = Sort @ Select[sampleRootsRes["QByPhi"][key], NumericQ];
        mcRejectedOutsideQRangeCount += Count[candidates, q_ /; !(qPlotRange[[1]] <= q <= qPlotRange[[2]])];
        candidates = Select[candidates, qPlotRange[[1]] <= # <= qPlotRange[[2]] &];
        If[!KeyExistsQ[centralTargetQByKey, key],
          mcRejectedOutsidePhiSupportCount += Length[candidates];
          Continue[]
        ];
        If[candidates === {}, Continue[]];
        centralQ = centralTargetQByKey[key];
        bestQ = First @ SortBy[candidates, Abs[# - centralQ] &];
        mcRejectedNotClosestCount += Max[Length[candidates] - 1, 0];
        mcRejectedNonTargetRootsCount += Max[Length[candidates] - 1, 0];
        If[useDQMaxCut && Abs[bestQ - centralQ] >= dqMax,
          mcRejectedDistanceCutCount++;
          Continue[]
        ];
        mcQByPhi[key] = Join[mcQByPhi[key], {bestQ}];
        AppendTo[samplePts, {phiDeg, bestQ}];
        ,
        {phiDeg, mcPhiGridDeg}
      ];

      If[samplePts === {},
        mcSampleRejectedCount++,
        mcSampleAcceptedCount++;
        AppendTo[mcTargetSampleCurves, samplePts]
      ]
    ],
    {i, Length[sampledParamSets]}
  ];

  mcBandRows = Table[
    Module[{phiDeg, key, vals, lower, upper, qmin, qmax},
      phiDeg = mcPhiGridPlotDeg[[i]];
      key = PhiKey[phiDeg];
      vals = mcQByPhi[key];
      If[vals =!= {} && !VectorQ[vals, NumericQ],
        Print["[ERROR] Non-numeric q-values at phi bin ", phiDeg, ": ", vals];
        Abort[]
      ];
      If[Length[vals] >= minAcceptedPerPhiBin,
        lower = N @ Quantile[vals, 0.16];
        upper = N @ Quantile[vals, 0.84];
        qmin = N @ Min[vals];
        qmax = N @ Max[vals];
        <|
          "PhiDeg" -> phiDeg,
          "Count" -> Length[vals],
          "Q16" -> lower,
          "Q84" -> upper,
          "QMin" -> qmin,
          "QMax" -> qmax
        |>,
        <|
          "PhiDeg" -> phiDeg,
          "Count" -> Length[vals],
          "Q16" -> Missing["InsufficientSamples"],
          "Q84" -> Missing["InsufficientSamples"],
          "QMin" -> If[vals === {}, Missing["NoSamples"], Min[vals]],
          "QMax" -> If[vals === {}, Missing["NoSamples"], Max[vals]]
        |>
      ]
    ],
    {i, Length[mcPhiGridPlotDeg]}
  ];

  insufficientPhiBins = Count[mcBandRows, row_ /; row["Count"] < minAcceptedPerPhiBin];
  envelopeLower = Cases[mcBandRows, row_ /; NumericQ[row["Q16"]] :> {row["PhiDeg"], row["Q16"]}];
  envelopeUpper = Cases[mcBandRows, row_ /; NumericQ[row["Q84"]] :> {row["PhiDeg"], row["Q84"]}];
  envelopeMin = Cases[mcBandRows, row_ /; NumericQ[row["QMin"]] :> {row["PhiDeg"], row["QMin"]}];
  envelopeMax = Cases[mcBandRows, row_ /; NumericQ[row["QMax"]] :> {row["PhiDeg"], row["QMax"]}];
  envelopeMedian = Table[
    Module[{vals},
      vals = mcQByPhi[PhiKey[mcPhiGridPlotDeg[[i]]]];
      If[vals =!= {} && !VectorQ[vals, NumericQ],
        Print["[ERROR] Non-numeric q-values in median bin ", mcPhiGridPlotDeg[[i]], ": ", vals];
        Abort[]
      ];
      If[Length[vals] >= minAcceptedPerPhiBin,
        {mcPhiGridPlotDeg[[i]], N @ Quantile[vals, 0.50]},
        Nothing
      ]
    ],
    {i, Length[mcPhiGridPlotDeg]}
  ];
  envelopeMedian = DeleteCases[envelopeMedian, Nothing];

  envelopeRows = Cases[
    mcBandRows,
    row_ /; NumericQ[row["Q16"]] && NumericQ[row["Q84"]] :>
      <|
        "PhiDeg" -> N @ row["PhiDeg"],
        "QLow16" -> N @ row["Q16"],
        "QMedian" -> N @ Lookup[Association[Rule @@@ envelopeMedian], row["PhiDeg"], Missing["NoMedian"]],
        "QHigh84" -> N @ row["Q84"],
        "NAcceptedAtPhi" -> row["Count"],
        "QMin" -> N @ row["QMin"],
        "QMax" -> N @ row["QMax"],
        "HighGreaterThanLow" -> TrueQ[row["Q84"] > row["Q16"]]
      |>
  ];
  If[!And @@ (NumericQ /@ envelopeRows[[All, "QLow16"]]) ||
     !And @@ (NumericQ /@ envelopeRows[[All, "QMedian"]]) ||
     !And @@ (NumericQ /@ envelopeRows[[All, "QHigh84"]]),
    Print["[ERROR] Non-numeric envelope row detected: ", Select[
      envelopeRows,
      !(NumericQ[#["QLow16"]] && NumericQ[#["QMedian"]] && NumericQ[#["QHigh84"]]) &
    ]];
    Abort[]
  ];

  validEnvelopeBinCount = Length[envelopeRows];
  If[validEnvelopeBinCount > 0,
    qLowMinMax = MinMax[envelopeRows[[All, "QLow16"]]];
    qHighMinMax = MinMax[envelopeRows[[All, "QHigh84"]]];
    envelopeWidths = envelopeRows[[All, "QHigh84"]] - envelopeRows[[All, "QLow16"]];
    maxEnvelopeWidth = Max[envelopeWidths];
    meanEnvelopeWidth = Mean[envelopeWidths];
    allBinsOrderedQ = And @@ envelopeRows[[All, "HighGreaterThanLow"]],
    qLowMinMax = Missing["NoEnvelope"];
    qHighMinMax = Missing["NoEnvelope"];
    maxEnvelopeWidth = Missing["NoEnvelope"];
    meanEnvelopeWidth = Missing["NoEnvelope"];
    allBinsOrderedQ = False
  ];

  bandPolygonPts = If[validEnvelopeBinCount >= 2,
    Join[
      Transpose[{envelopeRows[[All, "PhiDeg"]], envelopeRows[[All, "QHigh84"]]}],
      Reverse @ Transpose[{envelopeRows[[All, "PhiDeg"]], envelopeRows[[All, "QLow16"]]}]
    ],
    {}
  ];

  Print["[MC] samples attempted = ", nSamples];
  Print["[MC] samples with matched target branch = ", mcSampleAcceptedCount];
  Print["[MC] samples rejected (no matched target branch) = ", mcSampleRejectedCount];
  Print["[MC] rejected non-target roots = ", mcRejectedNonTargetRootsCount];
  Print["[MC] rejected roots outside target phi support = ", mcRejectedOutsidePhiSupportCount];
  Print["[MC] rejected roots outside q range = ", mcRejectedOutsideQRangeCount];
  Print["[MC] rejected roots not closest to central target = ", mcRejectedNotClosestCount];
  Print["[MC] rejected roots by dqMax cut = ", mcRejectedDistanceCutCount];
  Print["[MC] dqMax cut active = ", useDQMaxCut];
  Print["[MC] dqMax = ", dqMax];
  Print["[MC] phi bins with insufficient accepted points = ", insufficientPhiBins];
  Print["[MC] final target envelope bins = ", validEnvelopeBinCount];
  Print["[MC] first 10 envelope rows = ", Take[envelopeRows, UpTo[10]]];
  Print["[MC] qLow16 min/max = ", qLowMinMax];
  Print["[MC] qHigh84 min/max = ", qHighMinMax];
  Print["[MC] max band width = ", maxEnvelopeWidth];
  Print["[MC] mean band width = ", meanEnvelopeWidth];
  Print["[MC] qHigh84 > qLow16 for all valid bins = ", allBinsOrderedQ];

  (* -------- output directory -------- *)
  outRoot = FileNameJoin[{projectDir, "results", "q_phi_scp"}];
  If[!DirectoryQ[outRoot],
    CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
  ];

  csvPath = FileNameJoin[{outRoot, "root_diagnostics_mc_helper.csv"}];
  txtPath = FileNameJoin[{outRoot, "root_diagnostics_mc_helper.txt"}];

  csvHeader = Keys[First[diagnosticRows]];
  csvData = Prepend[(Values /@ diagnosticRows), csvHeader];
  Export[csvPath, csvData];

  plotPath = FileNameJoin[{outRoot, "purple_band_mc_helper_plot.png"}];
  plotFig = If[acceptedPlotPoints === {},
    Graphics[
      Text[Style["No plot-range accepted roots", 16, Black], {Mean[phiPlotRangeDeg], Mean[qPlotRange]}],
      PlotRange -> {phiPlotRangeDeg, qPlotRange},
      Frame -> True
    ],
    Show[
      If[Length[bandPolygonPts] >= 4,
        Graphics[{
          Directive[RGBColor[0.76, 0.62, 0.88], Opacity[0.35], EdgeForm[None]],
          Polygon[bandPolygonPts]
        }],
        Graphics[{}]
      ],
      Sequence @@ DeleteCases[{
        If[acceptedPlotPoints === {},
          Nothing,
          ListLinePlot[
            {acceptedPlotPoints},
            PlotStyle -> Directive[RGBColor[0.50, 0.12, 0.70], Thick],
            InterpolationOrder -> 1
          ]
        ]
      },
        Nothing
      ],
      Frame -> True,
      Axes -> False,
      FrameStyle -> Directive[Black, 14],
      FrameTicksStyle -> Directive[Black, 12],
      FrameLabel -> {
        Style["\[Phi] [Degree]", 16, Black],
        Style["q", 16, Black]
      },
      PlotRange -> {phiPlotRangeDeg, qPlotRange},
      PlotRangePadding -> Scaled[0.02],
      ImageSize -> 650,
      AspectRatio -> 0.72,
      Background -> White
    ]
  ];
  Export[plotPath, plotFig, ImageResolution -> 300];

  mcMxPath = FileNameJoin[{outRoot, "purple_band_data_mc_helper.mx"}];
  mcTxtPath = FileNameJoin[{outRoot, "purple_band_data.txt"}];
  Export[
    mcMxPath,
    <|
      "Settings" -> <|
        "FastMode" -> fastMode,
        "NSamples" -> nSamples,
        "DQMax" -> dqMax,
        "UseDQMaxCut" -> useDQMaxCut,
        "PhiScanRangeDeg" -> phiScanRangeDeg,
        "PhiStepDeg" -> phiStepDeg,
        "PhiPlotRangeDeg" -> phiPlotRangeDeg,
        "QPlotRange" -> qPlotRange,
        "QMinCut" -> qMinCut,
        "MinAcceptedPerPhiBin" -> minAcceptedPerPhiBin
      |>,
        "Diagnostics" -> <|
        "SamplesAttempted" -> nSamples,
        "SamplesMatchedTargetBranch" -> mcSampleAcceptedCount,
        "SamplesRejectedNoMatchedTargetBranch" -> mcSampleRejectedCount,
        "RejectedNonTargetRoots" -> mcRejectedNonTargetRootsCount,
        "RejectedOutsidePhiSupportRoots" -> mcRejectedOutsidePhiSupportCount,
        "RejectedOutsideQRangeRoots" -> mcRejectedOutsideQRangeCount,
        "RejectedNotClosestRoots" -> mcRejectedNotClosestCount,
        "RejectedDistanceCutRoots" -> mcRejectedDistanceCutCount,
        "InsufficientPhiBins" -> insufficientPhiBins,
        "MaxBandWidth" -> maxEnvelopeWidth,
        "MeanBandWidth" -> meanEnvelopeWidth
      |>,
      "CentralTargetBranch" -> acceptedPlotPoints,
      "BandRows" -> mcBandRows,
      "EnvelopeRows" -> envelopeRows,
      "EnvelopeLower" -> envelopeLower,
      "EnvelopeMedian" -> envelopeMedian,
      "EnvelopeUpper" -> envelopeUpper,
      "EnvelopeMin" -> envelopeMin,
      "EnvelopeMax" -> envelopeMax,
      "BandPolygonPoints" -> bandPolygonPts
    |>
  ];
  mcCsvPath = FileNameJoin[{outRoot, "purple_target_branch_band_points.csv"}];
  Export[
    mcCsvPath,
    Prepend[
      (Values /@ envelopeRows[[All, {"PhiDeg", "QLow16", "QMedian", "QHigh84", "NAcceptedAtPhi"}]]),
      {"PhiDeg", "QLow16", "QMedian", "QHigh84", "NAcceptedAtPhi"}
    ]
  ];
  Export[
    mcTxtPath,
    StringRiffle[
      Flatten @ {
        {
        "Purple target-branch percentile envelope",
        "fastMode = " <> ToString[fastMode],
        "samplesAttempted = " <> ToString[nSamples],
        "samplesMatchedTargetBranch = " <> ToString[mcSampleAcceptedCount],
        "samplesRejectedNoMatchedTargetBranch = " <> ToString[mcSampleRejectedCount],
        "rejectedNonTargetRoots = " <> ToString[mcRejectedNonTargetRootsCount],
        "rejectedOutsidePhiSupportRoots = " <> ToString[mcRejectedOutsidePhiSupportCount],
        "rejectedOutsideQRangeRoots = " <> ToString[mcRejectedOutsideQRangeCount],
        "rejectedNotClosestRoots = " <> ToString[mcRejectedNotClosestCount],
        "rejectedDistanceCutRoots = " <> ToString[mcRejectedDistanceCutCount],
        "dqMaxCutActive = " <> ToString[useDQMaxCut],
        "dqMax = " <> ToString[dqMax, InputForm],
        "phiBinsInsufficient = " <> ToString[insufficientPhiBins],
        "centralTargetPhiRangeDeg = " <> ToString[MinMax[centralTargetBranch[[All, 1]]], InputForm],
        "centralTargetQRange = " <> ToString[MinMax[centralTargetBranch[[All, 2]]], InputForm],
        "maxBandWidth = " <> ToString[maxEnvelopeWidth, InputForm],
        "meanBandWidth = " <> ToString[meanEnvelopeWidth, InputForm],
        "validEnvelopeBins = " <> ToString[validEnvelopeBinCount],
        "qLow16MinMax = " <> ToString[qLowMinMax, InputForm],
        "qHigh84MinMax = " <> ToString[qHighMinMax, InputForm],
        "allBinsHighGreaterThanLow = " <> ToString[allBinsOrderedQ],
        "phiScanRangeDeg = " <> ToString[phiScanRangeDeg, InputForm],
        "phiStepDeg = " <> ToString[phiStepDeg],
        "phiPlotRangeDeg = " <> ToString[phiPlotRangeDeg, InputForm],
        "qPlotRange = " <> ToString[qPlotRange, InputForm],
        "qMinCut = " <> ToString[qMinCut]
        },
        {"", "Per-phi envelope table:"},
        Table[
          "phiDeg = " <> ToString[row["PhiDeg"], InputForm] <>
          ", qLow16 = " <> ToString[row["QLow16"], InputForm] <>
          ", qMedian = " <> ToString[row["QMedian"], InputForm] <>
          ", qHigh84 = " <> ToString[row["QHigh84"], InputForm] <>
          ", nAcceptedAtPhi = " <> ToString[row["NAcceptedAtPhi"], InputForm],
          {row, envelopeRows}
        ]
      },
      "\n"
    ],
    "String"
  ];

  txtLines = Flatten @ {
    {
      "Purple root diagnostics",
      "R = " <> ToString @ N[rObs["Value"], 12] <> " +/- " <> ToString @ N[rObs["Error"], 8],
      "phiScanDegList = " <> ToString[phiScanDegList, InputForm],
      "phiScanRangeDeg = " <> ToString[phiScanRangeDeg, InputForm],
      "phiStepDeg = " <> ToString[phiStepDeg],
      "phiPlotRangeDeg = " <> ToString[phiPlotRangeDeg, InputForm],
      "qPlotRange = " <> ToString[qPlotRange, InputForm],
      "qMinCut = " <> ToString[qMinCut],
      "residualTolerance = " <> ToString[residualTolerance],
      "fastMode = " <> ToString[fastMode],
      "samplesAttempted = " <> ToString[nSamples],
      "samplesMatchedTargetBranch = " <> ToString[mcSampleAcceptedCount],
      "samplesRejectedNoMatchedTargetBranch = " <> ToString[mcSampleRejectedCount],
      "rejectedNonTargetRoots = " <> ToString[mcRejectedNonTargetRootsCount],
      "rejectedOutsidePhiSupportRoots = " <> ToString[mcRejectedOutsidePhiSupportCount],
      "rejectedOutsideQRangeRoots = " <> ToString[mcRejectedOutsideQRangeCount],
      "rejectedNotClosestRoots = " <> ToString[mcRejectedNotClosestCount],
      "rejectedDistanceCutRoots = " <> ToString[mcRejectedDistanceCutCount],
      "dqMaxCutActive = " <> ToString[useDQMaxCut],
      "dqMax = " <> ToString[dqMax, InputForm],
      "phiBinsInsufficient = " <> ToString[insufficientPhiBins],
      "centralTargetPhiRangeDeg = " <> ToString[MinMax[centralTargetBranch[[All, 1]]], InputForm],
      "centralTargetQRange = " <> ToString[MinMax[centralTargetBranch[[All, 2]]], InputForm],
      "maxBandWidth = " <> ToString[maxEnvelopeWidth, InputForm],
      "meanBandWidth = " <> ToString[meanEnvelopeWidth, InputForm]
    },
    Table[
      {
        "",
        "[Branch " <> ToString[row["phi00BranchIndex"]] <> ", phi=" <> ToString[row["phiScanDeg"]] <> " deg]",
        "phi00Deg = " <> ToString[row["phi00Deg"]],
        "nRawRoots = " <> ToString[row["nRawRoots"]],
        "nAcceptedRoots = " <> ToString[row["nAcceptedRoots"]],
        "nearZeroRoots = " <> row["nearZeroRoots"],
        "normalPositiveRootsBeforePlotRange = " <> row["normalPositiveRootsBeforePlotRange"],
        "plotRangeAcceptedRoots = " <> row["plotRangeAcceptedRoots"],
        "rootsRejectedByResidual = " <> row["rootsRejectedByResidual"],
        "rootsRejectedByDenominator = " <> row["rootsRejectedByDenominator"],
        "rootsRejectedByPlotRange = " <> row["rootsRejectedByPlotRange"]
      },
      {row, summaryRows}
    ]
  };
  Export[txtPath, StringRiffle[txtLines, "\n"], "String"];

  Print["[OK] Saved: ", csvPath];
  Print["[OK] Saved: ", txtPath];
  Print["[OK] Saved: ", plotPath];
  Print["[OK] Saved: ", mcMxPath];
  Print["[OK] Saved: ", mcCsvPath];
  Print["[OK] Saved: ", mcTxtPath];

  Exit[0];

,
  Print["[ERROR] Diagnostic script failed."];
  Exit[1];
];
