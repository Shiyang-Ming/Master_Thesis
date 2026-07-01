(* ::Package:: *)

(*============== Script for final combined q-\[Phi] plot ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "uncertainty_MC.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "q_phi_scan.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "q_phi_scp_constraint.wl"}]];

  (* -------- settings -------- *)
  plotRange = {{-80, 80}, {0, 3.05}};
  imageSize = 650;
  aspectRatio = 0.72;

  blueGreenBranches = {
    {2, 1},  (* green *)
    {3, 1}   (* blue *)
  };

  branchColors = QPhiDefaultBranchColors[];
  purpleFillColor = RGBColor[0.76, 0.62, 0.88];
  purpleLineColor = RGBColor[0.50, 0.12, 0.70];

  phiSM = 0.;
  qSM = 0.68;
  qSMErr = If[ValueQ[qErr] && NumericQ[qErr], qErr, 0.15];

  (* -------- q-phi charged-band settings -------- *)
  qphiNGrid = 1350;
  qphiNSamples = 1000;

  (* -------- purple-band settings -------- *)
  phiPlotRangeDeg = {-75, 78};
  qPlotRange = {0, 3.05};
  phiScanRangeDeg = {-85, 85};
  phiStepDeg = 1.0;
  qMinCut = 0.05;
  residualTolerance = 10^-7;
  denominatorEps = 10^-8;
  fastMode = False;
  minAcceptedPerPhiBin = 5;
  purpleNSamples = If[fastMode, 150, 1500];
  dqMax = Infinity;
  useDQMaxCut = NumericQ[dqMax] && dqMax < Infinity;

  (* -------- local helpers for purple workflow -------- *)
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

  ClearAll[MakePurpleBandGraphic];
  MakePurpleBandGraphic[purpleEnvelopeRows_List] :=
   Module[{upper, lower, upperPieces, lowerPieces, nPieces},
    upper = Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QHigh84"]]}];
    lower = Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QLow16"]]}];
    upperPieces = SplitCurveByJump[upper, 6];
    lowerPieces = SplitCurveByJump[lower, 6];
    nPieces = Min[Length[upperPieces], Length[lowerPieces]];
    Graphics[{
      Directive[purpleFillColor, Opacity[0.35], EdgeForm[None]],
      Table[
        If[
          Length[upperPieces[[i]]] >= 2 && Length[lowerPieces[[i]]] >= 2,
          Polygon[Join[upperPieces[[i]], Reverse[lowerPieces[[i]]]]],
          Nothing
        ],
        {i, nPieces}
      ]
    }]
   ];

  ClearAll[MakePurpleBoundaryPlots];
  MakePurpleBoundaryPlots[purpleEnvelopeRows_List] :=
   Module[{upper, lower},
    upper = Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QHigh84"]]}];
    lower = Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QLow16"]]}];
    DeleteCases[
      {
        If[upper === {}, Nothing,
          ListLinePlot[
            SplitCurveByJump[upper, 6],
            PlotStyle -> Directive[purpleLineColor, Thick],
            InterpolationOrder -> 1
          ]
        ],
        If[lower === {}, Nothing,
          ListLinePlot[
            SplitCurveByJump[lower, 6],
            PlotStyle -> Directive[purpleLineColor, Thick],
            InterpolationOrder -> 1
          ]
        ]
      },
      Nothing
    ]
   ];

  ClearAll[PrintPurpleEnvelopeDiagnostics];
  PrintPurpleEnvelopeDiagnostics[purpleCentralCurve_List, purpleEnvelopeRows_List] :=
   Module[{phiVals, qLowVals, qHighVals, widths, gapPos, gapRows},
    Print["[DEBUG] Purple central branch point count = ", Length[purpleCentralCurve]];
    If[purpleCentralCurve =!= {},
      Print["[DEBUG] Purple central branch phi range = ", MinMax[purpleCentralCurve[[All, 1]]]];
      Print["[DEBUG] Purple central branch q range = ", MinMax[purpleCentralCurve[[All, 2]]]]
    ];

    Print["[DEBUG] Purple envelope row count = ", Length[purpleEnvelopeRows]];
    If[purpleEnvelopeRows === {}, Return[]];

    phiVals = purpleEnvelopeRows[[All, "PhiDeg"]];
    qLowVals = purpleEnvelopeRows[[All, "QLow16"]];
    qHighVals = purpleEnvelopeRows[[All, "QHigh84"]];
    widths = (qHighVals - qLowVals)/2;

    Print["[DEBUG] Purple envelope phi range = ", MinMax[phiVals]];
    Print["[DEBUG] Purple envelope width range = ", MinMax[widths]];
    Print["[DEBUG] Purple envelope sorted-by-phi Q = ", OrderedQ[phiVals]];
    Print["[DEBUG] Purple envelope any negative width Q = ", AnyTrue[widths, # < 0 &]];
    Print["[DEBUG] Purple envelope any non-numeric width Q = ", !VectorQ[widths, NumericQ]];

    gapPos = Flatten @ Position[Differences[phiVals], d_ /; d > 6];
    Print["[DEBUG] Purple envelope gap indices (>6 deg) = ", gapPos];
    If[gapPos =!= {},
      gapRows = Flatten[
        Table[
          {
            <|
              "Side" -> "BeforeGap",
              "PhiDeg" -> phiVals[[i]],
              "QLow16" -> qLowVals[[i]],
              "QHigh84" -> qHighVals[[i]],
              "SigmaHalfWidth" -> widths[[i]]
            |>,
            <|
              "Side" -> "AfterGap",
              "PhiDeg" -> phiVals[[i + 1]],
              "QLow16" -> qLowVals[[i + 1]],
              "QHigh84" -> qHighVals[[i + 1]],
              "SigmaHalfWidth" -> widths[[i + 1]]
            |>
          },
          {i, gapPos}
        ],
        1
      ];
      Print["[DEBUG] Purple envelope rows around gaps = ", gapRows]
    ];
   ];

  ClearAll[MakeSMGraphic];
  MakeSMGraphic[] :=
   Graphics[{
     Black,
     AbsoluteThickness[2.0],
     Line[{{phiSM, qSM - qSMErr}, {phiSM, qSM + qSMErr}}],
     AbsolutePointSize[8],
     Point[{phiSM, qSM}],
     Text[Style["SM", 15, Black], {phiSM + 6, qSM + 0.16}]
   }];

  (* -------- blue/green workflow from scripts/q_phi_scan.wl -------- *)
  centralPars = <|
    "gamma" -> gamma0,
    "RTC" -> RTC0,
    "Vus" -> Vus0,
    "Vud" -> Vud0,
    "rhoc" -> rhoc0,
    "thetac" -> thetac0,
    "rc" -> rc0,
    "deltac" -> deltac0,
    "AcpPiPlusK0" -> AcppipKs0,
    "AcpPi0KPlus" -> Acppi0Kp0,
    "BrPiPlusK0" -> BrBppipK00,
    "BrPi0KPlus" -> BrBppi0Kp0,
    "BrPiPlusPi0" -> BrBppippi00
  |>;

  errorPars = <|
    "gamma" -> gammaErr,
    "RTC" -> RTCErr,
    "Vus" -> VusErr,
    "Vud" -> VudErr,
    "rhoc" -> rhocErr,
    "thetac" -> thetacErr,
    "rc" -> rcErr,
    "deltac" -> deltacErr,
    "AcpPiPlusK0" -> AcppipKsErr,
    "AcpPi0KPlus" -> Acppi0KpErr,
    "BrPiPlusK0" -> BrBppipK0Err,
    "BrPi0KPlus" -> BrBppi0KpErr,
    "BrPiPlusPi0" -> BrBppippi0Err
  |>;

  qphiParamSpecs = {
    {"gamma", centralPars["gamma"], errorPars["gamma"], "Any"},
    {"RTC", centralPars["RTC"], errorPars["RTC"], "Positive"},
    {"Vus", centralPars["Vus"], errorPars["Vus"], "Positive"},
    {"Vud", centralPars["Vud"], errorPars["Vud"], "Positive"},
    {"rhoc", centralPars["rhoc"], errorPars["rhoc"], "Positive"},
    {"thetac", centralPars["thetac"], errorPars["thetac"], "Any"},
    {"rc", centralPars["rc"], errorPars["rc"], "Positive"},
    {"deltac", centralPars["deltac"], errorPars["deltac"], "Any"},
    {"AcpPiPlusK0", centralPars["AcpPiPlusK0"], errorPars["AcpPiPlusK0"], "Any"},
    {"AcpPi0KPlus", centralPars["AcpPi0KPlus"], errorPars["AcpPi0KPlus"], "Any"},
    {"BrPiPlusK0", centralPars["BrPiPlusK0"], errorPars["BrPiPlusK0"], "Positive"},
    {"BrPi0KPlus", centralPars["BrPi0KPlus"], errorPars["BrPi0KPlus"], "Positive"},
    {"BrPiPlusPi0", centralPars["BrPiPlusPi0"], errorPars["BrPiPlusPi0"], "Positive"}
  };

  Print["[INFO] Building blue/green q-phi bands with nGrid = ", qphiNGrid];

  centerRes = ComputeQPhiBranchData[centralPars, qphiNGrid, MKp0];
  If[centerRes === $Failed,
    Print["[ERROR] Central q-phi branch construction failed."];
    Abort[];
  ];
  branchData = centerRes["BranchData"];

  qphiCompute = Function[in, ComputeQPhiBranchData[in, qphiNGrid, MKp0]];
  qphiSampler = Function[{}, SampleParametersMC[qphiParamSpecs]];

  Print["[INFO] Running blue/green Monte Carlo with nSamples = ", qphiNSamples];

  qphiSampleResults = GenerateMCSamples[
    qphiCompute,
    qphiSampler,
    qphiNSamples,
    Function[res, res["BranchData"]]
  ];
  qphiValidSamples = Length[qphiSampleResults];
  If[qphiValidSamples == 0,
    Print["[ERROR] No valid blue/green Monte Carlo samples were produced."];
    Abort[];
  ];

  blueGreenBandGraphics = MakeQPhiBandGraphics[
    branchData,
    qphiSampleResults,
    branchColors,
    "KeepBranches" -> blueGreenBranches,
    "Jump" -> 40,
    "Opacity" -> 0.22
  ];

  blueGreenCurvePlots = MakeQPhiBranchPlots[
    branchData,
    branchColors,
    "KeepBranches" -> blueGreenBranches,
    "Jump" -> 40
  ];

  (* -------- purple workflow from scripts/diagnose_q_phi_scp_constraint.wl -------- *)
  purpleMcParamSpecs = {
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

  purplePars = <|
    "r" -> r0,
    "delta" -> delta0,
    "rc" -> rc0,
    "deltac" -> deltac0,
    "gamma" -> gamma0,
    "rhoc" -> rhoc0,
    "thetac" -> thetac0
  |>;

  purpleObs = <|
    "R" -> ComputeRPiKFromInputs[],
    "AcpPimKp" -> AcppimKp0,
    "ScpPi0Ks" -> Scppi0Ks0,
    "AcpPi0Ks" -> Acppi0Ks0,
    "phid" -> phid0
  |>;

  purplePhiGridDeg = Range[phiScanRangeDeg[[1]], phiScanRangeDeg[[2]], phiStepDeg];
  purplePhiGridPlotDegFull = Select[purplePhiGridDeg, phiPlotRangeDeg[[1]] <= # <= phiPlotRangeDeg[[2]] &];

  centralAcceptedRootsByPhiRes = ComputeAcceptedRootsByPhi[purplePars, purpleObs, purplePhiGridDeg];
  If[centralAcceptedRootsByPhiRes === $Failed,
    Print["[ERROR] Failed to build central accepted roots for purple band."];
    Abort[]
  ];
  purpleCentralCurve = BuildCentralTargetBranch[
    centralAcceptedRootsByPhiRes["QByPhi"],
    purplePhiGridPlotDegFull
  ];
  If[purpleCentralCurve === {},
    Print["[ERROR] Failed to identify central purple target branch."];
    Abort[]
  ];

  purpleCentralQByKey = Association @ Map[PhiKey[#[[1]]] -> #[[2]] &, purpleCentralCurve];
  purplePhiSupportDeg = purpleCentralCurve[[All, 1]];
  purpleQByPhi = AssociationMap[{} &, PhiKey /@ purplePhiSupportDeg];

  purpleSampler = Function[{}, SampleParametersMC[purpleMcParamSpecs]];
  sampledPurpleParamSets = Table[purpleSampler[], {purpleNSamples}];

  Print["[INFO] Running purple Monte Carlo with nSamples = ", purpleNSamples];

  purpleMatchedSamples = 0;
  Do[
    Module[{sampleState, sampleRootsRes, key, phiDeg, candidates, bestQ, centralQ, matchedAnyQ},
      sampleState = BuildPurpleSampleState[sampledPurpleParamSets[[i]]];
      sampleRootsRes = ComputeAcceptedRootsByPhi[sampleState["Pars"], sampleState["Obs"], purplePhiGridDeg];
      If[sampleRootsRes === $Failed, Continue[]];

      matchedAnyQ = False;
      Do[
        key = PhiKey[phiDeg];
        candidates = Sort @ Select[sampleRootsRes["QByPhi"][key], NumericQ];
        candidates = Select[candidates, qPlotRange[[1]] <= # <= qPlotRange[[2]] &];
        If[!KeyExistsQ[purpleCentralQByKey, key] || candidates === {}, Continue[]];
        centralQ = purpleCentralQByKey[key];
        bestQ = First @ SortBy[candidates, Abs[# - centralQ] &];
        If[useDQMaxCut && Abs[bestQ - centralQ] >= dqMax, Continue[]];
        purpleQByPhi[key] = Join[purpleQByPhi[key], {bestQ}];
        matchedAnyQ = True;
        ,
        {phiDeg, purplePhiGridDeg}
      ];
      If[matchedAnyQ, purpleMatchedSamples++];
    ],
    {i, Length[sampledPurpleParamSets]}
  ];

  purpleEnvelopeRows = Cases[
    Table[
      Module[{phiDeg, vals},
        phiDeg = purplePhiSupportDeg[[i]];
        vals = purpleQByPhi[PhiKey[phiDeg]];
        If[vals =!= {} && !VectorQ[vals, NumericQ],
          Print["[ERROR] Non-numeric purple q-values at phi bin ", phiDeg, ": ", vals];
          Abort[]
        ];
        If[Length[vals] >= minAcceptedPerPhiBin,
          <|
            "PhiDeg" -> N @ phiDeg,
            "QLow16" -> N @ Quantile[vals, 0.16],
            "QMedian" -> N @ Quantile[vals, 0.50],
            "QHigh84" -> N @ Quantile[vals, 0.84],
            "NAcceptedAtPhi" -> Length[vals]
          |>,
          Nothing
        ]
      ],
      {i, Length[purplePhiSupportDeg]}
    ],
    Except[Nothing]
  ];
  If[purpleEnvelopeRows === {},
    Print["[ERROR] Purple envelope construction produced no valid bins."];
    Abort[]
  ];

  PrintPurpleEnvelopeDiagnostics[purpleCentralCurve, purpleEnvelopeRows];

  purpleBandGraphic = MakePurpleBandGraphic[purpleEnvelopeRows];
  purpleBoundaryPlots = MakePurpleBoundaryPlots[purpleEnvelopeRows];

  Print["[INFO] Purple matched samples = ", purpleMatchedSamples, "/", purpleNSamples];

  (* -------- annotations -------- *)
  smGraphic = MakeSMGraphic[];
  scpLabelGraphic = Graphics[{
    Text[
      Style[
        Column[
          {
            Row[{Subscript["S", "CP"], " current data"}],
            Row[{Superscript[Subscript["S", "CP"], "\[Pi]0 Ks"], " = 0.64"}]
          },
          Spacings -> 0.2
        ],
        17,
        purpleLineColor
      ],
      {12, 2.12}
    ]
  }];

  (* -------- final combined figure -------- *)
  fig = Show[
    blueGreenBandGraphics,
    purpleBandGraphic,
    Sequence @@ blueGreenCurvePlots,
    Sequence @@ purpleBoundaryPlots,
    ListLinePlot[
      {purpleCentralCurve},
      PlotStyle -> Directive[purpleLineColor, Dashed, Thick],
      InterpolationOrder -> 1
    ],
    smGraphic,
    scpLabelGraphic,
    Frame -> True,
    Axes -> False,
    FrameStyle -> Directive[Black, 14],
    FrameTicksStyle -> Directive[Black, 12],
    FrameLabel -> {
      Style["\[Phi] [deg]", 16, Black],
      Style["q", 16, Black]
    },
    PlotRange -> plotRange,
    PlotRangeClipping -> True,
    PlotRangePadding -> Scaled[0.02],
    ImageSize -> imageSize,
    AspectRatio -> aspectRatio,
    Background -> White
  ];

  (* -------- export -------- *)
  outRoot = FileNameJoin[{projectDir, "results", "q_phi"}];
  If[!DirectoryQ[outRoot],
    CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
  ];

  runTag = DateString[{"Year", "Month", "Day", "_", "Hour", "Minute", "Second"}];
  outDir = FileNameJoin[{outRoot, "run_" <> runTag}];
  If[!DirectoryQ[outDir],
    CreateDirectory[outDir, CreateIntermediateDirectories -> True]
  ];

  base = FileNameJoin[{outDir, "q_phi_scan_combined"}];
  Export[base <> ".png", fig, ImageResolution -> 300];
  Export[base <> ".pdf", fig];
  Export[
    base <> ".txt",
    StringRiffle[
      {
        "combined q-phi plot export",
        "blueGreenNSamples = " <> ToString[qphiNSamples],
        "blueGreenValidSamples = " <> ToString[qphiValidSamples],
        "purpleNSamples = " <> ToString[purpleNSamples],
        "purpleMatchedSamples = " <> ToString[purpleMatchedSamples],
        "blueGreenBranches = " <> ToString[blueGreenBranches, InputForm]
      },
      "\n"
    ],
    "String"
  ];

  Print["[OK] Saved: ", base <> ".png"];
  Print["[OK] Saved: ", base <> ".pdf"];
  Print["[OK] Saved: ", base <> ".txt"];

  Exit[0];

,
  Print["[ERROR] Combined q-phi plotting script failed."];
  Exit[1];
];
