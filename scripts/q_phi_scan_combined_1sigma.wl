(* ::Package:: *)

(*============== Script for final combined q-\[Phi] plot (1-sigma FD) ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
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

  (* -------- purple-band settings -------- *)
  phiPlotRangeDeg = {-75, 78};
  qPlotRange = {0, 3.05};
  phiScanRangeDeg = {-85, 85};
  phiStepDeg = 1.0;
  qMinCut = 0.05;
  residualTolerance = 10^-7;
  denominatorEps = 10^-8;

  (* -------- generic geometry helpers -------- *)
  ClearAll[ValidPoint2DQ];
  ValidPoint2DQ[pt_] := MatchQ[pt, {_?NumericQ, _?NumericQ}];

  ClearAll[SplitCurveByJump];
  SplitCurveByJump[data_List, jump_ : 35] :=
   Module[{clean, pieces = {}, current, dphi},
    clean = Cases[data, {_?NumericQ, _?NumericQ}];
    If[clean === {}, Return[{}]];

    current = {First[clean]};
    Do[
      dphi = Abs[clean[[i, 1]] - clean[[i - 1, 1]]];
      If[dphi > jump,
        AppendTo[pieces, current];
        current = {clean[[i]]},
        AppendTo[current, clean[[i]]]
      ],
      {i, 2, Length[clean]}
    ];

    AppendTo[pieces, current];
    pieces
   ];

  ClearAll[PhiKey];
  PhiKey[phiDeg_] := Round[1000 phiDeg];

  ClearAll[BranchKey];
  BranchKey[{i_, s_}] := "B" <> ToString[i] <> "_" <> ToString[s];

  ClearAll[MakeBandPolygonsFromUpperLower];
  MakeBandPolygonsFromUpperLower[upper_List, lower_List, jump_ : 40] :=
   Module[{upperPieces, lowerPieces, nPieces},
    upperPieces = SplitCurveByJump[upper, jump];
    lowerPieces = SplitCurveByJump[lower, jump];
    nPieces = Min[Length[upperPieces], Length[lowerPieces]];
    DeleteCases[
      Table[
        If[
          Length[upperPieces[[i]]] >= 2 && Length[lowerPieces[[i]]] >= 2,
          Polygon[Join[upperPieces[[i]], Reverse[lowerPieces[[i]]]]],
          Nothing
        ],
        {i, nPieces}
      ],
      Nothing
    ]
   ];

  ClearAll[MakeQPhiFDBandGraphics];
  Options[MakeQPhiFDBandGraphics] = {
    "KeepBranches" -> Automatic,
    "Jump" -> 40,
    "Opacity" -> 0.22
  };

  MakeQPhiFDBandGraphics[fdBands_Association, branchColors_List, OptionsPattern[]] :=
   Module[{keepBranches, jump, alpha},
    keepBranches = Replace[OptionValue["KeepBranches"], Automatic -> QPhiDefaultKeepBranches[]];
    jump = OptionValue["Jump"];
    alpha = OptionValue["Opacity"];

    Graphics[
      Flatten @ Table[
        Module[{brKey, i, bandAssoc, polys},
          i = br[[1]];
          brKey = BranchKey[br];
          If[!KeyExistsQ[fdBands, brKey],
            {},
            bandAssoc = fdBands[brKey];
            polys = MakeBandPolygonsFromUpperLower[bandAssoc["Upper"], bandAssoc["Lower"], jump];
            If[polys === {}, {}, {Directive[branchColors[[i]], Opacity[alpha], EdgeForm[None]], polys}]
          ]
        ],
        {br, keepBranches}
      ]
    ]
   ];

  ClearAll[PrintRepresentativeCurvePoints];
  PrintRepresentativeCurvePoints[label_, curve_List] :=
   Module[{idxs},
    If[Length[curve] == 0,
      Print["[DEBUG] ", label, " has no points."];
      Return[];
    ];
    idxs = DeleteDuplicates @ Round /@ {1, Max[1, Length[curve]/2], Length[curve]};
    Print["[DEBUG] ", label, " representative central points = ", N[curve[[idxs]], 8]];
   ];

  ClearAll[PrintRepresentativeBandWidths];
  PrintRepresentativeBandWidths[label_, upper_List, lower_List] :=
   Module[{n, idxs, rows},
    n = Min[Length[upper], Length[lower]];
    If[n == 0,
      Print["[DEBUG] ", label, " has no band points."];
      Return[];
    ];
    idxs = DeleteDuplicates @ Round /@ {1, Max[1, n/2], n};
    rows = Table[
      <|
        "PhiDeg" -> N[upper[[i, 1]], 8],
        "CenterQ" -> N[(upper[[i, 2]] + lower[[i, 2]])/2, 8],
        "SigmaQ" -> N[(upper[[i, 2]] - lower[[i, 2]])/2, 8]
      |>,
      {i, idxs}
    ];
    Print["[DEBUG] ", label, " representative widths = ", rows];
   ];

  (* -------- local helpers for purple workflow -------- *)
  ClearAll[BuildPurpleInputState];
  BuildPurpleInputState[s_Association] :=
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

  ClearAll[BuildPurpleMatchedCurve];
  BuildPurpleMatchedCurve[inputAssoc_Association, phiGridDeg_List, phiSupportDeg_List, centralQByKey_Association] :=
   Module[{state, sampleRootsRes, curve},
    state = BuildPurpleInputState[inputAssoc];
    sampleRootsRes = ComputeAcceptedRootsByPhi[state["Pars"], state["Obs"], phiGridDeg];
    If[sampleRootsRes === $Failed, Return[$Failed]];

    curve = Table[
      Module[{key, candidates, centralQ, bestQ},
        key = PhiKey[phiDeg];
        If[!KeyExistsQ[centralQByKey, key], Return[Missing["NoCentralReference"]]];
        candidates = Sort @ Select[sampleRootsRes["QByPhi"][key], NumericQ];
        candidates = Select[candidates, qPlotRange[[1]] <= # <= qPlotRange[[2]] &];
        If[candidates === {}, Return[Missing["NoCandidate"]]];
        centralQ = centralQByKey[key];
        bestQ = First @ SortBy[candidates, Abs[# - centralQ] &];
        {phiDeg, bestQ}
      ],
      {phiDeg, phiSupportDeg}
    ];

    If[!VectorQ[curve, ValidPoint2DQ], Return[$Failed]];
    curve
   ];

  ClearAll[BuildPurpleBandRowsFD];
  BuildPurpleBandRowsFD[centerInputs_Association, paramSpecs_List, phiGridDeg_List, phiSupportDeg_List, centralQByKey_Association] :=
   Module[{centerCurve, names, variations, rows},
    centerCurve = BuildPurpleMatchedCurve[centerInputs, phiGridDeg, phiSupportDeg, centralQByKey];
    If[centerCurve === $Failed, Return[$Failed]];

    names = paramSpecs[[All, 1]];
    variations = Association @ Table[
      Module[{name, sig, upIn, dnIn, upCurve, dnCurve},
        name = paramSpecs[[i, 1]];
        sig = paramSpecs[[i, 3]];
        upIn = centerInputs;
        dnIn = centerInputs;
        upIn[name] = upIn[name] + sig;
        dnIn[name] = dnIn[name] - sig;
        upCurve = BuildPurpleMatchedCurve[upIn, phiGridDeg, phiSupportDeg, centralQByKey];
        dnCurve = BuildPurpleMatchedCurve[dnIn, phiGridDeg, phiSupportDeg, centralQByKey];
        name -> <|"Up" -> upCurve, "Down" -> dnCurve|>
      ],
      {i, Length[paramSpecs]}
    ];

    rows = Cases[
      Table[
        Module[{phiDeg, qCenter, halfDiffs, sigmaQ},
          phiDeg = centerCurve[[i, 1]];
          qCenter = centerCurve[[i, 2]];
          halfDiffs = Table[
            If[
              AssociationQ[variations[name]] &&
              ListQ[variations[name]["Up"]] &&
              ListQ[variations[name]["Down"]] &&
              Length[variations[name]["Up"]] >= i &&
              Length[variations[name]["Down"]] >= i &&
              ValidPoint2DQ[variations[name]["Up"][[i]]] &&
              ValidPoint2DQ[variations[name]["Down"][[i]]],
              (variations[name]["Up"][[i, 2]] - variations[name]["Down"][[i, 2]])/2,
              Missing["NoFiniteDifference"]
            ],
            {name, names}
          ];
          If[MemberQ[halfDiffs, _Missing], Return[Nothing]];
          sigmaQ = Sqrt[Total[halfDiffs^2]];
          <|
            "PhiDeg" -> N @ phiDeg,
            "QCentral" -> N @ qCenter,
            "QLow16" -> N @ (qCenter - sigmaQ),
            "QMedian" -> N @ qCenter,
            "QHigh84" -> N @ (qCenter + sigmaQ)
          |>
        ],
        {i, Length[centerCurve]}
      ],
      Except[Nothing]
    ];

    <|
      "CenterCurve" -> centerCurve,
      "BandRows" -> rows,
      "Variations" -> variations
    |>
   ];

  ClearAll[MakePurpleBandGraphic];
  MakePurpleBandGraphic[purpleEnvelopeRows_List] :=
   Module[{upper, lower},
    upper = Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QHigh84"]]}];
    lower = Reverse @ Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QLow16"]]}];
    Graphics[{
      Directive[purpleFillColor, Opacity[0.35], EdgeForm[None]],
      Polygon[Join[upper, lower]]
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

  (* -------- blue/green workflow with FD 1-sigma -------- *)
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

  Print["[INFO] Building blue/green q-phi bands with finite-difference 1-sigma and nGrid = ", qphiNGrid];

  centerRes = ComputeQPhiBranchData[centralPars, qphiNGrid, MKp0];
  If[centerRes === $Failed,
    Print["[ERROR] Central q-phi branch construction failed."];
    Abort[];
  ];
  branchData = centerRes["BranchData"];

  qphiFDBands = Association @ Flatten @ Table[
    Module[{brKey, singleFDRes},
      brKey = BranchKey[br];
      singleFDRes = BuildCurveBandsFD[
        Function[in, ComputeQPhiBranchData[in, qphiNGrid, MKp0]],
        Function[res,
          <|
            brKey -> Cases[res["BranchData"][[br[[1]], br[[2]]]], {_?NumericQ, _?NumericQ}]
          |>
        ],
        centralPars,
        qphiParamSpecs
      ];
      {brKey -> singleFDRes["Bands"][brKey]}
    ],
    {br, blueGreenBranches}
  ];

  blueGreenBandGraphics = MakeQPhiFDBandGraphics[
    qphiFDBands,
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

  PrintRepresentativeCurvePoints["Green central curve", Cases[branchData[[2, 1]], {_?NumericQ, _?NumericQ}]];
  PrintRepresentativeCurvePoints["Blue central curve", Cases[branchData[[3, 1]], {_?NumericQ, _?NumericQ}]];
  PrintRepresentativeBandWidths["Green band", qphiFDBands[BranchKey[{2, 1}]]["Upper"], qphiFDBands[BranchKey[{2, 1}]]["Lower"]];
  PrintRepresentativeBandWidths["Blue band", qphiFDBands[BranchKey[{3, 1}]]["Upper"], qphiFDBands[BranchKey[{3, 1}]]["Lower"]];

  (* -------- purple workflow with FD 1-sigma -------- *)
  purpleParamSpecs = {
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

  purpleCentralInputs = <|
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

  Print["[INFO] Building purple 1-sigma finite-difference envelope..."];

  purpleFDRes = BuildPurpleBandRowsFD[
    purpleCentralInputs,
    purpleParamSpecs,
    purplePhiGridDeg,
    purplePhiSupportDeg,
    purpleCentralQByKey
  ];
  If[purpleFDRes === $Failed,
    Print["[ERROR] Purple finite-difference envelope construction failed."];
    Abort[]
  ];

  purpleEnvelopeRows = purpleFDRes["BandRows"];
  If[purpleEnvelopeRows === {},
    Print["[ERROR] Purple envelope construction produced no valid bins."];
    Abort[]
  ];

  purpleBandGraphic = MakePurpleBandGraphic[purpleEnvelopeRows];
  purpleBoundaryPlots = MakePurpleBoundaryPlots[purpleEnvelopeRows];

  PrintRepresentativeCurvePoints["Purple central curve", purpleCentralCurve];
  PrintRepresentativeBandWidths[
    "Purple band",
    Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QHigh84"]]}],
    Transpose[{purpleEnvelopeRows[[All, "PhiDeg"]], purpleEnvelopeRows[[All, "QLow16"]]}]
  ];

  Print["[INFO] Purple central support points = ", Length[purpleCentralCurve]];
  Print["[INFO] Purple 1-sigma envelope bins = ", Length[purpleEnvelopeRows]];

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
      {24, 2.12}
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

  base = FileNameJoin[{outDir, "q_phi_scan_combined_1sigma"}];
  Export[base <> ".png", fig, ImageResolution -> 300];
  Export[base <> ".pdf", fig];
  Export[
    base <> ".txt",
    StringRiffle[
      {
        "combined q-phi plot export (1-sigma finite difference)",
        "blueGreenMethod = BuildCurveBandsFD",
        "purpleMethod = finite-difference target-branch propagation",
        "blueGreenBranches = " <> ToString[blueGreenBranches, InputForm],
        "purpleCentralSupportPoints = " <> ToString[Length[purpleCentralCurve]],
        "purpleEnvelopeBins = " <> ToString[Length[purpleEnvelopeRows]]
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
  Print["[ERROR] Combined q-phi 1-sigma plotting script failed."];
  Exit[1];
];
