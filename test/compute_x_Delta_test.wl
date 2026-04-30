(* ::Package:: *)

(*============== Test script for x and \[CapitalDelta] computation (MC uncertainty) ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "uncertainty_MC.wl"}]];

  (* -------- local MC helpers -------- *)
  ClearAll[MCQuantileSigma];
  MCQuantileSigma[vals_List, qlo_ : 0.16, qhi_ : 0.84] :=
   Module[{clean},
    clean = Cases[vals, _?NumericQ];
    If[Length[clean] < 2,
      Indeterminate,
      (Quantile[clean, qhi] - Quantile[clean, qlo])/2
    ]
   ];

  ClearAll[BuildCurveBandMC];
  BuildCurveBandMC[centerCurve_List, sampleCurves_List, qlo_ : 0.16, qhi_ : 0.84] :=
   Module[{xGrid, yCenter, yUpper, yLower, yVals},
    xGrid = centerCurve[[All, 1]];
    yCenter = centerCurve[[All, 2]];

    yUpper = Table[
      yVals = Cases[sampleCurves[[All, k, 2]], _?NumericQ];
      If[Length[yVals] < 2, yCenter[[k]], Quantile[yVals, qhi]],
      {k, Length[centerCurve]}
    ];

    yLower = Table[
      yVals = Cases[sampleCurves[[All, k, 2]], _?NumericQ];
      If[Length[yVals] < 2, yCenter[[k]], Quantile[yVals, qlo]],
      {k, Length[centerCurve]}
    ];

    <|
      "Upper" -> Transpose[{xGrid, yUpper}],
      "Lower" -> Transpose[{xGrid, yLower}]
    |>
   ];

  (* -------- load cached R00 / Rpm -------- *)
  {R000, R00Err} =
    Import[FileNameJoin[{projectDir, "results", "B_pipi_R00_Rpm", "R00.mx"}]];

  {Rpm0, RpmErr} =
    Import[FileNameJoin[{projectDir, "results", "B_pipi_R00_Rpm", "Rpm.mx"}]];

  (* -------- inputs -------- *)
  inputs = <|
    "A"     -> Acppippim0,
    "S"     -> Scppippim0,
    "gamma" -> gamma0,
    "phid"  -> phid0,
    "Rpm"   -> Rpm0,
    "R00"   -> R000,
    "Acp00" -> Acppi0pi00
  |>;

  paramSpecs = {
    {"A",     Acppippim0, AcppippimErr, "Any"},
    {"S",     Scppippim0, ScppippimErr, "Any"},
    {"gamma", gamma0,     gammaErr,     "Any"},
    {"phid",  phid0,      phidErr,      "Any"},
    {"Rpm",   Rpm0,       RpmErr,       "Positive"},
    {"R00",   R000,       R00Err,       "Positive"}
  };

  deltaStepDeg = 1;
  nSamples = 1000;
  qloMC = 0.16;
  qhiMC = 0.84;

  (* -------- compute -------- *)
  centerRes = ComputeXDeltaCenter[inputs, "DeltaStepDeg" -> deltaStepDeg];

  compute = Function[in, ComputeXDeltaCenter[in, "DeltaStepDeg" -> deltaStepDeg]];
  sampler = Function[{}, Join[inputs, SampleParametersMC[paramSpecs]]];

  sampleResults = GenerateMCSamples[
    compute,
    sampler,
    nSamples
  ];

  nValidSamples = Length[sampleResults];

  sampleCenters = sampleResults[[All, "Center"]];
  sampleCurvesRpm = sampleResults[[All, "Curves", "Rpm"]];
  sampleCurvesR00 = sampleResults[[All, "Curves", "R00"]];

  uncRes = <|
    "CenterResult" -> centerRes,
    "Scalars" -> <|
      "Center" -> centerRes["Center"],
      "Sigma" -> <|
        "x" -> MCQuantileSigma[Lookup[sampleCenters, "x"], qloMC, qhiMC],
        "Delta" -> MCQuantileSigma[Lookup[sampleCenters, "Delta"], qloMC, qhiMC],
        "d" -> MCQuantileSigma[Lookup[sampleCenters, "d"], qloMC, qhiMC],
        "theta" -> MCQuantileSigma[Lookup[sampleCenters, "theta"], qloMC, qhiMC]
      |>
    |>,
    "Curves" -> <|
      "Center" -> centerRes["Curves"],
      "Bands" -> <|
        "Rpm" -> BuildCurveBandMC[centerRes["Curves"]["Rpm"], sampleCurvesRpm, qloMC, qhiMC],
        "R00" -> BuildCurveBandMC[centerRes["Curves"]["R00"], sampleCurvesR00, qloMC, qhiMC]
      |>,
      "Samples" -> <|
        "Rpm" -> sampleCurvesRpm,
        "R00" -> sampleCurvesR00
      |>
    |>,
    "MCSettings" -> <|
      "NSamplesRequested" -> nSamples,
      "NSamplesValid" -> nValidSamples,
      "LowerQuantile" -> qloMC,
      "UpperQuantile" -> qhiMC
    |>
  |>;

  (* -------- numbers -------- *)
  x0 = centerRes["Center"]["x"];
  dx = uncRes["Scalars"]["Sigma"]["x"];

  Delta0deg = centerRes["Center"]["Delta"]/Degree;
  dDeltadeg = uncRes["Scalars"]["Sigma"]["Delta"]/Degree;

  Print["x = ", N[x0, 8], " \[PlusMinus] ", N[dx, 4]];
  Print["\[CapitalDelta] = ", N[Delta0deg, 8], "\[Degree] \[PlusMinus] ", N[dDeltadeg, 4], "\[Degree]"];

  (* -------- plot -------- *)
  fig = PlotXDelta[centerRes, uncRes, "Title" -> "x ~ \[CapitalDelta]"];

  outRoot = FileNameJoin[{projectDir, "test", "results_x_Delta_MC"}];
  If[!DirectoryQ[outRoot],
    CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
  ];

  runTag = DateString[{"Year", "Month", "Day", "_", "Hour", "Minute", "Second"}];
  outDir = FileNameJoin[{outRoot, "run_" <> runTag}];

  If[!DirectoryQ[outDir],
    CreateDirectory[outDir, CreateIntermediateDirectories -> True]
  ];

  base = FileNameJoin[{outDir, "xDelta_MC"}];

  (* -------- export figure (PNG only) -------- *)
  Export[base <> ".png", fig, ImageResolution -> 300];

  (* -------- export numbers -------- *)

  (* 1) .mx : machine-readable *)
  Export[
    base <> ".mx",
    <|
      "x" -> x0,
      "dx" -> dx,
      "DeltaDeg" -> Delta0deg,
      "dDeltaDeg" -> dDeltadeg,
      "DeltaStepDeg" -> deltaStepDeg,
      "MCSettings" -> uncRes["MCSettings"]
    |>
  ];

  (* 2) .txt : human-readable *)
  Export[
    base <> ".txt",
    StringRiffle[
      {
        "x = " <> ToString @ N[x0, 12] <>
          " +/- " <> ToString @ N[dx, 8],

        "Delta(deg) = " <> ToString @ N[Delta0deg, 12] <>
          " +/- " <> ToString @ N[dDeltadeg, 8],

        "DeltaStepDeg = " <> ToString[deltaStepDeg],
        "MC samples requested = " <> ToString[nSamples],
        "MC samples valid = " <> ToString[nValidSamples],
        "MC quantiles = {" <> ToString[qloMC] <> ", " <> ToString[qhiMC] <> "}"
      },
      "\n"
    ],
    "String"
  ];

  Print["[OK] Saved: ", base <> ".png"];
  Print["[OK] Saved: ", base <> ".mx"];
  Print["[OK] Saved: ", base <> ".txt"];
  Print["[MC] Settings: nSamples = ", nSamples, ", validSamples = ", nValidSamples,
    ", quantiles = {", qloMC, ", ", qhiMC, "}"];
  Print["[MC] Output directory: ", outDir];

  Exit[0];

,
  Print["[ERROR] Test script failed."];
  Exit[1];
];
