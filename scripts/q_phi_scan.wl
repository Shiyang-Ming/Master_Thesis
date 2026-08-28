(* ::Package:: *)

(*============== Script for q-\[Phi] scan computation ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "uncertainty_MC.wl"}]];
  Get[FileNameJoin[{projectDir, "src", "q_phi_scan.wl"}]];

  (* -------- inputs -------- *)
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

  paramSpecs = {
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

  nGrid = 1350;
  nSamples = 1000;
  branchColors = QPhiDefaultBranchColors[];
  keepBranches = QPhiDefaultKeepBranches[];
  rcColor = RGBColor[0.0, 0.65, 0.68];
  phiGridRc = Subdivide[-180, 180, 1600];

  Print["[INFO] Building central q-phi branches with nGrid = ", nGrid];

  (* -------- central branch data -------- *)
  centerRes = ComputeQPhiBranchData[centralPars, nGrid, MKp0];
  If[centerRes === $Failed,
    Print["[ERROR] Central branch construction failed."];
    Abort[];
  ];

  branchData = centerRes["BranchData"];

  (* -------- Monte Carlo sampling -------- *)
  compute = Function[in, ComputeQPhiBranchData[in, nGrid, MKp0]];
  sampler = Function[{}, SampleParametersMC[paramSpecs]];

  Print["[INFO] Running Monte Carlo with nSamples = ", nSamples];

  sampleResults = GenerateMCSamples[
    compute,
    sampler,
    nSamples,
    Function[res, res["BranchData"]]
  ];

  nValidSamples = Length[sampleResults];
  Print["[INFO] Valid MC samples = ", nValidSamples, "/", nSamples];

  If[nValidSamples == 0,
    Print["[ERROR] No valid Monte Carlo samples were produced."];
    Abort[];
  ];

  (* -------- Rc constraint -------- *)
  rcRes = ComputeRcObservable[
    BrBppi0Kp0, BrBppi0KpErr,
    BrBppipK00, BrBppipK0Err
  ];

  Rc0 = rcRes["Value"];
  RcErr = rcRes["Error"];

  Print["[INFO] Rc = ", N[Rc0, 6], " +/- ", N[RcErr, 6]];

  rcGraphics = MakeQPhiRcPlots[
    phiGridRc, Rc0, RcErr,
    gamma0, rhoc0, thetac0, rc0, deltac0,
    rcColor
  ];

  (* -------- plot -------- *)
  bandGraphics = MakeQPhiBandGraphics[
    branchData,
    sampleResults,
    branchColors,
    "KeepBranches" -> keepBranches,
    "Jump" -> 40,
    "Opacity" -> 0.22
  ];

  branchPlots = MakeQPhiBranchPlots[
    branchData,
    branchColors,
    "KeepBranches" -> keepBranches,
    "Jump" -> 40
  ];

  fig = PlotQPhiScan[
    bandGraphics,
    branchPlots,
    rcGraphics["Plots"],
    rcGraphics["LabelGraphic"],
    "PlotRange" -> {{-180, 180}, {0, 3}}
  ];

  (* -------- output directory -------- *)
  outRoot = FileNameJoin[{projectDir, "results", "q_phi"}];
  If[!DirectoryQ[outRoot],
    CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
  ];

  runTag = DateString[{"Year", "Month", "Day", "_", "Hour", "Minute", "Second"}];
  outDir = FileNameJoin[{outRoot, "run_" <> runTag}];
  If[!DirectoryQ[outDir],
    CreateDirectory[outDir, CreateIntermediateDirectories -> True]
  ];

  base = FileNameJoin[{outDir, "q_phi_scan"}];

  (* -------- export outputs -------- *)
  Export[base <> ".png", fig, ImageResolution -> 300];

  Export[
    base <> ".mx",
    <|
      "BranchData" -> branchData,
      "SampleResults" -> sampleResults,
      "Rc" -> <|"Value" -> Rc0, "Error" -> RcErr|>,
      "Settings" -> <|
        "NGrid" -> nGrid,
        "NSamplesRequested" -> nSamples,
        "NSamplesValid" -> nValidSamples,
        "KeepBranches" -> keepBranches
      |>
    |>
  ];

  Export[
    base <> ".txt",
    StringRiffle[
      {
        "q-phi scan export",
        "Rc = " <> ToString @ N[Rc0, 12] <> " +/- " <> ToString @ N[RcErr, 8],
        "nGrid = " <> ToString[nGrid],
        "MC samples requested = " <> ToString[nSamples],
        "MC samples valid = " <> ToString[nValidSamples],
        "KeepBranches = " <> ToString[keepBranches, InputForm]
      },
      "\n"
    ],
    "String"
  ];

  Print["[OK] Saved: ", base <> ".png"];
  Print["[OK] Saved: ", base <> ".mx"];
  Print["[OK] Saved: ", base <> ".txt"];
  Print["[INFO] Output directory: ", outDir];

  Exit[0];

,
  Print["[ERROR] Script failed."];
  Exit[1];
];
