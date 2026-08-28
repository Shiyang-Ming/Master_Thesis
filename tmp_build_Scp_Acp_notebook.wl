codeCell[s_String] := Cell[s, "Input"];

initializationCode = "ClearAll[\"Global`*\"];
SetDirectory[NotebookDirectory[]];
projectRoot = FileNameJoin[{NotebookDirectory[], \"..\"}];";

loadInputsCode = "Get[FileNameJoin[{\"..\", \"src\", \"inputs_exp.wl\"}]];
Get[FileNameJoin[{\"..\", \"src\", \"parameters.wl\"}]];";

definitionsCode = "ClearAll[realClean];
realClean[x_, tol_ : 10^-10] := Module[{xc = Chop[x, tol]},
  If[NumericQ[xc] && Chop[Im[xc], tol] == 0, Re[xc], Indeterminate]
];

ClearAll[safeArcCos];
safeArcCos[x_] := ArcCos[Clip[realClean[x], {-1, 1}]];

ClearAll[PropagateScalarsFD];
Options[PropagateScalarsFD] = {
  \"StepScale\" -> 1,
  \"Combine\" -> (Sqrt[Total[#^2]] &)
};

PropagateScalarsFD[
  compute_,
  scalarExtractor_,
  inputs_Association,
  paramSpecs_List,
  opts : OptionsPattern[]
] := Module[
  {centerRes, center, names, stepScale, vars, combine, sigmas},
  stepScale = OptionValue[\"StepScale\"];
  combine = OptionValue[\"Combine\"];
  centerRes = compute[inputs];
  center = scalarExtractor[centerRes];
  names = paramSpecs[[All, 1]];
  vars = Association@Table[
    Module[{name = paramSpecs[[i, 1]], sig = paramSpecs[[i, 3]]*stepScale, upIn, dnIn, upC, dnC},
      upIn = inputs;
      dnIn = inputs;
      upIn[name] = upIn[name] + sig;
      dnIn[name] = dnIn[name] - sig;
      upC = scalarExtractor @ compute[upIn];
      dnC = scalarExtractor @ compute[dnIn];
      name -> <|\"Up\" -> upC, \"Down\" -> dnC, \"HalfDiff\" -> (upC - dnC)/2|>
    ],
    {i, Length[paramSpecs]}
  ];
  sigmas = Association @ Map[
    Function[key, key -> combine@Table[vars[n][\"HalfDiff\"][key], {n, names}]],
    Keys[center]
  ];
  <|
    \"CenterResult\" -> centerRes,
    \"Center\" -> center,
    \"Sigma\" -> sigmas,
    \"Variations\" -> vars
  |>
];

ClearAll[BuildCurveBandsFD];
Options[BuildCurveBandsFD] = {
  \"StepScale\" -> 1,
  \"Combine\" -> (Sqrt[Total[#^2]] &)
};

BuildCurveBandsFD[
  compute_,
  curveExtractor_,
  inputs_Association,
  paramSpecs_List,
  opts : OptionsPattern[]
] := Module[
  {centerRes, curves, curveNames, xGrid, yCenter, names, stepScale, combine, dY, sigY, bands},
  stepScale = OptionValue[\"StepScale\"];
  combine = OptionValue[\"Combine\"];
  centerRes = compute[inputs];
  curves = curveExtractor[centerRes];
  curveNames = Keys[curves];
  names = paramSpecs[[All, 1]];
  xGrid = curves[curveNames[[1]]][[All, 1]];
  yCenter = AssociationMap[curves[#][[All, 2]] &, curveNames];
  dY = Association@Table[
    Module[{pname = paramSpecs[[i, 1]], sig = paramSpecs[[i, 3]]*stepScale, upIn, dnIn, upCur, dnCur},
      upIn = inputs;
      dnIn = inputs;
      upIn[pname] = upIn[pname] + sig;
      dnIn[pname] = dnIn[pname] - sig;
      upCur = curveExtractor @ compute[upIn];
      dnCur = curveExtractor @ compute[dnIn];
      pname -> AssociationMap[(upCur[#][[All, 2]] - dnCur[#][[All, 2]])/2 &, curveNames]
    ],
    {i, Length[paramSpecs]}
  ];
  sigY = AssociationMap[
    Function[cname, Sqrt @ Total @ (Table[dY[p][cname], {p, names}]^2)],
    curveNames
  ];
  bands = AssociationMap[
    Function[cname,
      <|
        \"Upper\" -> Transpose[{xGrid, yCenter[cname] + sigY[cname]}],
        \"Lower\" -> Transpose[{xGrid, yCenter[cname] - sigY[cname]}]
      |>
    ],
    curveNames
  ];
  <|
    \"CenterResult\" -> centerRes,
    \"Curves\" -> curves,
    \"Bands\" -> bands,
    \"SigmaY\" -> sigY,
    \"dY\" -> dY
  |>
];

ClearAll[ComputeUncertaintyFD];
Options[ComputeUncertaintyFD] = Join[Options[PropagateScalarsFD], Options[BuildCurveBandsFD]];

ComputeUncertaintyFD[
  compute_,
  scalarExtractor_,
  curveExtractor_,
  inputs_Association,
  paramSpecs_List,
  opts : OptionsPattern[]
] := Module[{scalar, band},
  scalar = PropagateScalarsFD[
    compute, scalarExtractor, inputs, paramSpecs,
    \"StepScale\" -> OptionValue[\"StepScale\"],
    \"Combine\" -> OptionValue[\"Combine\"]
  ];
  band = BuildCurveBandsFD[
    compute, curveExtractor, inputs, paramSpecs,
    \"StepScale\" -> OptionValue[\"StepScale\"],
    \"Combine\" -> OptionValue[\"Combine\"]
  ];
  <|
    \"CenterResult\" -> scalar[\"CenterResult\"],
    \"Scalars\" -> <|
      \"Center\" -> scalar[\"Center\"],
      \"Sigma\" -> scalar[\"Sigma\"],
      \"Variations\" -> scalar[\"Variations\"]
    |>,
    \"Curves\" -> <|
      \"Center\" -> band[\"Curves\"],
      \"Bands\" -> band[\"Bands\"],
      \"SigmaY\" -> band[\"SigmaY\"],
      \"dY\" -> band[\"dY\"]
    |>
  |>
];

ClearAll[phi00sol];
phi00sol[input_Association, A_?NumericQ] := Module[
  {A00, Amp, A32, A00bar, Ampbar, A32bar, phi32, phi1, phi2, phi001, phi002, phi003, phi004},
  A00 = Sqrt[2 input[\"BrB0pi0K0\"] (1 - A)];
  Amp = Sqrt[input[\"BrB0pimKp\"] (1 - input[\"AcpPimKp\"])];
  A32 =
    input[\"RTC\"]*(input[\"Vus\"]/input[\"Vud\"])*
    Sqrt[2 input[\"BrBppippi0\"]]*
    (Exp[I input[\"gamma\"]] - input[\"q\"] Exp[I input[\"phi\"]]);
  A00bar = Sqrt[2 input[\"BrB0pi0K0\"] (1 + A)];
  Ampbar = Sqrt[input[\"BrB0pimKp\"] (1 + input[\"AcpPimKp\"])];
  A32bar =
    input[\"RTC\"]*(input[\"Vus\"]/input[\"Vud\"])*
    Sqrt[2 input[\"BrBppippi0\"]]*
    (Exp[-I input[\"gamma\"]] - input[\"q\"] Exp[-I input[\"phi\"]]);
  phi32 = Arg[A32]/Degree;
  phi1 = safeArcCos[
    (Abs[A00]^2 + Abs[A32]^2 - Abs[Amp]^2)/
    (2 Abs[A00] Abs[A32])
  ]/Degree;
  phi2 = safeArcCos[
    (Abs[A00bar]^2 + Abs[A32bar]^2 - Abs[Ampbar]^2)/
    (2 Abs[A00bar] Abs[A32bar])
  ]/Degree;
  phi001 = phi1 + phi2 - 2 phi32;
  phi002 = 360 - (phi2 + 2 phi32 - phi1);
  phi003 = 360 - phi1 - phi2 - 2 phi32;
  phi004 = 360 - (2 phi32 - phi2 + phi1);
  {phi001, phi002, phi003, phi004}
];

ClearAll[ScpBranch];
ScpBranch[input_Association, A_?NumericQ, k_Integer] := Module[{phi00k},
  phi00k = phi00sol[input, A][[k]];
  Sqrt[1 - A^2] * Sin[input[\"phid\"] - phi00k Degree]
];

ClearAll[ComputeScpAcpCenter];
Options[ComputeScpAcpCenter] = {
  \"ARange\" -> {-0.4, 0.30},
  \"AStep\" -> 0.001
};

ComputeScpAcpCenter[input_Association, OptionsPattern[]] := Module[
  {aRange, aStep, as, curves},
  aRange = OptionValue[\"ARange\"];
  aStep = OptionValue[\"AStep\"];
  as = Range[aRange[[1]], aRange[[2]], aStep];
  curves = Association@Table[
    With[{branch = k},
      \"Branch\" <> ToString[branch] ->
        Table[{A, realClean @ ScpBranch[input, A, branch]}, {A, as}]
    ],
    {k, 1, 4}
  ];
  <|\"GridA\" -> as, \"Curves\" -> curves|>
];

ClearAll[BuildScpAcpFigure];
Options[BuildScpAcpFigure] = {
  \"PlotRange\" -> {{-0.4, 0.30}, {-1.0, 1.0}}
};

BuildScpAcpFigure[centerRes_Association, unc_Association, inputs_Association, OptionsPattern[]] := Module[
  {
    curves, bands, pts1, pts2, pts3, pts4, sig1, sig2, sig3, sig4,
    plotRange, amin, amax, ymin, ymax, xMajor, xMinor, yMajor, yMinor,
    ticksX, ticksY, frameTicksTarget, xLabel, yLabel, axisLikeTarget,
    xL, xR, sumRuleBandEpilog, sumRuleTextEpilog,
    currentDataEpilog, currentDataTextEpilog,
    colBlue, colGreen, colOrange, colGray, leftPlot, rightPlot, imgSingle, imgRow
  },
  curves = centerRes[\"Curves\"];
  bands = unc[\"Curves\"][\"Bands\"];
  pts1 = curves[\"Branch1\"];
  pts2 = curves[\"Branch2\"];
  pts3 = curves[\"Branch3\"];
  pts4 = curves[\"Branch4\"];
  sig1 = bands[\"Branch1\"];
  sig2 = bands[\"Branch2\"];
  sig3 = bands[\"Branch3\"];
  sig4 = bands[\"Branch4\"];
  plotRange = OptionValue[\"PlotRange\"];
  {{amin, amax}, {ymin, ymax}} = plotRange;
  xMajor = Range[-0.3, 0.2, 0.1];
  xMinor = Range[amin, amax, 0.02];
  yMajor = Range[-1, 1, 0.5];
  yMinor = Range[ymin, ymax, 0.1];
  ClearAll[mkFrameTicks];
  mkFrameTicks[major_, minor_, labfmt_, majLen_, minLen_] := Module[{minorOnly, clean},
    clean[x_] := Chop[x, 10^-10];
    minorOnly = Complement[minor, major];
    Join[
      ({clean[#], labfmt[clean[#]], {majLen, 0}} &) /@ major,
      ({clean[#], \"\", {minLen, 0}} &) /@ minorOnly
    ]
  ];
  ticksX = mkFrameTicks[xMajor, xMinor, (NumberForm[#, {2, 1}] &), 0.020, 0.010];
  ticksY = mkFrameTicks[yMajor, yMinor, (NumberForm[#, {2, 1}] &), 0.020, 0.010];
  frameTicksTarget = {
    {ticksY, ({#[[1]], \"\", #[[3]]} &) /@ ticksY},
    {ticksX, ({#[[1]], \"\", #[[3]]} &) /@ ticksX}
  };
  xLabel = Style[
    TraditionalForm @ Superscript[Subscript[\"A\", \"CP\"], \"\[Pi]0 Ks\"],
    20, Black, FontFamily -> \"Times\"
  ];
  yLabel = Style[
    TraditionalForm @ Superscript[Subscript[\"S\", \"CP\"], \"\[Pi]0 Ks\"],
    20, Black, FontFamily -> \"Times\"
  ];
  axisLikeTarget = {
    Background -> White,
    Axes -> False,
    Frame -> True,
    FrameStyle -> Directive[Black, AbsoluteThickness[1.0]],
    FrameTicks -> frameTicksTarget,
    FrameTicksStyle -> Directive[Black, 14, FontFamily -> \"Times\"],
    TicksStyle -> Directive[Black, 14, FontFamily -> \"Times\"],
    PlotRange -> plotRange,
    PlotRangePadding -> 0,
    AspectRatio -> 0.70,
    ImagePadding -> {{70, 25}, {55, 20}},
    FrameLabel -> {xLabel, yLabel}
  };
  xL = Acppi0KsSM0 - Acppi0KsSMErr;
  xR = Acppi0KsSM0 + Acppi0KsSMErr;
  sumRuleBandEpilog = {
    {Directive[Red, Opacity[0.18]], Rectangle[{xL, ymin}, {xR, ymax}]},
    {Directive[Red, AbsoluteThickness[2.2]], Line[{{xL, ymin}, {xL, ymax}}]},
    {Directive[Red, AbsoluteThickness[2.2]], Line[{{xR, ymin}, {xR, ymax}}]},
    {Directive[Red, AbsoluteThickness[2.2], Dashing[{0.012, 0.012}]],
      Line[{{Acppi0KsSM0, ymin}, {Acppi0KsSM0, ymax}}]}
  };
  sumRuleTextEpilog = {
    Inset[Style[\"Sum rule\\nprediction\", 20, Red, FontFamily -> \"Times\"], {Acppi0KsSM0 + 0.1, 0.0}]
  };
  currentDataEpilog = {
    {Black, AbsolutePointSize[6], Point[{Acppi0Ks0, Scppi0Ks0}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0 - Acppi0KsErr, Scppi0Ks0}, {Acppi0Ks0 + Acppi0KsErr, Scppi0Ks0}}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0, Scppi0Ks0 - Scppi0KsErr}, {Acppi0Ks0, Scppi0Ks0 + Scppi0KsErr}}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0 - Acppi0KsErr, Scppi0Ks0 - 0.02}, {Acppi0Ks0 - Acppi0KsErr, Scppi0Ks0 + 0.02}}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0 + Acppi0KsErr, Scppi0Ks0 - 0.02}, {Acppi0Ks0 + Acppi0KsErr, Scppi0Ks0 + 0.02}}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0 - 0.005, Scppi0Ks0 - Scppi0KsErr}, {Acppi0Ks0 + 0.005, Scppi0Ks0 - Scppi0KsErr}}]},
    {Black, AbsoluteThickness[2.2],
      Line[{{Acppi0Ks0 - 0.005, Scppi0Ks0 + Scppi0KsErr}, {Acppi0Ks0 + 0.005, Scppi0Ks0 + Scppi0KsErr}}]}
  };
  currentDataTextEpilog = {
    Inset[Style[\"Current\\ndata\", 18, Black, FontFamily -> \"Times\"], {Acppi0Ks0 + 0.02, Scppi0Ks0 - 0.15}, {Left, Center}]
  };
  colBlue = RGBColor[0, 0, 1];
  colGreen = RGBColor[0, 1, 0];
  colOrange = RGBColor[1, 0.5, 0];
  colGray = GrayLevel[0.45];
  ClearAll[bandPoly];
  bandPoly[centerPts_, sigmaBand_] := Module[{pairs, up, dn},
    pairs = Select[
      Transpose[{centerPts, sigmaBand[\"Upper\"], sigmaBand[\"Lower\"]}],
      NumericQ[#[[1, 2]]] && NumericQ[#[[2, 2]]] && NumericQ[#[[3, 2]]] &&
        #[[1, 2]] =!= Indeterminate && #[[2, 2]] =!= Indeterminate && #[[3, 2]] =!= Indeterminate &
    ];
    up = #[[2]] & /@ pairs;
    dn = #[[3]] & /@ pairs;
    Polygon[Join[up, Reverse[dn]]]
  ];
  imgSingle = 700;
  imgRow = 1300;
  leftPlot = ListLinePlot[
    {pts1, pts2},
    PlotStyle -> {
      Directive[colGreen, Thick, Dashing[{0.015, 0.015}]],
      Directive[colBlue, Thick, Dashing[{0.015, 0.015}]]
    },
    Evaluate @ axisLikeTarget,
    Prolog -> {
      Directive[colGreen, Opacity[0.15]], bandPoly[pts1, sig1],
      Directive[colBlue, Opacity[0.15]], bandPoly[pts2, sig2]
    },
    Epilog -> Join[sumRuleBandEpilog, sumRuleTextEpilog, currentDataEpilog, currentDataTextEpilog],
    ImageSize -> imgSingle
  ];
  rightPlot = ListLinePlot[
    {pts3, pts4},
    PlotStyle -> {
      Directive[colGray, Thick, Dashing[{0.015, 0.015}]],
      Directive[colOrange, Thick, Dashing[{0.015, 0.015}]]
    },
    Evaluate @ axisLikeTarget,
    Prolog -> {
      Directive[colGray, Opacity[0.15]], bandPoly[pts3, sig3],
      Directive[colOrange, Opacity[0.15]], bandPoly[pts4, sig4]
    },
    Epilog -> Join[sumRuleBandEpilog, sumRuleTextEpilog, currentDataEpilog, currentDataTextEpilog],
    ImageSize -> imgSingle
  ];
  GraphicsRow[{leftPlot, rightPlot}, Spacings -> 0.35, ImageSize -> imgRow, Background -> White]
];";

experimentalInputsCode = "centralInputs = <|
  \"BrB0pi0K0\" -> BrB0pi0K00,
  \"BrB0pimKp\" -> BrB0pimKp0,
  \"BrBppippi0\" -> BrBppippi00,
  \"RTC\" -> RTC0,
  \"Vus\" -> Vus0,
  \"Vud\" -> Vud0,
  \"gamma\" -> gamma0,
  \"phid\" -> phid0,
  \"AcpPimKp\" -> AcppimKp0,
  \"q\" -> q0,
  \"phi\" -> phi0
|>;

paramSpecs = {
  {\"BrB0pi0K0\", centralInputs[\"BrB0pi0K0\"], BrB0pi0K0Err},
  {\"BrB0pimKp\", centralInputs[\"BrB0pimKp\"], BrB0pimKpErr},
  {\"BrBppippi0\", centralInputs[\"BrBppippi0\"], BrBppippi0Err},
  {\"RTC\", centralInputs[\"RTC\"], RTCErr},
  {\"Vus\", centralInputs[\"Vus\"], VusErr},
  {\"Vud\", centralInputs[\"Vud\"], VudErr},
  {\"gamma\", centralInputs[\"gamma\"], gammaErr},
  {\"phid\", centralInputs[\"phid\"], phidErr},
  {\"AcpPimKp\", centralInputs[\"AcpPimKp\"], AcppimKpErr},
  {\"q\", centralInputs[\"q\"], qErr}
};

aRange = {-0.4, 0.30};
aStep = 0.001;
plotRange = {aRange, {-1.0, 1.0}};";

plotConstructionCode = "compute = Function[in, ComputeScpAcpCenter[in, \"ARange\" -> aRange, \"AStep\" -> aStep]];
scalarExtractor = Function[res, <|\"NPoints\" -> Length[res[\"GridA\"]]|>];
curveExtractor = Function[res, res[\"Curves\"]];

centerRes = compute[centralInputs];
uncRes = ComputeUncertaintyFD[compute, scalarExtractor, curveExtractor, centralInputs, paramSpecs];
fig = BuildScpAcpFigure[centerRes, uncRes, centralInputs, \"PlotRange\" -> plotRange];";

exportCode = "outRoot = FileNameJoin[{\"..\", \"results\", \"Scp_Acp\"}];
If[!DirectoryQ[outRoot],
  CreateDirectory[outRoot, CreateIntermediateDirectories -> True]
];

runTag = DateString[{\"Year\", \"Month\", \"Day\", \"_\", \"Hour\", \"Minute\", \"Second\"}];
outDir = FileNameJoin[{outRoot, \"run_\" <> runTag}];
If[!DirectoryQ[outDir],
  CreateDirectory[outDir, CreateIntermediateDirectories -> True]
];

exportBase = FileNameJoin[{outDir, \"Scp_Acp\"}];
exportFigurePath = exportBase <> \".png\";
Export[exportFigurePath, fig, ImageResolution -> 300];
exportFigurePath";

SetDirectory[FileNameJoin[{Directory[], "notebooks"}]];
projectRoot = FileNameJoin[{Directory[], ".."}];
ToExpression[loadInputsCode];
ToExpression[definitionsCode];
ToExpression[experimentalInputsCode];
ToExpression[plotConstructionCode];
Print["plot-built"];

figBoxes = ToBoxes[fig, StandardForm];
Print["plot-boxed"];

exportFigurePath = ToExpression[exportCode];
exportPathBoxes = ToBoxes[exportFigurePath, StandardForm];
Print["figure-exported"];

cells = {
  Cell["Scp_Acp", "Title"],
  Cell["Initialization", "Section"],
  codeCell[initializationCode],
  Cell["Load inputs", "Section"],
  codeCell[loadInputsCode],
  Cell["Definitions for S_CP and A_CP", "Section"],
  codeCell[definitionsCode],
  Cell["Experimental inputs", "Section"],
  codeCell[experimentalInputsCode],
  Cell["Plot construction", "Section"],
  codeCell[plotConstructionCode],
  Cell["Generate final figure", "Section"],
  codeCell["fig"],
  Cell[BoxData[figBoxes], "Output", GeneratedCell -> True, CellAutoOverwrite -> True],
  Cell["Export figure", "Section"],
  codeCell[exportCode],
  Cell[BoxData[exportPathBoxes], "Output", GeneratedCell -> True, CellAutoOverwrite -> True]
};

nb = Notebook[
  cells,
  WindowSize -> {1440, 1080},
  WindowMargins -> {{Automatic, 80}, {Automatic, 40}},
  StyleDefinitions -> "Default.nb"
];

notebookPath = FileNameJoin[{Directory[], "Scp_Acp.nb"}];
Put[nb, notebookPath];

Print[notebookPath];
Print[exportFigurePath];
