(* ::Package:: *)

(*============== Script for x and \[CapitalDelta] computation ==============*)

Quiet @ Check[

  (* -------- paths -------- *)
  scriptDir  = DirectoryName[$InputFileName];
  projectDir = ParentDirectory[scriptDir];

  Get[FileNameJoin[{projectDir, "src", "init.wl"}]];

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
    "Acp00" -> Acppi0pi00,
    "Scp00" -> Scppi0pi00
  |>;

  paramSpecs = {
    {"A",     Acppippim0, AcppippimErr},
    {"S",     Scppippim0, ScppippimErr},
    {"gamma", gamma0,     gammaErr},
    {"phid",  phid0,      phidErr},
    {"Rpm",   Rpm0,       RpmErr},
    {"R00",   R000,       R00Err}
  };

  deltaStepDeg = 1;

  (* -------- compute -------- *)
  centerRes = ComputeXDeltaCenter[inputs, "DeltaStepDeg" -> deltaStepDeg];

  compute = Function[in, ComputeXDeltaCenter[in, "DeltaStepDeg" -> deltaStepDeg]];
  scalarExtractor = Function[res, res["Center"]];
  curveExtractor  = Function[res, res["Curves"]];

  uncRes = ComputeUncertaintyFD[
    compute, scalarExtractor, curveExtractor,
    inputs, paramSpecs
  ];

  (* -------- numbers -------- *)
  x0 = centerRes["Center"]["x"];
  dx = uncRes["Scalars"]["Sigma"]["x"];

  Delta0deg = centerRes["Center"]["Delta"]/Degree;
  dDeltadeg = uncRes["Scalars"]["Sigma"]["Delta"]/Degree;

  Print["x = ", N[x0, 8], " \[PlusMinus] ", N[dx, 4]];
  Print["\[CapitalDelta] = ", N[Delta0deg, 8], "\[Degree] \[PlusMinus] ", N[dDeltadeg, 4], "\[Degree]"];

  (* -------- plot -------- *)
  fig = PlotXDelta[centerRes, uncRes, "Title" -> "x ~ \[CapitalDelta]"];

  outDir = FileNameJoin[{projectDir, "results", "B_pipi_x_Delta"}];
  If[!DirectoryQ[outDir],
    CreateDirectory[outDir, CreateIntermediateDirectories -> True]
  ];

  base = FileNameJoin[{outDir, "xDelta"}];

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
      "DeltaStepDeg" -> deltaStepDeg
    |>
  ];

  (* 2) .txt : human-readable *)
  Export[
    base <> ".txt",
    StringRiffle[
      {
        "x = " <> ToString@N[x0, 12] <>
          " +/- " <> ToString@N[dx, 8],

        "Delta(deg) = " <> ToString@N[Delta0deg, 12] <>
          " +/- " <> ToString@N[dDeltadeg, 8],

        "DeltaStepDeg = " <> ToString[deltaStepDeg]
      },
      "\n"
    ],
    "String"
  ];

  Print["[OK] Saved: results/B_pipi_x_Delta/xDelta.png"];
  Print["[OK] Saved: results/B_pipi_x_Delta/xDelta.mx"];
  Print["[OK] Saved: results/B_pipi_x_Delta/xDelta.txt"];

  Exit[0];

,
  Print["[ERROR] Script failed."];
  Exit[1];
];
