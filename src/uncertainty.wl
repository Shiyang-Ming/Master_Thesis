(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  uncertainty.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Implements uncertainty propagation strategies for
    derived observables.

  Content:
    - One-by-one \[PlusMinus]1\[Sigma] parameter variation
    - Quadrature combination
    - Curve band construction
    - (Optional future) Monte Carlo sampling methods

  Design Principle:
    - Completely separated from physics formulae
    - Accepts computation functions as arguments
    - Returns structured uncertainty information

===========================================================
*)


(*=========================================================*)
(*------ Uncertainty propagation (finite difference) ------*)
(*=========================================================*)
ClearAll[PropagateScalarsFD];
Options[PropagateScalarsFD] = {
  "StepScale" -> 1,                    (* \:7528 sig * StepScale \:505a\:6b65\:957f *)
  "Combine" -> (Sqrt[Total[#^2]] &)    (* \:9ed8\:8ba4 quadrature combine *)
};

PropagateScalarsFD[
  compute_,                             (* compute[inputs_] := result *)
  scalarExtractor_,                     (* scalarExtractor[result_] := <|...|> *)
  inputs_Association,
  paramSpecs_List,                      (* {{name, central, sigma}...} \:6216\:4f60\:53ea\:7528\:5230 {name,_,sigma} *)
  opts : OptionsPattern[]
] :=
 Module[
  {centerRes, center, names, stepScale, vars, combine, sigmas},

  stepScale = OptionValue["StepScale"];
  combine = OptionValue["Combine"];

  centerRes = compute[inputs];
  center = scalarExtractor[centerRes];
  names = paramSpecs[[All, 1]];

  (* \:6bcf\:4e2a\:53c2\:6570\:6270\:52a8\:540e\:7684\:6807\:91cf\:7ed3\:679c *)
  vars = Association @ Table[
     Module[{name = paramSpecs[[i, 1]], sig = paramSpecs[[i, 3]]*stepScale, upIn, dnIn, upC, dnC},
      upIn = inputs; dnIn = inputs;
      upIn[name] = upIn[name] + sig;
      dnIn[name] = dnIn[name] - sig;

      upC = scalarExtractor @ compute[upIn];
      dnC = scalarExtractor @ compute[dnIn];

      name -> <|"Up" -> upC, "Down" -> dnC, "HalfDiff" -> (upC - dnC)/2|>
     ],
     {i, Length[paramSpecs]}
   ];

  (* \:5bf9\:6bcf\:4e2a\:6807\:91cf key\:ff1a\:628a\:6240\:6709\:53c2\:6570\:8d21\:732e\:5408\:6210 *)
  sigmas = Association @ Map[
     Function[key,
      key -> combine@Table[ vars[n]["HalfDiff"][key], {n, names} ]
     ],
     Keys[center]
   ];

  <|
    "CenterResult" -> centerRes,
    "Center" -> center,
    "Sigma" -> sigmas,
    "Variations" -> vars
  |>
];

ClearAll[BuildCurveBandsFD];
Options[BuildCurveBandsFD] = {
  "StepScale" -> 1,
  "Combine" -> (Sqrt[Total[#^2]] &)
};

BuildCurveBandsFD[
  compute_,
  curveExtractor_,
  inputs_Association,
  paramSpecs_List,
  opts : OptionsPattern[]
] :=
 Module[
  {centerRes, curves, curveNames, xGrid, yCenter, names, stepScale, combine,
   dY, sigY, bands},

  stepScale = OptionValue["StepScale"];
  combine = OptionValue["Combine"];

  centerRes = compute[inputs];
  curves = curveExtractor[centerRes];
  curveNames = Keys[curves];
  names = paramSpecs[[All, 1]];

  (* \:9ed8\:8ba4\:4f7f\:7528\:7b2c\:4e00\:6761\:66f2\:7ebf\:7684 x \:4f5c\:4e3a\:516c\:5171\:7f51\:683c *)
  xGrid = curves[curveNames[[1]]][[All, 1]];
  yCenter = AssociationMap[curves[#][[All, 2]] &, curveNames];

  (* dY[name][curve] = (up - down)/2 \:7684 y-array *)
  dY = Association@Table[
     Module[{pname = paramSpecs[[i, 1]], sig = paramSpecs[[i, 3]]*stepScale, upIn, dnIn, upCur, dnCur},
      upIn = inputs; dnIn = inputs;
      upIn[pname] = upIn[pname] + sig;
      dnIn[pname] = dnIn[pname] - sig;

      upCur = curveExtractor @ compute[upIn];
      dnCur = curveExtractor @ compute[dnIn];

      pname -> AssociationMap[(upCur[#][[All, 2]] - dnCur[#][[All, 2]])/2 &, curveNames]
     ],
     {i, Length[paramSpecs]}
   ];

  (* \:5148\:53d6\:6bcf\:6761\:66f2\:7ebf\:3001\:6bcf\:4e2a\:53c2\:6570\:7684 half-diff y-array\:ff0c\:7136\:540e\:9010\:70b9 quadrature \:5408\:6210 *)
  sigY = AssociationMap[
  Function[cname,
    Sqrt @ Total @ (Table[dY[p][cname], {p, names}]^2)
  ],
  curveNames
  ];

  bands = AssociationMap[
    Function[cname,
      <|
        "Upper" -> Transpose[{xGrid, yCenter[cname] + sigY[cname]}],
        "Lower" -> Transpose[{xGrid, yCenter[cname] - sigY[cname]}]
      |>
    ],
    curveNames
  ];

  <|
    "CenterResult" -> centerRes,
    "Curves" -> curves,
    "Bands" -> bands,
    "SigmaY" -> sigY,
    "dY" -> dY
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
] :=
 Module[{scalar, band},
  scalar = PropagateScalarsFD[compute, scalarExtractor, inputs, paramSpecs,
    "StepScale" -> OptionValue["StepScale"], "Combine" -> OptionValue["Combine"]
  ];
  band = BuildCurveBandsFD[compute, curveExtractor, inputs, paramSpecs,
    "StepScale" -> OptionValue["StepScale"], "Combine" -> OptionValue["Combine"]
  ];

  <|
    "CenterResult" -> scalar["CenterResult"],  (* \:6216 band["CenterResult"]\:ff0c\:5e94\:4e00\:81f4 *)
    "Scalars" -> <|"Center" -> scalar["Center"], "Sigma" -> scalar["Sigma"], "Variations" -> scalar["Variations"]|>,
    "Curves" -> <|"Center" -> band["Curves"], "Bands" -> band["Bands"], "SigmaY" -> band["SigmaY"], "dY" -> band["dY"]|>
  |>
];

