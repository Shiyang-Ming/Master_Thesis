(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  uncertainty_MC.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Implements reusable Monte Carlo uncertainty utilities
    for sampled scans and curve-band construction.

  Content:
    - Gaussian random sampling helpers
    - Generic parameter sampling from specifications
    - Monte Carlo sample collection
    - Curve splitting and local normal-vector helpers
    - Quantile-based band polygon construction

  Design Principle:
    - Keeps the numerical Monte Carlo logic close to the
      notebook implementation
    - Stays independent of analysis-specific observables
    - Returns simple data structures for scripts to plot

===========================================================
*)


(*=========================================================*)
(*------------- Monte Carlo sampling helpers --------------*)
(*=========================================================*)
ClearAll[DrawNormalPositive];
DrawNormalPositive[mu_, sigma_] :=
 Module[{x = RandomVariate[NormalDistribution[mu, sigma]]},
  While[x <= 0, x = RandomVariate[NormalDistribution[mu, sigma]]];
  x
 ];

ClearAll[DrawNormalAny];
DrawNormalAny[mu_, sigma_] := RandomVariate[NormalDistribution[mu, sigma]];

ClearAll[SampleParametersMC];
SampleParametersMC[paramSpecs_List] :=
 Association @ Table[
   Module[{name, mu, sigma, mode},
    name = paramSpecs[[i, 1]];
    mu = paramSpecs[[i, 2]];
    sigma = paramSpecs[[i, 3]];
    mode = If[Length[paramSpecs[[i]]] >= 4, paramSpecs[[i, 4]], "Any"];

    name -> Switch[
      mode,
      "Positive", DrawNormalPositive[mu, sigma],
      "Any", DrawNormalAny[mu, sigma],
      _, DrawNormalAny[mu, sigma]
    ]
   ],
   {i, Length[paramSpecs]}
 ];

ClearAll[GenerateMCSamples];
GenerateMCSamples[
  compute_,
  sampler_,
  nSamples_Integer?Positive,
  resultExtractor_ : Identity
] :=
 DeleteCases[
  Table[
   Module[{res},
    res = compute[sampler[]];
    If[res === $Failed, Nothing, resultExtractor[res]]
   ],
   {nSamples}
  ],
  Nothing
 ];


(*=========================================================*)
(*---------------- Curve geometry helpers -----------------*)
(*=========================================================*)
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

ClearAll[CurvePiecesFromIndices];
CurvePiecesFromIndices[curve_List, jump_ : 40] :=
 Module[{validIdx, pieces = {}, current = {}, dphi},
  validIdx = Flatten @ Position[curve, {_?NumericQ, _?NumericQ}];
  If[validIdx === {}, Return[{}]];

  current = {First[validIdx]};
  Do[
   dphi = Abs[curve[[validIdx[[i]], 1]] - curve[[validIdx[[i - 1]], 1]]];
   If[dphi > jump,
    AppendTo[pieces, current];
    current = {validIdx[[i]]},
    AppendTo[current, validIdx[[i]]]
   ],
   {i, 2, Length[validIdx]}
  ];

  AppendTo[pieces, current];
  pieces
 ];

ClearAll[Tangent2D];
Tangent2D[data_, k_] :=
 Module[{n = Length[data], v},
  If[n < 2, Return[{1., 0.}]];

  Which[
   k <= 1, v = data[[Min[2, n]]] - data[[1]],
   k >= n, v = data[[n]] - data[[Max[1, n - 1]]],
   True, v = data[[k + 1]] - data[[k - 1]]
  ];

  v = Chop @ Re @ N[v, 20];
  If[
   !VectorQ[v, NumericQ] || Length[v] != 2 || Norm[v] == 0,
   {1., 0.},
   v/Norm[v]
  ]
 ];

ClearAll[Normal2D];
Normal2D[data_, k_] :=
 Module[{t},
  t = Tangent2D[data, k];
  Chop @ Re @ N[{-t[[2]], t[[1]]}, 20]
 ];


(*=========================================================*)
(*------------- Monte Carlo curve-band helpers ------------*)
(*=========================================================*)
ClearAll[MakeBandPolygonMC];
Options[MakeBandPolygonMC] = {
  "LowerQuantile" -> 0.16,
  "UpperQuantile" -> 0.84,
  "MinPointsPerSlice" -> 10
};

MakeBandPolygonMC[
  centralCurve_List,
  sampleCurves_List,
  idxPiece_List,
  OptionsPattern[]
] :=
 Module[
  {
   upperPts = {}, lowerPts = {}, pos, k, center, subcurve, nn, pts, proj,
   qlo, qhi, qloVal, qhiVal, minPts
  },

  qloVal = OptionValue["LowerQuantile"];
  qhiVal = OptionValue["UpperQuantile"];
  minPts = OptionValue["MinPointsPerSlice"];

  subcurve = centralCurve[[idxPiece]];
  subcurve = Cases[subcurve, {_?NumericQ, _?NumericQ}];
  subcurve = Chop @ Re @ N[subcurve, 20];
  If[Length[subcurve] < 3, Return[Nothing]];

  Do[
   k = idxPiece[[pos]];
   center = centralCurve[[k]];
   If[!ValidPoint2DQ[center], Continue[]];
   center = Chop @ Re @ N[center, 20];

   nn = Normal2D[subcurve, pos];
   If[!VectorQ[nn, NumericQ] || Length[nn] != 2, Continue[]];

   pts = Cases[sampleCurves[[All, k]], {_?NumericQ, _?NumericQ}];
   pts = Chop @ Re @ N[pts, 20];
   pts = Select[pts, ValidPoint2DQ];
   If[Length[pts] < minPts, Continue[]];

   proj = ((# - center) . nn) & /@ pts;
   proj = Chop @ Re @ N[proj, 20];
   proj = Select[proj, NumericQ];
   If[Length[proj] < minPts, Continue[]];

   qlo = Quantile[proj, qloVal];
   qhi = Quantile[proj, qhiVal];

   AppendTo[upperPts, center + qhi nn];
   AppendTo[lowerPts, center + qlo nn];
   ,
   {pos, Length[idxPiece]}
  ];

  upperPts = Select[upperPts, ValidPoint2DQ];
  lowerPts = Select[lowerPts, ValidPoint2DQ];

  If[
   Length[upperPts] < 3 || Length[lowerPts] < 3,
   Nothing,
   Polygon @ Join[upperPts, Reverse[lowerPts]]
  ]
 ];

ClearAll[BuildBandPolygonsMC];
Options[BuildBandPolygonsMC] = Join[
  {"Jump" -> 40},
  Options[MakeBandPolygonMC]
];

BuildBandPolygonsMC[
  centralCurve_List,
  sampleCurves_List,
  OptionsPattern[]
] :=
 Module[{pieces, polys},
  pieces = CurvePiecesFromIndices[centralCurve, OptionValue["Jump"]];
  polys = DeleteCases[
    MakeBandPolygonMC[
      centralCurve,
      sampleCurves,
      #,
      "LowerQuantile" -> OptionValue["LowerQuantile"],
      "UpperQuantile" -> OptionValue["UpperQuantile"],
      "MinPointsPerSlice" -> OptionValue["MinPointsPerSlice"]
    ] & /@ pieces,
    Nothing
  ];
  polys
 ];
