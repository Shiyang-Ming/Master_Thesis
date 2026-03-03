(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  plotting.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Provides visualization utilities for computed observables.

  Content:
    - Plotting routines
    - Uncertainty bands visualization
    - Formatting, legends, and export utilities

  Design Principle:
    - No physics calculations inside this file
    - Consumes structured outputs from computations.wl
    - Purely responsible for graphical representation

===========================================================
*)


(*===============================================*)
(*------------- Plot x\[Dash]Delta figure -------------*)
(*===============================================*)

(*ClearAll[PlotXDelta];
Options[PlotXDelta] = {
  "PlotRange" -> {{0, 360}, {0, 3.5}},
  "LegendPos" -> Scaled[{0.20, 0.86}],
  "Title" -> "x ~ \[CapitalDelta]"
};

PlotXDelta[centerRes_Association, unc_Association : Missing["NoUncertainty"], OptionsPattern[]] :=
 Module[
  {
   curves, bands, hasBands,
   centerRpm, centerR00, acp1, acp2,
   upperRpm, lowerRpm, upperR00, lowerR00
  },

  curves = centerRes["Curves"];
  centerRpm = curves["Rpm"];
  centerR00 = curves["R00"];
  acp1 = curves["Acp1"];
  acp2 = curves["Acp2"];

  hasBands = AssociationQ[unc] && KeyExistsQ[unc, "CurveBands"];
  If[hasBands,
   bands = unc["CurveBands"];
   upperRpm = bands["RpmUpper"]; lowerRpm = bands["RpmLower"];
   upperR00 = bands["R00Upper"]; lowerR00 = bands["R00Lower"];
  ];

  Show[
   {
    (* center dashed curves *)
    ListLinePlot[{centerRpm, centerR00},
     PlotStyle -> {{Gray, Dashed, Thick}, {Red, Dashed, Thick}}
    ],

    (* error bands if provided *)
    If[hasBands,
     {
      ListLinePlot[{upperRpm, lowerRpm}, PlotStyle -> {{Gray, Thick}, {Gray, Thick}}],
      ListLinePlot[{upperR00, lowerR00}, PlotStyle -> {{Red, Thick}, {Red, Thick}}]
     },
     {}
    ],

    (* blue dashed Acp00 branches if provided *)
    If[ListQ[acp1] && acp1 =!= Missing["NotProvided"], ListLinePlot[{acp1}, PlotStyle -> {{Blue, Dashed, Thick}}], {}],
    If[ListQ[acp2] && acp2 =!= Missing["NotProvided"], ListLinePlot[{acp2}, PlotStyle -> {{Blue, Dashed, Thick}}], {}]
   },

   Frame -> True,
   FrameLabel -> {"\[CapitalDelta] (\[Degree])", "x"},
   PlotLabel -> Style[OptionValue["Title"], 14, Black],
   ImageSize -> Large,
   LabelStyle -> Directive[FontSize -> 13],
   FrameStyle -> Black,
   TicksStyle -> Black,
   PlotRange -> OptionValue["PlotRange"],
   Background -> White,

   Epilog -> {
     Inset[
      LineLegend[
       {GrayLevel[.5], Red, Directive[Blue, Dashed]},
       {
        Style["R^{\[Pi]\[Pi]}_{+-}", 16, FontFamily -> "Times", Black],
        Style["R^{\[Pi]\[Pi]}_{00}", 16, FontFamily -> "Times", Black],
        Style["A^{\[Pi]^0 \[Pi]^0}_{CP}", 16, FontFamily -> "Times", Black]
       },
       LegendMarkers -> None,
       LegendLayout -> "Column",
       LabelStyle -> {FontFamily -> "Times", FontSize -> 16},
       Spacings -> .7,
       LegendFunction -> (Framed[#, FrameStyle -> None, Background -> None] &)
      ],
      OptionValue["LegendPos"]
     ]
   }
  ]
 ];*)
 
 ClearAll[PlotXDelta];
Options[PlotXDelta] = {
  "PlotRange" -> {{0, 360}, {0, 3.5}},
  "LegendPos" -> Scaled[{0.20, 0.86}],
  "Title" -> "x ~ \[CapitalDelta]"
};

PlotXDelta[centerRes_Association, unc_Association, OptionsPattern[]] :=
 Module[
  {
   curves,
   centerRpm, centerR00, acp1, acp2,
   bands, upperRpm, lowerRpm, upperR00, lowerR00
  },

  curves = centerRes["Curves"];
  centerRpm = curves["Rpm"];
  centerR00 = curves["R00"];
  acp1 = Lookup[curves, "Acp1", Missing["NotProvided"]];
  acp2 = Lookup[curves, "Acp2", Missing["NotProvided"]];

  (* NEW format only: unc["Curves"]["Bands"][curveName]["Upper"/"Lower"] *)
  bands = unc["Curves"]["Bands"];

  upperRpm = bands["Rpm"]["Upper"];  lowerRpm = bands["Rpm"]["Lower"];
  upperR00 = bands["R00"]["Upper"];  lowerR00 = bands["R00"]["Lower"];

  Show[
   {
    (* center dashed curves *)
    ListLinePlot[{centerRpm, centerR00},
     PlotStyle -> {{Gray, Dashed, Thick}, {Red, Dashed, Thick}}
    ],

    (* error bands *)
    {
     ListLinePlot[{upperRpm, lowerRpm}, PlotStyle -> {{Gray, Thick}, {Gray, Thick}}],
     ListLinePlot[{upperR00, lowerR00}, PlotStyle -> {{Red, Thick}, {Red, Thick}}]
    },

    (* blue dashed Acp00 branches if provided *)
    If[ListQ[acp1] && acp1 =!= Missing["NotProvided"],
      ListLinePlot[{acp1}, PlotStyle -> {{Blue, Dashed, Thick}}],
      {}
    ],
    If[ListQ[acp2] && acp2 =!= Missing["NotProvided"],
      ListLinePlot[{acp2}, PlotStyle -> {{Blue, Dashed, Thick}}],
      {}
    ]
   },

   Frame -> True,
   FrameLabel -> {"\[CapitalDelta] (\[Degree])", "x"},
   PlotLabel -> Style[OptionValue["Title"], 14, Black],
   ImageSize -> Large,
   LabelStyle -> Directive[FontSize -> 13],
   FrameStyle -> Black,
   TicksStyle -> Black,
   PlotRange -> OptionValue["PlotRange"],
   Background -> White,

   Epilog -> {
     Inset[
      LineLegend[
       {GrayLevel[.5], Red, Directive[Blue, Dashed]},
       {
        Style["R^{\[Pi]\[Pi]}_{+-}", 16, FontFamily -> "Times", Black],
        Style["R^{\[Pi]\[Pi]}_{00}", 16, FontFamily -> "Times", Black],
        Style["A^{\[Pi]^0 \[Pi]^0}_{CP}", 16, FontFamily -> "Times", Black]
       },
       LegendMarkers -> None,
       LegendLayout -> "Column",
       LabelStyle -> {FontFamily -> "Times", FontSize -> 16},
       Spacings -> .7,
       LegendFunction -> (Framed[#, FrameStyle -> None, Background -> None] &)
      ],
      OptionValue["LegendPos"]
     ]
   }
  ]
];

