(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  q_phi_scan.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Implements reusable functions for the q-\[Phi] scan
    analysis and its Monte Carlo uncertainty bands.

  Content:
    - Amplitude-building helpers
    - q-\[Phi] branch construction
    - Rc constraint curves
    - Plot assembly helpers

  Design Principle:
    - Keeps the numerical logic close to the working
      test implementation
    - Reuses generic Monte Carlo helpers from
      uncertainty_MC.wl
    - Leaves script-level execution and exports to
      scripts/q_phi_scan.wl

===========================================================
*)


(*=========================================================*)
(*----------------- Low-level q-\[Phi] helpers -----------------*)
(*=========================================================*)
ClearAll[QPhiPhaseSpaceMomentum];
QPhiPhaseSpaceMomentum[mB_, m1_, m2_] :=
  Sqrt[(mB^2 - (m1 + m2)^2) (mB^2 - (m1 - m2)^2)]/(2 mB);

ClearAll[QPhiAmplitudeMagnitude];
QPhiAmplitudeMagnitude[br_, tau_, mB_, m1_, m2_] :=
  Sqrt[N[br/(tau*QPhiPhaseSpaceMomentum[mB, m1, m2]), 20]];

ClearAll[QPhiTriangleAngle];
QPhiTriangleAngle[a_, b_, x_, sign_ : 1] :=
 Module[{cosa},
  cosa = N[(x^2 + a^2 - b^2)/(2 x a), 20];
  cosa = Chop @ Re @ N[cosa, 20];
  If[NumericQ[cosa] && -1 <= cosa <= 1,
    sign ArcCos[cosa],
    Indeterminate
  ]
 ];


(*=========================================================*)
(*---------------- q-\[Phi] branch construction ----------------*)
(*=========================================================*)
ClearAll[ComputeQPhiBranchData];
ComputeQPhiBranchData[pars_Association, nGrid_ : 350, mKPlus_ : Automatic] :=
 Module[
  {
   gammaVal, RTCVal, epsVal, AcpPiPlusK0, AcpPi0KPlus, BrPiPlusK0,
   BrPi0KPlus, BrPiPlusPi0, Aplus0Mag, Abarplus0Mag, A0plusMag,
   Abar0plusMag, ApiPiMag, TpCpMag, Aplus0Model, Abarplus0Model,
   phiCVal, xMin1, xMax1, xMin2, xMax2, xMin, xMax, tGrid, xGrid,
   signPairs, qPhiPointLocal, deltaPhi32Local, branchData, mKPlusVal
  },

  mKPlusVal = If[mKPlus === Automatic, MKp0, mKPlus];

  gammaVal = pars["gamma"];
  RTCVal = pars["RTC"];
  epsVal = pars["Vus"]/pars["Vud"];
  AcpPiPlusK0 = pars["AcpPiPlusK0"];
  AcpPi0KPlus = pars["AcpPi0KPlus"];
  BrPiPlusK0 = pars["BrPiPlusK0"];
  BrPi0KPlus = pars["BrPi0KPlus"];
  BrPiPlusPi0 = pars["BrPiPlusPi0"];

  Aplus0Mag =
    QPhiAmplitudeMagnitude[BrPiPlusK0 (1 - AcpPiPlusK0), TBp0, MBp0, Mpip0, MK00];
  Abarplus0Mag =
    QPhiAmplitudeMagnitude[BrPiPlusK0 (1 + AcpPiPlusK0), TBp0, MBp0, Mpip0, MK00];

  A0plusMag =
    Sqrt[2] QPhiAmplitudeMagnitude[BrPi0KPlus (1 - AcpPi0KPlus), TBp0, MBp0, Mpi00, mKPlusVal];
  Abar0plusMag =
    Sqrt[2] QPhiAmplitudeMagnitude[BrPi0KPlus (1 + AcpPi0KPlus), TBp0, MBp0, Mpi00, mKPlusVal];

  ApiPiMag =
    Sqrt[2] QPhiAmplitudeMagnitude[BrPiPlusPi0, TBp0, MBp0, Mpip0, Mpi00];

  TpCpMag = RTCVal*epsVal*ApiPiMag;
  TpCpMag = Chop @ Re @ N[TpCpMag, 20];
  If[!NumericQ[TpCpMag] || TpCpMag <= 0, Return[$Failed]];

  Aplus0Model =
    1 + pars["rhoc"] Exp[I pars["thetac"]] Exp[I gammaVal];
  Abarplus0Model =
    1 + pars["rhoc"] Exp[I pars["thetac"]] Exp[-I gammaVal];

  phiCVal = Arg[Abarplus0Model Conjugate[Aplus0Model]];
  phiCVal = Chop @ Re @ N[phiCVal, 20];

  xMin1 = Abs[Aplus0Mag - A0plusMag];
  xMax1 = Aplus0Mag + A0plusMag;
  xMin2 = Abs[Abarplus0Mag - Abar0plusMag];
  xMax2 = Abarplus0Mag + Abar0plusMag;

  xMin = Max[xMin1, xMin2];
  xMax = Min[xMax1, xMax2];
  xMin = Chop @ Re @ N[xMin, 20];
  xMax = Chop @ Re @ N[xMax, 20];

  If[!(NumericQ[xMin] && NumericQ[xMax] && xMin < xMax), Return[$Failed]];

  tGrid = Subdivide[0., 1., nGrid];
  xGrid = xMin + tGrid (xMax - xMin);

  signPairs = Tuples[{-1, 1}, 2];

  deltaPhi32Local[x_, {s1_, s2_}] :=
   Module[{alpha, beta, val},
    alpha = QPhiTriangleAngle[Aplus0Mag, A0plusMag, x, s1];
    beta = phiCVal + QPhiTriangleAngle[Abarplus0Mag, Abar0plusMag, x, s2];
    val = If[NumericQ[alpha] && NumericQ[beta], alpha - beta, Indeterminate];
    If[NumericQ[val], Chop @ Re @ N[val, 20], Indeterminate]
   ];

  qPhiPointLocal[x_, dphi_, eta_] :=
   Module[{Nval, cval, sval, qval, phival, pt},
    Nval = x/TpCpMag;
    cval = eta Nval Cos[dphi/2];
    sval = eta Nval Sin[dphi/2];
    qval = Sqrt[(Cos[gammaVal] - cval)^2 + (Sin[gammaVal] - sval)^2];
    phival = Arg[(Cos[gammaVal] - cval) + I (Sin[gammaVal] - sval)];
    pt = Chop @ Re @ N[{phival/Degree, qval}, 20];
    If[ValidPoint2DQ[pt], pt, Missing["InvalidPoint"]]
   ];

  branchData = Table[
    Module[{dphiList, curvePlus, curveMinus},
      dphiList = Table[
        deltaPhi32Local[xGrid[[k]], signPairs[[i]]],
        {k, Length[xGrid]}
      ];

      curvePlus = Table[
        If[NumericQ[dphiList[[k]]],
          qPhiPointLocal[xGrid[[k]], dphiList[[k]], 1],
          Missing["InvalidPoint"]
        ],
        {k, Length[xGrid]}
      ];

      curveMinus = Table[
        If[NumericQ[dphiList[[k]]],
          qPhiPointLocal[xGrid[[k]], dphiList[[k]], -1],
          Missing["InvalidPoint"]
        ],
        {k, Length[xGrid]}
      ];

      {curvePlus, curveMinus}
    ],
    {i, Length[signPairs]}
  ];

  <|"BranchData" -> branchData, "TGrid" -> tGrid|>
 ];


(*=========================================================*)
(*---------------- Rc-constraint utilities ----------------*)
(*=========================================================*)
ClearAll[ComputeRcObservable];
ComputeRcObservable[brPi0KPlus_, brPi0KPlusErr_, brPiPlusK0_, brPiPlusK0Err_] :=
 Module[{Rc0, RcErr},
  Rc0 = 2 brPi0KPlus/brPiPlusK0;
  RcErr = Rc0 Sqrt[
    (brPi0KPlusErr/brPi0KPlus)^2 +
    (brPiPlusK0Err/brPiPlusK0)^2
  ];
  <|"Value" -> Rc0, "Error" -> RcErr|>
 ];

ClearAll[ComputeQPhiRcCurves];
ComputeQPhiRcCurves[phiDegGrid_, RcVal_, gammaVal_, rhocVal_, thetacVal_, rcVal_, deltacVal_] :=
 Module[{phi, A, B, C, disc, qPlus, qMinus},
  Table[
    phi = phiDeg Degree;

    A = rcVal^2;

    B = 2 rcVal (
      Cos[deltacVal] Cos[phi] -
      (rcVal - rhocVal Cos[thetacVal - deltacVal]) Cos[gammaVal - phi]
    );

    C = (1 + 2 rhocVal Cos[thetacVal] Cos[gammaVal] + rhocVal^2) (1 - RcVal)
        - 2 rhocVal rcVal Cos[thetacVal - deltacVal]
        - 2 rcVal Cos[deltacVal] Cos[gammaVal]
        + rcVal^2;

    disc = Chop @ Re @ N[B^2 - 4 A C, 20];

    If[NumericQ[disc] && disc >= 0,
      qPlus = Chop @ Re @ N[(-B + Sqrt[disc])/(2 A), 20];
      qMinus = Chop @ Re @ N[(-B - Sqrt[disc])/(2 A), 20];
      {
        If[NumericQ[qPlus] && qPlus > 0, {phiDeg, qPlus}, Missing["Invalid"]],
        If[NumericQ[qMinus] && qMinus > 0, {phiDeg, qMinus}, Missing["Invalid"]]
      },
      {Missing["Invalid"], Missing["Invalid"]}
    ],
    {phiDeg, phiDegGrid}
  ]
 ];


(*=========================================================*)
(*------------------- Plotting utilities ------------------*)
(*=========================================================*)
ClearAll[QPhiDefaultBranchColors];
QPhiDefaultBranchColors[] := {
  RGBColor[0.90, 0.40, 0.10],
  RGBColor[0.15, 0.65, 0.20],
  RGBColor[0.02, 0.12, 0.72],
  GrayLevel[0.40]
};

ClearAll[QPhiDefaultKeepBranches];
QPhiDefaultKeepBranches[] := {
  {1, 1}, {1, 2},
  {2, 1},
  {3, 1},
  {4, 1}, {4, 2}
};

ClearAll[MakeQPhiBandGraphics];
Options[MakeQPhiBandGraphics] = {
  "KeepBranches" -> Automatic,
  "Jump" -> 40,
  "Opacity" -> 0.22
};

MakeQPhiBandGraphics[branchData_, sampleResults_, branchColors_List, OptionsPattern[]] :=
 Module[{keepBranches, jump, alpha},
  keepBranches = Replace[OptionValue["KeepBranches"], Automatic -> QPhiDefaultKeepBranches[]];
  jump = OptionValue["Jump"];
  alpha = OptionValue["Opacity"];

  Graphics[
    Flatten @ Table[
      Module[{i, s, centralCurve, sampleCurves, polys},
        i = br[[1]];
        s = br[[2]];
        centralCurve = branchData[[i, s]];
        sampleCurves = sampleResults[[All, i, s]];
        polys = BuildBandPolygonsMC[centralCurve, sampleCurves, "Jump" -> jump];
        If[polys === {}, {}, {Directive[branchColors[[i]], Opacity[alpha], EdgeForm[None]], polys}]
      ],
      {br, keepBranches}
    ]
  ]
 ];

ClearAll[MakeQPhiBranchPlots];
Options[MakeQPhiBranchPlots] = {
  "KeepBranches" -> Automatic,
  "Jump" -> 40
};

MakeQPhiBranchPlots[branchData_, branchColors_List, OptionsPattern[]] :=
 Module[{keepBranches, jump},
  keepBranches = Replace[OptionValue["KeepBranches"], Automatic -> QPhiDefaultKeepBranches[]];
  jump = OptionValue["Jump"];

  DeleteCases[
    Table[
      Module[{i, s, pieces},
        i = br[[1]];
        s = br[[2]];
        pieces = SplitCurveByJump[branchData[[i, s]], jump];
        If[pieces === {},
          Nothing,
          ListLinePlot[
            pieces,
            PlotStyle -> Directive[branchColors[[i]], Thick],
            InterpolationOrder -> 1
          ]
        ]
      ],
      {br, keepBranches}
    ],
    Nothing
  ]
 ];

ClearAll[MakeQPhiRcPlots];
MakeQPhiRcPlots[phiGridRc_, Rc0_, RcErr_, gammaVal_, rhocVal_, thetacVal_, rcVal_, deltacVal_, rcColor_] :=
 Module[
  {
   rcLowData, rcHighData, rcLowPlus, rcLowMinus, rcHighPlus, rcHighMinus,
   rcPlots, rcLabelGraphic
  },

  rcLowData = ComputeQPhiRcCurves[
    phiGridRc, Rc0 - RcErr, gammaVal, rhocVal, thetacVal, rcVal, deltacVal
  ];
  rcHighData = ComputeQPhiRcCurves[
    phiGridRc, Rc0 + RcErr, gammaVal, rhocVal, thetacVal, rcVal, deltacVal
  ];

  rcLowPlus = Cases[rcLowData[[All, 1]], {_?NumericQ, _?NumericQ}];
  rcLowMinus = Cases[rcLowData[[All, 2]], {_?NumericQ, _?NumericQ}];
  rcHighPlus = Cases[rcHighData[[All, 1]], {_?NumericQ, _?NumericQ}];
  rcHighMinus = Cases[rcHighData[[All, 2]], {_?NumericQ, _?NumericQ}];

  rcPlots = DeleteCases[
    {
      If[rcLowPlus === {}, Nothing,
        ListLinePlot[
          SplitCurveByJump[rcLowPlus, 8],
          PlotStyle -> Directive[rcColor, Dashed, Thick],
          InterpolationOrder -> 1
        ]
      ],
      If[rcLowMinus === {}, Nothing,
        ListLinePlot[
          SplitCurveByJump[rcLowMinus, 8],
          PlotStyle -> Directive[rcColor, Dashed, Thick],
          InterpolationOrder -> 1
        ]
      ],
      If[rcHighPlus === {}, Nothing,
        ListLinePlot[
          SplitCurveByJump[rcHighPlus, 8],
          PlotStyle -> Directive[rcColor, Dashed, Thick],
          InterpolationOrder -> 1
        ]
      ],
      If[rcHighMinus === {}, Nothing,
        ListLinePlot[
          SplitCurveByJump[rcHighMinus, 8],
          PlotStyle -> Directive[rcColor, Dashed, Thick],
          InterpolationOrder -> 1
        ]
      ]
    },
    Nothing
  ];

  rcLabelGraphic = Graphics[
    Text[
      Style[Subscript["R", "c"], 18, rcColor, Bold],
      {-62, 2.18}
    ]
  ];

  <|"Plots" -> rcPlots, "LabelGraphic" -> rcLabelGraphic|>
 ];

ClearAll[PlotQPhiScan];
Options[PlotQPhiScan] = {
  "PlotRange" -> {{-180, 180}, {0, 3}},
  "ImageSize" -> 650,
  "AspectRatio" -> 0.72
};

PlotQPhiScan[bandGraphics_, branchPlots_List, rcPlots_List : {}, rcLabelGraphic_ : Graphics[{}], OptionsPattern[]] :=
 Show[
  bandGraphics,
  Sequence @@ branchPlots,
  Sequence @@ rcPlots,
  rcLabelGraphic,
  Frame -> True,
  Axes -> False,
  FrameStyle -> Directive[Black, 14],
  FrameTicksStyle -> Directive[Black, 12],
  FrameLabel -> {
    Style["\[Phi] [Degree]", 16, Black],
    Style["q", 16, Black]
  },
  PlotRange -> OptionValue["PlotRange"],
  PlotRangePadding -> Scaled[0.02],
  ImageSize -> OptionValue["ImageSize"],
  AspectRatio -> OptionValue["AspectRatio"],
  Background -> White
 ];
