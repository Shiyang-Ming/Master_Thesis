(* ::Package:: *)

(* ::Package:: *)
(**)


(*
===========================================================
  init.wl
  Master Thesis Project
-----------------------------------------------------------

  Purpose:
    Central loader for all core modules of the project.
    This file does NOT perform any computation.
    It only loads source files in the correct order.

  Structure:
    inputs_exp.wl        \[Dash] Experimental inputs
    functions.wl         \[Dash] Core physics/helper functions
    computations.wl      \[Dash] Derived quantities and main computations
    uncertainty.wl       \[Dash] Error propagation utilities
    plotting.wl          \[Dash] Plotting utilities

  Usage:
    In scripts:
        Get["../src/init.wl"];

    In notebooks:
        SetDirectory[NotebookDirectory[]];
        Get["../src/init.wl"];

===========================================================
*)


srcDir = DirectoryName[$InputFileName];
SetDirectory[srcDir];

Get["inputs_exp.wl"];
Get["functions.wl"];
Get["computations.wl"];
Get["uncertainty.wl"];
Get["plotting.wl"];
Get["parameters.wl"];

