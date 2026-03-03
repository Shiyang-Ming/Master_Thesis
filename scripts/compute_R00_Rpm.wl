(* ::Package:: *)

(*============== Script for R00 and Rpm computation ==============*)
scriptDir = DirectoryName[$InputFileName];
projectDir = ParentDirectory[scriptDir];

Get[FileNameJoin[{projectDir, "src", "init.wl"}]];

outDir = FileNameJoin[{projectDir, "results", "B_pipi_R00_Rpm"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

res = ComputeR00Rpm[];

Print["R00 = ", res["R00"]["Value"], " \[PlusMinus] ", res["R00"]["Error"]];
Print["Rpm = ", res["Rpm"]["Value"], " \[PlusMinus] ", res["Rpm"]["Error"]];

Export[FileNameJoin[{outDir, "R00.mx"}], {res["R00"]["Value"], res["R00"]["Error"]}];
Export[FileNameJoin[{outDir, "Rpm.mx"}], {res["Rpm"]["Value"], res["Rpm"]["Error"]}];

