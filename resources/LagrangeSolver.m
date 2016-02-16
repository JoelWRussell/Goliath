(* ::Package:: *)

BeginPackage["LagrangeSolver`"]



PrepareData::usage = "PrepareData[sz, deltat, df]"

Begin["`Private`"]

makeControlTrajectory[vars_, tTrain_, tStep_] := Flatten[{
     symbolAppend[#, ""] -> Table[0.1 t, {t, 0, tTrain, tStep}],
     symbolAppend[#, "d"] -> Table[0.1 , {t, 0, tTrain, tStep}],
     symbolAppend[#, "dd"] -> Table[0, {t, 0, tTrain, tStep}]} & /@ 
   vars]

symbolAppend[symbol_, postfix_] := 
  Symbol[SymbolName[symbol] <> postfix];
dotVar[v_] := symbolAppend[v, "dot"];


PrepareData[sz_, deltat_, df_] := Module[
  {experimentalData, times, vars, expTh, expW, expA, expD, 
myExperimentalData, myControlData},
  
experimentalData = 
 Transpose[ 
  Partition[ ReadList[ sz, Number, RecordSeparators -> {","}], df]];
times =  Range[0, (Length[experimentalData[[1]]] - 1)]*deltat;
vars = symbolAppend[th, ToString[#]] & /@ Range[df];
expTh = Transpose[{times, # + 1*^-20}] & /@ experimentalData;
expW = (Interpolation[#])'[times] & /@ expTh;
expA = (Interpolation[#])''[times] & /@ expTh;
expD = Transpose[{ vars, 
   experimentalData, expW, expA}];
 myExperimentalData = 
 Flatten[{#[[1]] -> #[[2]], symbolAppend[#[[1]], "d"] -> #[[3]], 
     symbolAppend[#[[1]], "dd"] -> #[[4]]} & /@ expD];
myControlData = 
 makeControlTrajectory[vars, times[[Length[times]]], deltat];

{myExperimentalData, myControlData}

]




End[]
EndPackage[]