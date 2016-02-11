(* ::Package:: *)

BeginPackage["LagrangeSolver`"]



PrepareData::usage = "PrepareData[sz, deltat]"

Begin["`Private`"]

makeControlTrajectory[vars_, tTrain_, tStep_] := Flatten[{
     symbolAppend[#, ""] -> Table[0.1 t, {t, 0, tTrain, tStep}],
     symbolAppend[#, "d"] -> Table[0.1 , {t, 0, tTrain, tStep}],
     symbolAppend[#, "dd"] -> Table[0, {t, 0, tTrain, tStep}]} & /@ 
   vars]

symbolAppend[symbol_, postfix_] := 
  Symbol[SymbolName[symbol] <> postfix];
dotVar[v_] := symbolAppend[v, "dot"];


PrepareData[sz_, deltat_] := Module[
  {experimentalData, times, expth1, expth2, expw1, expw2, expa1, expa2},
  
  experimentalData = 
   Transpose[ 
    Partition[ ReadList[ sz, Number, RecordSeparators -> {","}], 2]];
times =  Range[0, (Length[experimentalData[[1]]] - 1)]*deltat;
expth1 = Transpose[{times, experimentalData[[1]] + 1*^-15}];
expth2 = Transpose[{times, experimentalData[[2]] + 1*^-15}];
expw1 = (Interpolation[expth1])'[times];
expw2 = (Interpolation[expth2])'[times];
expa1 = (Interpolation[expth1])''[times];
expa2 = (Interpolation[expth2])''[times];
expth1 = experimentalData[[1]];
expth2 = experimentalData[[2]];
myExperimentalData = 
  Flatten[{ {symbolAppend[th1,""] -> expth1, symbolAppend[th1, "d"] -> expw1, 
     symbolAppend[th1, "dd"] -> expa1},
    {symbolAppend[th2,""] -> expth2, symbolAppend[th2, "d"] -> expw2, 
     symbolAppend[th2, "dd"] -> expa2}}];
myControlData = 
  makeControlTrajectory[{th1, th2}, times[[Length[times]]], deltat];
{myExperimentalData, myControlData}
]




End[]
EndPackage[]
