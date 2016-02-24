(* ::Package:: *)

makeControlTrajectory[vars_, tTrain_, tStep_] := Flatten[{
     symbolAppend[#, ""] -> Table[0.1 t, {t, 0, tTrain, tStep}],
     symbolAppend[#, "d"] -> Table[0.1 , {t, 0, tTrain, tStep}],
     symbolAppend[#, "dd"] -> Table[0, {t, 0, tTrain, tStep}]} & /@
   vars]

symbolAppend[symbol_, postfix_] :=
  Symbol[SymbolName[symbol] <> postfix];

dotVar[v_] := symbolAppend[v, "dot"];

prepareData[sz_, deltat_, df_] := Module[
  {experimentalData, times, vars, expTh, expW, expA, expD, myExperimentalData, myControlData},
experimentalData = Transpose[Partition[ ReadList[ sz, Number, RecordSeparators -> {","}], df]];
times =  Range[0, (Length[experimentalData[[1]]] - 1)]*deltat;
vars = symbolAppend[th, ToString[#]] & /@ Range[df];
expTh = Transpose[{times, # + 1*^-20}] & /@ experimentalData;
expW = (Interpolation[#])'[times] & /@ expTh;
expA = (Interpolation[#])''[times] & /@ expTh;
expD = Transpose[{ vars, experimentalData, expW, expA}];
myExperimentalData = Flatten[{#[[1]] -> #[[2]], symbolAppend[#[[1]], "d"] -> #[[3]], symbolAppend[#[[1]], "dd"] -> #[[4]]} & /@ expD];
myControlData = makeControlTrajectory[vars, times[[Length[times]]], deltat];
{myExperimentalData, myControlData}
]

prepareDataSim[sz_, deltat_, df_] := Module[
  {expD, times, vars, expTh, expW, expA,  myExperimentalData, myControlData},
expD = Transpose[Partition[ ReadList[ sz, Number, RecordSeparators -> {","}], 3 df]];
vars = symbolAppend[th, ToString[#]] & /@ Range[df];
times =  Range[0, (Length[expD[[1]]] - 1)]*deltat;
myExperimentalData = Flatten[{vars[[#]] -> expD[[#]], symbolAppend[vars[[#]], "d"]->expD[[#+df]], symbolAppend[vars[[#]], "dd"]->expD[[#+2 df]] } & /@ Range[df]];
myControlData = makeControlTrajectory[vars, times[[Length[times]]], deltat];
{myExperimentalData, myControlData}
]

generateTransformations[vars_]:=Module[{dotVars,toTimeDependent,toDerivatives},
dotVars=dotVar/@vars;
toTimeDependent=#->#[t]&/@Join[vars,dotVars];
toDerivatives=dotVar[#][t]->#'[t]&/@vars;
{toTimeDependent,toDerivatives}
]

generateELScore[vars_]:=Module[{toTimeDependent,toDerivatives},
{toTimeDependent,toDerivatives}=generateTransformations[vars];
Function[l,
Plus@@((D[D[l/.toTimeDependent,dotVar[#][t]]/.toDerivatives,t]-D[l/.toTimeDependent,#[t]]/.toDerivatives)^2&/@vars)
]]

generateNScore[vars_]:=Module[{toTimeDependent,toDerivatives},
{toTimeDependent,toDerivatives}=generateTransformations[vars];
Function[l,
Plus@@((D[D[l/.toTimeDependent,dotVar[#][t]]/.toDerivatives,t])^2+(D[l/.toTimeDependent,#[t]]/.toDerivatives)^2&/@vars)
]]

targetUnity[x_]:=Log[10^-25+x]^2+1

pathSum[vars_,score_,p_]:=Module[{toData},
toData=Flatten[{#[t]->#,#'[t]->symbolAppend[#,"d"],#''[t]->symbolAppend[#,"dd"]}&/@vars];
Plus@@((score/.toData)/.p)
]

generateScore[vars_,l_,trajectory_,controlTrj_]:=Log[(10^-10+Expand[pathSum[vars,generateELScore[vars][l],trajectory]])/(10^-10+Expand[pathSum[vars,generateELScore[vars][l],controlTrj]]) targetUnity[Expand[pathSum[vars,generateNScore[vars][l],controlTrj]]]targetUnity[Expand[pathSum[vars,generateELScore[vars][l],controlTrj]]]]

nmSolve[scoreExpression_,coeffs_,model_]:=Module[{sol},
err={};step=0;eval=0;dimensions=Length[coeffs];bestModel=0;
sol=FindMinimum[scoreExpression,
coeffs,
Method->"QuasiNewton",
AccuracyGoal->10,
PrecisionGoal->10,
MaxIterations->50 10^5,
StepMonitor:>(step++;If[Mod[step,500]==0,err=Append[err,scoreExpression];bestModel=model]),
EvaluationMonitor:>(eval++)];
{sol[[1]],sol[[2]]}
]

polynomialiser2[vars_, indv_]:= Module[{nVars, fullVarList, monoList, coeffs},
nVars = 2*Length[vars];
fullVarList=Flatten[{vars, dotVar[#]&/@vars}];
monoList=Times@@(fullVarList^#)& /@ indv;
coeffs = Symbol["c"<>ToString[#]]&/@Range[Length[monoList]];
{coeffs, monoList.coeffs}
];

scoreAndGetCoefficients[polyData_, expData_, controlData_, df_]:=Module[{vars, myPolynomial, scoreFn, result},
vars = symbolAppend[th, ToString[#]] & /@ Range[df];
myPolynomial = polynomialiser2[vars, polyData];
scoreFn = generateScore[vars,myPolynomial[[2]], expData, controlData];
result = nmSolve[scoreFn, myPolynomial[[1]], myPolynomial[[2]]];
Flatten[{result[[1]], myPolynomial[[1]]/.result[[2]]}]
]


