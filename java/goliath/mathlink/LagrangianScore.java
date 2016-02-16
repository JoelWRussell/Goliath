/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package goliath.mathlink;

import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

/**
 * @author User
 */
public class LagrangianScore {

    private static String szMath = "-linkmode launch -linkname ";

    ///Put commands between OpenLink and CloseLink; there is no exception
    //handling; the exception is just bounced out of the program
    //exeCmd discards the return value
    //exeCmdDouble expects to return a double
        /* examples....
   OpenLink();
   exeCmd("x=5");
   exeCmd("x=10");
   exeCmd("f[x_]:=x^2");
   System.out.println(exeCmdDouble("f[2]"));
   System.out.println(exeCmdDoubleArray("{1.0, 2.0}")[1]);
   CloseLink();
                */
    //recall that \ is the escape character in Java so " \" "  is the string which is just "
    //also \ escapes itself so "\\" is the string \
    //also copying and pasting from Mathematica to Java seems to automatically convert the
    //escape characters etc...!
    private static KernelLink ml = null;

    private static void OpenLink(String mathPath) throws MathLinkException {
        ml = MathLinkFactory.createKernelLink(szMath + "\"" + mathPath + "\"");
        ml.discardAnswer();
    }

    private static void CloseLink() {
        if (ml != null) ml.close();
    }

    private static void exeCmd(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.discardAnswer();
    }

    private static double exeCmdDouble(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getDouble());
    }

    private static int exeCmdInt(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getInteger());
    }

    private static Boolean exeCmdBoolean(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getBoolean());
    }

    public static String exeCmdString(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getString());
    }

    private static double[] exeCmdDoubleArray(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getDoubleArray1());
    }


    // input a list of integers and it will score to the data etc
    public static synchronized double[] GetScore(int[] poly, int df) throws MathLinkException {
        //clojure will send an array of doubles
        //package into groups of 4
        String szPoly = "{{";
        for (int f = 0; f < poly.length; f++) {
            if (f % (2*df) == 0 && f > 0) {
                szPoly += "},{";
            }
            szPoly += poly[f] + (((f - (2*df-1)) % (2*df) == 0) ? "" : ",");
        }
        szPoly += "}}";
        System.out.println("Poly received: " + szPoly);
        
        return (exeCmdDoubleArray("ScoreAndGetCoefficients[" + szPoly + ", data[[1]], data[[2]],"+df+"]"));

    }

    private static void TestScore() throws MathLinkException {

        //Make sure InitFunctions has been called to load the data
        double[] a = exeCmdDoubleArray("ScoreAndGetCoefficients[{{1, 1, 1, 1}, {2, 2, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}}, data[[1]], data[[2]]]");
        for (double i : a) System.out.println(i);
    }

    public static void InitFunctions(String mathPath, String packagePath, String dataName, double deltat, int df) throws MathLinkException {
        OpenLink(mathPath);
        exeCmd("Get[\"" + packagePath + "LagrangeSolver.m\"]");
        exeCmd("data = PrepareData[\"" + packagePath + "/" + dataName + "\","+ deltat + "," + df +"];");
        exeCmd("symbolAppend[symbol_, postfix_] := Symbol[SymbolName[symbol] <> postfix]; "
                + "dotVar[v_] := symbolAppend[v, \"dot\"];");
        exeCmd("generateTransformations[vars_] := Module[{dotVars, toTimeDependent, toDerivatives}," +
                "dotVars = dotVar /@ vars;" +
                "toTimeDependent = # -> #[t] & /@ Join[vars, dotVars];" +
                "toDerivatives = dotVar[#][t] -> #'[t] & /@ vars;" +
                "{toTimeDependent, toDerivatives}]");
        exeCmd("polynomialiser2[vars_, indv_] :=" +
                "Module[{nVars, fullVarList, monoList, coeffs}," +
                "nVars = 2*Length[vars];" +
                "fullVarList = Flatten[ {vars, dotVar[#] & /@ vars}];" +
                "monoList = Times @@ (fullVarList^#) & /@ indv;" +
                "coeffs = Symbol[\"c\" <> ToString[#]] & /@ Range[Length[monoList]];" +
                "{coeffs, monoList.coeffs}" +
                "];");
        exeCmd("generateELScore[vars_] := Module[{toTimeDependent, toDerivatives}," +
                "{toTimeDependent, toDerivatives} = generateTransformations[vars];" +
                "Function[l," +
                " Plus @@ ((D[D[l /. toTimeDependent, dotVar[#][t]] /. toDerivatives," +
                "t] - D[l /. toTimeDependent, #[t]] /. toDerivatives)^2 & /@" +
                "vars)" +
                "]]");
        exeCmd("generateNScore[vars_] := Module[{toTimeDependent, toDerivatives}," +
                "{toTimeDependent, toDerivatives} = generateTransformations[vars];" +
                "Function[l," +
                "Plus @@ ((D[D[l /. toTimeDependent, dotVar[#][t]] /. toDerivatives," +
                "t])^2 + (D[l /. toTimeDependent, #[t]] /. " +
                "toDerivatives)^2 & /@ vars)" +
                "]]");
        exeCmd("targetUnity[x_] := Log[10^-25 + x]^2 + 1");
        exeCmd("pathSum[vars_, score_, p_] := Module[{toData}," +
                "toData = " +
                "Flatten[{#[t] -> #, #'[t] -> symbolAppend[#, \"d\"], #''[t] -> " +
                "symbolAppend[#, \"dd\"]} & /@ vars];" +
                "Plus @@ ((score /. toData) /. p)" +
                "]");
        exeCmd("generateScore[vars_, l_, trajectory_, controlTrj_] :=" +
                "Log[((10^-10 +" +
                "Expand[pathSum[vars, generateELScore[vars][l]," +
                "trajectory]])/(10^-10 +" +
                "Expand[pathSum[vars, generateELScore[vars][l]," +
                "controlTrj]])) targetUnity[" +
                "Expand[pathSum[vars, generateNScore[vars][l]," +
                "controlTrj]]] targetUnity[" +
                "Expand[pathSum[vars, generateELScore[vars][l], controlTrj]]]]");
//        exeCmd("nmSolve[scoreExpression_, coeffs_, model_] := Module[{sol}," +
//                "err = {}; step = 0; eval = 0; dimensions = Length[coeffs];" +
//                "bestModel = 0;" +
//                "sol = NMinimize[" +
//                "Join[{scoreExpression}, Thread[-1.0 < # < 1.0 & /@ coeffs]]," +
//                "coeffs," +
//                "AccuracyGoal -> 10," +
//                "PrecisionGoal -> 10," +
//                "Method -> {" +
//                "\"NelderMead\"," +
//                "\"PostProcess\" -> False," +
//                "\"ExpandRatio\" -> 1 + (2/Length[coeffs])," +
//                "\"ContractRatio\" -> (3/4) - (1/(2 Length[coeffs]))," +
//                "\"ShrinkRatio\" -> 1 - (1/Length[coeffs])}," +
//                "MaxIterations -> 50 10^5," +
//                "StepMonitor :> (step++; " +
//                "If[Mod[step, 500] == 0, err = Append[err, scoreExpression];" +
//                "bestModel = model])," +
//                "EvaluationMonitor :> (eval++)];" +
//                "{\"steps\" -> step, \"bestScore\" -> sol[[1]], \"solution\" -> sol[[2]]," +
//                "\"model\" -> model /. sol[[2]]}" +
//                "]");
        exeCmd("nmSolve[scoreExpression_, coeffs_, model_] := Module[{sol},\n" +
                "  err = {}; step = 0; eval = 0; dimensions = Length[coeffs]; \n" +
                "  bestModel = 0;\n" +
                "  sol = FindMinimum[scoreExpression,\n" +
                "    coeffs,\n" +
                "    Method -> \"QuasiNewton\",\n" +
                "    AccuracyGoal -> 10,\n" +
                "    PrecisionGoal -> 10,\n" +
                "    MaxIterations -> 50 10^5,\n" +
                "    StepMonitor :> (step++; \n" +
                "      If[Mod[step, 500] == 0, err = Append[err, scoreExpression]; \n" +
                "       bestModel = model]),\n" +
                "    EvaluationMonitor :> (eval++)];\n" +
                "  {\"steps\" -> step, \"bestScore\" -> sol[[1]], \"solution\" -> sol[[2]], \n" +
                "   \"model\" -> model /. sol[[2]]}\n" +
                "  ]");
        exeCmd("ScoreAndGetCoefficients[polyData_, expData_, controlData_, df_] := " +
                "Module[{vars, myPolynomial, scoreFn, result}," +
                "vars = symbolAppend[th, ToString[#]] & /@ Range[df];" +
                "myPolynomial = polynomialiser2[vars, polyData];" +
                "scoreFn = " +
                "generateScore[vars, myPolynomial[[2]], expData, controlData];" +
                "result = nmSolve[scoreFn, myPolynomial[[1]], myPolynomial[[2]]];" +
                "Flatten[{\"bestScore\" /. result[[2]]," +
                " myPolynomial[[1]] /. \"solution\" /. result[[3]]}]" +
                "]");
        

    }

    public static void Shutdown() {
        CloseLink();
    }


    public static void main(String[] args) throws MathLinkException {
        InitFunctions(args[0], args[1], args[2], Double.valueOf(args[3]), Integer.valueOf(args[4]));
        int[] a = {1, 1, 1, 1, 2, 2, 2, 2};
        double[] res = GetScore(a,2);
        for (double i : res) System.out.println(i);
        Shutdown();
    }
}