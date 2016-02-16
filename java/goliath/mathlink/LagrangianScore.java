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
            if (f % 4 == 0 && f > 0) {
                szPoly += "},{";
            }
            szPoly += poly[f] + (((f - 3) % 4 == 0) ? "" : ",");
        }
        szPoly += "}}";
        System.out.println("Poly received: " + szPoly);
        
        return (exeCmdDoubleArray("scoreAndGetCoefficients[" + szPoly + ", data[[1]], data[[2]],"+df+"]"));

    }

    private static void TestScore() throws MathLinkException {

        //Make sure InitFunctions has been called to load the data
        double[] a = exeCmdDoubleArray("scoreAndGetCoefficients[{{1, 1, 1, 1}, {2, 2, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}}, data[[1]], data[[2]]]");
        for (double i : a) System.out.println(i);
    }

    public static void InitFunctions(String mathPath, String packagePath, String dataName, double deltat, int df) throws MathLinkException {
        OpenLink(mathPath);
        exeCmd("Get[\"" + packagePath + "LagrangeSolver.m\"]");
        exeCmd("data = prepareData[\"" + packagePath + "/" + dataName + "\","+ deltat + "," + df +"];");
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
