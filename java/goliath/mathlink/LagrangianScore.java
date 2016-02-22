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

	private static double bestScore = 1000000.0;
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
        //System.out.println(szPoly);
        double[] sc = exeCmdDoubleArray("scoreAndGetCoefficients[" + szPoly + ", data[[1]], data[[2]],"+df+"]");

	bestScore = (sc[0] < bestScore)? sc[0] : bestScore;
        System.out.println("" + bestScore);

        return sc;

    }

 
    public static void InitFunctions(String mathPath, String packagePath, String dataName, double deltat, int df) throws MathLinkException {
        OpenLink(mathPath);
        exeCmd("Get[\"" + packagePath + "" + "LagrangeSolverCompact.m\"]");
        exeCmd("data = prepareData[\"" + packagePath + "" + dataName + "\","+ deltat + "," + df +"]");
    }

    public static void Shutdown() {
        
        CloseLink();
    }


    public static void main(String[] args) throws MathLinkException {
        ;
    }
}