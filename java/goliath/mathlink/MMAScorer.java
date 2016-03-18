/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import com.wolfram.jlink.*;
/**
 *
 * @author User
 */
public class MMAScorer {
    public boolean DEBUG = true;
    public MMAScorer(String szMath, String szMathPath, String szPackagePath){
        this.szMath = szMath;
        this.szMathPath = szMathPath;
        this.szPackagePath = szPackagePath;
    }
    public void InitKernel() throws MathLinkException{
        if (DEBUG) cout("MMScorer:InitKernel");
        OpenLink(szMathPath);
    }
    public void DestroyKernel(){
        if (DEBUG) cout("MMScorer:DestroyKernel");  
        CloseLink();
    }
    public void PrepareData(String dataName, DataFormat format) throws MathLinkException {
        if (DEBUG) cout("MMScorer:PrepareData");
        this.bSim = format.bSim;
        this.df = format.df;
        this.deltat = format.deltat;

        if (bSim){
            InitFunctionsSim(szMathPath, szPackagePath, dataName,deltat, df );
        }
        else{
            InitFunctions(szMathPath, szPackagePath, dataName,deltat, df );
        }
    }
    public Population ProcessPopulation(Population population) throws MathLinkException {
        if (DEBUG) cout("MMScorer:ProcessPopulation");
        for (int f=0; f<population.polys.size(); f++){
            Poly pl = population.polys.get(f);

                double[] res = GetScore(pl.powers, df);
               

                for (int ff=0; ff<res.length; ff++){
                if (ff==0){
                    pl.score = res[0];
                }
                else{
                    pl.coefficients[ff-1] = res[ff];
                }
            }
 

            population.polys.set(f, pl);          
            

        }
        cout(population.toString());
        return population;
     } 
    
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    private  KernelLink ml = null;
    private double[] GetScore(int[] poly, int df) throws MathLinkException {
        //poly is an array that represents a single poly
        //it returns a score and a series of coefficients

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
        return sc;

    }
    private void InitFunctions(String mathPath, String packagePath, String dataName, double deltat, int df) throws MathLinkException {
        exeCmd("Get[\"" + packagePath + "" + "LagrangeSolverCompact.m\"]");
        exeCmd("data = prepareData[\"" + packagePath + "" + dataName + "\","+ deltat + "," + df +"]");
    }  
    private  void InitFunctionsSim(String mathPath, String packagePath, String dataName, double deltat, int df) throws MathLinkException {
        exeCmd("Get[\"" + packagePath + "" + "LagrangeSolverCompact.m\"]");
        exeCmd("data = prepareDataSim[\"" + packagePath + "" + dataName + "\","+ deltat + "," + df +"]");
        //processes the data with vels and accels
     
    }
    ///////utility functions//////////////////////////////////////////////////////////
    private  void OpenLink(String mathPath) throws MathLinkException {
        ml = MathLinkFactory.createKernelLink(szMath + "\"" + mathPath + "\"");
        ml.discardAnswer();
    }
    private  void CloseLink() {
        if (ml != null) ml.close();
        ml = null;
    }
    private  void exeCmd(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.discardAnswer();
    }
    private  double exeCmdDouble(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getDouble());
    }
    private  int exeCmdInt(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getInteger());
    }
    private  Boolean exeCmdBoolean(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getBoolean());
    }
    private  String exeCmdString(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getString());
    }
    private  double[] exeCmdDoubleArray(String cmd) throws MathLinkException {
        ml.evaluate(cmd);
        ml.waitForAnswer();
        return (ml.getDoubleArray1());
    }
    private String szMath, szMathPath, szPackagePath;
    private double deltat = 0.1;
    private int df = 1;
    private boolean bSim = false;
}
