/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author User
 */
public class Zeitgeist implements Serializable {
    public String szData;
    //so that we know what this zg is about
    public static void main(String[] Args){
        int[] p1 = {4,6,2,2,2,2,2,2};
        int[] p2 = {1,1,1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8};
        Zeitgeist zg = new Zeitgeist();
        zg.SetData(8, p1, p2, 1, 0.1);
        zg.SetNumberClients(5);
        System.out.println(zg.toString());
        double[] d = zg.GetCoefficientsList();
        for (double dd: d){
            System.out.print(dd + " ");
        }
        try {
            zg.SaveToFile("TestZeitgeist.txt");
        } catch (IOException ex) {
            Logger.getLogger(Zeitgeist.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    public Zeitgeist(){
  
         polyList = new ArrayList<Poly>(); 
         popList = new ArrayList<Population>();
    }
    @Override public String toString(){
        String sz = "POLY LIST*****\n";
        for (Poly p : polyList){
            sz += p.toString();
        }
        sz += "*********\n";
        sz += "POP LIST*****\n";
        for (Population p : popList){
            sz += p.toString();
        }
        sz += "********\n";
        return sz;
    }
    //SetData---reads the incomming data into an ArrayList polyList
    public void SetData(int numPolys, int[] sizePolyList, int[] pList, int df, double deltat){
        this.df = df;
        this.deltat = deltat;
        int c = 0;
        //create a new polyList each time that this is run.
        polyList = null;
        polyList = new ArrayList<Poly>();
        for (int f=0; f<numPolys; f++){
            int[] ccc = new int[sizePolyList[f]];
            for (int c1 = 0; c1<sizePolyList[f]; c1++){
                ccc[c1] = pList[c1 + c];
            }
   
            polyList.add(new Poly(ccc, df));
            c += sizePolyList[f];
        }
    }
    public void SetPopulation(int i, Population pop){
        popList.set(i, pop);
    }
    public void FinalizeZeitgeist(){
        polyList = new ArrayList<Poly>();
        int total = 0;
        for (int f=0; f<numClients; f++){
            total += popList.get(f).polys.size();
        } 
        int[] cnt = new int[numClients];
        for (int f=0; f<numClients; f++){
            cnt[f] = 0;
        }
        for (int f=0; f<total; f++){
            polyList.add(popList.get(f%numClients).polys.get(cnt[f%numClients]));
            cnt[f%numClients]++;
        }
    }
    public void SaveToFile(String szFile) throws FileNotFoundException, IOException{
        File fout = new File(szFile);
	FileOutputStream fos = new FileOutputStream(fout);
	BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
        for (int f=0; f<polyList.size(); f++){
            bw.write(polyList.get(f).toString());

        }
        bw.close();
    }
    public int GetNumPolys(){
        return polyList.size();
    }
    //SetNumberClients---partitions the polyList into an ArrayList of Populations which is each Serializable
    public void SetNumberClients(int num){
        numClients = num;
        popList = new ArrayList<Population>();
        //set up the arraylist of populations
        for (int f=0; f<numClients; f++){
            Population p = new Population();
            p.deltat = deltat;
            p.df = df;
            popList.add(p);
            
        }
        //sort all of the polys into populations-order of the zeitgeist does
        //not matter
       
        for (int f=0; f<polyList.size(); f++){
           
           popList.get(f%numClients).polys.add(polyList.get(f));
          
        }
        

      
    }
    //GetScoresList, GetCoefficientsList gives the results of the clients
    public double[] GetScoresList(){
        double[] res = new double[polyList.size()];
        for (int f=0; f<polyList.size(); f++){
           res[f] = polyList.get(f).score;
        }
        return res;
    }
    public double[] GetCoefficientsList(){
     
        int arraySize = 0;
        for (int f=0; f<polyList.size(); f++){
            arraySize += polyList.get(f).coefficients.length;
        }
        double[] res = new double[arraySize];
        int cnt = 0;
        for (int f=0; f<polyList.size(); f++){
            for (int ff=0; ff<polyList.get(f).coefficients.length; ff++){
                res[cnt] = polyList.get(f).coefficients[ff];
                cnt++;
            }
        }
        return res;
    }
    public Population GetPopulation(int f){
        return popList.get(f);
    }
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    private int df;
    private double deltat;
    private int numClients;
    private ArrayList<Poly> polyList; 
    private ArrayList< Population> popList;
   
}
