/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import java.io.Serializable;

/**
 *
 * @author User
 */
//static variables can not be serialised
public class Poly implements Serializable {
    //This is a polynomial it packages information about the coefficients,
    //structure of the polynomial and the score.
    public static void main(String[] Args){
        int[] pp = {1,1,2,2};
        Poly p = new Poly(pp,1);
        System.out.println(p.toString());
        
        
    }
    public Poly(int[] p, int df){
        powers = p;
        coefficients = new double[p.length/(2*df)];
     
        score = 10000;
    }
    @Override public String toString(){
        String sz="";
        for (int f : powers){
            sz = sz + f+" ";
        }
       

        sz += "   "+score+"\n";
        return sz;
    }
    public int[] powers;//[2 0 0 0] [0 2 0 0]  as a flat list
    public double[] coefficients;//c0, c1, c2 etc for each monomial
    public double score;//score of the poly returned from  the remote clients
    
    
    
    
}
