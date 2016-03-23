/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author User
 */
//Population --- is just a struct to hold data
//### use as a proper class later on
//it is serializable which means it can be sent down the wire
//
//the population knows the data it is supposed to work on
public class Population implements Serializable {
    public String szData;
    public int df;
    public double deltat;
    public ArrayList<Poly> polys = new ArrayList<Poly>();
    @Override public String toString(){
        String sz = "Population##\n";
        for (Poly p : polys){
           if (p.score<0) sz += p.toString();
        }
    
        return sz;
    }
    
}
