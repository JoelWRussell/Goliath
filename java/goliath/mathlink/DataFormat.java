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
public class DataFormat implements Serializable {
    public DataFormat(boolean bSim, int df, double deltat ){
        this.bSim = bSim;
        this.df = df;
        this.deltat = deltat;
    }
    public int df;
    public boolean bSim;
    public double deltat;
}
