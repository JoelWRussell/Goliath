/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author User
 */

public class Data implements Serializable {
    //almost a struct. Other than gather info about the data ; df, deltat, bSim 
    //it also expands the data into a byte[] and can be inflated and deflated via
    //serialize - to be sent along the wire.
    public static void main(String[] Args){
        try {
            Data data = new Data("mma_double.csv");

        } catch (IOException ex) {
            cout("cant find data file?");
        }
    }
    public Data(String szCsv) throws IOException{
     SetBytes(szCsv);
  
        
    }

    
    public byte[] GetBytes(){return rawData;}
    
    public void SetBytes(String szCsv) throws IOException{
        rawData = null; //to be reentrant
        Path path = Paths.get(szCsv);
        rawData = Files.readAllBytes(path);
    }
    
 /////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////////////////

    private byte[] rawData;
}
