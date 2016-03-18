/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 *
 * @author User
 */
public class Utility {
    static public String website = "http://www.lagrangianmining.com";
    static public String getServer = "/cgi-bin/Get_Server.pl";
    static public String registerServer = "/cgi-bin/Register_Server.pl";

    static public String szMath = "-linkmode launch -linkname ";
    static public String szMathPath = "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe";
    static public String szPackagePath = "C:/Users/User/Documents/NetBeansProjects/com.lagrangianmining/";//might need to look at this
    static public String szExpData = "experimental.txt";
//I am so used to c++ that I need cin anc cout
      public static String cin() throws IOException{
        return  (new BufferedReader(new InputStreamReader(System.in))).readLine();    
        }
       public static void cout(String sz){
        System.out.println(sz);
    }  
}
