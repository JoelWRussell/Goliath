/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import static com.lagrangianmining.Utility.getServer;
import static com.lagrangianmining.Utility.website;
import java.io.IOException;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

/**
 *
 * @author User
 */
public class GorillaClient {
    public boolean DEBUG = true;
    public static void main(String[] Args){
        try {
  
            GorillaClient gor = new GorillaClient();
            gor.InitMathKernel();//start the kernel
            gor.NewData("mma_double.csv");//copy the data across to all hosts
            gor.PrepareData(false, 2, 0.1);//load worksheet in mma and the correct data
            
            //for 2 polynomials
            int[] polyLength = {4, 4};// ie  (2 0 0 2) would be 4
            int[] poly = {1,1,1,1,   2,2,2,2};
            gor.NewZeitgeist(2, polyLength, poly, 2, 0.1);
             double[] r1 = gor.GetScores();
             double[] r2 = gor.GetCoefficients();
           
            for (double ddd : r1){
                System.out.print(ddd + " ");
            }
            System.out.println("");
            for (double ddd : r2){
                System.out.print(ddd + " ");
            }

          gor.DestroyMathKernel();
           gor.DisconnectServer();
//            
        } catch (Exception ex) {
             System.out.println("GorillaClient: error");
        }
    }
    
    public GorillaClient() throws Exception{ 
        this("localHost");
    }
    public GorillaClient(String serverName) throws Exception{
        cout("GorillaClient: connects between a clj file and the LagrangeServer");
        cout("#################################################################");
        cout("The clj document with the GA code uses this GorillaClient - so calls methods \non it to send the data and zeitgeist to the LagrangeServer (which needs to be running).");
        FindServerFromWeb(serverName);
        GetServer();
    }
    public void DisconnectServer() throws IOException{
        socket.DisconnectSocket();
    }
    public boolean NewData(String szCSV) throws IOException{
        //send NEW_DATA and a Data
        Data data = new Data(szCSV);
        //if (DEBUG) System.out.print(data.toString());
        socket.WriteString("NEW_DATA");
        socket.WriteObject(data);
       
        String szMsg = socket.ReadString(); // receives an OK message to prevent it collapsing the socket
        if (DEBUG) System.out.println(szMsg);
        if (szMsg.equals("OK")){
            return true;
        }
        else return false;
        
    }
    public boolean NewZeitgeist(int numPolynomials, int[] sizeOfEachPolynomial, int[] polynomials, int df, double deltat) throws IOException, ClassNotFoundException{
        //send NEW_ZEITGEIST and Zeitgeist
        zg = new Zeitgeist();
        zg.SetData(numPolynomials, sizeOfEachPolynomial, polynomials, df, deltat);
        socket.WriteString("NEW_ZEITGEIST");
        socket.WriteObject(zg);
        String szResult = socket.ReadString();
        if (szResult.equals("RESULTS")){
            zg = (Zeitgeist)socket.ReadObject();
            if (DEBUG) System.out.println("GorillaClient: Got results zeitgeist");
           //if (DEBUG) System.out.println(zg.toString());
          
            return true;
           
        }
        
        return false;

    }
    public boolean InitMathKernel() throws IOException{
        socket.WriteString("INIT_MATH_KERNEL");
        String szMsg = socket.ReadString(); 
        if (DEBUG) System.out.println(szMsg);
        if (szMsg.equals("OK")){
            return true;
        }
        else return false;
    }
    public boolean DestroyMathKernel() throws IOException{
        socket.WriteString("DESTROY_MATH_KERNEL");
        String szMsg = socket.ReadString(); 
        if (DEBUG) System.out.println(szMsg);
        if (szMsg.equals("OK")){
            return true;
        }
        else return false;
    }
    public boolean PrepareData(boolean bSim, int df, double deltat) throws IOException{
                //send NEW_DATA and a Data
        DataFormat format = new DataFormat(bSim, df, deltat);
        //if (DEBUG) System.out.print(data.toString());
        socket.WriteString("PREPARE_DATA");
        socket.WriteObject(format);
       
        String szMsg = socket.ReadString(); // receives an OK message to prevent it collapsing the socket
        if (DEBUG) System.out.println(szMsg);
        if (szMsg.equals("OK")){
            return true;
        }
        else return false;
    }
    public double[] GetScores(){
       
        return zg.GetScoresList();
        
        
    }
    public double[] GetCoefficients(){
         return zg.GetCoefficientsList();
        
        
    }

 ///////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////
    private void  GetServer(){
        try {
           socket = new ObjectSocket(new Socket(ip, port));
           socket.InitSocket();
           if (DEBUG) cout("GorillaClient: found server");
       } catch (IOException ex) {
           if (DEBUG) cout("GorillaClient: unable to find server");
           GetServer();
            try {
                Thread.sleep(1000);
            } catch (InterruptedException ex1) {
                Logger.getLogger(GorillaClient.class.getName()).log(Level.SEVERE, null, ex1);
            }
       }
      
    }
    private void FindServerFromWeb(String serverName) throws IOException, Exception {
        //gets the ip:port form the www given a short name from the website
            XMLHelper xml = new XMLHelper();
            
            CloseableHttpClient httpclient = HttpClients.createDefault();
            HttpGet httpget = new HttpGet(website+getServer+"?"+"server=" + serverName);
            CloseableHttpResponse response = httpclient.execute(httpget); 
               
                HttpEntity entity = response.getEntity();
                if (entity != null){
                    String sz = EntityUtils.toString(entity);
                    xml.setDocument(sz);
                    /////###set permissions on Network Address Translation
                    ////##NAT, there is no TCP hole punching
                    
                    ip = xml.EvaluateXPath("/server/ip");
                    port = Integer.parseInt(xml.EvaluateXPath("/server/port1")); 


                }
                
            if (DEBUG) cout("Connecting to LagrangeServer at: ip="+ip+" port="+port);
    }

    private String ip;
    private int port; 
    private ObjectSocket socket; 
    private Zeitgeist zg = null; 
    
}
