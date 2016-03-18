/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import static com.lagrangianmining.Utility.getServer;
import static com.lagrangianmining.Utility.szExpData;
import static com.lagrangianmining.Utility.szMath;
import static com.lagrangianmining.Utility.szMathPath;
import static com.lagrangianmining.Utility.szPackagePath;
import static com.lagrangianmining.Utility.website;
import com.wolfram.jlink.MathLinkException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
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
public class WorkerClient {
    public boolean DEBUG = true;
    public static void main (String[] Args){
        try {
            WorkerClient client = new WorkerClient ();
            client.GetServer();
            client.Start();
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
            cout("cant find server");
        }
        
        cout("exit");
        
    }
    public WorkerClient(){
        this("localHost");
    }
    public WorkerClient(String serverName){
        cout("Starting WorkerClient");
        cout("#####################");
        cout("Worker client will init Mathematica kernel and process the population.");
        scorer = new MMAScorer(szMath, szMathPath, szPackagePath);
        try {
            FindServerFromWeb(serverName);
        } catch (Exception ex) {
            cout("WorkerClient: unable to contact www");
        }
    }
    public void GetServer() throws InterruptedException{
        try {
           socket = new ObjectSocket(new Socket(ip, port));
           socket.InitSocket();
           if (DEBUG) cout("WorkerClient: found LagrangeServer");
           Start();
       } catch (Exception ex) {
           System.out.println("WorkerClient: cant find LagrangeServer");
           Thread.sleep(10000);
           GetServer();
       }
    }
    public void Start() throws Exception{
 
        try{
        while (true){
      
                String szMessage = socket.ReadString();
                System.out.println(szMessage);
                switch (szMessage) {
                    case "INIT_MATH_KERNEL"://creates a new mma kernel
                        InitMathKernel();
                        break;
                    case "DESTROY_MATH_KERNEL"://destroys the mma kernel
                        DestroyMathKernel();
                        break;
                    case "NEW_DATA"://send new data and save it to file also getdf, bSim, delat
                        NewData();
                        break;
                    case "PREPARE_DATA"://load lagraniansolver.m into mma and also load data
                        PrepareData();
                        break;
                    case "NEW_POPULATION"://process part of a zeitgeist
                        System.out.println("WorkerClient:NewPopulation");
                        NewPopulation();
                        break;
                    default:
                        //probably an error
                        break;

                }
        
        }
        } catch (IOException ex) {
                if (DEBUG) cout("WorkerClient: IOException");
            } catch (MathLinkException ex) {
                if (DEBUG) cout("WorkerClient: MathLinkError");
            } catch (ClassNotFoundException ex) {
                ;
            }
        //if the execution gets to this point then there has been a problem
        socket.DisconnectSocket();
        socket = null;
        throw new Exception();
    }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
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
                    port = Integer.parseInt(xml.EvaluateXPath("/server/port2")); 


                }
                
            if (DEBUG) cout("Looking for LagrangeServer at ip="+ip+" port="+port);
    }
    private void InitMathKernel() throws MathLinkException{
        if (DEBUG) cout("WorkerClient: InitMathKernel");
        scorer.InitKernel();
    }
    private void NewData() throws IOException, ClassNotFoundException, MathLinkException{
        Data data = (Data) socket.ReadObject();
        //write out this file        
        FileOutputStream fos = new FileOutputStream(szExpData);
        fos.write(data.GetBytes());
        fos.close();
        if (DEBUG) cout("WorkerClient: Downloaded experimental data from server");
        
    }
    private void PrepareData() throws MathLinkException, IOException, ClassNotFoundException{
        if (DEBUG) cout("WorkerClient: PrepareData");
        DataFormat format = (DataFormat) socket.ReadObject();
        scorer.PrepareData(szExpData, format);
    }
    private void NewPopulation() throws IOException, ClassNotFoundException, MathLinkException{
            if (DEBUG) System.out.println("LagrangeServer:NewPopulation");
            //can only run if the experimental data is available
            Population population = (Population) socket.ReadObject();
            if (DEBUG) cout(population.toString());
            population = scorer.ProcessPopulation(population);
            socket.WriteObject(population);
  
        
    }
    private void DestroyMathKernel(){
        if (DEBUG) cout("WorkerClient: DestroyMathKernel");
        scorer.DestroyKernel();
    }
    private boolean bSim;
    private double deltat;
    private int df;
    private String ip;
    private int port;    
    private ObjectSocket socket;    
    private MMAScorer scorer;
}
