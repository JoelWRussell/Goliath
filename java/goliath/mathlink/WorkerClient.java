/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import static com.lagrangianmining.Utility.getServer;
import static com.lagrangianmining.Utility.szExpData1;
import static com.lagrangianmining.Utility.szMath;
import static com.lagrangianmining.Utility.szMathPath1;
import static com.lagrangianmining.Utility.szPackagePath1;
import static com.lagrangianmining.Utility.szWorkerConfig;
import static com.lagrangianmining.Utility.website;
import com.wolfram.jlink.MathLinkException;
import java.awt.Toolkit;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;
import static javafx.application.Platform.exit;
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
//keep an error file so that we know
//what the mathematica script failed at 
public class WorkerClient {
    public boolean DEBUG = true;
    public static void main (String[] Args){
        try{
        WorkerClient client;
        if (Args.length !=0){
           client = new WorkerClient(Args[0]); 
        }
        else client = new WorkerClient();

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
        szMathPath = szMathPath1;
        szPackagePath = szPackagePath1;
        szExpData = szExpData1;
        
        cout("Starting WorkerClient");
        cout("#####################");
        cout("Worker client will init Mathematica kernel and process the population.");
        cout("");
        ReadConfigFile(szWorkerConfig);
        cout("Read configuration file:"+szWorkerConfig+"(defaults in Utility.java)");
        cout("MMAKernel="+szMathPath+"\nPackagePath=" + szPackagePath+"\nExperimentalData="+szExpData);
        cout("");
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
                    case "NEW_DATA"://
                        NewData();
                        break;
                    case "PREPARE_DATA"://load lagraniansolver.m into mma and also load data
                        PrepareData();
                        break;
                    case "NEW_POPULATION"://process part of a zeitgeist
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
    public void ReadConfigFile(String sz){
        File file = new File(sz);
        try(BufferedReader br = new BufferedReader(new FileReader(file))) {
            for(String line; (line = br.readLine()) != null; ) {
                String[] splited = line.split(",");
                if (splited.length != 2) continue;
                switch (splited[0]) {
                    case "MathKernelPath":
                        szMathPath = splited[1];
                        break;
                    case "PackagePath":
                        szPackagePath = splited[1];
                        break;
                    case "ExperimentalData":
                        szExpData = splited[1];
                        break;
                }
            }

        }catch (IOException ex) {
             if (DEBUG) cout("WorkerClient: Cant find the configuration file, this contains info like the location \nof the math kernel, using default values in Utility.java.");
    
        }
        
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    private String szMathPath, szPackagePath, szExpData;
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
        FileOutputStream fos = new FileOutputStream(szPackagePath + szExpData);
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
            if (DEBUG) System.out.println("WorkerClient:NewPopulation");
            //can only run if the experimental data is available
            Population population = (Population) socket.ReadObject();
            scorer.ProcessPopulation(population);
            if (DEBUG) System.out.println("Population size: "+population.polys.size());
            //cout("WORKERCLIENT:START_SLEEP");
            //Toolkit.getDefaultToolkit().beep();
            socket.WriteObject(population);
            cout("WORKERCLIENT:PROCESSED_POPULATION");
  
        
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
