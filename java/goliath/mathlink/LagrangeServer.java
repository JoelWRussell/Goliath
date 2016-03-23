/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import static com.lagrangianmining.Utility.registerServer;
import static com.lagrangianmining.Utility.website;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
//keeps an up to date check on the number of clients
/**
 *
 * @author User
 */
public class LagrangeServer {
    public boolean DEBUG = true;
    public static void main(String[] Args){
        LagrangeServer ls;
        if (Args.length !=0){
           ls = new LagrangeServer(Args[0], Integer.parseInt(Args[1]), Integer.parseInt(Args[2])); 
        }
        else ls = new LagrangeServer();
   
    }
    public LagrangeServer(){
        this(8000, 8001);
    }
    public LagrangeServer(int gorillaPort, int workerPort){
        this("localHost", gorillaPort, workerPort);
 
    }
    public LagrangeServer(String name, int gorillaPort, int workerPort){
 
        cout("Starting LagrangeServer");
        cout("#######################");
        cout("This server routes the connections between a machine that is running the GA (GorillaClient) \nand worker computers (WorkerClient) which share out the MMA computation.");
        cout("This servers {name, GorillaClient port, WorkerClient port} are registered at www.lagrangianmining.com. \nThis way port information and ip addresses etc can be retrieved from the website as a lookup service \nfor convenience.");
        cout("name = "+name);
        cout("gorillaPort = "+gorillaPort);
        cout("workerPort = "+workerPort);
        System.out.println("");

        try {
            RegisterServerWithWebsite(name, gorillaPort, workerPort);
        } catch (Exception ex) {
            cout("Server registered with www");
        }
        
        try {
            synchServer = new SynchServer(this, gorillaPort);
            Thread th = new Thread(synchServer);
            th.start();
        } catch (Exception ex) {
            if (DEBUG) cout("LagrangeServer: cannot start the server at:"+gorillaPort);
        }
        try {
            asynchServer = new AsynchServer(this, workerPort);
            Thread th = new Thread(asynchServer);
            th.start();
        } catch (Exception ex) {
            if (DEBUG) cout("LagrangeServer: cannot start the server at:"+workerPort);
        }
    }
    public boolean ResetServer() throws InterruptedException{
        if (DEBUG) cout("LagrangeServer:ResetServer");
        boolean bOK = DestroyMathKernel();
        asynchServer.Reset();
        return true;
        
        
    }
    public boolean InitMathKernel() throws InterruptedException{
        //This function inits all of the attached client computers
        //so it can fail if during the execution some are lost, or gained
        //it can also fail if at least one does not have a successful init
        //
        //to return true means that by the end of the function all of the attached
        //computers have init
        //
        if (DEBUG) System.out.println("LagrangeServer:InitMathKernel");
        int numClients = asynchServer.GetNumClients();
        //GetNumClients brings
        if (numClients == 0) return false;
        Thread[] threads = new Thread[numClients];
        InitMathKernelRemoteHost[] remotes = new InitMathKernelRemoteHost[numClients];
        for (int f=0; f<numClients; f++){
            
            remotes[f] = new InitMathKernelRemoteHost(asynchServer.GetSocket(f));
        }
        for (int f=0; f<numClients; f++){
            threads[f] = new Thread(remotes[f]);
            threads[f].start();
        }
        for (int f=0; f<numClients; f++){
            threads[f].join();
        }
        boolean bOK = true;
        for (int f=0; f<numClients; f++){
            bOK = bOK && remotes[f].GetStatus();
        }
        return bOK;
    }
    public boolean DestroyMathKernel() throws InterruptedException{
        if (DEBUG) System.out.println("LagrangeServer:DestroyMathKernel");
        int numClients = asynchServer.GetNumClients();
        if (numClients == 0) return false;
        Thread[] threads = new Thread[numClients];
        DestroyMathKernelRemoteHost[] remotes = new DestroyMathKernelRemoteHost[numClients];
        for (int f=0; f<numClients; f++){
            remotes[f] = new DestroyMathKernelRemoteHost(asynchServer.GetSocket(f));
        }
        for (int f=0; f<numClients; f++){
            threads[f] = new Thread(remotes[f]);
            threads[f].start();
        }
        for (int f=0; f<numClients; f++){
            threads[f].join();
        }
        boolean bOK = true;
        for (int f=0; f<numClients; f++){
            bOK = bOK && remotes[f].GetStatus();
        }
        return bOK;
    }
    public boolean NewData(Data data) throws IOException, InterruptedException{
        if (DEBUG) System.out.println("LagrangeServer:NewData");
        int numClients = asynchServer.GetNumClients();
        if (numClients == 0 || data==null) return false;
        Thread[] threads = new Thread[numClients];
        NewDataToRemoteHost[] remotes = new NewDataToRemoteHost[numClients];
        for (int f=0; f<numClients; f++){
            remotes[f] = new NewDataToRemoteHost(asynchServer.GetSocket(f), data);
        }
        for (int f=0; f<numClients; f++){
            threads[f] = new Thread(remotes[f]);
            threads[f].start();
        }
        for (int f=0; f<numClients; f++){
            threads[f].join();
        }
        boolean bOK = true;
        for (int f=0; f<numClients; f++){
            bOK = bOK && remotes[f].GetStatus();
        }
        return bOK;
    }
    public boolean PrepareData(DataFormat format) throws IOException, InterruptedException{
        if (DEBUG) System.out.println("LagrangeServer:PrepareData");
        int numClients = asynchServer.GetNumClients();
        if (numClients == 0 || format==null) return false;
        Thread[] threads = new Thread[numClients];
        PrepareDataRemoteHost[] remotes = new PrepareDataRemoteHost[numClients];
        for (int f=0; f<numClients; f++){
            remotes[f] = new PrepareDataRemoteHost(asynchServer.GetSocket(f), format);
        }
        for (int f=0; f<numClients; f++){
            threads[f] = new Thread(remotes[f]);
            threads[f].start();
        }
        for (int f=0; f<numClients; f++){
            threads[f].join();
        }
        boolean bOK = true;
        for (int f=0; f<numClients; f++){
            bOK = bOK && remotes[f].GetStatus();
        }
        return bOK;
    }
    public boolean NewZeitgeist(Zeitgeist zg) throws InterruptedException{
        if (DEBUG) System.out.println("LagrangeServer:NewZeitgeist");
        //cout("#####################################################");
        //cout("STARTING ZEITGEIST");
        //cout(zg.toString());
        //cout("#####################################################");
  
        int numClients = asynchServer.GetNumClients();
        if (numClients == 0 || zg==null) return false;
        if (DEBUG) System.out.println("LagrangeServer: numClients = " + numClients);
        zg.SetNumberClients(numClients);  
        //cout("#####################################################");
        //cout("REORGANISED ZEITGEIST WITH NUMCLIENTS");
        //cout(zg.toString());
        //cout("#####################################################");
        Thread[] threads = new Thread[numClients];
        ScoreRemoteHost[] remotes = new ScoreRemoteHost[numClients];
        for (int f=0; f<numClients; f++){
            remotes[f] = new ScoreRemoteHost(asynchServer.GetSocket(f), zg.GetPopulation(f));
        }
        //System.out.println("LAGRANGESERVER-START ALL THREADS");
        for (int f=0; f<numClients; f++){
            threads[f] = new Thread(remotes[f]);
            threads[f].start();
        }
        for (int f=0; f<numClients; f++){
            threads[f].join();
        }
        boolean bOK = true;
        for (int f=0; f<numClients; f++){
            bOK = bOK && remotes[f].GetStatus();
        }
        //System.out.println("LAGRANGESERVER-FINISHED ALL THREADS");
        if (bOK){
            for (int f=0; f<numClients; f++){
                zg.SetPopulation(f, remotes[f].GetPopulation());
            }
            zg.FinalizeZeitgeist();//copies from Population to Poly
        }
        //cout("#####################################################");
        //cout("FINAL ZEITGEIST");
        //cout(zg.toString());
        //cout("#####################################################");
        //try {
        //    zg.SaveToFile("zeit3.txt");
        //} catch (Exception ex) {
        //    cout("save to file error");
        //}
        cout("SIZE OF ZEITGEIST "+zg.GetNumPolys());
        return bOK;
        
    }
    private void RegisterServerWithWebsite(String name, int gorillaPort, int workerPort) throws IOException, Exception{
        //should I use try-with resources to close these 
        //connections or does it happen automatically?
            XMLHelper xml = new XMLHelper();          
            CloseableHttpClient httpclient = HttpClients.createDefault();
            String szAddr = website+registerServer+"?"+"server=" + name + "&port1="+gorillaPort + "&port2="+workerPort;
            HttpGet httpget = new HttpGet(szAddr);
            CloseableHttpResponse response = httpclient.execute(httpget); 
               
                HttpEntity entity = response.getEntity();
                if (entity != null){
                    String sz = EntityUtils.toString(entity);
                    xml.setDocument(sz);
                    String status = xml.EvaluateXPath("/status");
                    if (!status.equals("OK")){
                        httpclient.close();
                      
                    }


                }
            httpclient.close();
           
    }
    private AsynchServer asynchServer;
    private SynchServer synchServer;

    private class PrepareDataRemoteHost implements Runnable{
    
        public PrepareDataRemoteHost(ObjectSocket socket, DataFormat format){
            this.socket = socket;
            this.format = format;
            bStatus = true;
        }
        @Override public void run(){
            try {

                //send the population off to the client
                socket.WriteString("PREPARE_DATA");
                socket.WriteObject(format);

            } catch (Exception ex){
                bStatus = false;
            }
        }
        public boolean GetStatus(){ return bStatus;}
 
        
        private ObjectSocket socket;
        private DataFormat format;
        private boolean bStatus;
    
} 
    private class ScoreRemoteHost implements Runnable{
    
        public ScoreRemoteHost(ObjectSocket socket, Population population){
            this.socket = socket;
            this.population = population;
            bStatus = true;
        }
        @Override public void run(){
            try {
                socket.WriteString("NEW_POPULATION");
                socket.WriteObject(population);
                population = (Population)socket.ReadObject();
               // if (DEBUG) cout(population.toString());
            } catch (Exception ex){
                bStatus = false;
            }
        }
        public boolean GetStatus(){ return bStatus;}
        public Population GetPopulation(){return population;}
        
        private ObjectSocket socket;
        private Population population;
        private boolean bStatus;
    
}
    private class NewDataToRemoteHost implements Runnable{
        public NewDataToRemoteHost(ObjectSocket socket, Data data){
            sock = socket;
            this.data = data;
            bStatus = true;
        }
        @Override public void run(){
           try {
               sock.WriteString("NEW_DATA");
               sock.WriteObject(data);
              
               }
             catch (IOException ex) {
                cout("Error with socket");
                sock.DisconnectSocket();
            } 
        }
        public boolean GetStatus(){
            return bStatus;
        }
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        private ObjectSocket sock;
        private Data data;
        private boolean bStatus;
       
    } 
    private class InitMathKernelRemoteHost implements Runnable{
        public InitMathKernelRemoteHost(ObjectSocket socket){
            sock = socket; 
            bStatus = true;
        }
        @Override public void run(){
           try {
               sock.WriteString("INIT_MATH_KERNEL");
              
               }
             catch (IOException ex) {
                cout("Error with socket");
                sock.DisconnectSocket();
            } 
        }
        public boolean GetStatus(){
            return bStatus;
        }
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        private ObjectSocket sock;
        private boolean bStatus;

       
    } 
    private class DestroyMathKernelRemoteHost implements Runnable{
        public DestroyMathKernelRemoteHost(ObjectSocket socket){
            sock = socket;
            bStatus = true;
        }
        @Override public void run(){
           try {
               sock.WriteString("DESTROY_MATH_KERNEL");
              
               }
             catch (IOException ex) {
                cout("Error with socket");
                sock.DisconnectSocket();
            } 
        }
        public boolean GetStatus(){
            return bStatus;
        }
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        private ObjectSocket sock;
        private boolean bStatus;

       
    } 
 

}
