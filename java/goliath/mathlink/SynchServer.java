/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author User
 */
//This is part of the manager class and it wont really register
//itself
public class SynchServer implements Runnable {
    public boolean DEBUG = true;
    public SynchServer(LagrangeServer ls, int port){
        this.port = port;
        lagrangeServer = ls;
    }
    @Override
    public void finalize() throws IOException{
        if (serverSocket != null) {
            serverSocket.close();
            if (DEBUG) cout("SychServer: closed serverSocket");
        }
    }
    @Override
    public void run(){
        try {
            serverSocket = new ServerSocket(port);
            if (DEBUG) System.out.println("SychServer:Listening at port "+port);
            while (true){             
                Socket socket = serverSocket.accept();
                ObjectSocket sock = new ObjectSocket(socket);
                sock.InitSocket();
                if (DEBUG) System.out.println("SychServer:Incomming Transmission");
                try{
                while(socket.isConnected())//this will throw an exception when the socket is closed
                {
                    String szMsg = sock.ReadString();
                    if (DEBUG) System.out.println("SychServer: Received Message="+szMsg);
                    switch (szMsg) {
                        case "RESET_SERVER":
                        {
                            boolean bOK = lagrangeServer.ResetServer();
                            if (bOK)   sock.WriteString("OK");
                            else sock.WriteString("FAIL");                            
                        }
                        break;
                        case "INIT_MATH_KERNEL":
                        {                         
                            boolean bOK = lagrangeServer.InitMathKernel();
                            if (bOK)   sock.WriteString("OK");
                            else sock.WriteString("FAIL"); 
                        
                        }
                        break;
                        case "DESTROY_MATH_KERNEL":
                        {
                            boolean bOK = lagrangeServer.DestroyMathKernel();
                            if (bOK)   sock.WriteString("OK");
                            else sock.WriteString("FAIL");                         
                        }
                        break;

                        case "NEW_DATA":
                        {
                           
                            Data data = null;
                            try {
                            data = (Data)sock.ReadObject();
                            boolean bOK =  lagrangeServer.NewData(data);
                            if (bOK)   sock.WriteString("OK");
                            else sock.WriteString("FAIL");                 
                            } catch (ClassNotFoundException ex) {
                            cout("SychServer: ClassNotFoundException in run() -- some sort of unknown class type has been sent through the ObjectSocket.");
                        }
                        }
                            break;
                        
                        case "PREPARE_DATA":
                        {
                            DataFormat format = null;
                            try {
                            format = (DataFormat)sock.ReadObject();
                            boolean bOK =  lagrangeServer.PrepareData(format);
                            if (bOK)   sock.WriteString("OK");
                            else sock.WriteString("FAIL");                 
                            } catch (ClassNotFoundException ex) {
                            cout("SychServer: ClassNotFoundException in run() -- some sort of unknown class type has been sent through the ObjectSocket.");
                            }
                        }
                            
                            
                            
                        
                        break;
                        case "NEW_ZEITGEIST":
                        {
                           
                            Zeitgeist zg = null;
                            try {
                            zg = (Zeitgeist)sock.ReadObject();

                           
                            boolean bOK = lagrangeServer.NewZeitgeist(zg);
                             if (DEBUG) System.out.println("SychServer: processed zeitgeist");
                            if (bOK){
                                sock.WriteString("RESULTS");
                                sock.WriteObject(zg);
                            }
                            else {
                                sock.WriteString("FAIL");
                            }
                            } catch (ClassNotFoundException ex) {
                            cout("SychServer: ClassNotFoundException in run() -- some sort of unknown class type has been sent through the ObjectSocket.");
                        
                        }
                            
                        }
                            break;
                    }
                    
                }
                }catch (IOException ex){
                    //this exception means that the socket is closed
                   sock.DisconnectSocket();
                   sock = null;
                   if (DEBUG) System.out.println("SychServer:connection with GorillaClient closed");
                } catch (InterruptedException ex) {
                   if (DEBUG) System.out.println("SynchServer:unable to close ObjectSocket after the transmission with GorillaClient finished");
                }    
            }
        } catch (IOException ex) {
            System.out.println("SychServer: Unable to make a ServerSocket at gorillaPort. This is probably because another process is currently using that port. Is an instance of this server already running on that port?");
        }

    }
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    private LagrangeServer lagrangeServer;
    private ServerSocket serverSocket;
    private int port;
 


}
