/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import static com.lagrangianmining.Utility.cout;
import java.io.IOException;
import java.net.ServerSocket;
import java.util.ArrayList;

/**
 *
 * @author User
 */
//lagrangeserver is called when sockets is and is checked for change in numClients
public class AsynchServer implements Runnable{
    public boolean DEBUG = true;
    //Server---registers name:port with website to advertise for clients
    //throws exception is this contact with website fails
    public AsynchServer(LagrangeServer ls, int port){
        lagrangeServer = ls;
        this.port = port;  
        bGotData = false;
        sockets = new ArrayList<ObjectSocket>();
    }
    public void finalize(){
        try {
            //close all sockets
            for(ObjectSocket s:sockets){
                s.DisconnectSocket();//throw if cant close-so we know that something went wrong
            }
            serverSocket.close();
        } catch (IOException ex) {
            if (DEBUG) cout("AsychServer: Was unable to close the ServerSocket or the ObjectSockets in Finalize()");
        }
    }
    public void Reset(){
        for (ObjectSocket sock : sockets){
            sock.DisconnectSocket();
        }
        sockets = new ArrayList<ObjectSocket>();
    }
    @Override public void run(){

        try {

            serverSocket = new ServerSocket(port);
            if (DEBUG) cout("AsychServer:Listening at port "+port);
            while (true){
                try{
                ObjectSocket socket = new ObjectSocket(serverSocket.accept());
                socket.InitSocket();
                if (DEBUG) cout ("AsychServer:received an incomming connection");
                InitRemoteHost hst = new InitRemoteHost(socket);
                Thread thd = new Thread(hst);
                thd.run();
                }catch (IOException ex){
                   if (DEBUG) cout("AsychServer: ObjectSocket exception in run()"); 
                }
            }
        } catch (IOException ex) {
            if (DEBUG) cout("AsychServer: serverSocket exception in run()");
        }

    }
    public ObjectSocket GetSocket(int f){
        return sockets.get(f);

    }
    public void RemoveClient(ObjectSocket obj){
        sockets.remove(obj);
    }
    public int GetNumClients(){
        return sockets.size();
    }
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    private int port;//port to listen on for client requests
    private boolean bGotData;
    private ServerSocket serverSocket;
    private ArrayList<ObjectSocket> sockets;
    private LagrangeServer lagrangeServer;
    private class InitRemoteHost implements Runnable{
        public InitRemoteHost(ObjectSocket socket){
          sock = socket;  
        }
        @Override public void run(){

          sockets.add(sock);
          if (DEBUG) cout("AsynchServer: Added workerClient, numClients=" + sockets.size());

        }
        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
        private ObjectSocket sock;
       
    } 
    
}

