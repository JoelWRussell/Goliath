/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.lagrangianmining;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author User
 */
public class ObjectSocket {
    private boolean DEBUG = true;
    public static void main(String[] Args){
        try {
            try (Socket socket = new Socket("127.0.0.1", 8000)) {

                ObjectSocket sock = new ObjectSocket(socket);
                sock.InitSocket();
            }
        } catch (IOException ex) {
            System.out.println("error");
        }
    }
    public ObjectSocket(Socket socket){
        try {
            //this should be given an active socket
            this.socket = socket;
            oos = new ObjectOutputStream(socket.getOutputStream());
            oos.flush();
 
           

        } catch (IOException ex) {
            if (DEBUG) System.out.println("ObjectSocket: Error getting oos, ois");
        }
    }
    public void InitSocket() throws IOException{
        ois = new ObjectInputStream(socket.getInputStream());
    }
    public void DisconnectSocket(){
        try {
            socket.close();
        } catch (IOException ex) {
            ;
        }
    }
    public boolean IsConnected(){
        return socket.isConnected();
    }
    @Override public void finalize(){
        DisconnectSocket();
    }
    public String ReadString() throws IOException{
        return ois.readUTF();
    }
    public void WriteString(String sz) throws IOException{
        oos.writeUTF(sz);
        oos.flush();
    }
    public Object ReadObject() throws IOException, ClassNotFoundException{
        return (ois.readObject());
    }
    public void WriteObject(Object obj) throws IOException{
        oos.writeObject(obj);     
        oos.flush();
       
    }
    
    private ObjectInputStream ois;
    private ObjectOutputStream oos;
    private Socket socket;
}
