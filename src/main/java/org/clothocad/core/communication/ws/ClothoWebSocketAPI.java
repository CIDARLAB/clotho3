package org.clothocad.core.communication.ws;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.net.URI;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketClient;
import org.eclipse.jetty.websocket.WebSocketClientFactory;

@Slf4j
public class ClothoWebSocketAPI {
    // simplistic serverAPI methods done here
    // singleton for the client websocket to be implemented
    
    public static void getFromSocket(String uuid, String uri){
        //establish websocket connection --> decompose from the uuid
        
        System.out.println("The getFromSocket is being called");
        try{
              //this should be a thread 
            ClothoWebSocketServer ws = ClothoWebSocketServer.getInstance();
//            System.out.println("Is ClothoWebsocket open? " + ws.isOpen());
//            System.out.println("The Id of the ClothoWebsocket: " + ws.getId());
            
            WebSocket.Connection connect = ws.getConnection(uri); 
            String getCommand = "{\"channel\":\"get\", \"data\":\""+ uuid + "\"}";
//            System.out.println("Command sent to the websocket: " + getCommand);
            connect.sendMessage(getCommand);    
            
        }catch(Throwable t){
            t.printStackTrace();
        }
        
    }
    public static String buildUri(String rawUri){
        String output = "wss://";
        Boolean changePrefix = rawUri.contains("http://") || rawUri.contains("https://");
        if(changePrefix){
            int position = rawUri.indexOf("//")+2;
            output = output.concat(rawUri.substring(position));
        }else{
            output = output.concat(rawUri);
        }
        output += "/websocket";
        return output;
    }

   public static Map<String,Object> receiveFromSocket(String uri){ //return Map<String, Object>
       ClothoWebSocketServer ws = ClothoWebSocketServer.getInstance();
       //set time limit
       int count = 0;
       while(!ws.gotMessage(uri) || count <10){
           count++;
           try {
               Thread.sleep(1000);
           } catch (InterruptedException ex) {
               Logger.getLogger(ClothoWebSocketAPI.class.getName()).log(Level.SEVERE, null, ex);
           }
       }
       if(ws.gotMessage(uri)){
           System.out.println(ws.getMessage(uri));
           return new HashMap<String, Object>();
       }else{
           //message not received
           return null;
       }
       
   }
   
   public static void terminateConnection(String uri){
       
   }
}
