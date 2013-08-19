/*
 * 
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.aspects.Ambassador;

import org.clothocad.core.aspects.Hopper;
import org.clothocad.core.aspects.Ambassador.Sister.TrustLevel;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.objbases.Badge;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashMap;
import java.util.Map;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.util.JSON;

/**
 * The Connector and Communicator have somewhat similar functions but not identical.
 * 
 * Communicator is for speaking to a client device, like an ipad, smartphone, or browser.
 * As such, it involves a single user at a time and the communication is all the API stuff
 * but also a bunch of visual stuff like GUI tracking.  It is a high-resolution channel
 * of communication which assumes that both the client and server have a fully-entrusted
 * relationship.
 * 
 * Ambassador is for speaking server to server or clothocore to clothocore.  As such, it is
 * just for CRUD operations, relaying of search requests, and confirmation of badges (and possibly
 * a few other things).
 * But, it doesn't involve GUIs.  The critical difference is that the Connector does not assume
 * its partner is trustworthy, and it is a many-to-many people interaction rather than a one-to-many
 * relationship like you get with a client.
 * 
 * Currently this is connecting via TCP/IP, but eventually this should do something like Router
 * to abstract server-t0-server communication to allow more mechanisms of communication.
 * 
 * This is for clotho-to-clotho as well as clotho-to-other-type-of-server communications.  
 * 
 * @author John Christopher Anderson
 * @author Kevin Meng
 */

@Slf4j
public class Ambassador implements Aspect {
    
    Persistor persistor;
    
    private Ambassador() {
        System.out.println("The ambassador has been initiated");
        
        //Initiate the singular Thread that is in an endless loop to intercept TCP/IP messages
        try {
            serverThread = new Thread(
                new Runnable() {
                    @Override
                    public void run() {
                        while(true) {
                            listen();
                        }
                    }
                } 
                );
            serverThread.start();
            
            final Thread counterThread = new Thread();
            Thread runnerThread = new Thread(
                new Runnable() {
                    @Override
                    public void run() {
                        for(;;) {
                            try{
                                counterThread.sleep(1000);
                                count++;
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }
                );

            int startPortRange = 49152;
            int stopPortRange = 49162;

            for (int i = startPortRange; i <= stopPortRange; i++) {
                try {
                    Socket ServerSok = new Socket("127.0.0.1", i);
                    System.out.println("Port usable: " + i);
                    port = i;
                    ServerSok.close();
                    break;
                } catch (Exception e) {
//                    e.printStackTrace();
                    System.out.println("Port cannot be used: " + i);
                }
            }
            
            server = new ServerSocket(port);
            System.out.println(server.getLocalSocketAddress());

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * The listen method is called within an endless while loop
     * within a single thread.  That thread is stalled until a message
     * is received, at which point it reads the message, responds, then
     * goes back to waiting.
     * 
     * TO DO:  this is the asynchronous receiving thing that is supposed to do
     * a Node-like pattern - the non-blocking bit.  In this case, I think that 
     * means it should do this stalling business as it is, but when it gets a 
     * message it does a Process.nextTick() deal to put everything on one singular
     * Executor thread, that is owned by Executor and controls all tasks received
     * by the order they are received, and executes those tasks sequentially and
     * synchronized.  The nextTick() of an action listener adds that action to the
     * top of the heap.  That action listener object should be a doo, so this is
     * a Doo hopper scenario.
     */
    
    public String getIP() {
        try {
            URL whatismyip = new URL("http://automation.whatismyip.com/n09230945.asp");
            URLConnection connection = whatismyip.openConnection();
            connection.addRequestProperty("Protocol", "Http/1.1");
            connection.addRequestProperty("Connection", "keep-alive");
            connection.addRequestProperty("Keep-Alive", "1000");
            connection.addRequestProperty("User-Agent", "Web-Agent");
            InputStreamReader input = new InputStreamReader(whatismyip.openStream());

            BufferedReader in = new BufferedReader(input);
            String ip = in.readLine(); //you get the IP as a String
            System.out.println(ip);
            return ip;
         } catch (IOException e) {
            e.printStackTrace();
            return "No IP";
         }
    }
    private void listen() {
        boolean a = false;
    	try {
            Socket socket = server.accept();
            ObjectInputStream ois = new ObjectInputStream(socket.getInputStream());
            String response = (String) ois.readObject();
            receiveClothoMessage(socket.getInetAddress().getHostName(), response);
            ois.close();
            a = true;
      	} catch (Exception e) {
      		e.printStackTrace();
    	}
        if (count % 30 == 0) {
            if (a) {
                sendClothoMessage(null, null);     // Try to send to the port it is currently listening to the new IP address
            }
        }
    }
    
    /**
     * Send an outgoing json message
     * 
     * @param url  the internet address of the Clotho you are calling
     * @param port the port at the address of the Clotho you are calling
     * @param message a json message
     */
    public void sendClothoMessage(String url,
                                  String message) {
    	Socket socket = null;
    	try {
          //
          // Create a connection to the server socket on the server application
          //
          socket = new Socket(url, port);
          //
          // Send a message to the client application
          //
          ObjectOutputStream oos = new ObjectOutputStream(socket.getOutputStream());
          oos.writeObject(message);
          System.out.println("Message: " + message);
          
          oos.close();
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	//
        // to do:  relay the message onto the url/port.  I don't know if port
        //needs to be included, or if it can always be the same.
    }
    
    /**
     * The interpretation and response to an incoming message.
     * You probably shouldn't be calling this.
     * 
     * @param url the URL
     * @param message a json message
     */
    private void receiveClothoMessage(String url,
                                     String message) {
            	
        try {
            //
            //Read the message and pull any waiting Doos
            //
            Map<String, Object> jsonObj = JSON.deserializeObject(message);
            String callerId = null;
            Doo waitingDoo = null;
            if(jsonObj.containsKey("calling_doo_id")) {
                callerId = (String) jsonObj.get("calling_doo_id");
                waitingDoo = Hopper.get().extract(callerId);
            }
            
            //
            //Instantiate a Doo to manage this request
            //
            AmbassadorDoo doo = new AmbassadorDoo(waitingDoo, url, message);
            doo.jsonObj = jsonObj;
            doo.callerDooId = callerId;
            doo.url = url;
            
            //
            //Do some bookeeping about the Sister
            //
            doo.sisterId = doo.jsonObj.get("sister_id").toString();
            doo.port = Integer.parseInt(doo.jsonObj.get("port").toString());
            Sister sister = persistor.get(Sister.class, new ObjectId(doo.sisterId));
            if(sister==null) {
                //Create a new Sister to represent this newcomer
                sister = new Sister(doo.sisterId, doo.url, doo.port);
                doo.securityCode = sister.getSecurityCode();
            } else {
                //If I've bloacked that clothocore, terminate the doo and abort
                if(sister.getTrust().equals(TrustLevel.BLOCKED)) {
                    doo.setMessage("Aborted request since the Sister is blocked");
                    doo.terminate();
                    return;
                }
                
                sister.setPort(doo.port);
                sister.setUrl(doo.url);
            }
            
            persistor.save(sister);
            
            //
            //Extract any security codes, and the cmd token
            //
            doo.cmd = (String) doo.jsonObj.get("cmd");
            if(doo.jsonObj.containsKey("security_code")) {
                doo.securityCode = (String) doo.jsonObj.get("security_code");
            }
            
            //
            //Redirect the message using the cmd token
   //******** this is effectively the server-to-server API ********//
            //
            if(doo.cmd.equals("responding")) {
                reviveDooFromHopper(doo);
            } else if(doo.cmd.equals("get")) {
                respondToGet(doo);
            } else if(doo.cmd.equals("has_badge")) {
                respondToHasBadge(doo);
            }
            
            //
            //etcetera, this all defines the api and its form
            //the sent message object needs to be built up in a case-specific way, so
            //it probably can't be meaningfully abstracted, so things just call
            //the sendClothoMsg(...) directly.  So, that should be a public method
            //of the Aspect
            //
            
        } catch (Exception ex) {
            System.out.println("Ambassador (Exception ex) Response to: " + message);
//            Logger.getLogger(Ambassador.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void reviveDooFromHopper(AmbassadorDoo doo) {
        try {
            String waitingDooId = doo.callerDooId;
            AmbassadorDoo waitingDoo = (AmbassadorDoo) Hopper.get().extract(waitingDooId);
            if(waitingDoo==null) {
                log.error("THIS SHOULD NEVER EVER HAPPEN!");
                return;
            }

            //The waiting doo's commands should now take control
            doo.terminate();
            waitingDoo.command.run();
            
            //Say thank you
            //TO BE ADDED
            waitingDoo.terminate();
        
        } catch(Exception err) {
            doo.setMessage("failed revival from the Hopper");
            doo.terminate();
        }
        
    }
    
    /**
     * Encodes the response to the question:
     * 
               Can I get the Sharable with this Id?
     
     * @param doo 
     */
    private void respondToGet(AmbassadorDoo doo) {
        try {
            String serverId = doo.securityCode;  //The uuid of the server--the uuid that this clotho gave it
        //so, that will be null or blank if this clotho is a stranger, and it should be in a "protected" sort of mode and should be tracked
        
        //if that server_id is null/blank, then the caller cannot be trusted at all, so only public stuff gets included
        //you can't know that it is honestly someone's name calling so, yeah, totally just public stuff if that's blank
        
        //if that server_id is one I know, I probably have a rating for it, like, fully trusted, somewhat, and then blocked
        
            //
            //If I don't like this clotho, ping back.
            //
            if(serverId!=null && !serverId.equals("") ) {
                singleton.sendClothoMessage(doo.url, "rejected");
            }
        
            //
            //Go fetch the sharable, wrap it into the response, then send it back
            //
            doo.sharableJSON = persistor.getAsJSON(new ObjectId(doo.itemId));
            doo.response = new HashMap<>();
            doo.response.put("sharable_item", doo.sharableJSON);
            
              //NEED TO ALSO put in a doo routing number, maybe some secutity in here
            sendClothoMessage(doo.url, doo.response.toString());
            
            //Put it in the Hopper to await a callback to confirm receipt
            Hopper.get().add(doo);
        } catch(Exception ex) {
            log.error("respondToGet:", ex);
        }
        
    }
    
    /**
     * Encodes the response to the question:
     * 
               Does this person have this badge?
     
     Doo in format:
        int port: 7777
        String rawMsg:      N/A but present
        Map<String, Object> jsonObj; N/A but present
        String cmd;         has_badge
        String url;         http://andersonlab.qb3.berkeley.edu/clotho
        String item_id:     static_badge_isuuid
        String caller_id;   static_admin_isuuid
     * 
     * 
     * @param doo 
     */
    private void respondToHasBadge(AmbassadorDoo doo) {
        try {
            //
            //Go fetch the Badge in question and check if person has that badge
            //
            String badgeId = doo.itemId;
            Badge badge = persistor.get(Badge.class, new ObjectId(badgeId));
            boolean out = badge.hasBadge(doo.callerDooId);
            
            //
            //Wrap it up and return it
            //
            doo.response = new HashMap<>();
            doo.response.put("has_badge", out);
            this.sendClothoMessage(doo.url, doo.response.toString());
            doo.terminate();
        } catch(Exception ex) {
            log.error("respondToHasBadge:", ex);
        }
    }

    
    private class AmbassadorDoo extends Doo {

        private AmbassadorDoo(Doo parent, String url, String message) {
            super(parent, true); //Always log these communications since this is where all trouble emenates, so 'true'
            this.url = url;
            this.rawMsg = message;
        }
        
        int port;               //the port it came in on, not sure this is important
        String rawMsg;          //the original String sent in
        Map<String, Object> jsonObj;     //the parsing of the message
        Map<String, Object> response;    //this doo's response message
        String cmd;             //the command, like "get" "set" "check_badge", etc.
        String url;             //Th url that called me, I should have it always update that in a registry
        String sisterId;        //The clotho-wide UUID of the calling ClothoCore
        String securityCode;    //That clotho's uuid in Clotho land, not its url which could change
        String callingDooId;    //The uuid of the Doo making the call
        String itemId;          //for gets and sets, etc, the uuid of the Sharable in question, or a badge Id
        String callerDooId;        //for gets and sets, etc, the uuid of the Person making the request
        Map<String, Object> sharableJSON;//the json of the sharable that was sent as that id
        Runnable command; //the actions that should be taken when the Doo is revived from the Hopper
        
    }
    
    public void setCount(int newCount) {
        count = newCount;
    }

    public static Ambassador get() {
        return singleton;
    }
    
    public static void main(String[] args) {
    	get();
//        try {
//        	InetAddress host = InetAddress.getLocalHost();
//                System.out.println("your host is " + host.getHostAddress());
//        	get().sendClothoMessage("136.152.166.251", "hi there");
//    	} catch (UnknownHostException e) {
//            e.printStackTrace();
//    	}
    }
    
    
    private static final Ambassador singleton = new Ambassador();
    private Thread serverThread;
    private ServerSocket server;
    private int port = 1933;
    private int count = 0;

}