package org.clothocad.core.layers.communication.apollo;

import java.util.UUID;

import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.Session;
import org.clothocad.core.layers.communication.Communicator;

import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.clothocad.core.layers.communication.mind.Mind;
import org.fusesource.stomp.jms.message.StompJmsMessage;
import org.json.JSONException;
import org.json.JSONObject;

public class ClothoMessageConsumer 
	implements Runnable {
	
	private Session session;
	private StompJmsMessage message;
	
	public ClothoMessageConsumer(Session session, Message message) 
			throws Exception {
		this.session = session;
		
		if(message!=null && message instanceof StompJmsMessage) {			
			this.message = (StompJmsMessage)message;
		} else {
			throw new Exception("INVALID MESSAGE!");
		}
	}
	
	@Override
	public void run() {
            try {
                if(this.message.propertyExists("request")) {

                    JSONObject json = new JSONObject(
                            message.getStringProperty("request"));

                    // get the message's correlation id
                    String sCorrelationID = this.message.getJMSCorrelationID();
                    if(null == sCorrelationID) {
                        sCorrelationID = UUID.randomUUID().toString();
                        this.message.setJMSCorrelationID(sCorrelationID);
                    }

				
                    // HERE, we need to contact the mind about the connection... 
                    String auth_key = null;
                    try {
                         auth_key = json.getString("auth_key");
                    } catch(Exception error) {
                        error.printStackTrace();
                        JSONObject responseJSON = new JSONObject();
                        responseJSON.put("CLOTHO-ERROR", error.getMessage());
                        return;
                    }
                    
                    if(auth_key != null) {
                        ApolloConnection connection = null;
                        Mind mind = Communicator.get().getMind(auth_key);
                        if(null == mind.getClientConnection()) {
                            connection = new ApolloConnection(
                                                    session);
                        } else {
                            ClientConnection conn = mind.getClientConnection();
                            if(conn instanceof ApolloConnection) {
                                connection = (ApolloConnection)conn;
                            } else {
                                Router.get().receiveMessage(
                                        connection,
                                        json.getString("channel"), 
                                        json);
                                return;
                            }
                        }
                        
                        // store the callback-handler in the callback-handler table
                        CallbackHandlerTable.put(
                                connection.getId(), 
                                new CallbackHandler(this.session, this.message));
                        
                        // route the message
                        Router.get().receiveMessage(
                            connection,
                            json.getString("channel"), 
                            json);
                    } else {
                        // this needs to be clarified...
                    }
                }
        } catch (JMSException e) {
            e.printStackTrace();
        } catch (JSONException e) {
            e.printStackTrace();
        }
    } 
}
