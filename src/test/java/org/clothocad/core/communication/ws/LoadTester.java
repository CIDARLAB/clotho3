package org.clothocad.core.communication.ws;

import java.util.Random;
import org.eclipse.jetty.websocket.client.WebSocketClient;

public class LoadTester {

	private static final int NR_OF_CLIENTS = 1000;
	private static final int NR_OF_MESSAGES = 10000;
	
	private static final String HOST = "localhost";
	private static final int PORT = 8080;
	
	public static void main(String[] args) 
			throws Exception {
		
		WebSocketClient client = new WebSocketClient();
		client.setConnectTimeout(30000);
		client.setMaxTextMessageBufferSize(1024 * 64);
                client.start();
		
		ClothoLoadClient[] clotho = new ClothoLoadClient[NR_OF_CLIENTS];
		for (int i=0; i<NR_OF_CLIENTS; i++) {
			clotho[i]= new ClothoLoadClient("client"+i, client, HOST, PORT);
		}
		
		// Send messages
		Random random = new Random();
		for (int i=0;i<NR_OF_MESSAGES;i++) {
		    ClothoLoadClient c = clotho[random.nextInt(NR_OF_CLIENTS)];
		    
			String msg = "{\"channel\":\"say\",\"data\":{\"text\":\"This is a test message\",\"from\":\"client-"+i+"\",\"class\":\"\",\"timestamp\":1369177796313}}";
		    c.send(msg);
		}		
		
		// close all connections
		for (int i=0; i<NR_OF_CLIENTS; i++) {
		    clotho[i].disconnect();
		}
	}
}
