package org.clothocad.core.communication.ws;

import java.net.URI;
import java.util.concurrent.TimeUnit;

import org.clothocad.core.util.FileUtils;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketClient;
import org.eclipse.jetty.websocket.WebSocketClientFactory;

public class ClothoWebSocketClient {

	private final String CLOTHO_WEBSOCKET_LOCATION = "ws://localhost:9090/websocket";
	
	public ClothoWebSocketClient(String jsonFile) 
			throws Exception {
		
		WebSocketClientFactory factory = new WebSocketClientFactory();
		factory.start();		
		
		WebSocketClient client = factory.newWebSocketClient();
		
		WebSocket.Connection connection = client.open(
			new URI(CLOTHO_WEBSOCKET_LOCATION), 
			new WebSocket.OnTextMessage() {					
                                public void onOpen(Connection connection) {
                                        // open notification
                                }

                                public void onClose(int closeCode, String message) {
                                        // close notification
                                }

                                public void onMessage(String data) {
                                        System.out.println("[ClothoWebSocketClient.onMessage] RECEIVED MESSAGE "+data);
                                }
			}
		).get(5, TimeUnit.SECONDS);
				
		String json = FileUtils.readFile(jsonFile);
		System.out.println("sending...");
		System.out.println(json);
		
		connection.sendMessage(json);
	}
	
	public static void main(String[] args) 
			throws Exception {
		System.out.println("**** SAY ****");
		ClothoWebSocketClient wsClient = new ClothoWebSocketClient("./clotho3-web/json-examples/say.json");
		
		// 1. connect to the websocket
		System.out.println("**** CREATE ****");
		wsClient = new ClothoWebSocketClient("./clotho3-web/json-examples/create-schema.json");		
	}

}
