package org.clothocad.core.communication.ws;

import java.net.URI;
import java.util.concurrent.TimeUnit;

import org.clothocad.core.util.FileUtils;
import org.eclipse.jetty.util.ssl.SslContextFactory;
import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.*;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;
import org.eclipse.jetty.websocket.client.WebSocketClient;
import java.util.concurrent.Future;

public class ClothoWebSocketClient {

	private final String CLOTHO_WEBSOCKET_LOCATION = "ws://localhost:9090/websocket";
        
        @WebSocket
        private class ClothoWebSocket {
            @OnWebSocketConnect
            public void onOpen(Session session) {
                    // open notification
            }
            @OnWebSocketClose
            public void onClose(int closeCode, String message) {
                    // close notification
            }
            @OnWebSocketMessage
            public void onMessage(String data) {
                    System.out.println("[ClothoWebSocketClient.onMessage] RECEIVED MESSAGE "+data);
            }
        }
	
	public ClothoWebSocketClient(String jsonFile) 
			throws Exception {
		
		SslContextFactory factory = new SslContextFactory(true);
		WebSocketClient client = new WebSocketClient(factory);
                client.setMaxBinaryMessageBufferSize(999999);
                client.setMaxTextMessageBufferSize(999999);
		Future<Session> fut = client.connect(new ClothoWebSocket(), new URI(CLOTHO_WEBSOCKET_LOCATION));
//		WebSocket.Connection connection = client.open(
//			new URI(CLOTHO_WEBSOCKET_LOCATION), 
//			new WebSocket.OnTextMessage() {					
//                                public void onOpen(Connection connection) {
//                                        // open notification
//                                }
//
//                                public void onClose(int closeCode, String message) {
//                                        // close notification
//                                }
//
//                                public void onMessage(String data) {
//                                        System.out.println("[ClothoWebSocketClient.onMessage] RECEIVED MESSAGE "+data);
//                                }
//			}
//		).get(5, TimeUnit.SECONDS);
				
		String json = FileUtils.readFile(jsonFile);
		System.out.println("sending...");
		System.out.println(json);
                
		fut.get(5, TimeUnit.SECONDS).getRemote().sendString(json);
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
