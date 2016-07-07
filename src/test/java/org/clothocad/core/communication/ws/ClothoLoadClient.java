package org.clothocad.core.communication.ws;

import java.io.IOException;
import java.net.URI;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.atomic.AtomicLong;
import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.*;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;
import org.eclipse.jetty.websocket.client.WebSocketClient;

@WebSocket
public class ClothoLoadClient{
	
	private static final AtomicLong sent = new AtomicLong(0);
	private static final AtomicLong received = new AtomicLong(0);
	private static final Set<ClothoLoadClient> members = new CopyOnWriteArraySet<ClothoLoadClient>();
	private final String name;
	private final Session session;
 
	public ClothoLoadClient(String username, WebSocketClient client,String host, int port)
				  throws Exception {
		name=username;
		session = client.connect(this, new URI("ws://"+host+":"+port+"/websocket")).get();
	}
 
	public void send(String message) 
			throws IOException {
		session.getRemote().sendString(message);
	}
 
        @OnWebSocketConnect
	public void onOpen(Session session) {
		members.add(this);
	}
 
        @OnWebSocketClose
	public void onClose(int closeCode, String message) {
		members.remove(this);
	}
 
        
	public void onMessage(String data) {
		received.incrementAndGet();
	}
 
	public void disconnect() 
		throws IOException {
		session.close();
	}
}