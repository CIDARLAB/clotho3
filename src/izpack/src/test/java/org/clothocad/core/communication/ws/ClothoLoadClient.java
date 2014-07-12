package org.clothocad.core.communication.ws;

import java.io.IOException;
import java.net.URI;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.atomic.AtomicLong;

import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketClient;

public class ClothoLoadClient 
	implements WebSocket.OnTextMessage {
	
	private static final AtomicLong sent = new AtomicLong(0);
	private static final AtomicLong received = new AtomicLong(0);
	private static final Set<ClothoLoadClient> members = new CopyOnWriteArraySet<ClothoLoadClient>();
	private final String name;
	private final Connection connection;
 
	public ClothoLoadClient(
		  String username, WebSocketClient client,String host, int port)
				  throws Exception {
		name=username;
		connection=client.open(new URI("ws://"+host+":"+port+"/websocket"),this).get();
	}
 
	public void send(String message) 
			throws IOException {
		connection.sendMessage(message);
	}
 
	public void onOpen(Connection connection) {
		members.add(this);
	}
 
	public void onClose(int closeCode, String message) {
		members.remove(this);
	}
 
	public void onMessage(String data) {
		received.incrementAndGet();
	}
 
	public void disconnect() 
		throws IOException {
		connection.disconnect();
	}
}