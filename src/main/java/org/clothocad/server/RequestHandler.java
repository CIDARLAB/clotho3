package org.clothocad.server;

import java.net.Socket;

public class RequestHandler {

	private Socket clientSocket;
	
	public RequestHandler(Socket clientSocket) {
		this.clientSocket = clientSocket;
	}
}
