package org.clothocad.core.layers.communication.connection.apollo;

import java.util.UUID;

import javax.jms.Session;

import org.clothocad.core.layers.communication.connection.ClientConnection;

public class ApolloConnection 
	extends ClientConnection {

	private Session session;
	
	public ApolloConnection(String id, Session session) {
		super(id);
		this.session = session;
	}
	
	public ApolloConnection(Session session) {
		super(UUID.randomUUID().toString());
		this.session = session;
	}

	public Session getSession() {
		return this.session;
	}
}
