package org.clothocad.core.layers.communication.connection.apollo;

import java.util.UUID;

import javax.jms.Session;

import org.clothocad.core.layers.communication.connection.ClientConnection;

public class ApolloConnection 
	extends ClientConnection {

	private Session session;
	private String correlationId;
	
	public ApolloConnection(String id, Session session, String correlationId) {
		super(id);
		this.session = session;
		this.correlationId = correlationId;
	}
	
	public ApolloConnection(Session session, String correlationId) {
		super(UUID.randomUUID().toString());
		this.session = session;
		this.correlationId = correlationId;
	}

	public Session getSession() {
		return this.session;
	}
	
	public String getCorrelationId() {
		return this.correlationId;
	}
}
