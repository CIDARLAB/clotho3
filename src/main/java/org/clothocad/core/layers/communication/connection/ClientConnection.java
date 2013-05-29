package org.clothocad.core.layers.communication.connection;

public abstract class ClientConnection {

	private String id;
	public ClientConnection(String id) {
		this.id = id;
	}
	public String getId() {
		return this.id;
	}
}
