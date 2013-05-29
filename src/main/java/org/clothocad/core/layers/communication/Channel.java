package org.clothocad.core.layers.communication;

public enum Channel {
	EXECUTION,
	AUTOCOMPLETION,
	NOTIFICATION,
	ACCESS,
	UPDATES,   // a channel where the server can push updates to subscribed clients
	RESPONSE
}
