package org.clothocad;

import org.clothocad.core.jetty.ClothoWebServer;
import org.clothocad.server.ClothoServer;

public class ClothoServerStarter {
	public static void main(String[] args) {
		try {
			// first, we start the Clotho core
			new ClothoServer().start();
			
			// second, we start the web server
			new ClothoWebServer().start();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
