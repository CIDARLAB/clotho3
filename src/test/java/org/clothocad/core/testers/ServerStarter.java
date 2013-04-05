package org.clothocad.core.testers;

import org.clothocad.server.ClothoCore;

public class ServerStarter {

	public static void main(String[] args) {

    	/**
    	int nPort = 8161;
    	
    	if(args.length > 1) {
    		System.err.println("Usage: java -jar clotho-server.jar [<port>]");
    		System.err.println("       The default port is 8161.");
    		System.exit(1);
    	} else if (args.length == 1) {
	    	try {
	    		nPort = Integer.parseInt(args[0]);
	    	} catch(Exception e) {
	    		System.err.println(args[0]+" is an invalid port number!");
	    		System.exit(1);    		
	    	}
    	}
    	**/

		try {
			new ClothoCore();
		} catch(Exception e) {
			e.printStackTrace();
		}
	}

}
