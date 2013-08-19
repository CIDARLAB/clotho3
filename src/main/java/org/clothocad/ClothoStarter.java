package org.clothocad;

import java.util.Arrays;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.webserver.jetty.ClothoWebserver;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter
        implements Daemon {
    
    private static final Logger logger = LoggerFactory.getLogger(ClothoStarter.class);
    private static ClothoWebserver server;
    
	public static void main(String[] args) 
			throws Exception {
		
		int nPort = 8080;
		if(args.length == 1) {
			nPort = Integer.parseInt(args[0]);
		}
		
		// start the Jetty webserver
		new ClothoWebserver(nPort);
		logger.info("The Clotho Webserver is running...");
	}
        
    private DaemonContext context;
        
    @Override
    public void init(DaemonContext dc){
        context = dc;
    }

    @Override
    public void start() throws Exception{
        System.out.println("starting with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }

    @Override
    public void stop() throws Exception {
        
        server.getServer().stop();
    }

    @Override
    public void destroy() {
        System.out.println("done.");
    }
}
