package org.clothocad;

import java.util.Arrays;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.broker.ClothoBroker;
import org.clothocad.core.ClothoCore;
import org.clothocad.webserver.jetty.ClothoWebserver;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter
        implements Daemon {
    
    private static final Logger logger = LoggerFactory.getLogger(ClothoStarter.class);
    private static ClothoWebserver server;
    private static ClothoBroker broker;
    private static ClothoCore core;
    
	public static void main(String[] args) 
			throws Exception {
		int nPort = 9090;

		if(args.length == 1) {
			nPort = Integer.parseInt(args[0]);
		}
		
                try{
                    broker = new ClothoBroker();
                    // wait a bit until the broker is running
                    Thread.sleep(4000);
		
                    // now, start the core (i.e. the server)
                    core = new ClothoCore();
                    core.start();
                } catch (Exception e){
                    logger.error("Exception while initializing JMS Broker and Listener", e);
                }

		// start the Jetty webserver
                logger.debug("Starting server on port {}", nPort);
		server = new ClothoWebserver(nPort);

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
        System.out.println("stopping ...");
        if (core != null){
            core.shutdown();
        }
        
        if (broker != null){
            broker.stop();
        }
        
        server.getServer().stop();
    }

    @Override
    public void destroy() {
        System.out.println("done.");
    }
}

