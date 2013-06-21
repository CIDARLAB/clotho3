package org.clothocad;

import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.broker.ClothoBroker;
import org.clothocad.core.ClothoCore;
import org.clothocad.webserver.jetty.ClothoWebserver;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter implements Daemon {
    
	public static void main(String[] args) 
			throws Exception {
		

		if(args.length == 1) {
			nPort = Integer.parseInt(args[0]);
		}
		
		//LogManager.getLogManager().reset();
		
            System.out.println("Ernst, ClothoBroker crashes on my machine, but if I silence this line I can run everything else fine");
           /**
            Exception in thread "main" javax.jms.JMSException: Could not connect: Virtual host stopped
            at org.fusesource.stomp.jms.StompJmsExceptionSupport.create(StompJmsExceptionSupport.java:59)
            */
                    
		// start the message broker
		new ClothoBroker();
		
		// wait a bit until the broker is running
		Thread.sleep(4000);
		
		// now, start the core (i.e. the server)
		new ClothoCore().start();
		
		// wait a bit until the broker is running
		Thread.sleep(4000);
		
		// start the Jetty webserver
		new ClothoWebserver(nPort);
		//System.out.println("The Clotho Webserver is running...");
		
		//Object lock = new Object();
        //synchronized (lock) {
        //    lock.wait();
        //}

	}
        
        
    private DaemonContext context;
        
    @Override
    public void init(DaemonContext dc){
        context = dc;
    }

    @Override
    public void start() throws Exception{
        System.out.println("starting ...");
        main(context.getArguments());
    }

    @Override
    public void stop() throws Exception {
        System.out.println("stopping ...");
    }

    @Override
    public void destroy() {
        System.out.println("done.");
    }
}
