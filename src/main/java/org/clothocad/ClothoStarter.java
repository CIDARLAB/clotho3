package org.clothocad;

import com.google.inject.Guice;
import com.google.inject.Injector;
import org.clothocad.broker.ClothoBroker;
import org.clothocad.core.ClothoCore;
import org.clothocad.core.ClothoModule;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.webserver.jetty.ClothoWebserver;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter {

    public static void main(String[] args)
            throws Exception {

        int nPort = 8080;
        if (args.length == 1) {
            nPort = Integer.parseInt(args[0]);
        }

        //LogManager.getLogManager().reset();

        System.out.println("Ernst, ClothoBroker crashes on my machine, but if I silence this line I can run everything else fine");
        /**
         * Exception in thread "main" javax.jms.JMSException: Could not connect:
         * Virtual host stopped at
         * org.fusesource.stomp.jms.StompJmsExceptionSupport.create(StompJmsExceptionSupport.java:59)
         */
        /*	// start the message broker
         new ClothoBroker();
		
         // wait a bit until the broker is running
         Thread.sleep(4000);
		
         // now, start the core (i.e. the server)
         new ClothoCore().start();
		
         // wait a bit until the broker is running
         Thread.sleep(4000);*/
        // start the Jetty webserver
        Injector injector = Guice.createInjector(new ClothoModule(), new MongoDBModule());

        new ClothoWebserver(nPort);
        //System.out.println("The Clotho Webserver is running...");

        //Object lock = new Object();
        //synchronized (lock) {
        //    lock.wait();
        //}

    }
}
