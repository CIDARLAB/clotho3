package org.clothocad.broker;

import org.apache.activemq.apollo.broker.Broker;
import org.apache.activemq.apollo.broker.store.leveldb.dto.*;
import org.apache.activemq.apollo.dto.*;
import org.iq80.leveldb.util.FileUtils;

import java.io.File;

public class ClothoBroker {

	private Broker broker = null;
	
	public ClothoBroker() {
        //
        // Creating and initially configuring the broker.
        this.broker = new Broker();
        
        this.broker.setTmp(new File("./tmp"));
        this.broker.setConfig(createConfig());

        this.broker.start(new Runnable() {
            public void run() {
            	System.out.println("The Broker is running...");
            }
        });
        
        Runtime.getRuntime().addShutdownHook(new Thread() {
            public void run() {
            	// delete the data directory
            	FileUtils.deleteDirectoryContents(new File("./data"));
            }
        });
	}
	
	public void stop() {
        // The broker stops asynchronously. The runnable is invoked once
        // the broker if fully stopped.
        broker.stop(new Runnable(){
            public void run() {
                System.out.println("The broker has now stopped.");
            }
        });
	}
	
    /**
     * Builds a simple configuration model with just plain Java.  Corresponds 1 to 1 with
     * the XML configuration model.  See the Apollo user guide for more details.
     * @return
     */
    private BrokerDTO createConfig() {
        BrokerDTO broker = new BrokerDTO();

        // Brokers support multiple virtual hosts.
        VirtualHostDTO host = new VirtualHostDTO();
        host.id = "localhost";
        host.host_names.add("localhost");
        host.host_names.add("127.0.0.1");

        // The message store is configured on the virtual host.
        LevelDBStoreDTO store = new LevelDBStoreDTO();
        store.directory = new File("./data");
        host.store = store;

        broker.virtual_hosts.add(host);

        /**
  <connector id="tcp" bind="tcp://0.0.0.0:61613" connection_limit="2000"/>
  <connector id="tls" bind="tls://0.0.0.0:61614" connection_limit="2000"/>
  <connector id="ws"  bind="ws://0.0.0.0:61623"  connection_limit="2000"/>
  <connector id="wss" bind="wss://0.0.0.0:61624" connection_limit="2000"/>
         **/
        
        // Control which ports and protocols the broker binds and accepts
        AcceptingConnectorDTO tcpConnector = new AcceptingConnectorDTO();
        tcpConnector.id = "tcp";
        tcpConnector.bind = "tcp://0.0.0.0:61613";
        broker.connectors.add(tcpConnector);

        AcceptingConnectorDTO tlsConnector = new AcceptingConnectorDTO();
        tlsConnector.id = "tls";
        tlsConnector.bind = "tls://0.0.0.0:61614";
        broker.connectors.add(tlsConnector);

        AcceptingConnectorDTO wsConnector = new AcceptingConnectorDTO();
        wsConnector.id = "ws";
        wsConnector.bind = "ws://0.0.0.0:61623";
        broker.connectors.add(wsConnector);
        
        AcceptingConnectorDTO wssConnector = new AcceptingConnectorDTO();
        wssConnector.id = "wss";
        wssConnector.bind = "wss://0.0.0.0:61624";
        broker.connectors.add(wssConnector);

        //
        // Fires up the web admin console on HTTP.
        //WebAdminDTO webadmin = new WebAdminDTO();
        //webadmin.bind = "http://0.0.0.0:8080";
        //broker.web_admins.add(webadmin);

        return broker;
    }

	public static void main(String[] args) 
			 throws Exception {

		ClothoBroker broker = new ClothoBroker();
		
		Object lock = new Object();
        synchronized (lock) {
            lock.wait();
        }
        
	}
}
