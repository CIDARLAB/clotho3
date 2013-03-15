package org.clothocad.core.jetty;

import org.clothocad.core.settings.Settings;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.bio.SocketConnector;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;

public class ClothoWebServer {
	
	private boolean running = false;
	
	public ClothoWebServer() {
		running = false;
	}
	
	public boolean isRunning() {
		return running;
	}
	
	public void start() 
			throws Exception {
        
		
		/**
        ServletContextHandler context0 = createContext();
        if (context0 == null) {
            Logger.log(Logger.Level.FATAL, "Cannot create context");
            return;
        }
        handlers.setHandlers(new Handler[] {context0});
		**/
        
		SocketConnector coreConnector = new SocketConnector();
        coreConnector.setHost(Settings.getHost());
        coreConnector.setPort(Settings.getPort());

        ResourceHandler resourceHandler = new ResourceHandler();
        resourceHandler.setDirectoriesListed(false);
        resourceHandler.setWelcomeFiles(new String[]{ "login.html" }); 
        resourceHandler.setResourceBase("./web/");

        HandlerList handlers = new HandlerList();
        handlers.setHandlers(new Handler[] { resourceHandler, new RESTHandler() });        

        // now, start the server
		Server server = new Server();
        server.setConnectors(new Connector[] {coreConnector});
        server.setHandler(handlers);
        server.start();
        
        System.out.println("Clotho's web server is running on 8080...");
    }

	/***
    private static ServletContextHandler createContext() {
        ServletContextHandler context =
            new ServletContextHandler(ServletContextHandler.SESSIONS);
        
        
        context.setContextPath("/web");

        context.addServlet(new ServletHolder(WSServlet.get()),
                           "/servlet/websocket");
        context.addServlet(new ServletHolder(XHRServlet.get()),
                           "/servlet/xhr");
        context.addServlet(new ServletHolder(new PageServlet()),
                           "/servlet/page");
        context.addServlet(new ServletHolder(new StaticServlet()), "/");

        List<EchoServlet> echo_servlets =
            createEchoServlets("./web/includes/");
        if (echo_servlets == null)
            return null;
        for (EchoServlet s : echo_servlets) {
            context.addServlet(new ServletHolder(s),
                               "/servlet/include/" + s.name());
        }
        return context;
    }
	***/
	
	/***
    private static List<EchoServlet> createEchoServlets(String path) {
        File includes_dir = new File(path);
        File[] includes_files = includes_dir.listFiles();
        if (includes_files == null) {
            Logger.log(Logger.Level.FATAL,
                       "Cannot find directory " + path);
            return null;
        }

        List<EchoServlet> out = new ArrayList<EchoServlet>();
        for (File f : includes_files) {
            if (!f.isFile())
                continue;
            FileInputStream fis;
            try {
                fis = new FileInputStream(f);
            } catch (FileNotFoundException e) {
                continue;
            }

            byte[] contents;
            try {
                contents = FileUtils.dumpInputStream(fis);
            } catch (IOException e) {
                continue;
            }
            out.add(new EchoServlet(f.getName(), new String(contents)));
        }
        return out;
    }
    ***/
	
}
