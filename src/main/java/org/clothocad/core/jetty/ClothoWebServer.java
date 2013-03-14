package org.clothocad.core.jetty;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.clothocad.core.settings.Settings;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

public class ClothoWebServer {
	public static void main(String[] args) 
			throws Exception {
		
        HandlerList handlers = new HandlerList();
        ServletContextHandler context0 = createContext();
        if (context0 == null) {
            Logger.log(Logger.Level.FATAL, "Cannot create context");
            return;
        }
        handlers.setHandlers(new Handler[] {context0});

        SelectChannelConnector connector0 = new SelectChannelConnector();
        connector0.setHost(Settings.getHost());
        connector0.setPort(Settings.getPort());

		Server server = new Server(8080);
        server.setConnectors(new Connector[] {connector0});
        server.setHandler(handlers);
        try {
            server.start();
        } catch (Exception e) {
            Logger.log(Logger.Level.FATAL, "Cannot start", e);
            return;
        }
        Logger.log(Logger.Level.INFO, "Running");
    }

    private static ServletContextHandler createContext() {
        ServletContextHandler context =
            new ServletContextHandler(ServletContextHandler.SESSIONS);
        context.setContextPath("/");

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
}
