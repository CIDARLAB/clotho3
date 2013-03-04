/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.aspects.Router;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.ContextHandler;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.clothocad.core.settings.Settings;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;

/**
 * This class exists merely to keep `Router.java` short and sweet
 * Note that nobody keeps permanent references to this class,
 * so it will be GC'ed soon after instantiation.
 *
 * Apache Tomcat project provides Server Side Includes (SSI) functionality
 * Jetty framework provides web server and WS/XHR services
 *      http://eclipse.org/jetty
 *
 * @author Kelvin Li
 */
class RouterServer {
    static void run() {
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

        Server server = new Server();
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
            createEchoServlets("./trunk/web/includes/");
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
