package org.clothocad.core;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import java.util.Properties;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.webserver.jetty.ClothoWebserver;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter
        implements Daemon {

    private static final Logger logger = LoggerFactory.getLogger(ClothoStarter.class);
    protected static ClothoWebserver server;

    public static void main(String[] args)
            throws Exception {

        CommandLine cmd = parseArgs(args);

        Injector injector = Guice.createInjector(new ClothoModule(commandToProperties(cmd)), new MongoDBModule());

        server = injector.getInstance(ClothoWebserver.class);
        
        server.start();
    }
    
    protected static CommandLine parseArgs(String[] args) throws ParseException{
        Options options = ConfigOption.getOptions();
        CommandLineParser parser = new PosixParser();
        
        return parser.parse(options, args);
    }

    protected static Properties commandToProperties(CommandLine cmd) {
        Properties properties = new Properties();
        for (ConfigOption configOption : ConfigOption.values()){
            if (cmd.hasOption(configOption.abbreviation))
                properties.setProperty(configOption.name(), cmd.getOptionValue(configOption.abbreviation));
        }
        return properties;
    }
    protected DaemonContext context;

    @Override
    public void init(DaemonContext dc) {
        context = dc;
    }

    @Override
    public void start() throws Exception {
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