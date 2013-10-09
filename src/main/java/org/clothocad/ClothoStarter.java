package org.clothocad;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import java.util.Properties;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.core.ClothoModule;
import org.clothocad.core.ConfigOption;
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

        Options options = ConfigOption.getOptions();
        CommandLineParser parser = new PosixParser();
        
        CommandLine cmd = parser.parse(options, args);
        

        Injector injector = Guice.createInjector(new ClothoModule(commandToProperties(cmd)), new MongoDBModule());

        server = injector.getInstance(ClothoWebserver.class);
        
        server.start();
    }

    private static Properties commandToProperties(CommandLine cmd) {
        Properties properties = new Properties();
        for (ConfigOption configOption : ConfigOption.values()){
            if (cmd.hasOption(configOption.abbreviation))
                properties.setProperty(configOption.abbreviation, cmd.getOptionValue(null));
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
