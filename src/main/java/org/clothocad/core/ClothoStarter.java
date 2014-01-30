package org.clothocad.core;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Properties;
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.util.JSON;
import org.clothocad.webserver.jetty.ClothoWebserver;

//Start then navigate to:  http://localhost:8080/#/
@Slf4j
public class ClothoStarter
        implements Daemon {

    protected static ClothoWebserver server;

    public static void main(String[] args)
            throws Exception {

        try {
            CommandLine cmd = parseArgs(args);
            
            if (cmd.hasOption("help")){
                printHelp();
                return;
            }
            //TODO: if keystorepass option passed w/o arg, prompt for password 
            Injector injector = Guice.createInjector(new ClothoModule(commandToProperties(cmd)), new JongoModule());

            Persistor persistor = injector.getInstance(Persistor.class);
        
            ensureMinimalObjects(persistor);
        
            server = injector.getInstance(ClothoWebserver.class);
        
            server.start();
        } catch (ParseException e){
            //TODO: customise message to include default values
            System.out.println(e.getMessage());
            printHelp();
        }
    }
    
    protected static void printHelp(){
        new HelpFormatter().printHelp(" ", ConfigOption.getOptions());
    }
    
    protected static void ensureMinimalObjects(Persistor p){
        JSON.importTestJSON(Paths.get("src", "main", "resources", "json", "essential").toString(), p, false);
    }
    
    protected static CommandLine parseArgs(String[] args) throws ParseException{
        Options options = ConfigOption.getOptions();
        CommandLineParser parser = new PosixParser();
        
        return parser.parse(options, args);
    }

    protected static Properties commandToProperties(CommandLine cmd) {
        Properties properties = new Properties();
        
        //properties in properties file are overwritten by command-line arguments
        if (cmd.hasOption(ConfigOption.propfile.name())){
            String path = cmd.getOptionValue(ConfigOption.propfile.name());
            try {
                properties.load(Files.newInputStream(Paths.get(path)));
                log.info("Properties loaded from file at {}", path);
            } catch (IOException ex) {
                log.warn("Could not load properties file at {}", path);
            }
        } else {
            log.debug("No property file specified.");
        }
        
        for (ConfigOption configOption : ConfigOption.values()){
            if (configOption.equals(ConfigOption.propfile)) continue;
            if (cmd.hasOption(configOption.name()))
                properties.setProperty(configOption.name(), cmd.getOptionValue(configOption.name()));
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
