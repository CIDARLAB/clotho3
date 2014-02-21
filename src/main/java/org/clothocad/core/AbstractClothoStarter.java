package org.clothocad.core;

import com.google.inject.Injector;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.clothocad.webserver.jetty.ClothoWebserver;

@Slf4j
abstract public class AbstractClothoStarter implements Daemon {
    /** This is part of a mechanism to factor out common code
      * across multiple main classes.
      */
    protected static interface MainHook {
        Injector getInjector(CommandLine cmd);
        void call(Injector injector);
    }

    private static ClothoWebserver server;
    protected DaemonContext context;

    /** Generic program entry point.
      * Customized setup is passed in the hook parameter.
      */
    protected static void
    baseMain(String[] args, MainHook hook) throws Exception {
        final CommandLine cmd;
        try {
            cmd = parseArgs(args);
        } catch (ParseException e){
            //TODO: customise message to include default values
            System.out.println(e.getMessage());
            printHelp();
            return;
        }
        if (cmd.hasOption("help")){
            printHelp();
            return;
        }
        /* Do custom setup */
        /* TODO: if keystorepass option passed without arg,
         * then prompt for password
         */
        Injector injector = hook.getInjector(cmd);
        hook.call(injector);

        server = injector.getInstance(ClothoWebserver.class);
        server.start();
    }

    private static void printHelp() {
        new HelpFormatter().printHelp(" ", ConfigOption.getOptions());
    }

    static CommandLine parseArgs(String[] args) throws ParseException {
        return new PosixParser().parse(ConfigOption.getOptions(), args);
    }

    /** Construct a Properties object from a properties file and command line
      * arguments. Command line arguments take precedence over the properties
      * file.
      */
    protected static Properties commandToProperties(CommandLine cmd) {
        Properties properties = new Properties();

        if (cmd.hasOption(ConfigOption.propfile.name())) {
            String path = cmd.getOptionValue(ConfigOption.propfile.name());
            try {
                properties.load(Files.newInputStream(Paths.get(path)));
            } catch (IOException ex) {
                log.warn("Could not load properties file at {}", path);
            }
            log.info("Properties loaded from file at {}", path);
        } else {
            log.debug("No property file specified.");
        }

        for (ConfigOption configOption : ConfigOption.values()){
            if (configOption.equals(ConfigOption.propfile))
                continue;
            String key = configOption.name();
            String value = cmd.getOptionValue(key);
            if (value != null)
                properties.setProperty(key, value);
        }
        return properties;
    }

    @Override
    public void init(DaemonContext dc) {
        context = dc;
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
