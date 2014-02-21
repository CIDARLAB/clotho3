package org.clothocad.core;

import java.nio.file.Paths;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 *
 * @author spaige
 */
public enum ConfigOption {
    port("HTTP listening port", "8080", "port"),
    confidentialport("HTTPS listening port", "8443", "port"),
    dbname("database name", "clotho", "name"),
    dbhost("database url hostname", "localhost", "hostname"),
    dbport("database url port", "27017", "port"),
    //loglevel,
    keystorepath("SSL keystore path",
                 Paths.get(System.getProperty("java.home"))
                      .resolve(Paths.get("lib", "security", "cacerts"))
                      .toString(),
                 "path"),
    keystorepass("SSL keystore password", "", "password"),
    propfile("path to properties file", "", "path"),
    clientdirectory("path to client files directory",
                    Paths.get("clotho3-web", "dist").toString(),
                    "path");
    
    final String description;
    final String defaultValue;
    final String argName;
    private ConfigOption(final String description,
                         final String defaultValue,
                         final String argName) {
        this.description = description;
        this.defaultValue = defaultValue;
        this.argName = argName;
    }

    private Option
    toOption() {
        final String desc =
            String.format("%s (default: %s)", description, defaultValue);
        final Option out = new Option(this.name(), true, desc);
        out.setArgName(argName);
        return out;
    }

    public static Options getOptions() {
        Options options = new Options();
        for (ConfigOption opt : ConfigOption.values()) {
            options.addOption(opt.toOption());
        }
        options.addOption(new Option("help", false, "print this message"));
        return options;
    }
}
