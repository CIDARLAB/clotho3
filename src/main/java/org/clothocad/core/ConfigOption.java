/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import java.io.File;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 *
 * @author spaige
 */
public enum ConfigOption {
    port("http listening port", "8080", "port"),
    confidentialport("https listening port", "8443", "port"),
    dbname("database name", "clotho", "name"),
    dbhost("database url hostname", "localhost", "hostname"),
    dbport("database url port", "27017", "port"),
    //loglevel,
    keystorepath("ssl keystore path", System.getProperty("java.home") + "/lib/security/cacerts".replace('/', File.separatorChar), "path"),
    keystorepass( "ssl keystore password", "", "password"),
    propfile( "path to properties file", "", "path"),
    clientdirectory( "path to client files directory", "clotho3-web/dist", "path");
    
    public final String description;
    public final String defaultValue;
    public final String argName;
    public final boolean hasarg;
    private ConfigOption(final String description, final String defaultValue, final String argName){
        hasarg = true;
        this.description = description;
        this.defaultValue = defaultValue;
        this.argName = argName;
    }
    
    public static Options getOptions() {
        Options options = new Options();
        for (ConfigOption option : ConfigOption.values()){
            Option cliOption = new Option(option.name(), option.hasarg, option.description + String.format(" (default: %s)", option.defaultValue));
            cliOption.setArgName(option.argName);
            options.addOption(cliOption);
        }
        options.addOption(new Option("help", false, "print this message"));
        return options;
    }
}
