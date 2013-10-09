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
    
    port("port", "the port the Clotho accepts connections on", "8080"),
    confidentialport("cport", "the port Clotho accepts ssl connections on", "8443"),
    dbname("dbn", "the name of the mongodb database to use", "clotho"),
    dbhost("dbh", "the hostname of the machine the database is on", "localhost"),
    dbport("dbp", "mongodb port", "27017"),
    //loglevel,
    keystorepath("ks", "path to the keystore", System.getProperty("java.home") + "/lib/security/cacerts".replace('/', File.separatorChar)),
    keystorepass("kspass", "keystore password", ""),
    keymanagerpass("ksmpass", "keymanager password", "");
    
    public final String abbreviation;
    public final String description;
    public final String defaultValue;
    public final boolean hasarg;
    private ConfigOption(final String abbreviation, final String description, final String defaultValue){
        hasarg = true;
        this.abbreviation = abbreviation;
        this.description = description;
        this.defaultValue = defaultValue;
    }
    
    public static Options getOptions() {
        Options options = new Options();
        for (ConfigOption option : ConfigOption.values()){
            options.addOption(new Option(option.abbreviation, option.name(), option.hasarg, option.description));
        }
        return options;
    }
}
