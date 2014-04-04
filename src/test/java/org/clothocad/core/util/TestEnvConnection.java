/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.net.UnknownHostException;
import javax.inject.Inject;
import javax.inject.Named;
import org.clothocad.core.persistence.mongodb.ClothoMapper;
import org.clothocad.core.persistence.mongodb.MongoDBConnection;

/**
 * A simple way to guarantee that the database will be flushed completely
 * in the test environment before startup
 * 
 * @author spaige
 */
public class TestEnvConnection extends MongoDBConnection {

    @Inject
    public TestEnvConnection(@Named("dbport") int port, @Named("dbhost")String host, @Named("dbname") String dbName, ClothoMapper mapper) {
        super(port, host, dbName, mapper);
        try {
            this.connect();
        } catch (UnknownHostException ex) {
            throw new RuntimeException(ex);
        }
        this.deleteAll();
    }
    
}
