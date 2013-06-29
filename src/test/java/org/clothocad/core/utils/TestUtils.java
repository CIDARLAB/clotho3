/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.utils;

import com.google.inject.Guice;
import com.google.inject.Injector;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.testers.ClothoTestModule;

/**
 *
 * @author spaige
 */
public class TestUtils {
    //not sure if this should be static 
    
    private static Injector injector;

    static {
        Injector injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());

    }

    public static <T> T getA(Class<T> type) {
        return injector.getInstance(type);
    }
}
