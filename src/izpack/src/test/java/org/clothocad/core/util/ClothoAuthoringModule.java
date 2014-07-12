/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.name.Names;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.dataauthoring.FileHookPersistor;
import org.clothocad.core.testers.ClothoTestModule;

/**
 *
 * @author spaige
 */
class ClothoAuthoringModule extends ClothoTestModule {
    public ClothoAuthoringModule(Properties config) {
        super(config);
    }
    
    public ClothoAuthoringModule(){
        this(null);
    }

    @Override
    protected void configure() {
        super.configure(); //To change body of generated methods, choose Tools | Templates.
        bind(Path.class).annotatedWith(Names.named("storagefolder")).toInstance(Paths.get("src", "test", "resources", "authoredJSON"));
        bind(Persistor.class).to(FileHookPersistor.class);
    }
    
    
}
