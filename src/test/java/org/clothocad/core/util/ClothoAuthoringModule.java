package org.clothocad.core.util;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.dataauthoring.FileHookPersistor;
import org.clothocad.core.testers.ClothoTestModule;

import com.google.inject.name.Names;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;

/**
 *
 * @author spaige
 * @author billcao
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
        super.configure();
        bind(Path.class).annotatedWith(Names.named("storagefolder")).toInstance(Paths.get("src", "test", "resources", "authoredJSON"));
        bind(Path.class).annotatedWith(Names.named("resourcefolder")).toInstance(Paths.get("src", "test", "resources", "authoredResources"));
        bind(Persistor.class).to(FileHookPersistor.class);
    }

}
