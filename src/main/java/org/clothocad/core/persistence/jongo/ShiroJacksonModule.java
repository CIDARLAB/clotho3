package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.databind.module.SimpleModule;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.SimpleByteSource;

/**
 *
 * @author spaige
 */
public class ShiroJacksonModule extends SimpleModule {

    public ShiroJacksonModule() {
        super("ShiroModule", new Version(0,0,1, null, "org.clothocad", "clotho"));
    }
    
    @Override
    public void setupModule(SetupContext context) {
        context.setMixInAnnotations(SimpleHash.class, SimpleHashMixin.class);
        context.setMixInAnnotations(SimpleByteSource.class, SimpleByteSourceMixin.class);
        
        super.setupModule(context);
    }

    static abstract class SimpleHashMixin {
        
        @JsonCreator
        public SimpleHashMixin(@JsonProperty("algorithmName") String algorithmName){}
    }

    static abstract class SimpleByteSourceMixin {

        @JsonCreator
        public SimpleByteSourceMixin(@JsonProperty("bytes") byte[] bytes) {}
    }
}
