/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import static com.fasterxml.jackson.databind.SerializationFeature.FAIL_ON_EMPTY_BEANS;
import org.jongo.Mapper;
import org.jongo.ObjectIdUpdater;
import org.jongo.marshall.Marshaller;
import org.jongo.marshall.jackson.JacksonIdFieldSelector;
import org.jongo.marshall.jackson.configuration.AbstractMappingBuilder;
import org.jongo.marshall.jackson.configuration.MapperModifier;
import org.jongo.marshall.jackson.configuration.Mapping;
import org.jongo.query.BsonQueryFactory;
import org.jongo.query.QueryFactory;

/**
 *
 * @author spaige
 */
public class ClothoMapper implements Mapper {
    
    private final RefResolvingJacksonEngine engine;
    private final ObjectIdUpdater objectIdUpdater;
    private final QueryFactory queryFactory;

    public ClothoMapper(){
        //behold my mistreatment of the builder pattern and weep
        ClothoMapperBuilder builder = new ClothoMapperBuilder();
        Mapping mapping = builder.createMapping();
        engine = new RefResolvingJacksonEngine(mapping);
        objectIdUpdater = new ClothoObjectIdUpdater(new JacksonIdFieldSelector());
        queryFactory = new BsonQueryFactory(engine);
    }

    @Override
    public Marshaller getMarshaller() {
        return engine;
    }

    @Override
    public ExtendedUnmarshaller getUnmarshaller() {
        return engine;
    }

    @Override
    public ObjectIdUpdater getObjectIdUpdater() {
        return objectIdUpdater;
    }

    @Override
    public QueryFactory getQueryFactory() {
        return queryFactory;
    }
    
    //exposing some protected utility methods
    public class ClothoMapperBuilder extends AbstractMappingBuilder<ClothoMapperBuilder>{

        public ClothoMapperBuilder() {
            super();
        }
        
        public Mapping createMapping(){
            //add customizations
            addModifier(new MapperModifier() {
                @Override
                public void modify(ObjectMapper mapper) {
                        mapper.disable(FAIL_ON_EMPTY_BEANS);
                        //write types into serialized objects
                        mapper.enableDefaultTyping(ObjectMapper.DefaultTyping.OBJECT_AND_NON_CONCRETE);
                        //mapper.setDefaultTyping(typer);
                        
                        // databind module handles ids and circular structures
                        //mapper.registerModule(new ClothoDatabindModule());
                    }
                });
            return super.createMapping();
        }

        @Override
        protected ClothoMapperBuilder getBuilderInstance() {
            return this;
        }
    }
}
