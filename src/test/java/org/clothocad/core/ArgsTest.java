package org.clothocad.core;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import org.apache.commons.cli.ParseException;
import org.junit.Assert;
import org.junit.Test;

/**
 *
 * @author spaige
 */
public class ArgsTest {
    private static ClothoModule
    moduleFromArgs(String[] args) throws ParseException {
        return new ClothoModule(
            AbstractClothoStarter.commandToProperties(
                AbstractClothoStarter.parseArgs(args)));
    }
    
    @Test
    public void testDefaults() throws ParseException{
        ClothoModule clothoModule = moduleFromArgs(new String[]{});
        Properties properties = clothoModule.properties;
        for (ConfigOption configOption : ConfigOption.values()) {
            if (configOption.equals(ConfigOption.propfile)) continue;
            Assert.assertEquals(configOption.defaultValue, clothoModule.properties.getProperty(configOption.name()));
        }
    }
   
    @Test 
    public void testArgs() throws ParseException{
        Map<String, String> values = new HashMap<>();
        values.put("port", "8181");
        values.put("confidentialport", "8444");
        values.put("dbname", "customname");
        values.put("dbhost", "customhost");
        values.put("dbport", "1");
        values.put("keystorepath", "kpath");
        
        List<String> args = new ArrayList<>();
        for (String key : values.keySet()){
            args.add("--"+key);
            args.add(values.get(key));
        }
        
        ClothoModule clothoModule = moduleFromArgs(args.toArray(new String[0]));
        for (String key : values.keySet()){
            Assert.assertEquals(values.get(key), clothoModule.properties.getProperty(key));
        }
    }
    
    @Test
    public void testPropFile() throws ParseException {
        List<String> args = new ArrayList<>();
        args.add("--propfile");
        args.add(Paths.get("src", "test", "resources", "testproperties.properties").toString());
        ClothoModule clothoModule = moduleFromArgs(args.toArray(new String[0]));
        Assert.assertEquals("loaded from property file", clothoModule.properties.getProperty("dbname"));
    }
}
