package org.clothocad.core.util;

import static org.junit.Assert.assertEquals;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import org.apache.commons.cli.ParseException;
import org.clothocad.core.ConfigOption;
import org.junit.Test;

/**
 *
 * @author spaige
 */
public class ConfigTest {
    @Test
    public void testPrecedence() throws ParseException {
        final Properties subject = Config.get(Config.parseArgs(new String[] {
            "--configfile",
                Paths.get("src", "test", "resources", "clothoconfig")
                     .toString(),
            "--port", "123",
            "--port", "9876",
            "--dbhost", "goodhost"
        }));
        assertEquals("8443", subject.getProperty("confidentialport"));
        assertEquals("config file DB", subject.getProperty("dbname"));
        assertEquals("9876", subject.getProperty("port"));
        assertEquals("goodhost", subject.getProperty("dbhost"));
    }
}
