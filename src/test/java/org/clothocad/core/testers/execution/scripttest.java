/* Simple javax.script.ScriptEngine demonstration
 *
 * Expects the file 'scripttest.js' to be in the current working
 * directory of the running Java process.
 */
package org.clothocad.core.testers.execution;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class scripttest {
    static final Logger logger = LoggerFactory.getLogger(scripttest.class);
    public static void main(String[] args) {
        try {
            test();
        } catch (Exception e) {
            logger.error("",e);
        }
    }

    private static void test() throws ScriptException,
                                      FileNotFoundException {
        ScriptEngineManager sem = new ScriptEngineManager();
        ScriptEngine engine = sem.getEngineByName("JavaScript");

        /* Construct input */
        Map<String, Object> j_spam = new HashMap<>();
        j_spam.put("name", "My Spam");
        j_spam.put("age", 10);
        engine.put("spam", j_spam);

        /* Construct skeleton output (dynamically modified by script) */
        Map<String, Object> j_eggs = new HashMap<>();
        engine.put("eggs", j_eggs);

        /* Feed file handle of script into engine and execute */
        File f = new File("test/org/clothocad/core/scripting/scripttest.js");
        engine.eval(new FileReader(f));

        /* Get output */
        String j_eggs_name = j_eggs.get("name").toString();
        Double j_eggs_age = Double.parseDouble(j_eggs.get("age").toString());

        /* Print output */
        logger.info(
                   "j_eggs is called '{}'", j_eggs_name);
        logger.info(
                   "j_eggs is {} days old", String.valueOf(j_eggs_age));
    }
}
