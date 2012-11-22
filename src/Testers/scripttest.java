/* Simple javax.script.ScriptEngine demonstration
 *
 * Expects the file 'scripttest.js' to be in the current working
 * directory of the running Java process.
 */
package Testers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.json.JSONException;
import org.json.JSONObject;
import org.clothocad.core.aspects.Logger;

public class scripttest {
    public static void main(String[] args) {
        try {
            test();
        } catch (Exception e) {
            Logger.log(Logger.Level.FATAL, "", e);
        }
    }

    private static void test() throws ScriptException,
                                      JSONException,
                                      FileNotFoundException {
        ScriptEngineManager sem = new ScriptEngineManager();
        ScriptEngine engine = sem.getEngineByName("JavaScript");

        /* Construct input */
        JSONObject j_spam = new JSONObject();
        j_spam.put("name", "My Spam");
        j_spam.put("age", 10);
        engine.put("spam", j_spam);

        /* Construct skeleton output (dynamically modified by script) */
        JSONObject j_eggs = new JSONObject();
        engine.put("eggs", j_eggs);

        /* Feed file handle of script into engine and execute */
        File f = new File("scripttest.js");
        engine.eval(new FileReader(f));

        /* Get output */
        String j_eggs_name = j_eggs.getString("name");
        Double j_eggs_age = j_eggs.getDouble("age");

        /* Print output */
        Logger.log(Logger.Level.INFO,
                   "j_eggs is called '" + j_eggs_name + "'");
        Logger.log(Logger.Level.INFO,
                   "j_eggs is " + String.valueOf(j_eggs_age) + " days old");
    }
}
