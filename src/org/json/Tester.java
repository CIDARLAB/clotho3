
package org.json;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author John Christopher Anderson
 */


public class Tester {
    public static void main(String[] args) {
        try {
            JSONElement json = new JSONObject("{\"email\":\"cindy@su.com\",\"name\":\"cindysu\"}");
            System.out.println(json.getEscaped());
            System.out.println(json.toString());
        } catch (JSONException ex) {
            Logger.getLogger(Tester.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}


