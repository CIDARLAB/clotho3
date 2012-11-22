/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.json;

/**
 *
 * @author jcanderson
 */
public interface JSONElement {
    
    /**
     * This returns the string representation of a the JSONObject or JSONArray
     * properly escaped so it can be re-parsed into an object later
     * 
     * The implementation as of 7/1 just deals with quotes, not other special
     * characters that likely need to be added.  There's probably a library
     * somewhere that can be adapted for this.
     * @return 
     */
    public String getEscaped();
}
