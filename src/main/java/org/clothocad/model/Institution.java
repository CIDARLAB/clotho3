/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.JCAClothoSchema;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *org.clothocad.model.Institution
 * @author jcanderson
 */
@NoArgsConstructor
public class Institution extends ObjBase implements JCAClothoSchema {

    /**
     * Constructor from raw data
     * @param name
     * @param city
     * @param state
     * @param country
     */
    
    @Setter
    @Getter
    private String city, state, country;
    
    //TODO:unique name constraint
    public Institution( String name, String city, String state, String country ) {
        super(name);
        this.city = city;
        this.state = state;
        this.country = country;
    }


    /** drag/drop handler
    public boolean addObject( ObjBase dropObject ) {
        switch ( dropObject.getType() ) {
            case LAB:
                Lab theLab = (Lab) dropObject;
                theLab.changeInstitition( this );
                return true;
            default:
                return false;
        }
    }*/

    /* GETTERS
     * */

    public static Institution retrieveByName( String name ) {
        throw new UnsupportedOperationException("not implemented yet.");
    }

    @Override
    public boolean validate(JSONObject obj) {
        //Stephanie:  I believe this can all be expressed with @notNull type things, so this is redundant
        //This validate method covers the per-object verification, of which there is none for Institution unless it confirms that the cities and states are existant entities
        try {
            if(obj.getString("name")==null || obj.getString("name").equals("")) {
                return false;
            }
            if(obj.getString("city")==null || obj.getString("city").equals("")) {
                return false;
            }
            if(obj.getString("state")==null || obj.getString("state").equals("")) {
                return false;
            }
            if(obj.getString("country")==null || obj.getString("country").equals("")) {
                return false;
            }
            return true;
        } catch (JSONException ex) {
            return false;
        }
    }
}
