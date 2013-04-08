/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.dm.clotho2;

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public class Institution extends ObjBase {

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
}