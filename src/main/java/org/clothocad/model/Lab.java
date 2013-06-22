    /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import com.github.jmkgreen.morphia.annotations.Reference;
import javax.validation.constraints.NotNull;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.JCAClothoSchema;

import org.clothocad.core.datums.ObjBase;
import org.hibernate.validator.constraints.URL;
import org.json.JSONObject;


/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public class Lab extends ObjBase implements JCAClothoSchema {

    @Getter
    @Setter
    private String department, address;
          
    @Getter
    @Setter
    @NotNull
    private String website;
    @Getter
    @Setter
    @Reference
    private Person PI;
    @Getter
    @Setter
    @Reference
    private Institution institution;
    //TODO: validation: unique name criterion

    public Lab( Institution inst, Person PI, String name, String department, String address ) {
        super(name);
        this.department = department;
        this.address = address;
        institution = inst;
        this.PI = PI;
    }





    /*(@Override
    public boolean addObject( ObjBase dropObject ) {
        switch ( dropObject.getType() ) {
            case PERSON:
                //PUT THE PERSON IN this LAB, or change their affiliation
                return true;
            default:
                return false;
        }

    }*/

   


    public static Lab retrieveByName( String name ) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean validate(JSONObject obj) {
        return true;
    }

}
