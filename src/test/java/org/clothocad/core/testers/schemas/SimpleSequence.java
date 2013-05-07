/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.schemas;

import javax.validation.constraints.Pattern;
import lombok.Getter;
import lombok.Setter;

/**
 *
 * @author spaige
 */
public class SimpleSequence {
    
    public SimpleSequence(String name, String sequence){
        this.name = name;
        this.sequence = sequence;
    }
    
    @Getter
    String name;
    
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    String sequence;
}
