/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import javax.validation.constraints.Pattern;
import lombok.Getter;
import lombok.Setter;

/**
 *
 * @author prashantvaidyanathan
 */
public class TestSeq {
    
    public TestSeq(String sequence){
        this.sequence = sequence;
    }
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    String sequence;
    
    
}
