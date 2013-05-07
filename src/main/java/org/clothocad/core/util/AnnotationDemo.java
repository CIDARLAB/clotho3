/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;
import javax.validation.constraints.Size;
import org.objectweb.asm.util.ASMifier;

/**
 *
 * @author spaige
 */
public class AnnotationDemo {
    
    @Pattern(regexp="[A-Z]", flags={Pattern.Flag.CASE_INSENSITIVE})
    protected String pattern;
    
    @NotNull
    @Size(min=1, max=16)
    private String firstname;
    
    public static void main(String[] args) throws Exception{
        ASMifier.main(new String[]{"org.clothocad.core.util.AnnotationDemo"});
    }
}
