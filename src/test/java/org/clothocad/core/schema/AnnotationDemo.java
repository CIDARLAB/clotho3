/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.io.PrintWriter;
import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;
import javax.validation.constraints.Size;
import org.clothocad.core.schema.Schema;
import org.objectweb.asm.ClassReader;
import org.objectweb.asm.util.ASMifier;
import org.objectweb.asm.util.TraceClassVisitor;

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
        ASMifier.main(new String[]{"org.clothocad.core.testers.schemas.AnnotationDemo"});
                ClothoSchemaTest.setUpClass();
        Schema feature = ClothoSchemaTest.createFeatureSchema();
        ClassReader cr = new ClassReader(feature.getClassData());

                cr.accept(new TraceClassVisitor(null,
                new ASMifier(),
                new PrintWriter(System.out)), 0);
    }
}
