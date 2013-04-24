/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.Set;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.model.Person;
import org.objectweb.asm.ClassWriter;
import org.objectweb.asm.FieldVisitor;
import static org.objectweb.asm.Opcodes.*;
import org.objectweb.asm.Type;

/**
 *
 * @author spaige
 */
public class ClothoSchema extends Schema {
    
    protected Set<ClothoField> fields;
    protected Set<Function> methods;
    protected Schema superClass;
    
    public ClothoSchema(String name, String description, Person author, Schema superClass, Set<ClothoField> fields){
        super(name, description, author);
        this.fields = fields;
        this.superClass = superClass;
    }

    @Override
    public Language getLanguage() {
        return Language.JSONSCHEMA;
    }

    @Override
    public void setSource(String source) {
        throw new UnsupportedOperationException("Cannot set source on JSON-based schemas"); //TODO: should be equivalent to setJSON?
    }
    
    @Override 
    public byte[] getClassData(){
        if (classData == null){
            classData = generateClassData();
        }
        return classData;
    }
    
    private byte[] generateClassData(){
        ClassWriter cw = new ClassWriter(0);
        String superClassName = this.superClass == null ? "org/clothocad/core/datums/ObjBase" : this.superClass.getInternalName();
        
        cw.visit(V1_7, ACC_PUBLIC, this.getInternalName(), null, superClassName, new String[]{});
        //store original name
        cw.visitField(ACC_PUBLIC + ACC_FINAL + ACC_STATIC, "SCHEMA_NAME", Type.getType(String.class).getInternalName(), null, this.getName()).visitEnd();
        
        //fields
        for (ClothoField field : fields){
            FieldVisitor fv = cw.visitField(field.getAccess(), field.getName(), Type.getType(field.getType()).getInternalName(), null, null);
            //TODO: annotating embed vs reference
            if (field.isReference()){
                fv.visitAnnotation(Type.getType(Reference.class).getInternalName(), true).visitEnd();
            }
            
            fv.visitEnd();
        }
        //TODO: methods

        
        /*for (Function method : methods){
            mv = cw.visitMethod(method.getAccess(), method.getName(), Type.getType(method.getReturnType()).getInternalName(), method.getSignature(), null)
            //if method name is isValid, attach validator annotation
            //also append all field-level validators
            //bind handle to function call
        }*/

        cw.visitEnd();
        return cw.toByteArray();
    }
    
}
