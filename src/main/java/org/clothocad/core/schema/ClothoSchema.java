/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.lang.reflect.Array;
import java.util.Set;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.model.Person;
import org.objectweb.asm.AnnotationVisitor;
import org.objectweb.asm.ClassWriter;
import org.objectweb.asm.FieldVisitor;
import org.objectweb.asm.MethodVisitor;
import static org.objectweb.asm.Opcodes.*;
import org.objectweb.asm.Type;

/**
 *
 * @author spaige
 */
public class ClothoSchema extends Schema {
    public ClothoSchema() {}

    
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
        cw.visitField(ACC_PUBLIC + ACC_FINAL + ACC_STATIC, "SCHEMA_NAME", Type.getType(String.class).getDescriptor(), null, this.getName()).visitEnd();
        
        //no-args constructor
        MethodVisitor constructorVisitor = cw.visitMethod(ACC_PUBLIC, "<init>", Type.getMethodDescriptor(Type.VOID_TYPE, new Type[] {}), null, null);
        constructorVisitor.visitCode();
        constructorVisitor.visitVarInsn(ALOAD, 0);
        constructorVisitor.visitMethodInsn(INVOKESPECIAL, superClassName, "<init>", Type.getMethodDescriptor(Type.VOID_TYPE, new Type[] {}));
        constructorVisitor.visitInsn(RETURN);
        constructorVisitor.visitMaxs(1,1);
        constructorVisitor.visitEnd();
        
        //fields
        for (ClothoField field : fields){
            FieldVisitor fv = cw.visitField(accessToOpcode(field.getAccess()), field.getName(), Type.getType(field.getType()).getDescriptor(), null, null);
            //TODO: annotating embed vs reference
            if (field.isReference()){
                fv.visitAnnotation(Type.getType(Reference.class).getInternalName(), true).visitEnd();
            }
            for (Constraint c : field.getConstraints()){
               AnnotationVisitor av =  fv.visitAnnotation(Type.getDescriptor(c.getAnnotation()), true);
               for (String value : c.getValues()){
                   handleAnnotationValue(av,value,c.getValue(value));
               }
               av.visitEnd();
            }
            fv.visitEnd();
            //getters and setters
            if (field.getAccess() != Access.PRIVATE){
                if (field.getAccess() != Access.READONLY){
                    //setter
                    MethodVisitor mv = cw.visitMethod(ACC_PUBLIC, field.getSetterName(), Type.getMethodDescriptor(Type.VOID_TYPE, new Type[] {Type.getType(String.class)}),
                            null, null);
                    mv.visitCode();
                    mv.visitVarInsn(ALOAD, 0);
                    mv.visitVarInsn(ALOAD, 1);
                    mv.visitFieldInsn(PUTFIELD, this.getInternalName(), field.getName(), Type.getDescriptor(field.getType()));
                    mv.visitInsn(RETURN);
                    mv.visitMaxs(2, 2);
                    mv.visitEnd();
                }
                // getter
                MethodVisitor mv = cw.visitMethod(ACC_PUBLIC, field.getGetterName(), Type.getMethodDescriptor(Type.getType(field.getType()), new Type[]{}),
                        null, null);
                mv.visitCode();
                mv.visitVarInsn(ALOAD, 0);
                mv.visitFieldInsn(GETFIELD, this.getInternalName(), field.getName(), Type.getDescriptor(field.getType()));
                mv.visitInsn(ARETURN);
                mv.visitMaxs(1,1);
                mv.visitEnd();
            }
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
    
    public static int accessToOpcode(Access access){
        if (access == Access.PUBLIC) return ACC_PUBLIC;
        else return ACC_PROTECTED;
    }
    
    private void handleAnnotationArray(AnnotationVisitor arrayVisitor, Object[] a){
        for (Object value : a){
            handleAnnotationValue(arrayVisitor, null, value);
        }
    }
    
    private void handleAnnotationValue(AnnotationVisitor visitor, String name, Object value){
        if (value.getClass().isArray()){
            handleAnnotationArray(visitor.visitArray(name), (Object[]) value);
        }
        else if (value instanceof Enum){
            visitor.visitEnum(name, Type.getDescriptor(value.getClass()), value.toString());
        }
        //XXX: not handling the annotation case right now; no use case at the moment
        else{
            visitor.visit(name, value);
        }
    }
}
