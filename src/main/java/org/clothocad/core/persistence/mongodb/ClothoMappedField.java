/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.annotations.AlsoLoad;
import com.github.jmkgreen.morphia.annotations.ConstructorArgs;
import com.github.jmkgreen.morphia.annotations.Reference;
import com.github.jmkgreen.morphia.mapping.MappedField;
import com.github.jmkgreen.morphia.mapping.MappingException;
import com.github.jmkgreen.morphia.utils.ReflectionUtils;
import java.lang.annotation.Annotation;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.GenericArrayType;
import java.lang.reflect.GenericDeclaration;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.persistence.Add;
import org.clothocad.core.persistence.Replace;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author spaige
 */
public class ClothoMappedField extends MappedField{
    Logger log = LoggerFactory.getLogger(ClothoMappedField.class);
    
    protected Object value = null;
    protected Method provider = null;
    protected String fieldName;

    public static final String VIRTUAL_PREFIX = "__";
    
     /**
     * The constructor.
     */
    ClothoMappedField(Field f, Class<?> clazz) {
        f.setAccessible(true);
        field = f;
        persistedClass = clazz;
        discover();
    }
    
    ClothoMappedField(Add annotation, Class<?> clazz){
        this(annotation.provider(), annotation.value(), annotation.name(), annotation.concreteClass(), clazz);
        if (annotation.isReference()) this.foundAnnotations.put(Reference.class, new FakeReference());
    }
    
    ClothoMappedField(Replace annotation, MappedField f, Class<?> clazz){
        this(annotation.encoder(), annotation.value(), f.getJavaFieldName(), annotation.concreteClass(), clazz);
        this.foundAnnotations = f.getAnnotations();
    }
    
    
    ClothoMappedField(String provider, String value, String fieldName, Class concreteClass, Class<?> clazz) {
        this.fieldName = fieldName;
        persistedClass = clazz;
        realType = discoverVirtualFieldType(provider, value);
        ctor = discoverCTor(concreteClass);
        discoverMultivalued();
        
                // check the main type
        isMongoType = ReflectionUtils.isPropertyType(realType);

        // if the main type isn't supported by the Mongo, see if the subtype is.
        // works for T[], List<T>, Map<?, T>, where T is Long/String/etc.
        if (!isMongoType && subType != null)
            isMongoType = ReflectionUtils.isPropertyType(subType);

        if (!isMongoType && !isSingleValue && (subType == null || subType.equals(Object.class))) {
            log.warn("The multi-valued field '{}' is a possible heterogenous collection. It cannot be verified. Please declare a valid type to get rid of this warning. {}", getFullName(), subType);
            isMongoType = true;
        }
    }

    @Override
    public void addAnnotation(Class<? extends Annotation> clazz) {
        if (field == null){
            return;
        }
        super.addAnnotation(clazz); 
    }
    
    
    
    /*public static final List<Class<? extends Annotation>> interestingAnnotations = new ArrayList<>();
    static {
        interestingAnnotations.addAll(MappedField.interestingAnnotations);
        interestingAnnotations.add(Replace.class);
        interestingAnnotations.add(DBOnly.class);
    }*/

    @Override
    public void setFieldValue(Object classInst, Object value) throws IllegalArgumentException {
        if (field == null) return;
        super.setFieldValue(classInst, value); 
    }
    
   
    


    private Class discoverVirtualFieldType(String provider, String value) {
        Class type;
        Type gType;

        
        if (!"".equals(provider)){
            try {
                Method methodProvider = persistedClass.getMethod(provider, new Class[]{});
                type = methodProvider.getReturnType();
                gType = methodProvider.getGenericReturnType();
                this.provider = methodProvider;
            }
            catch (NoSuchMethodException me){
                    log.error("Could not find public method or field named {} on class {}", provider, persistedClass.getCanonicalName());
                    throw new IllegalArgumentException("Could not find named method.");
        
            }                
        } else {
            throw new UnsupportedOperationException("Add by value not implemented yet.");
        }
        
        TypeVariable<GenericDeclaration> tv = null;
        ParameterizedType pt = null;
        if (gType instanceof TypeVariable)
            tv = (TypeVariable<GenericDeclaration>) gType;
        else if (gType instanceof ParameterizedType)
            pt = (ParameterizedType) gType;

        if (tv != null) { 
//			type = ReflectionUtils.getTypeArgument(persistedClass, tv);
            Class typeArgument = ReflectionUtils.getTypeArgument(persistedClass, tv);
            if (typeArgument != null)
                type = typeArgument;
        } else if (pt != null) {
            if (log.isDebugEnabled())
                log.debug("found instance of ParameterizedType : {}", pt);
        }

        if (Object.class.equals(realType) && (tv != null || pt != null))
            log.warn("Parameterized types are treated as untyped Objects. See field '{}' on {}" , field.getName(), field.getDeclaringClass());

        if (type == null)
            throw new MappingException("A type could not be found for " + this.field);

        return type;
    }

    private Constructor discoverCTor(Class ctorType) {
        Constructor returnCtor = null;

        if (!ctorType.equals(Object.class))
            try {
                returnCtor = ctorType.getDeclaredConstructor();
                returnCtor.setAccessible(true);
            } catch (NoSuchMethodException e) {
                if (!hasAnnotation(ConstructorArgs.class))
                    log.warn("No usable constructor for {}" + ctorType.getName());
            }
        else {
            // see if we can create instances of the type used for declaration
            ctorType = getType();
            if (ctorType != null)
                try {
                    returnCtor = ctorType.getDeclaredConstructor();
                    returnCtor.setAccessible(true);
                } catch (NoSuchMethodException e) {
                    // never mind.
                } catch (SecurityException e) {
                    // never mind.
                }
        }
        return returnCtor;
    }

    private void discoverMultivalued() {
        //copied straight from MappedField
        if (realType.isArray() ||
                Collection.class.isAssignableFrom(realType) ||
                Map.class.isAssignableFrom(realType)) {

            isSingleValue = false;

            isList = List.class.isAssignableFrom(realType);
            isMap = Map.class.isAssignableFrom(realType);
            isSet = Set.class.isAssignableFrom(realType);
            //for debugging
            isCollection = Collection.class.isAssignableFrom(realType);
            isArray = realType.isArray();

            //for debugging with issue
            if (!isMap && !isSet && !isCollection && !isArray)
                throw new MappingException("type is not a map/set/collection/array : " + realType);

            // get the subtype T, T[]/List<T>/Map<?,T>; subtype of Long[], List<Long> is Long
            subType = (realType.isArray()) ? realType.getComponentType() : getParameterizedType((isMap) ? 1 : 0);

            if (isMap)
                mapKeyType = getParameterizedType(0);
        }
    }

    protected Type getParameterizedType(int i){
        if (field != null) return ReflectionUtils.getParameterizedType(field, i);
        else if (provider != null) return getParameterizedType(provider.getGenericReturnType(), i);
        //else if (value != null) return getParameterizedType(value.getClass().)
        else return null;
    }
    
    protected static Type getParameterizedType(final Type type, final int index) {
        if (type instanceof ParameterizedType) {
            ParameterizedType ptype = (ParameterizedType) type;
            if ((ptype.getActualTypeArguments() != null) && (ptype.getActualTypeArguments().length <= index)) {
                return null;
            }
            Type paramType = ptype.getActualTypeArguments()[index];
            if (paramType instanceof GenericArrayType) {
                return ((GenericArrayType) paramType).getGenericComponentType();
            } else {
                if (paramType instanceof ParameterizedType) {
                    return paramType;
                } else {
                    if (paramType instanceof TypeVariable) {
                        // TODO: Figure out what to do... Walk back up the to
                        // the parent class and try to get the variable type
                        // from the T/V/X
//						throw new MappingException("Generic Typed Class not supported:  <" + ((TypeVariable) paramType).getName() + "> = " + ((TypeVariable) paramType).getBounds()[0]);
                        return paramType;
                    } else if (paramType instanceof Class) {
                        return (Class) paramType;
                    } else {
                        throw new MappingException("Unknown type... pretty bad... call for help, wave your hands... yeah!");
                    }
                }
            }
        }
        return null;
    }
    
    @Override
    public Object getFieldValue(Object classInst) throws IllegalArgumentException {
        if (this.value != null) return value;
        if (this.provider != null) try {
            return provider.invoke(classInst);
        } catch (IllegalAccessException e) {
            log.error("Error accessing provider method", e);
            throw new RuntimeException(e);
        } catch (InvocationTargetException ex) {
            log.error("Exception in provider method", ex);
            return null;
        }
        return super.getFieldValue(classInst); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String getJavaFieldName() {
        if (fieldName != null) return fieldName;
        return super.getJavaFieldName(); //To change body of generated methods, choose Tools | Templates.
    }
   
    @Override
    public List<String> getLoadNames() {
        ArrayList<String> names = new ArrayList<String>();
        names.add(getNameToStore());

        AlsoLoad al = (AlsoLoad) this.foundAnnotations.get(AlsoLoad.class);
        if (al != null && al.value() != null && al.value().length > 0)
            names.addAll(Arrays.asList(al.value()));

        return names;
    }

    @Override
    public String getNameToStore() {
        String mappedName;
        try{
            mappedName = super.getNameToStore();
        } catch (NullPointerException e){
            mappedName = this.getJavaFieldName();
        }
        
        if (this.hasAnnotation(Add.class)) return VIRTUAL_PREFIX + mappedName;
        return mappedName;
    }

    @Override
    public String getFullName() {
        return getDeclaringClass().getName() +"."+ getJavaFieldName();
    }

    @Override
    public Class getDeclaringClass() {
        if (field == null) return persistedClass;
        return field.getDeclaringClass();
    }
    
    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getNameToStore()).append(" (");
        sb.append(" type:").append(realType.getSimpleName()).append(",");

        if (isSingleValue())
            sb.append(" single:true,");
        else {
            sb.append(" multiple:true,");
            sb.append(" subtype:").append(getSubClass()).append(",");
        }
        if (isMap()) {
            sb.append(" map:true,");
            if (getMapKeyClass() != null)
                sb.append(" map-key:").append(getMapKeyClass().getSimpleName());
            else
                sb.append(" map-key: class unknown! ");
        }

        if (isSet())
            sb.append(" set:true,");
        if (isCollection)
            sb.append(" collection:true,");
        if (isArray)
            sb.append(" array:true,");

        //remove last comma
        if (sb.charAt(sb.length() - 1) == ',')
            sb.setLength(sb.length() - 1);

        sb.append("); ").append(this.foundAnnotations.toString());
        return sb.toString();
    }
    
    
    
}
