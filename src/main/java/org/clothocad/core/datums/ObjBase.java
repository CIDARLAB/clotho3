package org.clothocad.core.datums;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import lombok.AccessLevel;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.bson.types.ObjectId;

import java.util.Arrays;
import java.util.Date;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.core.persistence.Reference;
import org.clothocad.core.security.Visibility;

/**
 *
 * @author spaige
 */
@EqualsAndHashCode(exclude = {"dateCreated", "lastModified", "lastAccessed", "isDeleted"})
@Data()
@NoArgsConstructor
@Slf4j
public abstract class ObjBase {

    //add schema
    //remove schema
    //can only manipulate schema set if you have write privs
    public ObjBase(String name) {
        this.name = name;
    }
    private ObjectId id;
    private String name;
    @DBOnly
    private boolean isDeleted;
    @Setter(AccessLevel.NONE)
    private Date dateCreated;
    @DBOnly
    private Date lastModified, lastAccessed;
    
    @Getter
    @Setter
    private Visibility visibility;

    public void onUpdate() {

        System.out.println("[ObjBase.onUpdate] ERNST's Task!! -> Push object via pubsub.");


        // here we need to call the client-side API
        // which forwards the update message 
        // to ``subscribed'' clients

        //so, do all setters need to check to see if the value changed and then call onUpdate?

    }

    public List<ObjBase> getChildren() {
        ArrayList<ObjBase> children = new ArrayList<>();

        for (Field f : getAllReferences(this.getClass())) {
            boolean accessible = f.isAccessible();
            try {
                f.setAccessible(true);
                Object value = f.get(this);
                if (value == null) continue;
                //reference might be a collection of references
                if (java.util.Collection.class.isInstance(value)) {
                    //TODO: not typesafe
                    children.addAll((java.util.Collection) value);

                } 
                else if (value.getClass().isArray()){
                    List list = Arrays.asList((Object[])value);
                    children.addAll(list);
                }
                else {
                    children.add((ObjBase) value);
                }
            } catch (IllegalArgumentException | IllegalAccessException ex) {
                log.error("getChildren: ", ex);
            } finally {
                f.setAccessible(accessible);
            }
        }

        return children;
    }

    private static List<Field> getAllReferences(Class c) {
        ArrayList<Field> output = new ArrayList<>();
        while (c != null && c != Object.class) {
            for (Field f : c.getDeclaredFields()) {
                if (f.getAnnotation(Reference.class) != null) {
                    output.add(f);
                }
            }
            c = c.getSuperclass();
        }
        return output;
    }
}
