/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.InjectableValues;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class InstantiatedReferencesCache extends InjectableValues.Std{
    
    List<ObjBase> undone = new ArrayList<>();
    
    public void addObjBase(Object valueId, ObjBase value){
        addValue((String) valueId, value);
        undone.add(value);
    }    

    @Override
    public Object findInjectableValue(Object valueId, DeserializationContext ctxt, BeanProperty forProperty, Object beanInstance) {
        
        if (ObjBase.class.isAssignableFrom(forProperty.getType().getRawClass())){
            Class<ObjBase> c = (Class<ObjBase>) forProperty.getType().getRawClass();
            try {
                //already made an object for this id - return it
                return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
            } catch (IllegalArgumentException e){
                // else, create it, cache it, return it
                //XXX: just assuming there's a lombok no-arg constructor for now
                try{
                    ObjBase newValue = c.newInstance();
                    newValue.setId(new ObjectId(valueId));
                    this.addObjBase(valueId, newValue);
                    return newValue;
                } catch (InstantiationException|IllegalAccessException ie){
                    throw new RuntimeException(ie);
                }
            }
            
        } else {
            return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
        }
    }
    
    public Boolean done(){
        return undone.isEmpty();
    }
    
    public ObjBase getNextUndone(){
        return undone.remove(0);
    }

    boolean has(String id) {
        return _values.containsKey(id);
    }
    
    public Object get(String id){
        return _values.get(id);
    }
}
