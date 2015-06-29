package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public abstract class Module extends SharableObjBase {

    @NotNull
    @Getter
    @Setter
    protected ModuleRole role;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Influence> influences;

    @Getter
    @Setter
    @Reference
    protected Module parentModule;

    public Module(String name, ModuleRole role, Person author) {
        super(name, author);
        this.role = role;
    }

    public Module(String name, String description, ModuleRole role, Person author) {
        super(name, author, description);
        this.role = role;
    }

    // Feel free to add more of these
    public static enum ModuleRole {
        TRANSCRIPTION, TRANSLATION, EXPRESSION, COMPARTMENTALIZATION, LOCALIZATION, SENSOR, REPORTER, ACTIVATION, REPRESSION;
    }
    
    public void addInfluence(Influence influence) {
    	if (influences == null) {
    		influences = new HashSet<Influence>();
    	}
    	influences.add(influence);
    }

}
