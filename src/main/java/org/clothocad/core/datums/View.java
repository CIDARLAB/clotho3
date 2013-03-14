package org.clothocad.core.datums;

import java.util.List;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Permissions;
import org.clothocad.core.datums.util.ServerScript;

import flexjson.JSONDeserializer;


/**
 * A View is the factory for HTML widgets
 * 
 * It is one of the 4 primary Sharables
 * 
 * @author John Christopher Anderson
 */
public class View 
		extends Sharable {

    private View(Person author, 
            String name, 
            String description,
            List<ClothoField> inputArguments,
            ServerScript canUpdate, 
            String graphicsScript,
            String onShowScript,
            String onUpdateScript) {
    	
    	super(author, SharableType.VIEW);
        
        this.inputArguments = inputArguments;
        this.canUpdate = canUpdate;
        this.graphicsScript = graphicsScript;
        this.onShowScript = onShowScript;
        this.onUpdateScript = onUpdateScript;
        this.name = name;
        this.description = description;
        
    }

    public static View create(
            Person author, 
            String name, 
            String description,
            List<ClothoField> inputArgs,
            ServerScript canUpdate, 
            String graphicsScript,
            String onShowScript,
            String onUpdateScript) {
        
    	View view = new View(author, name, description, 
    			inputArgs, canUpdate, graphicsScript, 
    			onShowScript, onUpdateScript);
    	
        Collector.get().add(view);
        return view;
    }
    
    /***
    @Override
    public JSONObject toJSON() {
        try {
            JSONSerializer serializer = new JSONSerializer().exclude("*.class");
            serializer.prettyPrint(true);
            String serial = serializer.deepSerialize( this );
            return new JSONObject(serial);
        } catch (JSONException ex) {
            return null;
        }
    }
	***/
    
    public static View deserialize(String json) {
        View out = new JSONDeserializer<View>().deserialize(json, View.class);
        return out;
    }


    public String getDescription() {
        return description;
    }

    public int getInstanceCount() {
        return instanceCount;
    }

    public String getLargeIconURL() {
        return largeIconURL;
    }

    public String getName() {
        return name;
    }

    public Permissions getPermissions() {
        return permissions;
    }

    public String getSmallIconURL() {
        return smallIconURL;
    }

    public List<ClothoField> getInputArguments() {
        return inputArguments;
    }

    public ServerScript getCanUpdate() {
        return canUpdate;
    }

    public String getOnUpdateScript() {
        return onUpdateScript;
    }
    
    public String getGraphicsScript() {
        return this.graphicsScript;
    }
    
    public String getOnShowScript() {
        return this.onShowScript;
    }

    //Still need to implement this:
    private List<ClothoField> inputArguments;
    private ServerScript canUpdate;
    private String onUpdateScript;
    private String graphicsScript;
    private String onShowScript;
    
    //Permissions
    private Permissions permissions = new Permissions();

    //Metadata
    private String id;
    private String name;
    private String description;
    
    private String authorId;
    private String smallIconURL;
    private String largeIconURL;
    
    private int instanceCount = 0;
}
