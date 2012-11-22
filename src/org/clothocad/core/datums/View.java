package org.clothocad.core.datums;

import flexjson.JSONDeserializer;
import flexjson.JSONSerializer;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.util.ServerScript;
import org.json.JSONObject;
import org.json.JSONException;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Permissions;


/**
 * A View is the factory for HTML widgets
 * 
 * It is one of the 4 primary Sharables
 * 
 * @author John Christopher Anderson
 */
public class View implements Sharable {

    private View() {
        
    }

    public static View create(
            Person author, 
            String name, 
            String description,
            List<ClothoField> inputArguments,
            ServerScript canUpdate, 
            String graphicsScript,
            String onShowScript,
            String onUpdateScript) {
        
        View out = new View();
        if(author!=null) {
            out.authorId = author.getId();
        }
        out.inputArguments = inputArguments;
        out.canUpdate = canUpdate;
        out.graphicsScript = graphicsScript;
        out.onShowScript = onShowScript;
        out.onUpdateScript = onUpdateScript;
        out.name = name;
        out.description = description;
        
        Collector.get().add(out);
        return out;
    }
    
    @Override
    public Person extractAuthor() {
        Person out = (Person) Collector.get().getDatum(authorId);
        return out;
    }

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

    public static View deserialize(String json) {
        View out = new JSONDeserializer<View>().deserialize(json, View.class);
        return out;
    }


    @Override
    public String getId() {
        return id;
    }

    public String getAuthorId() {
        return authorId;
    }

    public ClothoDate getDateCreated() {
        return dateCreated;
    }

    public ClothoDate getDateLastAccessed() {
        return dateLastAccessed;
    }

    public ClothoDate getDateLastModified() {
        return dateLastModified;
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

    @Override
    public SharableType type() {
        return SharableType.VIEW;
    }

    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        //Check that the requestor has set permissions on this object, if not return false
        
        //Run validation on newvalue, if it fails, return false
        
        
        //set the new data
        View newView = deserialize(newvalue.toString());
        
        //Save things
        Persistor.get().persistDatum(newView);
        return true;
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
    
    private ClothoDate dateCreated = new ClothoDate();
    private ClothoDate dateLastModified = new ClothoDate();
    private ClothoDate dateLastAccessed = new ClothoDate();


}
