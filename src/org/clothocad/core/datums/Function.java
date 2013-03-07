/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.datums;

import flexjson.JSONDeserializer;
import flexjson.JSONSerializer;
import java.util.List;
import org.clothocad.core.datums.util.ServerScript;
import java.util.Map;
import java.util.UUID;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Permissions;
import org.json.JSONObject;

/**
 * TODO: *ADD ADD AND REMOVE METHODS FOR ARGUMENTS
 * @author John Christopher Anderson
 */
public class Function 
		extends Sharable {

    private Function(Person author,
            String name, 
            String description,
            ServerScript dooIt,
            ServerScript canDooIt,
            List<ClothoField> inputArguments,
            List<ClothoField> outputArguments) {
    	
    	super(author, SharableType.FUNCTION);
    	
        this.name = name;
        this.description = description;
        this.id = UUID.randomUUID().toString();
        this.dooIt = dooIt;
        this.canDooIt = canDooIt;
        this.inputArguments = inputArguments;
        this.outputArguments = outputArguments;
    }

    /**
     * An Assistant is the functional 'bean' that inputs DataFields like ObjBases, alters
     * them and outputs new ones.
     * 
     * @param name is a name for the Assistant
     * @param description a natural language description of the function of the assistant
     * @param dooIt a Script that inputs the inputArguments, runs its thing, and returns outputArguments
     * @param canDooIt a Script that inputs the inputArguments and returns true or false
     * @param inputArguments a map of tokens, like "my_sequence" to a Schema UUID (it remains loosely-coupled)
     * @param outputArguments a map of tokens, like "output_seq" to a Schema UUID (it remains loosely-coupled)
     */
    public static Function create(
            Person author,
            String name, 
            String description,
            ServerScript dooIt,
            ServerScript canDooIt,
            List<ClothoField> inputArgs,
            List<ClothoField> outputArgs ) {
        
        Function out = new Function(
        		author, name, description, dooIt, canDooIt, inputArgs, outputArgs);
        
        Collector.get().add(out);
        return out;
    }


    public boolean validate(JSONObject jsondata) {
        try {
            String resultStr = this.canDooIt.run(jsondata);
            Boolean bool = Boolean.parseBoolean(resultStr);
            return bool;
        } catch (Exception ex) {
            System.out.println("Assistant: This needs to realy a developer message as the validate method failed to run properly");
            ex.printStackTrace();
            return false;
        }

    }


    /***
    @Override
    public JSONObject toJSON() {
        try {
            JSONSerializer serializer = new JSONSerializer().exclude("*.class");
            serializer.prettyPrint(true);
            String serial = serializer.deepSerialize( this );
            return new JSONObject(serial);
        } catch (Exception ex) {
            return null;
        }
    }
    **/
    
    public static Function deserialize(String json) {
        Function out = new JSONDeserializer<Function>().deserialize(json, Function.class);
        return out;
    }
    
    @Override
    public String getId() {
        return id;
    }

    public String getAuthorId() {
        return authorId;
    }

    public ServerScript getCanDooIt() {
        return canDooIt;
    }

    public String getDescription() {
        return description;
    }

    public ServerScript getDooIt() {
        return dooIt;
    }

    public List<ClothoField> getInputArguments() {
        return inputArguments;
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

    public List<ClothoField> getOutputArguments() {
        return outputArguments;
    }

    public Permissions getPermissions() {
        return permissions;
    }

    public String getSmallIconURL() {
        return smallIconURL;
    }

    /***
    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public SharableType type() {
        return SharableType.FUNCTION;
    }
    ***/
    
    private ServerScript dooIt;
    private ServerScript canDooIt;

    //The inputs and outputs are namespaced within the script by these String fieldtokens
    private List<ClothoField> inputArguments;
    private List<ClothoField> outputArguments;

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
