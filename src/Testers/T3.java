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

package Testers;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.aspects.Logger;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.json.JSONException;
import org.json.JSONObject;


/**
 * This test creates the makeabuddy assisant, which I've already renamed accordingly the uuid
 * @author John Christopher Anderson
 */

public class T3 {
    public static void main(String[] args) {
        createMakeABuddyFunction();
        createCindyPersonInstance();
        System.out.println("T3:  All is good!");
    }
    
    
    public static Function createMakeABuddyFunction() {
        try {
            //Start constructing the bits to describe an Function
            ServerScript dooit = new ServerScript(
                      "clotho.log(\"INFO\", '+++running Dooit in js from makeabuddy');\n"
                    + "clotho.log(\"INFO\", '+++cindy of makeabuddy: ' + cindy);\n"
                    + "println('cindy has middlename ' + cindy.middlename);\n"
                    + "return cindy;",
                       Language.JavaScript);
            
            ServerScript candooit = new ServerScript("return true;", Language.JavaScript);
            
            Schema personSchema = Person.getSchema();

            
            //Create input and output arguments
            List<ClothoField> inputArgs = new ArrayList<ClothoField>();
            ClothoField cindy = new ClothoField("cindy", FieldType.SCHEMA, personSchema.getId(), 1);
            inputArgs.add(cindy);
                
            List<ClothoField> outputArgs = new ArrayList<ClothoField>();
            ClothoField herbuddy = new ClothoField("herbuddy", FieldType.SCHEMA, personSchema.getId(), 1);
            outputArgs.add(herbuddy);

            //Create an Function MakeABuddy
            Function asst = Function.create(
                    Person.getAdmin(),
                    "MakeABuddy",
                    "Inputs a personSchema and returns a buddy with her name included in his",
                    dooit, 
                    candooit, 
                    inputArgs, 
                    outputArgs);
            
            
            //Change the Id
            JSONObject obj = asst.toJSON();
            obj.put("id", "specific-makeabuddy-is-uuid");
            asst = Function.deserialize(obj.toString());
            Persistor.get().persistDatum(asst);
            Logger.log(Logger.Level.INFO, asst.getId() + "   " + asst.getName() + "\n" + asst.getDescription() + "\n...was created successfully, all good!");
            return asst;
        } catch (JSONException ex) {
            ex.printStackTrace();
            System.exit(0);
            return null;
        }
    }
    
    public static Person createCindyPersonInstance() {
        try {
            JSONObject fields = new JSONObject();
            
            fields.put("email",  "cindy@cindysu.com" );
            fields.put("displayname",  "cindysu" );
            fields.put("givenname",  "Cindy" );
            fields.put("middlename",  "Meredith" );
            fields.put("surname",  "Su" );
            fields.put("nickname",  "Cindy" );

            String str = fields.toString();

            Instance obj = Instance.create(Person.getAdmin(), Person.getSchema(), str);
            Person cindy = Person.create(obj);
            
            //Change the Id
            JSONObject json = cindy.toJSON();
            json.put("id", "specific-cindysu-is-uuid");
            obj = Instance.deserialize(json.toString());
            cindy = Person.create(obj);
            Persistor.get().persistDatum(cindy);
            Logger.log(Logger.Level.INFO, cindy.getId() + "   " + cindy.getName() + "\n...was created successfully, all good!");
            return cindy;
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(0);
            return null;
        }
}
}
