/*
 * Copyright 2005-2012 Roger Kapsi, Sam Berlin
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

/**
 * @usedby Michael Lin
 */

package org.clothocad.core.aspects.Interpreter;

import groovy.lang.Tuple;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.aspects.Interpreter.RadixTrie.PatriciaTrie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.Trie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.StringKeyAnalyzer;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.JSON;

public class GlobalTrie {
    /* AutoComplete Contructor */
    /*This constructor is never used except in the InterpreterAC test*/
    public GlobalTrie () {
        
        trie = new PatriciaTrie<String, Object> (StringKeyAnalyzer.INSTANCE);
        //Load up the words from the word bank into the Trie


    }
    
    /*
     * New constructor to accept the persigstor from the ServerSideAPI
     */
    public GlobalTrie(List<Map> data) {
        trie = new PatriciaTrie<String, Object> (StringKeyAnalyzer.INSTANCE);
        List<Map> temp = data;
        for(Map map : temp){
            try {
                String name = (String) map.get("name");
                trie.put(name.toLowerCase(), map);
            } catch(Exception err) {
                System.out.println("Couldn't add to autocomplete Trie: " + map.toString());
            }
        }
    }

    /**
     * JCA:  this is the important method -- when constructing the message
     * to send autocompletes to the Client, this is where those are extracted.
     * 
     * Extracting options from the Trie
     */
    public List<Map> getCompletions(String subString) {
        SortedMap<String, Object> subTrie = trie.prefixMap(subString.toLowerCase());
        //System.out.println("Size of subtrie: " + subTrie.size());
        List<Map> options = new ArrayList<>();
        for (Map.Entry<String, Object> entry : subTrie.entrySet()) {
            HashMap tempMap = (HashMap) entry.getValue();
            options.add(tempMap);
        }
        return options;
    }

    /**
     * Adds new options into the trie
     */
    public void put(Map data) {
        try {
            Map map = new HashMap();
            map.put("name", data.get("name"));
            map.put("id", data.get("id").toString());
            map.put("schema", data.get("schema"));
            if(data.containsKey("description")) {
                map.put("description", data.get("description"));
            } else if(data.containsKey("shortDescription")) {
                map.put("description", data.get("shortDescription"));
            }
            
            String name = (String) data.get("name");
            trie.put(name.toLowerCase(), map);
        } catch(Exception err) {
            System.out.println("Couldn't add new object to autocomplete Trie: " + data.toString());
        }
    }

    private Trie<String, Object> trie;
}
