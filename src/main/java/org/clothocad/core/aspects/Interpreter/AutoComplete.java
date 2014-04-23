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

public class AutoComplete {
    Persistor persistor;
    
    /* AutoComplete Contructor */
    public AutoComplete () {
        
        trie = new PatriciaTrie<String, Object> (StringKeyAnalyzer.INSTANCE);
        //Load up the words from the word bank into the Trie
        for(String word : getWordBank()) {
            trie.put(word, map.get(word));
        }

    }
    /*
     * New constructor to accept the persistor from the ServerSideAPI
     */
    public AutoComplete(Persistor persistorMongo) {
        persistor = persistorMongo;
        trie = new PatriciaTrie<String, Object> (StringKeyAnalyzer.INSTANCE);
        map = new HashMap();
        for(String word: getWordBank()){
            HashMap temp = new HashMap();
            temp.put("name",word);
            temp.put("uuid", map.get(word));
            temp.put("text", word);
            temp.put("type", "phrase");
            trie.put(word.toLowerCase(),temp);
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
    public void put(String word) {
        try {
            System.out.println("Trie is going to store: " + word);
            trie.put(word, word);
            wordBank.add(word);
            //FIXME: persistor.save(wordBank); (make a wordbank class)
            System.out.println("Stephanie, change persistance to db instead of flatfile");
            FileUtils.writeFile(wordBank.toString(), "wordbank.txt");
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    /*
     * Creates a wordbank using the list of tuples from the persistor
     */
    private Set<String> getWordBank() {
        try {
            if(wordBank==null) {
                Tuple[] temp = persistor.getTuplesMongo();
                //System.out.println("Temp size: " + temp.length);
                //String sfile = FileUtils.readFile("wordbank.txt");
                //List listy = JSON.deserializeList(sfile);
                //if (listy == null) return new HashSet<>(); //XXX
                wordBank = new HashSet<String>();
                //for(int i=0; i<listy.size(); i++) {
                //    String str = listy.get(i).toString();
                //    wordBank.add(str);
                //}
                for(int i = 0; i< temp.length;i++){
                    if (temp[i].get(0) != null){
                        String str = temp[i].get(0).toString();
                        Object uuid = temp[i].get(1);
                        wordBank.add(str);
                        map.put(str,uuid);
                    }
                }
                return wordBank;
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        return wordBank;
    }

    private Trie<String, Object> trie;
    transient private Set<String> wordBank;
    private HashMap map;
    
    public Object getUUID(String key){

        HashMap output = (HashMap) trie.selectValue(key);
        return output.get("uuid");
    }
    
}
