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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.aspects.Interpreter.RadixTrie.PatriciaTrie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.Trie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.StringKeyAnalyzer;
import org.clothocad.core.util.FileUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

public class AutoComplete {
    Persistor persistor;
    
    /* AutoComplete Contructor */
    public AutoComplete () {
        trie = new PatriciaTrie<String, String> (StringKeyAnalyzer.CHAR);
        
        //Load up the words from the word bank into the Trie
        for(String word : getWordBank()) {
            trie.put(word, word);
        }
    }

    /**
     * JCA:  this is the important method -- when constructing the message
     * to send autocompletes to the Client, this is where those are extracted.
     * 
     * Extracting options from the Trie
     */
    public ArrayList<String> getCompletions(String subString) {
        SortedMap<String, String> subTrie = trie.prefixMap(subString);
        ArrayList<String> options = new ArrayList<String>();
        for (Map.Entry<String, String> entry : subTrie.entrySet()) {
            options.add(entry.getValue());
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
    
    private Set<String> getWordBank() {
        try {
            if(wordBank==null) {
                String sfile = FileUtils.readFile("wordbank.txt");
                JSONArray listy = new JSONArray(sfile);
                wordBank = new HashSet<String>();
                for(int i=0; i<listy.length(); i++) {
                    String str = listy.getString(i);
                    wordBank.add(str);
                }
                return wordBank;
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
        return wordBank;
    }

    private Trie<String, String> trie;
    transient private Set<String> wordBank;
}
