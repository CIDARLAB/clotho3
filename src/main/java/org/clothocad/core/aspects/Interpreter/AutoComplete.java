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
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.aspects.Interpreter.RadixTrie.PatriciaTrie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.Trie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.StringKeyAnalyzer;

public class AutoComplete {
    /* AutoComplete Contructor */
    public AutoComplete () {
        trie = new PatriciaTrie<String, String> (StringKeyAnalyzer.CHAR);
        
        //Load up the words from the word bank into the Trie
        for (String word : getWordBank()) {
            //JCA:  the Trie holds mappings of phrases to displayed results, so the wordBank really should be a mapping of these things
            //ehhh, no, there are two separate things going on.  This is just autocomplete, and it's holding all phrases verbatim, so it's
            //not necessarily the same as the disambiguation bit.  These are just exact phrases people have put in before
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
        System.out.println("Trie is going to store: " + word);
        trie.put(word, word);
        if(!wordBank.contains(word)) {
            wordBank.add(word);
            //Persistor.get().persistWordBank(wordBank);
        }
    }
    
    private List<String> getWordBank() {
        if(wordBank==null) {
            //wordBank = Persistor.get().loadWordBank();
            if(wordBank==null) {
                wordBank = new ArrayList<String>();
            }
        }
        return wordBank;
    }

    private Trie<String, String> trie;
    transient private List<String> wordBank;
}
