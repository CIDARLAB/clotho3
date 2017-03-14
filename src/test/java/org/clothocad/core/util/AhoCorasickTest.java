/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;


import java.util.Collection;
import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.junit.Test;

/**
 *
 * @author David
 */
public class AhoCorasickTest {
    
    
    @Test
    public void howToTest()
    {
        String speech = "I have a lot of words and words are really nice sometimes but you know at the end of the day fuck words";
        
        Trie trie = Trie.builder().removeOverlaps().caseInsensitive()
                .addKeyword("word")
                .addKeyword("time")
                .addKeyword("end")
                .build();
        
        Collection<Emit> parseText = trie.parseText(speech);
        
        for(Emit e : parseText)
        {
            System.out.println(e.getKeyword());
            System.out.println("  Start : " + e.getStart());
            System.out.println("  End : " + e.getEnd());
        }
        
    }
    
}
