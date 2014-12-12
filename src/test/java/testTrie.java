/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import org.clothocad.core.aspects.Interpreter.RadixTrie.PatriciaTrie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.StringKeyAnalyzer;
/**
 *
 * @author jingeunlee
 */
public class testTrie {
    public static void main(String[] args){
        PatriciaTrie<String,String> p = new PatriciaTrie<String,String>(StringKeyAnalyzer.CHAR);
        p.put("ab", "Hello");
        p.put("ac", "Bye");
        String x = p.select("a").getValue();
        System.out.println(x);
    }
}
