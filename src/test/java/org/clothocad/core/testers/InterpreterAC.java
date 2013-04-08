package org.clothocad.core.testers;

import java.util.ArrayList;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.util.Logger;

class InterpreterAC {
    public static AutoComplete completer1() {
        String[] word_bank1 = {
            "walk my cat",
            "sail on my boat",
            "samba",
            "walk my dog",
            "walk my friend",
            "walk too much",
            "walking",
            "walker",
            "Walter",
        };
//        return new AutoComplete(word_bank1);
        return new AutoComplete();
    }

    public static AutoComplete completer2() {
        String[] word_bank2 = {
            "complete this sentence",
            "Complete this sentence",
            "",
            "complete thissentence",
            "complete this sentences",
            "complete this sentence",
            "COMPLETE THIS SENTENCE"
        };
//        return new AutoComplete(word_bank2);
        return new AutoComplete();
    }

    public static void test1() {
        AutoComplete completer1 = completer1();
        ArrayList<String> results = completer1.getCompletions("wa");
        check(results.get(0), "walk my cat");
        check(results.get(1), "walk my dog");
        check(results.get(2), "walk my friend");
        check(results.size(), 6);
    }

    public static void test2() {
        AutoComplete completer1 = completer1();
        ArrayList<String> results = completer1.getCompletions("sa");
        check(results.get(0), "sail on my boat");
        check(results.get(1), "samba");
        check(results.size(), 2);
    }

    public static void test3() {
        AutoComplete completer2 = completer2();
        ArrayList<String> results = completer2.getCompletions("com");
        check(results.get(0), "complete this sentence");
        check(results.get(2), "complete thissentence");
        check(results.size(), 3);
    }

    public static void check(String x, String y) {
        if (!x.equals(y)) {
            Logger.log(Logger.Level.WARN, x + " is not " + y);
        }
    }
     public static void check(int actual, int expected) {
         if (actual != expected) {
             Logger.log(Logger.Level.WARN, "Expected: " + expected +
                                            "\nGot: " + actual);
         }
     }

    public static void main(String[] args) {
        test1();
        test2();
        test3();
    }
}