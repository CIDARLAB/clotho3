package org.clothocad.core.testers;

import java.util.Set;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


class IntprT1 {
    static final Logger logger = LoggerFactory.getLogger(IntprT1.class);
    public static void test1() {
        Interpreter.get().learnNative("walk my cat", "run(assistantWalker, Cat)");
        Interpreter.get().learnNative("walk my dog", "run(assistantWalker, Dog)");
        Set<String> results = Interpreter.get().receiveNative("walk my cat");
        String result0 = (String) (results.toArray())[0];
        String result1 = (String) (results.toArray())[1];
        check(result0, "run(assistantWalker, Cat)");
        check(result1, "run(assistantWalker, Dog)");
        Interpreter.get().forgetNative("walk my cat", "run(assistantWalker, Cat)");
        Interpreter.get().forgetNative("walk my dog", "run(assistantWalker, Dog)");
    }

    public static void test2() {
        multiTrain("walk my cat", "run(assistantWalker, Cat)", 40);
        multiTrain("walk my dog", "run(assistantWalker, Dog)", 70);
        Set<String> resultsA = Interpreter.get().receiveNative("walk my cat");
        Set<String> resultsB = Interpreter.get().receiveNative("walk my dog");
        Set<String> resultsC = Interpreter.get().receiveNative("walk dog");
        Set<String> resultsD = Interpreter.get().receiveNative("walk cat");
        Set<String> resultsE = Interpreter.get().receiveNative("my cat");
        Set<String> resultsF = Interpreter.get().receiveNative("my dog");

        String result0 = (String) (resultsA.toArray())[0];
        String result1 = (String) (resultsA.toArray())[1];
        String result2 = (String) (resultsB.toArray())[0];
        String result3 = (String) (resultsC.toArray())[0];
        String result4 = (String) (resultsD.toArray())[0];
        String result5 = (String) (resultsE.toArray())[0];
        String result6 = (String) (resultsF.toArray())[0];

        check(result0, "run(assistantWalker, Cat)");
        check(result1, "run(assistantWalker, Dog)");
        check(result2, "run(assistantWalker, Dog)");
        check(result3, "run(assistantWalker, Dog)");
        check(result4, "run(assistantWalker, Cat)");
        check(result5, "run(assistantWalker, Cat)");
        check(result6, "run(assistantWalker, Dog)");
        multiTrain("walk my cat", "run(assistantWalker, Cat)", -40);
        multiTrain("walk my dog", "run(assistantWalker, Dog)", -70);
    }

    public static void check(String act, String pred) {
        //TODO: convert to Assert.assertEqual
        if (!act.equals(pred)) {
            logger.warn( act + " is not " + pred);
        }
    }

    public static void multiTrain(String cmd, String action, int reps) {
        int temp = Math.abs(reps);
        if (reps > 0) {
            while (temp > 0) {
                Interpreter.get().learnNative(cmd, action);
                temp --;
            }
        } else {
            while (temp > 0) {
                Interpreter.get().forgetNative(cmd, action);
                temp --;
            }
        }
    }

    public static void main(String[] args) {
        test1();
        test2();
    }
}
