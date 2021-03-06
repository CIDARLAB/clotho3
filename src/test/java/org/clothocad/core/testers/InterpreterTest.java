package org.clothocad.core.testers;

import com.fasterxml.jackson.core.JsonParseException;
import java.util.Set;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.aspects.Interpreter.Trainer;
import org.clothocad.core.aspects.Interpreter.StdIn;

class InterpreterTest {
    /**
     * Simulating the client input with a StdIn.
     * Set a threshold of when to execute directly and when to ask.
     * @param args
     */
    public static void main(String[] args) throws JsonParseException{
        StdIn inputReader = new StdIn();
        while (true) {
            System.out.println("Type your command: ");
            String cmd = inputReader.readString();
            if (cmd.equalsIgnoreCase("Train")) {
                Trainer.inputTrainingData();
            } else if (cmd.equalsIgnoreCase("Stop")) {
                break;   
            } else {
                Set<String> results = Interpreter.get().receiveNative(cmd);
                for (String str : results) {
                    System.out.println(str);
                }
            }
        }
    }
}
