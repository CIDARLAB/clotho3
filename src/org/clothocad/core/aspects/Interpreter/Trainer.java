package org.clothocad.core.aspects.Interpreter;

import org.clothocad.core.aspects.Logger;

/**
 * Singleton trainer used as developer interface to store all command-to-action 
 * arguments in .txt file.
 */
/* TODO: make not public. Finally intgrate some parts into Communicator which is 
 * the only entry point from the outside */
public class Trainer {

    /* Instantiating Trainers own scanner. This is temporary */
    private static StdIn inputReader = new StdIn();

    /**
     * Allows user to input and save all the command/action relationships.
     */
    public static void inputTrainingData(){
        while (true) {
            System.out.println("Type in relationship 'command/action': ");
            String input = inputReader.readString();
            switch (checkInput(input)) {
            case INVALID_COMMAND:
                System.out.println("Invalid Command");
                break;
            case ARGS_NOT_EXACTLY_TWO:
                System.out.println("Train with exactly two arguments");
                break;
            case BADLY_FORMED_MULTI_TRAIN:
                Logger.log(Logger.Level.WARN,
                       "Badly formed multi train, use 'Multi train:(num_rep)");
                break;
            case OK:
                break;
            default:
                Logger.log(Logger.Level.FATAL, "Badly formed enums");
                break;
            }

            if (input.equalsIgnoreCase("Stop Train")) {
                break;
            } else if(input.contains("Multi train")) {
                String[] inArr = input.split(":");
                multiTrain(Integer.parseInt(inArr[1].trim()));
            } else if (input.contains("/")) {
                String[] cmdAct = input.split("/");
                Interpreter.get().learnNative(cmdAct[0], cmdAct[1]);
            } 
        }
    }

    /**
     * Temporary hacky way to test the NGrams algorithm
     */
    private static void multiTrain (int rep) {
        int reps = rep;
        System.out.println("MultiTrain~ Type in relationship 'command/action': ");
        String input = inputReader.readString();
        if (input.contains("/")) {
            while (reps > 0) {
                String[] cmdAct = input.split("/");
                Interpreter.get().learnNative(cmdAct[0], cmdAct[1]);
                reps -= 1;
            }
        } else {
            System.out.println("Train with exactly two arguments");
        }
    }

    private static ErrorBag checkInput(String in) {
        if (in.contains("/") && (in.split("/").length != 2)) {
            return ErrorBag.ARGS_NOT_EXACTLY_TWO;

        } else if (in.startsWith("Multi train") &&
                   !(in.contains(":") && in.split(":").length == 2)) {
            return ErrorBag.BADLY_FORMED_MULTI_TRAIN;

        } else if (!in.contains("/") &&
                   !in.equalsIgnoreCase("Stop Train") &&
                   !in.equalsIgnoreCase("View Train Set") &&
                   !in.startsWith("Multi train")) {
            return ErrorBag.INVALID_COMMAND;

        } else {
            return ErrorBag.OK;
        }
    }

    public static enum ErrorBag {
        INVALID_COMMAND,
        BADLY_FORMED_MULTI_TRAIN,
        ARGS_NOT_EXACTLY_TWO,
        OK,
    }
}
