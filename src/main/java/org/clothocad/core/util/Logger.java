/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.util;

/**
 * @author Kelvin Li
 */
public class Logger {
    public static void log(Level level, Object mesg_obj) {
        StackTraceElement stackelement = getArtificialStackTrace();
        real_log(level, mesg_obj, null, getArtificialStackTrace());
    }

    public static void log(Level level, Object mesg_obj, Throwable ex) {
        StackTraceElement stackelement = getArtificialStackTrace();
        real_log(level, mesg_obj, ex, getArtificialStackTrace());
    }

    public static enum Level {
        INFO,
        WARN,
        FATAL,
    }

    private static StackTraceElement getArtificialStackTrace() {
        Throwable myex = new Exception();
        StackTraceElement stackelement = null;
        try {
            stackelement = myex.getStackTrace()[2];
        } catch (ArrayIndexOutOfBoundsException e) {
            log(Level.FATAL, "Could not obtain caller information", e);
        }
        return stackelement;
    }

    private static void real_log(Level level,
                                 Object mesg_obj,
                                 Throwable ex,
                                 StackTraceElement stackelement) {
        if (!doOrigin(formatMsgOrigin(stackelement))) {
            return;
        }
        String out_string = formatMsg(level,
                                      String.valueOf(mesg_obj),
                                      stackelement);
        logPrint(out_string);
        if (ex != null) {
            logPrint("-------Begin Stack Trace--------\n");
            ex.printStackTrace(System.out);
            logPrint("-------End of Stack Trace-------\n");
        }
        logPrint("\n");
    }

    private static String formatMsg(Logger.Level level,
                                    String message,
                                    StackTraceElement stackelement) {
        String out = "";

        switch (level) {
        case INFO:
            out = "INFO";
        break;
        case WARN:
            out = "WARN";
        break;
        case FATAL:
            out = "FATAL";
        break;
        default:
            return "ILLEGAL LOG LEVEL";
        }

        out += ": " + message + "\n";
        out += "    (" + formatMsgOrigin(stackelement) + ")\n";
        return out;
    }

    private static String formatMsgOrigin(StackTraceElement ste) {
        String origin = ste.getClassName() + "$" + ste.getMethodName() + "~";
        origin += ste.getFileName() + "@" + ste.getLineNumber();
        return origin;
    }

    private static boolean doOrigin(String origin) {
        return (!hasPrefix(ignored_origins, origin) && hasPrefix(allowed_origins, origin));
    }

    private static boolean hasPrefix(String[] prefixes, String msg) {
        /* TODO: need to check for null? */
        for (String p : prefixes) {
            if (msg.startsWith(p)) {
                return true;
            }
        }
        return false;
    }

    private static void logPrint(String s) {
        System.out.print(s);
    }

    /* Logger suppresses every message whose origin matches any prefix below. */
    private static String[] ignored_origins = {
        /* EXAMPLES
        "org.clothocad.core.aspects",
        "org.clothocad.core.aspects.Router.RouterServer",
        "org.clothocad.core.aspects.Router.RouterServer$run",
        "org.clothocad.core.aspects.Router.RouterServer$run~RouterServer.java@73",
         */

        /* Uncomment to shut-up Logger completely */
        /* "", */
    };

    /* Logger prints only the messages whose origins match a prefix below. */
    private static String[] allowed_origins = {
        "", /* allow all */
    };
}
