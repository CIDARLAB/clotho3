package org.clothocad.core.execution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.execution.subprocess.SubprocessExec;
import org.clothocad.core.execution.Mind;

public class ExecTimingTest {
    private static final SubprocessExec.EventHandler noopHandler =
    new SubprocessExec.EventHandler () {
        @Override public void onFail(final String x) {}
        @Override public void onSuccess(final String x) {}
    };

    public static void
    main(String[] unused) throws Exception
    {
        testJS();
        testPython();
    }

    private static void
    testJS() throws Exception
    {
        final Mind mind = new Mind();
        final long startTime = System.currentTimeMillis();
        for (int i = 0; i < 1000; i++) {
            final Object out = mind.evalFunction(
                "function tester () {return 2;}",
                "tester",
                new ArrayList<Object>(),
                null
            );
            if (!new Integer(2).equals(out))
                throw new RuntimeException("bad JS result");
        }
        System.out.print("JS execution took: ");
        System.out.print(System.currentTimeMillis() - startTime);
        System.out.println(" milliseconds");
    }

    private static void
    testPython() {
        final Map<String, Object> sourceJSON = new HashMap<>();
        sourceJSON.put("code", "def run(): return 2");
        sourceJSON.put("name", "tester");
        sourceJSON.put("language", "python2");
        final List<Object> args = new ArrayList<>();
        final long startTime = System.currentTimeMillis();
        for (int i = 0; i < 1000; i++) {
            final Object out =
                SubprocessExec.run(null, sourceJSON, args, noopHandler);
            if (!new Integer(2).equals(out))
                throw new RuntimeException("bad Python result");
        }
        System.out.print("Python execution took: ");
        System.out.print(System.currentTimeMillis() - startTime);
        System.out.println(" milliseconds");
    }
}
