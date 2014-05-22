package org.clothocad.core.execution.subprocess;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.util.CloseableProcess;
import org.clothocad.core.util.CloseableRunnable;
import org.clothocad.core.util.CloseableThread;

/** Executes function in an external language interpreter subprocess
 *
 * Our Java process (the "host") communicates with the subprocess over
 * the following channels:
 *   - subprocess standard input:  host --> subprocess messages
 *   - temporary file:             subprocess --> host messages
 *   - subprocess standard output: copied to host; not interpreted
 *   - subprocess standard error:  copied to host; not interpreted
 *
 * The temporary file is created by the host and the file name is passed
 * to the subprocess via its standard input. This enables communication
 * from subprocess to host.
 *
 * Subprocess standard output and error is captured in byte arrays and
 * is not interpreted.
 *
 * Messages are UTF-8 encoded, null-byte terminated JSON values.
 * The types of messages that can be sent to the subprocess are:
 *   function_init:
 *     {"type": "func", "tmpfile": <string>, "code": <string>, "args": <array>}
 *
 *   api_return:
 *     {"type": "api", "return": <value>}
 *
 *   api_error:
 *     {"type": "api_error", "message": <string>}
 *
 * The messages that can be received from the subprocess are:
 *   function_return:
 *     {"type": "func", "return": <value>}
 *
 *   api_call:
 *     {"type": "api", "name": <string>, "args": <list>}
 *
 * The message protocol and threads of execution are
 * summarized in the following diagram:
 *
 *  Host                    Subprocess
 *main thread
 *  ...
 *   |                       [is born]
 *   |                          |
 *   |                          X (blocks)
 *   |      function_init
 *   | -----------------------> \ (resumes)
 *   X                          |
 *                              |
 *  ...                        ...
 *                              |
 *            api_call          |
 *   / <----------------------- |
 *   |                          X
 *   | api_return or api_error
 *   | -----------------------> \
 *   X                          |
 *                              |
 *  ...                        ...
 *                              |
 *         function_return      |
 *   / <----------------------- |
 *   |                        [dies]
 *  ...
 *
 * The host creates helper threads for reading from the standard output and
 * standard error streams of the subprocess and from the temporary file. When
 * the subprocess dies, so do the helper threads. The temporary file is also
 * deleted. The diagram above does not show the helper threads.
 */
public class SubprocessExec {
    public static Object
    run(final ServerSideAPI api,
        final Map<String, Object> sourceJSON,
        final List<Object> args,
        final EventHandler eventHandler) {

        try (final CloseableProcess proc =
             ProcessProvider.newProcess((String) sourceJSON.get("language"));
             final CloseableTmpFile tmpfile = new CloseableTmpFile()
        ) {
            final Object returnValue;
            final Condition procCond = new Condition() {
                @Override public boolean poll() {
                    try {
                        proc.getProcess().exitValue();
                    } catch (Exception e) {
                        return true;
                    }
                    return false;
                }
            };
            final Path tmpfile_path = tmpfile.setup();
            final PipeDumper outDumper =
                new PipeDumper(proc.getProcess().getInputStream(), procCond);
            final PipeDumper errDumper =
                new PipeDumper(proc.getProcess().getErrorStream(), procCond);
            final JSONStreamReader reader =
                new JSONStreamReader(tmpfile.getInputStream(), procCond);
            try {
                returnValue = initThreadsAndRun(
                    api,
                    reader,
                    new JSONStreamWriter(proc.getProcess().getOutputStream()),
                    new CloseableRunnable[] {reader, outDumper, errDumper},
                    tmpfile_path,
                    (String) sourceJSON.get("code"),
                    args
                );
            } catch (final Exception e) {
                eventHandler.onFail(outDumper.getBytes(),
                                    errDumper.getBytes());
                throw e;
            }
            eventHandler.onSuccess(outDumper.getBytes(), errDumper.getBytes());
            return returnValue;
        }
    }

    private static Object
    initThreadsAndRun(final ServerSideAPI api,
                      final JSONStreamReader r,
                      final JSONStreamWriter w,
                      final CloseableRunnable[] runnables,
                      final Path tmpfile_path,
                      final String code,
                      final List<Object> args) {
        final CloseableThread[] threads =
            new CloseableThread[runnables.length];

        for (int i = 0; i < runnables.length; i++)
            threads[i] = new CloseableThread(runnables[i]);

        try (final TryAll tryall = new TryAll(threads)) {
            for (final CloseableThread thread : threads)
                thread.getThread().start();
            return new ExecutionContext(api, r, w)
                   .start(tmpfile_path, code, args);
        }
    }

    public static interface EventHandler {
        void onFail(byte[] standardOutput, byte[] standardError);
        void onSuccess(byte[] standardOutput, byte[] standardError);
    }
}
