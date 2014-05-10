package org.clothocad.core.execution.subprocess;

import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.util.CloseableProcess;
import org.clothocad.core.util.CloseableThread;

/** Executes function in an external language interpreter subprocess
 *
 * Communication with the subprocess is done over pipes connected to
 * the standard input/output/error of the subprocess. Messages over standard
 * input/output are UTF-8 encoded, null-byte terminated JSON values. Standard
 * error is captured in a byte array and is not interpreted.
 *
 * The types of messages that can be sent to the subprocess are:
 *   function_defcall ("define and call"):
 *     {"type": "func", "code": <string>, "args": <array>}
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
 *   |    function_defcall
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
 * standard error streams of the subprocess. When the subprocess dies, so do
 * the helper threads. The diagram above does not show the helper threads.
 */
public class SubprocessExec {
    public static Object
    run(final ServerSideAPI api,
        final Map<String, Object> sourceJSON,
        final List<Object> args,
        final EventHandler eventHandler) {

        try (final CloseableProcess proc =
             ProcessProvider.newProcess((String) sourceJSON.get("language"))
        ) {
            final JSONStreamReader reader = new JSONStreamReader(
                proc.getProcess().getInputStream()
            );
            final ErrorDumper errorDumper =
                new ErrorDumper(proc.getProcess().getErrorStream());
            final Object returnValue;
            try {
                try (final CloseableThread readerThread =
                     new CloseableThread(reader);
                     final CloseableThread errorDumperThread =
                     new CloseableThread(errorDumper)
                ) {
                    readerThread.getThread().start();
                    errorDumperThread.getThread().start();
                    returnValue = new ExecutionContext(
                        api,
                        reader,
                        proc.getProcess().getOutputStream(),
                        (String) sourceJSON.get("code"),
                        args
                    ).start();
                }
            } catch (final Exception e) {
                eventHandler.onFail(errorDumper.getBytes());
                throw e;
            }
            eventHandler.onSuccess(errorDumper.getBytes());
            return returnValue;
        }
    }

    public static interface EventHandler {
        void onFail(final byte[] standardError);
        void onSuccess(final byte[] standardError);
    }
}
