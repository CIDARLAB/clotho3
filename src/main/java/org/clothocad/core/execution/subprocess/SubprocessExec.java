package org.clothocad.core.execution.subprocess;

import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.util.CloseableProcess;
import org.clothocad.core.util.CloseableThread;

/** Executes function in an external language interpreter subprocess
 *
 * Communication with the subprocess is done over pipes connected to
 * the standard input/output of the subprocess. Messages are encoded
 * as null-byte terminated, UTF-8 JSON values.
 *
 * The types of messages that can be sent to the subprocess are:
 *   function_defcall ("define and call"):
 *     {"type": "func", "code": <string>,
 *      "name": <string>, "args": <array>}
 *
 *   api_return: (send only)
 *     {"type": "api", "return": <value>}
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
 *  Host             Subprocess
 *main thread
 *  ...
 *   |                 [born]
 *   |                   |
 *   |                   X (blocks)
 *   | function_defcall
 *   | ----------------> \ (resumes)
 *   X                   |
 *                       |
 *                       |
 *         api_call      |
 *   / <---------------- |
 *   |                   X
 *   |     api_return
 *   | ----------------> \
 *   X                   |
 *                       |
 *  ...                 ...
 *                       |
 *      function_return  |
 *   / <---------------- |
 *   |                 [dies]
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
            } catch (Exception e) {
                eventHandler.onFail(errorDumper.getString());
                throw e;
            }
            eventHandler.onSuccess(errorDumper.getString());
            return returnValue;
        }
    }

    public static interface EventHandler {
        void onFail(final String standardError);
        void onSuccess(final String standardError);
    }
}
