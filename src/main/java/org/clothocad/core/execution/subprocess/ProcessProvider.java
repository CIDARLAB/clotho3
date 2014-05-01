package org.clothocad.core.execution.subprocess;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.util.CloseableProcess;

class ProcessProvider {
    private static final Map<String, List<String>> COMMANDS = new HashMap<>();
    static {
        COMMANDS.put(
            "python2",
            Arrays.asList(new String[] {
                "python",
                Paths.get("src", "main", "python", "runner.py").toString()
            })
        );
    }

    static CloseableProcess newProcess(final String language) {
        final List<String> cmd = COMMANDS.get(language);
        if (cmd == null)
            throw new IllegalArgumentException(language);
        try {
            return new CloseableProcess(new ProcessBuilder(cmd).start());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
