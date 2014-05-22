package org.clothocad.core.execution.subprocess;

import java.io.InputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

class CloseableTmpFile implements AutoCloseable {
    private Path path;
    private InputStream instream;

    CloseableTmpFile() {}

    InputStream
    getInputStream() {return instream;}

    Path
    setup() {
        if (path != null || instream != null)
            throw new IllegalStateException();
        try {
            path = Files.createTempFile(null, ".clothotmp");
            instream = Files.newInputStream(path);
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
        return path;
    }

    @Override public void
    close() {
        try {
            instream.close();
        } catch (IOException e) {
        }
        try {
            /* we don't care whether file really existed; but we avoid
             * Files.delete(...) because it raises more exceptions and we
             * don't want to die here.
             */
            Files.deleteIfExists(path);
        } catch (IOException | RuntimeException e) {
        }
    }
}
