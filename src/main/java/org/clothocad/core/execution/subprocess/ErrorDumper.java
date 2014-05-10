package org.clothocad.core.execution.subprocess;

import java.io.InputStream;
import java.io.IOException;
import org.clothocad.core.util.ByteArray;

class ErrorDumper implements Runnable {
    private final InputStream errorStream;
    private final ByteArray buffer = new ByteArray();

    ErrorDumper(final InputStream errorStream) {
        this.errorStream = errorStream;
    }

    byte[] getBytes() {
        return buffer.getArray();
    }

    @Override public void
    run() {
        while (!Thread.interrupted()) {
            final int b;
            try {
                b = errorStream.read();
            } catch (IOException e) {
                break;
            }
            if (b < 0)
                break;
            buffer.add((byte) b);
        }
    }
}
