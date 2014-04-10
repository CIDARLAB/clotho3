package org.clothocad.core.execution.subprocess;

import java.io.InputStream;
import java.io.IOException;

class ErrorDumper implements Runnable {
    private final InputStream errorStream;
    private final StringBuffer buffer = new StringBuffer(0);

    ErrorDumper(final InputStream errorStream) {
        this.errorStream = errorStream;
    }

    String getString() {
        return buffer.toString();
    }

    @Override public void
    run() {
        while (!Thread.interrupted()) {
            final int c;
            try {
                c = errorStream.read();
            } catch (IOException e) {
                break;
            }
            if (c < 0)
                break;
            buffer.append((char) c);
        }
    }
}
