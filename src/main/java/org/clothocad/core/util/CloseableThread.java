package org.clothocad.core.util;

import java.io.InputStream;

public class CloseableThread implements AutoCloseable {
    private final Thread thread;

    public
    CloseableThread(final Runnable runnable) {
        thread = new Thread(runnable);
    }

    public Thread
    getThread() {return thread;}

    @Override public void
    close() {
        thread.interrupt();
        try {
            thread.join();
        } catch (InterruptedException e) {
            /* give up joining */
            thread.interrupt();
        }
    }
}
