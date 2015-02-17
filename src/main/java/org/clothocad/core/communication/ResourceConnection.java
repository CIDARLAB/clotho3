package org.clothocad.core.communication;

import java.io.IOException;
import java.io.InputStream;

/**
 *
 * @author billcao
 */
public class ResourceConnection extends ClientConnection {

    private boolean done = false;
    private InputStream fileStream;

    public ResourceConnection(String id) {
        super(id);
    }

    @Override
    public void send(Message msg) {
        try {
            fileStream = (InputStream) msg.getData();
            this.done = true;
        } catch (Exception e) {
            System.out.println("Got here (failed to find file, exception)");
            fileStream = null;
        }
    }

    public InputStream getResult() {
        return fileStream;
    }
}
