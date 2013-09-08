/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.utils;

import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.layers.communication.Message;
import org.clothocad.core.layers.communication.ServerSideAPI;
import org.clothocad.core.persistence.Persistor;

/**
 *
 * @author spaige
 */
public class DummyAPI extends ServerSideAPI {

    public DummyAPI(Persistor persistor) {
        super(null, persistor, null);
    }

    @Override
    protected void say(String message, Severity severity, String recipients, boolean isUser) {
    }

    @Override
    protected void send(Message message) {
    }

    
    
}
