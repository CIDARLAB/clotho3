/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.io.IOException;
import java.util.List;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import static org.junit.Assert.assertEquals;
import org.junit.Before;

/**
 *
 * @author spaige
 */
public class AbstractRouterTest extends AuthorizedShiroTest {
    protected static Router router;
    protected static List<ObjectId> ids;

    public AbstractRouterTest() {
        router = injector.getInstance(Router.class);

    }

    @Before
    public void setUp() {
        injector.getInstance(Persistor.class).deleteAll();
        ids = TestUtils.setupTestData(injector.getInstance(Persistor.class), injector.getInstance(ClothoRealm.class));
    }
    
    static void assertMatch(Message m1, Message m2) {
        assertEquals(m1.getChannel(), m2.getChannel());
        assertEquals(m1.getRequestId(), m2.getRequestId());
    }
    
    protected void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
    
}
