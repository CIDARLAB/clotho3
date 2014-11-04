/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.util.List;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.core.util.TestUtils;
import org.junit.Before;

/**
 *
 * @author spaige
 */
public class AbstractServerAPITest extends AuthorizedShiroTest {
    protected static ServerSideAPI api;
    protected static ServerSideAPI unprivilegedUser;
    protected static Persistor persistor;
    protected static List<ObjectId> ids;
    protected static Mind mind;
    protected static Router router;
    protected static ClothoRealm realm;

    public AbstractServerAPITest() {
        super();
        persistor = injector.getInstance(Persistor.class);
        router = injector.getInstance(Router.class);
        realm  = injector.getInstance(ClothoRealm.class);
        mind = new Mind();
        api = new ServerSideAPI(mind, persistor, router, realm, null);
        mind.setConnection(new TestConnection("test"));
    }
    
    @Before
    public void setUp() {
        persistor.deleteAll();
        ids = TestUtils.setupTestData(persistor, realm);
    }
}
