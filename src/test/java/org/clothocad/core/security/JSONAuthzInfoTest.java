package org.clothocad.core.security;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.TestConnection;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Institution;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class JSONAuthzInfoTest extends AbstractSecurityTest {

    private static String permissions = "$$permissions";
    protected Router router;

    protected Message getMessage() {
        ObjectId id = util.getPrivate().getId();
        return new Message(Channel.get, id, "");
    }

    protected Message findMessage() {
        Map<String, Object> query = new HashMap<>();
        query.put("schema", Institution.class.getCanonicalName());

        return new Message(Channel.query, query, "");
    }

    @Override
    public void beforeTest() {
        super.beforeTest();
        router = injector.getInstance(Router.class);
        initAPI("");
    }

    protected void assertMyPermissions(Set<ClothoAction> actions, Map<String, Object> permissionData) {
        assertNotNull(permissionData);

        //contains 'mine' permissions
        Set<ClothoAction> myPermissions = (Set) permissionData.get("mine");
        assertNotNull(myPermissions);

        //contains correct 'mine' permissions
        assertEquals(actions, myPermissions);
    }

    protected void assertNoAdvancedPermissions(Map<String, Object> permissionData) {
        //does not contain 'user' or 'group' permissions        
        assertNull(permissionData.get("user"));
        assertNull(permissionData.get("group"));
    }

    protected void assertGroupPermissionsExist(Map<String,Object> permissionData){
        assertNotNull(permissionData.get("group"));
    }
    
    protected void assertUserPermissions(Map<String,Object> permissionData){
        Map<String,Object> userPermissionData = (Map) permissionData.get("user");
        assertNotNull(userPermissionData);
        assertEquals(userPermissionData.get("owner"), ClothoPermission.OWN.actions);
        assertEquals(userPermissionData.get("write"), ClothoPermission.WRITE.actions);
    }
    
    @Test
    public void testGet() throws IOException {
        login("write");
        TestConnection connection = new TestConnection("");
        sendMessage(getMessage(), connection);
        Map<String, Object> result = (Map) connection.messageDataByChannelAndId.get("get");

        Map<String, Object> permissionData = (Map) result.get(permissions);

        assertMyPermissions(ClothoPermission.WRITE.actions, permissionData);
        assertNoAdvancedPermissions(permissionData);
    }

    @Test
    public void testFind() throws IOException {
        login("write");
        TestConnection connection = new TestConnection("");
        sendMessage(findMessage(), connection);
        List<Map<String, Object>> results = (List) connection.messageDataByChannelAndId.get("query");

        for (Map<String, Object> result : results) {
            Map<String, Object> permissionData = (Map) result.get(permissions);

            assertMyPermissions(ClothoPermission.WRITE.actions, permissionData);
        }
    }

    @Test
    public void testOwnerGet() throws IOException {
        login("owner");
        TestConnection connection = new TestConnection("");
        sendMessage(getMessage(), connection);
        Map<String, Object> result = (Map) connection.messageDataByChannelAndId.get("get");

        Map<String, Object> permissionData = (Map) result.get(permissions);

        assertMyPermissions(ClothoPermission.OWN.actions, permissionData);
        assertGroupPermissionsExist(permissionData);
        assertUserPermissions(permissionData);
    }

    @Test
    public void testOwnerQuery() throws IOException {
        login("owner");
        TestConnection connection = new TestConnection("");
        sendMessage(findMessage(), connection);
        List<Map<String, Object>> results = (List) connection.messageDataByChannelAndId.get("query");

        for (Map<String, Object> result : results) {
            Map<String, Object> permissionData = (Map) result.get(permissions);

            assertMyPermissions(ClothoPermission.OWN.actions, permissionData);
            assertGroupPermissionsExist(permissionData);
            assertUserPermissions(permissionData);
        }
    }

    protected void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
}
