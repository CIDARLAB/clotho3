package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.ResourceConnection;
import org.clothocad.core.communication.Router;

import java.io.InputStream;
import java.io.IOException;
import java.io.OutputStream;

import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.ServletException;

/**
 *
 * @author billcao
 */
public class ResourceApi extends HttpServlet {

    private static Router router;
    private static Message m;
    private static final ResourceConnection rc = new ResourceConnection("ResourceConnection");

    public ResourceApi(Router router) {
        ResourceApi.router = router;
    }

    // TODO: For directory paths, serve directory path + index.html
    @Override
    protected void doGet(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {

        // Allows GET calls from domains other than this Clotho's domain
        response.addHeader("Access-Control-Allow-Origin", "*");

        String[] pathID = request.getPathInfo().split("/");
        String resourcePath = request.getPathInfo();
        String fileName = pathID[pathID.length - 1];

        m = new Message(Channel.resource, resourcePath, null, null);

        ResourceApi.router.receiveMessage(ResourceApi.rc, m);

        InputStream fileStream = ResourceApi.rc.getResult();

        if (fileStream == null) {
            response.setStatus(HttpServletResponse.SC_NOT_FOUND);
            return;
        }

        OutputStream responseOutputStream = response.getOutputStream();

        response.setStatus(HttpServletResponse.SC_OK);
        response.setHeader("Content-Disposition", "inline; filename='" + fileName + "'");

        int bytes;
        while ((bytes = fileStream.read()) != -1) {
            responseOutputStream.write(bytes);
        }
    }
}
