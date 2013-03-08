package org.clothocad.core.jetty;

import java.io.IOException;
import java.net.URL;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.apache.catalina.ssi.SSIServlet;
import org.clothocad.core.layers.communication.mind.PageMode;
import org.clothocad.core.settings.Settings;

class PageServlet extends SSIServlet {
    @Override
    public void doGet(HttpServletRequest request,
                      HttpServletResponse response)
                     throws IOException, ServletException {
        String mode_name = request.getQueryString();
        if (mode_name == null) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND);
            return;
        }

        PageMode mode;
        try {
            mode = Enum.valueOf(PageMode.class, mode_name);
        } catch (IllegalArgumentException e) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND);
            return;
        }

        processSSI(request,
                   response,
                   new URL(Settings.getRootURL() + "/pages/" + mode.name() + ".shtml"));
    }

    String getPath() {
        return "/servlet/page";
    }
}
