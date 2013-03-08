/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.jetty;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.FileNotFoundException;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;

/* TODO: special-case the root path "/" */
class StaticServlet extends HttpServlet {
    @Override
    public void doGet(HttpServletRequest request,
                      HttpServletResponse response) {
        String resource_path = getResourcePath(request, response);
        if (resource_path == null)
            return;

        byte[] contents = getContents(request, response, resource_path);
        if (contents == null)
            return;

        sendResponse(response, contents, resource_path);
    }

    private String getResourcePath(HttpServletRequest request,
                                   HttpServletResponse response) {
        String resource_path = request.getRequestURI();
        if (resource_path == null || !resource_path.startsWith("/")) {
            doError(response);
            return null;
        }

        /* Special-case root path */
        if ("/".equals(resource_path)) {
            response.setStatus(HttpServletResponse.SC_SEE_OTHER);
            response.setHeader("Location", "/servlet/page?HOMEPAGE");
            try {
                response.flushBuffer();
            } catch (IOException e) {
                Logger.log(Logger.Level.FATAL, "Cannot send response", e);
            }
            return null;
        }

        return "./web" + resource_path;
    }

    private byte[] getContents(HttpServletRequest request,
                               HttpServletResponse response,
                               String resource_path) {
        File f = new File(resource_path);
        if (!f.isFile()) {
            doError(response);
            return null;
        }

        FileInputStream fis;
        try {
            fis = new FileInputStream(f);
        } catch (FileNotFoundException e) {
            Logger.log(Logger.Level.WARN, "cannot find " + resource_path);
            return null;
        }
        byte[] contents;
        try {
            contents = FileUtils.dumpInputStream(fis);
        } catch (IOException e) {
            doError(response);
            return null;
        }

        return contents;
    }

    private void sendResponse(HttpServletResponse response,
                              byte[] contents,
                              String resource_path) {
        response.setHeader("Content-Type", getMimeFromExt(resource_path));
        response.setStatus(HttpServletResponse.SC_OK);
        try {
            response.getOutputStream().write(contents);
        } catch (IOException e) {
            Logger.log(Logger.Level.FATAL, "cannot send response", e);
        }
    }

    private String getMimeFromExt(String path) {
             if (path.endsWith(".js"))  return "application/javascript";
        else if (path.endsWith(".css")) return "text/css";
        else if (path.endsWith(".png")) return "image/png";
        else if (path.endsWith(".jpg")) return "image/jpeg";
        else if (path.endsWith(".gif")) return "image/gif";
        else                            return "text/html";
    }

    private void doError(HttpServletResponse response) {
        try {
            response.sendError(HttpServletResponse.SC_NOT_FOUND);
        } catch (IOException e) {
            Logger.log(Logger.Level.FATAL, "cannot send response", e);
        }
    }
}
