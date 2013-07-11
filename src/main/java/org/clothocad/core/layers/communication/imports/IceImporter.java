/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication.imports;

import com.gargoylesoftware.htmlunit.BrowserVersion;
import com.gargoylesoftware.htmlunit.FailingHttpStatusCodeException;
import com.gargoylesoftware.htmlunit.HttpMethod;
import com.gargoylesoftware.htmlunit.Page;
import com.gargoylesoftware.htmlunit.WebClient;
import com.gargoylesoftware.htmlunit.WebRequest;
import com.gargoylesoftware.htmlunit.WebResponse;
import com.gargoylesoftware.htmlunit.html.HtmlButton;
import com.gargoylesoftware.htmlunit.html.HtmlPage;
import com.gargoylesoftware.htmlunit.html.HtmlPasswordInput;
import com.gargoylesoftware.htmlunit.html.HtmlTextInput;
import com.gargoylesoftware.htmlunit.util.NameValuePair;
import com.gargoylesoftware.htmlunit.xml.XmlPage;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.bson.types.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.JSON;
import org.clothocad.model.FreeForm;
import org.clothocad.model.NucSeq;
import org.jbei.ice.lib.entry.model.Entry;
import org.jbei.ice.lib.entry.model.Part;
import org.jbei.ice.lib.models.Sequence;
import org.jbei.ice.lib.utils.IceXmlSerializer;
import org.jbei.ice.lib.utils.UtilityException;
import org.jbei.ice.lib.vo.CompleteEntry;
import org.python.google.common.base.Joiner;

/**
 *
 * @author spaige
 */
public class IceImporter {

    static final IceXmlSerializer serializer = new IceXmlSerializer();
    static final String protocol = "https://";
    static final String host = "public-registry.jbei.org/";
    static final String worked = "https://public-registry.jbei.org/export?type=XML&entries=745";
    
    
    static final String usernamePath = "//*[@id=\"user_login_input\"]/input";
    static final String pwdPath = "//*[@id=\"user_password_input\"]/input";
    static final String loginBtnPath = "//*[@id=\"login_button\"]/button";
    static final String searchBarPath = "//*[@id=\"search_panel\"]/table/tbody/tr/td[1]/table/tbody/tr/td[2]/input";

    public static String convertStreamToString(java.io.InputStream is) {
        java.util.Scanner s = new java.util.Scanner(is).useDelimiter("\\A");
        return s.hasNext() ? s.next() : "";
    }

    public static void doImport(Persistor p, List<Integer> entries, int attempt) {
        String file = String.format("export?type=XML&entries=%s", Joiner.on(", ").join(entries));

        //navigate to jbei
        WebClient client = new WebClient(BrowserVersion.FIREFOX_17);
        client.waitForBackgroundJavaScriptStartingBefore(1000);
        try {
            //login
            HtmlPage page = client.getPage("https://public-registry.jbei.org");
            HtmlTextInput uname = page.getFirstByXPath(usernamePath);
            uname.setValueAttribute("stmpaige@gmail.com");
            HtmlPasswordInput pwd = page.getFirstByXPath(pwdPath);
            pwd.setValueAttribute("blue j 5");
            HtmlButton button = page.getFirstByXPath(loginBtnPath);
            page = button.click();
            Thread.sleep(5000);
            
            String address = protocol+host+file;
            
            XmlPage export = client.getPage(address);
            String xml = export.getContent();


            //parse entries
            List<CompleteEntry> jbeiEntries = serializer.deserializeJbeiXml(xml);

            //parse classes to BSON
            //recordId -> uuid


            for (CompleteEntry entry : jbeiEntries) {
                Entry realEntry = entry.getEntry();
                Map<String, Object> data;
                if (realEntry instanceof Part){
                    Sequence sequence = entry.getSequence();
                    NucSeq seq = new NucSeq(sequence.getSequence(), true, false);
                    
                    org.clothocad.model.BasicPart clothoPart = new org.clothocad.model.BasicPart(realEntry.getOneName().toString(),
                            realEntry.getShortDescription(), entry.getSequence().getSequence(), new FreeForm(), null);
                    //clothoPart.setUUID(new ObjectId(Joiner.on("").join(realEntry.getRecordId().split("-"))));
                    p.save(clothoPart);
                }
            }



        } catch (MalformedURLException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        } catch (FailingHttpStatusCodeException ex) {
            ex.printStackTrace();
        } catch (UtilityException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

    }
}
