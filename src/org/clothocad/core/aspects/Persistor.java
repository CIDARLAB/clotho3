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

package org.clothocad.core.aspects;

import java.io.File;
import java.util.HashMap;
import com.thoughtworks.xstream.XStream;
import java.util.List;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.Datum;

/**
 * @author jcanderson
 */
public class Persistor implements Aspect {
    public Datum loadDatum (String uuid) {
        String xml = FileUtils.readFile(dir + "Datum/" + uuid);
        if (xml.equals("")) {
            Logger.log(Logger.Level.WARN,
                       "There is no data with that uuid");
            return null;
        } else {
            return (Datum) xstream.fromXML(xml);
        }
    }

    public HashMap<String, Integer> loadFeature (String feature) {
        String xml = FileUtils.readFile(dir + "FeaturesDB/" + feature);
        if (xml == "") {
            return null;
        } else {
            @SuppressWarnings("unchecked")
            HashMap<String, Integer> out = (HashMap<String, Integer>)
                                                    xstream.fromXML(xml);
            return out;
        }
    }

    public List<String> loadWordBank () {
        String xml = FileUtils.readFile(dir + "WordBank/words");
        if (xml == "") {
            Logger.log(Logger.Level.WARN, 
                        "No word bank existing.");
            return null;
        } else {
            return (List<String>) xstream.fromXML(xml);
        }
    }

    public void persistDatum (Datum obj) {
        String path = FileUtils.getFilePath(obj.getId(), dir + "Datum");
        try {
            String xml = xstream.toXML(obj);
            FileUtils.writeFile(xml, path);
        } catch (Exception err) {
            err.printStackTrace();
        }
    }

    public void persistFeature (Object obj, String filePath) {
        String path = FileUtils.getFilePath(filePath, dir + "FeaturesDB");
        try {
            String xml = xstream.toXML(obj);
            FileUtils.writeFile(xml, path);
        } catch (Exception err) {
            err.printStackTrace();
        }
    }

    public void persistWordBank (List<String> wordBank) {
        String path = FileUtils.getFilePath("words", dir + "WordBank/");
        try {
            String xml = xstream.toXML(wordBank);
            FileUtils.writeFile(xml, path);
        } catch (Exception err) {
            err.printStackTrace();
        }
    }

    public static Persistor get() {
        return singleton;
    }
    
    private static final Persistor singleton = new Persistor();
    public static final String dir = "database/";
    private XStream xstream = new XStream();
}
