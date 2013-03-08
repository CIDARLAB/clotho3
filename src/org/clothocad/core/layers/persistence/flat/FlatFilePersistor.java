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

package org.clothocad.core.layers.persistence.flat;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;

import com.thoughtworks.xstream.XStream;

/**
 * @author jcanderson
 */
public class FlatFilePersistor 
		extends Persistor {
	
	public boolean persist(Collection<ObjBase> col) {
		for(Datum datum : col) {
			if(!this.persistDatum(datum)) {
				return false;
			}
		}
		return true;
	}
	
	public boolean persist(ObjBase obj) {
		// here, we forward the datum to the MongoDBPersistor
		return obj.save();
	}

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

    public boolean persistDatum (Datum datum) {
        String path = FileUtils.getFilePath(datum.getId(), dir + "Datum");
        try {
            String xml = xstream.toXML(datum);
            FileUtils.writeFile(xml, path);
            return true;
        } catch (Exception err) {
            err.printStackTrace();
        }
        return false;
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

    public static final String dir = "database/";
    private XStream xstream = new XStream();
	
    
	public Datum get(Class c, ObjBase obj) {
        throw new UnsupportedOperationException("Not supported yet.");
	}

	public Datum get(Class c, ObjectId id) {
        throw new UnsupportedOperationException("Not supported yet.");
	}
    
	@Override
	public ObjBase[] get(ObjBase obj) {
        throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public <T> T[] get(Class<T> t) {
        throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean clearDB() {
		// TODO Auto-generated method stub
		return false;
	}
}
	
