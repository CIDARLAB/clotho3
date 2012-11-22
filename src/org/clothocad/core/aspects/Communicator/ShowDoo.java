/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Communicator;

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Sharable;
import org.json.JSONArray;

/**
 *
 * @author Tao
 */
public class ShowDoo extends Doo {
    
    public ShowDoo(Doo doo) {
        super(doo, false);
    }
    
    public String widgetId;
    public JSONArray commandMessageArray;
    public String viewId;
    public List<String> sharableIds = new ArrayList<String>();

    void collectSharables(List<Sharable> shareList) {
        for(Sharable sharable : shareList) {
            if(sharable!=null) {
                sharableIds.add(sharable.getId());
            }
        }
    }
    

}
