/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import java.util.ArrayList;
import java.util.List;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Sharable;

/**
 *
 * @author Tao
 */
public class ShowDoo extends Doo {
    
    public ShowDoo(Doo doo) {
        super(doo, false);
    }
    
    public String widgetId;
    public List commandMessageArray;
    public String viewId;
    public List<ObjectId> sharableIds = new ArrayList<>();

    void collectSharables(List<Sharable> shareList) {
        for(Sharable sharable : shareList) {
            if(sharable!=null) {
                sharableIds.add(sharable.getId());
            }
        }
    }
    

}
