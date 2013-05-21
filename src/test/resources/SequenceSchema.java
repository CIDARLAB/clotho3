package org.clothocad.core.testers.schemas

import lombok.Getter;
import lombok.Setter;
/**
 * A simple and sloppy representation of a Feature or other DNA sequence
 */
public class SimpleFeature {

    @Getter
    @Setter
    /**
     * The sequence of the feature
     */
    @Example("ATACCGGA")
    private String sequence;

}

