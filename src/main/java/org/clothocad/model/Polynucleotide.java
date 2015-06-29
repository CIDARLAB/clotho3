package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.io.Serializable;
import java.util.Date;

import javax.validation.constraints.NotNull;

@NoArgsConstructor
public class Polynucleotide extends SharableObjBase implements Serializable {

    @Getter
    @Setter
    protected String accession;

    @NotNull
    @Getter
    @Setter
    protected boolean isLinear, isSingleStranded;

    @Getter
    @Setter
    protected Date submissionDate;

    @Getter
    @Setter
    @Reference
    protected Sequence sequence;

    @Getter
    @Setter
    @Reference
    protected Polynucleotide parentPolynucleotide;

    public Polynucleotide(String name, boolean isLinear, boolean isSingleStranded, Person author) {
        super(name, author);
        this.isLinear = isLinear;
        this.isSingleStranded = isSingleStranded;
    }

    public Polynucleotide(String name, String description, boolean isLinear, boolean isSingleStranded,
            Person author) {
        super(name, author, description);
        this.isLinear = isLinear;
        this.isSingleStranded = isSingleStranded;
    }

}
