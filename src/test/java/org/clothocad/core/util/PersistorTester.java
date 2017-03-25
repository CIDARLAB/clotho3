/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;
import org.junit.Test;

/**
 *
 * @author David
 */
public class PersistorTester extends AuthorizedShiroTest {
    
    private Persistor persistor;

    public PersistorTester() {
        this.persistor = injector.getInstance(Persistor.class);
    }

    @Test
    public void timeToBulkCreate()
    {
        System.out.println("Testing Bulk Create");
        Person roxanne = new Person("Roxanne Shank");       
        long start = System.currentTimeMillis();
        for (int i = 0; i < 5000; i++)
        {           

            Sequence seqK249002 = new Sequence("K249" + i + " Sequence", "atgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgctttaaccatgcgc"
				+ "tggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaaattccggtgg"
				+ "cggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagcgaagtgagca"
				+ "aaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcaccaaacatggca"
				+ "acaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctag", roxanne);
            persistor.save(seqK249002);
        }
        long end = System.currentTimeMillis();
        System.out.println("Bulk Create in Persistor took " + (end - start) + " MilliSeconds");      
    }
}
