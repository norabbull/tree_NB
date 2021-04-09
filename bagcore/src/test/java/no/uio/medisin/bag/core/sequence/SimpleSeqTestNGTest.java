/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;

import java.util.ArrayList;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author simonray
 */
public class SimpleSeqTestNGTest {
    
    public SimpleSeqTestNGTest() {
    }

    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    @Test
    public static void testConstructors(){
        
        SimpleSeq sseq1 = new SimpleSeq();
        Assert.assertEquals(0, sseq1.getLength());
        Assert.assertEquals("", sseq1.getName());
        Assert.assertNull(sseq1.getAccessionNumber());
        
        
        SimpleSeq sseq2 = new SimpleSeq("ID", "ACGT");
        Assert.assertEquals("ID", sseq2.getId());
        Assert.assertEquals("ID", sseq2.getName());
        Assert.assertEquals(4, sseq2.getLength());
        Assert.assertEquals("ACGT", sseq2.getSeq());
        
        Assert.assertEquals("ID_1-3", sseq2.splitSequence(1, 3).get(0).getName());
        Assert.assertEquals("ID_2-4", sseq2.splitSequence(1, 3).get(1).getName());
        
        
        SimpleSeq sseq3 = new SimpleSeq(sseq2);
        sseq3.setChromosome("1");
        sseq3.setName("name3");
        sseq3.setID("ID3");
        sseq3.setStrand("+");
        Assert.assertEquals("name3", sseq3.getName());
        Assert.assertEquals("ID3", sseq3.getId());
        Assert.assertEquals("1", sseq3.getChromosome());
        Assert.assertEquals("+", sseq3.getStrand());
        

        sseq3.setSeq("ACCGGGTTTT");
        Assert.assertEquals(10, sseq3.getLength());
        sseq3.setStartPos(21);
        Assert.assertEquals(21, sseq3.getStartPos());
        Assert.assertEquals(30, sseq3.getEndPos());

        sseq3.setEndPos(40);
        sseq3.setStartPos(sseq3.getEndPos()-sseq3.getLength()+1);
        Assert.assertEquals(31, sseq3.getStartPos());
                
    }
    
    @Test
    public static void testSequenceOperations(){
        SimpleSeq sseq4 = new SimpleSeq("Seq4", "ACCGGGTTTT");
        
        Assert.assertEquals(1, sseq4.NTcount('A'));
        Assert.assertEquals(0, sseq4.NTcount('a'));
        Assert.assertEquals(2, sseq4.NTcount('C'));
        Assert.assertEquals(0, sseq4.NTcount('c'));
        Assert.assertEquals(3, sseq4.NTcount('G'));
        Assert.assertEquals(0, sseq4.NTcount('g'));
        Assert.assertEquals(4, sseq4.NTcount('T'));
        Assert.assertEquals(0, sseq4.NTcount('t'));
        
        Assert.assertEquals("AAAACCCGGT", SimpleSeq.complement("ACCGGGTTTT"));
        Assert.assertEquals("ACCGGGUUUU", SimpleSeq.dna2rna("ACCGGGTTTT"));
        Assert.assertEquals("ACCGGGTTTT", SimpleSeq.rna2dna("ACCGGGUUUU"));
        
        Assert.assertEquals("AaAaCcCgGtT", SimpleSeq.complement("AaCcGgGtTtT"));
        Assert.assertEquals("AaCcGgGuUuU", SimpleSeq.dna2rna("AaCcGgGtTtT"));
        Assert.assertEquals("AaCcGgGtTtT", SimpleSeq.rna2dna("AaCcGgGtUuU"));
        
        Assert.assertEquals((double)0.50, SimpleSeq.GCfraction("AaCcGgGtTt"));

        Assert.assertEquals(2, SimpleSeq.NTcount("AaCcGgGtTt".toUpperCase(), 'A'));
        Assert.assertEquals((double)0.20, SimpleSeq.Afraction("AaCcGgGtTt"));

        Assert.assertEquals(4, SimpleSeq.NTcount("CcCcGgGtTt".toUpperCase(), 'C'));
        Assert.assertEquals((double)0.40, SimpleSeq.Cfraction("CcCcGgGtTt"));

        Assert.assertEquals(7, SimpleSeq.NTcount("GgGgGgGtTt".toUpperCase(), 'G'));
        Assert.assertEquals((double)0.70, SimpleSeq.Gfraction("GgGgGgGtTt"));

        Assert.assertEquals(6, SimpleSeq.NTcount("TtTtGgGgTt".toUpperCase(), 'T'));
        Assert.assertEquals((double)0.60, SimpleSeq.Tfraction("TtTtGgGgTt"));
        

    }
    
    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    @AfterMethod
    public void tearDownMethod() throws Exception {
    }
}
