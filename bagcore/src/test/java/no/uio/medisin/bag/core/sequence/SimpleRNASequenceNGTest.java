/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;

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
public class SimpleRNASequenceNGTest {
    
    public SimpleRNASequenceNGTest() {
    }

    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}
    
    @Test
    public static void testConstructors(){
        SimpleRNASequence srseq1 = new SimpleRNASequence();
        Assert.assertEquals(0, srseq1.getLength());
        Assert.assertEquals("", srseq1.getName());
        Assert.assertNull(srseq1.getAccessionNumber());
        
        
        SimpleRNASequence srseq2 = new SimpleRNASequence("ID", "ACGT");
        Assert.assertEquals("ID", srseq2.getId());
        Assert.assertEquals("ID", srseq2.getName());
        Assert.assertEquals(4, srseq2.getLength());
        Assert.assertEquals("ACGT", srseq2.getSeq());
        srseq2.setSeq("ACCGGGTTTT");
        Assert.assertEquals(10, srseq2.getLength());
        srseq2.setStartPos(21);
        Assert.assertEquals(21, srseq2.getStartPos());
        Assert.assertEquals(30, srseq2.getEndPos());

        srseq2.setEndPos(40);
        srseq2.setStartPos(srseq2.getEndPos()-srseq2.getLength()+1);
        Assert.assertEquals(31, srseq2.getStartPos());
        
        srseq2.setChromosome("1");
        
                
        SimpleRNASequence srseq3 = new SimpleRNASequence(srseq2);
        Assert.assertEquals("ID", srseq3.getId());
        Assert.assertEquals("ID", srseq3.getName());
        Assert.assertEquals(10, srseq3.getLength());
        Assert.assertEquals("ACCGGGTTTT", srseq3.getSeq());
        Assert.assertEquals("1", srseq3.getChromosome());
        

        
        
        SimpleSeq sseq1 = new SimpleSeq("ID", "ACGT");
        sseq1.setChromosome("1");
        sseq1.setName("name3");
        sseq1.setID("ID3");
        sseq1.setStrand("+");

        
        
        SimpleRNASequence srseq4 = new SimpleRNASequence(sseq1);
        Assert.assertEquals("name3", srseq4.getName());
        Assert.assertEquals("ID3", srseq4.getId());
        Assert.assertEquals("1", srseq4.getChromosome());
        Assert.assertEquals("+", srseq4.getStrand());

        
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
