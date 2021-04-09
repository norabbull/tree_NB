/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mapping;

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
public class SAMEntryNGTest {
    
    public SAMEntryNGTest() {
    }

    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}
    @Test
    public static void headerLineTest(){
        SAMEntry headerSAM = new SAMEntry("@SQ	SN:gi|23898|emb|X12811.1|	LN:2231");
        Assert.assertEquals(true, (boolean) headerSAM.isHeaderLine());
        
        /*
617-36  16      gi|555853|gb|U13369.1|HSU13369  1232    255     18M     *       0       0       GAGCGCACCGGACGGGCC      IIIIIIIIIIIIIIIIII      XA:i:2  MD:Z:1C9C6      NM:i:2

        */
        SAMEntry unmappedRead = new SAMEntry(
                "618-36\t4\t*\t0\t0\t"
                +   "*\t*\t0\t0\t"
                +   "TAGCTTATCAGACTGGTGTTGT\tIIIIIIIIIIIIIIIIIIIIII\tXM:i:0\n");
        Assert.assertEquals(false, (boolean) unmappedRead.isHeaderLine());
        
    }
            
    
    @Test 
    public static void testLineParsing(){
        SAMEntry mappedRead = new SAMEntry(
                "617-36\t16\tgi|555853|gb|U13369.1|HSU13369\t1232\t255\t18M\t"
                        + "*\t0\t0\tGAGCGCACCGGACGGGCC\tIIIIIIIIIIIIIIIIII\t"
                        + "XA:i:2\t MD:Z:1C9C6\t NM:i:2");
        Assert.assertEquals("617-36", mappedRead.getqName());
        Assert.assertEquals(16, mappedRead.getFlags());
        Assert.assertEquals("gi|555853|gb|U13369.1|HSU13369", mappedRead.getrName());        
        Assert.assertEquals(1232, mappedRead.getStartPos());        
        Assert.assertEquals(255, mappedRead.getMapQ());        
        Assert.assertEquals("18M", mappedRead.getCigar());        
        Assert.assertEquals("*", mappedRead.getrNext());        
        Assert.assertEquals(0, mappedRead.getpNext());        
        Assert.assertEquals(0, mappedRead.gettLen());        
        Assert.assertEquals("GAGCGCACCGGACGGGCC", mappedRead.getSeq());        
        Assert.assertEquals("IIIIIIIIIIIIIIIIII", mappedRead.getQual());        
        Assert.assertEquals("XA:i:2\t MD:Z:1C9C6\t NM:i:2", mappedRead.getBowtieTagString());        
        
        Assert.assertEquals(true, (boolean)mappedRead.isMappedRead());
        
        SAMEntry unmappedRead = new SAMEntry(
                "618-36\t4\t*\t0\t0\t"
                +   "*\t*\t0\t0\t"
                +   "TAGCTTATCAGACTGGTGTTGT\tIIIIIIIIIIIIIIIIIIIIII\tXM:i:0\n");
        Assert.assertEquals(false, (boolean)unmappedRead.isMappedRead());

        
    }
    
    
    @Test
    public static void testFlags(){
        SAMEntry mappedRead = new SAMEntry(
                "617-36\t16\tgi|555853|gb|U13369.1|HSU13369\t1232\t255\t18M\t"
                        + "*\t0\t0\tGAGCGCACCGGACGGGCC\tIIIIIIIIIIIIIIIIII\t"
                        + "XA:i:2\t MD:Z:1C9C6\t NM:i:2");
        Assert.assertEquals(true, (boolean)mappedRead.isThisFlagSet(SAMEntry.FLAG_REVCOMP));

        SAMEntry mappedReadRev = new SAMEntry(
                "617-36\t0\tgi|555853|gb|U13369.1|HSU13369\t1232\t255\t18M\t"
                        + "*\t0\t0\tGAGCGCACCGGACGGGCC\tIIIIIIIIIIIIIIIIII\t"
                        + "XA:i:2\t MD:Z:1C9C6\t NM:i:2");
        Assert.assertEquals(false, (boolean)mappedReadRev.isThisFlagSet(SAMEntry.FLAG_REVCOMP));
        
    }
    
    
    @Test
    public static void testTagExtraction(){
         SAMEntry mappedRead = new SAMEntry(
                "617-36\t16\tgi|555853|gb|U13369.1|HSU13369\t1232\t255\t18M\t"
                        + "*\t0\t0\tGAGCGCACCGGACGGGCC\tIIIIIIIIIIIIIIIIII\t"
                        + "XA:i:2\t MD:Z:1C9C6\t NM:i:3");
         Assert.assertEquals("2", mappedRead.getTagValue("XA"));
         Assert.assertEquals("1C9C6", mappedRead.getTagValue("MD"));
         Assert.assertEquals("3", mappedRead.getTagValue("NM"));
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
