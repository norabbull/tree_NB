/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

import org.testng.Assert;
import static no.uio.medisin.bag.core.mirna.MiRNA.logger;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 *
 * @author simonray
 */
public class MiRNATestNGTest {
    
    public MiRNATestNGTest() {
    }

    
    @Test(description="test the class constructors")
    public static void testConstructors(){
        
        MiRNA mirseq1 = new MiRNA();
        Assert.assertEquals(0, mirseq1.getLength());
        Assert.assertEquals("", mirseq1.getName());
        Assert.assertNull(mirseq1.getAccessionNumber());
        
        
//      public MiRNA(String cName, String cChrom, int cStart, int cEnd, String cStrand){
        MiRNA mirseq2 = new MiRNA("miR1", "1", 1, 2, "+");
        Assert.assertEquals(0, mirseq2.getLength());
        Assert.assertEquals("miR1", mirseq2.getName());
        Assert.assertNull(mirseq2.getAccessionNumber());
        
        Assert.assertEquals("1", mirseq2.getChromosome());
        Assert.assertEquals("+", mirseq2.getStrand());
        Assert.assertEquals(1, mirseq2.getStartPos());
        Assert.assertEquals(2, mirseq2.getEndPos());
        
        
//      public MiRNA(String cName, String cChrom, int cStart, int cEnd, String cStrand, String cMIMAT, String cParentID, String cSeq)        
        MiRNA mirseq3 = new MiRNA("miR1", "1", 1, 2, "+", "MIMAT1", "PARENT", "ACGT");
        Assert.assertEquals(4, mirseq3.getLength());
        Assert.assertEquals("miR1", mirseq3.getName());
        Assert.assertEquals("ACGT", mirseq3.getSeq());
        Assert.assertEquals("PARENT", mirseq3.getParent());
        Assert.assertNull(mirseq3.getAccessionNumber());
        
        Assert.assertEquals("1", mirseq3.getChromosome());
        Assert.assertEquals("+", mirseq3.getStrand());
        Assert.assertEquals(1, mirseq3.getStartPos());
        Assert.assertEquals(2, mirseq3.getEndPos());
        
        mirseq3.setSeq("ACCGGGTTTT");
        Assert.assertEquals(10, mirseq3.getLength());
        mirseq3.setStartPos(21);
        Assert.assertEquals(21, mirseq3.getStartPos());
        Assert.assertEquals(30, mirseq3.getEndPos());

        mirseq3.setEndPos(40);
        mirseq3.setStartPos(mirseq3.getEndPos()-mirseq3.getLength()+1);
        Assert.assertEquals(31, mirseq3.getStartPos());
        
    }

    
    @Test(description = "test the functions for parsing miRBase entries")
    public static void testMiRBaseParsing()
    {
        String miRBaseLine = "chr1\t.\tmiRNA\t3580220\t3580240\t.\t+\t.\tID=MIMAT0011792;"
                + "Alias=MIMAT0011792;Name=bta-miR-2284i;Derives_from=MI0011294";
        MiRNA mirseq3 = new MiRNA();
        mirseq3.parseMiRBaseGFFEntry(miRBaseLine);
        Assert.assertEquals("1", mirseq3.getChromosome());
        Assert.assertEquals(3580220,mirseq3.getStartPos());
        Assert.assertEquals("+",mirseq3.getStrand());
        
        Assert.assertEquals("bta-miR-2284i", mirseq3.getName());
        Assert.assertEquals("MIMAT0011792", mirseq3.getMimatID());
        Assert.assertEquals("MIMAT0011792", mirseq3.getMimatAlias());
        Assert.assertEquals("MI0011294", mirseq3.getParent());
        
    }
    
    @Test(description="test feature push/pull")
    public static void testFeatureMethods(){
        MiRNA mirseq4 = new MiRNA("miR1", "1", 1, 2, "+", "MIMAT1", "PARENT", "ACGT");

        mirseq4.addFeature("miRNA_Feature_1",                  "miRNA_feature_value_1");
        mirseq4.addFeature("miRNA_Feature_2",                  "miRNA_feature_value_2");
        mirseq4.addFeature("miRNA_Feature_3",                  "miRNA_feature_value_3");
        
        Assert.assertEquals("miRNA_feature_value_1",mirseq4.getFeatureSet().get("miRNA_Feature_1"));
        
    }
    
    
    
    @Test(description="test isomiR processing function")
    public static void testIsomiRProcessing(){

        String miRBaseLine = "chr1\t.\tmiRNA\t3580220\t3580240\t.\t+\t.\tID=MIMAT0011792;"
                        + "Alias=MIMAT0011792;Name=bta-miR-2284i;Derives_from=MI0011294";
        MiRNA mirseq = new MiRNA();
        mirseq.parseMiRBaseGFFEntry(miRBaseLine);
        
        mirseq.addIsomiR("isomiR1-1000", 3580220, "18M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "GAGCGCACCGGACGGGCC");
        mirseq.addIsomiR("isomiR2-100",  3580220, "16M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "GTGCGCACCGGACGGGCC");
        mirseq.addIsomiR("isomiR3-60",   3580221, "18M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2",  "AGCGCACCGGACGGGCCT");
        mirseq.addIsomiR("isomiR4-25",   3580221, "19M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "AAGCGCACCGGACGGGCC");
        
        logger.error(mirseq.prettyReportIsomiRs(5, 200));
        
        // this should report 2 isomiRs for the cutoff of 5 percent + header and footer lines
        Assert.assertEquals(4, mirseq.prettyReportIsomiRs(6, 200).split(System.getProperty("line.separator")).length);
        
        // this should report 3 isomiRs for the cutoff of 4 percent + header and footer lines
        Assert.assertEquals(5, mirseq.prettyReportIsomiRs(5, 200).split(System.getProperty("line.separator")).length);
        
        // this should report 0 isomirs because the miRNA doesn't have sufficient counts
        Assert.assertEquals(1, mirseq.prettyReportIsomiRs(4, 1186).split(System.getProperty("line.separator")).length);
        
        
        mirseq.removeIsomiRs();

        // this should report 0 isomirs because we just removed them
        Assert.assertEquals(1, mirseq.prettyReportIsomiRs(4, 1).split(System.getProperty("line.separator")).length);
        
        
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
