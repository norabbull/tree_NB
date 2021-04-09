/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mapping;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
public class MappedReadNGTest {

    public MappedReadNGTest() {
    }

    /**
     * this tests the compare method for two mapped reads
     * A read is be compared by location, strand and reference feature 
     * 
     * (i.e. chromosome)
     */
    @Test
    public void testMappedReadsEquivalence(){
        //(int s, int e, String c, String t, int n
        MappedRead m1 = new MappedRead(1, 2, "1", "+", 0);
        MappedRead m2 = new MappedRead(1, 2, "1", "+", 0);
        Assert.assertEquals(Boolean.TRUE,  (Boolean) m1.equals(m2));
        
        /* different read counts shouldnt affect equivalence */
        m2 = new MappedRead(1, 2, "1", "+", 1);
        Assert.assertEquals(Boolean.TRUE,  (Boolean) m1.equals(m2));        

        /* all of these should be unequal */
        MappedRead m3 = new MappedRead(2, 3, "1", "+", 0);
        Assert.assertEquals(Boolean.FALSE,  (Boolean) m1.equals(m3));
        
        m3 = new MappedRead(2, 3, "1", "+", 0);
        Assert.assertEquals(Boolean.FALSE,  (Boolean) m1.equals(m3));        
        
        m3 = new MappedRead(1, 2, "1", "-", 0);
        Assert.assertEquals(Boolean.FALSE,  (Boolean) m1.equals(m3));
        
        m3 = new MappedRead(1, 2, "2", "+", 0);
        Assert.assertEquals(Boolean.FALSE,  (Boolean) m1.equals(m3));
        
        
    }
    
    
    @Test 
    public void testCompare(){
        
        MappedRead m1 = new MappedRead(1, 2, "1", "+", 0);
        MappedRead m2 = new MappedRead(1, 3, "1", "+", 0);
        MappedRead m3 = new MappedRead(1, 2, "2", "+", 0);
        MappedRead m4 = new MappedRead(1, 2, "2", "-", 0);
        
        ArrayList<MappedRead> readList = new ArrayList<MappedRead>();
        readList.add(m4);
        readList.add(m3);
        readList.add(m2);
        readList.add(m1);
        Collections.sort(readList);
        
        Assert.assertEquals(readList.get(0), m3);
        Assert.assertEquals(readList.get(1), m2);
        Assert.assertEquals(readList.get(2), m1);
        Assert.assertEquals(readList.get(3), m4);
        
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
