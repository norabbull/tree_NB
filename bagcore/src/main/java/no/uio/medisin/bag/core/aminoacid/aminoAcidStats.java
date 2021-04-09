/*
 * aminoAcidStats.java
 *
 * Created on January 4, 2009, 3:37 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

/**
 *
 * @author sr
 */
  public class aminoAcidStats extends aminoAcidProperties
          implements java.io.Serializable
  {
     int count;
     
     public aminoAcidStats()
     {
         count = 0;
     }
     
     public aminoAcidStats(String lName, String mName, String sName, 
             String propString)
     {
         super(lName, mName, sName, propString);
         count = 0;
     }
     
     public void incCount()
     {
        count++;
     }

     public int getCount()
     {
        return count;
     }
     
     public void setCount(int c)
     {
        count = c;
     }
     
     public void copy(aminoAcidStats aas)
     {
        super.copy(aas);
        this.setCount(aas.getCount());
     }
     
  }
