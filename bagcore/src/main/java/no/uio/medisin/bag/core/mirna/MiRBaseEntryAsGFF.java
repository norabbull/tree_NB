/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

/**
 * convert a pre-miRNA/miRNA entry to miRBase GFF format
 * 
 * @author simonray
 */
public class MiRBaseEntryAsGFF {
  
    public static final int                 COL_CHR                 = 0;
    public static final int                 COL_SOURCE              = 1;
    public static final int                 COL_FEATURETYPE         = 2;
    public static final int                 COL_STARTPOS            = 3;
    public static final int                 COL_ENDPOS              = 4;
    public static final int                 COL_SCORE               = 5;
    public static final int                 COL_STRAND              = 6;
    public static final int                 COL_FRAME               = 7;
    public static final int                 COL_ATTRIBUTES          = 8;
  
    
    /**
     * convert pre-miR to GFF String
     * 
     * @param preMiR
     * @return 
     */
    public static String  preMiR2GFFEntry(PreMiRNA preMiR){
      
      String gffString = preMiR.getChromosome() + "\t.\tmiRNA_primary_transcript\t";
      gffString = gffString.concat(Integer.toString(preMiR.getStartPos())+"\t");
      gffString = gffString.concat(Integer.toString(preMiR.getEndPos())+"\t");
      gffString = gffString.concat(".\t");
      gffString = gffString.concat(preMiR.getStrand() + "\t.\t");
      gffString = gffString.concat("ID=" + preMiR.getMiID() + ";");
      gffString = gffString.concat("Alias=" + preMiR.getMiID() + ";");
      gffString = gffString.concat("Name=" + preMiR.getName() + "");
      return gffString;
    }
}
