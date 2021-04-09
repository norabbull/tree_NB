/*
 * aminoAcids.java
 *
 * Created on January 4, 2009, 10:16 AM
 *
 * contains properties of all amino acids
 */

package no.uio.medisin.bag.core.aminoacid;

import java.util.*;
/**
 *
 * @author sr
 */
public class aminoAcids {
   
   /** Creates a new instance of aminoAcids */
   private static Set<aminoAcidProperties> aaProperties 
           = new HashSet<aminoAcidProperties>();
   
   public aminoAcids() {
      aaProperties.add(new aminoAcidProperties("ISOLEUCINE",    "ILE", "I", "X....X...."));
      aaProperties.add(new aminoAcidProperties("LEUCINE",       "LEU", "L", "X....X...."));
      aaProperties.add(new aminoAcidProperties("VALINE",        "VAL", "V", "X.X..X...."));
      aaProperties.add(new aminoAcidProperties("CYSTEINE",      "CYS", "C", "X.X......."));
      aaProperties.add(new aminoAcidProperties("ALANINE",       "ALA", "A", "X.X.X....."));
      aaProperties.add(new aminoAcidProperties("GLYCINE",       "GLY", "G", "X.X.X....."));
      aaProperties.add(new aminoAcidProperties("METHIONINE",    "MET", "M", "X........."));
      aaProperties.add(new aminoAcidProperties("PHENYLALANINE", "PHY", "F", "X.....X..."));
      aaProperties.add(new aminoAcidProperties("TYROSINE",      "TYR", "Y", "XX....X..."));
      aaProperties.add(new aminoAcidProperties("TRYPTOPHAN",    "TRP", "W", "XX....X..."));
      aaProperties.add(new aminoAcidProperties("HISTIDINE",     "HIS", "H", "XX....XX.X"));
      aaProperties.add(new aminoAcidProperties("LYSINE",        "LYS", "K", ".X.....X.X"));
      aaProperties.add(new aminoAcidProperties("ARGININE",      "ARG", "R", ".X.....X.X"));
      aaProperties.add(new aminoAcidProperties("GLUTAMATE",     "GLU", "E", ".X......XX"));
      aaProperties.add(new aminoAcidProperties("GLUTAMINE",     "GLN", "Q", ".X........"));
      aaProperties.add(new aminoAcidProperties("ASPARTATE",     "ASP", "D", ".XX.....XX"));
      aaProperties.add(new aminoAcidProperties("ASPARAGINE",    "ASN", "N", ".XX......."));
      aaProperties.add(new aminoAcidProperties("SERINE",        "SER", "S", ".XX.X....."));
      aaProperties.add(new aminoAcidProperties("THREONINE",     "THR", "T", "X.X......."));
      aaProperties.add(new aminoAcidProperties("PROLINE",       "PRO", "P", ".XX.X....."));
      aaProperties.add(new aminoAcidProperties("ASPARTICACID",  "ASX", "B", ".X........"));
      aaProperties.add(new aminoAcidProperties("GLUTAMICACID",  "GLX", "Z", ".X........"));
      aaProperties.add(new aminoAcidProperties("UNKNOWN",       "XAA", "X", "XXXXXXXXXX"));
      aaProperties.add(new aminoAcidProperties("GAP",           "DAS", "-", "XXXXXXXXXX"));
      aaProperties.add(new aminoAcidProperties("ASTERIK",       "ASK", "*", "XXXXXXXXXX"));
      aaProperties.add(new aminoAcidProperties("QUESTION",      "DUH", "?", "XXXXXXXXXX"));
      aaProperties.add(new aminoAcidProperties("SQUIGGLY",      "SQL", "~", "XXXXXXXXXX"));
      
   }
 
   public aminoAcidProperties getProperties(String property)
   {
      aminoAcidProperties currAA = new aminoAcidProperties();
      switch(property.length())
      {
         // single letter
         case 1:
            for (Iterator<aminoAcidProperties> it = aaProperties.iterator(); it.hasNext();)
            {
               currAA = it.next();
               if(currAA.getLetter().equalsIgnoreCase(property))
               {
                  return currAA;
               }
            }
            currAA.setPropertyString("X");
            break;   
         
         case 3:
            for (Iterator<aminoAcidProperties> it = aaProperties.iterator(); it.hasNext();)
            {
               currAA = it.next();
               if(currAA.getShortName().equalsIgnoreCase(property))
               {
                  return currAA;
               }
            }
            currAA.setPropertyString("X");
            break;

         default:
            for (Iterator<aminoAcidProperties> it = aaProperties.iterator(); it.hasNext();)
            {
               currAA = it.next();
               if(currAA.getLongName().equalsIgnoreCase(property))
               {
                  return currAA;
               }
            }
            currAA.setPropertyString("X");
            break;
      }
      return currAA;
      
   }
   
   
}
