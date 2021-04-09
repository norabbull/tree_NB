/*
 * aminoAcidList.java
 *
 * Created on January 4, 2009, 12:52 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
import java.io.BufferedWriter;
/**
 * Stores a list of amino acids and the number of times they occurred.
 * Useful for calculating entropy and mutual information
 * The list comprises a set of <code>Biojava</code>
 * <a href="http://www.biojava.org/docs/api/org/biojava/bio/symbol/SimpleAtomicSymbol.html">
 * SimpleSymbols</a> stored in a
 * <a href="http://www.biojava.org/docs/api/org/biojava/bio/symbol/SimpleSymbolList.html"
 * SimpleSymbolList</a>
 * since these store all kinds of useful information about the properties
 * of the amino acid.
 * The entropy can also be calculated and stored.
 * stored.
 *
 * @author sr
 */
public class aminoAcidListx implements java.io.Serializable {
   
  //Namespace wivaalns;

  private SimpleSymbolList aaList;
  private String listName;
  private java.util.Vector <Integer> symbolCount;
  private double entropy;

  /** Creates a new instance of aminoAcidList */
  public aminoAcidListx(String name, Symbol startAA)
  {
    listName = name;
    //wivaalns = (Namespace) org.biojavax.RichObjectFactory
    //        .getObject(SimpleNamespace.class, new Object[]{"wivaalNS"});
    aaList = new SimpleSymbolList(ProteinTools.getAlphabet());
    symbolCount = new java.util.Vector<Integer>();
    entropy = -1.0;
    try{
      aaList.addSymbol(startAA);
    }
    catch(BioException bex)
    {
      System.err.print("error adding first symbol to Symbol List\n");
      System.err.print(bex);
    }
    symbolCount.add(0);
//    symbolCount.set(0, 1);
  }
  public aminoAcidListx(String name)
  {
    listName = name;
//    wivaalns = (Namespace) org.biojavax.RichObjectFactory
//            .getObject(SimpleNamespace.class, new Object[]{"wivaalNS"});
    aaList = new SimpleSymbolList(ProteinTools.getAlphabet());
    symbolCount = new java.util.Vector<Integer>();
  }

  /**
   * Returns the amino acid at position <code>pos</code> within the list
   * @param pos 
   * @return the String-ified AminoAcid
   */
  public String getStringSymbolAt(int pos)
  {
    return aaList.subStr(pos, pos);
  }

  /**
   * increments the amino acid count of the specifed amino acid
   * if there is no match, then nothing happens
   * @param String aa - the single letter code for the amino acid
   * @return void
   */
  public void inc(String aa)
  {
    for(int a=0; a<aaList.length(); a++)
    {
      if(aaList.symbolAt(a+1).toString().equalsIgnoreCase(aa))
      {
        symbolCount.set(a, symbolCount.get(a)+1);
      }
    }
  }

  /**
   * increments the amino acid count of the specifed amino acid
   * if there is no match, then nothing happens
   * @param Symbol aa - the amino acid symbol
   * @return void
   */
  public void inc(Symbol aa)
  {
    for(int a=0; a<aaList.length(); a++)
    {
      if(aaList.symbolAt(a+1) == aa)
      {
        symbolCount.set(a, symbolCount.get(a)+1);
      }
    }
  }

  /**
   * checks whether the amino acid is present in the list
   * @param Symbol aa - the amino acid symbol
   * @return true if found, false if not
   */
  public boolean hasAA(String aa)
  {
    for(int a=0; a<aaList.length(); a++)
    {
       if(aaList.symbolAt(a+1).toString().equalsIgnoreCase(aa))
       {
          return true;
       }
    }
    return false;
  }

  /**
   * checks whether the amino acid is present in the list
   * @param String aa - the single letter code for the amino acid
   * @return true if found, false if not
   */
  public boolean hasAA(Symbol aa)
  {
    for(int a=0; a<aaList.length(); a++)
    {
       if(aaList.symbolAt(a+1) == aa)
       {
          return true;
       }
    }
    return false;
  }

  /**
   * writes a list summary to the specified BufferedWriter
   * @param BufferedWriter for output
   * @return void
   */
  public void reportList(BufferedWriter bw)
  {
    try{
      for(int a=0; a<aaList.length(); a++)
      {
          Symbol s = aaList.symbolAt(a+1);
          bw.write(s.getName() + "\t|"
                  + symbolCount.get(a) + "\n");
  /*        bw.write(aaStats.getLongName() + "\t| "
                     + aaStats.getLetter() + "|\t"
                     + aaStats.getPropertyString() + "|\t"
                     + aaStats.getCount() + "\n");*/
      }
    }
    catch(java.io.IOException iex)
    {
      System.err.print("error writing Site report to file\n");
      System.err.print(iex);
    }
  }

  /**
   * calculates the entropy of the list
   * @param none
   * @return entropy as a double
   */
  public double calcEntropy()
  {
    int noOfSeqs = 0;
    for(int a=0; a<symbolCount.size(); a++)
    {
      noOfSeqs += symbolCount.elementAt(a);
    }
    entropy = 0;
    for(int a=0;  a<symbolCount.size(); a++)
    {
      double prob = (double)symbolCount.elementAt(a)/(double)noOfSeqs;
      double logProb = Math.log(prob)/Math.log(2);
      double i = -1.0d*prob*logProb;
      entropy+= i;
    }
    return entropy;
  }

  /**
   * Adds an amino acid to the list.
   * @param aa. The Symbol of the new amino acid.
   */
  public void addAA(Symbol aa)
  {
    try{
      aaList.addSymbol(aa);
      symbolCount.add(1);
    }
    catch(BioException bex)
    {
      System.err.print("Error adding Symbol to Site\n");
      System.err.print(bex);
    }
  }
  
  // --------------------- accessor methods

  /**
   * returns the count for the specified amino acid symbol
   * @param Symbol for the amino acid query
   * @return count for the amino acid as an integer
   */
  public int getCount(Symbol qSym)
  {
    for(int a=0; a<aaList.length(); a++)
    {
       if(aaList.symbolAt(a+1) == qSym)
       {
          return 0;
       }
    }
    return -1;
  }

  /**
   * returns the amino acid list
   * @return amino acid list as a <a href="http://www.biojava.org/docs/api/org/biojava/bio/symbol/SimpleSymbolList.html"
   * SimpleSymbolList</a>
   */
  public SimpleSymbolList getAAList()
  {
    return aaList;
  }

  /**
   * returns the name of the AminoAcidList
   * @return name of the amino acid list as a String
   */
  public String getName()
  {
    return listName;
  }

  // --------------------- mutator methods

  /**
   * Sets the name of the AminoAcidList
   * @param newName AminoAcidList as a String
   */
  public void setName(String newName)
  {
    listName = newName;
  }

}
