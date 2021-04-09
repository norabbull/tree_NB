/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.alignment;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
//import org.biojava.bio.seq.*;
/**
 * A feeble object for storing amino acid sequences.
 * It has the sole advantage that it is Serializable.
 *
 * @author sr
 */
public class PoorAASequence implements java.io.Serializable{

  private String  sequenceName;
  private SimpleSymbolList sequence;
  private long position;

  /**
   * creates empty PoorAASequence
   */
  public PoorAASequence()
  {
    sequence = new SimpleSymbolList(ProteinTools.getAlphabet());
  }

  /**
   * creates PoorAASequence with Name and Sequence
   * @param name: Sequence Name
   * @param seq:  Amino Acid Sequence
   * 
   */
  public PoorAASequence(String name, String seq)
  {
    try{
       CharacterTokenization ct
               = new CharacterTokenization(
                      ProteinTools.getAlphabet(), false);
      sequence = new SimpleSymbolList(ProteinTools.getAlphabet());
      Sequence s = ProteinTools.createProteinSequence(seq, name);
      for(int i=1;i<=s.length(); i++)
      {
        sequence.addSymbol(s.symbolAt(i));
      }
      sequenceName = name;
    }
    catch(IllegalSymbolException exBJ)
    {
      System.err.println("couldn't create Symbol List to store sequence");
      exBJ.printStackTrace();
      System.err.println(name);
      System.err.println(seq);
    }
  }

  /************************************************************************/
  /*                            ACCESSOR METHODS                          */
  /************************************************************************/

  /**
   * get sequence name
   * @return sequence name as String
   *
   */
  public String getName()
  {
    return sequenceName;
  }

  /**
   * get Sequence Length
   *
   * @return sequence length as integer
   */
  public int getLength()
  {
    return sequence.length();
  }

  /**
   * return the sequence
   * 
   * @return
   */
  public String getSequence()
  {
    return sequence.seqString();
  }


  /**
   * return sub sequence.  first position in sequence is 1
   * If start is greater than stop assumes a circularized sequence
   * and returns a sequence tied together across the join.
   * 
   * @param start first position
   * @param end   last position
   * @return  Subsequence as String
   */
  public String subSeq(int start, int end)
  {
    try{
      if(start > sequence.length())
      {
        System.out.println("warning: start pos is beyond the end of alignment");
        System.out.println("         subtracting alignment length");
        System.out.println("         start pos is     " + start);
        while(start > sequence.length())
        {
          start -= sequence.length();
          System.out.println("         ---------------> " + start);
        }
        System.out.println("         start pos is now " + start);

      }
      if(start > end)
      {
        return sequence.subStr(start, sequence.length())
                + sequence.subStr(1, end);
      }
      else
      {
        return sequence.subStr(start, end);
      }
    }
    catch(java.lang.Exception exLang)
    {
      System.err.println(start + " --> " + end);
      System.err.println("error in sequence coords.  " +
              "Note: first position in sequence is 1");
      exLang.printStackTrace();
    }
    return "";
  }

  /**
   * Returns Amino acid at position b in the sequence
   * @param b base position
   * @return String representation of amino acid
   * 
   */
  public String baseAt(int b)
  {
    return sequence.subStr(b, b);
  }

  /**
   * Returns Amino acid at position b in the sequence
   * @param b base position
   * @return Symbol representation of amino acid
   *
   */
  public Symbol symbolAt(int b)
  {
    return sequence.symbolAt(b);
  }
}
