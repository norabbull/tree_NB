/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;
import no.uio.medisin.bag.core.aminoacid.aminoAcidListx;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import java.io.*;

/**
 *
 * @author sr
 */
public class wivSimpleSymbolList implements java.io.Serializable{

  public wivSimpleSymbolList()
  {
    try
    {
    sdx = new aminoAcidListx("nameOfList");
    sdx2 = new aminoAcidListx("nameOfList");
    Sequence seq = ProteinTools.createProteinSequence("futnt", "seqname");
    Sequence seq2 = ProteinTools.createProteinSequence("slgg", "seqname2");
    sdx.addAA(seq.symbolAt(7));
    sdx.addAA(seq.symbolAt(6));
    sdx.addAA(seq.symbolAt(5));
    sdx.addAA(seq.symbolAt(4));
    sdx.addAA(seq.symbolAt(3));
    sdx.addAA(seq.symbolAt(2));
    sdx.addAA(seq.symbolAt(1));
    sdx2.addAA(seq2.symbolAt(1));
    sdx2.addAA(seq2.symbolAt(2));
    sdx2.addAA(seq2.symbolAt(3));
    sdx2.addAA(seq2.symbolAt(4));
    ssl = new SimpleSymbolList(seq);
    ssl2 = new SimpleSymbolList(seq2);
    }
    catch(org.biojava.bio.symbol.IllegalSymbolException exBJ)
    {
      System.err.print("couldn't create sequence");
      exBJ.printStackTrace();
    }
  }
  public void writeOut()
  {
    try{
      ObjectOutputStream osS =
            new ObjectOutputStream(new FileOutputStream("data.txt"));
      osS.writeObject(ssl);
      osS.writeObject(ssl2);
      osS.writeObject(sdx);
      osS.writeObject(sdx2);
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing data");
      exIO.printStackTrace();
    }
    
  }
  
  public void ReadIn()
  {
    try{
      ObjectInputStream osI = 
            new ObjectInputStream(new FileInputStream("data.txt"));
      ssl2 = (SimpleSymbolList) osI.readObject();
      ssl = (SimpleSymbolList) osI.readObject();
      sdx2 = (aminoAcidListx) osI.readObject();
      sdx = (aminoAcidListx) osI.readObject();

    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing data");
      exIO.printStackTrace();
    }
    catch(java.lang.ClassNotFoundException exLang)
    {
      System.err.println("problem casting object");
      exLang.printStackTrace();
    }
  }



  public SimpleSymbolList ssl;
  public SimpleSymbolList ssl2;
  aminoAcidListx sdx;
  aminoAcidListx sdx2;
}
