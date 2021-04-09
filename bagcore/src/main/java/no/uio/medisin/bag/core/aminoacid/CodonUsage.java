/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

import java.io.*;

/**
 * Simple hash table to get information about codon usage for each amino acid
 *
 * @author sr
 */
public class CodonUsage implements Serializable {
  private static java.util.ArrayList<AAProperties> codon2AA;

  public CodonUsage()
  {
    codon2AA = new java.util.ArrayList<AAProperties>();
    initTable();
  }

  private void initTable()
  {

    getCodon2AA().add(new AAProperties("TTT", "Phe",  "F", "Phenylalanine"));
    getCodon2AA().add(new AAProperties("TTC", "Phe",  "F", "Phenylalanine"));
    getCodon2AA().add(new AAProperties("TTA", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("TTG", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("CTT", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("CTC", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("CTA", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("CTG", "Leu",  "L", "Leucine"));
    getCodon2AA().add(new AAProperties("ATT", "Ile",  "I", "Isoleucine"));
    getCodon2AA().add(new AAProperties("ATC", "Ile",  "I", "Isoleucine"));
    getCodon2AA().add(new AAProperties("ATA", "Ile",  "I", "Isoleucine"));
    getCodon2AA().add(new AAProperties("ATG[A]", "Met",  "M", "Methionine"));
    getCodon2AA().add(new AAProperties("GTT", "Val",  "V", "Valine"));
    getCodon2AA().add(new AAProperties("GTC", "Val",  "V", "Valine"));
    getCodon2AA().add(new AAProperties("GTA", "Val",  "V", "Valine"));
    getCodon2AA().add(new AAProperties("GTG", "Val",  "V", "Valine"));
    getCodon2AA().add(new AAProperties("TCT", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("TCC", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("TCA", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("TCG", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("CCT", "Pro",  "P", "Proline"));
    getCodon2AA().add(new AAProperties("CCC", "Pro",  "P", "Proline"));
    getCodon2AA().add(new AAProperties("CCA", "Pro",  "P", "Proline"));
    getCodon2AA().add(new AAProperties("CCG", "Pro",  "P", "Proline"));
    getCodon2AA().add(new AAProperties("ACT", "Thr",  "T", "Threonine"));
    getCodon2AA().add(new AAProperties("ACC", "Thr",  "T", "Threonine"));
    getCodon2AA().add(new AAProperties("ACA", "Thr",  "T", "Threonine"));
    getCodon2AA().add(new AAProperties("ACG", "Thr",  "T", "Threonine"));
    getCodon2AA().add(new AAProperties("GCT", "Ala",  "A", "Alanine"));
    getCodon2AA().add(new AAProperties("GCC", "Ala",  "A", "Alanine"));
    getCodon2AA().add(new AAProperties("GCA", "Ala",  "A", "Alanine"));
    getCodon2AA().add(new AAProperties("GCG", "Ala",  "A", "Alanine"));
    getCodon2AA().add(new AAProperties("TAT", "Tyr",  "Y", "Tyrosine"));
    getCodon2AA().add(new AAProperties("TAC", "Tyr",  "Y", "Tyrosine"));
    getCodon2AA().add(new AAProperties("TAA", "Och",  ".", "Ochre Stop"));
    getCodon2AA().add(new AAProperties("TAG", "Amb",  ".", "Amber Stop"));
    getCodon2AA().add(new AAProperties("CAT", "His",  "H", "Histidine"));
    getCodon2AA().add(new AAProperties("CAC", "His",  "H", "Histidine"));
    getCodon2AA().add(new AAProperties("CAA", "Gln",  "Q", "Glutamine"));
    getCodon2AA().add(new AAProperties("CAG", "Gln",  "Q", "Glutamine"));
    getCodon2AA().add(new AAProperties("AAT", "Asn",  "N", "Asparagine"));
    getCodon2AA().add(new AAProperties("AAC", "Asn",  "N", "Asparagine"));
    getCodon2AA().add(new AAProperties("AAA", "Lys",  "K", "Lysine"));
    getCodon2AA().add(new AAProperties("AAG", "Lys",  "K", "Lysine"));
    getCodon2AA().add(new AAProperties("GAT", "Asp",  "D", "Aspartic acid"));
    getCodon2AA().add(new AAProperties("GAC", "Asp",  "D", "Aspartic acid"));
    getCodon2AA().add(new AAProperties("GAA", "Glu",  "E", "Glutamic acid"));
    getCodon2AA().add(new AAProperties("GAG", "Glu",  "E", "Glutamic acid"));
    getCodon2AA().add(new AAProperties("TGT", "Cys",  "C", "Cysteine"));
    getCodon2AA().add(new AAProperties("TGC", "Cys",  "C", "Cysteine"));
    getCodon2AA().add(new AAProperties("TGA", "Opl",  ".", "Opal Stop"));
    getCodon2AA().add(new AAProperties("TGG", "Trp",  "W", "Tryptophan"));
    getCodon2AA().add(new AAProperties("CGT", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("CGC", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("CGA", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("CGG", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("AGT", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("AGC", "Ser",  "S", "Serine"));
    getCodon2AA().add(new AAProperties("AGA", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("AGG", "Arg",  "R", "Arginine"));
    getCodon2AA().add(new AAProperties("GGT", "Gly",  "G", "Glycine"));
    getCodon2AA().add(new AAProperties("GGC", "Gly",  "G", "Glycine"));
    getCodon2AA().add(new AAProperties("GGA", "Gly",  "G", "Glycine"));
    getCodon2AA().add(new AAProperties("GGG", "Gly",  "G", "Glycine"));


  }

  public String getShort(String codon)
  {
    java.util.ListIterator itAA = getCodon2AA().listIterator();
    while(itAA.hasNext())
    {
      AAProperties currentAAP = (AAProperties) itAA.next();
      if(currentAAP.getCodon().equals(codon.trim().toUpperCase()))
        return currentAAP.getShortName();
    }
    return "noMatch";

  }

  /**
   * @return the codon2AA
   */
  public java.util.ArrayList<AAProperties> getCodon2AA() {
    return codon2AA;
  }
}
