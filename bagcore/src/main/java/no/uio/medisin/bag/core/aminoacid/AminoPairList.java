/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;


import no.uio.medisin.bag.core.epistasis.EpistaticAAPair;


//import no.uio.medisin.bag.core.epistasis.EpistaticAAPair;


/**
 *
 * @author sr
 */
public class AminoPairList {
  private java.util.ArrayList<EpistaticAAPair> aaPairs;

  public AminoPairList()
  {
    aaPairs = new java.util.ArrayList<EpistaticAAPair>();
  }

  /**
   * see whether the list contains the provided AApair
   *
   * @param aacid1
   * @param aacid2
   * @return
   */
  public int contains(String aacid1, String aacid2)
  {
    int i=0;
    java.util.ListIterator itAA = getAaPairs().listIterator();
    while(itAA.hasNext())
    {
      EpistaticAAPair aaPair = (EpistaticAAPair) itAA.next();
      if(aacid1.equalsIgnoreCase(aaPair.getAa1()) &&
              aacid2.equalsIgnoreCase(aaPair.getAa2()))
        return i;
      i++;
    }
    return -1;
  }

  /**
   * increment the count of the specified AApair
   *
   * @param aacid1
   * @param aacid2
   * @return
   */
  public int incCount(String aacid1, String aacid2)
  {
    int i=0;
    java.util.ListIterator itAA = getAaPairs().listIterator();
    while(itAA.hasNext())
    {
      EpistaticAAPair aaPair = (EpistaticAAPair) itAA.next();
      if(aacid1.equalsIgnoreCase(aaPair.getAa1()) &&
              aacid2.equalsIgnoreCase(aaPair.getAa2()))
      {
        aaPair.setCount(aaPair.getCount()+1);
        return i;
      }
      i++;
    }
    return -1;
  }

  /**
   * increment the count of the specified element
   * 
   * @param index
   * @return
   */
  public int incCount(int index)
  {
    if(index < getAaPairs().size())
    {
      EpistaticAAPair aaPair = (EpistaticAAPair)getAaPairs().get(index);
      aaPair.setCount(aaPair.getCount()+1);
        return aaPair.getCount();
    }
    return -1;
  }

  /**
   * get the element at index i
   *
   * @param index
   * @return
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public EpistaticAAPair getAApair(int index) throws java.lang.ArrayIndexOutOfBoundsException
  {
    return (EpistaticAAPair) aaPairs.get(index);
  }

  /**
   * set the element at index i
   *
   * @param index
   * @param aaPair
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void setAApair(int index, EpistaticAAPair aaPair)
          throws java.lang.ArrayIndexOutOfBoundsException
  {
    aaPairs.set(index, aaPair);
  }

  /**
   * add new element at end of list
   * 
   * @param aaPair
   */
  public void add(EpistaticAAPair aaPair)
  {
    aaPairs.add(aaPair);
  }

  /**
   * @return the aaPairs
   */
  public java.util.ArrayList<EpistaticAAPair> getAaPairs() {
    return aaPairs;
  }

  /**
   * @param aaPairs the aaPairs to set
   */
  public void setAaPairs(java.util.ArrayList<EpistaticAAPair> aaPairs) {
    this.aaPairs = aaPairs;
  }

}
