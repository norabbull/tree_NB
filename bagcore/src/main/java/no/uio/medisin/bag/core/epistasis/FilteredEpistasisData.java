/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import no.uio.medisin.bag.core.aminoacid.CodonUsage;
import no.uio.medisin.bag.core.alignment.AlignPropCount;
import java.io.*;
import no.uio.medisin.bag.core.alignment.NumericalAlignment;
/**
 *
 * @author sr
 */
public class FilteredEpistasisData implements Serializable{

  private int maxCount;
  private int minCount;

  private int counts[][];
  private int sites[];
  private java.util.ArrayList<Integer> uniqueSites;
  private CodonUsage cu;


  public FilteredEpistasisData()
  {
    uniqueSites = new java.util.ArrayList<Integer>();
    cu = new CodonUsage();
  }

  public FilteredEpistasisData(int n)
  {
    counts = new int [n][n];
    //consensusSet = new int [n][n];
    uniqueSites = new java.util.ArrayList<Integer>();
    cu = new CodonUsage();
  }

  /**
   * find the range of values in the count array
   *
   * @param acp
   */
  public void findRange(AlignPropCount acp)
  {
    maxCount = 0;
    int ref[][] = acp.getCount();
    for(int i=0; i< ref.length;i++)
      for(int j=0; j<ref[i].length; j++)
      {
        if(ref[i][j] > maxCount)
          maxCount = ref[i][j];
      }

    minCount = maxCount;
   for(int i=0; i< ref.length;i++)
    for(int j=0; j<ref[i].length; j++)
    {
      if(ref[i][j] < minCount)
        minCount = ref[i][j];
    }

  }


  /**
   * Count up how many elements are above the cutoff
   *
   * @param acp
   * @param cutoff
   * @return
   */
  public int countByCeiling(int ref[][], int cutoff)
  {

   int fCount = 0;

   for(int i=0; i< ref.length;i++)
    for(int j=0; j<ref[i].length; j++)
      if(counts[i][j] > 0 && ref[i][j] >= cutoff)
        fCount++;
    return fCount;
  }


  /**
   * remove any elements that are below the ceiling
   *
   *
   * @param ref : the array to be filtered
   * @param ceiling : cutoff value - all elements below this value are set to 0
   */
  public int[][] filterByCeiling(int ref[][] , int cutoff)
  {
    //int [][] filteredCounts = new int[ref.length][ref.length];

    for(int i=0; i< ref.length;i++)
      for(int j=0; j<ref[i].length; j++)
        if(counts[i][j] > 0 && ref[i][j] >= cutoff)
          getCounts()[i][j] = ref[i][j];
        else
          getCounts()[i][j] = 0;


    return getCounts();
  }


  /**
   * find how many elements have sites with a conserved sequence greater than
   * the cutoff value
   *
   * @param acp
   * @param cutoff
   * @return
   */
  public int countByProbability(NumericalAlignment ntAl, double cutoff)
  {
    int fCount = 0;
    for(int i=0; i< ntAl.getEpiDataNum().length;i++)
      for(int j=0; j<ntAl.getEpiDataNum()[i].length; j++)
        if(counts[i][j] > 0)
          if(((ntAl.getConsensusProb()[i] > cutoff)
              || (ntAl.getConsensusProb()[j] > cutoff)))
            fCount++;

    return fCount;
  }


  /**
   * remove any elements that are likely to occur by chance, based on the
   * uncertainty in the consensus sequence
   * 
   * @param ref
   * @param counts
   * @param prob EpistasisFilterType eType = EpistasisFilterType.PAIR;
   * @param nSeqs
   * @return filteredCounts[][]  integer array of filtered results
   */
  public int[][] filterByProbability(NumericalAlignment ntAl, double cutoff)
  {

    for(int i=0; i< ntAl.getEpiDataNum().length;i++)
      for(int j=0; j<ntAl.getEpiDataNum()[i].length; j++)
        if(counts[i][j] > 0)
          if(((ntAl.getConsensusProb()[i] < cutoff)
              || (ntAl.getConsensusProb()[j] < cutoff)))
            getCounts()[i][j] = 0;
    return getCounts();
  }


  /**
   * find a list of all the sites in the array
   *
   * @param acb
   * @param ceiling
   */
  public int getUniqueSiteList(AlignPropCount acb, int ceiling)
  {
    int fCount = 0;
    for(int i=0; i< acb.getCount().length;i++)
    {
      for(int j=0; j<acb.getCount()[i].length; j++)
      {
        if(acb.getCount()[i][j] >= ceiling
                && getUniqueSites().contains(i)==false)
        {
          getUniqueSites().add(i);
          fCount++;
        }
      }
    }
    return fCount;

  }

  /**
   * this is a botched way to set the sites[] array in the object
   * because of a problem in how the AlignCoProb object is implemented
   *
   * @param na
   * @return
   */
  public int setSites(NumericalAlignment na)
  {
    sites = na.getSites();
    return getSites().length;
  }


  /**
   * write out the counts
   *
   * @param bw
   * @throws java.io.IOException
   */
  public void writeCounts(BufferedWriter bw) throws java.io.IOException
  {
    for(int s=0; s<getSites().length; s++)
      bw.write(getSites()[s] + "\t");
    bw.write("\n");

    for(int i=0; i< getCounts().length;i++)
    {
      for(int j=0; j<getCounts()[i].length; j++)
        bw.write(getCounts()[i][j] + "\t");
      bw.write("\n");
    }

  }

  /**
   * write out the list of unique sites in the array.
   *
   * @param bw
   * @throws java.io.IOException
   */
  public void writeUniqueSites(BufferedWriter bw) throws java.io.IOException
  {
    java.util.ListIterator itUS = uniqueSites.listIterator();
    while(itUS.hasNext())
    {
      bw.write((Integer)itUS.next() + "\n");
    }
  }

  /**
   * @return the maxCount
   */
  public int getMaxCount() {
    return maxCount;
  }

  /**
   * @param maxCount the maxCount to set
   */
  public void setMaxCount(int maxCount) {
    this.maxCount = maxCount;
  }

  /**
   * @return the minCount
   */
  public int getMinCount() {
    return minCount;
  }

  /**
   * @param minCount the minCount to set
   */
  public void setMinCount(int minCount) {
    this.minCount = minCount;
  }

  /**
   * @return the uniqueSites
   */
  public java.util.ArrayList<Integer> getUniqueSites() {
    return uniqueSites;
  }


  /**
   * @return the sites
   */
  public int[] getSites() {
    return sites;
  }

  /**
   * @return the filteredCounts
   */
  public int[][] getFilteredCounts() {
    return getCounts();
  }

  /**
   * @param filteredCounts the filteredCounts to set
   */
  public void setFilteredCounts(int[][] filteredCounts) {
    this.setCounts(filteredCounts);
  }

  /**
   * @return the counts
   */
  public int[][] getCounts() {
    return counts;
  }

  /**
   * @param counts the counts to set
   */
  public void setCounts(int[][] counts) {
    this.counts = counts;
  }

}
