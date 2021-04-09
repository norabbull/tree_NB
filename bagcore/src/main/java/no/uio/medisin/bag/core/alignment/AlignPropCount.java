/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.alignment;

import java.io.*;
import no.uio.medisin.bag.core.epistasis.FilteredEpistasisData;
/**
 *
 * @author sr
 */
public class AlignPropCount implements Serializable{

  private int noOfDatasets;
  private double prob[][];
  private int count[][];
  private int noOfSigElements;
  private int sites[];

  public AlignPropCount(int len, int nS, int s[])
  {
    noOfDatasets = 0;
    prob = new double[len][len];
    count = new int[len][len];
    noOfSigElements = 0;
    sites = s;
  }

  public AlignPropCount(AlignPropCount acp, FilteredEpistasisData fepi)
  {
    noOfDatasets = acp.noOfDatasets;
    prob = acp.prob;
    count = fepi.getCounts();
    noOfSigElements = acp.noOfSigElements;
    sites = acp.sites;
  }

  public AlignPropCount()
  {
    
  }

  /**
   * get probability element (i,j)
   * @param i
   * @param j
   * @return
   */
  public double getProbAt(int i, int j)
  {
    return prob[i][j];
  }

  /**
   * get count element (i, j)
   * @param i
   * @param j
   * @return
   */
  public double getCountAt(int i, int j)
  {
    return getCount()[i][j];
  }

  /**
   * increment element in count array
   * @param i
   * @param j
   */
  private void incCount(int i, int j)
  {
    getCount()[i][j]++;
  }

  /**
   * compare two arrays element by element.
   * One is a reference array, one is a query array.  If the ref array element value
   * is greater than the query element then increment the corresponding element
   * in the count array.  This is used to calculate probability for bootstrapping
   *
   * @param ref
   * @param query
   * @throws java.lang.Exception
   */
  public void compare(int ref[][], int query[][]) throws java.lang.Exception
  {
    if(ref.length != getSites().length
            || ref.length != query.length
            || ref[0].length != query[0].length)
      throw new java.lang.ArrayIndexOutOfBoundsException(
              "the reference and query arrays are different sizes");

    for(int i=0; i< ref.length;i++)
      for(int j=0; j<ref[i].length; j++)
        if(ref[i][j] > query[i][j])
          getCount()[i][j]++;
    noOfDatasets++;
  }

  /**
   * adds the counts from the supplied array to the counts stored in this instance
   * also updates the noOfDatasets
   * @param c
   * @param cNoData
   * @throws java.lang.Exception
   */
  public void addCounts(int c[][]) throws java.lang.Exception
  {
    if(getCount().length != c.length || getCount()[0].length != c[0].length)
      throw new java.lang.ArrayIndexOutOfBoundsException(
              "the specified and local arrays are different sizes");
    for(int i=0; i< getCount().length;i++)
      for(int j=0; j<getCount()[i].length; j++)
          count[i][j]+=c[i][j];
  }

  /**
   * find which elements are above the defined pValue
   *
   * @param pValue
   */
  public void calcProb(double pValue, int noOfData)
  {
    for(int i=0; i< getCount().length;i++)
      for(int j=0; j<getCount()[i].length; j++)
      {
        prob[i][j]=0.0;
        if(getCount()[i][j] >= noOfData*(1-pValue))
        {
          prob[i][j] = (double)getCount()[i][j]/(double)noOfData;
          noOfSigElements++;
        }
      }
  }

  /**
   * write the calculated probabilities
   *
   * @param bw
   * @throws java.io.IOException
   */
  public void writeProb(BufferedWriter bw) throws java.io.IOException
  {
    for(int s=0; s<getSites().length; s++)
      bw.write(getSites()[s] + "\t");
    bw.write("\n");

    for(int i=0; i< prob.length;i++)
    {
      for(int j=0; j<prob[i].length; j++)
        bw.write(prob[i][j] + "\t");
      bw.write("\n");
    }
  }

  /**
   * write the counts
   *
   * @param bw
   * @throws java.io.IOException
   */
  public void writeCount(BufferedWriter bw ) throws java.io.IOException
  {
    for(int s=0; s<getSites().length; s++)
      bw.write(getSites()[s] + "\t");
    bw.write("\n");

    for(int i=0; i< getCount().length;i++)
    {
      for(int j=0; j<getCount()[i].length; j++)
        bw.write(getCount()[i][j] + "\t");
      bw.write("\n");
    }

  }

  /**
   * write points as a list
   * 
   * @param bw
   * @throws java.io.IOException
   */
  public void writePointsList(BufferedWriter bw) throws java.io.IOException
  {
    for(int i=0; i<getCount().length; i++)
      for(int j=0; j<getCount().length; j++)
        if(getCount()[i][j] > 0)
          bw.write(i +"\t" + j + "\n");

  }
  /**
   * increment the number of datasets;
   */
  public void incNoOfDatasets()
  {
    noOfDatasets++;
  }
  /**
   * @return the noOfSeqs
   */
  public int getnoOfDatasets() {
    return noOfDatasets;
  }

  /**
   * @param noOfSeqs the noOfSeqs to set
   */
  public void setnoOfDatasets(int noOfData) {
    this.noOfDatasets = noOfData;
  }

  /**
   * @return the noOfSigElements
   */
  public int getNoOfSigElements() {
    return noOfSigElements;
  }

  /**
   * @return the count
   */
  public int[][] getCount() {
    return count;
  }

  /**
   * @return the sites
   */
  public int[] getSites() {
    return sites;
  }

  /**
   * @param sites the sites to set
   */
  public void setSites(int[] sites) {
    this.sites = sites;
  }

}
