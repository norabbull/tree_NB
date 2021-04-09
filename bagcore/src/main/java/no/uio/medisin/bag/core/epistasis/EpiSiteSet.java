/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import no.uio.medisin.bag.core.aminoacid.CodonUsage;
import no.uio.medisin.bag.core.alignment.AlignPropCount;
import no.uio.medisin.bag.core.tree.BranchMutualPair;
import java.io.*;
import java.util.*;
import no.uio.medisin.bag.core.InformationLevel;
import static no.uio.medisin.bag.core.InformationLevel.LONG;
import static no.uio.medisin.bag.core.InformationLevel.MEDIUM;
import static no.uio.medisin.bag.core.InformationLevel.SHORT;
import no.uio.medisin.bag.core.alignment.NumericalAlignment;
import no.uio.medisin.bag.core.tree.SiteMutData;
import no.uio.medisin.bag.core.tree.SubstitutionType;
import no.uio.medisin.bag.core.tree.WivTree;
/**
 *
 * @author sr
 */
public class EpiSiteSet implements Serializable {
  private WivTree wTree;
  private no.uio.medisin.bag.core.alignment.AlignPropCount epiData;
  private NumericalAlignment ntAlFull;
  private CodonUsage cu;
  private java.util.ArrayList<EpiDataPoint> epiDataPointList;
  private java.util.ArrayList<BranchMutualPair> branchMutPairList;
  private EpiFilterSet range;



  public EpiSiteSet()
  {
    cu = new CodonUsage();
    epiDataPointList = new java.util.ArrayList<EpiDataPoint>();
    wTree = new WivTree();
    epiData = new AlignPropCount();
    ntAlFull = new NumericalAlignment();
    branchMutPairList = new java.util.ArrayList<BranchMutualPair>();
    range = new EpiFilterSet();
  }

  /**
   * read newark tree into wivtree instance
   *
   * @param filename
   * @throws java.io.IOException
   * @throws java.lang.ClassNotFoundException
   */
  public void readTree(String filename) throws java.io.IOException
  {
    //BufferedReader brWT = new BufferedReader(new FileReader(filename));
      wTree.readTree(filename);
    //brWT.close();
  }

  /**
   * read in epistatic data properties file. This file has a .epi.proc extension
   * and is created by running ProcessEpiData.java
   *
   * @param filename
   * @throws java.io.IOException
   */
  public void readEpiDataProps(String filename) throws java.io.IOException,
          java.lang.ClassNotFoundException
  {
    ObjectInputStream isEA =
          new ObjectInputStream(
            new FileInputStream(filename));
      epiData = (AlignPropCount)isEA.readObject();
    isEA.close();
  }

  /**
   * read a full alignment including conserved sites in WNA format.
   * A PHYLIP or FASTA alignment can be converted to WNA by running aln2wna.java
   *
   * @param filename
   * @throws java.io.IOException
   * @throws java.lang.ClassNotFoundException
   */
  public void readFullAlign(String filename) throws java.io.IOException,
          java.lang.ClassNotFoundException
  {
    ObjectInputStream isNA =
          new ObjectInputStream(
            new FileInputStream(filename));
      ntAlFull = (NumericalAlignment)isNA.readObject();
      ntAlFull.calcConsensus();
    isNA.close();
  }


  /**
   * calculate the properties of each data point
   *
   *
   * @param na
   */
  public void getDetailedEpiProperties()
  {
   for(int i=0; i< epiData.getCount().length;i++)
    for(int j=0; j<epiData.getCount()[i].length; j++)
    {
      if(epiData.getCount()[i][j] > 0) // get the data for this pair
      {
        EpiDataPoint epiDP = new EpiDataPoint(i, j, epiData.getCount().length);
        int iReal = epiData.getSites()[i];
        int jReal = epiData.getSites()[j];
        epiDP.calcEpiDataPtProperties(ntAlFull, wTree, iReal, jReal);
        epiDataPointList.add(epiDP);
      }
    }
    System.out.println("found " + epiDataPointList.size() + " epistatic data points");
  }



  /**
   * This is the same as the @ref getDetailedEpiProperties() except it
   * works through a list of site pairs rather than the entire count array
   *
   * first index is the sequence, second index is the site
   * i.e. numData[seq][site]
   *
   * @param na
   */
  public void getDetailedEpiProperties(java.util.ArrayList<java.awt.Point> list)
  {
    int p=0;
    java.util.ListIterator itPt = list.listIterator();
    while(itPt.hasNext())
    {
      java.awt.Point pt = (java.awt.Point) itPt.next();
      int i=pt.x;
      int j=pt.y;
      if(epiData.getCount()[i][j] > 0) // get the data for this pair
      {
        EpiDataPoint epiDP = new EpiDataPoint(i, j, epiData.getCount().length);
        int iReal = epiData.getSites()[i];
        int jReal = epiData.getSites()[j];
        epiDP.calcEpiDataPtProperties(ntAlFull, wTree, iReal, jReal);
        epiDataPointList.add(epiDP);
      }
      p++;
      if(p%10==0)
        System.out.println("-" + p);
    }
    System.out.println("found " + epiDataPointList.size() + " epistatic data points");
  }

  /**
   * get the range of values
   */
  public void findRange()
  {
    getRange().resetToZero();

    java.util.ListIterator itEDP = epiDataPointList.listIterator();

    EpiDataPoint edp = (EpiDataPoint) itEDP.next();
    getRange().setMinSMetric(edp.getSMetric());
    getRange().setMaxSMetric(edp.getSMetric());
    getRange().setMinEpiSeqCount(edp.getEpiSeqCount());
    getRange().setMaxEpiSeqCount(edp.getEpiSeqCount());
    getRange().setMinMutSeqCount(edp.getOneMutSeqCount());
    getRange().setMaxMutSeqCount(edp.getOneMutSeqCount());
    getRange().setMinAlignDist(edp.getAlignDist());
    getRange().setMaxAlignDist(edp.getAlignDist());
    getRange().setMinBLcount(edp.getBranchLengthList().size());
    getRange().setMaxBLcount(edp.getBranchLengthList().size());
    getRange().setMinBL(edp.getBlMin());
    getRange().setMaxBL(edp.getBlMax());
    getRange().setMinBLspread(edp.getBlMax()-edp.getBlMin());
    getRange().setMinBLspread(edp.getBlMax()-edp.getBlMin());
    getRange().setMinBLmean(edp.getBlMean());
    getRange().setMaxBLmean(edp.getBlMean());
    getRange().setMinBLsd(edp.getBlSD());
    getRange().setMaxBLsd(edp.getBlSD());
    getRange().setMinBLmoment(edp.getBlMoment());
    getRange().setMaxBLmoment(edp.getBlMoment());


    while(itEDP.hasNext())
    {
      edp = (EpiDataPoint) itEDP.next();
      if(getRange().getMinSMetric() > edp.getSMetric())
        getRange().setMinSMetric(edp.getSMetric());
      if(getRange().getMaxSMetric() < edp.getSMetric())
        getRange().setMaxSMetric(edp.getSMetric());

      if(getRange().getMinEpiSeqCount() > edp.getEpiSeqCount())
        getRange().setMinEpiSeqCount(edp.getEpiSeqCount());
      if(getRange().getMaxEpiSeqCount() < edp.getEpiSeqCount())
        getRange().setMaxEpiSeqCount(edp.getEpiSeqCount());

      if(getRange().getMinMutSeqCount() > edp.getOneMutSeqCount())
        getRange().setMinMutSeqCount(edp.getOneMutSeqCount());
      if(getRange().getMaxMutSeqCount() < edp.getOneMutSeqCount())
        getRange().setMaxMutSeqCount(edp.getOneMutSeqCount());

      if(getRange().getMinAlignDist() > edp.getAlignDist())
        getRange().setMinAlignDist(edp.getAlignDist());
      if(getRange().getMaxAlignDist() < edp.getAlignDist())
        getRange().setMaxAlignDist(edp.getAlignDist());

      if(getRange().getMinBLcount() > edp.getBranchLengthList().size())
        getRange().setMinBLcount(edp.getBranchLengthList().size());
      if(getRange().getMaxBLcount() < edp.getBranchLengthList().size())
        getRange().setMaxBLcount(edp.getBranchLengthList().size());

      if(getRange().getMinBL() > edp.getBlMin())
        getRange().setMinBL(edp.getBlMin());
      if(getRange().getMaxBL() < edp.getBlMax())
        getRange().setMaxBL(edp.getBlMax());

      if(getRange().getMinBLspread() > (edp.getBlMax()-edp.getBlMin()))
        getRange().setMinBLspread(edp.getBlMax()-edp.getBlMin());
      if(getRange().getMaxBLspread() < (edp.getBlMax()-edp.getBlMin()))
        getRange().setMaxBLspread(edp.getBlMax()-edp.getBlMin());

      if(getRange().getMinBLmean() > edp.getBlMean())
        getRange().setMinBLmean(edp.getBlMean());
      if(getRange().getMaxBLmean() < edp.getBlMean())
        getRange().setMaxBLmean(edp.getBlMean());

      if(getRange().getMinBLsd() > edp.getBlSD())
        getRange().setMinBLsd(edp.getBlSD());
      if(getRange().getMaxBLsd() < edp.getBlSD())
        getRange().setMaxBLsd(edp.getBlSD());

      if(getRange().getMinBLmoment() > edp.getBlMoment())
        getRange().setMinBLmoment(edp.getBlMoment());
      if(getRange().getMaxBLmoment() < edp.getBlMoment())
        getRange().setMaxBLmoment(edp.getBlMoment());
    }

  }

  public String printRange()
  {
    String s =
          "SMetric:\t minimum = " + getRange().getMinSMetric()
            + "\tmaximum = " + getRange().getMaxSMetric() + "\n"
        + "EpiSeqCount:\t minimum = " + getRange().getMinEpiSeqCount()
            + "\tmaximum = " + getRange().getMaxEpiSeqCount() + "\n"
        + "MutSeqCount:\t minimum = " + getRange().getMinEpiSeqCount()
            + "\tmaximum = " + getRange().getMaxEpiSeqCount() + "\n"
        + "BLcount:\t minimum = " + getRange().getMinBLcount()
            + "\tmaximum = " + getRange().getMaxBLcount() + "\n"
        + "BLspread:\t minimum = " + getRange().getMinBLspread()
            + "\tmaximum = " + getRange().getMaxBLspread() + "\n"
        + "BL:\t minimum = " + getRange().getMinBL()
            + "\tmaximum = " + getRange().getMaxBL() + "\n"
        + "BLmean:\t minimum = " + getRange().getMinBLmean()
            + "\tmaximum = " + getRange().getMaxBLmean() + "\n"
        + "BLSD:\t minimum = " + getRange().getMinBLsd()
            + "\tmaximum = " + getRange().getMaxBLsd() + "\n"
        + "BLmoment:\t minimum = " + getRange().getMinBLmoment()
            + "\tmaximum = " + getRange().getMaxBLmoment() + "\n";
    return s;
  }

  /**
   *  for each epistatic data point
   *   (i)  calculate the mutual information
   *   (ii) generate a double [] array of the branch lengths and calculate
   *        the standard deviation
   *  Add these to a list
   *
   */
  public void getBranchLengthMIpairs()
  {
    java.util.ListIterator itPL = epiDataPointList.listIterator();
    while(itPL.hasNext())
    {
      EpiDataPoint epiDataPoint = (EpiDataPoint) itPL.next();
      branchMutPairList.add(
              new BranchMutualPair(
              epiDataPoint.getSd(),
              epiDataPoint.getMutualIJ()));
    }
  }


  /**
   * find how many points have Mutual Information values within the specified
   * window
   *
   * @param low
   * @param high
   * @return
   */
  public int countMIbyWindow(double low, double high)
  {
    int count = 0;
    java.util.ListIterator itPL = epiDataPointList.listIterator();
    while(itPL.hasNext())
    {
      EpiDataPoint epiDataPoint = (EpiDataPoint) itPL.next();
      if(epiDataPoint.getMutualIJ() >= low && epiDataPoint.getMutualIJ() <= high)
        count++;
    }
    return count;
  }


  /**
   * count how many points are within the specified window
   *
   * @param epiFilter
   * @return
   */
  public int countbyWindow(no.uio.medisin.bag.core.epistasis.EpiFilterSet epiFilter)
  {
    int count = 0;
    java.util.ListIterator itPL = epiDataPointList.listIterator();
    while(itPL.hasNext())
    {
      EpiDataPoint epiDataPoint = (EpiDataPoint) itPL.next();
      if(epiDataPoint.isPointWithinWindow(epiFilter))
        count++;
    }
    return count;

  }

  /**
   * 
   * write points that are within the specified window
   *
   * @param epiFilter
   * @param level
   * @param bw
   * @return
   * @throws java.io.IOException
   */
  public int writebyWindow(
          no.uio.medisin.bag.core.epistasis.EpiFilterSet epiFilter,
          InformationLevel level,
          BufferedWriter bw) throws java.io.IOException
  {
    int count = 0;
    bw.write(EpiDataPoint.printShortHeader());
    java.util.ListIterator itPL = epiDataPointList.listIterator();
    while(itPL.hasNext())
    {
      EpiDataPoint epiDataPoint = (EpiDataPoint) itPL.next();
      if(epiDataPoint.isPointWithinWindow(epiFilter))
      {
        switch(level)
        {
          case SHORT:
            bw.write(epiDataPoint.printShortSummary());
            break;
          case MEDIUM:
            bw.write(epiDataPoint.printMediumSummary());
            break;
          case LONG:
            bw.write(epiDataPoint.printDetailedSummary());
            break;
          default:
            bw.write("unrecognized detail level for the print summary");
            break;
        }
        count++;
      }
    }
    return count;

  }



  /**
   * write BranchLength / MutualInformation pairs
   * @param bw
   * @throws java.io.IOException
   */
  public void writeBranchMIData(BufferedWriter bw) throws java.io.IOException
  {
    ListIterator itBM = branchMutPairList.listIterator();
    while(itBM.hasNext())
    {
      BranchMutualPair bmPair = (BranchMutualPair) itBM.next();
      bw.write(bmPair.getBranchLength() + "\t" + bmPair.getMutualInfo() + "\n");
    }
  }



  public void writeShortSummary(BufferedWriter bw) throws java.io.IOException
  {
    bw.write(EpiDataPoint.printShortHeader());
    ListIterator itEDP = this.epiDataPointList.listIterator();
    while(itEDP.hasNext())
    {
      EpiDataPoint edp = (EpiDataPoint)itEDP.next();
      bw.write(edp.printShortSummary());
    }

  }

  
  public EpiSiteMutData getSeqSiteMutData(int s1, int iReal, int jReal)
  {

    SiteMutData siteMDA = new SiteMutData();

    /******************************************************
     * SEQUENCE SITE A
     *****************************************************/

    //*** nucleotide data
    siteMDA.setNtPos(iReal);
    if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("A"))
      siteMDA.setNt(1);
    else
      if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("C"))
        siteMDA.setNt(2);
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("G"))
          siteMDA.setNt(3);
        else
          siteMDA.setNt(4);
    siteMDA.setConservation(ntAlFull.getConsensusProb()[iReal-1]);

    //*** nucleotide data Consensus sequence
    if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[iReal-1]).equalsIgnoreCase("A"))
      siteMDA.setConsNT(1);
    else
      if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[iReal-1]).equalsIgnoreCase("C"))
        siteMDA.setConsNT(2);
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[iReal-1]).equalsIgnoreCase("G"))
          siteMDA.setConsNT(3);
        else
          siteMDA.setConsNT(4);

    //*** codon data
    int codonStartA = iReal%3+3;
    siteMDA.setCodonPos(codonStartA);
    codonStartA = iReal - codonStartA; // remember, this is the pos starting from 0
    siteMDA.setCodonStart(codonStartA);


    String codonSeqA =
            ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartA))
          + ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartA+1))
          + ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartA+2));
    siteMDA.setCodon(codonSeqA);
    siteMDA.setaAcid(cu.getShort(codonSeqA));

    String codonConA =
            ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartA])
          + ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartA+1])
          + ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartA+2]);
    siteMDA.setConsCodon(codonConA);
    siteMDA.setConsAAcid(cu.getShort(codonConA));


    if(cu.getShort(codonSeqA).equals(cu.getShort(codonConA)))
      siteMDA.setSubType(SubstitutionType.SYNONYMOUS);
    else
      siteMDA.setSubType(SubstitutionType.NONSYNONYMOUS);


    /******************************************************
     * SEQUENCE SITE B
     *****************************************************/
    SiteMutData siteMDB = new SiteMutData();


    //*** nucleotide data
    siteMDB.setNtPos(jReal);
    if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("A"))
      siteMDB.setNt(1);
    else
      if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("C"))
        siteMDB.setNt(2);
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("G"))
          siteMDB.setNt(3);
        else
          siteMDB.setNt(4);
    siteMDB.setConservation(ntAlFull.getConsensusProb()[jReal-1]);


    //*** nucleotide data Consensus sequence
    if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[jReal-1]).equalsIgnoreCase("A"))
      siteMDB.setConsNT(1);
    else
      if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[jReal-1]).equalsIgnoreCase("C"))
        siteMDB.setConsNT(2);
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[jReal-1]).equalsIgnoreCase("G"))
          siteMDB.setConsNT(3);
        else
          siteMDB.setConsNT(4);

    //*** codon data
    int codonStartB = jReal%3+3;
    siteMDB.setCodonPos(codonStartB);
    codonStartB = jReal - codonStartB;
    siteMDB.setCodonStart(codonStartB);


    String codonSeqB =
            ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartB))
          + ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartB+1))
          + ntAlFull.getNum2nt().get(ntAlFull.getNumData(s1, codonStartB+2));
    siteMDB.setCodon(codonSeqB);
    siteMDB.setaAcid(cu.getShort(codonSeqB));

    String codonConB =
            ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartB])
          + ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartB+1])
          + ntAlFull.getNum2nt().get(ntAlFull.getConsensusData()[codonStartB+2]);
    siteMDB.setConsCodon(codonConB);
    siteMDB.setConsAAcid(cu.getShort(codonConB));


    if(cu.getShort(codonSeqB).equals(cu.getShort(codonConB)))
      siteMDB.setSubType(SubstitutionType.SYNONYMOUS);
    else
      siteMDB.setSubType(SubstitutionType.NONSYNONYMOUS);


    return new EpiSiteMutData(siteMDA, siteMDB);

  }

  /**
   * @return the wTree
   */
  public WivTree getwTree() {
    return wTree;
  }

  /**
   * @return the epiProps
   */
  public no.uio.medisin.bag.core.alignment.AlignPropCount getEpiProps() {
    return epiData;
  }

  /**
   * @return the ntAlFull
   */
  public NumericalAlignment getNtAlFull() {
    return ntAlFull;
  }

  /**
   * @return the branchMutPairList
   */
  public java.util.ArrayList<BranchMutualPair> getBranchMutPairList() {
    return branchMutPairList;
  }

  /**
   * @return the range
   */
  public EpiFilterSet getRange() {
    return range;
  }


}