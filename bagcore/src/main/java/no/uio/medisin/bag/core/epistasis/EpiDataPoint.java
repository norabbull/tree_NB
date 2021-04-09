/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import no.uio.medisin.bag.core.alignment.NumericalAlignment;
import no.uio.medisin.bag.core.aminoacid.CodonUsage;
import no.uio.medisin.bag.core.tree.NodeType;
import no.uio.medisin.bag.core.tree.SiteMutData;
import no.uio.medisin.bag.core.tree.SubstitutionType;
import no.uio.medisin.bag.core.tree.WivTree;
import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.apache.commons.math3.stat.descriptive.rank.Min;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;


/**
 *
 * @author sr
 */
public class EpiDataPoint {

  private int iSite;
  private int jSite;
  private int iReal;
  private int jReal;
  private java.util.ArrayList<EpiSiteMutData> epiEpiMutDataList;
  private java.util.ArrayList<Double> branchLengthList;
  private double evDistArray[][];
  private int ntPairFreq[][];
  private int alignDist;
  private double entropyI;
  private double entropyJ;
  private double entropyIJ;
  private double mutualIJ;
  private double blSD;
  private double blMean;
  private double blMoment;
  private double blMedian;
  private double blMax;
  private double blMin;
  private double blSpread;
  private int synCountA;
  private int synCountB;
  private int epiSeqCount; // number of epistatic sequences
  private int oneMutSeqCount; // number of sequences with one mut
  private CodonUsage cu;
  private EpiMutationStats mutStats;


  // also need to calculate the MI for the dataset

  public EpiDataPoint(int x, int y, int numSeqs)
  {
    iSite = x;
    jSite = y;
    epiEpiMutDataList = new java.util.ArrayList<EpiSiteMutData>();
    branchLengthList = new java.util.ArrayList<Double>();
    ntPairFreq = new int[4][4];
    evDistArray = new double[numSeqs][numSeqs];
    cu = new CodonUsage();
    mutStats = new EpiMutationStats();
  }

  /**
   *
   * this retrieves a number of properties that describe an epistatic pair
   * specifically, for each sequence it finds"
   * 1. nucleotide
   * 2. position
   * 3. tip/internal node
   * 4. codon
   * 5. symmetric/non-symmetric substitution
   *
   * for 4 & 5, need the whole sequence, rather than the one that was generated
   * to only contain the invariant sites
   *
   * this is a little more fiddly because we have to use the site
   * data to map consensusSet (which only contains invariant sites)
   * on to ntFull which contains all sites.
   *
   * This may seem really dumb, but we remove the invariant sites to
   * speed things up; however, we need the full alignment to calculate anything
   * to do with the codons.
   *
   * first index is the sequence, second index is the site
   * i.e. numData[seq][site]
   *
   * @param na
   */
  public void calcEpiDataPtProperties(
          NumericalAlignment ntAlFull,
          WivTree wTree,
          int iReal, int jReal)
  {
    //System.out.println(iReal + "|" + jReal);

   this.setiReal(iReal);
   this.setjReal(jReal);
   this.setAlignDist(java.lang.Math.abs(iReal-jReal));

    // Run through all the sequences and collect the mutation information
    for(int s1=0; s1<ntAlFull.getNumSeqs(); s1++)
    {
      int iNT = 0;
      if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("A"))
        iNT = 0;
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("C"))
          iNT = 1;
        else
          if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][iReal-1]).equalsIgnoreCase("G"))
            iNT = 2;
          else
            iNT = 3;
      int jNT = 0;
      if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("A"))
        jNT = 0;
      else
        if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("C"))
          jNT = 1;
        else
          if(ntAlFull.getNum2nt().get(ntAlFull.getNumData()[s1][jReal-1]).equalsIgnoreCase("G"))
            jNT = 2;
          else
            jNT = 3;
      this.incCount(iNT, jNT);
      EpiSiteMutData seqEpiMutData = this.getSeqSiteMutData(ntAlFull, s1, iReal, jReal);
      this.getEpiPairDataList().add(seqEpiMutData);
    }

    // now cycle through the sequence pairs and calculate branch length
    // and get mutation pairs.
    int e = 0;
    for(int s1=0; s1<ntAlFull.getNumSeqs(); s1++)
    {
      // we only want to get the branch lengths for epistatic sequences
      // so skip it if it isn't
      if(this.getEpiPairDataList().get(s1).isSeqEpistatic() == false)
        continue;

      String seqNameA = ntAlFull.getSeqNames()[s1];
      this.getEpiPairDataList().get(s1).setSeqName(seqNameA);
      // check whether the sequence is internal node or tip
      if(seqNameA.indexOf("NODE") == 0)
        this.getEpiPairDataList().get(s1).setNodeLocation(NodeType.INTERNAL);
      else
        this.getEpiPairDataList().get(s1).setNodeLocation(NodeType.EXTERNAL);

      // cycle through the sequence pairs and calculate the branch lengths
      // between the epistatic ones
      for(int s2=s1+1; s2<ntAlFull.getNumSeqs(); s2++)
      {
        if(this.getEpiPairDataList().get(s2).isSeqEpistatic() == false)
          continue;
        e++;
        String seqNameB = ntAlFull.getSeqNames()[s2];
        double branchLength = wTree.findPairDistance(
                seqNameA,
                seqNameB);
        this.getEvDistArray()[s1][s2] = branchLength;
        this.getEvDistArray()[s2][s1] = branchLength;
        if(branchLength > 0)
          this.getBranchLengthList().add(new Double(branchLength));

      }

    }
    //System.out.println("no of branch lengths = " + epiDP.getBranchLengthList().size());
    //System.out.println("epi="+ e);
    this.calculateMutualInfo();
    this.countEpiSeqs();
    this.calcBranchLengthStats();
    this.calculateSiteStats();
  }

  public EpiSiteMutData getSeqSiteMutData(NumericalAlignment ntAlFull,
          int s1, int iReal, int jReal)
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
    int codonStartA = iReal%3;
    if(codonStartA == 0)
      codonStartA = 3;
    else
      if(codonStartA == 1)
        codonStartA = 3;
      else
        codonStartA = 2;

    if(codonStartA > 3)
      System.out.println("weird codon pos" + codonStartA);
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
    int codonStartB = iReal%3;
    if(codonStartB == 0)
      codonStartB = 3;
    else
      if(codonStartB == 1)
        codonStartB = 3;
      else
        codonStartB = 2;

    if(codonStartB > 3)
      System.out.println("weird codon pos" + codonStartB);
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


  public String printSummary()
  {
    String s = iSite + "\t" + jSite + "\t"
            + alignDist + "\t" + entropyI + "\t" + entropyI + "\t"
            + entropyIJ + "\t" + mutualIJ + "\t"
            + getBlSD() + "\n";
    return s;
  }



  /**
   * SMetric = (double)epiSeqCount/[(double)oneMutSeqCount+1] is supposed
   * to give some measure of how much of the dinucleotide variation is
   * associated with epistatic sequences, as opposed to sequences with
   * a mutation at only one of the sites.  We use this method instead of
   * the more commmonly used Mutual Information because the epistatic sequences
   * in generally only represent a small fraction of the total number of
   * sequences (10 - 20% is typical) so the MI value is dominated by the
   * contribution from the conserved sequences.
   * 
   * @return
   */
  public double getSMetric()
  {
    return (double)epiSeqCount/((double)oneMutSeqCount+1);
  }


  /**
   * header line for short information
   * @return
   */
  public static String printShortHeader()
  {
    String s = "iSite" + "\t" + "jSite" + "\t"
            + "blCount" +"\t"
            + "blSD" + "\t" + "blMean"  + "\t" + "blMoment" + "\t"
            + "blMedian" + "\t "+ "blMin" + "\t" + "blMax" + "\t"
            + "[blMax-blMin]" + "\t"
            + "alignDist" + "\t"
            + "epiSeqC" + "\t" + "MutSeqC" + "\t"
            + "epiSeqC/[MutSeqC + 1]" + "\t"
            + "entropyI" + "\t" + "entropyJ" + "\t"
            + "entropyIJ" + "\t" + "mutualIJ"
            + "\n";
    return s;
  }

  /**
   * short information string for the data point
   * only lists the calculated variables 
   * @return
   */
  public String printShortSummary()
  {
    NumberFormat formatter = new DecimalFormat("0.00000");
    String s = iReal + "\t" + jReal + "\t"
            + branchLengthList.size() +"\t"
            + formatter.format(getBlSD()) + "\t" + formatter.format(getBlMean())
            + "\t" + formatter.format(getBlMoment()) + "\t"
            + formatter.format(getBlMedian()) + "\t "
            + formatter.format(getBlMin()) + "\t" + formatter.format(getBlMax()) + "\t"
            + formatter.format(getBlMax()-getBlMin()) + "\t"
            + alignDist + "\t"
            + epiSeqCount + "\t" + oneMutSeqCount + "\t"
//            + epiSeqCount + "\t" + oneMutSeqCount + "\t"
            + formatter.format((double)epiSeqCount/((double)oneMutSeqCount + 1.0)) + "\t"
            + formatter.format(entropyI) + "\t" + formatter.format(entropyJ) + "\t"
            + formatter.format(entropyIJ) + "\t" + formatter.format(mutualIJ)
            + "\n";
    return s;
  }

  /**
   * provides a more detailed description of the data points including
   * 
   *   short summary
   *   the epistatic sequences at this site
   *   the dinucleotide pair frequencies
   *   a summary of the nt mutations, codons and amino acid changes at the site
   *   all the branch lengths
   * 
   * @return
   */
  public String printMediumSummary()
  {
    String s="";
    s = s.concat(printShortHeader());
    s = s.concat(this.printShortSummary() + "\n\n");
    s = s.concat(this.printEpiSequenceNames() + "\n\n");
    s = s.concat(this.printPairFrequency() + "\n\n");
    s = s.concat(this.printMutationSummary() + "\n\n");
    s = s.concat(this.printBranchLengths(10) + "\n\n");
    s = s.concat("------------------------------------------" +
            "------------------------------------------\n\n\n\n");

    return s;
  }

  public String printDetailedSummary()
  {
    String s="";
    return s;
  }
  /**
   * prints a tab delimited string of all the branch lengths
   *
   * @param  lineCount - number of entries / line
   * @return
   */
  public String printBranchLengths(int lineCount)
  {
    NumberFormat formatter = new DecimalFormat("0.00000");
    String s = "\tBranch Lengths\n";
    int b=0;
    java.util.ListIterator itBL = branchLengthList.listIterator();
    while(itBL.hasNext())
    {
      if(b%10==0) s = s.concat("\n");
      Double d = (Double) itBL.next();
      s = s.concat("\t" + formatter.format(d).toString());
      b++;
    }
    return s;
  }


  /**
   * get a list of the names of all the epistatic sequences
   *
   * @return
   */
  public String printEpiSequenceNames()
  {
    java.util.ListIterator itMu = epiEpiMutDataList.listIterator();
    String s = "";
    while(itMu.hasNext())
    {
      EpiSiteMutData mutData = (EpiSiteMutData) itMu.next();
      if(mutData.getSeqName()!=null)
        s = s.concat("\t" + mutData.getSeqName() + "\n");
    }
    return s;
  }


  /**
   * print a summary of the mutations that occurred for this site.
   * Summarizes codon usage, synonymous/non-synonymous counts and lists
   * all codon substitutions that occurred.
   * 
   * @return
   */
  public String printMutationSummary()
  {
    String s;
    s = "\t\t\tSiteA\t\t\t\tSiteB\n"
        +    "\tCodon Pos:\t1\t2\t3\t\t1\t2\t3\n"
        +    "\tcount\t"
        + mutStats.getSiteAcodon1count() + "\t"
        + mutStats.getSiteAcodon2count() + "\t"
        + mutStats.getSiteAcodon3count() + "\t"
        + "\t"
        + mutStats.getSiteBcodon1count() + "\t"
        + mutStats.getSiteBcodon2count() + "\t"
        + mutStats.getSiteBcodon3count() + "\t"
        + "\n";


    String SiteAconsAA
            = ((SiteMutData)epiEpiMutDataList.get(0).getSeqAEpiData()).getConsAAcid();
    String SiteBconsAA
            = ((SiteMutData)epiEpiMutDataList.get(0).getSeqBEpiData()).getConsAAcid();

    s = s.concat("\n\t\t\t\tAA1\t\t\t\tAA2\t\tcount:\n");
    java.util.ListIterator itAA = mutStats.getAaPairs().getAaPairs().listIterator();
    while(itAA.hasNext())
    {
      EpistaticAAPair aaPair = (EpistaticAAPair)itAA.next();
      s = s.concat("\t\t\t" + SiteAconsAA + "\t->\t" + aaPair.getAa1()
              + "\t"
              + "\t" + SiteBconsAA + "\t->\t" + aaPair.getAa2()
              + "\t" + aaPair.getCount() + "\n");
    }
    return s;
  }

  /**
   * prints the dinucleotide frequency array for the point
   * in a 4x4 format
   * @return
   */
  public String printPairFrequency()
  {

    String s = "\tdinucleotide frequencies\n\t\tA\tC\tG\tT\n";

      int i=0;
      s = s.concat("\tA");
      s = s.concat("\t" + ntPairFreq[i][0]
              + "\t" + ntPairFreq[i][1]
              + "\t" + ntPairFreq[i][2]
              + "\t" + ntPairFreq[i][3] + "\n");
      i=1;
      s = s.concat("\tC");
      s = s.concat("\t" + ntPairFreq[i][0]
              + "\t" + ntPairFreq[i][1]
              + "\t" + ntPairFreq[i][2]
              + "\t" + ntPairFreq[i][3] + "\n");
      i=2;
      s = s.concat("\tG");
      s = s.concat("\t" + ntPairFreq[i][0]
              + "\t" + ntPairFreq[i][1]
              + "\t" + ntPairFreq[i][2]
              + "\t" + ntPairFreq[i][3] + "\n");
      i=3;
      s = s.concat("\tT");
      s = s.concat("\t" + ntPairFreq[i][0]
              + "\t" + ntPairFreq[i][1]
              + "\t" + ntPairFreq[i][2]
              + "\t" + ntPairFreq[i][3] + "\n");
    return s;
  }


  /**
   * check whether this point passes the filter
   * 
   * @param epiFilter
   * @return
   */
  public Boolean isPointWithinWindow(no.uio.medisin.bag.core.epistasis.EpiFilterSet epiFilter)
  {
    /*
     * yeah, i know, i could write this as a single 'if' statement
     * but it's more readable this way.
     */

    if(epiFilter.getMinSMetric() <= this.getSMetric() 
            && epiFilter.getMaxSMetric() >= this.getSMetric())
    {
      if(epiFilter.getMinEpiSeqCount() <= this.getEpiSeqCount()
              && epiFilter.getMaxEpiSeqCount() >= this.getEpiSeqCount())
      {
        if(epiFilter.getMinMutSeqCount() <= this.getOneMutSeqCount()
                && epiFilter.getMaxMutSeqCount() >= this.getOneMutSeqCount())
        {
          if(epiFilter.getMinAlignDist() <= this.getAlignDist()
                  && epiFilter.getMaxAlignDist() >= this.getAlignDist())
          {
            if(epiFilter.getMinBLcount() <= this.getBranchLengthList().size()
                    && epiFilter.getMaxBLcount() >= this.getBranchLengthList().size())
            {
              if(epiFilter.getMinBL() <= this.getBlMin()
                      && epiFilter.getMaxBL() >= this.getBlMax())
              {
                if(epiFilter.getMinBLspread() <= this.getBlSpread()
                        && epiFilter.getMaxBLspread() >= this.getBlSpread())
                {
                  if(epiFilter.getMinBLmean() <= this.getBlMean()
                          && epiFilter.getMaxBLmean() >= this.getBlMean())
                  {
                    if(epiFilter.getMinBLsd() <= this.getBlSD()
                            && epiFilter.getMaxBLsd() >= this.getBlSD())
                    {
                      if(epiFilter.getMinBLmoment() <= this.getBlMoment()
                              && epiFilter.getMaxBLmoment() >= this.getBlMoment())
                      {
                        return true;
                      }
                      else
                        return false;
                    }
                    else
                      return false;
                  }
                  else
                    return false;
                }
                else
                  return false;
              }
              else
                return false;
            }
            else
              return false;
          }
          else
            return false;
        }
        else
          return false;
      }
      else
        return false;
    }
    else
      return false;
  }



  /**
   * count up how many epistatic sequences have synonymous substitutions
   *
   */
  public void countSynSubstitutions()
  {
    java.util.ListIterator itMu = epiEpiMutDataList.listIterator();
    synCountA = 0;
    synCountB = 0;
    while(itMu.hasNext())
    {
      EpiSiteMutData currEpiMutData = (EpiSiteMutData)itMu.next();
      if(currEpiMutData.getSeqAEpiData().getSubType() == SubstitutionType.SYNONYMOUS)
        synCountA++;
      if(currEpiMutData.getSeqBEpiData().getSubType() == SubstitutionType.SYNONYMOUS)
        synCountB++;

    }
    System.out.println("of the " + epiEpiMutDataList + " epistatic sequences");
    System.out.println("there are " + synCountA + " at A &  "+ synCountB + " at B");

  }


  public void calcBranchLengthStats()
  {
    double eDist[] = new double[branchLengthList.size()];
    java.util.ListIterator itBL = branchLengthList.listIterator();
    int d=0;
    while(itBL.hasNext())
    {
      eDist[d++] = (Double) itBL.next();
    }
    StandardDeviation sd = new StandardDeviation();
    blSD = sd.evaluate(eDist);

    Mean mean = new Mean();
    blMean = mean.evaluate(eDist);

    Median median = new Median();
    blMedian = median.evaluate(eDist);

    Max max = new Max();
    blMax = max.evaluate(eDist);

    Min min = new Min();
    blMin = min.evaluate(eDist);

    if(getBlMin() < 0)
    {
      System.out.println("why");
    }

  }
  
  


  /**
   * calculate the mutual information and joint entropy based on the di-nucleotide
   * distribution
   * 
   */
  public void calculateMutualInfo()
  {

    // first calculate the entropy
    java.util.ListIterator itMu = epiEpiMutDataList.listIterator();
    int [] ntAcount = new int[4];
    int [] ntBcount = new int[4];
    while(itMu.hasNext())
    {
      EpiSiteMutData currEpiMutData = (EpiSiteMutData)itMu.next();
      ntAcount[currEpiMutData.getSeqAEpiData().getNt()-1]++;
      ntBcount[currEpiMutData.getSeqBEpiData().getNt()-1]++;
    }
    for(int i=0; i<4; i++)
    {
      double pA = (double)ntAcount[i]/(double)epiEpiMutDataList.size();
      if(pA == 0)
        entropyI -= pA /java.lang.Math.log10(2);
      else
        entropyI -= pA * java.lang.Math.log10(pA) /java.lang.Math.log10(2);

      double pB = (double)ntBcount[i]/(double)epiEpiMutDataList.size();
      if(pB == 0)
        entropyJ -= pB /java.lang.Math.log10(2);
      else
        entropyJ -= pB * java.lang.Math.log10(pB) /java.lang.Math.log10(2);
    }

    // now calculate joint entropy H(i,j)
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
      {
        double pAB = (double)ntPairFreq[i][j]/(double)epiEpiMutDataList.size();
        if(pAB==0)
          entropyIJ -= pAB /java.lang.Math.log10(2);
        else
          entropyIJ -= pAB*java.lang.Math.log10(pAB) /java.lang.Math.log10(2);
      }
    }

    // finally, calculate the MutualInfo
    mutualIJ = entropyI + entropyJ - entropyIJ;

  }


  /**
   * get the stats for the codon and amino acid usage for this site
   * 
   */
  public void calculateSiteStats()
  {
    java.util.ListIterator itMu = epiEpiMutDataList.listIterator();

    while(itMu.hasNext())
    {
      EpiSiteMutData mutData = (EpiSiteMutData) itMu.next();
      if(mutData.isSeqEpistatic())
      {
        /*************************************************
         * Site A - Non-Syn/Syn substitutions
         ************************************************/
        if(mutData.getSeqAEpiData().getaAcid().equalsIgnoreCase(
                mutData.getSeqBEpiData().getConsAAcid()))
          mutStats.incSiteASynCount();
        else
        {
          mutStats.incSiteANonSynCount();

         // Amino Acid Transitions
        int index = mutStats.getAaPairs().contains(
                mutData.getSeqAEpiData().getaAcid(),
                mutData.getSeqBEpiData().getaAcid());
        if(index >= 0)
          mutStats.getAaPairs().incCount(index);
        else
          mutStats.addPair(new EpistaticAAPair(mutData.getSeqAEpiData().getaAcid(),
                  mutData.getSeqBEpiData().getaAcid()));
        }

        /*************************************************
         * Site A - Codon Usage
         ************************************************/
        switch(mutData.getSeqAEpiData().getCodonPos())
        {
          case 1:
            mutStats.incSiteACodon1Count();
            break;
          case 2:
            mutStats.incSiteACodon2Count();
            break;
          case 3:
            mutStats.incSiteACodon3Count();
            break;
        }
        
        /*************************************************
         * Site B - Non-Syn/Syn substitutions
         ************************************************/
        if(mutData.getSeqBEpiData().getaAcid().equalsIgnoreCase(
                mutData.getSeqBEpiData().getConsAAcid()))
          mutStats.incSiteBSynCount();
        else
        {
          mutStats.incSiteBNonSynCount();

         // Amino Acid Transitions
        int index = mutStats.getAaPairs().contains(
                mutData.getSeqBEpiData().getaAcid(),
                mutData.getSeqBEpiData().getaAcid());
        if(index >= 0)
          mutStats.getAaPairs().incCount(index);
        else
          mutStats.addPair(new EpistaticAAPair(mutData.getSeqBEpiData().getaAcid(),
                  mutData.getSeqBEpiData().getaAcid()));
        }

        /*************************************************
         * Site B - Codon Usage
         ************************************************/
        switch(mutData.getSeqBEpiData().getCodonPos())
        {
          case 1:
            mutStats.incSiteBCodon1Count();
            break;
          case 2:
            mutStats.incSiteBCodon2Count();
            break;
          case 3:
            mutStats.incSiteBCodon3Count();
            break;
        }
      }
    }

  }

  /**
   * find how many of the sequences are epistatic, and how many contain a single
   * mutation
   * 
   */
  public void countEpiSeqs()
  {

    int max = 0;
    int iMax = 0;
    int jMax = 0;

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        if(ntPairFreq[i][j]> max)
        {
          max = ntPairFreq[i][j];
          iMax = i;
          jMax = j;
        }


    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        if(i == iMax && j == jMax)
          continue;
        else
        {
          if(i == iMax || j==jMax)
            oneMutSeqCount += ntPairFreq[i][j];
          else
            epiSeqCount+= ntPairFreq[i][j];
        }
   // System.out.println(oneMutSeqCount + "|" + epiSeqCount);
  }

  public void incCount(int i, int j)
  {
    ntPairFreq[i][j]++;
  }

  /**
   * @return the i
   */
  public int getI() {
    return iSite;
  }

  /**
   * @return the j
   */
  public int getJ() {
    return jSite;
  }

  /**
   * @return the epiPairDataList
   */
  public java.util.ArrayList<EpiSiteMutData> getEpiPairDataList() {
    return epiEpiMutDataList;
  }

  /**
   * @return the evolDistList
   */
  public java.util.ArrayList<Double> getBranchLengthList() {
    return branchLengthList;
  }

  /**
   * @return the ntPairFreq
   */
  public int[][] getNtPairFreq() {
    return ntPairFreq;
  }

  /**
   * @return the alignDist
   */
  public int getAlignDist() {
    return alignDist;
  }

  /**
   * @param alignDist the alignDist to set
   */
  public void setAlignDist(int alignDist) {
    this.alignDist = alignDist;
  }

  /**
   * @return the evDistArray
   */
  public double[][] getEvDistArray() {
    return evDistArray;
  }

  /**
   * @return the entropyI
   */
  public double getEntropyI() {
    return entropyI;
  }

  /**
   * @return the entropyJ
   */
  public double getEntropyJ() {
    return entropyJ;
  }

  /**
   * @return the entropyIJ
   */
  public double getEntropyIJ() {
    return entropyIJ;
  }

  /**
   * @return the mutualIJ
   */
  public double getMutualIJ() {
    return mutualIJ;
  }

  /**
   * @param sd the sd to set
   */
  public void setSd(double sd) {
    this.blSD = sd;
  }

  /**
   * @return the sd
   */
  public double getSd() {
    return getBlSD();
  }

  /**
   * @return the synCountA
   */
  public int getSynCountA() {
    return synCountA;
  }

  /**
   * @return the synCountB
   */
  public int getSynCountB() {
    return synCountB;
  }

  /**
   * @return the epiSeqCount
   */
  public int getEpiSeqCount() {
    return epiSeqCount;
  }

  /**
   * @return the oneMutSeqCount
   */
  public int getOneMutSeqCount() {
    return oneMutSeqCount;
  }

  /**
   * @return the blSD
   */
  public double getBlSD() {
    return blSD;
  }

  /**
   * @return the blMean
   */
  public double getBlMean() {
    return blMean;
  }

  /**
   * @return the blMoment
   */
  public double getBlMoment() {
    return blMoment;
  }

  /**
   * @return the blMedian
   */
  public double getBlMedian() {
    return blMedian;
  }

  /**
   * @return the blMax
   */
  public double getBlMax() {
    return blMax;
  }

  /**
   * @return the blMin
   */
  public double getBlMin() {
    return blMin;
  }

  /**
   * @return the blSpread
   */
  public double getBlSpread() {
    return blMax-blMin;
  }

  /**
   * @return the iReal
   */
  public int getiReal() {
    return iReal;
  }

  /**
   * @param iReal the iReal to set
   */
  public void setiReal(int iReal) {
    this.iReal = iReal;
  }

  /**
   * @return the jReal
   */
  public int getjReal() {
    return jReal;
  }

  /**
   * @param jReal the jReal to set
   */
  public void setjReal(int jReal) {
    this.jReal = jReal;
  }

}
