/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.alignment;

/**
 * This provides a completely numerical representation of an alignment
 * permitting more rapid calculation of alignment properties.
 * The class stores the alignment as a 2D N x L array where N is
 * the number of sequences and L is the length of the alignment.
 * Since the principal purpose of the class is to calculate properties
 * such as Entropy and Mutual Information, including invariant sites serves
 * no purpose since they provide no information and merely extend the
 * calculation time unnecessarily.  Thus a Sites[L] can also be included
 * which stores the site position in the original alignment.
 * Finally, structural information about the alignment can also be defined
 * in the Structure [L][] array which stores information about which sites
 * are in close physical proximity to each site and should therefore be
 * incorporated into the MutualInformation calculation for that site.
 *
 * So, for example, if site 5 has Structure[5][0] = 10 & Structure[5][1] = 15
 * this means that changes occuring at site 10 and 15 should be considered
 * as equivalent mutations since site 10 & 15 are known to be in close
 * physical proximity in the known structure.
 * @author sr
 */

import no.uio.medisin.bag.core.epistasis.EpistasisColumn;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

public class NumericalAlignment implements Serializable{
  private int alignLength;
  private int numSequences;
  private int maxStructures;
  private int numDistinctMols; // How many different types of molecules are
                                // in the alignment.
  private int numData[][];
  private int consensusData[];
  private double consensusProb[];
  private int sites[];
  private int structData[][];
  private double dist;
  private String seqNames[];
  private double distMat[][];
  private Hashtable<String, Integer> nt2num;
  private Hashtable<Integer, String> num2nt;

  private double entropy[];
  private double mutualInfo[][];
  private double sharedEntropy[][];
  private double averagePosCorr[][];
  private double averageColMI[];
  private double averageEntropy;
  private double averageMI;
  private double prob[][];
  private java.util.ArrayList<Integer> acidList;// stores a list of all AAs
                                                // in the alignment
//  private java.util.ArrayList<MutationPair> epiSites;
  private ArrayList<EpistasisColumn> epiData;
  private int epiDataNum[][];
  private int minEpi;
  private int maxEpi;

  private boolean haveAlignment;
  private boolean haveAlignData;
  private boolean haveEntropy;
  private boolean haveMutualInfo;
  private boolean haveBootstraps;
  private boolean haveTree;

  // need to add map for AA <-> numerical representation

  /**
   * Basic Constructor
   */
  public NumericalAlignment()
  {
    alignLength = 0;
    numSequences = 0;
    maxStructures = 0;
    acidList = new java.util.ArrayList<Integer>();
    epiData = new ArrayList<EpistasisColumn>();
  }

  /**
   * This constructor initializes the data arrays
   *
   * @param aL  Length of alignment
   * @param nS  Number of Sequences
   * @param mS  Maximum number of atoms mapped to any one site.
   */
  public NumericalAlignment(int aL, int nS, int mS)
  {
    alignLength = aL;
    numSequences = nS;
    maxStructures = mS;
    dist = -1.0;

    numData = new int[nS][aL];
    consensusData = new int[aL];
    consensusProb = new double[aL];
    sites = new int[aL];
    //structData = new int[aL][mS];
    seqNames = new String[nS];
    //distMat = new double[aL][aL];
    nt2num = new Hashtable<String, Integer>();
    num2nt = new Hashtable<Integer, String>();

    //entropy = new double[aL];
    //mutualInfo = new double[aL][aL];
    //sharedEntropy = new double[aL][aL];
    //prob = new double[aL][aL];
    acidList = new java.util.ArrayList<Integer>();

    epiData = new ArrayList<EpistasisColumn>();
    epiDataNum = new int[aL][aL];
    //averagePosCorr = new double [aL][aL];
    //averageColMI = new double [aL];

  }

  public NumericalAlignment(
          int [][] seq,
          String[] names,
          int[] s,
          int[][] struct,
          Hashtable<String, Integer> aaMap)
  {
    alignLength = seq[0].length;
    numSequences = names.length;
    dist = -1.0;

    numData = new int[numSequences][alignLength];
    consensusData = new int[alignLength];
    consensusProb = new double[alignLength];
    sites = new int[alignLength];
    //structData = new int[alignLength][numSequences];
    seqNames = new String[numSequences];

    //entropy = new double[alignLength];
    //mutualInfo = new double[alignLength][alignLength];
    //sharedEntropy = new double[alignLength][alignLength];
    //averagePosCorr = new double [alignLength][alignLength];
    //averageColMI = new double [alignLength];
    //prob = new double[alignLength][alignLength];
    acidList = new java.util.ArrayList<Integer>();
    epiData = new ArrayList<EpistasisColumn>();
    nt2num = new Hashtable<String, Integer>();
    num2nt = new Hashtable<Integer, String>();

    int i=0;
    while(i<numSequences)
    {
      int j=0;
      while(j<alignLength)
      {
        numData[i][j] = seq[i][j];
        structData[j][i] = struct[j][i];
        j++;
      }
      seqNames[i] = names[i];
      i++;
    }

    int k=0;
    while(k<alignLength)
    {
      sites[k] = s[k];
      k++;
    }

    nt2num.putAll(aaMap);
    findDistinctMols();
  }

  public NumericalAlignment(
          int [][] seq,
          String[] names,
          int[] s,
          Hashtable<String, Integer> aaMap)
  {
    alignLength = seq[0].length;
    numSequences = names.length;
    dist = -1.0;

    numData = new int[numSequences][alignLength];
    consensusData = new int[alignLength];
    consensusProb = new double[alignLength];
    sites = new int[alignLength];
    //structData = new int[alignLength][numSequences];
    seqNames = new String[numSequences];

    entropy = new double[alignLength];
    mutualInfo = new double[alignLength][alignLength];
    sharedEntropy = new double[alignLength][alignLength];
    averagePosCorr = new double [alignLength][alignLength];
    averageColMI = new double [alignLength];
    prob = new double[alignLength][alignLength];
    acidList = new java.util.ArrayList<Integer>();
    epiData = new ArrayList<EpistasisColumn>();
    nt2num = new Hashtable<String, Integer>();
    num2nt = new Hashtable<Integer, String>();

    int i=0;
    while(i<numSequences)
    {
      int j=0;
      while(j<alignLength)
      {
        numData[i][j] = seq[i][j];
        j++;
      }
      seqNames[i] = names[i];
      i++;
    }

    int k=0;
    while(k<alignLength)
    {
      sites[k] = s[k];
      k++;
    }

    nt2num.putAll(aaMap);
    findDistinctMols();
  }

  /**
   * create new NumericalAlignment from a WivAAAlignment
   * @param wAlign - a WivAlignment
   */
  public NumericalAlignment(WivAAAlignment wAAAlign)
  {

    alignLength = wAAAlign.getAlignmentLength();
    numSequences = wAAAlign.getNoOfSeqs();
    maxStructures = wAAAlign.getMaxStruct();
    dist = -1.0;

    numData = new int[numSequences][alignLength];
    consensusData = new int[alignLength];
    consensusProb = new double[alignLength];
    java.util.ArrayList<Integer> wSites = wAAAlign.getSitesList();
    sites = new int[alignLength];
    java.util.ListIterator itS = wAAAlign.getSitesList().listIterator();
    int s=0;
    while(itS.hasNext())
    {
      sites[s++] = (Integer)itS.next();
    }
    //structData = new int[alignLength][numSequences];
    seqNames = new String[numSequences];

    entropy = new double[alignLength];
    mutualInfo = new double[alignLength][alignLength];
    sharedEntropy = new double [alignLength][alignLength];
    prob = new double[alignLength][alignLength];
    averagePosCorr = new double [alignLength][alignLength];
    averageColMI = new double [alignLength];
    acidList = new java.util.ArrayList<Integer>();
    epiData = new ArrayList<EpistasisColumn>();
    nt2num = new Hashtable<String, Integer>();
    num2nt = new Hashtable<Integer, String>();
    aaAlignment2Numbers(wAAAlign);
    findDistinctMols();
  }

  /**
   * create new NumericalAlignment from a WivAAAlignment
   * @param wAlign - a WivAlignment
   */
  public NumericalAlignment(WivNTAlignment wNTAlign)
  {

    alignLength = wNTAlign.getAlignmentLength();
    numSequences = wNTAlign.getNoOfSeqs();
    maxStructures = wNTAlign.getMaxStruct();
    dist = -1.0;

    numData = new int[numSequences][alignLength];
    consensusData = new int[alignLength];
    consensusProb = new double[alignLength];

    java.util.ArrayList<Integer> wSites = wNTAlign.getSitesList();
    sites = new int[alignLength];
    java.util.ListIterator itS = wNTAlign.getSitesList().listIterator();
    int s=0;
    while(itS.hasNext())
    {
      sites[s++] = (Integer)itS.next();
    }
    //structData = new int[alignLength][numSequences];
    seqNames = new String[numSequences];

    //entropy = new double[alignLength];
    //mutualInfo = new double[alignLength][alignLength];
    //sharedEntropy = new double [alignLength][alignLength];
    //prob = new double[alignLength][alignLength];
    //averagePosCorr = new double [alignLength][alignLength];
    //averageColMI = new double [alignLength];
    acidList = new java.util.ArrayList<Integer>();
    epiData = new ArrayList<EpistasisColumn>();
    epiDataNum = new int[alignLength][alignLength];
    nt2num = new Hashtable<String, Integer>();
    num2nt = new Hashtable<Integer, String>();
    ntAlignment2Numbers(wNTAlign);
    copyNames(wNTAlign);
    findDistinctMols();
  }

  /**
   * create new Numerical Alignment from a WivAlignment
   * and an AminoAcid distance matrix obtained from the
   * corresponding PDB protein structure file.
   *
   * @param wAlign  a WivAlignment
   * @param distMat double [][] containing distances.
   */
  public NumericalAlignment(WivAAAlignment wAlign, double dMat[][])
  {
    this(wAlign);
    distMat = dMat;
    aaAlignment2Numbers(wAlign);
    findDistinctMols();
  }

  public void copyNames(WivNTAlignment ntAlign)
  {
    java.util.ListIterator itNm =  ntAlign.getNames().listIterator();
    int s=0;
    while(itNm.hasNext())
    {
      this.seqNames[s++] = (String)itNm.next();
    }
  }
  /**
   * This does the actual translation into numbers from a text alignment.
   * It generates a numerical representation of the array plus a hashTable
   * so we can translate back to a text alignment to interpret any analysis.
   *
   * @param wA
   * @return void
   */
  public void aaAlignment2Numbers(WivAAAlignment wA)
  {
    int uniqueAA = 0;
    java.util.ListIterator itSeq = wA.getSeq().listIterator();
    int s=0;
    PoorAASequence currSeq = new PoorAASequence();
    while(itSeq.hasNext())
    {
      currSeq = (PoorAASequence)itSeq.next();
      for(int b=1; b<=currSeq.getLength(); b++)
      {
        if(getNt2num().containsKey(currSeq.subSeq(b, b)))
        {
          numData[s][b-1] = getNt2num().get(currSeq.subSeq(b, b));
        }
        else
        {
          getNt2num().put(currSeq.subSeq(b, b), ++uniqueAA);
          getNum2nt().put(getNt2num().get(currSeq.subSeq(b, b)), currSeq.subSeq(b, b));
          numData[s][b-1] = getNt2num().get(currSeq.subSeq(b, b));
        }
      }
      s++;
    }

    PoorAASequence aaSeq = wA.getConsensus();
    for(int b=1; b<=currSeq.getLength(); b++)
    {
      consensusData[b-1] = getNt2num().get(aaSeq.subSeq(b, b));
    }
  }

  /**
   * This does the actual translation into numbers from a text alignment.
   * It generates a numerical representation of the array plus a hashTable
   * so we can translate back to a text alignment to interpret any analysis.
   *
   * @param wA
   * @return void
   */
  public void ntAlignment2Numbers(WivNTAlignment wA)
  {
    int uniqueNT = 0;
    java.util.ListIterator itSeq = wA.getSeq().listIterator();
    int s=0;
    PoorNTSequence currSeq = new PoorNTSequence();
    while(itSeq.hasNext())
    {
      currSeq = (PoorNTSequence)itSeq.next();
      //System.out.println(currSeq.getSequence());
      for(int b=1; b<=currSeq.getLength(); b++)
      {
        if(getNt2num().containsKey(currSeq.subSeq(b, b)))
        {
          numData[s][b-1] = getNt2num().get(currSeq.subSeq(b, b));
        }
        else
        {
          getNt2num().put(currSeq.subSeq(b, b), ++uniqueNT);
          getNum2nt().put(getNt2num().get(currSeq.subSeq(b, b)), currSeq.subSeq(b, b));
          numData[s][b-1] = getNt2num().get(currSeq.subSeq(b, b));
        }
      }
      s++;
    }
    PoorNTSequence ntSeq = wA.getConsensus();
    //System.out.println(ntSeq.getSequence());
    for(int b=1; b<=ntSeq.getLength(); b++)
    {
      consensusData[b-1] = getNt2num().get(ntSeq.subSeq(b, b));
    }
  }

  /**
   * calculates which neighboring amino acidList are within a distance 'd'
   * from each amino acid. Results are stored in the StructureData array
   *
   * @param d
   */
  public void calcAAneighbours(double dist)
  {
    for(int i=0; i<alignLength; i++)
    {
      int currSite = i;
      int pairedSiteCount = 0;
      for(int j=0; j<alignLength; j++)
      {
        if(distMat[currSite][j] < dist && java.lang.Math.abs(currSite - j) > 2)
        {
          structData[currSite][pairedSiteCount] = j;
        }
        pairedSiteCount++;
      }
      if(pairedSiteCount > 0)
        structData[currSite][0]=1;
    }
  }



  /**
   * calculates Mutual Information of alignment according to
   *
   *  H(i) = Sum(s_i) [P(i) log(P(i))]
   *  and
   *  P(i) is siteProb of amino acid s_i @ site i
   *
   * but also takes into account structural data.
   * any sites within a distance 'd' are considered to represent
   * the same site.
   *
   * @param d
   */
  public void calcMutualInfo(double d)
  {

    dist = d;
    // calculate the entropy as before
    calcEntropy();

    // now pool the sites according to the structure data
    for(int i=0; i<alignLength; i++)
    {
      if(structData[i][0] != 0)
      {

      }
    }

    haveEntropy = true;

  }

  /**
   * find the minimum epistasis value in the array
   */
  public void findEpistasisRange()
  {
    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
      for(int siteJ=siteI+1; siteJ<getAlignLength();siteJ++)
      {
        if(getEpiDataNum()[siteI][siteJ] > getMaxEpi())
          maxEpi = getEpiDataNum()[siteI][siteJ];
      }
    }

    minEpi = getMaxEpi();
    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
      for(int siteJ=siteI+1; siteJ<getAlignLength();siteJ++)
      {
        if(getEpiDataNum()[siteI][siteJ] < getMinEpi()
                && getEpiDataNum()[siteI][siteJ] > 0)
          minEpi = getEpiDataNum()[siteI][siteJ];
      }
    }
  }

  /**
   * subtract a value from the array.  useful for making plots
   * @param baseline
   */
  public void subtractEpiBaseline(int baseline)
  {
    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
      for(int siteJ=siteI+1; siteJ<getAlignLength();siteJ++)
      {
        epiDataNum[siteI][siteJ] -= baseline;
        epiDataNum[siteJ][siteI] -= baseline;
      }
      epiDataNum[siteI][siteI] -= getMinEpi() - baseline;
    }
  }

  /**
   * find site pairs which differ from the consensus sequence
   */
  public void findEpistaticSites()
  {
    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
//      EpistasisColumn epiCol = new EpistasisColumn();
      for(int siteJ=siteI+1; siteJ<getAlignLength();siteJ++)
      {
        int epi=0;
        int nonepi=0;
//        EpistasisSet currES = new EpistasisSet();
        for(int s=0; s<this.getNumSeqs(); s++)
        {
/*          System.out.println(numData[s][siteI] + ":"
                  + consensusData[siteI] + ":(" + consensusProb[siteI] + ")"
                  + numData[s][siteJ]+ ":"
                  + consensusData[siteJ] + ":(" + consensusProb[siteJ] + ")" 
                  + (numData[s][siteI] != consensusData[siteI]
                  && numData[s][siteJ] != consensusData[siteJ])
                  );
 * 
 */
          if(numData[s][siteI] != getConsensusData()[siteI]
                  && numData[s][siteJ] != getConsensusData()[siteJ])
          {
            // we have epistatic sites
            epi++;
            getEpiDataNum()[siteI][siteJ]++;
            getEpiDataNum()[siteJ][siteI]++;

          }
          else
          {
            nonepi++;
          }
        }
        //System.out.println("epi:non" + epi + ":" + nonepi);

      }

    }
      System.out.println("done");

  }

  /**
   * This just prints out the counts, it doesn't export the individual mutation
   * data.
   *
   * @param bw
   * @throws java.io.IOException
   */
  public void writeEpistaticCount(BufferedWriter bw) throws java.io.IOException
  {
    for(int i=0; i<getSites().length; i++)
    {
      bw.write(getSites()[i] + "\t");
    }
    bw.write("\n");
    for(int i=0;i<alignLength;i++)
    {
      for(int j=0;j<alignLength;j++)
      {
        bw.write(getEpiDataNum()[i][j] + "\t");
      }
      bw.write("\n");
    }
/*    ListIterator itEpD = epiData.listIterator();
    while(itEpD.hasNext())
    {
      EpistasisColumn epiCol = (EpistasisColumn) itEpD.next();
      ListIterator itEpC = epiCol.getEpistasisCol().listIterator();
      while(itEpC.hasNext())
      {
        EpistasisSet currES = (EpistasisSet)itEpC.next();
        bw.write(currES.getEpistasisSet().size() + "\t");
      }
      bw.write("\n");
    }
 */

  }

  public void calcConsensus()
  {
    for(int site=0;  site<getAlignLength(); site++)
    {
      java.util.ArrayList<Integer> siteList = new java.util.ArrayList<Integer>();
      int atomCount[] = new int[getNumDistinctMols()];
      for(int a=0; a<this.acidList.size(); a++)
      {
          siteList.add(acidList.get(a));
      }
      int siteCount = 0;
      for(int seq=0; seq<getNumSeqs(); seq++)
      {
        int base = getNumData(seq, site);
        if(base >= atomCount.length)
          System.out.println("out of range " + base);
        else
         atomCount[base]++;
        /*
        int index = siteList.indexOf(base);
        if(siteList.indexOf(base) >= 0)
        {
          atomCount[index]++;
        }
        else
        {
          siteList.add(base);
          atomCount[siteCount]++;
          siteCount++;
        }
        */
      }
      int maxBaseCount = 0;
      int maxBase = -1;
      for(int s=0; s<siteList.size();s++)
      {
        if(atomCount[s]>maxBaseCount)
        {
          maxBaseCount = atomCount[s];
          maxBase = s;
        }
      }
      consensusData[site] = maxBase;
      consensusProb[site] = (double)maxBaseCount/(double)numSequences;
    }

    haveEntropy = true;
  }

  /**
   * Calculate entropy of alignment according to
   *
   *  H(i) = Sum(s_i) [P(i) log(P(i))]
   *  and
   *  P(i) is siteProb of amino acid s_i @ site i
   *
   */
  public void calcEntropy()
  {
    for(int site=0;  site<getAlignLength(); site++)
    {
      java.util.ArrayList<Integer> siteList = new java.util.ArrayList<Integer>();
      int atomCount[] = new int[getNumDistinctMols()];
      int siteCount = 0;
      double siteEntropy = 0.0;
      for(int seq=0; seq<getNumSeqs(); seq++)
      {
        int base = getNumData(seq, site);
        int index = siteList.indexOf(base);
        if(siteList.indexOf(base) >= 0)
        {
          atomCount[index]++;
        }
        else
        {
          siteList.add(base);
          atomCount[siteCount]++;
          siteCount++;
          //System.out.println("didnt find the molecule " + base);
        }
      }
      // now we have the mol count we can calc the entropy.
      for(int element: atomCount)
      {
        double siteProb = (double)(element)/(double)getNumSeqs();
        double logProb = 0.0;
        if(siteProb != 0.0)
        {
          logProb = Math.log(siteProb)/Math.log(2);
        }
        double i = -1.0d*siteProb*logProb;
        siteEntropy += i;
      }
      entropy[site] = siteEntropy;
    }

    haveEntropy = true;
  }

  /**
   *
   *  Calculate mutual information according to
   *  M(i,j) = H(i) + H(j) - H(i,j)
   *  where
   *  H(i,j) = Sum(s_i) Sum(s'_j) [P(i,j) log(P(i,j))]
   *  and
   *  P(i,j) is siteProb of s_i @ i and s'_j @ j
   *
   */
  public void calcMutualInfo()
  {
    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
      java.util.ArrayList<Integer> iSiteList
              = new java.util.ArrayList<Integer>();
      int iSiteCount = 0;
      for(int seq=0; seq<getNumSeqs(); seq++)
      {
        int index = iSiteList.indexOf(numData[seq][siteI]);
        if(index < 0)
        {
          iSiteList.add(numData[seq][siteI]);
          iSiteCount++;
        }
      }
      for(int siteJ=siteI; siteJ<getAlignLength();siteJ++)
      {
        // sum over all mols in the list.
          if(siteI != siteJ)
          {
          java.util.ArrayList<Integer> jSiteList
                  = new java.util.ArrayList<Integer>();
          int jSiteCount = 0;
          for(int seq=0; seq<getNumSeqs(); seq++)
          {
            int index = jSiteList.indexOf(numData[seq][siteJ]);
            if(index < 0)
            {
              jSiteList.add(numData[seq][siteJ]);
              jSiteCount++;
            }
          }
          double ijH = 0;
          java.util.Iterator iIter = iSiteList.listIterator();
          while(iIter.hasNext())
          {
            int currIMol = (Integer) iIter.next();
            java.util.Iterator jIter = jSiteList.listIterator();
            while(jIter.hasNext())
            {
              int s = 0;
              int count = 0;
              int currJMol = (Integer) jIter.next();
              while(s < getNumSeqs())
              {
                if((numData[s][siteI] == currIMol)
                  && (numData[s][siteJ] == currJMol))
                  count++;
                s++;
              }
              double ijdProb = (double)count
                      /(double)getNumSeqs();
              double dH = 0;
              if(count > 0)
              {
                double ijdLogProb = Math.log(ijdProb)/Math.log(2);
                dH = -1.0*ijdProb*ijdLogProb;
              }
              ijH += dH;
            }
          }
          sharedEntropy[siteI][siteJ] = ijH;
          mutualInfo[siteI][siteJ]
                  = getEntropy(siteI) + getEntropy(siteJ) - ijH;
          mutualInfo[siteJ][siteI] = mutualInfo[siteI][siteJ];
        }
      }
    }
  }

  /**
   * calculates Average Position Correction defined by Dunn S.D., Wahl, L.M. &
   * Gloor G.B. which estimates the background in the Mutual Information
   * caused by Entropy.
   *
   * APC = M(a, Avg(allSites))M(b, Avg(allSites)) / Avg(MI)
   *
   */
  public void calcAveragePositionCorrection()
  {
    /* First calculate average MI for each site */
    averageMI = 0;
    for(int siteI=0; siteI<getAlignLength(); siteI++)
    {
      averageColMI[siteI] = 0;
      for(int siteJ=0; siteJ<getAlignLength(); siteJ++)
      {
        averageColMI[siteI] += mutualInfo[siteI][siteJ];
        averageMI += mutualInfo[siteI][siteJ];
      }
      averageColMI[siteI] /= (float)(getAlignLength()-1);
    }
    averageMI /= (float)((getAlignLength()*(getAlignLength()-1)));

    for(int siteI=0; siteI<getAlignLength();siteI++)
    {
      for(int siteJ=siteI; siteJ<getAlignLength();siteJ++)
      {
        if(siteI == siteJ)
        {
          averagePosCorr[siteI][siteJ] = 0;
        }
        else
        {
          averagePosCorr[siteI][siteJ]
                  = averageColMI[siteI]*averageColMI[siteJ] / averageMI;
          averagePosCorr[siteJ][siteI] = averagePosCorr[siteI][siteJ];
        }
      }
    }

  }

  /**
   *
   *  Calculate mutual information according to
   *  M(i,j) = H(i) + H(j) - H(i,j)
   *  where
   *  H(i,j) = Sum(s_i) Sum(s'_j) [P(i,j) log(P(i,j))]
   *  and
   *  P(i,j) is siteProb of s_i @ i and s'_j @ j
   *
   */
  public void calcStructMutualInfo(double dist)
  {
    calcMutualInfo();
    calcAAneighbours(dist);
   // now we use the structure data and pool the sites
   for(int s=0; s<alignLength; s++)
   {
     if(structData[s][0] != 0)
     {
       int t=1;
       double sum=0.0;
       while(structData[s][t] != 0)
       {
         sum += mutualInfo[s][structData[s][t]];
         t++;
       }
       t=0;
       while(structData[s][t] != 0)
       {
         mutualInfo[s][structData[s][t]] = sum;
         t++;
       }
     }
   }
   haveMutualInfo = true;
  }

  /**
   * 
   * randomizes the alignment by shuffling the bases within each column
   * @param noOfShuffles - how many times to shuffle each column
   * 
   * @return  int[][] the shuffled alignment
   */
  public NumericalAlignment bootstrapAlignment(int noOfShuffles)
  {
    int bootAlign[][] = new int [numSequences][alignLength];
    
    bootAlign = numData;
    
    for(int c=0; c<alignLength; c++)
    {
      int r=0;
      while(r<noOfShuffles)
      {
        java.util.Random generator = new java.util.Random();
        int r1 = generator.nextInt(numSequences);
        int r2 = generator.nextInt(numSequences);
        if(r1 == r2) break;
        int temp = bootAlign[r1][c];
        bootAlign[r1][c] = bootAlign[r2][c];
        bootAlign[r2][c] = temp;
        r++;
      }
    }
    NumericalAlignment nA = new NumericalAlignment(
            bootAlign, getSeqNames(), getSites(),
            getNt2num());

    return nA;
  }

  /**
   * processes the specified bootstrap file to calculate probabilities
   *
   * @param dis
   */
  public void processBootstraps(DataInputStream dis)
  {
    int noOfBoots = 0;

    try{
      // we make no assumption the bootstrap file is complete
      // if the program crashed we might have incomplete data
      double bmut[][] ;
      int ii = 0;
      int jj = 0;
      bmut = new double[alignLength][alignLength];
      try {
        while (true)
        {
          // read next bootstrap array and find the max value
          double bmutMax=0.0;
          ii=0;
          while(ii<alignLength)
          {
            jj=0;
            while(jj<alignLength)
            {
             bmut[ii][jj] = dis.readDouble();
             if(bmut[ii][jj]>bmutMax)
             {
               bmutMax=bmut[ii][jj];
             }
             jj++;
            }
            ii++;
          }

          // check each site to see whether it exceeds bmutMax
          int iSite = 0;
          while(iSite<alignLength)
          {
            int jSite = 0;
            while (jSite < alignLength)
            {
              if(bmutMax>mutualInfo[iSite][jSite])
              {
                prob[iSite][jSite]++;
              }
              jSite++;
            }
            iSite++;
          }

          System.out.print(noOfBoots + "-");
          if(noOfBoots % 20 == 0)
          {
            System.out.print("\n");
          }
          noOfBoots++;
        }
      }
      catch (java.io.EOFException e)
      {
        if(ii*jj < alignLength*alignLength)
          System.err.println("caught unexpected end of file.  " + "\n" +
                  "likely due to incomplete bootstrap data due to aborted run."
                  + "\n" +
                  "If the number of processed bootstraps is correct. "
                  + "(" + noOfBoots + ")" +
                  "You can probably ignore this error");
        //e.printStackTrace();
      }
      dis.close();

      // finally, divide by the number of bootstraps to get the probability
      int iSite = 0;
      while(iSite<alignLength)
      {
        int jSite = 0;
        while (jSite < alignLength)
        {
          prob[iSite][jSite] = prob[iSite][jSite] / noOfBoots;
          jSite++;
        }
        iSite++;
      }
      iSite++;
    }
    catch(java.io.IOException ex)
    {
      System.out.print("bombed reading bootstrap file\n");
      System.err.print(ex + "\n");
    }


  }


  /**
   * Reads a NumericalAlignment in NEXUS format.
   * @param bw
   */
  public void readAlign(BufferedReader br)
          throws java.io.IOException,
          java.lang.ArrayIndexOutOfBoundsException
  {
    String headerLine = br.readLine();
    String params[] = headerLine.split("\\s+");
    setAlignLength(Integer.parseInt(params[0]));
    setNumSeqs(Integer.parseInt(params[1]));
    setMaxStructs(Integer.parseInt(params[2]));

    if (getNumData() == null)
    {
      initAlign();
    }
    String descriptor = br.readLine();
    while(descriptor != null)
    {
      // get next descriptor line
      if(descriptor.equalsIgnoreCase("[seqnames]"))
      {
        readSeqNames(br);
      }

      if(descriptor.equalsIgnoreCase("[alignment]"))
      {
        readAlignData(br);
        findDistinctMols();
      }

      if(descriptor.equalsIgnoreCase("[sites]"))
      {
        readSiteData(br);
      }

      if(descriptor.equalsIgnoreCase("[structure]"))
      {
        readStructureData(br);
      }
      descriptor = findNextDescriptor(br);
    }

   }

  /**
   * Reads the structure data associated with the alignment
   * @param br
   */
  public void readStructureData(BufferedReader br)
          throws java.io.IOException,
          java.lang.ArrayIndexOutOfBoundsException
  {
    int site = 0;
    String currentLine = "";

    while(!(currentLine = br.readLine()).isEmpty())
    {
      String vals[] = currentLine.split("\\s+");
      int i=0;
      for(String val: vals)
      {
        setStructData(site, i, Integer.parseInt(val));
        i++;
      }
      site++;
    }

  }

  /**
   * Reads the Sites data associated with the alignment
   * @param br
   * @throws java.io.IOException
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void readSiteData(BufferedReader br)
          throws java.io.IOException,
          java.lang.ArrayIndexOutOfBoundsException
  {
    String currentLine = br.readLine();
    String vals[] = currentLine.split("\\s+");
    int i=0;
    for(String val: vals)
    {
      setSiteData(i, Integer.parseInt(val));
      i++;
    }

  }
  /**
   * Reads the numerical data array section of a numerical alignment file.
   *
   * @param br  BufferedReader to the input file
   * @throws java.io.IOException
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void readAlignData(BufferedReader br)
          throws java.io.IOException,
          java.lang.ArrayIndexOutOfBoundsException
  {
    if(getAlignLength() == 0 || getNumSeqs() == 0 || getNumData() == null)
    {
      System.err.print("Data array not initialized");
      throw new java.lang.ArrayIndexOutOfBoundsException();
    }

    int line=0;
    while(line<getNumSeqs())
    {
      String currentLine = br.readLine();
      String vals[] = currentLine.split("\\s+");
      int i=0;
      for(String val: vals)
      {
        setNumData(i, line, Integer.parseInt(val));
        i++;
      }
      line++;
    }
    findDistinctMols();
  }

  /**
   * Writes entropy data with site position in ASCII format
   * @param bw
   */
  public void writeEntropy(BufferedWriter bw)
  {
    for(int i=0; i<entropy.length; i++)
    {
      try{
        bw.write(getSites()[i] + "\t" + entropy[i] + "\n");
      }
      catch(java.io.IOException exIO)
      {
        System.err.println("error writing entropy data");
        exIO.printStackTrace();
      }
    }
  }

  /**
   * writes mutual information in ASCII format
   * @param bw  BufferedWriter
   */
  public void writeMutualInfo( DataOutputStream dos)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          dos.writeDouble(mutualInfo[i][j]);
        }
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing mutual information");
      exIO.printStackTrace();
    }
  }

  /**
   * read Entropy from ASCII file
   * @param br BufferedReader
   */
  public void readEntropy(BufferedReader br)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        String nums[] = br.readLine().split("\\t");
        entropy[i] = Double.valueOf(nums[1]);
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error reading entropy data");
      exIO.printStackTrace();
    }
  }

  /**
   * read mutual information from ASCII file
   * @param br BufferedReader
   */
  public void readMutualInfo(BufferedReader br)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        String nums[] = br.readLine().split("\\t");
        for(int j=0; j<getAlignLength(); j++)
        {
          mutualInfo[i][j] = Double.valueOf(nums[j]);
        }
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error reading mutual information data");
      exIO.printStackTrace();
    }
  }

  /**
   * read shared entropy in ASCII format
   * @param bw  BufferedWriter
   */
  public void readSharedEntropy(BufferedReader br)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        String nums[] = br.readLine().split("\\t");
        for(int j=0; j<getAlignLength(); j++)
        {
          sharedEntropy[i][j] = Double.valueOf(nums[j]);
        }
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error reading mutual information data");
      exIO.printStackTrace();
    }
  }

  /**
   * writes mutual information in ASCII format
   * @param bw  BufferedWriter
   */
  public void writeMutualInfo(BufferedWriter bw)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          bw.write(mutualInfo[i][j] + "\t" );
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing mutual information data");
      exIO.printStackTrace();
    }
  }

  /**
   * writes shared entropy H(i,j) in ASCII format
   * @param bw  BufferedWriter
   */

  public void writeSharedEntropy(BufferedWriter bw)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          bw.write(sharedEntropy[i][j] + "\t" );
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing mutual information data");
      exIO.printStackTrace();
    }
  }

  /**
   * report all sites above the specified threshold
   * @param threshold:  Report
   * @param bw
   */
  public void reportMutualInfo(double miThreshold,
          double probThreshold, BufferedWriter bw)
  {
    int basesPerLine = 60;
    int numberOfSeqs = getSeqNames().length;
    int numberOfLines = numberOfSeqs/basesPerLine;
    DecimalFormat df = new DecimalFormat("0.00000");

    try{
      bw.write("MUTUAL INFORMATION REPORT:\n\n");
      for(int i=0; i<getAlignLength()-1; i++)
      {
        for(int j=i+1; j<getAlignLength(); j++)
        {
          java.util.ArrayList<String> pairList
                  = new java.util.ArrayList<String>();
          java.util.ArrayList<Integer> pairCount
                  = new ArrayList<Integer>();
//          if(mutualInfo[i][j] >= miThreshold && prob[i][j] <= probThreshold)
//          if(mutualInfo[i][j] >= miThreshold && prob[i][j] <= probThreshold)
          {
            double normalMI = 0;
            if(sharedEntropy[i][j] != 0)
              normalMI = mutualInfo[i][j]/sharedEntropy[i][j];
            bw.write("i\tj\tMI\tPROB\tAPC(i,h)\tNMI(i,j)\tH(i,j)\tH(i)\tH(j)\n");
            bw.write(i + "\t" + j + "\t"
                    + df.format(mutualInfo[i][j]) + "\t"
                    + df.format(prob[i][j]) + "\t"
                    + df.format(averagePosCorr[i][j]) + "\t"
                    + df.format(normalMI) + "\t"
                    + df.format(sharedEntropy[i][j]) + "\t"
                    + df.format(entropy[i]) + "\t"
                    + df.format(entropy[j])
                    + "\n\n");
            for(int l=0; l<numberOfLines; l++)
            {
              // print i
              // print matches
              // print j
              String colIstring = " " + (l*basesPerLine + 1);
              String matchString = "";
              String colJstring = "";
              for(int p=0; p<10; p++)
              {
                if(colIstring.length() < 10)
                  colIstring = colIstring.concat(" ");
                if(colJstring.length() < 10)
                  colJstring = colIstring.concat(" ");
                if(matchString.length() < 10)
                  matchString = colIstring.concat(" ");
              }
              for(int c=l*basesPerLine; c<(l+1)*basesPerLine; c++)
              {
                String baseI = getNum2nt().get(getNumData(c, i));
                String baseJ = getNum2nt().get(getNumData(c, j));
                String basePair = baseI + "-" + baseJ;
                if(pairList.contains(basePair))
                {
                  pairCount.set(pairList.indexOf(basePair),
                          pairCount.get(pairList.indexOf(basePair))+1);
                }
                else
                {
                  pairList.add(basePair);
                  pairCount.add(1);
                }
                colIstring = colIstring.concat(getNum2nt().get(getNumData(c, i)));
                colJstring = colJstring.concat(getNum2nt().get(getNumData(c, j)));
                if(getNumData(c, i) == getNumData(c, j))
                  matchString = matchString.concat(" ");
                else
                  matchString = matchString.concat("|");
              }
              bw.write("\t" + colIstring + "\n");
              bw.write("\t" + matchString + "\n");
              bw.write("\t" + colJstring + "\n");
              bw.write("\n\n");
            }
            // add the last fragment that doesn't fit onto a full line
            String colIstring = " " + (numberOfLines*basesPerLine + 1);
            String matchString = "";
            String colJstring = "";
            for(int p=0; p<10; p++)
            {
              if(colIstring.length() < 10)
                colIstring = colIstring.concat(" ");
              if(colJstring.length() < 10)
                colJstring = colIstring.concat(" ");
              if(matchString.length() < 10)
                matchString = colIstring.concat(" ");
            }
            for(int c=numberOfLines*basesPerLine; c<numberOfSeqs; c++)
            {
              String baseI = getNum2nt().get(getNumData(c, i));
              String baseJ = getNum2nt().get(getNumData(c, j));
              String basePair = baseI + "-" + baseJ;
              if(pairList.contains(basePair))
              {
                pairCount.set(pairList.indexOf(basePair),
                        pairCount.get(pairList.indexOf(basePair))+1);
              }
              else
              {
                pairList.add(basePair);
                pairCount.add(1);
              }
              colIstring = colIstring.concat(getNum2nt().get(getNumData(c, i)));
              colJstring = colJstring.concat(getNum2nt().get(getNumData(c, j)));
              if(getNumData(c, i) == getNumData(c, j))
                matchString = matchString.concat(" ");
              else
                matchString = matchString.concat("|");
            }
            bw.write("\t" + colIstring + "\n");
            bw.write("\t" + matchString + "\n");
            bw.write("\t" + colJstring + "\n");
            bw.write("\n");

            bw.write("bair pair distribution:\n");
            java.util.ListIterator itBPL = pairList.listIterator();
            java.util.ListIterator itBPC = pairCount.listIterator();
            bw.write("-------------------\n");
            bw.write("PAIR\tCOUNT\n");
            while(itBPL.hasNext())
            {
              bw.write(itBPL.next() + "\t" + itBPC.next() + "\n");
            }
            bw.write("\n\n" +
                    "-----------------------------------------------------" +
                    "\n\n\n" );
          }
        }
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing APC data");
      exIO.printStackTrace();
    }

  }

  /**
   * Object Serialization
   * @param dos DataOutputStream for the Serialization
   */
  public void write(DataOutputStream dos)
  {
    try{
      ObjectOutputStream osS =
            new ObjectOutputStream(new FileOutputStream("data.txt"));
        osS.writeObject(getSites());
        osS.writeObject(structData);
        osS.writeObject(dist);
        osS.writeObject(getSeqNames());
        osS.writeObject(distMat);
        osS.writeObject(getNt2num());
        osS.writeObject(getNum2nt());

        osS.writeObject( entropy);
        osS.writeObject(mutualInfo);
        osS.writeObject(sharedEntropy);
        osS.writeObject(averagePosCorr);
        osS.writeObject(averageColMI);
        osS.writeObject(averageEntropy);
        osS.writeObject(averageMI);
        osS.writeObject(prob);
        osS.writeObject(acidList);// stores a list of all AAs
        osS.writeObject(epiData);

        osS.writeBoolean(haveAlignment);
        osS.writeBoolean(haveAlignData);
        osS.writeBoolean(haveEntropy);
        osS.writeBoolean(haveMutualInfo);
        osS.writeBoolean(haveBootstraps);
        osS.writeBoolean(haveTree);
      osS.close();
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error serializing NumericalAlignment instance");
      exIO.printStackTrace();
    }
  }

  public void read(String fileName)
  {
    try{
    ObjectInputStream osS =
          new ObjectInputStream(new FileInputStream("data.txt"));
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error de-serializing NumericalAlignment instance");
      exIO.printStackTrace();
    }

  }
  /**
   * writes Average Position Correction Data in ASCII format
   * @param bw  BufferedWriter
   */

  public void writeAPC(BufferedWriter bw)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          bw.write(averagePosCorr[i][j] + "\t" );
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing APC data");
      exIO.printStackTrace();
    }
  }

  public void writeENormalizedMutualInfo(BufferedWriter bw)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          if(sharedEntropy[i][j] != 0)
            bw.write(mutualInfo[i][j]/sharedEntropy[i][j] + "\t" );
          else
            bw.write(0.0 + "\t");
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing Normalized MI data");
      exIO.printStackTrace();
    }
  }


  /**
   * write probability in ASCII format
   * @param bw  BufferedWriter
   */

  public void writeProbability(BufferedWriter bw)
  {
    try{
      for(int i=0; i<getAlignLength(); i++)
      {
        for(int j=0; j<getAlignLength(); j++)
        {
          bw.write(prob[i][j] + "\t" );
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error writing probability data");
      exIO.printStackTrace();
    }
  }

  /**
   * Scans the alignment to find how many distinct molecule types exist.
   * This can vary greatly from one alignment to another since particular
   * amino acidList may have been grouped by physical characteristics to form
   * a smaller set, or a single acid may have been subcategorized to generate
   * a series of subsets.
   */
  private void findDistinctMols()
  {
    for(int site=0; site<getAlignLength(); site++)
    {
      for(int seq=0; seq<getNumSeqs(); seq++)
      {
        //Integer mol = getNumData(seq, site);
        if(acidList.contains(getNumData(seq, site)) == false)
        {
          acidList.add(getNumData(seq, site));
        }
      }
    }
    numDistinctMols = acidList.size();
  }
  /**
   * Read in Sequence Names from supplied BufferedStream
   * @param br
   * @throws java.io.IOException
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void readSeqNames(BufferedReader br)
          throws java.io.IOException,
          java.lang.ArrayIndexOutOfBoundsException
  {
    int line=0;
    while(line<getNumSeqs())
    {
      setSeqName(line, br.readLine());
      line++;
    }

  }
  /**
   * searches for next descriptor line in the supplied BufferedStream
   * descriptor lines are of the form [DESCRIPTOR]
   * @param br BufferedReader Stream
   * @throws java.io.IOException
   * @return Descriptor Line
   */
  public String findNextDescriptor(BufferedReader br)
          throws java.io.IOException
  {
    String currentLine = "";

    while((currentLine = br.readLine()) != null)
    {
      currentLine = currentLine.trim();
      if(currentLine.startsWith("["))
      {
        return currentLine;
      }
    }
    return null;
  }

  /**
   * initializes all the data arrays using the given values for
   * Alignment Length, Number of Sequences & Max Struct.
   *
   */
  public void initAlign()
  {
    numData = new int[getAlignLength()][getNumSeqs()];
    sites = new int[getAlignLength()];
    structData = new int[getAlignLength()][getMaxStructs()];
    seqNames = new String[getNumSeqs()];
    entropy = new double[getAlignLength()];
    mutualInfo = new double[getAlignLength()][getAlignLength()];

  }

  /**
   * Returns Numerical Data Array
   * @return Alignment in Numerical Format;
   */
  public int[][] getNumData()
  {
    return numData;
  }

  /**
   * Returns single element of
   * @param s site
   * @param j sequence number
   * @return element at position (s,j)
   */
  public int getNumData(int s, int n)
  {
    return numData[s][n];
  }

  /**
   * set element of numerical alignment
   * @param site
   * @param seqNum
   * @param value
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void setNumData(int site, int seqNum, int value)
          throws java.lang.ArrayIndexOutOfBoundsException
  {
      numData[site][seqNum] = value;
  }

  public void setStructData(int site, int seqNum, int value)
          throws java.lang.ArrayIndexOutOfBoundsException
  {
      structData[site][seqNum] = value;
  }

  /**
   * set element of site data
   * @param site
   * @param value
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void setSiteData(int site, int value)
          throws java.lang.ArrayIndexOutOfBoundsException
  {
    sites[site] = value;
  }

  /**
   * set element of sequence nane array
   * @param site
   * @param name
   * @throws java.lang.ArrayIndexOutOfBoundsException
   */
  public void setSeqName(int site, String name)
          throws java.lang.ArrayIndexOutOfBoundsException
  {
    seqNames[site] = name;
  }


  /**
   * returns number of sequences in the alignment
   * @return numbersequences
   */
  public int getNumSeqs()
  {
    return numSequences;
  }

  /**
   * returns number of distinct molecule types in the alignment
   * @return
   */
  public int getNumDistinctMols()
  {
    return numDistinctMols;
  }

  /**
   * Set number of sequences in the alignment
   * @param n Number of sequences
   */
  public void setNumSeqs(int n)
  {
    numSequences = n;
  }

  /**
   * get alignment length
   * @return  Number of sites in the alignment
   */
  public int getAlignLength()
  {
    return alignLength;
  }

  /**
   * return distance value used in the mutual information calculation
   * returns a negative value if MIstruct has never been calculated.
   * 
   * @return
   */
  public double getDist()
  {
    return dist;
  }

  /**
   * returns Entropy for a single site.
   * @param site location
   * @return
   */
  public double getEntropy(int site)
  {
    return entropy[site];
  }

  /**
   * set alignment length
   * @param l Numbe of sites in the alignment.
   */
  public void setAlignLength(int l)
  {
    alignLength = l;
  }

  /**
   * The maximum number of sites that are associated with any one site
   * anywhere in the alignment.  See class documentation for more information
   * @return
   */
  public int getMaxStructs()
  {
    return maxStructures;
  }

  /**
   * sets the maximum number of sites that are associated with any one site
   * anywhere in the alignment.  See class documentation for more information
   * @param m Maximum no of Structures.
   */
  public void setMaxStructs(int m)
  {
    maxStructures = m;
  }

  /**
   * accessor method for the whether an Alignment exists for this instance
   * @return a boolean value
   */
  public boolean HaveAlignment()
  {
    return haveAlignment;
  }

  /**
   * accessor method for the whether an Alignment Data exists for this instance
   * @return a boolean value
   */
  public boolean HaveAlignData()
  {
    return haveAlignData;
  }

  /**
   * accessor method for the whether the Entropy exists for this instance
   * @return a boolean value
   */
  public boolean HaveEntropy()
  {
    return haveEntropy;
  }

  /**
   * accessor method for the whether the MutualInformation exists for this instance
   * @return a boolean value
   */
  public boolean HaveMutualInfo()
  {
    return haveMutualInfo;
  }

  /**
   * accessor method for the whether bootstrap data exists for this instance
   * @return a boolean value
   */
  public boolean HaveBootstraps()
  {
    return haveBootstraps;
  }

  /**
   * accessor method for the whether a NEXUS tree exists for this instance
   * @return a boolean value
  */
  public boolean HaveTree()
  {
    return haveTree;
  }

  public static void main(String args[])
  {
    NumericalAlignment n = new NumericalAlignment();
    try{
      BufferedReader br = new BufferedReader(new FileReader(args[0]));
      n.readAlign(br);
      n.calcEntropy();
      BufferedWriter bwE
              = new BufferedWriter(new FileWriter(args[0] + ".entropy"));
        n.writeEntropy(bwE);
      bwE.close();

      n.calcMutualInfo();
      BufferedWriter bwM
              = new BufferedWriter(new FileWriter(args[0] + ".mutinf"));
        n.writeMutualInfo(bwM);
      bwM.close();
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("problem reading file");
      exIO.printStackTrace();
    }
  }

  /**
   * find the range of the conservation in the alignment
   *
   * @return Point2D ( min, max )
   */
  public java.awt.geom.Point2D.Double findProbRange()
  {
    double pMax = 0;
    for (int p=0; p< this.consensusProb.length; p++)
    {
      if(consensusProb[p] > pMax)
        pMax = consensusProb[p];
    }
    double pMin = pMax;
    for (int p=0; p< this.consensusProb.length; p++)
    {
      if(consensusProb[p] < pMin)
        pMin = consensusProb[p];
    }
    return new java.awt.geom.Point2D.Double(pMin, pMax);
  }

  /**
   * @return the minEpi
   */
  public int getMinEpi() {
    return minEpi;
  }

  /**
   * @return the maxEpi
   */
  public int getMaxEpi() {
    return maxEpi;
  }

  /**
   * @return the epiDataNum
   */
  public int[][] getEpiDataNum() {
    return epiDataNum;
  }

  /**
   * @return the sites
   */
  public int[] getSites() {
    return sites;
  }

  /**
   * @return the consensusProb
   */
  public double[] getConsensusProb() {
    return consensusProb;
  }

  /**
   * @return the consensusData
   */
  public int[] getConsensusData() {
    return consensusData;
  }

  /**
   * @return the num2nt
   */
  public Hashtable<Integer, String> getNum2nt() {
    return num2nt;
  }

  /**
   * @return the seqNames
   */
  public String[] getSeqNames() {
    return seqNames;
  }

  /**
   * @return the nt2num
   */
  public Hashtable<String, Integer> getNt2num() {
    return nt2num;
  }
}
