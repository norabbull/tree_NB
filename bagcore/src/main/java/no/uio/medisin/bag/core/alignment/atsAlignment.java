/*
 * This is an extension of the wivAlignment class
 * It associates additional data that are associated with the entire
 * alignment (such as Structural info) and with specific sites (such as mutual
 * information and entropy)
 */

package no.uio.medisin.bag.core.alignment;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import no.uio.medisin.bag.core.aminoacid.aminoAcidListx;
import no.uio.medisin.bag.core.aminoacid.entropy;
import no.uio.medisin.bag.core.tree.MutualInfo;
import no.uio.medisin.bag.core.tree.SiteData;

/*import no.uio.medisin.bag.core.aminoacid.aminoAcidListx;
import no.uio.medisin.bag.core.alignment.WivAAAlignment;
import no.uio.medisin.bag.core.tree.SiteData;
import no.uio.medisin.bag.core.alignment.PoorAASequence;
import no.uio.medisin.bag.core.tree.MutualInfo;
import java.io.*;*/

/**
 *
 * @author sr
 */
public class atsAlignment extends WivAAAlignment
        implements java.io.Serializable {

  private boolean haveAlignment;
  private boolean haveAlignData;
  private boolean haveEntropy;
  private boolean haveMutualInfo;
  private boolean haveBootstraps;
  private boolean haveTree;

  private java.util.List<SiteData> siteInfo;

  public atsAlignment()
  {
    siteInfo = new java.util.ArrayList<SiteData>();
  }

  public atsAlignment(java.util.List<String> names)
  {
    super(names);
    siteInfo = new java.util.ArrayList<SiteData>();
  }

  public void reportAlignment(BufferedWriter bw)
  {
    // write results to output file alignFilename.results
    try
    {
      int base = 0;
      java.util.Iterator itSite = siteInfo.listIterator();
      while(itSite.hasNext())
      {
        aminoAcidListx currSite = (aminoAcidListx) itSite.next();
        bw.write("----------------------------------------\n");
        //System.out.print("----------------------------------------\n");
        bw.write("site " + base + " has "
                + currSite.getAAList().length());
        if(currSite.getAAList().length() == 1)
        {
          bw.write(" entry\n");
        }
        else
        {
          bw.write(" entries\n");
        }
        currSite.reportList(bw);
        bw.write("\n");
        base++;
      }
    }
    catch(java.io.IOException ex)
    {
      System.out.print("error writing out alignment analysis results\n");
      System.err.print(ex);
    }
  }

  public void calcEntropyAlign()
  {
    // write results to output file alignFilename.results
    int base = 0;
    java.util.Iterator itSite = siteInfo.listIterator();
    while(itSite.hasNext())
    {
      SiteData currData = (SiteData) itSite.next();
      aminoAcidListx currAAL = currData.getAAL();
      currData.getEntropy().setEntropy(currAAL.calcEntropy());
      base++;
    }
  }

  public void reportEntropyAlign(BufferedWriter bw)
  {
    // write results to output file alignFilename.results
    try
    {
      int base = 0;
      java.util.Iterator itSite = siteInfo.listIterator();
      while(itSite.hasNext())
      {
        SiteData currData = (SiteData) itSite.next();
        entropy e = currData.getEntropy();
        bw.write(base + "\t" + e.getEntropy() + "\n");
        base++;
      }
    }
    catch(java.io.IOException ex)
    {
      System.out.print("error writing out entropy results\n");
      System.err.print(ex);
    }
  }
  public void calcMutualInfo(double floor)
  {
    /*
     *  loop through all sites
     *  we calculate mutual information according to
     *  M(i,j) = H(i) + H(j) - H(i,j)
     *  where
     *  H(i,j) = Sum(s_i) Sum(s'_j) [P(i,j) log(P(i,j))]
     *  and
     *  P(i,j) is prob of s_i @ i and s'_j @ j
     *
     *  already have H(i)/H(j) so just need to calculate H(i,j) and M(i,j)
     */
    //printout("calculate mutual information..");
    java.util.List seqs = getSeq();
    int iBase = 0;
    java.util.Iterator itSite = siteInfo.listIterator();
    while(itSite.hasNext())
    {
      SiteData iCurrData = (SiteData) itSite.next();
      aminoAcidListx iCurrAAL = iCurrData.getAAL();
      double iH = iCurrData.getEntropy().getEntropy();

      // get alignment data for site i
      String iSlice[] = new String[getNoOfSeqs()];
      java.util.Iterator iSeqIter = seqs.listIterator();
      int iA = 0;

      while(iSeqIter.hasNext())
      {
        //String currSeq = (String) iSeqIter.next();
        PoorAASequence currSeq = (PoorAASequence) iSeqIter.next();
        iSlice[iA] = currSeq.subSeq(iBase+1, iBase+1);
        iA++;
      }

      // Site j
      int jBase = 0;
      java.util.Iterator jtSite = siteInfo.listIterator();
      while(jtSite.hasNext())
      {
        SiteData jCurrData = (SiteData) jtSite.next();
        if(iBase < jBase)
        {
          double jH = jCurrData.getEntropy().getEntropy();
          // get alignment data for site j
          String jSlice[] = new String[getNoOfSeqs()];
          java.util.Iterator jSeqIter = seqs.listIterator();
          int jA = 0;
          while(jSeqIter.hasNext())
          {
            PoorAASequence currSeq = (PoorAASequence) jSeqIter.next();
            jSlice[jA] = currSeq.subSeq(jBase+1, jBase+1);
            jA++;
          }

          aminoAcidListx jCurrAAL = jCurrData.getAAL();

          // now we have alignment slice for site i and j
          // and lists of which amino acids are present at i and j
          // so we can calc P(i,j) and H(i,j)
          float e = 0;
          int iE = 0;
          double ijH = 0;
          for(int ai=0; ai<iCurrAAL.getAAList().length(); ai++)
          {
            int jE = 0;
            String iCurrAA = iCurrAAL.getStringSymbolAt(ai+1);
            for(int aj=0; aj<jCurrAAL.getAAList().length(); aj++)
            {
              String jCurrAA = jCurrAAL.getStringSymbolAt(aj+1);
              int k=0;
              int count=0;
              while(k<getNoOfSeqs())
              {
                if(iCurrAA.equals(iSlice[k]) && jCurrAA.equals(jSlice[k]))
                {
                  count++;
                }
                k++;
              }

              double ijdProb = (double)count
                      /(double)getNoOfSeqs();
              double dH = 0;
              if(count > 0)
              {
                double ijdLogProb = Math.log(ijdProb)/Math.log(2);
                dH = -1.0*ijdProb*ijdLogProb;
              }
              else
              {
                dH = 0;
              }
              ijH += dH;
              jE++;
            }
            iE++;
          }
          if(iH + jH - ijH < floor)
          {
            iCurrData.getMutualInfo().setMInfo(jBase, 0.0);
          }
          else
          {
            iCurrData.getMutualInfo().setMInfo(jBase, iH + jH - ijH);
          }
          iCurrData.getMutualInfo().setEntropyIJ(jBase, ijH);
        }
        jBase++;
      }
      iBase++;
    }
//    printout("done\n");
  }

  public void reportMutualInfoAsc(BufferedWriter bw)
  {
    // write results to output file alignFilename.results
    try
    {
      java.util.Iterator itSite = siteInfo.listIterator();
      while(itSite.hasNext())
      {
        SiteData currData = (SiteData) itSite.next();
        MutualInfo mutInf = currData.getMutualInfo();
        int j=0;
        while(j<getAlignmentLength())
        {
           bw.write(Double.toString(mutInf.getMInfo(j)) + "\t");
          j++;
        }
        bw.write("\n");
      }
    }
    catch(java.io.IOException ex)
    {
      System.out.print("error writing out alignment analysis results\n");
      System.err.print(ex);
    }
  }

  public void reportMutualInfoBin(DataOutputStream dos)
  {
    // write results to output file alignFilename.results
    try
    {
      java.util.Iterator itSite = siteInfo.listIterator();
      while(itSite.hasNext())
      {
        SiteData currData = (SiteData) itSite.next();
        MutualInfo mutInf = currData.getMutualInfo();
        int j=0;
        while(j<getAlignmentLength())
        {
            dos.writeDouble(mutInf.getMInfo(j++));
        }
      }
    }
    catch(java.io.IOException ex)
    {
      System.out.print("error writing out mutual info data\n");
      System.err.print(ex);
    }
  }

  public void scanAlignment()
  {
    java.util.List seqs = getSeq();
    int base = 1;
    int ten = getAlignmentLength()/10;
    int hundred = (getAlignmentLength()/10)*10;

    while(base <= getAlignmentLength())
    {
      if((base+1)%ten == 0)
      {
        //printout( ((base+1-ten)*10)/ten + "%..");
      }
      //aminoAcidListx currentSite = new aminoAcidListx(Integer.toString(base));
      SiteData currentSiteInfo
              = new SiteData(getAlignmentLength(), Integer.toString(base));
      currentSiteInfo.setSitePosition(base);

      java.util.Iterator itSeqs = seqs.listIterator();
      while(itSeqs.hasNext())
      {
        PoorAASequence currSeq = (PoorAASequence) itSeqs.next();
        // now we check whether we already have the amino acid in the list
        if(currentSiteInfo.getAAL().hasAA(currSeq.symbolAt(base)) == false)
        {
          currentSiteInfo.getAAL().addAA(currSeq.symbolAt(base));
        }
        else
        {
          currentSiteInfo.getAAL().inc(currSeq.symbolAt(base));
        }
      }
      siteInfo.add(currentSiteInfo);
      base++;
    }
    haveAlignData = true;
  }

  public void reportSiteInfo(BufferedWriter bw)
  {
    // write results to output file  alignFilename.results
    try
    {
      int base = 0;
      java.util.Iterator itSite = siteInfo.listIterator();
      while(itSite.hasNext())
      {
        SiteData currData = (SiteData) itSite.next();
        aminoAcidListx currSite = currData.getAAL();
        bw.write("----------------------------------------\n");
        //System.out.print("----------------------------------------\n");
        bw.write("site " + base + " has "
                + currSite.getAAList().length());
        if(currSite.getAAList().length() == 1)
        {
          bw.write(" entry\n");
        }
        else
        {
          bw.write(" entries\n");
        }
        currSite.reportList(bw);
        bw.write("\n");
        base++;
      }
      bw.close();
    }
    catch(java.io.IOException ex)
    {
      System.out.print("error writing out alignment analysis results\n");
      System.err.print(ex);
    }
  }

  /**
   * writes the alignment in NumericalAlignment format
   * @param bw
   */
  public void exportAsNumericalAlignment(BufferedWriter bw)
  {

  }
  public void write(DataOutputStream dos)
  {
    try{
      ObjectOutputStream osS =
            new ObjectOutputStream(new FileOutputStream("data.txt"));
      osS.writeBoolean(haveAlignment);
      osS.writeBoolean(haveAlignData);
      osS.writeBoolean(haveEntropy);
      osS.writeBoolean(haveMutualInfo);
      osS.writeBoolean(haveBootstraps);
      osS.writeBoolean(haveTree);
      osS.writeObject(siteInfo);
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error serializing atsAlignment");
      exIO.printStackTrace();
    }
  }
  /***************************************************************************
   *                            accessor methods
   ***************************************************************************/

  public java.util.List<SiteData> getSiteInfo()
  {
    return siteInfo;
  }

  public java.util.ListIterator getSiteIter()
  {
    return siteInfo.listIterator();
  }
  /**
   * set whether an Alignment exists for this instance
   * @return void
   */
  public void setHaveAlignment(boolean b)
  {
    haveAlignment = b;
  }

  public void setHaveAlignData(boolean b)
  {
    haveAlignData = b;
  }

  /**
   * set whether the Entropy exists for this instance
   * @return void
   */
  public void setHaveEntropy(boolean b)
  {
    haveEntropy = b;
  }

  /**
   * set whether the MutualInformation exists for this instance
   * @return void
   */
  public void setHaveMutualInfo(boolean b)
  {
    haveMutualInfo = b;
  }

  /**
   * set whether bootstrap data exists for this instance
   * @return void
   */
  public void setHaveBootstraps(boolean b)
  {
    haveBootstraps = b;
  }

  /**
   * set whether a NEXUS tree exists for this instance
   * @return void
  */
  public void setHaveTree(boolean b)
  {
    haveTree = b;
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

}
