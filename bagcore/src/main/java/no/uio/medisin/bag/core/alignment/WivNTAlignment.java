
package no.uio.medisin.bag.core.alignment;

import no.uio.medisin.bag.core.alignment.PoorNTSequence;
import java.io.*;

import org.biojavax.bio.seq.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
/**
 * A basic class for storing an Amino Acid sequence alignment.
 * Sequences are stored as a RichSequence java.util.List
 * Alignments are created by loading a FASTA file
 *
 * @version 0.1
 * 2009-07-08 Replaced RichSequence with customized PoorNTSequence
 * to allow Serialization
 *
 *
 */
public class WivNTAlignment implements java.io.Serializable {

  java.util.ArrayList<PoorNTSequence> seqs;
  private int noOfSequences;
  private int alignmentLength;
  java.util.ArrayList<Integer> sites;
  private int longestName;
  PoorNTSequence consensus;
  private String alignmentFilename;

  /**
   *  Creates empty alignment
   */
  public WivNTAlignment()
  {
    seqs = new java.util.ArrayList<PoorNTSequence>();
    noOfSequences = 0;
    alignmentLength = 0;
    longestName = 0;
    sites = new java.util.ArrayList<Integer>();
  }

  /**
   * Initialize array with names and zero length sequences
   *
   * @param   java.util.ArrayList <String>  list of sequence names
   *
   */
  public WivNTAlignment(java.util.List<String> names)
  {
    seqs = new java.util.ArrayList<PoorNTSequence>();
    longestName = 0;
    java.util.ListIterator iN = names.listIterator();
    try{
      while (iN.hasNext())
      {
        String seqName = (String)iN.next();
        if(seqName.length() > longestName)
        {
          longestName = seqName.length();
        }
        seqs.add(new PoorNTSequence(seqName, ""));
      }
    }
    catch(java.lang.Exception exBioJ)
    {
      System.err.print("problem names only construction of alignment ");
      System.err.print(exBioJ);
    }

    noOfSequences = seqs.size();
    alignmentLength = seqs.get(0).getLength();
    sites = new java.util.ArrayList<Integer>();
    for(int s=1; s<=alignmentLength;s++)
    {
      sites.add(s);
    }
  }

  /**
   * Read FASTA alignment from file
   * @deprecated use @readFastaAlignment(BufferedReader br)
   * @param alignFile
   * @throws java.io.IOException
   */
  public void readFastaAlignment(File alignFile) throws java.io.IOException
  {
    BufferedReader br = new BufferedReader(new FileReader(alignFile));
      this.readFastaAlignment(br);
    br.close();
  }

  /**
   * read in a non-interleaved PHYLIP file
   * the file may contain multiple datasets and the BufferedReader is assumed
   * to be placed somewhere above the next data.
   *
   * @param br
   * @throws java.io.IOException
   */
  public void readPhylipAlignment(BufferedReader br) throws java.io.IOException
  {
    String line = "";
    String title = "";
    String seq = "";
    longestName = 0;
    alignmentLength = 0;
    noOfSequences = -1;
    seqs.clear();
    sites.clear();
    System.gc();
    
    // skip the blank lines and check for end of file
    title = br.readLine();
    while((title != null) && (title.equals("") == true))
    {
      title = br.readLine();
    }
    if(title==null)
      return;

    //title = br.readLine();
    String bits[] = title.split("\\s+");
    // if there is whitespace at the beginning of the line, need to skip the
    // first element.
    if(bits[0].length() == 0)
    {
      bits[0] = bits[1];
      bits[1] = bits[2];
    }
    // seqCount is only used to keep track of where we are in the file.
    int seqCount = Integer.valueOf(bits[0].trim());
    this.alignmentLength = Integer.valueOf(bits[1].trim());
    noOfSequences=0;  // if we are here, we have an alignment to read
    int noOfLines = 0;
    while(noOfLines < seqCount  && (line = br.readLine())!= null)
    {
      if(line.length() == 0)
        continue;
      String bits1[] = line.split("\\s");
      title = bits1[0].trim();
      //System.out.println(title);
      if(title.length() > longestName)
        longestName = title.length();
      for(int b=1; b<bits1.length;b++)
      {
        seq = seq.concat(bits1[b].trim());
        line = line.substring(1, line.length());
      }
      try
      {
        addSequence(title, seq);
      }
      catch(org.biojava.bio.BioException ex)
      {
        System.err.print("problem encountered adding sequence"
               + seq + "\n");
        System.err.print(ex);
      }
      if(this.getAlignmentLength() != 0 &&
             this.getAlignmentLength() != seq.length())
      {
        System.out.print("warning. " +
               "sequences appear to be different lengths\n");
      }
      this.setAlignmentLength(seq.length());
      title = line;
      seq = "";
      noOfLines++;
    }
    
    for(int s=1; s<=alignmentLength;s++)
    {
      sites.add(s);
    }
    this.findConsensus();
  }


  /**
   * Read FASTA alignment from file
   * @param alignFile
   * @throws java.io.IOException
   */
  public void readFastaAlignment(BufferedReader br) throws java.io.IOException
  {

    String line = "";
    String title = "";
    String seq = "";
    longestName = 0;
    alignmentLength = 0;
    seqs.clear();
    title = br.readLine();
    if(title.substring(0, 1).equals(">") == false)
    {
     throw new java.io.IOException("bad FASTA format on first line\n");
    }
    title = title.substring(1, title.length());
    while((line = br.readLine()) != null)
    {
      if(line.startsWith(">")==true)
      {
        line = line.substring(1, line.length());
        if(line.length() > longestName)
        {
          longestName = line.length();
        }
        try
        {
          addSequence(title, seq);
        }
        catch(org.biojava.bio.BioException ex)
        {
          System.err.print("problem encountered adding sequence"
                 + seq + "\n");
          System.err.print(ex);
        }
        if(this.getAlignmentLength() != 0 &&
               this.getAlignmentLength() != seq.length())
        {
          System.out.print("warning. " +
                 "sequences appear to be different lengths\n");
        }
        this.setAlignmentLength(seq.length());
        title = line;
        seq = "";
      }
      else
      {
        seq = seq.concat(line);
      }
    }
    try
    {
      this.addSequence(title, seq);
    }
    catch(org.biojava.bio.BioException ex)
    {
      System.err.print("problem encountered adding sequence"
             + seq + "\n");
      System.err.print(ex);
    }

    for(int s=1; s<=alignmentLength;s++)
    {
      sites.add(s);
    }
    this.findConsensus();
  }

  /**
   * write alignment in CLUSTAL format.
   *
   *
   * @param alignFile
   * @throws java.io.IOException
   *
   * writes 6 spaces after longest name
   * 50 sites / line
   */
  public void writeClustalAlignment(File alignFile) throws java.io.IOException
  {
    BufferedWriter bw
           = new BufferedWriter(new FileWriter(alignFile + ".aln"));
    bw.write("CLUSTAL 3.0.0 multiple sequence alignment FROM WIV\n\n\n");
    int numberOfBlocks = alignmentLength/50;
    int block = 0;
    while(block < numberOfBlocks)
    {
      java.util.ListIterator iSp = seqs.listIterator();
      while(iSp.hasNext())
      {
        PoorNTSequence ps = (PoorNTSequence) iSp.next();
        bw.write(ps.getName());
        for(int s=ps.getName().length(); s<this.getLongestName()+6; s++)
          bw.write(" ");
        bw.write(ps.subSeq(block*50+1, (block+1)*50));
        bw.write("\n");
      }
      bw.write("\n\n");
      block++;
    }

    java.util.ListIterator iSp = seqs.listIterator();
    while(iSp.hasNext())
    {
      PoorNTSequence ps = (PoorNTSequence) iSp.next();
      bw.write(ps.getName());
      for(int s=ps.getName().length(); s<this.getLongestName()+6; s++)
        bw.write(" ");
      bw.write(ps.subSeq(block*50+1, ps.getLength()));
      bw.write("\n");
    }
    bw.close();
  }


  /**
   * write alignment in PHYLIP format.
   * Not working properly...
   * 
   * @param alignFile
   * @throws java.io.IOException
   */
  public void writePhylipAlignment(File alignFile) throws java.io.IOException
   {
     // Needs to be finished
     BufferedWriter bw
             = new BufferedWriter(new FileWriter(alignFile + ".phylip"));
     java.util.ListIterator iSp = seqs.listIterator();

     int row = 0;
     while(row < alignmentLength%80 + 1)
     {
       while(iSp.hasNext())
       {
         PoorNTSequence ps = (PoorNTSequence) iSp.next();
         bw.write(">" + ps.getName());
         while(row < ps.getLength()%80 + 1)
         {
           bw.write(ps.subSeq(row*80, (row+1)*80));
           row++;
         }
       }
       row++;
     }
     row = 0;
     while(row < alignmentLength%80 + 1)
     {
       while(iSp.hasNext())
       {
         PoorNTSequence ps = (PoorNTSequence) iSp.next();
         bw.write(">" + ps.getName());
         while(row < ps.getLength()%80 + 1)
         {
           bw.write(ps.subSeq(row*80, (row+1)*80));
           row++;
         }
       }
       row++;
     }
     bw.close();
   }

  /**
   * remove sites that are invariant.
   * useful when calculating Mutual Information
   */
  public void removeInvariantSites()
  {
    // first of all mark the sites we want to keep
    java.util.ArrayList<Integer> variantSites = new java.util.ArrayList<Integer>();
    for(int b=1; b<=getAlignmentLength();b++)
    {
      java.util.ListIterator itPoor = getSeq().listIterator();
      String currBase = ((PoorNTSequence)getSeq().get(0)).subSeq(b, b);
      Boolean invariant = true;
      while(itPoor.hasNext())
      {
        PoorNTSequence currSeq = (PoorNTSequence) itPoor.next();
        if(currBase.equals(currSeq.subSeq(b, b)) == false)
        {
          //site not conserved.  drop it
          invariant = false;
          break;
        }
      }
      if(invariant == false)
      {
        variantSites.add(b);
      }
    }

    // now pull these sites out of the alignment
    java.util.ListIterator itS = getSeq().listIterator();
    int seqNo = 0;
    while(itS.hasNext())
    {
      String variantSeq = "";
      PoorNTSequence currSeq = (PoorNTSequence) itS.next();
      java.util.ListIterator itSites = variantSites.listIterator();
      while(itSites.hasNext())
      {
        variantSeq = variantSeq.concat(currSeq.baseAt((Integer)itSites.next()));
      }
      PoorNTSequence newSeq = new PoorNTSequence(currSeq.getName(), variantSeq);
      getSeq().set(seqNo, newSeq);
      seqNo++;
    }
    setAlignmentLength(this.getSeqAt(0).getLength());
    // update the consensus sequence
    java.util.ListIterator itSites = variantSites.listIterator();
    String newConsensus = "";
    //System.out.println(consensus.getSequence());
    while(itSites.hasNext())
    {
      int site = (Integer) itSites.next();
      newConsensus = newConsensus.concat(consensus.baseAt(site));
      //System.out.println(site + "," + consensus.baseAt(site));
    }
    //System.out.println(newConsensus);
    PoorNTSequence newSeq = new PoorNTSequence(consensus.getName(), newConsensus);
    consensus = newSeq;

    // finally replace the site info
    sites = variantSites;
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
      osS.writeInt(noOfSequences);
      osS.writeInt(alignmentLength);
      osS.writeInt(longestName);
      osS.writeObject(seqs);
    }
    catch(java.io.IOException exIO)
    {
      System.err.println("error serializing atsAlignment");
      exIO.printStackTrace();
    }
  }

  /**
   * Write out Consensus alignment for publication.
   *
   */
  /**
   *
   * @param start    start of region to write out
   * @param stop     end of region to write out
   * @param refSeq   name of sequence to use as Consensus reference
   * @param consChar which character to use to show consensus nucleotide
   * @param offset   Position to start numbering
   * @param numbers  show/hide nucleotide position
   * @param ruler    show/hide ruler
   *
   */
  public void writeConsensusAlignment(
          File alignFile,
          int start,
          int stop,
          String refName,
          String consChar,
          int offset,
          boolean numbers,
          boolean ruler) throws java.io.IOException
  {
    // first make sure the consensus sequence exists in the alignment
    int consSeqNum = -1;
    int s=0;
    java.util.ListIterator iS = seqs.listIterator();
    while(iS.hasNext())
    {
      PoorNTSequence ps = (PoorNTSequence) iS.next();
      if(ps.getName().equalsIgnoreCase(refName))
        consSeqNum = s;
      s++;
    }
    BufferedWriter bw
           = new BufferedWriter(new FileWriter(alignFile));
    if(consSeqNum < 0)
    {
      bw.write("couldn't find the consensus sequence in the alignment\n");
      bw.write("Reference sequence is: " + refName + "\n");
      iS = seqs.listIterator();
      bw.write("sequences in alignment are:\n");
      s = 0;
      while(iS.hasNext())
      {
        PoorNTSequence ps = (PoorNTSequence) iS.next();
        bw.write(s + "\t" + ps.getName() + "\n");
        s++;
      }
      bw.close();
      return;
    }

    // we found the reference sequence
    PoorNTSequence rS = seqs.get(consSeqNum);
    String refString = rS.subSeq(start, stop);

    // for now, can only write out the alignment as a single line

    // First write Number Line
    if(numbers)
    {
      /* get length of first number divisible by 10
       * s
       */
      if(offset == 0)
        offset = start;
      else
        offset = start + offset;
      
      int offset10 = offset;
      while(offset10%10 != 0)
      {
        offset10++;
      }

      // see whether we have enough space to fit the number string
      if(offset10 - offset <= String.valueOf(offset10).length())
        offset10 += 10;

      // calc how many numbers we need to display
      String numberLine = "";

      // add spacing for names
      for(int i=0; i<longestName+2;i++)
        numberLine = numberLine.concat(" ");

      //add first bit to the number line
      for(int t=start+offset-1; t<=offset10-Integer.toString(offset10).length();t++)
        numberLine = numberLine.concat(" ");
      numberLine = numberLine.concat(Integer.toString(offset10));

      int numOfNumbers = (int)(stop - start)/10 + 1;
      int n=1;
      while(n<numOfNumbers)
      {
        String numStr = Integer.toString(offset10 + n*10);
        int p=0;
        while(p + numStr.length() < 10)
        {
          numberLine = numberLine.concat(" ");
          p++;
        }
        numberLine = numberLine.concat(numStr);
        n++;
      }
      bw.write(numberLine + "\n");

      // write ruler
      String rulerString = "";
      int tenPos = numberLine.indexOf("0");

      // add spacing for names
      for(int i=0; i<longestName+2;i++)
        rulerString = rulerString.concat(" ");

      // add the first part of the ruler
      for(int i=offset; i<offset10; i++)
      {
        if(i%5 == 0)
          rulerString = rulerString.concat(":");
        else
          rulerString = rulerString.concat("-");
      }

      int t=offset10;
      while(t<stop+offset)  //subtract the bit we added at the front
      {
        if((t)%5 == 0)
          rulerString = rulerString.concat(":");
        else
          if((t-5)%10 == 0)
            rulerString = rulerString.concat("|");
          else
            rulerString = rulerString.concat("-");
        t++;
      }
      bw.write(rulerString + "\n");
    }

    iS = seqs.listIterator();
    while(iS.hasNext())
    {
      String outString = "";
      PoorNTSequence ps = (PoorNTSequence) iS.next();
      if(ps.getName().equals(refName))
      {
        outString = outString.concat(ps.getName());
        for(int i=0; i<longestName+2-ps.getName().length();i++)
        {
          outString = outString.concat(" ");
        }
        outString = outString.concat(refString);
        bw.write(outString + "\n");
      }
      else
      {
        outString = outString.concat(ps.getName());
        for(int i=0; i<longestName+2-ps.getName().length();i++)
        {
          outString = outString.concat(" ");
        }
        String seqString = ps.subSeq(start, stop);
        for(int t=0; t<seqString.length(); t++)
        {
          if(seqString.substring(t, t+1).equals(refString.substring(t,t+1)))
            outString = outString.concat(".");
          else
            outString = outString.concat(seqString.substring(t, t+1));
        }
        bw.write(outString + "\n");
      }
    }
    bw.close();
  }

  /**
   * Write Alignment in FASTA format
   * @param alignFile
   * @throws java.io.IOException
   */
  public void writeFastaAlignment(File alignFile) throws java.io.IOException
  {
   BufferedWriter bw
           = new BufferedWriter(new FileWriter(alignFile + ".fasta"));
   java.util.ListIterator iS = seqs.listIterator();

   while(iS.hasNext())
   {
     PoorNTSequence ps = (PoorNTSequence) iS.next();
     bw.write(">" + ps.getName() + "\n");
     int row = 0;
     while(row < alignmentLength/80)
     {
       bw.write(ps.subSeq(row*80+1, (row+1)*80) + "\n");
       row++;
     }
     bw.write(ps.subSeq(row*80+1, alignmentLength) + "\n");
   }
   bw.close();
  }

  /**
   * Method to allow extra sites to be added to an alignment
   * This may sound pointless, but useful for building sub-alignments
   * from regions of interest.
   * @param   site  Position in the alignment where the new sites will
   *          be inserted.
   * @param   rsList  List of the sites to be added.
   *          The sequences should be in the same order as the original
   *          alignment.  If they are, the routine will just kick up an
   *          exception.
   *
   * @return  an ArrayList of the sequence names
   */
   public void addSitesAt(int site, java.util.ArrayList<PoorNTSequence> colList)

   {
      try{
        if(getSeq().size() != 0 && colList.size() != getSeq().size())
        {
          System.err.print("the new sequences don't match the sequences" +
                "in the alignment - different sizes\n");
          throw new java.lang.Error();
        }
        java.util.ListIterator iAlign = getSeq().listIterator();
        java.util.ListIterator iCols = colList.listIterator();
        while(iCols.hasNext())
        {
          RichSequence currAlignSeq = (RichSequence) iAlign.next();
          RichSequence currColSeq = (RichSequence) iCols.next();
          if(currColSeq.getName()
                  .substring(0, currAlignSeq.getName().length())
                  .equals(currAlignSeq.getName()) == false)
          {
            System.err.print("the sequence names don't match\n");
            System.err.print("Aligment sequence name is\t"
                    + currAlignSeq.getName());
            System.err.print("Column sequence name is\t"
                    + currColSeq.getName());
            throw new java.lang.Error();
          }
          Edit e = new Edit(
                  currAlignSeq.length()+1,
                  0,
                  ProteinTools.createProtein(currColSeq.seqString()));
          currAlignSeq.edit(e);
          currAlignSeq.length();
        }

        setAlignmentLength(
                getAlignmentLength()
                + ((RichSequence) colList.get(0)).length()
                );
      }
      catch(java.lang.Error ex)
      {
        System.err.print("bombed trying to add sites to the alignment\n");
        System.err.print(ex);
      }
      catch(org.biojava.bio.symbol.IllegalSymbolException exSym)
      {
        System.err.print("bad symbol in add sites\n");
        System.err.print(exSym);
      }
      catch(org.biojava.bio.symbol.IllegalAlphabetException exSym)
      {
        System.err.print("bad symbol in add sites\n");
        System.err.print(exSym);
      }
   }
/**
   * Method to return a set of sites from the alignment
   *
   *
   * @param   start_site  First position in the alignment
   *
   * @param   end_site Last position in the alignment
   *          The sequences should be in the same order as the original
   *          alignment.  If they are, the routine will just kick up an
   *          exception.
   *
   * @return  a java.util.List<RichSequence>
   */
   public java.util.ArrayList<PoorNTSequence> getSites(int start_site, int end_site)
   {
     //java.util.ListIterator iSeqs = sequences.listIterator();
     java.util.ListIterator ipSeqs = seqs.listIterator();
     java.util.ArrayList<PoorNTSequence> colSeqs
             = new java.util.ArrayList<PoorNTSequence>();
     try{
       while(ipSeqs.hasNext())
       {
         PoorNTSequence currSeq = (PoorNTSequence) ipSeqs.next();
         colSeqs.add(new PoorNTSequence(
                 currSeq.getName(),
             currSeq.subSeq(start_site, end_site)
             ));
       }
     }
     catch(java.lang.Exception exBioJ)
     {
       System.err.print("error adding sub-sequence from alignment\n");
       System.err.print(exBioJ);
     }
     return colSeqs;
   }

   /**
    * returns a list of all the site positions in the alignment.
    * This is different from getSites, which returns a subalignment
    * of amino acids.  the Site positions are the positions of the amino acids
    * within some parent sequence.  The reason we do this is to allow
    * the alignment to represent a discontinous set of sites, such as variable
    * sites only within the original alignment.  This is useful when calculatin
    * Entropy or Mutual Information since in these cases the invariant sites
    * contain no useful information and confound efforts to analyze the results.
    *
    * @return List of sites in the alignment.
    * 
    */
   public java.util.ArrayList<Integer> getSitesList()
   {
     return sites;
   }

  /**
  * calculates consensus sequence from the aligmnent.
  *
  */
  public void findConsensus()
  {
    String cns = "";
    for(int b=1; b<=getAlignmentLength();b++)
    {
      java.util.ArrayList<Integer> count = new java.util.ArrayList<Integer>();
      java.util.ArrayList<String> AAlist = new java.util.ArrayList<String>();

      java.util.ListIterator itPoor = getSeq().listIterator();
      while(itPoor.hasNext())
      {
        PoorNTSequence currSeq = (PoorNTSequence) itPoor.next();
        String currBase = currSeq.subSeq(b, b);
        if(AAlist.contains(currBase))
        {
          count.set(AAlist.indexOf(currBase),
                  count.get(AAlist.indexOf(currBase)) + 1);

        }
        else
        {
          AAlist.add(currBase);
          count.add(1);
        }
      }
      java.util.ListIterator itC = count.listIterator();
      java.util.ListIterator itAA = AAlist.listIterator();
      int max = 0;
      String maxBase = "";
      while(itC.hasNext())
      {
        int currC = (Integer) itC.next();
        String base = (String) itAA.next();
        if(currC > max)
        {
          max = currC;
          maxBase = base;
        }
      }
      cns = cns.concat(maxBase);
    }
    consensus = new PoorNTSequence("CONSENSUS", cns);
  }

  /**
   * search for a particular sequence name in the alignment.
   * Returns row number if the alignment contains the sequence
   * -1 if the sequence name wasn't found
   * @param name
   * @return
   */
  public int findSequenceName(String name)
  {
    java.util.ListIterator itN = seqs.listIterator();
    int r=0;
    while(itN.hasNext())
    {
      PoorNTSequence pSeq = (PoorNTSequence)itN.next();
      if(pSeq.getName().equalsIgnoreCase(name))
        return r;
      else
        r++;
    }
    return -1;
  }

  
  /**
   * pulls out a subAlignment.  
   * If start is greater than stop it assumes a circularized genome and
   * ties sequence together across the join.
   * @param start
   * @param stop
   * @return
   */
  public WivNTAlignment getSubAlignment(int start, int stop)
  {
    WivNTAlignment subAlign = new WivNTAlignment();
    java.util.ListIterator itSeq = seqs.listIterator();
    while (itSeq.hasNext())
    {
      try{
        PoorNTSequence seq = (PoorNTSequence)itSeq.next();
        subAlign.addSequence(seq.getName() + "__" + start + "_" + stop,
                seq.subSeq(start, stop));
      }
      catch(org.biojava.bio.BioException exBJ)
      {
        System.err.println("error adding sequence to subAlignment");
        exBJ.printStackTrace();
      }
    }
    return subAlign;
  }

  /**
   * bootscan each sequence in the alignment with respect to the named reference
   * sequence
   * @param refName
   * @param offset
   * @param windowSize
   * @param step
   */
  public double[][] bootscanAlignment(
          String refName,
          int offset,
          int windowSize,
          int step
          )
    {
      double sim[][] = new double[seqs.size()+1][(alignmentLength-windowSize)/step];
      PoorNTSequence pa = (PoorNTSequence) seqs.get(findSequenceName(refName));
      java.util.ListIterator itSeq = seqs.listIterator();
      int s=0;
      while (itSeq.hasNext())
      {
        PoorNTSequence currPAASeq = (PoorNTSequence)itSeq.next();
        int site=0; int i=0;
        while(site + windowSize < getAlignmentLength())
        {
          int wSite = 1;
          int matches = 0;
          while(wSite <= windowSize)
          {
            if(currPAASeq.baseAt(site + wSite).equals(pa.baseAt(site + wSite)))
              matches++;
            wSite++;
          }
          double similarity = (double)matches/(double)windowSize;
          sim[s+1][i] = similarity;
          if(s==0)
            sim[0][i] = site + windowSize/2;
          site += step;
          i++;
        }
        s++;
      }
      return sim;
    }
  /**
   * doesn't do anything at the moment.
   * @return
   */
  public int getMaxStruct()
   {
     return 0;
   }
  /***************************************************************************
   *                            accessor methods
   ***************************************************************************/
  /**
   * accessor method for the names of the sequences in the alignment
   * @return an ArrayList of the sequence names
   */
public java.util.List getNames()
  {
    java.util.ListIterator itS = seqs.listIterator();
    java.util.List<String> names = new java.util.ArrayList <String> ();
    while(itS.hasNext())
    {
      names.add(((PoorNTSequence)itS.next()).getName());
    }

    return names;
  }

  /**
   * return consensus sequence for the alignment
   * @return
   */
  public PoorNTSequence getConsensus()
  {
    return consensus;
  }

  /**
   * returns a PoorNTSequence list of sequences in the alignment
   * @return
   */
  public java.util.List getSeq()
  {
    return seqs;
  }

  /**
   * returns the Sequence at specified row in the alignment.
   * @param row
   * @return
   */
  public PoorNTSequence getSeqAt(int row)
  {
    return (PoorNTSequence) seqs.get(row);
  }

  /**
   * returns the base/amino acid at specified row and column in the alignment
   * @param row
   * @param col
   * @return
   */
  public String getLetterAt(int row, int col)
  {
    return seqs.get(row).subSeq(col, col+1).toString();
  }

  /**
   * returns a String representation of the sequence at the specified row
   * and between the specified start and stop columns.
   * @param row
   * @param startCol
   * @param endCol
   * @return  String
   */
  public String getSeqAt(int row, int startCol, int endCol)
  {
    return seqs.get(row).subSeq(startCol, endCol).toString();
  }

  /**
   * returns number of sequences in the alignment
   * @return
   */
  public int getNoOfSeqs()
  {
    return noOfSequences;
  }

  /** 
   * Increment the sequence count.
   * Not sure why we would ever want to use this.
   */
  public void incNoOfSeqs()
  {
    noOfSequences++;
  }

  /**
   * Add a new sequence to the bottom of the alignment
   * @param name
   * @param seq
   * @throws org.biojava.bio.BioException
   */
  public void addSequence(String name, String seq) throws org.biojava.bio.BioException
  {
    seqs.add(new PoorNTSequence(name, seq));
    noOfSequences++;
    if(seqs.get(seqs.size()-1).getLength() > alignmentLength)
    {
      alignmentLength = seqs.get(seqs.size()-1).getLength();
    }
    if(seqs.get(seqs.size()-1).getName().length() > longestName)
    {
      longestName = seqs.get(seqs.size()-1).getName().length();
    }

  }

  /**
   * Add a new sequence to the bottom of the alignment
   * @param name
   * @param seq
   * @throws org.biojava.bio.BioException
   */
  public void addSequence(PoorNTSequence ps) throws org.biojava.bio.BioException
  {
    seqs.add(ps);
    noOfSequences++;
    if(seqs.get(seqs.size()-1).getLength() > alignmentLength)
    {
      alignmentLength = seqs.get(seqs.size()-1).getLength();
    }
    if(seqs.get(seqs.size()-1).getName().length() > longestName)
    {
      longestName = seqs.get(seqs.size()-1).getName().length();
    }

  }

  /**
   * set the length of the alignment.  There isn't any reason to use this method.
   * @param l
   */
  public void setAlignmentLength(int l)
  {
    alignmentLength = l;
  }

  public int getAlignmentLength()
  {
    return alignmentLength;
  }

  public int getLongestName()
  {
    return longestName;
  }

  public void setAlignFilename(String s)
  {
   alignmentFilename = s;
  }

  public String getAlignFilename()
  {
   return alignmentFilename;
  }


}

