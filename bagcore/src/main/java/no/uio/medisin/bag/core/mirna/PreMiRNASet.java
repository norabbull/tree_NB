/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import no.uio.medisin.bag.core.sequence.ChromosomeString;
import no.uio.medisin.bag.core.sequence.GenomeSeq;
import no.uio.medisin.bag.core.sequence.StrandString;

import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * Stores a set of pre-miRNAs.
 * This can be a miRBase release, or a subset of interesting pre-miRNAs
 * 
 * the class can load from miRBase GFF3, Fasta or from STR format (which contains the secondary structure information)
 * 
 * 
 * @author sr
 * 
 */
public class PreMiRNASet {
    
    static  Logger                      logger                          = LogManager.getLogger();

  	final static    String              INDENT1 = "--";
    final static    String              INDENT2 = "----";
  
    static  final   String              MIRBASEVERSION                  = "microRNAs:";
    static  final   String              GENOME_BUILD_ID                 = "genome-build-id:";
    static  final   String              GENOME_BUILD_ACC                = "genome-build-accession:";

    private int                         mirbaseVersion;
    private String                      buildID;
    private String                      buildAccession;
    private List<PreMiRNA>              preMiRNAList;
    private GenomeSeq                   genomeSeq;
    private HashMap <String, String>    MiRBase2GenomefeaturesMap ;
    private HashMap <String, String>    Genome2MiRBasefeaturesMap ;

    
    private String                      gffFilepath;
    private String                      matureFilepath;
    private String                      hairpinFilepath;
    private String                      structureFilepath;
    private String                      genomeFAFilepath;
    private String                      namesMapFilepath;
    private String                      organismCode;
    private String                      outputFolder;
    
    private Boolean                     singleSeqFile = false;

    private Boolean                     gffLoaded;
    private Boolean                     matureLoaded;
    private Boolean                     hairpinLoaded;
    private Boolean                     structureLoaded;
    
    private int                         noOfmiRNAs;

    public static void main(String args[]) throws IOException {
    	PreMiRNASet preMiRSet = new PreMiRNASet();
  		String gffFile = "/Users/simonray/Dropbox/data/mirbase/22.1/mmu.21pt1.gff3";
  		String strFile = "/Users/simonray/Dropbox/data/mirbase/22.1/miRNA.str";
      
  		String host = "mmu";
    	try {
    		logger.debug("reading from miRBase structure file <" + strFile + ">");
        preMiRSet.loadMiRBaseDataFromGFF(gffFile);
        preMiRSet.loadMiRBaseDataFromStructureData(host, strFile);   
      	  preMiRSet.trimBackHairpins();
      	  String trimmedStrFAFile = FilenameUtils.removeExtension(strFile) + "_" + host + "_trimmedStr.fa";
        preMiRSet.writePreMiRWithStructureInfosAsFA(trimmedStrFAFile);
    	}catch(IOException ex){
        logger.error("error reading miRBase gff reference file <" +  strFile + ">\n" + ex.toString());
        throw new IOException("error reading miRBase gff reference file <" +  strFile + ">");
		  }
    }
    
    public PreMiRNASet() {
    	preMiRNAList      = new ArrayList<>();
      genomeSeq         = new GenomeSeq();
      MiRBase2GenomefeaturesMap   = new HashMap<String, String>();
      Genome2MiRBasefeaturesMap   = new HashMap<String, String>();
     	gffLoaded         = false;
    	matureLoaded      = false;
    	hairpinLoaded     = false;
    	structureLoaded   = false;
    	
    	noOfmiRNAs = 0;
    }
    
    
    /**
     * load fasta sequence for miRNA entries. 
     * If the GFF file or Secondary Structure has been loaded, then match to existing entries.
     * at least one of these files must be loaded to get the pre-miRNA<->miRNA information
     * 
     * @param host
     * @param thisMatureFilePath
     * @param loadSequence
     * @throws IOException
     */
    public void loadMiRBaseDataFromMatureFA(String host, String thisMatureFilePath) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      	this.setMatureFilepath(thisMatureFilePath);
    	
    	if(this.getGffLoaded()==false && this.getSStructureLoaded()==false) {
    		logger.warn("a matching GFF/Secondary Structure file has not been loaded."
    				+ "You cannot load a mature file without matching pre-miRNA information");    		
    	}
    	int newEntries = 0;
    	this.setNoOfmiRNAs(0);
      try{
        String headerLine = null;
        String seqLine = null;
        
        logger.info("reading miRBase 'mature.fa' file <" +  thisMatureFilePath + ">");
        BufferedReader brMiR = new BufferedReader(new FileReader(new File(thisMatureFilePath)));
          headerLine = brMiR.readLine();            
          do{
            if(headerLine == null) break;
            seqLine = brMiR.readLine();
            if(headerLine.contains(host) && this.getGffLoaded()==true) {
            	
            	if(this.findMiRNAEntryByID(MiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine))>=0) {
            		MiRNA miR = this.getMiRNAEntryByID(MiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine));
            		miR.setSeq(seqLine.trim());
            		this.noOfmiRNAs++;
            	}else {
            		/* 
            		 * we shouldn't get here if we are parsing miRBase data as all entries should be present in all files
            		 */
            		logger.warn("didn't find miRNA <" + MiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine) + "> in"
            				+ " the miRNAs that were loaded from the GFF file <" + this.getGffFilepath() + ">");
            		newEntries++;
            	}

            }

          }while(true);
	      brMiR.close();
	      logger.info("read <" + getPreMiRNAList().size() + "> pre-miRNA entries");
	      logger.info("read <" + this.getNoOfmiRNAs() + "> miRNA entries");
	      logger.info("<" + newEntries + "> were not found in the existing list (i.e. were orphaned and not loaded)");
      
		  }
		  catch(IOException ex){
		      logger.error("error reading miRBase gff reference file <" +  thisMatureFilePath + ">\n" + ex.toString());
		      throw new IOException("error reading miRBase gff reference file <" +  thisMatureFilePath + ">");
		  }
      this.setMatureLoaded(true);
    	
    }

    /**
     * load fasta sequence for pre-miRNA entries. 
     * If the GFF file has been loaded, then match to existing entries. Unlike the mature.fa file, the GFF
     * file isn't required, because we can still create a pre-miRNA list from the hairpin sequences alone.
     * @param host
     * @param hairpinFilePath
     * @param loadSequence
     * @throws IOException
     */
    public void loadMiRBaseDataFromHairpinFA(String host, String hairpinFilePath) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
    	this.setHairpinFilepath(hairpinFilePath);
    	int newEntries = 0;
      try{
        String headerLine = null;
        String seqLine = null;
        
        logger.info("reading miRBase 'hairpin.fa' file <" +  hairpinFilePath + ">");
        BufferedReader brPre = new BufferedReader(new FileReader(new File(hairpinFilePath)));
          do{
            headerLine = brPre.readLine();            
            if(headerLine == null) 
              break;
            seqLine = brPre.readLine();
            if(headerLine.contains(host) && (this.getGffLoaded()==true || this.getSStructureLoaded())) {
            	
                if(this.findPreMiRNAEntryByID(PreMiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine))>=0) {
                  PreMiRNA preMiR = this.getPreMiRNAEntryByID(PreMiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine));
                  preMiR.setSeq(seqLine.trim());
                }else {
                  /* 
                   * we shouldn't get here if we are parsing miRBase data as all entries should be present in all files
                   */
                  logger.warn("didn't find pre-miRNA <" + MiRNA.parseIDFromMiRBaseFAHeaderLine(headerLine) + "> in"
                      + " the pre-miRNAs that were loaded from the GFF / SStructure file");
                  if(this.getGffLoaded())
                    logger.warn("GFF Filepath is <" + this.getGffFilepath() + ">");
                  if(this.getSStructureLoaded())
                    logger.warn("SStructure Filepath is <" + this.getStructureFilepath() + ">");
                  newEntries++;
                }

            }else{
              PreMiRNA preMiR = new PreMiRNA();
              preMiR.parseMiRBaseFAHeaderLine(headerLine);
              preMiR.setSeq(seqLine.trim()); 
              this.getPreMiRNAList().add(preMiR);
            }

          }while(true);
	      brPre.close();
	      logger.info("read <" + getPreMiRNAList().size() + "> pre-miRNA entries");
	      logger.info("<" + newEntries + "> were not found in the existing list");
      
		  }
		  catch(IOException ex){
		      logger.error("error reading miRBase hairpin.fa file <" +  hairpinFilePath + ">\n" + ex.toString());
		      throw new IOException("error reading miRBase hairpin.fa file <" +  hairpinFilePath + ">");
		  }
      this.setHairpinLoaded(true);
    	
    }

    /**
     * load secondary structures sequence for pre-miRNA miRNA entries. 
     * If the GFF file has been loaded, then match to existing entries.
     * There appears to be 8 lines / entry. For example:
     * 
     *  1.  >cel-let-7 (-42.90)   [cel-let-7-5p:17-38] [cel-let-7-3p:60-81]
     *  2.  
     *  3.  ------uaca    gga             U              ---  aaua 
     *  4.            cugu   uccggUGAGGUAG AGGUUGUAUAGUUu   gg    u
     *  5.            ||||   ||||||||||||| ||||||||||||||   ||     
     *  6.            gaca   aggCCAUUCCAUC UUUAACGUAUCaag   cc    u
     *  7.  agcuucucaa    --g             U              ugg  acca 
     *  8.
     *  
     * 
     * @param host
     * @param ssFilePath
     * @param loadSequence
     * @throws IOException
     */
    public void loadMiRBaseDataFromStructureData(String host, String ssFilePath) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      this.structureFilepath = ssFilePath;
      try{
      	int newEntries=0;
      	int oldEntries=0;
        String nextLine = null;
        
        ArrayList<String> lines = new ArrayList<>();
        logger.info("reading miRBase structure file <" +  ssFilePath + ">");
        BufferedReader brMiR = new BufferedReader(new FileReader(new File(ssFilePath)));
          nextLine = brMiR.readLine();            
          do{
            if(nextLine == null) break;
            if(nextLine.contains(">")) {
              do {
                  lines.add(nextLine);  
                  nextLine = brMiR.readLine();
                  if(nextLine == null || nextLine.contains(">")) break;
                  
              }while(true);
              /**
               * is the pre-miRNA already in the list?
               */
              String headerLine = lines.get(0);
              if(lines.get(0).contains(host)) {
//                        if(headerLine.contains("cfa-mir-371")){
//                          logger.debug("time for a break");
//                        }
              	if(this.findPreMiRNAEntryByName(PreMiRNA.parseNameFromMiRBaseSSHeaderLine(headerLine))>=0) {
              		PreMiRNA preMiR = this.getPreMiRNAEntryByName(PreMiRNA.parseNameFromMiRBaseSSHeaderLine(headerLine));
  	              preMiR.parseMiRBaseSStructureEntry(lines);
              		oldEntries++;
              	}else {
              		/* 
              		 * we shouldn't get here if we are parsing miRBase data as all entries should be present in all files
              		 */
//              		if(PreMiRNA.parseNameFromMiRBaseFAHeaderLine(headerLine).contentEquals("cel-mir-62"))
//              			logger.info("time to break");
              		logger.warn("didn't find pre-miRNA <" + PreMiRNA.parseNameFromMiRBaseSSHeaderLine(headerLine) + "> in"
              				+ " the pre-miRNAs that were loaded from the GFF / hairpin file");
              		if(this.getGffLoaded())
              			logger.warn("GFF Filepath is <" + this.getGffFilepath() + ">");
              		if(this.getSStructureLoaded())
              			logger.warn("SStructure Filepath is <" + this.getStructureFilepath() + ">");
  	              PreMiRNA preMiR = new PreMiRNA();
  	              preMiR.setName(PreMiRNA.parseNameFromMiRBaseFAHeaderLine(headerLine));
  	              preMiR.parseMiRBaseSStructureEntry(lines);
  	              //preMiR.trimBackHairpinStructure();
  	              this.getPreMiRNAList().add(preMiR);
              		newEntries++;
              	}
              	
              }
              lines.clear();
            }
          }while(true);
	      brMiR.close();
	      logger.info("read <" + (oldEntries + newEntries) + "> pre-miRNA entries");
	      logger.info("     <" + oldEntries + "> were already in the pre-miRNA list");
	      logger.info("     <" + newEntries + "> were new entries");
      
		  }
		  catch(IOException ex){
		      logger.error("error reading miRBase gff reference file <" +  ssFilePath + ">\n" + ex.toString());
		      throw new IOException("error reading miRBase gff reference file <" +  ssFilePath + ">");
		  }
      this.setStructureLoaded(true);
    	
    }
    
    /**
     * trim back the hairpin structure for this pre-miRNA
     */
    public void trimBackHairpins() {
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
    	for(PreMiRNA preMiR: this.getPreMiRNAList()) {
    		if(preMiR.hasStructureInfo()) {
    			preMiR.trimBackHairpinStructure();
    		}else {
    			logger.info("--no structure information for pre-miRNA [" + preMiR.getMiID() + "]");
    		}
    	}
    }
    
    

    /**
     * 1. 
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * 
     * This has been updated to remove loading the reference genome. This was to retrieve the sequence
     * of the miRNA and pre-miRNA based on the coordinates from the GFF3 file. However, since the sequence is 
     * specified in the 'mature.fa' and 'hairpin.fa' files, this seems a rather pointless exercise that introduces
     * unnecessary overhead.
     * 
     *
     * 
     * can also perform the following steps:
     *   2. Load miRNA sequence from Mature.fa file
     *   3. Load pre-miRNA sequence from hairpin.fa
     *   4. Load pre-miRNA secondary structure from miRNA.str
     * 
     * 
     * @param host String
     * @param gffFilePath String   
		 * 
		 * The loadSequence option has been removed from this method because, rather than retrieving
		 * the miRNA and pre-miRNA sequence from the reference genome, we can load from the 
		 * 'hairpin.fa' and 'mature.fa' files.
		 * The reference genome may be useful for probing the potential isomiRs, but we can still get this 
		 * information from the pre-mIRNA sequence.
		 * 
     * @throws IOException thrown if error encountered when reading the miRBase file
     * 
     */
    public void loadMiRBaseDataFromGFF(String gffFilePath) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        
        try{
            String nextLine = null;
            String thisLine = null;
            
            ArrayList<String> lines = new ArrayList<>();
            logger.info("reading miRBase gff reference file <" +  gffFilePath + ">");
            BufferedReader brMiR = new BufferedReader(new FileReader(new File(gffFilePath)));
                nextLine = brMiR.readLine();            
                do{

                    if(nextLine == null) break;
                    
                    if(nextLine.startsWith("#")){
                        this.parseGFFHeaderLine(nextLine);
                        nextLine = brMiR.readLine();                        
                        continue;
                    }
                    
                    if(nextLine.contains("miRNA_primary_transcript")) {
                        do {
                            lines.add(nextLine);  
                            nextLine = brMiR.readLine();
                            if(nextLine == null || nextLine.contains("miRNA_primary_transcript")) break;
                            
                        }while(true);
                        

                        
                        PreMiRNA preMiRNA = new PreMiRNA();
                        preMiRNA.parseMiRBaseGFFEntry(lines);
//                        if(preMiRNA.getMiID().equals("MI0007996")){
//                          logger.debug("time for a break");
//                        }
                        /* code to retrieve miRNA and pre-miRNA sequence from reference sequence
                         * 
                        if(loadSequence){
                            preMiRNA.setSeq(genomeFasta.getSubSeq(preMiRNA.getChromosome(), 
                                    StrandString.guessStrand(preMiRNA.getStrand()), 
                                    preMiRNA.getStartPos()-1, 
                                    preMiRNA.getEndPos()-1));
                            for(MiRNA miRNA: preMiRNA.getMiRNAList()){
                                miRNA.setSeq(genomeFasta.getSubSeq(miRNA.getChromosome(), 
                                    StrandString.guessStrand(miRNA.getStrand()), 
                                    miRNA.getStartPos()-1, 
                                    miRNA.getEndPos()-1));
                            }                            
                        }
                        */
                        this.getPreMiRNAList().add(preMiRNA);
                        
                        lines.clear();
                       
                    }

                }while(true);
            brMiR.close();
            logger.info("read " + getPreMiRNAList().size() + " pre-miRNA entries");
        }
        catch(IOException ex){
            logger.error("error reading miRBase gff reference file <" +  gffFilePath + ">\n" + ex.toString());
            throw new IOException("error reading miRBase gff reference file <" +  gffFilePath + ">");
        }
        this.setGffLoaded(true);
    }
    
    
    /**
     * write list in FastA format to a single FA file
     * @param fileName
     * @throws IOException 
     */
    public void writePreMiRAsSingleFA(String fileName) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
            
      String faString = "";
    	try {
	    	for(PreMiRNA preMiR: this.getPreMiRNAList()) {
          faString = faString.concat(preMiR.toFastA());
	    	}
    		BufferedWriter bwFA = new BufferedWriter(new FileWriter(new File(fileName)));
          bwFA.write(faString);
	    	bwFA.close();
    	}catch(IOException exIO) {
    		logger.error("error writing pre-miRNAs to FA file <" + fileName + ">");
    		exIO.printStackTrace();
    		throw new IOException("error writing pre-miRNAs to FA file <" + fileName + ">");
    	}

    }
  
  
    /**
     * write pre-miRNAs in FastA format, separate file for each feature.
     * @param fileName
     * @throws IOException 
     */
    public void writePreMiRAsMultipleFAs() throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      
      if((!new File(outputFolder).exists())){
        logger.info(INDENT2 + "creating output folder <" + getOutputFolder() + ">");
        new File(getOutputFolder()).mkdir();
      }else {
        logger.info(INDENT2 + "found output folder <" + getOutputFolder() + ">");
      }
      String preMiRFAFile = "";
    	try {
	    	for(PreMiRNA preMiR: this.getPreMiRNAList()) {
          preMiRFAFile  = FilenameUtils.concat(getOutputFolder(), preMiR.getMiID() + ".fa"); 
          BufferedWriter bwFA = new BufferedWriter(new FileWriter(new File(preMiRFAFile)));
            bwFA.write(preMiR.toFastA());
          bwFA.close();
        }

    	}catch(IOException exIO) {
    		logger.error("error writing pre-miRNA to FA file <" + preMiRFAFile + ">");
    		exIO.printStackTrace();
    		throw new IOException("error writing pre-miRNAs to FA file <" + preMiRFAFile + ">");
    	}

    }
    
    /**
     * write out the preMiRs in FastA format
     * 
     * @throws IOException 
     */
    public void writePreMirsAsFA() throws IOException{
      
      if(singleSeqFile){
        this.writePreMiRAsSingleFA(FilenameUtils.concat(
                FilenameUtils.getFullPath(this.getGffFilepath()), 
                FilenameUtils.getBaseName(this.getGffFilepath()) + ".fa"));
      }else{
        this.writePreMiRAsMultipleFAs();
      }
    }
    
    
    
    
    
    /**
     * output list in FastA format, but only pre-miRNAs with structure information
     */
    public void writePreMiRWithStructureInfosAsFA(String fileName) throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
    	try {
    		BufferedWriter bwFA = new BufferedWriter(new FileWriter(new File(fileName)));
	    	for(PreMiRNA preMiR: this.getPreMiRNAList()) {
	    		if(preMiR.hasStructureInfo()) {
	    			bwFA.write(preMiR.toFastA());
	    		}else {
	    			logger.info("--no structure information for pre-miRNA [" + preMiR.getName());
	    		}
	    	}
	    	bwFA.close();
    	}catch(IOException exIO) {
    		logger.error("error writing pre-miRNAs to FA file <" + fileName + ">");
    		exIO.printStackTrace();
    		throw new IOException("error writing pre-miRNAs to FA file <" + fileName + ">");
    	}

    }
    
    
    /**
     * read genome sequence from fasta file
     * 
     * @throws IOException 
     */
    public void loadGenomeFA() throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      if(getGenomeFAFilepath()!=null ){
        try{
          logger.info("reading genome file <" + getGenomeFAFilepath() + ">");
          genomeSeq.readFastaGenome(getGenomeFAFilepath());
          logger.info("finished ");
          logger.info("read " + genomeSeq.getNoOfBases() + " bases");
          logger.info("spanning " + genomeSeq.getNoOfChr() + " chromosomes");
        }
        catch(IOException exIO){
          logger.error("exception reading Genome reference file + <" + getGenomeFAFilepath() + ">");
          logger.error(exIO.toString());
          throw new IOException("exception reading Genome reference file + <" + getGenomeFAFilepath() + ">");
        }

      }else{
        logger.warn("no Genome FA file specified");
       
      }
      
    }
    
    /**
     * load name maps file that maps feature names in genome FA file to miRBase feature name used in s GFF file
     * 
     * @throws IOException 
     */
    public void loadNameMaps() throws IOException{
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      try{
        BufferedReader brNM = new BufferedReader(new FileReader(new File(this.getNamesMapFilepath())));
            String lineNM = null;
             
            while ((lineNM = brNM.readLine())!=null){
              if(lineNM.startsWith("#"))
                continue;
              MiRBase2GenomefeaturesMap.put(lineNM.split("\t")[0].trim(), lineNM.split("\t")[1].trim());
              this.Genome2MiRBasefeaturesMap.put(lineNM.split("\t")[1].trim(), lineNM.split("\t")[0].trim());
              logger.info(lineNM);
            }        
        brNM.close();
        logger.info("read <" + MiRBase2GenomefeaturesMap.size() + "> features");
      }
      catch(IOException ex){
        logger.error("error reading miRBase fasta reference file <" +  this.getNamesMapFilepath() + ">\n" + ex.toString());
        throw new IOException("error reading miRBase fasta reference file <" +  this.getNamesMapFilepath() + ">");
      }
    }
    
    
    /**
     * grab useful header information from header line. We are interested in the 
     * following lines:
     * 
     *   # microRNAs:               miRBase v21								
     *   # genome-build-id:         GRCh38								
     *   # genome-build-accession:  NCBI_Assembly:GCA_000001405.15	
     * 
     * @param headerLine 
     */
    private void parseGFFHeaderLine(String headerLine){
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());

        if(headerLine.contains(MIRBASEVERSION))
        {
            //# microRNAs:               miRBase v21	
            try{
                this.setMirbaseVersion(Integer.parseInt(headerLine.split(MIRBASEVERSION)[1].split("v")[1].trim()));
                logger.info("mirBase version is "  + this.getMirbaseVersion());
            }
            catch(NumberFormatException exNF){
                logger.info("error parsing miRBase version from miRBase header");
                logger.info(headerLine);
                logger.info("  -->  " + headerLine.split(MIRBASEVERSION)[1].split("v")[1].trim());
                logger.error("error parsing miRBase version from miRBase header");
                logger.error(headerLine);
                logger.error("  -->  " + headerLine.split(MIRBASEVERSION)[1].split("v")[1].trim());
                throw new NumberFormatException("error parsing miRBase version from miRBase header" + headerLine);
            }
            return;
        }
        
        if(headerLine.contains(GENOME_BUILD_ID))
        {
            this.setBuildID(headerLine.split(GENOME_BUILD_ID)[1].trim());
            logger.info("genome build ID is " + this.getBuildID());
            return;
        }
        
        if(headerLine.contains(GENOME_BUILD_ACC))
        {
            this.setBuildAccession(headerLine.split(GENOME_BUILD_ACC)[1]);
            logger.info("genome Accession Number is " + this.getBuildAccession());
        }
        
    }


    /**
     * cycle through the GFF entries and get the corresponding sequence from the
     * genome sequence
     */
    public void gff2seq(){
      logger.info(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      for(PreMiRNA preMiR: this.getPreMiRNAList()){
        logger.info(preMiR.getMiID()); 
        if(preMiR.getMiID().equals("MI0022143"))
          logger.info("time to rest");
        String chr = miRFeature2GenomeFeature(preMiR.getChromosome());
        if(chr==null) chr = preMiR.getChromosome();
        if(chr!=null){
//          preMiR.setSeq(this.genomeSeq.getSubSeq(miRFeature2GenomeFeature(preMiR.getChromosome()), StrandString.guessStrand(preMiR.getStrand()), preMiR.getStartPos(), preMiR.getEndPos()));
          logger.info(preMiR.getChromosome() + "from:\t" + preMiR.getStartPos() + " --->\t" + "to\t" + preMiR.getEndPos());
          preMiR.setSeq(this.genomeSeq.getSubSeq(chr, StrandString.guessStrand(preMiR.getStrand()), preMiR.getStartPos(), preMiR.getEndPos()));
        }else{
          preMiR.setSeq("");
        }
        logger.info(preMiR.getSeq());
      }

    }

    
    
    /**
     * get corresponding genomeFeatureName for specified miRBaseFeatureName
     * @param miRBaseFeatureName
     * @return 
     */
    public String miRFeature2GenomeFeature(String featureName){
      String genomeFeature = MiRBase2GenomefeaturesMap.get(featureName);
      return genomeFeature;
    }
    
    
    
    /**
     * get corresponding genomeFeatureName for specified miRBaseFeatureName
     * @param miRBaseFeatureName
     * @return 
     */
    public String GenomeFeature2miRFeature(String featureName){
      String genomeFeature = Genome2MiRBasefeaturesMap.get(featureName);
      return genomeFeature;
    }
    
    
    
    /**
     * search in the list for the specified miRNA miRBase name
     * returns the index of the parent pre-miRNA, or -1 if not found
     * @param miRNAname
     * @return
     */
    public int findMiRNAEntryByName(String miRNAname) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
      	for(MiRNA miR: preMiRNA.getMiRNAList()) {
      		if(miR.getMimatID().equals(miRNAname)){
            return i;
      		}
      	}
          if(preMiRNA.getMiID().equals(miRNAname)){
              return i;
          }
          i++;
      }
    	
    	return -1;
    }
    
    
    /**
     * search in the list for the specified miRNA miRBase ID
     * returns the index, or -1 if not found
     * @param miRNAname
     * @return
     */
    public int findMiRNAEntryByID(String miRNAID) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
        	for(MiRNA miR: preMiRNA.getMiRNAList()) {
        		if(miR.getMimatID().equals(miRNAID)){
            return i;
          }
        }
      }
    	
    	
    	return -1;
    }
    
    
    
    /**
     * return the miRNA in the list for the specified miRNA miRBase ID
     * @param miRNAname
     * @return MiRNA / null if not found
     */
    public MiRNA getMiRNAEntryByID(String miRNAID) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
        for(MiRNA miR: preMiRNA.getMiRNAList()) {
          if(miR.getMimatID().equals(miRNAID)){
            return miR;
          }
        }
      }
    	
    	
    	return null;
    }
    
    
    
    
    
    /**
     * search in the list for the specified miRNA miRBase ID
     * returns the index, or -1 if not found
     * @param miRNAname
     * @return
     */
    public int findPreMiRNAEntryByID(String preMiRNAID) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
    		if(preMiRNA.getMiID().equals(preMiRNAID)){
            return i;
      	  }
      }
    	
    	
    	return -1;
    }
    
    
    
    
    /**
     * search in the list for the specified miRNA miRBase Name
     * returns the index, or -1 if not found
     * @param miRNAname
     * @return
     */
    public int findPreMiRNAEntryByName(String preMiRNAName) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
    		if(preMiRNA.getName().equals(preMiRNAName)){
            return i;
          }
      }
    	
    	
    	return -1;
    }
    
    
    
    
    /**
     * search in the list for the specified pre-miRNA miRBase ID
     * returns the pre-miRNA, or null if not found
     * @param miRNAname
     * @return
     */
    public PreMiRNA getPreMiRNAEntryByID(String preMiRNAID) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());

      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
    		if(preMiRNA.getMiID().equals(preMiRNAID)){
            return preMiRNA;
        }
      }
    	
    	
    	return null;
    }
    
    
    
    
    /**
     * search in the list for the specified pre-miRNA miRBase Name
     * returns the pre-miRNA, or null if not found
     * @param miRNAname
     * @return
     */
    public PreMiRNA getPreMiRNAEntryByName(String preMiRNAName) {
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());

      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
    		if(preMiRNA.getName().equals(preMiRNAName)){
            return preMiRNA;
      	  }
      }
    	
    	
    	return null;
    }
    
    
    
    
    /**
     * empty the contents of the feature set
     */
    public void clearPreMiRNAList(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
          preMiRNA.getFeatureSet().clear();
          preMiRNA.getMiRNAList().clear();
          preMiRNA.getPubmedRefList().clear();
      }
        this.getPreMiRNAList().clear();
    }
    
    
    
    
    /**
     * remove miRNA entries that correspond to NGS reads. This leaves the
     * miRBase entries intact. MiRBase entries will have zero read count, 
     * NGS reads must have at least one read count. This is how we distinguish
     * the two groups.
     * The feature set and PubMed list remains untouched as these are describing 
     * features of the parent pre-miRNA
     * 
     */
    public void removeNGSReadEntries(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
          List<MiRNA> found = new ArrayList<>();
          for(MiRNA miRNA: preMiRNA.getMiRNAList()){
              if(miRNA.getReadCount()!=0)
                  found.add(miRNA);
          }
          preMiRNA.getMiRNAList().removeAll(found);
      }
        
    }
    
    
    
    
    /**
     * Search for this MiID in the feature list 
     * (the MiID should be unique, so it is reasonable to search against this parameter)
     * 
     * @param  qMiID String
     * @return int index of match
     */
    public int findEntryIndexByMiID(String qMiID){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      int i=0;
      for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
          if(preMiRNA.getMiID().equals(qMiID)){
              return i;
          }
          i++;
      }
      return -1;
    }
    
    
    /**
     * find an pre-miRNA entry by MiID
     * @param qMiID int
     * @return PreMiRNA instance match
     */
    public PreMiRNA getEntryByMiID(String qMiID){
        return this.getPreMiRNAList().get(this.findEntryIndexByMiID(qMiID));
    }
    
    
    /**
     * Is this MiID in the feature list?
     * @param qMiID int 
     * @return Boolean T/F match
     */
    public Boolean isMiIDInList(String qMiID){
        return this.findEntryIndexByMiID(qMiID) > -1;
    }
    
    /**
     * add entry to the pre-miRNA list
     * 
     * @param priMiRNA PreMiRNA to add
     * @return int length of the list
     * 
     */
    public int addPriMiRNA(PreMiRNA priMiRNA){
        this.getPreMiRNAList().add(priMiRNA);
        return this.getPreMiRNAList().size();
    }
    
    
    
    /**
     * sum the reads over all pre-miRNA entries
     * 
     * @return long total read count
     */
    public long countUpReadsInAllPreMiRNAs(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        long totalCounts = 0L;
        totalCounts 
                = this.getPreMiRNAList().stream().map(
                        (preMiRNA) -> preMiRNA.countReads()).reduce(totalCounts, (accumulator, _item) -> accumulator + _item);

        return totalCounts;
    }
    
    /**
     * add the miRNA to the pre-miRNA specified by the miID
     * 
     * @param miRNA MiRNA to add
     * @param miID String mimatID of the miRNA to be added
     */
    public void addMiRNAtoPreMiRNA(MiRNA miRNA, String miID){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        int i = this.findEntryIndexByMiID(miID);
        this.getPreMiRNAList().get(i).getMiRNAList().add(miRNA);
        this.getPreMiRNAList().get(i).addToReadCount(miRNA.getReadCount());

    }
    
    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start int 
     * @param stop int 
     * @param chr String
     * @param strand String
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return MiRNAFeature the MiRNA that the read overlaps
     */
    public MiRNA doesReadOverlapKnownMiRNA(int start, int stop, String chr, String strand, int bleed){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
      
        for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
            for(MiRNA miRRNAEntry: preMiRNA.getMiRNAList()){
                if (miRRNAEntry.chromosomeMatch(chr)){
                    if(strand.equals(miRRNAEntry.getStrand())){
                        if( java.lang.Math.abs(start - miRRNAEntry.getStartPos()) <= bleed){

                            if( java.lang.Math.abs(stop - miRRNAEntry.getEndPos()) <= bleed){
                                return miRRNAEntry;
                            }

                        }

                    }
                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start int
     * @param stop int 
     * @param chr String 
     * @param strand String
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return PreMiRNA the PreMiRNA that the read overlaps
     */
    public PreMiRNA doesReadOverlapKnownPreMiRNA(int start, int stop, String chr, String strand, int bleed){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        
        ChromosomeString chrString = new ChromosomeString(chr);
        for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
            if (ChromosomeString.GetChromNumber(preMiRNA.getChromosome())==chrString.getChrNo()){
                if(strand.equals(preMiRNA.getStrand())){
                    /*
                    this is different to the miRNA situation. Here, we never expect a read
                    to span the entire feature, we just specify how much it can "stick out" past 
                    the ends
                    */
                    if(preMiRNA.getStartPos() - start  <= bleed){

                        if(stop - preMiRNA.getEndPos() <= bleed){
                            return preMiRNA;
                        }

                    }

                }
            }
        }
        
        return null;
        
    }
    
    
    
    /**
     * does the pre-miRNA list contain the specified miRNA MIMATID?
     * @param qMimatID int 
     * @return PreMiRNA the pre-miRNA containing the miRNA
     */
    public PreMiRNA doesListContainMiRNA(String qMimatID){        
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        
        for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
            for(MiRNA miRRNAEntry: preMiRNA.getMiRNAList()){
                if (miRRNAEntry.getMimatID().equals(qMimatID)){
                    return preMiRNA;
                }
            }
        }
        return null;
    }
    
    
    
    /**
     * get the number of entries in this release
     * 
     * @return int the number of preMiRNAs in the list
     */
    public int getNumberOfPreMiRNAs(){
        return getPreMiRNAList().size();
    }

    
    
    /**
     * count up all the miRNA entries in the list
     * 
     * @return int the total number of miRNAs in the list
     */
    public int getNumberOfmiRNAs(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        
        //int miRCount = 0;
        int miRCount = this.getPreMiRNAList().stream().map((preMiRNA) -> preMiRNA.getMiRNAList().size()).reduce(0, Integer::sum);
        return miRCount;
    }
    
    
    
    /**
     * find the longest pre-miRNA in the set
     * 
     * @return int longest pre-miRNA
     */
    public int getLongestPreMiRNA(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        int longestPriMiRNA = 0;
        for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
            if(preMiRNA.getSeq().length() > longestPriMiRNA)
                longestPriMiRNA = preMiRNA.getSeq().length();
        }
        return longestPriMiRNA;
    }
    
    
    
    /**
     * find the longest miRNA in the set
     * @return int Longest miRNA
     */
    public int getLongestMiRNA(){
      logger.debug(INDENT1 + new Object(){}.getClass().getEnclosingMethod().getName());
        
        int longestMiRNA = 0;
        for(PreMiRNA preMiRNA: this.getPreMiRNAList()){
            for(MiRNA miRRNA: preMiRNA.getMiRNAList()){
            if(miRRNA.getSeq().length() > longestMiRNA)
                longestMiRNA = miRRNA.getSeq().length();
            }
        }
        return longestMiRNA;
    }
    /**
     * @return String buildID
     */
    public String getBuildID() {
        return buildID;
    }

    /**
     * @param buildID String the buildID to set
     */
    public void setBuildID(String buildID) {
        this.buildID = buildID;
    }

    /**
     * @return String buildAccession
     */
    public String getBuildAccession() {
        return buildAccession;
    }

    /**
     * @param buildAccession String the buildAccession value
     */
    public void setBuildAccession(String buildAccession) {
        this.buildAccession = buildAccession;
    }

    /**
     * @return int mirbaseVersion
     */
    public int getMirbaseVersion() {
        return mirbaseVersion;
    }

    /**
     * @param mirbaseVersion int the mirbaseVersion version
     */
    public void setMirbaseVersion(int mirbaseVersion) {
        this.mirbaseVersion = mirbaseVersion;
    }

    /**
     * @return List preMiRNAList
     */
    public List<PreMiRNA> getPreMiRNAList() {
        return preMiRNAList;
    }

    /**
     * 
     * @param refHost String
     * @param miRBaseGFFFile String
     * @param hairpinFileMirBase String
     * @throws IOException thrown when error encountered reading either the miRBase
     *                     GFF file or FA file
     */
    public void loadMiRBaseData(String refHost, String miRBaseGFFFile, String hairpinFileMirBase) throws IOException {
        HashMap <String, String> preMiRBaseSeq = new HashMap();
        try{
            logger.info("reading miRBase haripin fasta reference file <" +  hairpinFileMirBase + ">");
            BufferedReader brFA = new BufferedReader(new FileReader(new File(hairpinFileMirBase)));
            String lineFA = null;
            String preMimatID=null;
            String hairpinSeq = null;
            
            while ((lineFA = brFA.readLine())!=null){
                if (lineFA.startsWith(">")){
                    String entryHost = lineFA.split(" ")[0].substring(1).split("-")[0].trim();
                    if(entryHost.equals(refHost)){
                        logger.debug("premiRNA ID:Seq --- "+preMimatID + "\t"+hairpinSeq);
                        preMiRBaseSeq.put(preMimatID, hairpinSeq);
                        hairpinSeq = null;
                        preMimatID=null;
                        preMimatID = lineFA.split(" ")[1].trim();
                    }
                }
                else{
                    if (hairpinSeq != null){
                        hairpinSeq = hairpinSeq.concat(lineFA);
                    }
                    else{
                        hairpinSeq = lineFA;
                    }
                }

//                String seq = brFA.readLine().trim();
//                logger.debug("mydebug1 --- "+seq);
//                String entryHost = lineFA.split(" ")[0].substring(1).split("-")[0].trim();
//                if(entryHost.equals(refHost)){
//                    String preMimatID = lineFA.split(" ")[1].trim();
//                    logger.debug("mydebug2 --- "+preMimatID);
//                    preMiRBaseSeq.put(preMimatID, seq);
//                }

            }
            brFA.close();
            logger.info("read " + preMiRBaseSeq.size() + " entries");
            logger.info("done");

        }
        catch(IOException ex){
            logger.error("error reading miRBase fasta reference file <" +  hairpinFileMirBase + ">\n" + ex.toString());
            throw new IOException("error reading miRBase fasta reference file <" +  hairpinFileMirBase + ">");
        }
        
        
        
        try{

            String line = null;
            logger.info("reading miRBase gff reference file <" +  miRBaseGFFFile + ">");
            BufferedReader brMiR = new BufferedReader(new FileReader(new File(miRBaseGFFFile)));
                while((line = brMiR.readLine())!= null){

                    if(line.startsWith("#")) continue;
                    if(line.contains("miRNA_primary_transcript")) continue;
                    /*
                        chr1            chromosome
                        source          n/a here
                        miRNA           feature type (n/a)
                        start pos
                        end pos
                        score           n/a here               
                        strand          (+/-)
                        frame           n/a here
                        attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p;Derives_from=MI0022705

                    */
                    PreMiRNA preMiRNA = new PreMiRNA();
                  
                    preMiRNA.parseMiRBaseEntry(line);
                    
//                    MiRNA miRNA = new MiRNA();
//                    miRNA.parseMiRBaseEntry(line);
                    // this may need revising. In some cases we dont care if there is no sequence, we still want the feature
//                    String alias =preMiRNA.getMimatAlias();
                    String alias =preMiRNA.getMiID();

                    String seq = preMiRBaseSeq.get(alias);
                    logger.debug("to see preseq of mirna --- " + seq);
                    String preMirName = preMiRNA.getName();

                    preMiRNA.setSeq(seq);
                    boolean isDuplicate = checkForDuplicates(preMirName, this.preMiRNAList);
                    if (isDuplicate){
                        preMiRNA.setName(preMirName+".2");
                    }
                    this.preMiRNAList.add(preMiRNA);
                    if(seq == null) {
                        logger.warn("no sequence found for entry <" + alias + ">. Skipping");
                    }

                }
            brMiR.close();
            logger.info("read " + preMiRNAList.size() + " miRNA entries");
        }
        catch(IOException ex){
            logger.error("error reading miRBase gff reference file <" +  miRBaseGFFFile + ">\n" + ex.toString());
            throw new IOException("error reading miRBase gff reference file <" +  miRBaseGFFFile + ">");
        }
    }


    /**
     * check whether the query pre-MiRNA has duplicate entries in the loaded pre-miRNA list
     * @param preMirName        String query pre-miRNA
     * @param preMiRNAList  List of PreMiRNA entries
     * @return Boolean T/F is duplicate
     */
    private boolean checkForDuplicates(String preMirName, List<PreMiRNA> preMiRNAList) {
        boolean isDuplicate = false;
        for(PreMiRNA preMiRBaseEntry: preMiRNAList){
            if (preMirName.equals(preMiRBaseEntry.getName())) {
                isDuplicate = true;
            }
        }
        return isDuplicate;
    }

		public Boolean getGffLoaded() {
			return gffLoaded;
		}

		private void setGffLoaded(Boolean gffLoaded) {
			this.gffLoaded = gffLoaded;
		}

		public Boolean getMatureLoaded() {
			return matureLoaded;
		}

		private void setMatureLoaded(Boolean matureLoaded) {
			this.matureLoaded = matureLoaded;
		}

		public Boolean getHairpinLoaded() {
			return hairpinLoaded;
		}

		private void setHairpinLoaded(Boolean hairpinLoaded) {
			this.hairpinLoaded = hairpinLoaded;
		}

		private Boolean getSStructureLoaded() {
			return structureLoaded;
		}

		private void setStructureLoaded(Boolean structureLoaded) {
			this.structureLoaded = structureLoaded;
		}

		public String getGffFilepath() {
			return gffFilepath;
		}

		public String getMatureFilepath() {
			return matureFilepath;
		}

		public void setMatureFilepath(String matureFilepath) {
			this.matureFilepath = matureFilepath;
		}

		public String getHairpinFilepath() {
			return hairpinFilepath;
		}

		public void setHairpinFilepath(String hairpinFilepath) {
			this.hairpinFilepath = hairpinFilepath;
		}

		public String getStructureFilepath() {
			return structureFilepath;
		}

		public void setStructureFilepath(String structureFilepath) {
			this.structureFilepath = structureFilepath;
		}

		public int getNoOfmiRNAs() {
			return noOfmiRNAs;
		}

		public void setNoOfmiRNAs(int noOfmiRNAs) {
			this.noOfmiRNAs = noOfmiRNAs;
		}

  /**
   * @return the organismCode
   */
  public String getOrganismCode() {
    return organismCode;
  }

  /**
   * @param organismCode the organismCode to set
   */
  public void setOrganismCode(String organismCode) {
    this.organismCode = organismCode;
  }

  /**
   * @return the genomeFAFilepath
   */
  public String getGenomeFAFilepath() {
    return genomeFAFilepath;
  }

  /**
   * @param genomeFAFilepath the genomeFAFilepath to set
   */
  public void setGenomeFAFilepath(String genomeFAFilepath) {
    this.genomeFAFilepath = genomeFAFilepath;
  }

  /**
   * @param gffFilepath the gffFilepath to set
   */
  public void setGffFilepath(String gffFilepath) {
    this.gffFilepath = gffFilepath;
  }

  /**
   * @return the singleSeqFile
   */
  public Boolean getSingleSeqFile() {
    return singleSeqFile;
  }

  /**
   * @param singleSeqFile the singleSeqFile to set
   */
  public void setSingleSeqFile(Boolean singleSeqFile) {
    this.singleSeqFile = singleSeqFile;
  }

  /**
   * @return the namesMapFilepath
   */
  public String getNamesMapFilepath() {
    return namesMapFilepath;
  }

  /**
   * @param namesMapFilepath the namesMapFilepath to set
   */
  public void setNamesMapFilepath(String namesMapFilepath) {
    this.namesMapFilepath = namesMapFilepath;
  }

  /**
   * @return the outputFolder
   */
  public String getOutputFolder() {
    return outputFolder;
  }

  /**
   * @param outputFolder the outputFolder to set
   */
  public void setOutputFolder(String outputFolder) {
    this.outputFolder = outputFolder;
  }
}
