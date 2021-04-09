/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

import no.uio.medisin.bag.core.mirna.IsomirPoint;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * this stores the dispersion pattern for a set of isomiRs associated with a
 * miRNA.
 * 
 * This is kept separately from the miRNAFeature class to stop the class getting
 * too big and the overhead associated with creating an instance. Also, in many 
 * cases, an isomiR analysis may not be relevant, which is an additional reason to 
 * separate out this data
 * 
 * 
 * @author simon rayner
 */
public class IsomiRSet {
    private String mimatID;
    private String runID;
    private String note;
    
    private double totalPairwiseDist    = 0.0;
    private double total3pDist          = 0.0;
    private double total5pDist          = 0.0;
    private double totalPolyDist        = 0.0;
    private int    noOfIsomiRs          = 0;
    
    
    private ArrayList<HashMap> isomiRPts;
    

    
    
    /**
     * create from dispersion entry
     * @param  dispLine String
     */
    public IsomiRSet(String dispLine){
        
        String tokens[] = dispLine.split("\t");
        runID               = tokens[0].trim();
        note                = tokens[1].trim();
        mimatID             = tokens[2].trim();
        noOfIsomiRs         = Integer.parseInt(tokens[3].trim()); 
        totalPairwiseDist   = Double.parseDouble(tokens[4].trim());
        total5pDist         = Double.parseDouble(tokens[5].trim()); 
        total3pDist         = Double.parseDouble(tokens[6].trim()); 
        totalPolyDist       = Double.parseDouble(tokens[7].trim());        
        
    }
    
    /**
     * 
     * @param mID String the mimatID
     * @param n String Note
     * @param rID String runID associated with this isomiR
     * @param isoPts ArrayList of HashMap of isomiRs identified for this miRNA
     */
    public IsomiRSet(String mID, String n, String rID, ArrayList<HashMap> isoPts){
        mimatID = mID;
        note = n;
        runID = rID;
        isomiRPts = isoPts;
    }
    
    
    
    /**
     * report the isomiR dispersion for this miRNA and experiment
     * 
     * @return String tab limited summary of this instance
     * 
     */
    public String tabReportIsomiRSet(){
        String reportStr = "";
        
        reportStr = reportStr.concat(printSummary()); 
        reportStr = reportStr.concat("5p" + "\t" + "3p" + "poly" + "\t" + "fraction" + "\t" + "TPD" + "\n");
        for(HashMap hmIsomiR: getIsomiRPts()){
            reportStr = reportStr.concat(hmIsomiR.get("5p") + "\t" + hmIsomiR.get("3p") 
                    + "\t" + hmIsomiR.get("poly") + "\t" + hmIsomiR.get("fraction") + "\n" );
        }
        return reportStr;
    }
    
    
    
    /**
     * return header string for summary line
     * 
     * @return String tab delimited header for summary file
     */
    public static String printSummaryHeader(){
        
        return "runID\tnote\tmimatID\tno of isomiRs\ttotal dist\t5p dist\t3p dist\tpoly dist\n";
        
    }
    
    
    
    
    /**
     * return summary string for this isomiR
     * 
     * @return String summary string for this isomiR
     */
    public String printSummary(){
        
        String reportStr = "";
        
        reportStr = reportStr.concat(getRunID() + "\t" + getNote() + "\t" + getMimatID() + "\t" + getIsomiRPts().size() + "\t" + totalPairwiseDist 
          + "\t" + getTotal5pDist() + "\t" + getTotal3pDist() + "\t" + getTotalPolyDist() + "\n" ); 
        
        return reportStr;
    }
    
    /**
     * calculate sum of pairwise distances between all isomiRs, scaled by the fraction that each isomiR is present
     * A rather rough way of estimating dispersion, there has to be a better way
     * 
     * 
     */
    public void calcDistParameters(){
        
        totalPairwiseDist = 0;
        for(int i=0; i<getIsomiRPts().size(); i++){
            
            HashMap hmIsomiRi = getIsomiRPts().get(i);
            for(int j=i+1; j<getIsomiRPts().size(); j++){
                
                HashMap hmIsomiRj = getIsomiRPts().get(j);
                //(double) ((Integer) marks.get(i)).intValue();
                ((Integer) 1).doubleValue();
                double ff=((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"));

                total5pDist   += Math.abs(((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("5p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                total3pDist   += Math.abs(((Integer)hmIsomiRi.get("3p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("3p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                totalPolyDist += Math.abs(((Integer)hmIsomiRi.get("poly")).doubleValue() * ((Double)hmIsomiRi.get("fraction")) - ((Integer)hmIsomiRj.get("poly")).doubleValue() * ((Double)hmIsomiRj.get("fraction")));
                totalPairwiseDist += Math.sqrt(
                  Math.pow(((Integer)hmIsomiRi.get("5p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("5p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                + Math.pow(((Integer)hmIsomiRi.get("3p")).doubleValue() * ((Double)hmIsomiRi.get("fraction"))   - ((Integer)hmIsomiRj.get("3p")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                + Math.pow(((Integer)hmIsomiRi.get("poly")).doubleValue() * ((Double)hmIsomiRi.get("fraction")) - ((Integer)hmIsomiRj.get("poly")).doubleValue() * ((Double)hmIsomiRj.get("fraction")), 2.0)
                );
                
            }
            

        }

        int ii=0;
    }
    
    
    
    
    /**
     * @return String mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    
    
    /**
     * @param mimatID the mimatID to set
     */
    public void setMimatID(String mimatID) {
        this.mimatID = mimatID;
    }

    /**
     * @return the runID
     */
    public String getRunID() {
        return runID;
    }

    /**
     * @param runID the runID to set
     */
    public void setRunID(String runID) {
        this.runID = runID;
    }

    /**
     * @return the isomiRPts
     */
    public ArrayList<HashMap> getIsomiRPts() {
        return isomiRPts;
    }

    /**
     * @param isomiRPts the isomiRPts to set
     */
    public void setIsomiRPts(ArrayList<HashMap> isomiRPts) {
        this.isomiRPts = isomiRPts;
    }

    /**
     * @return the totPairwiseDist
     */
    public double getTotPairwiseDist() {
        return totalPairwiseDist;
    }

    /**
     * @return the total3pDist
     */
    public double getTotal3pDist() {
        return total3pDist;
    }

    /**
     * @return the total5pDist
     */
    public double getTotal5pDist() {
        return total5pDist;
    }

    /**
     * @return the totalPolyDist
     */
    public double getTotalPolyDist() {
        return totalPolyDist;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }

    /**
     * @return the noOfIsomiRs
     */
    public int getNoOfIsomiRs() {
        return noOfIsomiRs;
    }
    
    
}
