/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

/**
 *
 * @author sr
 */
public class MiRNADispAnalysisResult {
    
    private String          mimatID;
    private double          span;
    private double          pValue;
    private String          note;

    
    public MiRNADispAnalysisResult(String mID, double s, double p, String n){
        mimatID = mID;
        span = s;
        pValue = p;
        note = n;
    }
    
    
    /**
     * print header line with column names
     * 
     * @return String tab delimited header line
     */
    public static String printHeaderLine(){
        
        return "mimatID\tspan\tp-value\tnote\n";
        
    }
    
    
    /**
     * print summary line for this instance
     * 
     * @return String tab delimited summary line
     */
    public String printSummary(){
        return mimatID + "\t" + span + "\t" + pValue + "\t" + note + "\n";
    }
    
    /**
     * @return the mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    /**
     * @return the pValue
     */
    public double getpValue() {
        return pValue;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @return the span
     */
    public double getSpan() {
        return span;
    }
    
    
}
