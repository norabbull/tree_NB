/*
 * aminoAcidProperties.java
 *
 * Created on January 3, 2009, 5:53 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

/**
 *
 * @author sr
 */
public class aminoAcidProperties implements java.io.Serializable {
   private String longName;
   private String shortName;
   private String letter;
   private String properties;
   private double hydrophobic;
   private double polar;
   private double small;
   private double proline;
   private double tiny;
   private double aliphatic;
   private double aromatic;
   private double positive;
   private double negative;
   private double charged;
   /** Creates a new instance of aminoAcidProperties */
   
   public aminoAcidProperties() {
      longName = "";
      shortName = "";
      letter = "";
      properties = "";
      hydrophobic = 0.0;
      polar = 0.0;
      small = 0.0;
      proline = 0.0;
      tiny = 0.0;
      aliphatic = 0.0;
      aromatic = 0.0;
      positive = 0.0;
      negative = 0.0;
      charged = 0.0;
   }
   
   public aminoAcidProperties(String lName, 
           String mName, String sName, String propString)
   {
      longName = lName;
      shortName = sName;
      letter = sName;
      properties = propString;
      parsePropertyString();
   }
   
   public void copy(aminoAcidProperties aap)
   {
      this.setShortName(aap.getShortName());
      this.setLongName(aap.getLongName());
      this.setLetter(aap.getLetter());
      this.setPropertyString(aap.getPropertyString());
      this.parsePropertyString();
   }
   
   public void parsePropertyString()
   {
      // property string is a ten character string summarizing basic properties
      if(properties.substring(0,1).equals("X"))
      {
         hydrophobic = 1.0;
      }
      else
      {
         hydrophobic = 0.0;
      }
      if(properties.substring(1,2).equals("X"))
      {
         polar = 1.0;
      }
      else
      {
         polar = 0.0;
      }
      if(properties.substring(2,3).equals("X"))
      {
         small = 1.0;
      }
      else
      {
         small = 0.0;
      }
      if(properties.substring(3,4).equals("X"))
      {
         proline = 1.0;
      }
      else
      {
         proline = 0.0;
      }
      if(properties.substring(4,5).equals("X"))
      {
         tiny = 1.0;
      }
      else
      {
         tiny = 0.0;
      }
      if(properties.substring(5,6).equals("X"))
      {
         aliphatic = 1.0;
      }
      else
      {
         aliphatic = 0.0;
      }
      if(properties.substring(6,7).equals("X"))
      {
         aromatic = 1.0;
      }
      else
      {
         aromatic = 0.0;
      }
      if(properties.substring(7,8).equals("X"))
      {
         positive = 1.0;
      }
      else
      {
         positive = 0.0;
      }
      if(properties.substring(8,9).equals("X"))
      {
         negative = 1.0;
      }
      else
      {
         negative = 0.0;
      }
      if(properties.substring(8,9).equals("X"))
      {
         charged = 1.0;
      }
      else
      {
         charged = 0.0;
      }
   }
   
   public String getLongName()
   {
      return longName;
   }
   public void setLongName(String l)
   {
      longName = l;
   }
   
   public String getShortName()
   {
      return shortName;
   }
   public void setShortName(String s)
   {
      shortName = s;
   }
   
   public String getLetter()
   {
      return letter;
   }
   public void setLetter(String l)
   {
      letter = l;
   }
   
   public String getPropertyString()
   {
      return properties;
   }
   public void setPropertyString(String p)
   {
      properties = p;
   }
   
   public double getHydrophobicity()
   {
      return hydrophobic;
   }
   public void setHydrophobicity(double h)
   {
      hydrophobic = h;
   }
   
   public double getPolar()
   {
      return polar;
   }
   public void setPolar(double p)
   {
      polar = p;
   }
   
   public double getSmall()
   {
      return small;
   }
   public void setSmall(double s)
   {
      small = s;
   }
   
   public double getProline()
   {
      return proline;
   }
   public void setProline(double p)
   {
      proline = p;
   }
   
   public double getTiny()
   {
      return tiny;
   }
   public void setTiny(double t)
   {
      tiny = t;
   }

   public double getAliphatic()
   {
      return aliphatic;
   }
   public void setAliphatic(double a)
   {
      aliphatic = a;
   }
   
   public double getAromatic()
   {
      return aromatic;
   }
   public void setAromatic(double a)
   {
      aromatic = a;
   }
   
   public double getPositive()
   {
      return positive;
   }
   public void setPositive(double p)
   {
      positive = p;
   }
   
   public double getNegative()
   {
      return negative;
   }
   public void setNegative(double n)
   {
      negative = n;
   }
   
   public double getCharged()
   {
      return charged;
   }
   public void setCharged(double c)
   {
      charged = c;
   }
   
}
