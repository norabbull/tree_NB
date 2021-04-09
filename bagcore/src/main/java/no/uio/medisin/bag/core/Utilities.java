/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core;

import java.awt.Color;

/**
 *
 * @author simonrayner
 */
public class Utilities {

  /**
   * parse out the RGB values from the supplied String and create a new Color instance
   * @param colourString
   * @return : Color instance
   */
  public static Color rgbString2ColorString(String colourString) 
  {
    String bits[] = colourString.split(",");
  
    /*
    yeah, i know i should be able to do this using RegEx, but i couldn't get it to work
    */
    //Pattern c = Pattern.compile("rgb *\\( *([0-9]+), *([0-9]+), *([0-9]+) *\\)");
    //Matcher m = c.matcher(colourString);

    return new Color(Integer.valueOf(bits[0].split("\\[")[1].split("=")[1]),  // r
                     Integer.valueOf(bits[1].split("=")[1]),                  // g
                     Integer.valueOf(bits[2].split("\\]")[0].split("=")[1])); // b 
  }

  

  
}
