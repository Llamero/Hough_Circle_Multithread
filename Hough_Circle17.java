/** Hough_Circle.java:
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 @author Ben Smith (benjamin.smith@berkeley.edu)
 @Based on original plugin implementation by Hemerson Pistori (pistori@ec.ucdb.br) and Eduardo Rocha Costa
 @created February 4, 2017
 
 The Hough Transform implementation was based on 
 Mark A. Schulze applet (http://www.markschulze.net/)
 
*/

import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.awt.*;
import ij.plugin.HyperStackConverter;
import ij.plugin.frame.*;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.*;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Ben
 */
public class Hough_Circle implements PlugInFilter {

    //***GUI input variables***
    // <editor-fold desc="Initialize variables">
    //Search parameters
    private int radiusMin;  // Find circles with radius grater or equal radiusMin - argument syntax: "min=#"
    private int radiusMax;  // Find circles with radius less or equal radiusMax - argument syntax: "max=#"
    private int radiusInc;  // Increment used to go from radiusMin to radiusMax - argument syntax: "inc=#"
    private int minCircles; // Minumum number of circles to be found - argument syntax: "minCircles=#"    
    private int maxCircles; // Maximum number of circles to be found - argument syntax: "maxCircles=#"
    private int threshold = -1; // An alternative to maxCircles. All circles with a value in the hough space greater then threshold are marked. Higher thresholds result in fewer circles. - argument syntax: "threshold=#"
    private double thresholdRatio; //Ratio input from GUI that expresses threshold as ratio of resolution (highest possible # of votes)
    private int resolution; //The number of steps to use per transform (i.e. number of voting rounds)
    private double ratio; // Ratio of found circle radius to clear out surrounding neighbors
    private int searchBand = 0; //The +/- range of radii to search for relative to the last found radius - argument syntax: "bandwidth=#"
    private int searchRadius = 0; //The search radius to look for the next centroid relative to the last found centroid - argument syntax: "radius=#"
    private boolean reduce = false; //Cap the transform resolution by removeing redundant steps
    private boolean local = false; //Whether or not the search is going to be local
    
    //Output parameters
    private boolean houghSeries = false; //Contains whether the user wants the Hough series stack as an output - argument syntax: "show_raw"
    private boolean showCircles = false; //Contains whether the user wants the circles found as an output - argument syntax: "show_mask"
    private boolean showRadius = false; //Contains whether the user wants a map of centroids and radii outputed from search - argument syntax: "show_centroids"
    private boolean showScores = false; //Contains whether the user wants a map of centroids and Hough scores outputed from search - argument syntax: "show_scores"
    private boolean results = false; //Contains whether the user wants to export the measurements to a reuslts table 
    
    //Hough transform variables
    public ImagePlus imp; //Initalize the variable to hold the image
    private boolean isStack = false; //True if there is more than one slice in the input data
    private int stackSlices; //number of slices in the stack
    private int houghIndex; //Contains the current index (slice) position in the Hough radius series
    private int maxHough; //Contains the brights pixel in the entire Hough array
    private Point maxPoint; //Stores the location of the brightest pixel in the Hough array
    private int maxRadius; //Stores the radius of the brightest pixel in the Hough array
    private Rectangle r; //Stores the ROI on the original image
    private float imageValues[]; // Raw image (returned by ip.getPixels()) - float is used to allow 8, 16 or 32 bit images
    private int houghValues[][][]; // Hough Space Values [X coord][Y coord][radius index]
    private int localHoughValues[][][][]; //Local Hough space [circle#][X coord][Y coord][radius index]
    private int localHoughParameters[][]; //Array to pass local Hough space parameters to centroid search [circle#][parameter vector]
    private int width; // Hough Space width (depends on image width)
    private int height;  // Hough Space heigh (depends on image height)
    private int depth;  // Hough Space depth (depends on radius interval)
    private int offset; // Image Width
    private int offx;   // ROI x offset
    private int offy;   // ROI y offset
    private Point centerPoint[]; // Center Points of the Circles Found.
    private int centerRadii[]; //Corresponding radii of the cricles marked by the center points
    private int houghScores[]; //Corresponding Hough scores for each centroid
    private int circleID[]; //Corresponding ID # for each centroid
    private int circleIDcounter = 0; //Counter for keeping track of current max ID #
    private int lut[][][]; // LookUp Table for x and y tranform shifts in an octahedral manner
    private int lutSize; //Stores the actual number of transforms performed (<=selected resolution)
    private int nCircles = 0; //Number of circles found during search - <= maxCircles 
    
    //Variables for storing the results and exporting result images
    public ImagePlus houghPlus;
    public ImageStack houghStack;
    public ImagePlus circlePlus;
    public ImageStack circleStack;
    public ImagePlus radiusPlus;
    public ImageStack radiusStack;
    public ImagePlus scorePlus;
    public ImageStack scoreStack;
    public ImageProcessor circlesip;
    public ImageProcessor radiusip;
    public ImageProcessor scoresip;
    public ResultsTable rt; 
    
    //Variables for max Hough score search
    private int maxHoughArray[][]; //Matrix to store hough scores, radii, and points from multi-threaded max search
    private int ithread;
    // </editor-fold>
    
    @Override
    //Initialize the Help plugin menu
    public int setup(String arg, ImagePlus imp) {        
        
        if (arg.equals("about")) {
            showAbout();
            return DONE;
        }
        
        //Sends arduments to ImageJ that tells it how to run the plugin - tells it to accept all grayscale and supports selections
        return DOES_8G+DOES_16+DOES_32+SUPPORTS_MASKING; //Do not include DOES_STACKS, as this will call the GUI once for each slice
    }
    
    void showAbout() {
        IJ.showMessage("Hough Circle Transform v1.0.0",
                       "This plugin performs a Hough circle transform \n" +
                       "on an image or stack.  Hough circle transforms\n" +
                       "can be used to find the centroid and radius of\n" +
                       "circles embedded in complex images.  This plugin\n"+
                       "was inspired by the transform implementation\n"+
                       "of Hemerson Pistori (pistori@ec.ucdb.br)"
                      );
    }
    
    @Override
    //Start running the plugin - Grab ROI and then get transform parameters
    public void run(ImageProcessor ip) {

        //Initialize the results table
        rt = Analyzer.getResultsTable();
        if (rt == null) {
                rt = new ResultsTable();
                Analyzer.setResultsTable(rt);
        }
        
        //Show a Dialog Window for user input of parameters
        readParameters(); 
    }    

    //Create the GUI or extract parameters from macro arguments
    void readParameters() {
        // <editor-fold desc="Retrieve macro arguments">
        //See if any arguments have been passed to the plugin via a macro
        if (IJ.isMacro() && ij.Macro.getOptions() != null && !ij.Macro.getOptions().trim().isEmpty()){
            String[] arguments = ij.Macro.getOptions().split(",");

            //remove all spaces from array components
            for(int a = 0; a<arguments.length; a++){
                arguments[a] = arguments[a].trim();
            }
            
            //Parse the arguments to the corresponding variables
            for (String argument : arguments) { //passes each argment in the array to the variable "argument"
                if (argument.matches(".*minRadius.*=.*")) {
                    //Retrieve min radius
                    radiusMin = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*maxRadius.*=.*")) {
                    //Retrieve max radius
                    radiusMax = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*inc.*=.*")) {
                    //Retrieve radius increment
                    radiusInc = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*minCircles.*=.*")) {
                    //Retrieve number of circles
                    minCircles = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*maxCircles.*=.*")) {
                    //Retrieve number of circles
                    maxCircles = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*threshold.*=.*")) {
                    //Retrieve Hough score threshold
                    //Code from: http://www.itgo.me/a/x8683194173055835262/get-float-or-integer-value-from-the-string-in-java
                    Pattern pattern = Pattern.compile("\\d+(?:\\.\\d+)?"); // Match int or float
                    Matcher matcher = pattern.matcher(argument);
                    if(matcher.find()) thresholdRatio = Double.parseDouble(matcher.group());
                }
                else if (argument.matches(".*resolution.*=.*")) {
                    //Retrieve Hough score threshold
                    resolution = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*ratio.*=.*")) {
                    //Retrieve clearing neighbor radius ratio
                    //Code from: http://www.itgo.me/a/x8683194173055835262/get-float-or-integer-value-from-the-string-in-java
                    Pattern pattern = Pattern.compile("\\d+(?:\\.\\d+)?"); // Match int or float
                    Matcher matcher = pattern.matcher(argument);
                    if(matcher.find()) ratio = Double.parseDouble(matcher.group());
                }
                else if (argument.matches(".*bandwidth.*=.*")) {
                    //Retrieve Hough score threshold
                    searchBand = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*local_radius.*=.*")) {
                    //Retrieve Hough score threshold
                    searchRadius = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }

                //Retrieve checkbox status
                if (argument.matches(".*reduce.*")) reduce = true;
                if (argument.matches(".*local.*")) local = true;
                if (argument.matches(".*show_raw.*")) houghSeries = true;
                if (argument.matches(".*show_mask.*")) showCircles = true;
                if (argument.matches(".*show_centroids.*")) showRadius = true;
                if (argument.matches(".*show_scores.*")) showScores = true;
                if (argument.matches(".*results_table.*")) results = true;
            }
            
            //Start the Hough Transform
            startTransform();
          
        }
        // </editor-fold>
        //If arguments are not coming from macro, present user interface
        else{
            // <editor-fold desc="Swing GUI part 1.">
            //***Build GUI using Swing***
            //Initialize a panel (needed to house a frame)
            final JPanel guiPanel = new JPanel();

            //Initlialize a frame to hold the gui
            final JFrame guiFrame = new JFrame();

            //Set the frame to close when the window is closed
            guiFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            
            //Initialize Swing GUI components
            final JLabel guiTitle = new javax.swing.JLabel();
            final JLabel guiIntro1 = new javax.swing.JLabel();
            final JLabel guiIntro2 = new javax.swing.JLabel();
            final JLabel guiSearchLabel = new javax.swing.JLabel();
            final JLabel guiMinLabel = new javax.swing.JLabel();
            final JTextField guiMinText = new javax.swing.JTextField();
            final JLabel guiMaxLabel = new javax.swing.JLabel();
            final JTextField guiMaxText = new javax.swing.JTextField();
            final JLabel guiIncLabel = new javax.swing.JLabel();
            final JTextField guiIncText = new javax.swing.JTextField();
            final JLabel  guiMinNumLabel = new javax.swing.JLabel();
            final JTextField guiMinNumText = new javax.swing.JTextField();
            final JLabel guiMaxNumLabel = new javax.swing.JLabel();
            final JTextField guiMaxNumText = new javax.swing.JTextField();        
            final JLabel guiThreshLabel = new javax.swing.JLabel();
            final JTextField guiThreshText = new javax.swing.JTextField();
            final JLabel guiResLabel = new javax.swing.JLabel();
            final JTextField guiResText = new javax.swing.JTextField();
            final JLabel guiClearLabel = new javax.swing.JLabel();
            final JTextField guiClearText = new javax.swing.JTextField();
            final JLabel guiRadiusBandLabel = new javax.swing.JLabel();
            final JTextField guiRadiusBandText = new javax.swing.JTextField();
            final JLabel guiSearchRadLabel = new javax.swing.JLabel();
            final JTextField guiSearchRadText = new javax.swing.JTextField();
            final JCheckBox guiReduceBox = new javax.swing.JCheckBox();
            final JCheckBox guiLocalBox = new javax.swing.JCheckBox();
            final JLabel guiOutputLabel = new javax.swing.JLabel();
            final JCheckBox guiRawBox = new javax.swing.JCheckBox();
            final JCheckBox guiPointBox = new javax.swing.JCheckBox();
            final JCheckBox guiRadiusBox = new javax.swing.JCheckBox();
            final JCheckBox guiHoughBox = new javax.swing.JCheckBox();
            final JCheckBox guiResultsBox = new javax.swing.JCheckBox();
            final JButton guiOKButton = new javax.swing.JButton();
            
            //Format GUI text
            guiTitle.setFont(new java.awt.Font("Tahoma", 0, 24)); // NOI18N
            guiTitle.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiTitle.setText("Local Hough Circle Transform");

            guiIntro1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiIntro1.setText("This plugin performs a Hough circle transform on an image or stack. ");

            guiIntro2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiIntro2.setText("It can be used to find and measure circular objects within an image.");

            guiSearchLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
            guiSearchLabel.setText("Search parameters:");

            guiMinLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMinLabel.setText("Minimum search radius (in pixels):");

            guiMinText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMinText.setText("10");
            guiMinText.setToolTipText("Radius of smallest circle you expect to find");

            guiMaxLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMaxLabel.setText("Maximum search radius (in pixels):");

            guiMaxText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMaxText.setText("100");
            guiMaxText.setToolTipText("Radius of largest circle you expect to find");

            guiIncLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiIncLabel.setText("Radius search increment (in pixels):");

            guiIncText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiIncText.setText("2");
            guiIncText.setToolTipText("How much to increase the radius at each increment");

            guiMinNumLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMinNumLabel.setText("Minimum number of circles to be found:");

            guiMinNumText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMinNumText.setText("0");
            guiMinNumText.setToolTipText("The minimum number of circles you expect to find in the image");

            guiMaxNumLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiMaxNumLabel.setText("Maximum number of circles to be found:");

            guiMaxNumText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiMaxNumText.setText("1");
            guiMaxNumText.setToolTipText("The maximum number of circles you expect to find in the image");

            guiThreshLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiThreshLabel.setText("Hough score threshold (0.0 - 1.0):");
            guiThreshLabel.setToolTipText("This is the ratio of votes/resolution");

            guiThreshText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiThreshText.setText("0.5");
            guiThreshText.setToolTipText("The Hough score threshold above which circles are valid");

            guiResLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiResLabel.setText("Transform resolution (# of steps per transform):");
            guiResLabel.setToolTipText("The maximum number of transform steps per radius");

            guiResText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiResText.setText("1000");
            guiResText.setToolTipText("The number of steps to use per transform");

            guiClearLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiClearLabel.setText("Clear neighbors radius ratio:");
            guiClearLabel.setToolTipText("The ratio relative to the circle's radius");

            guiClearText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiClearText.setText("1.0");
            guiClearText.setToolTipText("The radius around a found circle to clear out the Hough search space");

            guiRadiusBandLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiRadiusBandLabel.setText("Local radius search bandwidth (+/- previous radius):");

            guiRadiusBandText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiRadiusBandText.setText("10");
            guiRadiusBandText.setToolTipText("The range of radii to be used in the next frame");

            guiSearchRadLabel.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
            guiSearchRadLabel.setText("Local radius search area (search radius near last centroid):");

            guiSearchRadText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiSearchRadText.setText("10");
            guiSearchRadText.setToolTipText("The area around a centroid to look for the next centroid");

            guiReduceBox.setSelected(true);
            guiReduceBox.setText("Reduce transform resolution by removing redundant transform steps");

            guiLocalBox.setSelected(true);
            guiLocalBox.setText("Use a local search space to speed up search (only applies to movies)");

            guiOutputLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
            guiOutputLabel.setText("Output options:");

            guiRawBox.setText("Raw Hough transform series");

            guiPointBox.setSelected(true);
            guiPointBox.setText("Circle centroid(s) marked on the original image");

            guiRadiusBox.setText("Map of circle radius at centroids (pixel intensity = circle radius)");

            guiHoughBox.setText("Map of circle score at centroids (pixel intensity = Hough score)");

            guiResultsBox.setSelected(true);
            guiResultsBox.setText("Export measurements to the results table");

            guiOKButton.setText("OK");
            guiOKButton.addActionListener((java.awt.event.ActionEvent evt) -> {
                // </editor-fold>
                // <editor-fold desc="Retrieve GUI arguments and calculate Hough parameters">
                //Retrive the numbers from the text boxes and combobox
                radiusMin = Integer.parseInt(guiMinText.getText());
                radiusMax = Integer.parseInt(guiMaxText.getText());
                radiusInc = Integer.parseInt(guiIncText.getText());
                minCircles = Integer.parseInt(guiMinNumText.getText()); 
                maxCircles = Integer.parseInt(guiMaxNumText.getText());
                thresholdRatio = Double.parseDouble(guiThreshText.getText());
                resolution = Integer.parseInt(guiResText.getText());
                ratio = Double.parseDouble(guiClearText.getText());
                searchBand = Integer.parseInt(guiRadiusBandText.getText());
                searchRadius = Integer.parseInt(guiSearchRadText.getText());
                reduce = guiReduceBox.isSelected();
                local = guiLocalBox.isSelected();
                
                //Retrieve the check box status
                houghSeries = guiRawBox.isSelected();
                showCircles = guiPointBox.isSelected();
                showRadius = guiRadiusBox.isSelected();
                showScores = guiHoughBox.isSelected();
                results = guiResultsBox.isSelected();
                
                //Remove the GUI frame now that it is no longer needed
                //guiFrame.dispose();
                
                //Start the Hough Transform
                startTransform();
            } //When push button is pushed, retrieve the state of
            );
            // </editor-fold>
            // <editor-fold desc="Swing GUI Part 2">
            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(guiFrame.getContentPane());
            guiFrame.getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(guiSearchLabel)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(guiOutputLabel)
                        .addGap(303, 303, 303))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(guiTitle, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiIntro1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiIntro2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiHoughBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiRadiusBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiRawBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiPointBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiResultsBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiLocalBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(guiReduceBox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(guiRadiusBandLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(guiRadiusBandText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(guiSearchRadLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(guiSearchRadText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addGap(0, 0, Short.MAX_VALUE)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                            .addComponent(guiMinLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(guiMaxLabel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(guiIncLabel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(guiMinNumLabel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(guiMinNumText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(guiIncText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(guiMaxText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(guiMinText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                    .addComponent(guiOKButton, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                            .addComponent(guiClearLabel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(guiResLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 237, Short.MAX_VALUE)
                                            .addComponent(guiThreshLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                .addComponent(guiThreshText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addComponent(guiResText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                                            .addComponent(guiClearText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(guiMaxNumLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(guiMaxNumText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addContainerGap())))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(guiTitle)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiIntro1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiIntro2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(guiSearchLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMinLabel)
                    .addComponent(guiMinText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMaxLabel)
                    .addComponent(guiMaxText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiIncText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiIncLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMinNumText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiMinNumLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiMaxNumText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiMaxNumLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiThreshText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(guiThreshLabel))
                .addGap(6, 6, 6)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiResLabel)
                    .addComponent(guiResText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiClearLabel)
                    .addComponent(guiClearText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiRadiusBandLabel)
                    .addComponent(guiRadiusBandText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(guiSearchRadLabel)
                    .addComponent(guiSearchRadText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiReduceBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiLocalBox)
                .addGap(18, 18, 18)
                .addComponent(guiOutputLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(guiRawBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiPointBox, javax.swing.GroupLayout.PREFERRED_SIZE, 23, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiRadiusBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiHoughBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiResultsBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(guiOKButton)
                .addContainerGap())
        );
                       
            //Use the Pack function to make the GUI frame shrink to the smallest size necassary
            guiFrame.pack();

            //Show the GUI
            guiFrame.setVisible(true); 
            // </editor-fold>
        }
    }
    
    public void startTransform(){
        
        //Calculate Hough parameters
        depth = ((radiusMax-radiusMin)/radiusInc)+1;
        
        //If radiusInc is not a divisor of radiusMax-radiusMin, return error
        if((radiusMax-radiusMin)%radiusInc != 0){
            IJ.showMessage("Radius increment must be a divisor of maximum radius - minimum radius.");
            IJ.showMessage("radiusMin=" + radiusMin + ", radiusMax=" + radiusMax + ", radiusInc=" + radiusInc);
            return;          
        }
        
        //Build the transform LUT (all necessary translations in Cartesian coordinates)
        //NOTE: This step must precede the calculatation of the threshold, as it can change the resolution from the original input value
        lutSize = buildLookUpTable();
        
        //Calculate the threshold based off the the actual resolution
        threshold = (int) Math.round(thresholdRatio * resolution);
        
        //If the threshold rounds to 0, set threshold to 1 as 0 is nonsensical
        if(threshold <= 0) {
            threshold = 1;
        }        
        // <editor-fold desc="Send arguments to record">
        //If the macro was being recorded, return the set values to the recorder
        if (Recorder.record){
            String Command  = "run(\"Local Hough Circle Transform\",\"minRadius=" + radiusMin + ", maxRadius=" + radiusMax + ", inc=" + radiusInc + 
                    ", minCircles=" + minCircles + ", maxCircles=" + maxCircles + ", threshold=" + thresholdRatio + ", resolution=" + resolution +
                     ", ratio=" + ratio + ", bandwidth=" + searchBand + ", local_radius=" + searchRadius + ", ";
            if(reduce) Command += " reduce";
            if(local) Command += " local";
            if(houghSeries) Command += " show_raw";
            if(showCircles) Command += " show_mask";
            if(showRadius) Command += " show_centroids";
            if(showScores) Command += " show_scores";
            if(results) Command += " results_table";
            Command += "\");\r\n";
            Recorder.recordString(Command);
        }
        // </editor-fold>

        //Create an ImagePlus instance of the currently active image
        imp = WindowManager.getCurrentImage();
        
        //Import the ImagePlus as a stack
        ImageStack stack = imp.getStack();
        
        //Get the ROI dimensions - no ROI = full image
        r = stack.getRoi();
        offx = r.x;
        offy = r.y;
        width = r.width;
        height = r.height;
        offset = stack.getWidth();

        //Convert the stack to float (allows 8, 16, and 32 stacks to all be processed as one type)
        ImageStack stackCopy = stack.duplicate();
        ImageStack floatStack = stackCopy.convertToFloat();

        if(houghSeries){
            //Frames is transform dimension
            houghStack = new ImageStack(width, height, depth*imp.getNSlices());
        }
        if(showCircles){
            circleStack = new ImageStack(width, height, imp.getNSlices());
        }
        if(showRadius){
            radiusStack = new ImageStack(width, height, imp.getNSlices());
        }
        if(showScores){
            scoreStack = new ImageStack(width, height, imp.getNSlices());
        }

        //See if input is stack
        stackSlices = imp.getStackSize();
        if(stackSlices > 1) isStack = true;
        
        //Build arrays for storing the circle result parameters
        houghScores = new int [maxCircles];
        centerRadii = new int [maxCircles];
        centerPoint = new Point[maxCircles];
        circleID = new int [maxCircles];
        Arrays.fill(circleID, -1); //Flag indeces as -1 = no circle found for this index
        
        //Start timer to delay log outputs
        long startTime = System.currentTimeMillis();
        
        
int totalTime = 0;    
        
        
        
        //Retrieve the current slice in the stack as an ImageProcessor
        for(int slice = 1; slice <= stackSlices; slice++){
         
            //Show tranform status if there is more than one slice
            if(isStack){
               if((System.currentTimeMillis() - startTime) > 500){
                   IJ.showStatus("Hough tranform - processing frame: " + slice + " of " + stackSlices + "");
                   
                   //Reset timer
                   startTime = System.currentTimeMillis();
               }
            }
            else{
               IJ.log("\\Clear");
               IJ.log("Running Hough tansfrom...\r\n");
            }
 
            ImageProcessor ip = floatStack.getProcessor(slice);
            imageValues = (float[]) ip.getPixels();

            // <editor-fold desc="Local search">
            //If the number of circles is greater than min circles - speed up the search by looking locally for each next circle
            //Multithread by giving each core a subset of the total number of circles
            if(local & slice > 1){
                //If there are still a sufficient number of circles, perform an exclusively local search
                if(nCircles>=minCircles){
                    
                    
long startSearch = System.currentTimeMillis();



                    startLocalTransform();
                    
                    
                    
totalTime += System.currentTimeMillis() - startSearch;
                    //Rebuild the local searches into a full transform if the Hough Series is desired
                    if(houghSeries){
                        convertLocaltoFullHough();
                        HoughSpaceSeries(slice, houghStack);
                    }
                }
                
                //Otherwise, if there still are valid circles, perform the full Hough transform, but perform a local search for the valid circles first before continuing onto a full search 
                else if(nCircles > 0){
                    houghTransform();
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                    startPartialLocalSearch();
                    getCenterPoints();
                }
                 //If there are no circles, and minCircles is 0, then return blank results
                else if(nCircles == 0 & minCircles == 0){
                    //Create an empty Hough array
                    houghValues = new int[width][height][depth];
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                } 

                //If no circles are found, then revert to a full Hough
                else{
                    houghTransform();
                    if (houghSeries){
                        //Create the hyperstach to put into the result if needed
                        HoughSpaceSeries(slice, houghStack);
                    }
                    // Mark the center of the found circles in a new image if user wants to find centers
                    if(showCircles || showRadius || showScores || results) getCenterPoints();
                }
            }
            //Otherwise, perform the full transform
            else{
                houghTransform();
                if (houghSeries){
                    //Create the hyperstach to put into the result if needed
                    HoughSpaceSeries(slice, houghStack);
                }
                // Mark the center of the found circles in a new image if user wants to find centers
                if(showCircles || showRadius || showScores || results) getCenterPoints();
            }
             // </editor-fold>
            
            if(!isStack) IJ.log("Searching Hough space for circles...\r\n");

            // Create image View for Marked Circles.
            if(showCircles) drawCircles(slice, circleStack, width, height, offx, offy, offset);
            
            //Create map of centroids where the intensity is the radius
            if (showRadius || showScores) drawCentroids(slice, radiusStack, scoreStack);
            
            //Export measurements to the results table
            if (results) resultsTable(slice);
            
        }
IJ.log("Total time: " + totalTime + ", Time per circle: " + totalTime/nCircles);
        //Draw the resulting stacks
         if(houghSeries){
            houghPlus = new ImagePlus("Hough Transform Series", houghStack);
            if(depth>1) houghPlus = HyperStackConverter.toHyperStack(houghPlus, 1, depth, imp.getNSlices(), "default", "grayscale");
            houghPlus.show();
         }
         if(showCircles) new ImagePlus("Centroid overlay", circleStack).show();
         if(showRadius) new ImagePlus("Centroid map", radiusStack).show();
         if(showScores) new ImagePlus("Score map", scoreStack).show(); 
         if(results) rt.show("Results");
         
         IJ.log("Hough transform complete!");
    }
    
    private void startLocalTransform(){
         //Initialize an array to store the local Hough transform from each thread
        localHoughValues = new int [nCircles][2*searchRadius + 1][2*searchRadius + 1][2*searchBand + 1];
        
        //Initialize an array to store the final dimensions and location of each local Hough Transform
        localHoughParameters = new int [nCircles][9];

        //Setup two processing pipelines to multithread only if there are more than 1 circles
        if(nCircles < 2){

            //Divide the task so that each core works on a subset of circles
            for(int circleNum = 0; circleNum < nCircles; circleNum++){
                localHoughTransform(circleNum);
                localGetCenterPoint(circleNum); //Make sure this returns -1 to index if no circle was found
            }
        }
        
        //Otherwise, multithread the local search
        else{
            //Build an array to store the result from each thread
            final Thread[] threads = newThreadArray();

            //Create an atomic integer counter that each thread can use to count through the radii
            final AtomicInteger ai = new AtomicInteger(0);  

            //Build a thread for as many CPUs as are available to the JVM 
            for (ithread = 0; ithread < threads.length; ithread++) {    

                // Concurrently run in as many threads as CPUs  
                threads[ithread] = new Thread() {  

                    { setPriority(Thread.NORM_PRIORITY); }  

                    @Override
                    public void run() {  

                        //Divide the task so that each core works on a subset of circles
                        for(int circleNum = ai.getAndAdd(1); circleNum < nCircles; circleNum = ai.getAndAdd(1)){
                            localHoughTransform(circleNum);
                            localGetCenterPoint(circleNum); //Make sure this returns -1 to index if no circle was found
                        }
                    }               
                };  
            }    
            startAndJoin(threads);
        }

        
        //Flag move all circles to the starting indexes, so there are no gaps between circles (i.e. if there were 3 circles, then they occupy indeces 0-2)
        collapseLocalResult();
    }
    
    private void startPartialLocalSearch(){
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();

        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(0);  

        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    

            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  

                { setPriority(Thread.NORM_PRIORITY); }  

                @Override
                public void run() {  

                    //Divide the task so that each core works on a subset of circles
                    for(int circleNum = ai.getAndAdd(1); circleNum < nCircles; circleNum = ai.getAndAdd(1)){
                        getLocalCenterPoint2(circleNum);
                    }
                }               
            };  
        }    
        startAndJoin(threads);
        
        //Flag move all circles to the starting indexes, so there are no gaps between circles (i.e. if there were 3 circles, then they occupy indeces 0-2)
        collapseLocalResult();
    }
    
    //Build an array of the cartesion transformations necessary for each increment at each radius
    private int buildLookUpTable() {
        //Build an array to store the X and Y coordinates for each angle increment (resolution) across each radius (depth)
        lut = new int[2][resolution][depth];
        
        //Initialize a variable that will allow for measuring the maximium LUT size
        int maxLUT = 0;
        
        //Step through all radii to be sampled
        for(int radius = radiusMax; radius>=radiusMin; radius -= radiusInc) {
            
            //Index counter that also tracks the largest actual LUT size (may be <= resolution)
            int i = 0; 
            for(int resStep = 0; resStep < resolution; resStep++) {
                
                //Calcualte the angle and corresponding X and Y displacements for the specified radius
                double angle = (2D*Math.PI * (double)resStep) / (double)resolution;
                int indexR = (radius-radiusMin)/radiusInc;
                int rcos = (int)Math.round ((double)radius * Math.cos (angle));
                int rsin = (int)Math.round ((double)radius * Math.sin (angle));
                
                //Test to make sure that the coordinate is a new coordinate
                //NOTE: A continuous circle is being discretized into pixels, it is possible that small angle steps for small radii circles will occupy the same pixel.
                //Since there is no point in making redundant calculations, these points are excluded from the LUT
                //NOTE: Using the minRadius as the transform cutoff results in a strong harmonic of the image forming near the min radius
                //threfore, using the max will push this harmonic outside the search range.
                if(radius == radiusMax && reduce){
                    if(i == 0) {
                        lut[0][i][indexR] = rcos;
                        lut[1][i][indexR] = rsin;
                        i++;
                    }
                    else if ((rcos != lut[0][i-1][indexR]) | (rsin != lut[1][i-1][indexR])){
                        lut[0][i][indexR] = rcos;
                        lut[1][i][indexR] = rsin;
                        i++;
                    }
                }
                else{
                    lut[0][i][indexR] = rcos;
                    lut[1][i][indexR] = rsin;
                    i++;
                }
            }
            
            //If this is the smallest radius, see how many transforms could be done, and set this as the new resolution
            if(radius == radiusMax){
                maxLUT = i;
                resolution = maxLUT;             
            }
        }
        return maxLUT;
    }

    
    private void houghTransform () {
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(radiusMin);  
        
        //Create an array to store the Hough values
        houghValues = new int[width][height][depth];
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  
  
                @Override
                public void run() {  
  
                //Divide the radius tasks across the cores available
                for (int radius = ai.getAndAdd(radiusInc); radius <= radiusMax; radius = ai.getAndAdd(radiusInc)) {  
                    int indexR=(radius-radiusMin)/radiusInc;
                    //For a given radius, transform each pixel in a circle, and add-up the votes 
                    for(int y = 1; y < height-1; y++) {
                        for(int x = 1; x < width-1; x++) {
                                if( imageValues[(x+offx)+(y+offy)*offset] != 0 )  {// Edge pixel found                                    
                                    for(int i = 0; i < lutSize; i++) {

                                        int a = x + lut[1][i][indexR]; 
                                        int b = y + lut[0][i][indexR]; 
                                        if((b >= 0) & (b < height) & (a >= 0) & (a < width)) {
                                            houghValues[a][b][indexR] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }  
            };  
        }    
        startAndJoin(threads);      
    }
    
    //The local Hough is an inversion of the inertial reference frame used in the full Hough
    //In the full Hough the image is seach for pixels with a value > 1, and if this is true
    //then the pixel is projected in a circle about that point.
    
    //To reduce the necessary transform space,
    private void localHoughTransform (int index) {
        //Initialize local search variables
        int startWidth = centerPoint[index].x - searchRadius;
        if (startWidth < 1) startWidth = 1; //Keep search area inside the image
        int endWidth = centerPoint[index].x + searchRadius;
        if (endWidth > width) endWidth = width; //Keep search area inside the image
        int lwidth = endWidth - startWidth + 1;
        
        int startHeight = centerPoint[index].y - searchRadius;
        if (startHeight < 1) startHeight = 1; //Keep search area inside the image
        int endHeight = centerPoint[index].y + searchRadius;
        if (endHeight > height) endHeight = height; //Keep search area inside the image
        int lheight = endHeight - startHeight + 1;
        
        int lradiusMin = centerRadii[index] - searchBand;
        if(lradiusMin < radiusMin) lradiusMin = radiusMin;
        int lradiusMax = centerRadii[index] + searchBand;
        if(lradiusMax > radiusMax) lradiusMax = radiusMax;
        int ldepth = ((lradiusMax-lradiusMin)/radiusInc)+1;
        
        //Store local Hough parameters into array
        localHoughParameters[index][0] = startWidth;
        localHoughParameters[index][1] = endWidth;
        localHoughParameters[index][2] = lwidth;
        localHoughParameters[index][3] = startHeight;
        localHoughParameters[index][4] = endHeight;
        localHoughParameters[index][5] = lheight;
        localHoughParameters[index][6] = lradiusMin;
        localHoughParameters[index][7] = lradiusMax;
        localHoughParameters[index][8] = ldepth;
        
        //Divide the radius tasks across the cores available
        for (int radius = lradiusMin; radius <= lradiusMax; radius += radiusInc) { 
            int indexR=(radius-lradiusMin)/radiusInc;
            //For a given radius, transform each pixel in a circle, and add-up the votes 
            for(int y = startHeight; y < endHeight-1; y++) {
                for(int x = startWidth; x < endWidth-1; x++) {
                    for(int i = 0; i < lutSize; i++){
                        int a = x + lut[1][i][indexR+((lradiusMin-radiusMin)/radiusInc)]; 
                        int b = y + lut[0][i][indexR+((lradiusMin-radiusMin)/radiusInc)];
                        
                        //Check to make sure pixel is within the image
                        if((a>1) & (a<width-1) & (b>1) & (b<height-1)){
                            //See if the pixel has an intensity > 1, if so, add 1 to the vote
                            if(imageValues[a + (b*width)] > 0 ){
                                localHoughValues[index][x-startWidth][y-startHeight][indexR] += 1;
                            } 
                        }
                    }
                }
            }
        }
    }  

    //Find the largest Hough pixel in the 3D Hough transform array to scale the 8-bit conversion
    private void houghMaximum () {
        //long startTime = System.currentTimeMillis(); //1337ms without multi, 319ms with multi, 175ms by writing to private variable per thread
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Build an array to store the max from each thread
        maxHoughArray = new int[threads.length][4];
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(0);
        
        //Create an integer for indexing the results array (one result per thread
        final AtomicInteger Index = new AtomicInteger(0);

        maxHough = 0;
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  

                //Search for the largest score with each thread
                @Override
                public void run() {
                    int maxHoughThread = -1;
                    int maxRadiusThread = -1;
                    Point maxPointThread = new Point (-1,-1);
                    for(int a=ai.getAndIncrement(); a<depth; a=ai.getAndIncrement()){
                        for(int j = 0; j < height; j++) {
                            for(int k = 0; k < width; k++){
                                if(houghValues[k][j][a] > maxHoughThread) {
                                    maxHoughThread = houghValues[k][j][a];
                                    maxPointThread = new Point(k,j);
                                    maxRadiusThread = a*radiusInc + radiusMin;
                                }
                            }
                        }
                    }
                    //Have each thread report the score to a common array
                    maxHoughArray[Index.getAndIncrement()] = new int[]{maxHoughThread, maxRadiusThread, maxPointThread.x, maxPointThread.y};                    
                }
            };
        }
        startAndJoin(threads);
        
        //Search common array for highest score
        for (int[] maxHoughArray1 : maxHoughArray) {
            if (maxHoughArray1[0] > maxHough) {
                maxHough = maxHoughArray1[0];
                maxRadius = maxHoughArray1[1];
                maxPoint = new Point((int) maxHoughArray1[2], (int) maxHoughArray1[3]);    
            }
        }
        //long stopTime = System.currentTimeMillis();
        //IJ.log("Elapsed time was " + (stopTime - startTime) + " miliseconds.");
    }
    
    private void convertLocaltoFullHough(){
        
        //Create an array to store the Hough values
        houghValues = new int[width][height][depth];
        
        for(int circle=0; circle<localHoughParameters.length; circle++){
            //Extract the local search parameters
            int startWidth = localHoughParameters[circle][0];
            int endWidth = localHoughParameters[circle][1];
            int lwidth = localHoughParameters[circle][2];
            int startHeight = localHoughParameters[circle][3];
            int endHeight = localHoughParameters[circle][4];
            int lheight = localHoughParameters[circle][5];
            int lradiusMin = localHoughParameters[circle][6];
            int lradiusMax = localHoughParameters[circle][7];
            int ldepth = localHoughParameters[circle][8];
            for(int w=0; w<lwidth; w++){
                for(int h=0; h<lheight; h++){
                    System.arraycopy(localHoughValues[circle][w][h], 0, houghValues[w+startWidth][h+startHeight], (lradiusMin-radiusMin)/radiusInc, ldepth);
                }
            }
        }
    }
    
    //Add transform series data to the hyperstack
    private void HoughSpaceSeries(int slice, ImageStack houghStack){
    
        //If transform was performed, find the brightest pixel in the whole Hough transform array to normalize stack intensity
        if(minCircles > 0 | nCircles > 0) houghMaximum();
        
        int startFrame = ((slice-1)*depth);

        for(int radius = radiusMin; radius<=radiusMax; radius += radiusInc) {
            //Create a new image for depositing the Hough pixels into
            ImageProcessor newip = new ByteProcessor(width, height);
            byte[] newpixels = (byte[])newip.getPixels();

            //Calculate the corresponding index
            houghIndex = (radius-radiusMin)/radiusInc;

            //If tansform was performed, retrieve the pixel array for the current Hough radius image
            if(minCircles > 0 | nCircles > 0) createHoughPixels(newpixels, houghIndex);
            
            //Deposit the array image into the corresponding slice in the stack
            houghStack.setPixels(newpixels, houghIndex+1+startFrame);

            //Give the current slice the appropriate radius label
            houghStack.setSliceLabel("Hough Space [r="+radius+", resolution="+resolution+"]", houghIndex+1+startFrame);
        }   
    }
    
    // Convert Values in Hough Space to an 8-Bit Image Space.
    private void createHoughPixels (byte houghPixels[], int index) {
	//Rescale all the Hough values to 8-bit to create the Hough image - 47ms to complete - single threading okay
        for(int l = 0; l < height; l++) {
            for(int i = 0; i < width; i++) {
                houghPixels[i + l * width] = (byte) Math.round ((houghValues[i][l][index] * 255D) / maxHough);
            }

        }
    }

    // Draw the circles found in the original image.
    private void drawCircles(int slice, ImageStack circleStack, int widthROI, int heightROI, int offsetX, int offsetY, int offsetROI) {
		
            // Copy original input pixels into output
            // circle location display image and
            // combine with saturation at 100
            circlesip = new ByteProcessor(widthROI, heightROI);
            byte[] circlespixels = (byte[])circlesip.getPixels();

            int roiaddr=0;
            for( int y = offsetY; y < offsetY+heightROI; y++) {
                    for(int x = offsetX; x < offsetX+widthROI; x++) {
                            // Copy;
                            circlespixels[roiaddr] = (byte) imageValues[x+offsetROI*y];
                            // Saturate
                            if(circlespixels[roiaddr] != 0 )
                                    circlespixels[roiaddr] = 100;
                            else
                                    circlespixels[roiaddr] = 0;
                            roiaddr++;
                    }
            }
            // Copy original image to the circlespixels image.
            // Changing pixels values to 100, so that the marked
            // circles appears more clear. Must be improved in
            // the future to show the resuls in a colored image.
            //for(int i = 0; i < width*height ;++i ) {
            //if(imageValues[i] != 0 )
            //if(circlespixels[i] != 0 )
            //circlespixels[i] = 100;
            //else
            //circlespixels[i] = 0;
            //}
            if(centerPoint == null) getCenterPoints();

            byte cor = -1;
            // Redefine these so refer to ROI coordinates exclusively
            offsetROI= widthROI;
            offsetX=0;
            offsetY=0;

            for(int l = 0; l < nCircles; l++) {
                    int i = centerPoint[l].x;
                    int j = centerPoint[l].y;
                    // Draw a gray cross marking the center of each circle.
                    for( int k = -10 ; k <= 10 ; ++k ) {
                            int p = (j+k+offsetY)*offsetROI+ (i+offsetX);
                            if(!outOfBounds(j+k+offsetY,i+offsetX))
                                    circlespixels[(j+k+offsetY)*offsetROI+ (i+offsetX)] = cor;
                            if(!outOfBounds(j+offsetY,i+k+offsetX))
                                    circlespixels[(j+offsetY)*offsetROI  + (i+k+offsetX)] = cor;
                    }
                    for( int k = -2 ; k <= 2 ; ++k ) {
                            if(!outOfBounds(j-2+offsetY,i+k+offsetX))
                                    circlespixels[(j-2+offsetY)*offsetROI+ (i+k+offsetX)] = cor;
                            if(!outOfBounds(j+2+offsetY,i+k+offsetX))
                                    circlespixels[(j+2+offsetY)*offsetROI+ (i+k+offsetX)] = cor;
                            if(!outOfBounds(j+k+offsetY,i-2+offsetX))
                                    circlespixels[(j+k+offsetY)*offsetROI+ (i-2+offsetX)] = cor;
                            if(!outOfBounds(j+k+offsetY,i+2+offsetX))
                                    circlespixels[(j+k+offsetY)*offsetROI+ (i+2+offsetX)] = cor;
                    }
            }
            circleStack.setPixels(circlespixels, slice);
            circleStack.setSliceLabel(nCircles + " circles found", slice);
    }
    
    // Draw the centroids found in the original image where intensity = radius.
    private void drawCentroids(int slice, ImageStack radiusStack, ImageStack scoreStack) {
            
        //Create arrays the same size as the images
        radiusip = new ShortProcessor(width, height);
        short[] radiuspixels = (short[])radiusip.getPixels();

        scoresip = new ShortProcessor(width, height);
        short[] scorepixels = (short[])scoresip.getPixels();
        
        for(int l = 0; l < nCircles; l++) {
                int i = centerPoint[l].x;
                int j = centerPoint[l].y;
                int radius = centerRadii[l];

                //Draw a point at the centroid and set the int to = radius
                if(showRadius) radiuspixels[j*offset + i] = (short) radius;

                //Record the Hough score on a 16-bit scale
                if(showScores) scorepixels[j*offset + i] = (short) Math.round(houghScores[l]/resolution*65535);
        }
        
        if(showRadius){
            radiusStack.setPixels(radiuspixels, slice);
            radiusStack.setSliceLabel(nCircles + " circles found", slice);
        }
        if(showScores){
            scoreStack.setPixels(scorepixels, slice);
            scoreStack.setSliceLabel(nCircles + " circles found", slice);
        }
    }
    
    //Export the results to the results table
    public void resultsTable(int frame){
        for(int a = 0; a < nCircles; a++){
            rt.incrementCounter();
            rt.addValue("ID", circleID[a]);
            rt.addValue("X", centerPoint[a].x);
            rt.addValue("Y", centerPoint[a].y);
            rt.addValue("Radius", centerRadii[a]);
            rt.addValue("Score", houghScores[a]/resolution);
            rt.addValue("nCircles", nCircles);
            rt.addValue("Resolution", lutSize);
            rt.addValue("Frame", frame);  
        }
    }
    

    private boolean outOfBounds(int y,int x) {
        if(x >= width)
            return(true);
        if(x <= 0)
            return(true);
        if(y >= height)
            return(true);
        if(y <= 0)
            return(true);
        return(false);
    }

    /** Search for a fixed number of circles.
    @param maxCircles The number of circles that should be found.  
    */
    private void getCenterPoints () {
   
        int countCircles = nCircles;
        maxHough = threshold;

        while(countCircles < maxCircles && maxHough >= threshold) {

            //Search for the highest remaining Hough score in the matrix
            houghMaximum();

            if(maxHough>=threshold){
                circleIDcounter++;
                circleID[countCircles] = circleIDcounter;
                centerPoint[countCircles] = maxPoint;
                centerRadii[countCircles] = maxRadius;
                houghScores[countCircles] = maxHough;
                countCircles++;
            }

            clearNeighbours(maxPoint.x,maxPoint.y, maxRadius);

        }

        nCircles = countCircles;
        
        //Clear the remainder of result arrays, since these results are from the previous round
        for(int i = nCircles; i<maxCircles; i++){
            circleID[i] = -1;
            centerPoint[i] = new Point(-1,-1);
            centerRadii[i] = -1;
            houghScores[i] = -1;
        }
    }
    
    private void getLocalCenterPoint2(int index){
        //Initialize local search variables
        int startWidth = centerPoint[index].x - searchRadius;
        if (startWidth < 1) startWidth = 1; //Keep search area inside the image
        int endWidth = centerPoint[index].x + searchRadius;
        if (endWidth > width) endWidth = width; //Keep search area inside the image
        int lwidth = endWidth - startWidth + 1;
        
        int startHeight = centerPoint[index].y - searchRadius;
        if (startHeight < 1) startHeight = 1; //Keep search area inside the image
        int endHeight = centerPoint[index].y + searchRadius;
        if (endHeight > height) endHeight = height; //Keep search area inside the image
        int lheight = endHeight - startHeight + 1;
        
        int lradiusMin = centerRadii[index] - searchBand;
        if(lradiusMin < radiusMin) lradiusMin = radiusMin;
        int lradiusMax = centerRadii[index] - searchBand;
        if(lradiusMax < radiusMax) lradiusMax = radiusMax;
        int ldepth = ((lradiusMax-lradiusMin)/radiusInc)+1;
        
        //Search locally for the highest hough score in the full Hough Space
        int maxLocalHough = -1;
        int maxLocalRadius = -1;
        Point maxLocalPoint = new Point (-1,-1);
        for(int radius = radiusMin; radius<=radiusMax; radius += radiusInc){
            int indexR=(radius-lradiusMin)/radiusInc;
            for(int j = startHeight; j < endHeight; j++) {
                for(int k = startWidth; k < endWidth; k++){
                    if(houghValues[k][j][indexR] > maxLocalHough) {
                        maxLocalHough = houghValues[k][j][indexR];
                        maxLocalPoint = new Point(k,j);
                        maxLocalRadius = radius;
                    }
                }
            }
        }
        
        //If the highest score is above the threshold, record the new circle
        if(maxHough>=threshold){
            centerPoint[index] = maxLocalPoint;
            centerRadii[index] = maxLocalRadius;
            houghScores[index] = maxLocalHough;
            clearNeighbours(maxLocalPoint.x,maxLocalPoint.y, maxLocalRadius);
        }
        //Otherwise, record that the circle was lost
        else{
            circleID[index] = -1;
            centerPoint[index] = new Point(-1,-1);
            centerRadii[index] = -1;
            houghScores[index] = -1;
        }  
    }
    
    private void localGetCenterPoint(int index){
        
        //Extract the local search parameters
        int startWidth = localHoughParameters[index][0];
        int endWidth = localHoughParameters[index][1];
        int lwidth = localHoughParameters[index][2];
        int startHeight = localHoughParameters[index][3];
        int endHeight = localHoughParameters[index][4];
        int lheight = localHoughParameters[index][5];
        int lradiusMin = localHoughParameters[index][6];
        int lradiusMax = localHoughParameters[index][7];
        int ldepth = localHoughParameters[index][8];
        
        //Search for the highest hough score in the local Hough Space
        int maxLocalHough = -1;
        int maxLocalRadius = -1;
        Point maxLocalPoint = new Point (-1,-1);
        for(int a=0; a<ldepth; a++){
            for(int j = 0; j < lheight; j++) {
                for(int k = 0; k < lwidth; k++){
                    if(localHoughValues[index][k][j][a] > maxLocalHough) {
                        maxLocalHough = localHoughValues[index][k][j][a];
                        maxLocalPoint = new Point(k + startWidth,j + startHeight);
                        maxLocalRadius = a*radiusInc + lradiusMin;
                    }
                }
            }
        }
        
        //If the highest score is above the threshold, record the new circle
        if(maxLocalHough>=threshold){            
            centerPoint[index] = maxLocalPoint;
            centerRadii[index] = maxLocalRadius;
            houghScores[index] = maxLocalHough;
        }
        //Otherwise, record that the circle was lost
        else{
            circleID[index] = -1;
            centerPoint[index] = new Point(-1,-1);
            centerRadii[index] = -1;
            houghScores[index] = -1;
        }     
    }

    private void collapseLocalResult(){
        //search for indeces containing circles (i.e. index != -1)
        int[] idIndeces = new int [circleID.length];
        int indexCounter = 0;
        for(int i=0; i<circleID.length; i++){
            //If a valid ID is found, record the index
            if(circleID[i] != -1){
                idIndeces[indexCounter] = i;
                indexCounter++;
            }
        }
       
        //Move all the found results to the starting indeces of the arrays if needed
        if(indexCounter < maxCircles){
            for(int i=0; i<indexCounter; i++){
                //If index is empty, then move the result
                if(circleID[i] == -1){
                    //Check to see index is empty
                    //Move results
                    circleID[i] = circleID[idIndeces[i]];
                    houghScores[i] = houghScores[idIndeces[i]];
                    centerRadii[i] = centerRadii[idIndeces[i]];
                    centerPoint[i] = centerPoint[idIndeces[i]];

                    //Clear original index
                    circleID[idIndeces[i]] = -1;
                    houghScores[idIndeces[i]] = -1;
                    centerRadii[idIndeces[i]] = -1;
                    centerPoint[idIndeces[i]] = new Point(-1,-1);
                }
            }
        }
        
        //Update nCircles to reflect the new found number of circles
        nCircles = indexCounter;
    };
    
    /** Clear, from the Hough Space, all the counter that are near (radius/2) a previously found circle C.
        
    @param x The x coordinate of the circle C found.
    @param x The y coordinate of the circle C found.
    @param x The radius of the circle C found.
    */
    private void clearNeighbours(int x,int y, int radius) {


        // The following code just clean the points around the center of the circle found.
	int radiusSquared = (int) Math.round(Math.pow(radius*ratio, 2D));
        //int radiusSquared = radius*radius;

        int y1 = (int)Math.floor (y - radius);
        int y2 = (int)Math.ceil (y + radius) + 1;
        int x1 = (int)Math.floor (x - radius);
        int x2 = (int)Math.ceil (x + radius) + 1;

        if(y1 < 0)
            y1 = 0;
        if(y2 > height)
            y2 = height;
        if(x1 < 0)
            x1 = 0;
        if(x2 > width)
            x2 = width;

        for(int indexR = 0; indexR<depth; indexR++){
            for(int i = y1; i < y2; i++) {
                for(int j = x1; j < x2; j++) {	      	     
                    if((int) (Math.pow (j - x, 2D) + Math.pow (i - y, 2D)) < radiusSquared) {
                        houghValues[j][i][indexR] = 0;
                    }
                }
            }
        }

    }
    
    /** Create a Thread[] array as large as the number of processors available. 
    * From Stephan Preibisch's Multithreading.java class. See: 
    * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
    */  
    private Thread[] newThreadArray() {  
        int n_cpus = Runtime.getRuntime().availableProcessors();  
        return new Thread[n_cpus];  
    }
    
    /** Start all given threads and wait on each of them until all are done. 
    * From Stephan Preibisch's Multithreading.java class. See: 
    * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
     * @param threads
    */  
    public static void startAndJoin(Thread[] threads)  
    {  
        for (int ithread = 0; ithread < threads.length; ++ithread)  
        {  
            threads[ithread].setPriority(Thread.NORM_PRIORITY);  
            threads[ithread].start();  
        }  
  
        try  
        {     
            for (int ithread = 0; ithread < threads.length; ++ithread)  
                threads[ithread].join();  
        } catch (InterruptedException ie)  
        {  
            throw new RuntimeException(ie);  
        }  
    } 

}
