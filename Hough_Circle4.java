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
 @created January 1, 2017
 
 The Hough Transform implementation was based on 
 Mark A. Schulze applet (http://www.markschulze.net/)
 
*/

import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.awt.*;
import ij.gui.*;
import ij.plugin.HyperStackConverter;
import ij.plugin.frame.*;
import java.awt.event.ActionEvent;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.*;
import java.util.regex.*;
import com.google.common.util.concurrent.AtomicDouble;

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
    private int nCircles; // Numbers of circles to be found - argument syntax: "nCircles=#"
    private int threshold = -1; // An alternative to nCircles. All circles with a value in the hough space greater then threshold are marked. Higher thresholds result in fewer circles. - argument syntax: "threshold=#"
    private int resolution; //The number of steps to use per transform (i.e. number of voting rounds)
    private boolean allowOverlap = false; //Allow circles to overlap on another - argument syntax: "overlap"
    private boolean reverseSearch = false; //Search from largest sphere to smallest sphere - argument syntax: "reverse"
    
    //Output parameters
    private boolean houghSeries = false; //Contains whether the user wants the Hough series stack as an output - argument syntax: "show_raw"
    private boolean showCircles = false; //Contains whether the user wants the circles found as an output - argument syntax: "show_select"
    private boolean showRadius = false; //Contains whether the user wants a map of centroids and radii outputed from search - argument syntax: "show_centroids"
    private boolean showScores = false; //Contains whether the user wants a map of centroids and Hough scores outputed from search - argument syntax: "show_scores"
    
    //Hough transform variables
    public ImagePlus imp; //Initalize the variable to hold the image
    private int houghIndex; //Contains the current index (slice) position in the Hough radius series
    private double maxPixel; //Contains the brights pixel in the entire Hough array

    private int radiusBegin; //Sets the starting radius for Hough search
    private int radiusStep; //Sets the increment used during the Hough Search
    private float imageValues[]; // Raw image (returned by ip.getPixels()) - byte is used to allow 8 or 16 bit images
    private double houghValues[][][]; // Hough Space Values [X coord][Y coord][radius index]
    private int width; // Hough Space width (depends on image width)
    private int height;  // Hough Space heigh (depends on image height)
    private int depth;  // Hough Space depth (depends on radius interval)
    private int offset; // Image Width
    private int offx;   // ROI x offset
    private int offy;   // ROI y offset
    private Point centerPoint[]; // Center Points of the Circles Found.
    private int centerRadii[]; //Corresponding radii of the cricles marked by the center points
    private int houghScores[]; //Corresponding Hough scores for each centroid
    private final static int VECTOR_SIZE_MAX = 500; //Maximum number of circles to be found during a search
    private boolean useThreshold = false; //Varible for whether to find a number of circles, or to use a hough score threshold.
    private int lut[][][]; // LookUp Table for x and y tranform shifts in an octahedral manner
    private int lutSize; //Stores the actual number of transforms performed (<=selected resolution)
    
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
    public ImageProcessor centroidsip;
    public ImageProcessor scoresip;
    
    //Variables for max Hough score search
    private double maxArray[];
    private int ithread;
    // </editor-fold>
    
    @Override
    //Initialize the Help plugin menu
    public int setup(String arg, ImagePlus imp) {        
        //Create an ImagePlus instance of the currently active image
        this.imp = imp;
        
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
        //Get the ROI dimensions - no ROI = full image
        Rectangle r = ip.getRoi();
        offx = r.x;
        offy = r.y;
        width = r.width;
        height = r.height;
        offset = ip.getWidth();
        IJ.log("" + imp.getStackSize() + "\r\n");
        
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
                if (argument.matches(".*min.*=.*")) {
                    //Retrieve min radius
                    radiusMin = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*max.*=.*")) {
                    //Retrieve max radius
                    radiusMax = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*inc.*=.*")) {
                    //Retrieve radius increment
                    radiusInc = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*nCircles.*=.*")) {
                    //Retrieve number of circles
                    nCircles = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*threshold.*=.*")) {
                    //Retrieve Hough score threshold
                    threshold = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                else if (argument.matches(".*resolution.*=.*")) {
                    //Retrieve Hough score threshold
                    resolution = Integer.parseInt(argument.replaceAll("\\D+", "")); //Remove all non digits
                }
                //Retrieve checkbox status
                if (argument.matches(".*overlap.*")) allowOverlap = true;
                if (argument.matches(".*reverse.*")) reverseSearch = true;
                if (argument.matches(".*show_raw.*")) houghSeries = true;
                if (argument.matches(".*show_select.*")) showCircles = true;
                if (argument.matches(".*show_centroids.*")) showRadius = true;
                if (argument.matches(".*show_scores.*")) showScores = true;
            }
          
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
            final JLabel guiCirNumLabel = new javax.swing.JLabel();
            final JTextField guiCirNumText = new javax.swing.JTextField();
            final JLabel guiThreshLabel = new javax.swing.JLabel();
            final JTextField guiThreshText = new javax.swing.JTextField();
            final JLabel guiResText = new javax.swing.JLabel();
            final JComboBox guiResBox = new javax.swing.JComboBox();
            final JCheckBox guiOverlapBox = new javax.swing.JCheckBox();
            final JCheckBox guiReverseBox = new javax.swing.JCheckBox();
            final JLabel guiOutputLabel = new javax.swing.JLabel();
            final JCheckBox guiRawBox = new javax.swing.JCheckBox();
            final JCheckBox guiPointBox = new javax.swing.JCheckBox();
            final JCheckBox guiRadiusBox = new javax.swing.JCheckBox();
            final JCheckBox guiHoughBox = new javax.swing.JCheckBox();
            final JButton guiOKButton = new javax.swing.JButton();
            
            //Format GUI text
            guiTitle.setFont(new java.awt.Font("Tahoma", 0, 18)); // NOI18N
            guiTitle.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            guiTitle.setText("Hough Circle Transform");
            guiIntro1.setText("This plugin performs a Hough circle transform on an image or stack. ");
            guiIntro2.setText("It can be used to find and measure circles embedded with an image.");
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
            guiCirNumLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiCirNumLabel.setText("Number of Circles (enter 0 to use threshold):");
            guiCirNumText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiCirNumText.setText("1");
            guiCirNumText.setToolTipText("The number of circles you expect to find in the image");
            guiThreshLabel.setHorizontalAlignment(javax.swing.SwingConstants.TRAILING);
            guiThreshLabel.setText("Hough score threshold (used if # of circles = 0):");
            guiThreshText.setHorizontalAlignment(javax.swing.JTextField.TRAILING);
            guiThreshText.setText("60");
            guiThreshText.setToolTipText("The Hough score threshold above which circles are valid");
            guiResText.setText("Hough transform resolution (# of steps in each transform)");
            guiResBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "4", "8", "16", "32", "64", "128", "256", "512", "1024" }));
            guiResBox.setSelectedIndex(2);
            guiOverlapBox.setText("Include overlapping circles in search");
            guiReverseBox.setText("Reverse search direction (large radius -> small radius)");
            guiOutputLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
            guiOutputLabel.setText("Output options:");
            guiRawBox.setText("Raw Hough transform series");
            guiPointBox.setText("Point selection(s) at circle centroid(s) overlaid on the original image");
            guiRadiusBox.setText("Map of circle centroids (pixel intensity = circle radius)");
            guiHoughBox.setText("Map of circle centroids (pixel intensity = Hough score)");
            guiOKButton.setText("OK");
            guiOKButton.addActionListener(new java.awt.event.ActionListener() {
            // </editor-fold>
                
            //When push button is pushed, retrieve the state of 
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                    // <editor-fold desc="Retrieve GUI arguments and calculate Hough parameters">
                    //Retrive the numbers from the text boxes and combobox
                    radiusMin = Integer.parseInt(guiMinText.getText());
                    radiusMax = Integer.parseInt(guiMaxText.getText());
                    radiusInc = Integer.parseInt(guiIncText.getText());
                    nCircles = Integer.parseInt(guiCirNumText.getText());
                    threshold = Integer.parseInt(guiThreshText.getText());
                    resolution = Integer.parseInt((String) guiResBox.getSelectedItem());
                    
                    //Retrieve the check box status
                    allowOverlap = guiOverlapBox.isSelected();
                    reverseSearch = guiReverseBox.isSelected();
                    houghSeries = guiRawBox.isSelected();
                    showCircles = guiPointBox.isSelected();
                    showRadius = guiRadiusBox.isSelected();
                    showScores = guiHoughBox.isSelected();
                    
                    //Remove the GUI frame now that it is no longer needed
                    guiFrame.dispose();
                    
                    //Start the Hough Transform
                    startTransform();
                    
                }    
            });
            // <editor-fold desc="Swing GUI Part 2">
            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(guiFrame.getContentPane());
            guiFrame.getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(guiTitle, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(layout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(guiMinLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(guiMinText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(guiMaxLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(guiMaxText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(guiIncLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(guiIncText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(guiCirNumLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(guiCirNumText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(guiThreshLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGroup(layout.createSequentialGroup()
                                    .addGap(0, 0, Short.MAX_VALUE)
                                    .addComponent(guiResText)))
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(guiResBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(guiThreshText, javax.swing.GroupLayout.PREFERRED_SIZE, 104, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGroup(layout.createSequentialGroup()
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(guiIntro1)
                                .addComponent(guiIntro2)
                                .addComponent(guiSearchLabel)
                                .addComponent(guiReverseBox)
                                .addComponent(guiOverlapBox)
                                .addComponent(guiHoughBox)
                                .addComponent(guiRadiusBox)
                                .addComponent(guiPointBox)
                                .addComponent(guiOutputLabel)
                                .addComponent(guiRawBox))
                            .addGap(0, 0, Short.MAX_VALUE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addGap(0, 0, Short.MAX_VALUE)
                            .addComponent(guiOKButton, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addContainerGap())
            );
            layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addComponent(guiTitle)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
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
                        .addComponent(guiCirNumText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(guiCirNumLabel))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(guiThreshText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(guiThreshLabel))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(guiResText)
                        .addComponent(guiResBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                    .addComponent(guiOverlapBox)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(guiReverseBox)
                    .addGap(18, 18, 18)
                    .addComponent(guiOutputLabel)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                    .addComponent(guiRawBox)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(guiPointBox)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(guiRadiusBox)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(guiHoughBox)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
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
    private void startTransform(){
        //Calculate Hough parameters
        depth = ((radiusMax-radiusMin)/radiusInc)+1;
        if (nCircles > 0) {
            useThreshold = false;
            threshold = -1;
        } else {
            useThreshold = true;
            if(threshold < 0) {
                IJ.showMessage("Threshold must be greater than 0");
                return;
            }
        }
        //If the reverse search is not selected then start at min, increment positively
        if(!reverseSearch){
            radiusBegin = radiusMin;
            radiusStep = radiusInc;
        }
        //If reverse search is selected then start at max, increment negatively
        else{
            radiusBegin = radiusMax;
            radiusStep = radiusInc * -1;
        }
        // </editor-fold>

        // <editor-fold desc="Send arguments to record">
        //If the macro was being recorded, return the set values to the recorder
        if (Recorder.record){
            String Command  = "run(\"Hough Circle Transform\",\"min=" + radiusMin + ", max=" + radiusMax + ", inc=" + radiusInc + 
                    ", nCircles=" + nCircles + ", threshold=" + threshold + ", resolution=" + resolution + "";
            if(allowOverlap) Command += " overlap";
            if(reverseSearch) Command += " reverse";
            if(houghSeries) Command += " show_raw";
            if(showCircles) Command += " show_select";
            if(showRadius) Command += " show_centroids";
            if(showScores) Command += " show_scores";
            Command += "\");\r\n";
            Recorder.recordString(Command);
        }
        // </editor-fold>

        //Import the ImagePlus as a stack
        ImageStack stack = imp.getStack();

        //Convert the stack to float (allows 8, 16, and 32 stacks to all be processed as one type)
        ImageStack stackCopy = stack.duplicate();
        ImageStack floatStack = stackCopy.convertToFloat();

        //Build the transform LUT (all necessary translations in Cartesian coordinates)
        lutSize = buildLookUpTable();


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

        //Retrieve the current slice in the stack as an ImageProcessor
        for(int slice = 1; slice <= imp.getStackSize(); slice++){
            ImageProcessor ip = floatStack.getProcessor(slice);
            imageValues = (float[]) ip.getPixels();

            //Calculate the Hough transform for the image
            houghTransform();

            if(houghSeries){
                //Create the hyperstach to put into the result if needed
                HoughSpaceSeries(slice, houghStack);
            }

            // Mark the center of the found circles in a new image if user wants to find centers
            if(showCircles || showRadius || showScores){
               if(useThreshold) getCenterPointsByThreshold(threshold);
               else getCenterPoints(nCircles);
            }


            // Create image View for Marked Circles.
            if(showCircles) drawCircles(slice, circleStack);
            
            //Create map of centroids where the intensity is the radius
            if (showRadius){
                centroidsip = new ByteProcessor(width, height);
                byte[] centroidpixels = (byte[])centroidsip.getPixels();
                drawCentroids(centroidpixels);
                new ImagePlus("Centroid Map contains " +nCircles+ " circles", centroidsip).show();
            }
            
            //Create a map of centroids where the intensity is the Hough score
            if(showScores){
                
            }
        }

        //Draw the resulting stacks
         if(houghSeries){
            houghPlus = new ImagePlus("Hough Transform Series", houghStack);
            houghPlus = HyperStackConverter.toHyperStack(houghPlus, 1, depth, imp.getNSlices(), "default", "grayscale");
            houghPlus.show();
         }
         if(showCircles){
             new ImagePlus(nCircles+" Circles Found", circleStack).show();
         }
         
               
    }
    //Build an array of the cartesion transformations necessary for each increment at each radius
    private int buildLookUpTable() {
        //Build an array to store the X and Y coordinates for each angle increment (resolution) across each radius (depth)
        lut = new int[2][resolution][depth];
        int radius = radiusBegin;
        
        //Initialize a variable that will allow for measuring the maximium LUT size
        int maxLUT = 0;
        
        
        //Step through all radii to be sampled
        while (radius >= radiusMin && radius <= radiusMax ) {
            
            //Index counter that also tracks the largest actual LUT size (may be <= resolution)
            int i = 0; 
            for(int resStep = 0; resStep < resolution; resStep++) {
                
                //Calcualte the angle and corresponding X and Y displacements for the specified radius
                double angle = (2*Math.PI * (double)resStep) / (double)resolution;
                int indexR = (radius-radiusMin)/radiusInc;
                int rcos = (int)Math.round ((double)radius * Math.cos (angle));
                int rsin = (int)Math.round ((double)radius * Math.sin (angle));
                
                //Test to make sure that the coordinate is a new coordinate
                //NOTE: A continuous circle is being discretized into pixels, it is possible that small angle steps for small radii circles will occupy the same pixel.
                //Since there is no point in making redundant calculations, these points are excluded from the LUT
                if((i == 0) || (rcos != lut[0][i][indexR]) && (rsin != lut[1][i][indexR])) {
                    lut[0][i][indexR] = rcos;
                    lut[1][i][indexR] = rsin;
                    i++;
                }
            }
            radius = radius + radiusStep; //Step the redius in the specified direction
            
            //Check to see if a new maximum LUT length was found
            if(i>maxLUT) maxLUT = i;
        }
        
        return maxLUT;
    }

    
    private void houghTransform () {
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(radiusMin);  
        
        //Create an array to store the Hough values
        houghValues = new double[width][height][depth];
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  
  
                @Override
                public void run() {  
  
                //Divide the radius tasks across the cores available
                for (int radius = ai.getAndAdd(radiusInc); radius <= radiusMax; radius = ai.getAndAdd(radiusInc)) {  
                    
                    //For a given radius, transform each pixel in a circle, and add-up the votes 
                    for(int y = 1; y < height-1; y++) {
                        for(int x = 1; x < width-1; x++) {
                                if( imageValues[(x+offx)+(y+offy)*offset] != 0 )  {// Edge pixel found
                                    int indexR=(radius-radiusMin)/radiusInc;
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
    

    //Find the largest Hough pixel in the 3D Hough transform array to scale the 8-bit conversion
    private void houghMaximum () {
        //long startTime = System.currentTimeMillis(); //1337ms without multi, 319ms with multi, 175ms by writing to private variable per thread
        //Build an array to store the result from each thread
        final Thread[] threads = newThreadArray();
        
        //Build an array to store the max from each thread
        maxArray = new double[threads.length];
        
        //Create an atomic integer counter that each thread can use to count through the radii
        final AtomicInteger ai = new AtomicInteger(0);
        
        //Create an integer for indexing the results array (one result per thread
        final AtomicInteger Index = new AtomicInteger(0);
        maxPixel = 0D;
        
        //Build a thread for as many CPUs as are available to the JVM 
        for (ithread = 0; ithread < threads.length; ithread++) {    
            
            // Concurrently run in as many threads as CPUs  
            threads[ithread] = new Thread() {  
                          
                { setPriority(Thread.NORM_PRIORITY); }  

                //Search for the largest score with each thread
                @Override
                public void run() {
                    double maxThread = 0D;
                    for(int a=ai.getAndIncrement(); a<depth; a=ai.getAndIncrement()){
                        for(int j = 0; j < height; j++) {
                            for(int k = 0; k < width; k++){
                                if(houghValues[k][j][a] > maxThread) {
                                    maxThread = houghValues[k][j][a];                                   
                                }
                            }
                        }
                    }
                    //Have each thread report the score to a common array
                    maxArray[Index.getAndIncrement()] = maxThread;
                }
            };
        }
        startAndJoin(threads);
        
        //Search common array for highest score
        for(int a = 0; a<maxArray.length; a++){
            if(maxArray[a]>maxPixel) maxPixel = maxArray[a];
        }
        
        //long stopTime = System.currentTimeMillis();
        //IJ.log("Elapsed time was " + (stopTime - startTime) + " miliseconds.");
    }
    
    //Add transform series data to the hyperstack
    private void HoughSpaceSeries(int slice, ImageStack houghStack){
    
        //Find the brightest pixel in the whole Hough transform array to normalize stack intensity
        houghMaximum();
        
        int startFrame = ((slice-1)*depth);
        
        int radius = radiusBegin;
        while (radius >= radiusMin && radius <= radiusMax ) {
            //Create a new image for depositing the Hough pixels into
            ImageProcessor newip = new ByteProcessor(width, height);
            byte[] newpixels = (byte[])newip.getPixels();

            //Calculate the corresponding index
            houghIndex = (radius-radiusMin)/radiusInc;
                
            //Retrieve the pixel array for the current Hough radius image
            createHoughPixels(newpixels, houghIndex);
                
            //Deposit the array image into the corresponding slice in the stack
            houghStack.setPixels(newpixels, houghIndex+1+startFrame);
                
            //Give the current slice the appropriate radius label
            houghStack.setSliceLabel("Hough Space [r="+radius+"]", houghIndex+1+startFrame);
            
            radius = radius + radiusStep;
        }
    }
    
    // Convert Values in Hough Space to an 8-Bit Image Space.
    private void createHoughPixels (byte houghPixels[], int index) {
	//Rescale all the Hough values to 8-bit to create the Hough image - 47ms to complete - single threading okay
        for(int l = 0; l < height; l++) {
            for(int i = 0; i < width; i++) {
                houghPixels[i + l * width] = (byte) Math.round ((houghValues[i][l][index] * 255D) / maxPixel);
            }

        }
    }

    // Draw the circles found in the original image.
    private void drawCircles(int slice, ImageStack circleStack) {
		
            // Copy original input pixels into output
            // circle location display image and
            // combine with saturation at 100
            circlesip = new ByteProcessor(width, height);
            byte[] circlespixels = (byte[])circlesip.getPixels();

            int roiaddr=0;
            for( int y = offy; y < offy+height; y++) {
                    for(int x = offx; x < offx+width; x++) {
                            // Copy;
                            circlespixels[roiaddr] = (byte) imageValues[x+offset*y];
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
            if(centerPoint == null) {
                    if(useThreshold)
                        getCenterPointsByThreshold(threshold);
                    else
                        getCenterPoints(nCircles);
            }
            byte cor = -1;
            // Redefine these so refer to ROI coordinates exclusively
            offset = width;
            offx=0;
            offy=0;

            for(int l = 0; l < nCircles; l++) {
                    int i = centerPoint[l].x;
                    int j = centerPoint[l].y;
                    // Draw a gray cross marking the center of each circle.
                    for( int k = -10 ; k <= 10 ; ++k ) {
                            int p = (j+k+offy)*offset + (i+offx);
                            if(!outOfBounds(j+k+offy,i+offx))
                                    circlespixels[(j+k+offy)*offset + (i+offx)] = cor;
                            if(!outOfBounds(j+offy,i+k+offx))
                                    circlespixels[(j+offy)*offset   + (i+k+offx)] = cor;
                    }
                    for( int k = -2 ; k <= 2 ; ++k ) {
                            if(!outOfBounds(j-2+offy,i+k+offx))
                                    circlespixels[(j-2+offy)*offset + (i+k+offx)] = cor;
                            if(!outOfBounds(j+2+offy,i+k+offx))
                                    circlespixels[(j+2+offy)*offset + (i+k+offx)] = cor;
                            if(!outOfBounds(j+k+offy,i-2+offx))
                                    circlespixels[(j+k+offy)*offset + (i-2+offx)] = cor;
                            if(!outOfBounds(j+k+offy,i+2+offx))
                                    circlespixels[(j+k+offy)*offset + (i+2+offx)] = cor;
                    }
            }
            circleStack.setPixels(circlespixels, slice);
    }
    
    // Draw the centroids found in the original image where intensity = radius.
    public void drawCentroids(byte[] centroidpixels) {
		
		for(int l = 0; l < nCircles; l++) {
			int i = centerPoint[l].x;
			int j = centerPoint[l].y;
                        byte radius = (byte) centerRadii[l];
                        			
                        //Draw a point at the centroid and set the int to = radius
                        centroidpixels[j*offset + i] = radius;
                        
                        //Record the Hough score on the next pixel below the centroid if desired
                        if(showScores){
                            centroidpixels[j*offset + i + 1] = (byte) houghScores[l];
                        }
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

    public Point nthMaxCenter (int i) {
        return centerPoint[i];
    }


    /** Search for a fixed number of circles.
    @param nCircles The number of circles that should be found.  
    */
    private void getCenterPoints (int nCircles) {
        long startTime = System.currentTimeMillis();
        houghScores = new int [nCircles];
        centerRadii = new int [nCircles];
        centerPoint = new Point[nCircles];
        int xMax = 0;
        int yMax = 0;
        int rMax = 0;




        for(int c = 0; c < nCircles; c++) {
            double counterMax = -1;
            int radius = radiusBegin;
            while (radius >= radiusMin && radius <= radiusMax ) {


                int indexR = (radius-radiusMin)/radiusInc;
                for(int y = 0; y < height; y++) {
                    for(int x = 0; x < width; x++) {
                        if(houghValues[x][y][indexR] > counterMax) {
                            counterMax = houghValues[x][y][indexR];
                            xMax = x;
                            yMax = y;
                            rMax = radius;
                        }
                    }

                }
                radius = radius + radiusStep;
            }

            centerPoint[c] = new Point (xMax, yMax);
            centerRadii[c] = rMax;
            houghScores[c] = (int) counterMax;
            

            clearNeighbours(xMax,yMax,rMax);
        }
        long stopTime = System.currentTimeMillis();
        IJ.log("Elapsed time was " + (stopTime - startTime) + " miliseconds.");
    }


    /** Search circles having values in the hough space higher than a threshold
    @param threshold The threshold used to select the higher point of Hough Space
    */
    private void getCenterPointsByThreshold (int threshold) {
        
        houghScores = new int [VECTOR_SIZE_MAX];
        centerRadii = new int [VECTOR_SIZE_MAX];
        centerPoint = new Point[VECTOR_SIZE_MAX];
        int xMax = 0;
        int yMax = 0;
        int countCircles = 0;

        int radius = radiusBegin;
        while (radius >= radiusMin && radius <= radiusMax ) {
            int indexR = (radius-radiusMin)/radiusInc;
            for(int y = 0; y < height; y++) {
                for(int x = 0; x < width; x++) {



                    if(houghValues[x][y][indexR] > threshold) {


                        if(countCircles < VECTOR_SIZE_MAX) {


                            centerPoint[countCircles] = new Point (x, y);
                            centerRadii[countCircles] = radius;
                            houghScores[countCircles] = (int) houghValues[x][y][indexR];

                            clearNeighbours(x,y,radius);

                            ++countCircles;
                        } else
                            break;
                    }
                }
            }
            radius = radius + radiusStep;
        }

        nCircles = countCircles;
    }

    /** Clear, from the Hough Space, all the counter that are near (radius/2) a previously found circle C.
        
    @param x The x coordinate of the circle C found.
    @param x The y coordinate of the circle C found.
    @param x The radius of the circle C found.
    */
    private void clearNeighbours(int x,int y, int radius) {


        // The following code just clean the points around the center of the circle found.
	double radiusSquared = radius*radius;


        int y1 = (int)Math.floor ((double)y - radius);
        int y2 = (int)Math.ceil ((double)y + radius) + 1;
        int x1 = (int)Math.floor ((double)x - radius);
        int x2 = (int)Math.ceil ((double)x + radius) + 1;



        if(y1 < 0)
            y1 = 0;
        if(y2 > height)
            y2 = height;
        if(x1 < 0)
            x1 = 0;
        if(x2 > width)
            x2 = width;



        int r = radius;
        while (r >= radiusMin && r <= radiusMax )  {
            int indexR = (r-radiusMin)/radiusInc;
            for(int i = y1; i < y2; i++) {
                for(int j = x1; j < x2; j++) {	      	     
                    if(Math.pow (j - x, 2D) + Math.pow (i - y, 2D) < radiusSquared) {
                        houghValues[j][i][indexR] = 0.0D;
                    }
                }
            }
            r = r + radiusStep;
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

 