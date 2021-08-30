import ij.*;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.gui.GenericDialog;
import java.awt.image.ColorModel;
//import ij.process.ImageProcessor;
import java.lang.Math;
import java.lang.Exception;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.GaussianBlur;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.*;
import java.util.*;

//import AutoThresholder;




public class LocalOrientation_orientationJ_20140803 implements PlugIn {

	ImagePlus experiment; // input stack of images
	ImageStack orientation_map; // stack of the orientation colour map
	ImageStack orientation_lines; // stack of the orientation lines at certain spacing
	int spacing_orientlines = 20; // spacing between the pts at which the orientation lines are plotted
	private String outtext="/home/mukund/Desktop/actinOP.txt"; // file to which order parameter values are written
	public boolean calculate_cylindrical_op = false;
	public boolean output_orientation_colourmap = true;
	public boolean output_orientation_directionlines = false;
	public float energy_threshold = 0.005f;
	private float coherency_threshold = 0.02f;

	int cylinder_centre_x = 720;
	int cylinder_centre_y = 562;
	float cylinder_radius = 30; 

	//********************************* run() function ********************************
	
	public void run(String arg) { 
		boolean consistent;
			
		do {
			consistent=true;
			if (!showDialog()) return;
			else {
				// obtain input image parameters
				int width=experiment.getWidth();
				int height=experiment.getHeight();
				int frames=experiment.getNSlices();

				Roi pr0 = experiment.getRoi();

				ImageStack imstack_input=experiment.getImageStack();
				if (output_orientation_colourmap) {orientation_map = new ImageStack(width, height);}
				if (output_orientation_directionlines) {orientation_lines = new ImageStack(width, height);}

				// calculating the order parameter for all the frames
				float [] order_parameter = new float [frames];
				for(int f=0;f<frames;f++) {
					ImageProcessor ip = imstack_input.getProcessor(f+1); // extracting the imageprocessor of the frame f
					FloatProcessor fip = new FloatProcessor(width, height); // floatprocessor to which ip would be converted
					fip = ip.toFloat(1,fip); // converting ip to float processor for local_gradient_orientation()
					// calling the local orient fn to get order parameter value and write output images
					order_parameter[f] = local_gradient_orientation(fip, pr0); 
					//IJ.log(""+order_parameter[f]);
					IJ.showProgress((double)f/(double)frames);
				}
				
				// showing output images
				if (output_orientation_colourmap) {
					ImagePlus orient_map = new ImagePlus("Orientation map ", orientation_map);
					orient_map.show();
				}
				if (output_orientation_directionlines) {
					ImagePlus orient_lines = new ImagePlus("Orientation field ", orientation_lines);
					orient_lines.show();
				}
				
						
				// writing order parameter to the output file 

				try {
					FileWriter pout=new FileWriter(outtext);
					PrintWriter printout=new PrintWriter(pout);
					for(int f=0; f<frames; f++) {
						printout.print(order_parameter[f]+"\n");	
					}
					printout.flush();
					printout.close();
					pout.close();
				}
				catch(FileNotFoundException fe)	{IJ.log("FileNotFoundException - writing failed : " + fe);}
				catch(IOException ioe){	IJ.log("IOException - writing failed  : " + ioe);}

			}		
		} while (!consistent);
		
		IJ.register(LocalOrientation_orientationJ_20140803.class);
	
	} // end of run()

	
	
	//***************************** Local orientation function *********************************************
	
	// this function returns the order parameter for the given frame
	// along with printing out the image with local orientation map

	public float local_gradient_orientation (FloatProcessor ip_orient, Roi pr0){
		
		int width = ip_orient.getWidth();
		int height = ip_orient.getHeight();
		
		FloatProcessor grad_x = (FloatProcessor)ip_orient.duplicate(); // initialized x derivative of image
		FloatProcessor grad_y = (FloatProcessor)ip_orient.duplicate(); // initialized y derivative of image
	
		// creating the DoG (Derivative of Gaussian) kernels

		float DoG_sigma = 1.0f;
		int DoG_size = (int)(3*DoG_sigma);//9

		
		final float[] DoG_x = new float [DoG_size*DoG_size];
		for (int x=0; x<DoG_size; x++) {
			for (int y=0; y<DoG_size; y++) {
				DoG_x[x+y*DoG_size] =  (x - DoG_size/2)*(float)Math.exp(-((x - DoG_size/2)*(x - DoG_size/2) + (y - DoG_size/2)*(y - DoG_size/2))/(2*DoG_sigma*DoG_sigma))/(float)(4*DoG_sigma*DoG_sigma*DoG_sigma*DoG_sigma*Math.PI);
			}
		}

		final float[] DoG_y = new float [DoG_size*DoG_size];
		for (int x=0; x<DoG_size; x++) {
			for (int y=0; y<DoG_size; y++) {
				DoG_y[x+y*DoG_size] =  (y - DoG_size/2)*(float)Math.exp(-((x - DoG_size/2)*(x - DoG_size/2) + (y - DoG_size/2)*(y - DoG_size/2))/(2*DoG_sigma*DoG_sigma))/(float)(4*DoG_sigma*DoG_sigma*DoG_sigma*DoG_sigma*Math.PI);
			}
		}

		// convolution with DoG filters to get derivatives in x and y direction
		grad_x.convolve(DoG_y, DoG_size, DoG_size); 
		grad_y.convolve(DoG_x, DoG_size, DoG_size); 
		
				
		//-------------------- calculating the elements of the structure tensor-----------------
		//---------------------------[ Sxx     Sxy ]---------------
		//---------------------------[ Sxy     Syy ]---------------
		
		// square of the derivatives
		FloatProcessor Ixx = (FloatProcessor)ip_orient.duplicate();
		FloatProcessor Iyy = (FloatProcessor)ip_orient.duplicate();
		FloatProcessor Ixy = (FloatProcessor)ip_orient.duplicate();
		
		for (int x=0;x<width;x++) {
			for (int y=0;y<height;y++) {
				Ixx.setf(x,y,grad_x.getf(x,y)*grad_x.getf(x,y));
				Iyy.setf(x,y,grad_y.getf(x,y)*grad_y.getf(x,y));
				Ixy.setf(x,y,grad_x.getf(x,y)*grad_y.getf(x,y));
			}
		}
		
		// the elements of the structure tensor
		FloatProcessor Sxx = (FloatProcessor)Ixx.duplicate();
		FloatProcessor Syy = (FloatProcessor)Iyy.duplicate();
		FloatProcessor Sxy = (FloatProcessor)Ixy.duplicate();

		float gauss_sigma = 1.0f;
		int gauss_size = (int)(7*gauss_sigma); //7
		final float[] kernel_gauss = new float [gauss_size*gauss_size];
		for (int x=0; x<gauss_size; x++) {
			for (int y=0; y<gauss_size; y++) {
				kernel_gauss[x+y*gauss_size] = (float)Math.exp(-((x - gauss_size/2)*(x - gauss_size/2) + (y - gauss_size/2)*(y - gauss_size/2))/(2*gauss_sigma*gauss_sigma))/(float)(2*Math.PI*gauss_sigma*gauss_sigma);
			}
		}
		Sxx.convolve(kernel_gauss, gauss_size, gauss_size);
		Sxy.convolve(kernel_gauss, gauss_size, gauss_size);
		Syy.convolve(kernel_gauss, gauss_size, gauss_size);

		// convolving with a gaussian of sigma 1
		//gaussianBlur.blurGaussian(Sxx,1.0,1.0,0.01);
		//gaussianBlur.blurGaussian(Syy,1.0,1.0,0.01);
		//gaussianBlur.blurGaussian(Sxy,1.0,1.0,0.01);

		//--------------- calculating the local theta, coherence, and energy---------------------
		
		// theta = 0.5*atan[(2Sxy/(Syy-Sxx)]
		// coherence = [(Sxx-Syy)/(Sxx+Syy)]^2
		// energy = Sxx + Syy
		
		FloatProcessor theta_image = new FloatProcessor(width, height); // image with local orientation in grayscale values
		FloatProcessor ortho_image = new FloatProcessor(width, height); // image for ortho coordinates
		FloatProcessor theta = new FloatProcessor(width, height); // local angles
		FloatProcessor coherency = new FloatProcessor(width, height); //  local coherence values
		FloatProcessor energy = new FloatProcessor(width, height); // local energy values
		float energy_max=0.0f;
		for (int x=0;x<width;x++) {
			for (int y=0;y<height;y++) {
				float angle_structtensor = (float)(0.5f*Math.atan2(2.0*(double)Sxy.getf(x,y),(double)(Syy.getf(x,y)-Sxx.getf(x,y)))); // radians, -PI/2 to + PI/2 ?
				float angle_grad = 0f;
				//if (grad_x.getf(x,y)!=0f) { angle_grad = (float)(Math.atan2(grad_y.getf(x,y)/grad_x.getf(x,y)));} // radians, -pi to +pi
				//else angle_grad = (float)Math.PI/2;
			
				theta.setf(x,y,angle_structtensor); // radians, -pi/2 to +pi/2 ?
				
				coherency.setf(x,y,(((Sxx.getf(x,y)-Syy.getf(x,y))*(Sxx.getf(x,y)-Syy.getf(x,y)))/((Sxx.getf(x,y)+Syy.getf(x,y))*(Sxx.getf(x,y)+Syy.getf(x,y)))));
				
				energy.setf(x,y,Sxx.getf(x,y) + Syy.getf(x,y));
				if (energy.getf(x,y)>energy_max) {energy_max=energy.getf(x,y);}
			}
		}

		// averaging the local theta values with a gaussian filter
		//theta.convolve(kernel_gauss, gauss_size, gauss_size);
		
		
		
		//---------------------calculating the Order parameter---------------------------------------
		//----------------------- = <Cos2(theta_mean - theta)> --------------------------------------
		
		Roi ppr = (PolygonRoi)pr0;
		
		// calculating the avg angle
		float num_pix = 0.0f;
		float sum_angles = 0.0f; // radians
		for (int x=0;x<width;x++) {
			for (int y=0;y<height;y++) {
				if ( ppr.contains(x,y) && energy.getf(x,y)/energy_max>energy_threshold){
					sum_angles = sum_angles + theta.getf(x,y);
					num_pix++;
				}
			}
		}
		float avg_angle = sum_angles/num_pix; // radians, -pi/2 to +pi/2
		//IJ.log(avg_angle+"");
		
		// calculating the order parameter
		
		float order_parameter_sum = 0; // sum of order parameter for all pixels 
		float num_pixel_op = 0;
		for (int x=0;x<width;x++) {
			for (int y=0;y<height;y++) {
				float angle_difference = avg_angle-theta.getf(x,y); // radians, 
				float pixel_op = (float)Math.cos(2*angle_difference); // -1 to +1
				if ( ppr.contains(x,y) && energy.getf(x,y)/energy_max>energy_threshold) {
					order_parameter_sum = order_parameter_sum + pixel_op;
					num_pixel_op++;
				}
			}
		}
		float order_parameter = Math.abs(order_parameter_sum/num_pixel_op); // -1 to +1 avg. order parameter

		
		//------------------------calculating the order parameter in cylindrical coordinates---------------------------

		if (calculate_cylindrical_op) {
			// find the centroid of the cell. use thresholding and centre of mass for it, or get user input for the position of the nucleus
			
			float[] op_cylinder_r = new float[(int)cylinder_radius];
			for (int r=0; r<cylinder_radius; r++) { 
				float op_sum = 0;
				float num_op_sum = 0;
				for (int angle=0; angle<180; angle++) {
					int x = cylinder_centre_x + (int)(r*Math.cos(angle));
					int y = cylinder_centre_y + (int)(r*Math.sin(angle));
					float angle_difference = angle - theta.getf(x,y); 
					if ( ppr.contains(x,y) && energy.getf(x,y)/energy_max>energy_threshold) {
						float pixel_op_rad = -(float)Math.cos(2*angle_difference); // -1 to +1
						op_sum = op_sum + pixel_op_rad;
						num_op_sum++;
					}
				}
				op_cylinder_r[r] = op_sum/num_op_sum;
				//IJ.log(op_cylinder_r[r]+" ");
			}
		}


		
		
		//------------------------output of images ------------------------------------
		
		
		// output image with local theta values
		//FloatProcessor theta_image = new FloatProcessor(width,height);
		for (int x=0;x<width;x++) {
			for (int y=0;y<height;y++) {
				if ( ppr.contains(x,y) && energy.getf(x,y)/energy_max>energy_threshold) {
					theta_image.setf(x,y,(int)((theta.getf(x,y))/Math.PI*180.0)); // converted theta to degrees, 0 to 180
				}
			else {theta_image.setf(x,y,0);}
			}
		}
		//ImagePlus orient_grayscale = new ImagePlus("orientation map grayscale",theta_image);
		//orient_grayscale.show();

		// output image with 2 colors for orthoradial coordinates
		if (calculate_cylindrical_op) {
			for (int x=0;x<width;x++) {
				for (int y=0;y<height;y++) {
					float angle_radial_vector = (float)Math.atan2((double)(y-cylinder_centre_y),(double)(x-cylinder_centre_x)); // angle of the radial vector from cylinder centre to (x,y)
					float angle_difference = angle_radial_vector-theta.getf(x,y); // radians, 
					if ( ppr.contains(x,y) && energy.getf(x,y)/energy_max>energy_threshold) {
						float pixel_op_rad = (float)Math.cos(2*angle_difference); // -1 to +1
						if (pixel_op_rad>0.3) {ortho_image.setf(x,y,255);}
						else if (pixel_op_rad<-0.3){ortho_image.setf(x,y,100);}
					}
					else {ortho_image.setf(x,y,0);}	
				}
			}
			ImagePlus orient_ortho = new ImagePlus("orientation map orthoradial",ortho_image);
			orient_ortho.show();
		}
		
		// image showing the local orientation as colors
		if (output_orientation_colourmap) {
			ColorProcessor cp = new ColorProcessor(width, height);
			byte[] H = (byte[]) theta.convertToByte(true).getPixels();
			byte[] S = (byte[]) theta.convertToByte(true).getPixels();
			byte[] B = (byte[]) theta.convertToByte(true).getPixels();
			for (int x = 0; x < width; x++) {
				for (int y = 0; y < height; y++) {
					H[x+y*width] = (byte) (255.0 * (theta.getf(x,y)+Math.PI/2)/Math.PI); // theta 0 to pi, and then normalized
					S[x+y*width] = (byte) (255.0);// * coherency[i] ) ; 
					if ((energy.getf(x,y)/energy_max<energy_threshold || coherency.getf(x,y)<coherency_threshold)) {
						B[x+y*width] = (byte) (0.0);
					}
					else B[x+y*width] = (byte) (255.0);
				}
			}		
			cp.setHSB(H, S, B);
			orientation_map.addSlice("sdfd",cp);
		}
		
		
		// image showing the local orientations as lines at certain spacing
		if (output_orientation_directionlines) {
			FloatProcessor ol = new FloatProcessor(width, height);
			for (int x = spacing_orientlines; x<width;x=x+spacing_orientlines) {
				for (int y = spacing_orientlines; y<height; y=y+spacing_orientlines) {
					float angle = theta.getf(x,y)+(float)Math.PI/2; // 0 to pi
					int x1 = x - (int)(0.75*spacing_orientlines*Math.cos(angle));
					int y1 = y - (int)(0.75*spacing_orientlines*Math.sin(angle)); 
					int x2 = x + (int)(0.75*spacing_orientlines*Math.cos(angle)); 
					int y2 = y + (int)(0.75*spacing_orientlines*Math.sin(angle)); 
					//if ((energy.getf(x,y)/energy_max>energy_threshold && coherency.getf(x,y)>coherency_threshold)){
						ol.drawLine(x1,y,x2,y2);
					//}
				}
			}
			orientation_lines.addSlice("",ol);
		}
		
		
		//------------ return the order parameter value for the frame-----------------------
		return order_parameter;
	
	}// end of local_gradient_orientation (FloatProcessor ip_orient)

	
	// ************************ function to show dialog box to get user input ****************************
	
	private boolean showDialog() {
		int[] wList = WindowManager.getIDList();
		if (wList==null || wList.length<1) {
			IJ.showMessage("Actin Local Orientation", "image required");
			return false;
		}
		String[] titles = new String[wList.length];
		for (int i=0; i<wList.length; i++) {
			ImagePlus imp = WindowManager.getImage(wList[i]);
			titles[i] = imp!=null?imp.getTitle():"";
		}
		
		GenericDialog gd = new GenericDialog("Actin Local Orientation");
		gd.addMessage("image");
		gd.addChoice("image or time series for cell", titles, titles[0]);
		gd.addMessage("output file");
		gd.addStringField("path and name of text file", outtext,40);
		gd.addCheckbox("calculate order parameter in cylindrical coordinates", calculate_cylindrical_op);
		gd.addMessage("output images");
		gd.addCheckbox("Local orientation Colour Map", output_orientation_colourmap);
		gd.addCheckbox("Local orintation direction lines", output_orientation_directionlines);
		gd.addNumericField("Threshold", energy_threshold,4);

		gd.showDialog();
		if (gd.wasCanceled()) return false;

		int index1 = gd.getNextChoiceIndex();
		experiment = WindowManager.getImage(wList[index1]);
		outtext=gd.getNextString();
		calculate_cylindrical_op = gd.getNextBoolean();
		output_orientation_colourmap = gd.getNextBoolean();
		output_orientation_directionlines = gd.getNextBoolean();
		energy_threshold = (float)gd.getNextNumber();
		
		
		return true;
	
	} // end of showDialog()


}// end of LocalOrientation_20130826
