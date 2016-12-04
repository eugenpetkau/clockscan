/** 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

import ij.*;
import ij.plugin.frame.RoiManager;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.measure.Calibration;

/**
 * This Plugin is used to measure the intensity of the brightness of the image
 * of the radial scan profiles in the region of interest (ROI). Reference:
 * <p>
 * "Clock-scan" protocol for image analysis Maxim Dobretsov, Dmitry Romanovsky
 * American Journal of Physiology - Cell Physiology Published 11 October 2006
 * Vol. 291 no. 5, C869-C879 DOI: 10.1152/ajpcell.00182.2006
 * 
 * @author Eugen Petkau
 */

public class Clock_Scan implements PlugIn {

	int widthInitial, heightInitial, widthTransform, heightTransform;
	double centerX, centerY;
	boolean isColor = false;
	ImageProcessor ip, ip2;
	ImagePlus imp, imp2;
	ImageProcessor ipTransform, ipInitial;
	ImagePlus iTransform, iInitial;
	boolean canceled = false;// aa
	private double X0;// a
	private double Y0;// a
	private FloatPolygon fp;
	private String btitle = ""; // need to file name
	private boolean realrad = false;// a
	private boolean subbackground = false;// a
	private boolean polar = false;// a
	private boolean paction = false;// a
	private float min = 255;
	private double limits = 1.2;// aa
	String center = "Center";
	static boolean clockWise = true; // true
	int[] rgbArray = new int[3];
	int[] xLyL = new int[3];
	int[] xLyH = new int[3];
	int[] xHyL = new int[3];
	int[] xHyH = new int[3];


	static boolean useCalibration = true;

	public void run(String arg) {

		ImagePlus imp = IJ.getImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}

		imp = WindowManager.getCurrentImage();

		Roi roi = imp.getRoi();
		RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();
		manager.addRoi(roi);
		manager.reset();
		

		// nur wenn eine selektion da ist
		if (roi != null) {
			doSetInterpolateRoi();
			setXYcenter();
		}

		ip = imp.getProcessor();
		if (roi == null)
			return;
		 doDialog() ;
		 if (!canceled){
		 doPicture(ip); 
		 }
		 WindowManager.setTempCurrentImage(imp);
		 
		 

	}

	public void doPicture(ImageProcessor ip) {

		imp = WindowManager.getCurrentImage();

		btitle = imp.getTitle();
		
		
		
		RoiManager rm = RoiManager.getInstance2();

		if (rm == null) {

			rm = new RoiManager();

		}

		Roi roi = imp.getRoi();
		RoiManager manager = RoiManager.getInstance();
		manager.reset();

		if (roi != null) {
			if (roi.getName() != "interpolate") {
				interpolate();
				roi = imp.getRoi();
				roi.setName("ROI interpolated");
				manager.addRoi(roi);
			}
		}


		setXYcenter();
		if (roi == null)
			return;

		if (canceled)
			return;
		doPaintCenter(ip2, X0, Y0);
		

		if (roi != null) {
			scaletocenter(roi, limits, limits, X0, Y0);
			Roi roi3 = imp.getRoi();
			roi3.setName("ROI with limit " + limits);
			manager.addRoi(roi3);
			// manager.runCommand(imp, "Measure");
		}

		doRadialClockScan(ip,roi);
		WindowManager.setTempCurrentImage(imp);

	}

	// Radial Clock Scan
	private void doRadialClockScan(ImageProcessor ip,Roi roi) {

		fp = imp.getRoi().getFloatPolygon();
		
		int xzahl = 0;
		int xmenge = 0;
		for (int i = 0; i < fp.npoints; i++) {
			double xmin = fp.getBounds().getCenterX() - 2;
			double xmax = fp.getBounds().getCenterX() + 2;
			double ycenter = fp.getBounds().getCenterY();
			double xx = (fp.xpoints[i]);
			double yy = (fp.ypoints[i]);
			if (yy > ycenter) {
				if (xmin <= xx && xx < xmax) {
					xzahl = i;
					xmenge = xmenge + 1;
				}
			}
		}
		int rechnung = xzahl + xmenge / 2;

		int entscheidung = 0;
		if (fp.xpoints[xzahl - xmenge] < fp.xpoints[xzahl + xmenge]) {
			entscheidung = 1;
		}
		
		
		double realradius = fp.npoints / (2 * Math.PI);//
		String pluntertitel = "scan length, ";
		if (!(realrad)) {
			realradius = 100;
			pluntertitel = "scan length, % of radius";
		}

		int nBins = (int) (realradius * limits); //

		float[][] Accpixel = new float[fp.npoints][nBins];
		float[][] Accpixelbild = new float[fp.npoints][nBins];
		float[][] Accumulator = new float[2][nBins];
				
		for (int i = rechnung; i < fp.npoints; i++) {

			double dx = (fp.xpoints[i] - X0);
			double dy = (fp.ypoints[i] - Y0);
			double radius = Math.sqrt(dx * dx + dy * dy);
			double sinA = dy / radius;
			double cosA = dx / radius;
			for (int j = 0; j < nBins; j++) {
				double newX = X0 + cosA * j * radius / nBins;
				double newY = Y0 + sinA * j * radius / nBins;
				Accpixelbild[i - rechnung][j] = ip.getPixel((int) newX, (int) newY);
				Accpixel[i - rechnung][j] = ip.getPixelValue((int) newX, (int) newY);
				Accumulator[0][j] = Accumulator[0][j] + 1;
				Accumulator[1][j] = Accumulator[1][j]
						+ ip.getPixelValue((int) newX, (int) newY);

			}
		}
		
		int rechnnerest = fp.npoints - rechnung;
		
		for (int i = 0; i < rechnung; i++) {

			double dx = (fp.xpoints[i] - X0);
			double dy = (fp.ypoints[i] - Y0);
			double radius = Math.sqrt(dx * dx + dy * dy);
			double sinA = dy / radius;
			double cosA = dx / radius;
			for (int j = 0; j < nBins; j++) {
				double newX = X0 + cosA * j * radius / nBins;
				double newY = Y0 + sinA * j * radius / nBins;
				Accpixelbild[i + rechnnerest][j] = ip.getPixel((int) newX, (int) newY);
				Accpixel[i + rechnnerest][j] = ip.getPixelValue((int) newX, (int) newY);
				Accumulator[0][j] = Accumulator[0][j] + 1;
				Accumulator[1][j] = Accumulator[1][j]+ ip.getPixelValue((int) newX, (int) newY);

			}
		}

		// real radius
		if (realrad) {

			// ///////////////////////////////////////////
			//IJ.log("you have real radius  ");
			if (polar) {
				Calibration cal = imp.getCalibration();
				
				ip = ip.createProcessor(fp.npoints, nBins);

				ip.setFloatArray(Accpixelbild);
				if (entscheidung == 1) {
					ip.flipHorizontal();
				}

				// cal
				imp2 = new ImagePlus(btitle + " linear transform", ip);
				imp2.setCalibration(cal);

				// imp2.show();
				ipInitial = imp2.getProcessor();
				int dept = imp2.getType();

				widthInitial = ipInitial.getWidth();
				heightInitial = ipInitial.getHeight();

				if (ipInitial instanceof ColorProcessor)
					isColor = true;
				cartesianto360();
				ipTransform.setMinAndMax(ipInitial.getMin(), ipInitial.getMax());
				ipTransform.setCalibrationTable(ipInitial.getCalibrationTable());
				iTransform = new ImagePlus(btitle + " polar transform ",ipTransform);
				iTransform.setCalibration(imp2.getCalibration());

				if (dept == 0) {
					ImageConverter ic = new ImageConverter(iTransform);

					ic.convertToGray8();
					iTransform.updateAndDraw();
				}

				if (dept == 3) {
					ImageConverter ic = new ImageConverter(iTransform);
					ic.convertToGray8();
					iTransform.updateAndDraw();
				}

				iTransform.show();

			}

			WindowManager.setTempCurrentImage(imp);

		}
		if (!realrad) {

			if (polar) {
				//IJ.log("you have radius 100% *limits ");
				ip = ip.createProcessor(fp.npoints, nBins);

				ip.setFloatArray(Accpixelbild);
				if (entscheidung == 1) {
					ip.flipHorizontal();
				}

				imp2 = new ImagePlus(btitle + " linear transform", ip);
				IJ.run(imp2, "Set Scale...","distance=1 known=1 pixel=1 unit=pixel");

				ipInitial = imp2.getProcessor();
				int dept = imp2.getType();

				widthInitial = ipInitial.getWidth();
				heightInitial = ipInitial.getHeight();

				if (ipInitial instanceof ColorProcessor)
					isColor = true;
				cartesianto360();
				ipTransform.setMinAndMax(ipInitial.getMin(), ipInitial.getMax());
				ipTransform.setCalibrationTable(ipInitial.getCalibrationTable());
				iTransform = new ImagePlus(btitle + " polar transform ",ipTransform);
				iTransform.setCalibration(imp2.getCalibration());

				if (dept == 0) {
					ImageConverter ic = new ImageConverter(iTransform);

					ic.convertToGray8();
					iTransform.updateAndDraw();
				}

				if (dept == 3) {
					ImageConverter ic = new ImageConverter(iTransform);
					ic.convertToGray8();
					iTransform.updateAndDraw();
				}

				iTransform.show();

			}

			WindowManager.setTempCurrentImage(imp);

		}

		// //////////////////////////////////////////////////////
        ///plot results

		Calibration cal = imp.getCalibration();

		// cal.setUnit("% radius");
		if (cal.getUnit() == "pixel")
			useCalibration = false;
		// cal.setUnit("% radius");
		// else {
		// normal

		min = Accumulator[1][0];
		for (int i = 0; i < nBins; i++) {
			Accumulator[1][i] = Accumulator[1][i] / Accumulator[0][i];
			Accumulator[0][i] = (float) (nBins * ((double) (i + 1) / nBins));
			if (Accumulator[1][i] < min) {
				min = Accumulator[1][i];
			}

		}
		// / bestimung std float[][] Accpixel = new
		// float[fp.npoints][nBins];
		float[][] accusum = new float[2][nBins];//

		// //varianz
		for (int i = 0; i < nBins; i++) {
			for (int j = 0; j < fp.npoints; j++) {
				accusum[1][i] += ((Accumulator[1][i] - Accpixel[j][i]) * (Accumulator[1][i] - Accpixel[j][i]))
						/ fp.npoints;
			}

		}
		// standard geviation
		for (int i = 0; i < nBins; i++) {
			if (fp.npoints > 1) {
				accusum[1][i] = (float) Math
						.sqrt((fp.npoints / (fp.npoints - 1.0)) * accusum[1][i]);
			} else {
				accusum[1][i] = (float) Math.sqrt(accusum[1][i]);
			}

		}

		// subbackground
		float[][] Accusubbackground = new float[1][nBins];

		for (int j = 0; j < nBins; j++) {
			Accusubbackground[0][j] = Accumulator[1][j] - min;// -Accumulator[1][1];
			Accumulator[0][j] = Accumulator[0][j];

		}
		// if (useCalibration) // in fale real radius
		if (realrad) {
			for (int i = 0; i < nBins; i++) {
				// Accumulator[1][i] = Accumulator[1][i] / Accumulator[0][i];
				Accumulator[0][i] = (float) (cal.pixelWidth * nBins * ((double) (i + 1) / nBins));

			}
		}

		Plot plot = null;
		// IJ.log("realradius " + realrad + "calibr ");
		String utitelzu = pluntertitel;
		if (realrad) {
			utitelzu = pluntertitel + " " + cal.getUnits() + " ";
		}
		// ""+ pluntertitel + " [" + cal.getUnits() + "]"

		if (!(subbackground)) {
			plot = new Plot("Clock Scan Profile Plot " + btitle, ""

			+ utitelzu + " ",

			"intensity, shades of grey", Accumulator[0], Accumulator[1]);
			if (paction) {
				plot.addErrorBars("SE", accusum[1]);
			}
		}
		// change nur subbackground
		if ((subbackground)) {
			plot = new Plot("Clock Scan Profile Plot " + btitle, ""

			+ utitelzu + "", "relative intensity, shades of grey",
					Accumulator[0], Accusubbackground[0]);
			if (paction) {
				// plot.addErrorBars(accusum[1]);
				plot.addErrorBars("SE", accusum[1]);
			}
		}

		plot.show();
	
	
		///// 
		WindowManager.setTempCurrentImage(imp);

	}// end of Radial Clock Scan


	
	// Roi Center
	private void setXYcenter() {
		fp = imp.getRoi().getFloatPolygon();
		X0 = fp.getBounds().getCenterX();
		Y0 = fp.getBounds().getCenterY();
	}

	// Paint Roi Center with radius 7
	private void doPaintCenter(ImageProcessor ip, double x, double y) {
		double mR = 7;
		IJ.makeOval((int) (x - mR), (int) (y - mR), (int) (2 * mR),
				(int) (2 * mR));
		Roi roi = imp.getRoi();
		String rname = "";
		roi.setName(center + " Y: " + (int) y + " X: " + (int) x);
		rname = center + " Y: " + (int) y + " X: " + (int) x;
		Selection sc = new Selection();
		sc.run("Fit Spline");
		RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();
		if (manager.getName() != rname)
			manager.addRoi(roi);

	}
	
	private void doSetInterpolateRoi() {
		imp = WindowManager.getCurrentImage();
		RoiManager rm = RoiManager.getInstance2();
		if (rm == null) {
			rm = new RoiManager();
		}
		Roi roi = imp.getRoi();
		RoiManager manager = RoiManager.getInstance();
		manager.reset();

		// roi interpolation
		if (roi != null) {
			if (roi.getName() != "interpolate") {
				interpolate();
				roi = imp.getRoi();

				roi.setName("ROI interpolated");
				manager.addRoi(roi);
			}
		}
	}

	private void doDialog() {
		canceled = false;
		GenericDialog gd = new GenericDialog("Clock Scan...", IJ.getInstance());
		fp = imp.getRoi().getFloatPolygon();
		// IJ.log("Anzahl point perimeter :"+ fp.npoints);
		gd.addNumericField("X center (pixels):", X0, 0);
		gd.addNumericField("Y center (pixels):", Y0, 0);
		gd.addMessage("scan limits (fraction of radius)");
		gd.addSlider(" ", 1, 2, 1.2);
		gd.addCheckbox("real radius ", realrad);
		gd.addCheckbox("subtract background", subbackground);
		gd.addCheckbox("polar transform ", polar);
		gd.addCheckbox("Plot with standard deviation ", paction);
		gd.addMessage("ROI selection length (pixels): " + fp.npoints);
		gd.addMessage("");
		String textcit = "If you use this plugin for your research,\n please cite the original article: \n Dobretsov & Romanovsky, \n Clock-scan protocol  for image analysis.AJP,\n Cell Physiology 291: C869-C879, 2006";
		gd.addMessage("");
		gd.addMessage(textcit);
		gd.showDialog();
		if (gd.wasCanceled()) {
			canceled = true;
			return;
		}
		X0 = gd.getNextNumber();
		Y0 = gd.getNextNumber();
		limits = gd.getNextNumber();
		realrad = gd.getNextBoolean();
		subbackground = gd.getNextBoolean();
		polar = gd.getNextBoolean();
		paction = gd.getNextBoolean();

		if (gd.invalidNumber()) {
			IJ.showMessage("Error", "Invalid input Number");
			canceled = true;
			return;
		}
	}


	
	
	// Roi interpolate
	public void interpolate() {
		Roi roi = imp.getRoi();
		if (roi == null) {
			return;
		}
		if (roi.getType() == Roi.POINT)
			return;
		double interval = 1;
		boolean smooth = true;
		Undo.setup(Undo.ROI, imp);
		boolean adjust = true;
		int sign = adjust ? -1 : 1;
		FloatPolygon poly = roi.getInterpolatedPolygon(sign * interval, smooth);
		int t = roi.getType();
		int type = roi.isLine() ? Roi.FREELINE : Roi.FREEROI;
		if (t == Roi.POLYGON && interval > 1.0)
			type = Roi.POLYGON;
		if ((t == Roi.RECTANGLE || t == Roi.OVAL || t == Roi.FREEROI)
				&& interval >= 8.0)
			type = Roi.POLYGON;
		if ((t == Roi.LINE || t == Roi.FREELINE) && interval >= 8.0)
			type = Roi.POLYLINE;
		if (t == Roi.POLYLINE && interval >= 8.0)
			type = Roi.POLYLINE;
		ImageCanvas ic = imp.getCanvas();
		if (poly.npoints <= 150 && ic != null && ic.getMagnification() >= 12.0)
			type = roi.isLine() ? Roi.POLYLINE : Roi.POLYGON;
		Roi p = new PolygonRoi(poly, type);
		if (roi.getStroke() != null)
			p.setStrokeWidth(roi.getStrokeWidth());
		p.setStrokeColor(roi.getStrokeColor());
		p.setName(roi.getName());
		transferProperties(roi, p);
		imp.setRoi(p);
	}

	// Roi transfer Properties
	private void transferProperties(Roi roi1, Roi roi2) {
		if (roi1 == null || roi2 == null)
			return;
		roi2.setStrokeColor(roi1.getStrokeColor());
		if (roi1.getStroke() != null)
			roi2.setStroke(roi1.getStroke());
		roi2.setDrawOffset(roi1.getDrawOffset());
	}

	// Roi scale to Center
	void scaletocenter(Roi roi, double xscale, double yscale, double x, double y) {

		FloatPolygon poly = roi.getFloatPolygon();
		int type = roi.getType();
		double xbase = x;
		double ybase = y;
		for (int i = 0; i < poly.npoints; i++) {
			poly.xpoints[i] = (float) ((poly.xpoints[i] - xbase) * xscale + xbase);
			poly.ypoints[i] = (float) ((poly.ypoints[i] - ybase) * yscale + ybase);
		}
		Roi roi2 = null;
		type = Roi.FREEROI;
		roi2 = new PolygonRoi(poly.xpoints, poly.ypoints, poly.npoints, type);
		roi2.setStrokeColor(roi.getStrokeColor());
		if (roi.getStroke() != null)
			roi2.setStroke(roi.getStroke());
		imp.setRoi(roi2);
	}

	public void cartesianto360() {

		heightTransform = heightInitial * 2;
		widthTransform = heightInitial * 2;
		getCartesianCenter();
		if (isColor)
			ipTransform = new ColorProcessor(widthTransform, heightTransform);		
		else
			ipTransform = new ShortProcessor(widthTransform, heightTransform);
		IJ.showStatus("Calculating...");
		for (int yy = 0; yy <= widthTransform; yy++) {
			for (int xx = 0; xx <= heightTransform; xx++) {
				double y = xx - centerY;
				double x = yy - centerX;
				double r = getRadius(x, y);
				double angle = getAngle(x, y);
				y = r;
				x = angle * (widthInitial / 360.5);
				if (isColor) {
					interpolateColorPixel(x, y);
					ipTransform.putPixel(xx, yy, rgbArray);
				} else {
					double newValue = ipInitial.getPixelInterpolated(x, y);
					ipTransform.putPixelValue(xx, yy, newValue);
				}
			}
			IJ.showProgress(yy, heightTransform);
		}
		IJ.showProgress(1.0);
	}

	void interpolateColorPixel(double x, double y) {

		int xL, yL;

		xL = (int) Math.floor(x);
		yL = (int) Math.floor(y);
		xLyL = ipInitial.getPixel(xL, yL, xLyL);
		xLyH = ipInitial.getPixel(xL, yL + 1, xLyH);
		xHyL = ipInitial.getPixel(xL + 1, yL, xHyL);
		xHyH = ipInitial.getPixel(xL + 1, yL + 1, xHyH);
		for (int rr = 0; rr < 3; rr++) {
			double newValue = (xL + 1 - x) * (yL + 1 - y) * xLyL[rr];
			newValue += (x - xL) * (yL + 1 - y) * xHyL[rr];
			newValue += (xL + 1 - x) * (y - yL) * xLyH[rr];
			newValue += (x - xL) * (y - yL) * xHyH[rr];
			rgbArray[rr] = (int) newValue;
		}
	}

	public void getCartesianCenter() {
		centerX = widthTransform / 2;
		centerY = heightTransform / 2;

	}

	double getRadius(double x, double y) {
		return Math.sqrt(x * x + y * y);
	}

	double getAngle(double x, double y) {

		double angle = Math.toDegrees(Math.atan2(y, x));
		if (angle <= 0) {
			angle += 360;
		}
		return clockWise ? 360 - angle : angle;
	}

	//
	void showAbout() {
		IJ.showMessage("Clock Skan ",
				"you need a selection for use the plugin.\n"
						+ "If you use this plugin for your research,\n"
						+ "please cite the original article:\n"
						+ "Dobretsov & Romanovsky, \n "
						+ "Clock-scan protocol  for image analysis.AJP,\n "
						+ "Cell Physiology 291: C869-C879, 2006");
	}

}


