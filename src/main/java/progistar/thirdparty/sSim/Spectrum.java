package progistar.thirdparty.sSim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class Spectrum {

	public static final double INVALID_PEAK = Double.MAX_VALUE;

	private String title = null;
	private int scanNum = -1;
	private int charge = 0;
	private int msLevel = 0;
	private int index = 0; // this is one-based.
	private double retentionTime = 0;
	private double precursorMz = 0;
	private double precursorInt = 0;
	private double[] basePeak = null;
	public ArrayList<double[]> peaks = null;

	// ad-hoc
	private String peptide = null;
	public void setPeptide(String peptide) {
		this.peptide = peptide;
	}

	public String getPeptide () {
		return this.peptide;
	}

	public Spectrum deepCopy() {
		ArrayList<double[]> newPeaks = new ArrayList<>();
		Spectrum s = new Spectrum(scanNum, charge, msLevel, precursorMz, newPeaks, retentionTime, index);
		s.peptide = this.peptide;
		for(double[] peak : this.peaks) {
			double[] newPeak = {peak[0], peak[1]};
			newPeaks.add(newPeak);
		}

		return s;
	}

	/**
	 * If you want to skip putting peak list, then just give null.<br>
	 *
	 * @param scanNum
	 * @param charge
	 * @param msLevel
	 * @param precursorMz
	 * @param peakList
	 * @param retentionTime
	 */
	public Spectrum (int scanNum, int charge, int msLevel, double precursorMz, ArrayList<double[]> peakList, double retentionTime, int index) {
		this.scanNum = scanNum;
		this.charge = charge;
		this.msLevel = msLevel;
		this.precursorMz = precursorMz;
		this.retentionTime = retentionTime;
		this.index = index;
		this.peaks = peakList;
		if(this.peaks == null) {
			this.peaks = new ArrayList<>();
		}
	}

	/**
	 * Constructor of mzXML.
	 *
	 *
	 * @param scanNum
	 * @param charge
	 * @param msLevel
	 * @param precursorMz
	 * @param peakList
	 * @param retentionTime
	 * @param index
	 */
	public Spectrum (int scanNum, int charge, int msLevel, double precursorMz, double[][] peakList, double retentionTime, int index) {
		this.scanNum = scanNum;
		this.charge = charge;
		this.msLevel = msLevel;
		this.precursorMz = precursorMz;
		this.retentionTime = retentionTime;
		this.index = index;
		this.loadPeaksFromMzXML(peakList);
	}

	public double getTotalIonChromatogram () {
		double TIC = .0;

		for(double[] peak : this.peaks) {
			TIC += peak[1];
		}

		return TIC;
	}

	public String getTitle () {
		return this.title;
	}

	public int getMsLevel () {
		return this.msLevel;
	}

	public double getMz (int index) {
		return this.peaks.get(index)[0];
	}

	public double getIntensity (int index) {
		return this.peaks.get(index)[1];
	}


	public int sizeOfPeaks () {
		return this.peaks.size();
	}

	public void addPeak (double[] peak) {
		if(this.peaks == null) {
			this.peaks = new ArrayList<>();
		}
		this.peaks.add(peak);
	}

	public double[] getPeak (int index) {
		return this.peaks.get(index);
	}

	public void removePeak (int index) {
		this.peaks.remove(index);
	}

	public ArrayList<double[]> getPeaks () {
		return this.peaks;
	}

	public double getPrecursorMz () {
		return this.precursorMz;
	}

	public double getPrecursorInt () {
		return this.precursorInt;
	}

	public int getCharge () {
		return this.charge;
	}

	public double getRT () {
		return this.retentionTime;
	}

	public void setTitle (String title) {
		this.title = title;
	}

	public void setRT (double RT) {
		this.retentionTime = RT;
	}

	public void setPrecursorInt (double precursorInt) {
		this.precursorInt = precursorInt;
	}

	/**
	 * loading peak information from mzXML.<br>
	 *
	 *
	 *
	 * @param peaks
	 */
	public void loadPeaksFromMzXML (double[][] peaks) {
		if(this.peaks == null) {
			this.peaks = new ArrayList<>();
		}
		int length = peaks[0].length;
		for(int i=0; i<length; i++) {
			double[] peak = new double[2];
			peak[0] = peaks[0][i];
			peak[1] = peaks[1][i];
			this.peaks.add(peak);
		}
	}
	/**
	 * mzOrInt: mz order = 0, int order = 1.<br>
	 *
	 * @param mzOrInt
	 * @param isIncreasingOrder
	 */
	public void sortPeaks (final int mzOrInt, final boolean isIncreasingOrder) {
		Collections.sort(this.peaks, new Comparator<double[]>() {
			@Override
			public int compare(double[] peak1, double[] peak2) {
				if(peak1[mzOrInt] > peak2[mzOrInt]) {
					return isIncreasingOrder ? 1 : -1;
				} else if(peak1[mzOrInt] < peak2[mzOrInt]) {
					return isIncreasingOrder ? -1 : 1;
				}
				return 0;
			}
		});
	}

	public int getScanNum () {
		return this.scanNum;
	}

	public int getIndex () {
		return this.index;
	}

	public double[] getBasePeak () {
		int sizeOfPeaks = this.sizeOfPeaks();
		this.basePeak = this.peaks.get(0);
		for(int i=1; i<sizeOfPeaks; i++) {
			if(this.peaks.get(i)[1] > this.basePeak[1]) {
				this.basePeak = this.peaks.get(i);
			}
		}
		return this.basePeak;
	}

	public void normalizationByValue(double value) {
		for(double[] peak : this.peaks) {
			peak[1] /= value;
		}
	}


	public void rankedWeightedSum (double tolerance) {
		this.sortPeaks(1, false);
		int sizeOfPeaks = this.sizeOfPeaks();
		double[] targetPeak = this.getPeak(0);
		for(int i=1; i<sizeOfPeaks; i++) {
			double[] pivotPeak = this.getPeak(i);
			if(pivotPeak[0] - targetPeak[0] < tolerance) {
				// weighted sum
				targetPeak[0] = (targetPeak[0] * targetPeak[1] + pivotPeak[0] * pivotPeak[1]) / (targetPeak[1] + pivotPeak[1]);
				targetPeak[1] += pivotPeak[1];
				pivotPeak[0] = INVALID_PEAK;
				pivotPeak[1] = INVALID_PEAK;
			} else {
				targetPeak = pivotPeak;
			}
		}

		this.sortPeaks(0, true);

		int index = sizeOfPeaks-1;
		while(index >= 0 && this.getMz(index) == INVALID_PEAK) {
			this.peaks.remove(index--);
		}
	}

	/**
	 * get subpeaks using binary search.
	 *
	 * @param fromMz
	 * @param toMz
	 * @return
	 */
	public ArrayList<double[]> getSubPeaks (double fromMz, double toMz) {
		ArrayList<double[]> subPeaks = new ArrayList<>();
		this.sortPeaks(0, true);

		int sizeOfPeaks = this.sizeOfPeaks();
		int leftBound = 0;
		int rightBound = sizeOfPeaks-1;
		int startIndex = 0;
		while(leftBound <= rightBound) {
			int mid = (leftBound+rightBound) / 2;
			startIndex = mid;
			if(this.getMz(mid) > fromMz) {
				rightBound = mid-1;
			} else if(this.getMz(mid) < fromMz) {
				leftBound = mid+1;
			} else {
				break;
			}
		}

		startIndex -= 1;
		if(startIndex < 0) {
			startIndex = 0;
		}

		while(startIndex < sizeOfPeaks) {
			double mz = this.getMz(startIndex);
			if(mz >= fromMz && mz <= toMz) {
				subPeaks.add(this.getPeak(startIndex));
			} else if(mz > toMz){
				break;
			}
			startIndex++;
		}

		return subPeaks;
	}
}