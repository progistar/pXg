package progistar.thirdparty.sSim;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Spectra {

	public static final int FILE_TYPE_MZXML = 0;
	public static final int FILE_TYPE_MGF = 1;
	
	private ArrayList<Spectrum> spectra = null;
	private TreeMap<Double, Spectrum> sortedRTSpectra = null;
	private Hashtable<Integer, Spectrum> scanIndexer = null;
	private double minRT = Double.MAX_VALUE;
	private double maxRT = 0;
	
	private double avgRTInterval = 0;
	private double varRTInterval = 0;
	
	private String fileName = null;
	private int fileType = -1;
	
	private Hashtable<String, String> targetedScans = null;
	
	public Spectra (String fileName, int fileType, Hashtable<String, String> targetedScans) {
		this.fileName = fileName;
		this.fileType = fileType;
		this.targetedScans = targetedScans;
		if(this.fileType == FILE_TYPE_MGF) readMGF (this.fileName, null);
	}
	
	
	public Spectra (String fileName, int fileType, ArrayList<Integer> selectiveScans) {
		this.fileName = fileName;
		this.fileType = fileType;
		if(this.fileType == FILE_TYPE_MGF) readMGF (this.fileName, selectiveScans);
	}
	
	//TODO: READ MGF
	private void readMGF (String fileName, ArrayList<Integer> selectiveScans) {
		spectra = new ArrayList<Spectrum>();
		sortedRTSpectra = new TreeMap<Double, Spectrum>();
		scanIndexer = new Hashtable<Integer, Spectrum>();
		
		try {
			BufferedReader BR = new BufferedReader(new FileReader(fileName));
			String line = null;
			int index = -1;
			double pepMass = 0;
			double precursorInt = 0;
			int charge = 0;
			double rt = 0;
			String title = null;
			int scanNum = -1;
			
			ArrayList<double[]> peaks = null;
			
			Pattern peakPattern = Pattern.compile("^[0-9]");
			Pattern scanPattern = Pattern.compile("[\\.]+[0-9]+[\\.]+");
			
			while((line = BR.readLine()) != null) {
				if(line.startsWith("BEGIN")) {
					peaks = new ArrayList<double[]>();
					index ++;
				}else if(line.startsWith("TITLE")) {
					title = line.split("\\=")[1];
					Matcher matcher = scanPattern.matcher(title);
					if(matcher.find()) {
						scanNum = Integer.parseInt(matcher.group().replace(".", ""));
					}
				}else if(line.startsWith("RTIN")) {
					rt = Double.parseDouble(line.split("\\=")[1]);
				}else if(line.startsWith("PEPMASS")) {
					pepMass = Double.parseDouble(line.split("\\s")[0].split("\\=")[1]);
					try {
						precursorInt = Double.parseDouble(line.split("\\s")[1]);
					}catch(Exception e) {
						System.out.println(line+"");
						e.printStackTrace();
					}
				}else if(line.startsWith("CHARGE")){
					charge = Integer.parseInt(line.split("\\=")[1].replace("+", ""));
				}else if(peakPattern.matcher(line).find()) {
					double[] peak = new double[2];
					String[] peakStr = line.split("\\s");
					peak[0] = Double.parseDouble(peakStr[0]);
					peak[1] = Double.parseDouble(peakStr[1]);
					peaks.add(peak);
				}else if(line.startsWith("END")) {
					if(targetedScans.get(title) == null) {
						continue;
					}
					Spectrum spectrum = new Spectrum(scanNum, charge, 2, pepMass, peaks, rt, index);
					spectrum.setPrecursorInt(precursorInt);
					spectrum.setTitle(title);
					spectra.add(spectrum);
					if(sortedRTSpectra.get(rt) != null) System.err.println("Duplicated RT was observed!");
					sortedRTSpectra.put(rt, spectrum);
					
					maxRT = maxRT < rt ? rt : maxRT;
					minRT = minRT > rt ? rt : minRT;
					
					// if scanNum is -1, then it is missing value.
					if(scanNum != -1) scanIndexer.put(scanNum, spectrum);
				}
			}
			
			BR.close();
		}catch(IOException ioe) {
			
		}
	}
	
	public ArrayList<Spectrum> getSpectraByRT (double from, double to) {
		ArrayList<Spectrum> rtSpectra = new ArrayList<Spectrum>();
		SortedMap<Double, Spectrum> list = sortedRTSpectra.subMap(from, to);
		if(list != null) {
			Iterator<Double> iter = (Iterator<Double>)list.keySet().iterator();
			while(iter.hasNext()) {
				rtSpectra.add(this.sortedRTSpectra.get(iter.next()));
			}
		}
		return rtSpectra;
	}
	
	public ArrayList<Spectrum> getMS1SpectraByRT (double from, double to) {
		ArrayList<Spectrum> rtSpectra = new ArrayList<Spectrum>();
		SortedMap<Double, Spectrum> list = sortedRTSpectra.subMap(from, to);
		if(list != null) {
			Iterator<Double> iter = (Iterator<Double>)list.keySet().iterator();
			while(iter.hasNext()) {
				Spectrum scan = this.sortedRTSpectra.get(iter.next());
				if(scan.getMsLevel() == 1) rtSpectra.add(scan);
			}
		}
		return rtSpectra;
	}
	
	public void normalizationByBasePeak () {
		
	}
	
	public int sizeOfSpectra () {
		return this.spectra.size();
	}
	
	/**
	 * 
	 * The return value is shallow copy of spectrum.
	 * 
	 * @param scanNum
	 * @return
	 */
	public Spectrum getSpectrumByScanNum (int scanNum) {
		Spectrum spectrum = this.scanIndexer.get(scanNum);
		if(spectrum == null) System.err.println("getSpectrumByScanNum: wrong scanNum! " + scanNum);
		return spectrum;
	}
	
	public int getFileType () {
		return this.fileType;
	}
	
	/**
	 * index is 0-based.
	 * 
	 * 
	 * @param index
	 * @return
	 */
	public Spectrum getSpectrumByIndex (int index) {
		return this.spectra.get(index);
	}
	
	public String getFileName () {
		return this.fileName;
	}
	
	public double getMaxRT () {
		return this.maxRT;
	}
	
	public double getMinRT () {
		return this.minRT;
	}
	/**
	 * 
	 * 
	 * 
	 * @param fileName
	 * @param msLevel
	 * @throws IOException
	 */
	public void writeToMGF (String fileName, int[] msLevels) throws IOException {
		Spectrum spectrum = null;
		int size = this.sizeOfSpectra();
		File outputFile = new File(fileName);
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		String titleHeader = outputFile.getName().substring(0, outputFile.getName().lastIndexOf("."));
		for(int i=0; i<size; i++) {
			spectrum = this.getSpectrumByIndex(i);
			if(spectrum != null) {
				int msLevel = spectrum.getMsLevel();
				
				boolean isTargeted = false;
				for(int targetMSLevel : msLevels) if(msLevel == targetMSLevel) isTargeted = true;
				
				if(isTargeted) {
					BW.append("BEGIN IONS");
					BW.newLine();
					BW.append("TITLE=").append(titleHeader).append(".").append(spectrum.getScanNum()+".").append(spectrum.getScanNum()+".").append(spectrum.getCharge()+"");
					BW.newLine();
					BW.append("RTINSECONDS=").append(spectrum.getRT()+"");
					BW.newLine();
					BW.append("PEPMASS=").append(spectrum.getPrecursorMz()+"\t").append(spectrum.getPrecursorInt()+"");
					BW.newLine();
					BW.append("CHARGE=").append(spectrum.getCharge()+"+");
					BW.newLine();
					int peakSize = spectrum.sizeOfPeaks();
					for(int j=0; j<peakSize; j++) {
						BW.append(spectrum.getPeak(j)[0]+"\t"+spectrum.getPeak(j)[1]);
						BW.newLine();
					}
					
					BW.append("END IONS");
					BW.newLine();
				}
			}
		}
		
		BW.close();
	}
	
	
}