package progistar.thirdparty.precursorswap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

class Spectrum implements Comparable<Spectrum> {
	public double pepMass = .0;
	public int charge = 0;
	public ArrayList<String> records = new ArrayList<String>();
	boolean isSwapped = false;
	@Override
	public int compareTo(Spectrum s) {
		if(this.pepMass < s.pepMass) {
			return -1;
		} else if(this.pepMass > s.pepMass) {
			return 1;
		}
		return 0;
	}
}

public class PSDG {

	public static double d = 14;
	public static String suffixName = "3";
	
	public static void main(String[] args) throws IOException {
		String filePath = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Cal_MGF_AddScanNum";
		String outputDir = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Cal_MGF_AddScanNum_PSDG/";
		File[] fileList = new File(filePath).listFiles();
		
		for(File file : fileList) {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String outputFileName = outputDir + file.getName().replace("mgf", d+"_"+suffixName+".mgf");
			ArrayList[] spectraPerCharge = new ArrayList[10];
			
			for(int i=0; i<10; i++) {
				spectraPerCharge[i] = new ArrayList<Spectrum>();
			}
			
			String line = null;
			
			Spectrum spectrum = null;
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			while((line = BR.readLine()) != null) {
				if(line.startsWith("BEGIN")) {
					spectrum = new Spectrum();
					spectra.add(spectrum);
				}else if(line.startsWith("CHARGE")) {
					int charge = Integer.parseInt(line.split("\\=")[1].charAt(0)+"");
					spectrum.charge = charge;
				}else if(line.startsWith("PEPMASS")) {
					double pepmass = Double.parseDouble(line.split("\\=")[1].split("\\s")[0]);
					spectrum.pepMass = pepmass;
				}
				spectrum.records.add(line);
			}
			
			BR.close();
			
			System.out.println("A total of spectra: "+spectra.size());
			for(int i=0; i<spectra.size(); i++) {
				spectrum = spectra.get(i);
				int charge = spectrum.charge;
				spectraPerCharge[charge].add(spectrum);
			}

			// precursor swap
			for(int i=0; i<spectraPerCharge.length; i++) {
				spectra = spectraPerCharge[i];
				long seed = System.currentTimeMillis();
				ArrayList<Spectrum> swappedSpectra = new ArrayList<Spectrum>();
				if(spectra.size() != 0) {
					System.out.println("Charge "+i+": "+spectra.size());
					// swap
					Random random = new Random(seed);
					
					for(int trial = 0; trial < 10; trial++) {
						for(int j=0; j< spectra.size(); j++) {
							Spectrum pivotSpectrum = spectra.get(j);
							if(pivotSpectrum.isSwapped) {
								continue;
							}
							for(int k=j+1; k<spectra.size(); k++) {
								Spectrum targetSpectrum = spectra.get(k);
								if(targetSpectrum.isSwapped) {
									continue;
								}
								if(targetSpectrum.pepMass - pivotSpectrum.pepMass > d) {
									int r = random.nextInt(2);
									if(r == 0) {
										pivotSpectrum.isSwapped = true;
										targetSpectrum.isSwapped = true;
										
										double pepMass = pivotSpectrum.pepMass;
										pivotSpectrum.pepMass = targetSpectrum.pepMass;
										targetSpectrum.pepMass = pepMass;
										
										swappedSpectra.add(pivotSpectrum);
										swappedSpectra.add(targetSpectrum);
										break;
									}
								}
							}
						}
					}
					
					int swapped = 0;
					int nonswapped = 0;
					for(int j=0; j<spectra.size(); j++) {
						Spectrum pivotSpectrum = spectra.get(j);
						if(!pivotSpectrum.isSwapped) {
							pivotSpectrum.pepMass += d;
							swappedSpectra.add(pivotSpectrum);
							nonswapped++;
						} else {
							swapped++;
						}
					}
					System.out.println(swapped+"\t"+nonswapped);
				}
				
				spectraPerCharge[i] = swappedSpectra;
			}
			
			// write file
			BufferedWriter BW = new BufferedWriter(new FileWriter(outputFileName));
			
			for(int i=0; i<spectraPerCharge.length; i++) {
				spectra = spectraPerCharge[i];
				for(Spectrum spec : spectra) {
					for(String record : spec.records) {
						if(record.startsWith("PEPMASS")) {
							BW.append("PEPMASS="+spec.pepMass+" 0.0");
						} else {
							BW.append(record);
						}
						BW.newLine();
					}
				}
			}
			
			BW.close();
			

		}
	}
}
