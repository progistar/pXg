package progistar.thirdparty.xCorr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class pXg2Input {

	public static final int PXG_SPEC_FILE_IDX = 1;
	public static final int PXG_SCAN_IDX = 4;
	public static final int PXG_CHARGE_IDX = 10;
	public static final int PXG_PEPTIDE_IDX = 20;
	public static final int PXG_IS_CANONICAL = 36;
	
	public static void main(String[] args) throws IOException {
		String outputName = "pXgXCorrInput.mgf";
		File[] pXgIDList = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg").listFiles();
		File[] mgfs = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Cal_MGF_AddScanNum").listFiles();
		
		Hashtable<String, String> pxgPSMs = new Hashtable<String, String>();
		// aggregate PSMs
		for(File file : pXgIDList) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			
			if(file.getName().endsWith(".pXg.BA")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				BR.readLine(); // skip header
				
				while((line = BR.readLine()) != null) {
					String[] fields = line.split("\t");
					//_uncalibrated
					String fileName = fields[PXG_SPEC_FILE_IDX].split("\\.")[0];
					fileName = fileName.replace("_uncalibrated", "");
					fileName = fileName.replace("_calibrated", "");
					
					String scanNum = fields[PXG_SCAN_IDX].split("\\:")[1];
					String charge = fields[PXG_CHARGE_IDX];
					String key = fileName+"."+scanNum+"."+scanNum+"."+charge;
					String peptide = fields[PXG_PEPTIDE_IDX];
					if(pxgPSMs.get(key) == null) {
						peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
						pxgPSMs.put(key, peptide);
					}
				}
				
				BR.close();
			}
		}
		
		Hashtable<String, ArrayList<String>> spectra = new Hashtable<String, ArrayList<String>>();
		for(File file : mgfs) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			
			if(file.getName().endsWith(".mgf")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				
				ArrayList<String> spectrum = new ArrayList<String>();
				String title = null;
				while((line = BR.readLine()) != null) {
					spectrum.add(line);
					if(line.startsWith("END")) {
						if(title != null) {
							spectra.put(title, spectrum);
							spectrum = new ArrayList<String>();
							title = null;
						}
					} else if(line.startsWith("TITLE=")) {
						title = line.split("\\=")[1];
					}
				}
				
				BR.close();
			}
		}
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
		pxgPSMs.forEach((title, peptide)->{
			ArrayList<String> spectrum = spectra.get(title);
			if(spectrum == null) {
				System.out.println(title+" is mssing");
			} else {
				for(String str : spectrum) {
					try {
						BW.append(str);
						BW.newLine();
						if(str.startsWith("PEPMASS=")) {
							BW.append("SEQ="+peptide);
							BW.newLine();
						}
					}catch(IOException ioe) {
						
					}
					
				}
			}
		});
		BW.close();
	}
}
