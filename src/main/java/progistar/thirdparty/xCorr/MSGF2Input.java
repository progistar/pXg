package progistar.thirdparty.xCorr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class MSGF2Input {

	public static final int MSGF_PEPTIDE_IDX = 9;
	public static final int MSGF_TITLE_IDX = 3;
	
	public static void main(String[] args) throws IOException {
		String outputName = "msgfXCorrInput.mgf";
		File[] msgfIDList = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/unmodified/fdr_q5").listFiles();
		File[] mgfs = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Cal_MGF_AddScanNum").listFiles();
		
		Hashtable<String, String> msgfPSMs = new Hashtable<String, String>();
		// aggregate PSMs
		for(File file : msgfIDList) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			
			if(file.getName().endsWith(".fdr.BA")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				BR.readLine(); // skip header
				
				while((line = BR.readLine()) != null) {
					String[] fields = line.split("\t");
					//_uncalibrated
					String key = fields[MSGF_TITLE_IDX];
					String peptide = fields[MSGF_PEPTIDE_IDX];
					if(msgfPSMs.get(key) == null) {
						peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
						msgfPSMs.put(key, peptide);
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
		msgfPSMs.forEach((title, peptide)->{
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
