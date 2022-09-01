package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class ProteomTools {
	
	public static int fileIdx = 0;
	public static int scanIdx = 1;
	public static int peptideIdx = 3;
	public static int chargeIdx = 12;
	
	public static ArrayList<String> msmsRecords = new ArrayList<String>();
	public static ArrayList<Spectra> spectraList = new ArrayList<Spectra>();

	public static void loadProteomeTools(String resFolder, String specFolder, Hashtable<String, String> targetPeptides) throws IOException {
		
		File[] fileList = new File(resFolder).listFiles();
		ArrayList<File> msmsFiles = new ArrayList<File>();
		
		// read msms file
		for(File file : fileList) {
			if(file.isDirectory()) {
				File[] subFiles = file.listFiles();
				for(File msFile : subFiles) {
					// msms.txt file identifier
					if(msFile.getName().equalsIgnoreCase("msms.txt")) {
						msmsFiles.add(msFile);
					}
				}
			}
		}
		
		Hashtable<String, String> idedScanTitleMapper = new Hashtable<String, String>();
		Hashtable<String, String> idedFileMapper = new Hashtable<String, String>();
		for(File file : msmsFiles) {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[peptideIdx];
				String title = fields[fileIdx]+"."+fields[scanIdx]+"."+fields[scanIdx]+"."+fields[chargeIdx];
				String fileName = fields[fileIdx]+".mgf";
				if(targetPeptides.get(peptide) != null) {
					msmsRecords.add(line);
					idedScanTitleMapper.put(title, "");
					idedFileMapper.put(fileName, "");
				}
			}
			
			BR.close();
		}
		
		fileList = new File(specFolder).listFiles();
		
		// read mgf
		for(File file : fileList) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".mgf")) continue;
			if(idedFileMapper.get(file.getName()) == null) continue;
			
			System.out.println("matched MGF file: "+file.getName());
			Spectra spectra = new Spectra(file.getAbsolutePath(), Spectra.FILE_TYPE_MGF, idedScanTitleMapper);
			spectraList.add(spectra);
		}
		
	}
}
