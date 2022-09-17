package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class SpectrumGen {
	
	public static int ProteomToolsFileIdx = 0;
	public static int ProteomToolsScanIdx = 1;
	public static int ProteomToolsPeptideIdx = 3;
	public static int ProteomToolsChargeIdx = 12;
	
	public static int pXgFileIdx = 1;
	public static int pXgScanIdx = 4;
	public static int pXgPeptideIdx = 20;
	public static int pXgChargeIdx = 10;
	
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
				String peptide = fields[ProteomToolsPeptideIdx];
				String title = fields[ProteomToolsFileIdx]+"."+fields[ProteomToolsScanIdx]+"."+fields[ProteomToolsScanIdx]+"."+fields[ProteomToolsChargeIdx];
				String fileName = fields[ProteomToolsFileIdx]+".mgf";
				if(targetPeptides.get(peptide) != null) {
					msmsRecords.add(line);
					idedScanTitleMapper.put(title, peptide);
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
	
	public static void loadpXg(String resFileName, String specFolder, Hashtable<String, String> targetPeptides) throws IOException {
		
		File resFile = new File(resFileName);
		ArrayList<File> msmsFiles = new ArrayList<File>();
		msmsFiles.add(resFile);
		
		
		Hashtable<String, String> idedScanTitleMapper = new Hashtable<String, String>();
		Hashtable<String, String> idedFileMapper = new Hashtable<String, String>();
		for(File file : msmsFiles) {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[pXgPeptideIdx];
				String scanNum = fields[pXgScanIdx].split("\\:")[1];
				String title = fields[pXgFileIdx].split("\\.")[0]+"."+scanNum+"."+scanNum+"."+fields[pXgChargeIdx];
				String fileName = fields[pXgFileIdx].split("\\.")[0]+".mgf";
				if(targetPeptides.get(peptide) != null) {
					msmsRecords.add(line);
					idedScanTitleMapper.put(title, peptide);
					idedFileMapper.put(fileName, "");
				}
			}
			
			BR.close();
		}
		
		File[] fileList = new File(specFolder).listFiles();
		
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
