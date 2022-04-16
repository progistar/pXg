package progistar.thirdparty.xCorr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class pXgAppender {
	public static final int PXG_SPEC_FILE_IDX = 1;
	public static final int PXG_SCAN_IDX = 4;
	public static final int PXG_CHARGE_IDX = 10;
	public static final int PXG_PEPTIDE_IDX = 20;
	public static final int PXG_IS_CANONICAL = 36;
	
	public static void main(String[] args) throws IOException {
		File xCorrTSV = new File("");
		File[] pXgIDList = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg").listFiles();
		
		
		BufferedReader BR_ = new BufferedReader(new FileReader(xCorrTSV));
		String line = null;
		
		BR_.readLine();// skip header
		Hashtable<String, String> titleToXcorr = new Hashtable<String, String>();
		while((line = BR_.readLine()) != null) {
			String[] fields = line.split("\t");
			String title = fields[0];
			String xCorr = fields[5];
			
			titleToXcorr.put(title, xCorr);
		}
		
		BR_.close();
		
		Hashtable<String, String> pxgPSMs = new Hashtable<String, String>();
		// aggregate PSMs
		for(File file : pXgIDList) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			
			if(file.getName().endsWith(".pXg.BA")) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".xCorr"));
				String header = BR.readLine(); // skip header
				header += "\txCorr";
				
				BW.append(header);
				BW.newLine();
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
						pxgPSMs.put(key, line);
						
						String xCorr = titleToXcorr.get(key);
						BW.append(line+"\t"+xCorr);
						BW.newLine();
					}
				}
				BW.close();
				BR.close();
			}
		}
		
	}
}
