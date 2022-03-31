package progistar.thirdparty.peaks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class PTMFrequency {
	
	public static int peptideWithPTMIdx = 3;
	public static int inferredPeptideIdx = 20;
	public static int isCannonicalIdx = 36;
	public static int baClass = 43;
	public static int bestHLAType = 44;

	// take pXg with BA prediction out.
	
	public static void main(String[] args) throws IOException {
		ArrayList<String[]> records = readPXG("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/5.withCalibrationAddScanNumWithoutDeami/pXg/S2.BA.fdr");
		
		showSiteFrequency(records, "C(+119.00)");
	}
	
	// return:
	// ArrayList of record information
	// Note that the first record is a header.
	public static ArrayList<String[]> readPXG (String fileName) throws IOException {
		File file = new File(fileName);
		ArrayList<String[]> records = new ArrayList<String[]>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			records.add(line.split("\t"));
		}
		
		BR.close();
		return records;
	}
	/**
	 * ptm = C(+119.00)
	 * 
	 * @param records
	 * @param ptm
	 * @throws IOException
	 */
	public static void showSiteFrequency (ArrayList<String[]> records, String ptm) throws IOException {
		
		Hashtable<String, String> list = new Hashtable<String, String>();
		
		for(int i=1; i<records.size(); i++) {
			String[] record = records.get(i);
			
			String peptide = record[peptideWithPTMIdx];
			String inferredPeptide = record[inferredPeptideIdx];
			String isCanonical = record[isCannonicalIdx];
			
			if(isCanonical.equalsIgnoreCase("TRUE")) {
				continue;
			}
			// pass not interesting
			if(!peptide.contains(ptm)) {
				continue;
			}
			
			StringBuilder inferredPeptideWithPTM = new StringBuilder();
			int pIdx = 0;
			int ipIdx = 0;
			
			for(pIdx=0; pIdx < peptide.length(); pIdx++) {
				// append modification information
				if(peptide.charAt(pIdx) == '(') {
					while(peptide.charAt(pIdx) != ')') {
						inferredPeptideWithPTM.append(peptide.charAt(pIdx++));
					}
					inferredPeptideWithPTM.append(peptide.charAt(pIdx));
				} else {
					inferredPeptideWithPTM.append(inferredPeptide.charAt(ipIdx++));
				}
			}
			
			list.put(inferredPeptideWithPTM.toString(), inferredPeptide);
		}
		
		// Length of interest
		for(int i=8; i<=15; i++) {
			int[] sites = new int[i+1];
			
			list.forEach((peptide, strip)->{
				if(strip.length() == sites.length-1) {

					int site = 0;
					for(int j=0; j<peptide.length(); j++) {
						if(peptide.charAt(j) == '('){
							String mod = peptide.charAt(j-1)+"";
							while(peptide.charAt(j) != ')') {
								mod += peptide.charAt(j);
								j++;
							}
							mod += peptide.charAt(j);
							if(ptm.equalsIgnoreCase(mod)) {
								sites[site]++;
							}
							
						} else if(peptide.charAt(j) == strip.charAt(site)) {
							site++;
						}
					}
				}
				
				
			});
			
			System.out.println(i);
			for(int site=1; site<=i; site++) {
				System.out.print(sites[site]+"\t");
			}
			System.out.println();
		}
	}
}
