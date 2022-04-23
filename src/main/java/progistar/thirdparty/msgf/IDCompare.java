package progistar.thirdparty.msgf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class IDCompare {

	public static final int MSGF_SPEC_FILE_IDX = 0;
	public static final int MSGF_SCAN_IDX = 2;
	public static final int MSGF_PEPTIDE_IDX = 9;
	
	public static final int PXG_SPEC_FILE_IDX = 1;
	public static final int PXG_SCAN_IDX = 4;
	public static final int PXG_PEPTIDE_IDX = 20;
	public static final int PXG_IS_CANONICAL = 36;
	
	public static Hashtable<String, String> readPXG (File file, boolean isCanonical) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		BR.readLine(); // skip header
		
		Hashtable<String, String> pxgPSMs = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[PXG_SPEC_FILE_IDX] +"_"+fields[PXG_SCAN_IDX].split("\\:")[1];
			String peptide = fields[PXG_PEPTIDE_IDX];
			String type = fields[PXG_IS_CANONICAL];
			
			if(isCanonical == type.equalsIgnoreCase("TRUE")) {
				if(pxgPSMs.get(key) == null) {
					peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
					pxgPSMs.put(key, peptide);
				}
			}
			
		}
		
		BR.close();
		
		return pxgPSMs;
	}
	
	public static Hashtable<String, String> readMSGF (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		Hashtable<String, String> msgfPSMs = new Hashtable<String, String>();
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[MSGF_SPEC_FILE_IDX] +"_"+fields[MSGF_SCAN_IDX];
			String peptide = fields[MSGF_PEPTIDE_IDX];
			peptide = peptide.replaceAll("[+-.0123456789*]", "");
			msgfPSMs.put(key, peptide);
		}
		
		BR.close();
		
		return msgfPSMs;
	}
	
	public static void main(String[] args) throws IOException {
		File[] msgfResSet = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/unmodified/fdr_q5").listFiles();
		File[] pxgResSet = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg").listFiles();
		
		Hashtable<String, String> msgfResHash = new Hashtable<String, String>();
		Hashtable<String, String> pXgResHash = new Hashtable<String, String>();
		
		for(File file : msgfResSet) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			if(file.getName().endsWith("BA.xCorr")) {
				msgfResHash.putAll(readMSGF(file));
			}
		}
		for(File file : pxgResSet) {
			if(file.getName().startsWith(".")) {
				continue;
			}
			if(file.getName().endsWith(".BA.xCorr")) {
				pXgResHash.putAll(readPXG(file, true));
			}
		}
		
		System.out.println("A total of MSGF IDs: "+msgfResHash.size());
		System.out.println("A total of pXg IDs : "+pXgResHash.size());
		
		int[] intersection = new int[1];
		pXgResHash.forEach((key, peptide) -> {
			String peptideMSGF = msgfResHash.get(key);
			if(peptideMSGF != null && peptideMSGF.equalsIgnoreCase(peptide)) {
				intersection[0]++;
			}
		});
		
		System.out.println("Intersected IDs: "+intersection[0]);
		
		
	}
}
