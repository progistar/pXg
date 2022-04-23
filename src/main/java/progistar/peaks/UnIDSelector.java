package progistar.peaks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;



public class UnIDSelector {

	public static final int PXG_SPEC_FILE_IDX = 1;
	public static final int PXG_SCAN_IDX = 4;
	public static final int PXG_PEPTIDE_IDX = 20;
	public static final int PXG_IS_CANONICAL = 36;
	
	public static Hashtable<String, String> readPXG (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		BR.readLine(); // skip header
		
		Hashtable<String, String> pxgPSMs = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[PXG_SPEC_FILE_IDX] +"_"+fields[PXG_SCAN_IDX].split("\\:")[1];
			String peptide = fields[PXG_PEPTIDE_IDX];
			String type = fields[PXG_IS_CANONICAL];
			if(pxgPSMs.get(key) == null) {
				peptide = peptide.replaceAll("[+-.\\(\\)0123456789*]", "");
				pxgPSMs.put(key, peptide);
			}
			
		}
		
		BR.close();
		
		return pxgPSMs;
	}
	
	
	public static void main(String[] args) throws IOException {
		String idedFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg/S4.PEAKS.UNMOD.rank10.pXg.BA.xCorr";
		String allFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/PEAKS/S4.UNMOD.PEAKS.csv";
		String unidedFileName = allFileName+".unided";
		
		Hashtable<String, String> idList = readPXG(new File(idedFileName));
		File allFile = new File(allFileName);
		
		BufferedReader BR = new BufferedReader(new FileReader(allFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(unidedFileName));
		
		String header = BR.readLine();
		BW.append(header);
		BW.newLine();
		String line = null;
		
		int remove = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split(",");
			String key = fields[PXG_SPEC_FILE_IDX] +"_"+fields[PXG_SCAN_IDX].split("\\:")[1];
			if(idList.get(key) == null) {
				BW.append(line);
				BW.newLine();
			} else {
				remove++;
			}
		}
		
		System.out.println(idList.size()+" were ided and "+remove+" were removed");
		BR.close();
		BW.close();
		
	}
}
