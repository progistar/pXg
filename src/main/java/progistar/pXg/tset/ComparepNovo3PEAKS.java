package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class ComparepNovo3PEAKS {

	public static void main(String[] args) throws IOException {
		File peaksFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pNovo3/DistinctIDAnalysis/S1.PEAKS.pXg");
		File pNovo3File = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pNovo3/DistinctIDAnalysis/S1.pNovo3.pXg");
		
		Hashtable<String, String> peaksRecordMap = new Hashtable<String, String>();
		Hashtable<String, String> peaksPeptideMap = new Hashtable<String, String>();
		Hashtable<String, String> pNovo3RecordMap = new Hashtable<String, String>();
		Hashtable<String, String> pNovo3PeptideMap = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(peaksFile));
		String line = null;
		
		// read PEAKS
		System.out.print(BR.readLine()); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String scan = fields[1].split("\\.")[0] +"."+ fields[4].split("\\:")[1];
			peaksRecordMap.put(scan, line);
			peaksPeptideMap.put(scan, fields[20]);
		}
		
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(pNovo3File));
		System.out.print("\t"+BR.readLine());
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String scan = fields[1].split("\\.")[0] + "." + fields[1].split("\\.")[1];
			pNovo3RecordMap.put(scan, line);
			pNovo3PeptideMap.put(scan, fields[6]);
		}
		
		BR.close();
		
		peaksPeptideMap.forEach((scan, peptide) -> {
			String pNovo3Peptide = pNovo3PeptideMap.get(scan);
			if(pNovo3Peptide != null && !pNovo3Peptide.equalsIgnoreCase(peptide)) {
				System.out.println(peaksRecordMap.get(scan)+"\t"+pNovo3RecordMap.get(scan));
			}
		});
		
	}
}
