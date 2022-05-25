package progistar.thirdparty.iedb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Mapping {

	public static void main(String[] args) throws IOException {
		File iedbFile = new File("/Users/gistar/projects/pXg/IEDB/IEDB.host");
		File pXgFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/IEDBInput/S4.PEAKS.UNMOD.rank10.pXg.BA.xCorr.fdr");
		
		Hashtable<String, String> iedbMapper = loadIEDB(iedbFile);
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgFile.getAbsolutePath()+".iedb"));
		String line = null;
		
		// skip header
		BW.append(BR.readLine());
		BW.append("\tIEDB");
		BW.newLine();
		
		Hashtable<String, String> duplicates = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[20];
			String isBA = fields[43];
			String isCanonical = fields[36];
			String eventCount = fields[32];
			String geneIDCount = fields[28];
			String genomicLoci = fields[19];
			
			if(!eventCount.equalsIgnoreCase("1")) {
				continue;
			}
			
			if(!(geneIDCount.equalsIgnoreCase("1") || geneIDCount.equalsIgnoreCase("0") || geneIDCount.equalsIgnoreCase("-"))) {
				continue;
			}
			
			if(!(genomicLoci.equalsIgnoreCase("1") || genomicLoci.equalsIgnoreCase("0") || genomicLoci.equalsIgnoreCase("-"))) {
				continue;
			}
			
			if(duplicates.get(peptide) != null) {
				continue;
			}
			
			duplicates.put(peptide, "");
			
			if(isCanonical.equalsIgnoreCase("TRUE")) {
				continue;
			}
			
			if(isBA.equalsIgnoreCase("NB")) {
				continue;
			}
			
			String iedbInfo = iedbMapper.get(peptide);
			if(iedbInfo == null) {
				iedbInfo = "NA";
			}
			
			BW.append(line+"\t"+iedbInfo);
			BW.newLine();
		}
		BW.close();
		BR.close();
	}
	
	public static Hashtable<String, String> loadIEDB (File file) throws IOException {
		Hashtable<String, String> iedb = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			iedb.put(peptide, line);
		}
		
		BR.close();
		
		return iedb;
	}
}
