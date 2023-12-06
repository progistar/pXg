package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Q2_1_NetMHCpanInput {
	
	public static void main(String[] args) throws IOException {
		
		File[] files = {new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.target.psm"),
				new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/4.ReadTranslation/B_LCL1.decoy.psm")};
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("B_LCL1.target.decoy.peptide"));
		Hashtable<String, String> map = new Hashtable<String, String>();
		for(File file : files) {

			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			BR.readLine();//
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[4].split("\\.")[1];
				if(map.get(peptide) == null) {
					BW.append(peptide);
					BW.newLine();
					map.put(peptide, "");
				}
			}
			
			BR.close();
		}
		
		BW.close();
		System.out.println(map.size() +" were detected.");
	}

}
