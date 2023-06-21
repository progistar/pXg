package progistar.thridparty.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class SummaryOfResultCategory {

	public static void main(String[] args) throws IOException {
		File[] files = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features").listFiles();
		System.out.println("Sample\tMatchedMSMS\tcTarget\tcDecoy\tncTarget\tncDecoy");
		for(File file : files) {
			if(!file.getName().endsWith("2.BA")) {
				continue;
			}
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String sampleName = file.getName().split("\\.")[0];
			String line = null;
			BR.readLine(); // header
			
			int cTarget = 0;
			int cDecoy = 0;
			int ncTarget = 0;
			int ncDecoy = 0;
			Hashtable<String, String> scans = new Hashtable<String, String>();
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String scanId = fields[0];
				String label = fields[2];
				String isCanonical = fields[41];
				Double binderRank = Double.parseDouble(fields[49]);
				
				if(binderRank < 0) {
					continue;
				}
				
				scans.put(scanId, "");
				
				if(label.equalsIgnoreCase("1")) {
					if(isCanonical.equalsIgnoreCase("true")) {
						cTarget++;
					} else {
						ncTarget++;
					}
				} else {
					if(isCanonical.equalsIgnoreCase("true")) {
						cDecoy++;
					} else {
						ncDecoy++;
					}
				}
			}
			
			System.out.println(sampleName+"\t"+scans.size()+"\t"+cTarget+"\t"+cDecoy+"\t"+ncTarget+"\t"+ncDecoy);
			
			BR.close();
		}
	}
}
