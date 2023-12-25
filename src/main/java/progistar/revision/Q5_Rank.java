package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

public class Q5_Rank {

	public static void main(String[] agrs) throws IOException {
		
		File[] pXgFiles = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut").listFiles();
		File[] pinFiles = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features").listFiles();
		
		String out = "/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features_rank/";
		
		Hashtable<String, File> pXgFileHash = new Hashtable<String, File>();
		Hashtable<String, File> pinFileHash = new Hashtable<String, File>();
		
		for(File file : pXgFiles) {
			if(file.getName().startsWith(".")) continue;
			if(file.getName().endsWith(".pxg")) {
				pXgFileHash.put(file.getName().split("\\.")[0], file);
			}
		}
		
		for(File file : pinFiles) {
			if(file.getName().startsWith(".")) continue;
			if(file.getName().endsWith(".feat2.pin")) {
				pinFileHash.put(file.getName().split("\\.")[0], file);
			}
		}
		
		
		Iterator<String> samples= (Iterator<String>) pXgFileHash.keys();
		while(samples.hasNext()) {
			String sample = samples.next();
			File pXgFile = pXgFileHash.get(sample);
			File pinFile = pinFileHash.get(sample);
			BufferedReader BR = new BufferedReader(new FileReader(pXgFile));
			String line = null;
			
			int rankIdx = 22;
			ArrayList<Integer> ranks = new ArrayList<Integer>();
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String rank = line.split("\t")[rankIdx];
				ranks.add(Integer.parseInt(rank));
			}
			
			BR.close();
			
			BR = new BufferedReader(new FileReader(pinFile));
			String header = BR.readLine();
			ArrayList<String> records = new ArrayList<String>();
			
			while((line = BR.readLine()) != null) {
				records.add(line);
			}
			
			BR.close();
			
			for(int i=1; i<=10; i++) {
				int rank = i;
				BufferedWriter BW = new BufferedWriter(new FileWriter(out+sample+".rank"+rank+".pin"));
				BW.append(header);
				BW.newLine();
				for(int j=0; j<records.size(); j++) {
					if(ranks.get(j) <= rank) {
						BW.append(records.get(j));
						BW.newLine();
					}
				}
				
				BW.close();
			}
			
			
		}
	}
	
	
}
