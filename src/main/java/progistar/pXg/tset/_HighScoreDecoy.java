package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class _HighScoreDecoy {

	public static void analysis() throws IOException {
		File pXgResFile = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/PEAKS_MHC_I_M009T.pxg");
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgResFile));
		String line = null;
		
		BR.readLine(); //skip header
		Hashtable<String, ArrayList<String>> collects = new Hashtable<String, ArrayList<String>>();
		Hashtable<String, Boolean> isDuplicated = new Hashtable<String, Boolean>();
		while((line = BR.readLine()) != null ) {
			String[] fields = line.split("\t");
			String fileName = fields[2]+"_"+fields[5];
			String uniqueKey = fileName+"_"+fields[19];
			if(isDuplicated.get(uniqueKey) != null) {
				continue;
			}
			isDuplicated.put(uniqueKey, true);
			ArrayList<String> collector = collects.get(fileName);
			if(collector == null) {
				collector = new ArrayList<String>();
			}
			collector.add(line);
			collects.put(fileName, collector);
		}
		
		BR.close();
		
		// enumerate
		ArrayList<String> outputs = new ArrayList<String>();
		
		collects.forEach((key, list)->{
			double maxTargetScore = 0;
			double maxDecoyScore = 0;
			double pXgTargetScore = 0;
			double pXgDecoyScore = 0;
			
			int targetCount = 0;
			int decoyCount = 0;
			for(String record : list) {
				String[] fields = record.split("\t");
				String label = fields[0];
				double score = Double.parseDouble(fields[8]);
				double rna = Math.log(Double.parseDouble(fields[36])+1)/Math.log(2);
				if(label.equalsIgnoreCase("1")) {
					maxTargetScore = Math.max(maxTargetScore, score);
					if(maxTargetScore == score) {
						pXgTargetScore = rna * score;
					}
					targetCount++;
				} else {
					maxDecoyScore = Math.max(maxDecoyScore, score);
					if(maxDecoyScore == score) {
						pXgDecoyScore = rna * score;
					}
					decoyCount++;
				}
			}
			
			outputs.add(list.get(0).split("\t")[5]+"\t"+maxTargetScore+"\t"+maxDecoyScore+"\t"+pXgTargetScore+"\t"+pXgDecoyScore+"\t"+targetCount+"\t"+decoyCount);
		});
		
		System.out.println("FileName\tMaxTargetScore\tMaxDecoyScore\tpXgTargetScore\tpXgDecoyScore\tTargetCount\tDecoyCount");
		for(int i=0; i<outputs.size(); i++) {
			System.out.println(outputs.get(i));
		}
	}
	
	public static void percolatorAnalysis () throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/pNovo3_MHC_I_M009T.target.tsv"));
		Hashtable<String, String> percolatorRes = new Hashtable<String, String>();
		String line = null;
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String psmId = fields[0];
			String peptide = fields[4];
			
			percolatorRes.put(psmId+"_"+peptide, line);
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/PEAKS_MHC_I_M009T.pxg"));
		System.out.println(BR.readLine()+"\tPSMId\tscore\tq-value\tpep\tpPept");
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[2]+fields[5]+"_"+fields[21];
			String pRecord = percolatorRes.get(key);
			if(pRecord != null) {
				System.out.println(line+"\t"+pRecord);
			}
		}
		
		BR.close();
	}
	
	public static void main(String[] args) throws IOException {
		percolatorAnalysis();
	}
}
