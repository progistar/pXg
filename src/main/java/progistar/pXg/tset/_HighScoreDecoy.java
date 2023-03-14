package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.thirdparty.netMHCpan.NetMHCpanParser;
import progistar.thirdparty.netMHCpan.NetMHCpanResult;

public class _HighScoreDecoy {

	public static void analysisTD(String fileName) throws IOException {
		File pXgResFile = new File(fileName);
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgResFile));
		String line = null;
		
		BR.readLine(); //skip header
		Hashtable<String, ArrayList<String>> collects = new Hashtable<String, ArrayList<String>>();
		while((line = BR.readLine()) != null ) {
			String[] fields = line.split("\t");
			String uniqueKey = fields[0];
			ArrayList<String> collector = collects.get(uniqueKey);
			if(collector == null) {
				collector = new ArrayList<String>();
			}
			collector.add(line);
			collects.put(uniqueKey, collector);
		}
		
		BR.close();
		
		// enumerate
		ArrayList<String> outputs = new ArrayList<String>();
		
		collects.forEach((key, list)->{
			double maxTargetScore = 0;
			double maxDecoyScore = 0;
			double targetRNA = 0;
			double decoyRNA = 0;
			
			for(String record : list) {
				String[] fields = record.split("\t");
				String label = fields[1];
				double score = Double.parseDouble(fields[9]);
				double rna = Math.log(Double.parseDouble(fields[37])+1)/Math.log(2);
				if(label.equalsIgnoreCase("1")) {
					maxTargetScore = Math.max(maxTargetScore, score);
					if(maxTargetScore == score) {
						targetRNA = Math.max(targetRNA, rna);
					}
				} else {
					maxDecoyScore = Math.max(maxDecoyScore, score);
					if(maxDecoyScore == score) {
						decoyRNA = Math.max(decoyRNA, rna);
					}
				}
			}
			String category = "";
			if(maxTargetScore >= maxDecoyScore) {
				if(maxDecoyScore == 0) {
					category = "Target only";
				} else {
					category = "Both but target win";
				}
			} else {
				if(maxTargetScore == 0) {
					category = "Decoy only";
				} else {
					category = "Both but decoy win";
				}
			}
			
			outputs.add(key+"\t"+maxTargetScore+"\t"+maxDecoyScore+"\t"+targetRNA+"\t"+decoyRNA+"\t"+category);
		});
		
		System.out.println("FileID\tMaxTargetScore\tMaxDecoyScore\tMaxTargetRNA\tMaxDecoyRNA\tCategory");
		for(int i=0; i<outputs.size(); i++) {
			System.out.println(outputs.get(i));
		}
	}
	
	public static void percolatorAnalysis () throws IOException {
		String sample = "4";
		BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S"+sample+".RAW.PEAKS.target.BA"));
		
		BufferedReader BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S"+sample+".RAW.PEAKS.target.tsv"));
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
		
		NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/tools/netMHCpan4.1/netMHCpan-4.1/nocut/peptide"+sample+".xls");
		
		BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S"+sample+".RAW.PEAKS.nocut.pxg"));
		BW.append(BR.readLine()+"\tpercolator_score\tq-value\tpep\t"+netMHCpanResult.getHeader());
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String psmId = fields[0];
			String peptide = fields[22];
			String key = psmId+"_"+peptide;
			String pRecord = percolatorRes.get(key);
			
			if(pRecord != null) {
				fields = pRecord.split("\t");
				pRecord = fields[1]+"\t"+fields[2]+"\t"+fields[3];
				BW.append(line).append("\t").append(pRecord).append("\t").append(netMHCpanResult.getHLATyping(peptide));
				BW.newLine();
			}
		}
		
		BR.close();
		BW.close();
	}
	
	public static void main(String[] args) throws IOException {
		//analysisTD("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S4.RAW.PEAKS.nocut.pxg");
		percolatorAnalysis();
	}
}
