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
		
		String feat = "feat1";
		String type = "psm";
		
		String[] sampleListDecoy = {
				"BLCL1."+feat+".decoy."+type,
				"BLCL2."+feat+".decoy."+type,
				"BLCL3."+feat+".decoy."+type,
				"BLCL4."+feat+".decoy."+type,
				"HBL1."+feat+".decoy."+type,
				"DOHH2."+feat+".decoy."+type,
				"SUDHL4."+feat+".decoy."+type,
				"THP1_1."+feat+".decoy."+type,
				"THP1_2."+feat+".decoy."+type,
				"THP1_3."+feat+".decoy."+type
		};
		
		String[] sampleList = {
				"BLCL1."+feat+".target."+type,
				"BLCL2."+feat+".target."+type,
				"BLCL3."+feat+".target."+type,
				"BLCL4."+feat+".target."+type,
				"HBL1."+feat+".target."+type,
				"DOHH2."+feat+".target."+type,
				"SUDHL4."+feat+".target."+type,
				"THP1_1."+feat+".target."+type,
				"THP1_2."+feat+".target."+type,
				"THP1_3."+feat+".target."+type
		};
		
		
		String[] pxgList = {
				"PEAKS_BLCL1.nocut.pxg",
				"PEAKS_BLCL2.nocut.pxg",
				"PEAKS_BLCL3.nocut.pxg",
				"PEAKS_BLCL4.nocut.pxg",
				"PEAKS_HBL1.nocut.pxg",
				"PEAKS_DOHH2.nocut.pxg",
				"PEAKS_SUDHL4.nocut.pxg",
				"PEAKS_THP1_1.nocut.pxg",
				"PEAKS_THP1_2.nocut.pxg",
				"PEAKS_THP1_3.nocut.pxg"
		};
		
		String[] BAList = {
				"peptide1.xls",
				"peptide2.xls",
				"peptide3.xls",
				"peptide4.xls",
				"HBL1.netmhcpan.xls",
				"DOHH2.netmhcpan.xls",
				"SUDHL4.netmhcpan.xls",
				"peptide_THP1_1.xls",
				"peptide_THP1_2.xls",
				"peptide_THP1_3.xls"
		};
		
		String[] sampleNames = {
				"B-LCL1",
				"B-LCL2",
				"B-LCL3",
				"B-LCL4",
				"HBL1",
				"DOHH2",
				"SUDHL4",
				"THP1-1",
				"THP1-2",
				"THP1-3"
		};
		
		String[] outList = {
				"_S1."+feat+".BA",
				"_S2."+feat+".BA",
				"_S3."+feat+".BA",
				"_S4."+feat+".BA",
				"_HBL1."+feat+".BA",
				"_DOHH2."+feat+".BA",
				"_SUDHL4."+feat+".BA",
				"_THP1_1."+feat+".BA",
				"_THP1_2."+feat+".BA",
				"_THP1_3."+feat+".BA"
		};
		
		String[] sampleList2 = {
				"HBL1."+feat+".target."+type,
				"DOHH2."+feat+".target."+type,
				"SUDHL4."+feat+".target."+type
		};
		
		String[] sampleList3 = {
				"HBL1."+feat+".decoy."+type,
				"DOHH2."+feat+".decoy."+type,
				"SUDHL4."+feat+".decoy."+type
		};
		
		
		String[] pxgList2 = {
				"PEAKS_HBL1.nocut.pxg",
				"PEAKS_DOHH2.nocut.pxg",
				"PEAKS_SUDHL4.nocut.pxg"
		};
		
		String[] BAList2 = {
				"HBL1.netmhcpan.xls",
				"DOHH2.netmhcpan.xls",
				"SUDHL4.netmhcpan.xls"
		};
		
		String[] sampleNames2 = {
				"HBL1",
				"DOHH2",
				"SUDHL4"
		};
		
		String[] outList2 = {
				"_HBL1."+feat+".BA",
				"_DOHH2."+feat+".BA",
				"_SUDHL4."+feat+".BA"
		};
		
		BufferedWriter BWAll = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads_features/_all."+feat+".BA"));
		
		int cnt = 0;
		
		for(int i=0; i<sampleList2.length; i++) {

			BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads_features/"+outList2[i]));
			
			BufferedReader BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads_features/"+sampleList2[i]));
			Hashtable<String, String> percolatorRes = new Hashtable<String, String>();
			String line = null;
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String psmId = fields[0];
				String peptide = fields[4];
				
				percolatorRes.put(psmId+"_"+peptide, line);
			}
			
			System.out.println("Percolator result: "+percolatorRes.size());
			BR.close();
			// for decoy
			BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads_features/"+sampleList3[i]));
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String psmId = fields[0];
				String peptide = fields[4];
				
				percolatorRes.put(psmId+"_"+peptide, line);
			}
			
			System.out.println("Percolator result: "+percolatorRes.size());
			BR.close();
			NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/"+BAList2[i]);
			
			BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_fullreads/"+pxgList2[i]));
			String header = BR.readLine()+"\tpercolator_score\tq-value\tpep\tEL_Rank\tBestHLAType\tSample";
			BW.append(header);
			BW.newLine();
			
			if(cnt == 0) {
				BWAll.append(header);
				BWAll.newLine();
			}
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String psmId = fields[0];
				String peptide = fields[23];
				String key = psmId+"_"+peptide;
				String pRecord = percolatorRes.get(key);
				
				if(pRecord != null) {
					fields = pRecord.split("\t");
					pRecord = fields[1]+"\t"+fields[2]+"\t"+fields[3];
					BW.append(line).append("\t").append(pRecord).append("\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()+"\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestHLAType()).append("\t")
					.append(sampleNames2[i]);
					BW.newLine();
					
					BWAll.append(line).append("\t").append(pRecord).append("\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()+"\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestHLAType()).append("\t")
					.append(sampleNames2[i]);
					BWAll.newLine();
					cnt++;
				}
			}
			
			BR.close();
			BW.close();
		}
		
		BWAll.close();
		
	}
	
	public static void main(String[] args) throws IOException {
		//analysisTD("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S4.RAW.PEAKS.nocut.pxg");
		percolatorAnalysis();
	}
}
