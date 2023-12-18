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
		
		String feat = "protonly";
		String type = "psm";
		
		boolean[] selectedSet = {
				true, // B-LCL1
				true, // B-LCL2
				true, // B-LCL3
				true, // B-LCL4
				true, // DOHH2
				true, // HBL1
				true, // SUDHL4
				true, // THP1-1
				true, // THP1-2
				true, // THP1-3
				
				false, // DI2T
				false, // DI5T
				false, // IN19T
				false, // IN26T
				false, // IN81T
				false, // R1_IN403T
				false, // R2_IN403T
				false, // IN407T
				false, // IN506T
				false, // IN524T
				false, // IN525T
				false, // IN526T
				false, // IN529T
				false, // R1_IN532T
				false, // R2_IN532T
				false, // M004T
				false, // M009T
				
				false, // DI2N
				false, // DI5N
				false, // IN19N
				false, // IN26N
				false, // IN81N
				false, // IN403N
				false, // IN407N
				false, // IN506N
				false, // IN524N
				false, // IN525N
				false, // IN526N
				false, // IN529N
				false, // IN532N
				false, // M004N
				false // M009N
		};
		
		
		String[] sampleListDecoy = {
				"B-LCL1."+feat+".decoy."+type,
				"B-LCL2."+feat+".decoy."+type,
				"B-LCL3."+feat+".decoy."+type,
				"B-LCL4."+feat+".decoy."+type,
				"HBL1."+feat+".decoy."+type,
				"DOHH2."+feat+".decoy."+type,
				"SUDHL4."+feat+".decoy."+type,
				"THP1-1."+feat+".decoy."+type,
				"THP1-2."+feat+".decoy."+type,
				"THP1-3."+feat+".decoy."+type,
				
				"DI2T."+feat+".decoy."+type, // DI2T
				"DI5T."+feat+".decoy."+type, // DI5T
				"IN19T."+feat+".decoy."+type, // IN19T
				"IN26T."+feat+".decoy."+type, // IN26T
				"IN81T."+feat+".decoy."+type, // IN81T
				"R1_IN403T."+feat+".decoy."+type, // R1_IN403T
				"R2_IN403T."+feat+".decoy."+type, // R2_IN403T
				"IN407T."+feat+".decoy."+type, // IN407T
				"IN506T."+feat+".decoy."+type, // IN506T
				"IN524T."+feat+".decoy."+type, // IN524T
				"IN525T."+feat+".decoy."+type, // IN525T
				"IN526T."+feat+".decoy."+type, // IN526T
				"IN529T."+feat+".decoy."+type, // IN529T
				"R1_IN532T."+feat+".decoy."+type, // R1_IN532T
				"R2_IN532T."+feat+".decoy."+type, // R2_IN532T
				"M004T."+feat+".decoy."+type, // M004T
				"M009T."+feat+".decoy."+type // M009T
		};
		
		String[] sampleListTarget = {
				"B-LCL1."+feat+".target."+type,
				"B-LCL2."+feat+".target."+type,
				"B-LCL3."+feat+".target."+type,
				"B-LCL4."+feat+".target."+type,
				"HBL1."+feat+".target."+type,
				"DOHH2."+feat+".target."+type,
				"SUDHL4."+feat+".target."+type,
				"THP1-1."+feat+".target."+type,
				"THP1-2."+feat+".target."+type,
				"THP1-3."+feat+".target."+type,
				
				"DI2T."+feat+".target."+type, // DI2T
				"DI5T."+feat+".target."+type, // DI5T
				"IN19T."+feat+".target."+type, // IN19T
				"IN26T."+feat+".target."+type, // IN26T
				"IN81T."+feat+".target."+type, // IN81T
				"R1_IN403T."+feat+".target."+type, // R1_IN403T
				"R2_IN403T."+feat+".target."+type, // R2_IN403T
				"IN407T."+feat+".target."+type, // IN407T
				"IN506T."+feat+".target."+type, // IN506T
				"IN524T."+feat+".target."+type, // IN524T
				"IN525T."+feat+".target."+type, // IN525T
				"IN526T."+feat+".target."+type, // IN526T
				"IN529T."+feat+".target."+type, // IN529T
				"R1_IN532T."+feat+".target."+type, // R1_IN532T
				"R2_IN532T."+feat+".target."+type, // R2_IN532T
				"M004T."+feat+".target."+type, // M004T
				"M009T."+feat+".target."+type // M009T
		};
		
		String[] sampleListTargetXML = {
				"B-LCL1."+feat+".xml",
				"B-LCL2."+feat+".xml",
				"B-LCL3."+feat+".xml",
				"B-LCL4."+feat+".xml",
				"HBL1."+feat+".xml",
				"DOHH2."+feat+".xml",
				"SUDHL4."+feat+".xml",
				"THP1-1."+feat+".xml",
				"THP1-2."+feat+".xml",
				"THP1-3."+feat+".xml",

				"DI2T."+feat+".xml", // DI2T
				"DI5T."+feat+".xml", // DI5T
				"IN19T."+feat+".xml", // IN19T
				"IN26T."+feat+".xml", // IN26T
				"IN81T."+feat+".xml", // IN81T
				"R1_IN403T."+feat+".xml", // R1_IN403T
				"R2_IN403T."+feat+".xml", // R2_IN403T
				"IN407T."+feat+".xml", // IN407T
				"IN506T."+feat+".xml", // IN506T
				"IN524T."+feat+".xml", // IN524T
				"IN525T."+feat+".xml", // IN525T
				"IN526T."+feat+".xml", // IN526T
				"IN529T."+feat+".xml", // IN529T
				"R1_IN532T."+feat+".xml", // R1_IN532T
				"R2_IN532T."+feat+".xml", // R2_IN532T
				"M004T."+feat+".xml", // M004T
				"M009T."+feat+".xml" // M009T
		};
		
		String[] pxgList = {
				"B-LCL1.pxg",
				"B-LCL2.pxg",
				"B-LCL3.pxg",
				"B-LCL4.pxg",
				"HBL1.pxg",
				"DOHH2.pxg",
				"SUDHL4.pxg",
				"THP1-1.pxg",
				"THP1-2.pxg",
				"THP1-3.pxg",
				
				"DI2T.pxg", // DI2T
				"DI5T.pxg", // DI5T
				"IN19T.pxg", // IN19T
				"IN26T.pxg", // IN26T
				"IN81T.pxg", // IN81T
				"R1_IN403T.pxg", // R1_IN403T
				"R2_IN403T.pxg", // R2_IN403T
				"IN407T.pxg", // IN407T
				"IN506T.pxg", // IN506T
				"IN524T.pxg", // IN524T
				"IN525T.pxg", // IN525T
				"IN526T.pxg", // IN526T
				"IN529T.pxg", // IN529T
				"R1_IN532T.pxg", // R1_IN532T
				"R2_IN532T.pxg", // R2_IN532T
				"M004T.pxg", // M004T
				"M009T.pxg", // M009T
				
				"DI2N.pxg",
				"DI5N.pxg",
				"IN19N.pxg",
				"IN26N.pxg",
				"IN81N.pxg",
				"IN403N.pxg",
				"IN407N.pxg",
				"IN506N.pxg",
				"IN524N.pxg",
				"IN525N.pxg",
				"IN526N.pxg",
				"IN529N.pxg",
				"IN532N.pxg",
				"M004N.pxg",
				"M009N.pxg"
		};
		
		String[] BAList = {
				"B-LCL1.netmhcpan.xls",
				"B-LCL2.netmhcpan.xls",
				"B-LCL3.netmhcpan.xls",
				"B-LCL4.netmhcpan.xls",
				"HBL1.netmhcpan.xls",
				"DOHH2.netmhcpan.xls",
				"SUDHL4.netmhcpan.xls",
				"THP1-1.netmhcpan.xls",
				"THP1-2.netmhcpan.xls",
				"THP1-3.netmhcpan.xls",
				
				"DI2T.netmhcpan.xls", // DI2T
				"DI5T.netmhcpan.xls", // DI5T
				"IN19T.netmhcpan.xls", // IN19T
				"IN26T.netmhcpan.xls", // IN26T
				"IN81T.netmhcpan.xls", // IN81T
				"R1_IN403T.netmhcpan.xls", // R1_IN403T
				"R2_IN403T.netmhcpan.xls", // R2_IN403T
				"IN407T.netmhcpan.xls", // IN407T
				"IN506T.netmhcpan.xls", // IN506T
				"IN524T.netmhcpan.xls", // IN524T
				"IN525T.netmhcpan.xls", // IN525T
				"IN526T.netmhcpan.xls", // IN526T
				"IN529T.netmhcpan.xls", // IN529T
				"R1_IN532T.netmhcpan.xls", // R1_IN532T
				"R2_IN532T.netmhcpan.xls", // R2_IN532T
				"M004T.netmhcpan.xls", // M004T
				"M009T.netmhcpan.xls", // M009T
				
				"DI2N.netmhcpan.xls", // DI2N
				"DI5N.netmhcpan.xls", // DI5N
				"IN19N.netmhcpan.xls", // IN19N
				"IN26N.netmhcpan.xls", // IN26N
				"IN81N.netmhcpan.xls", // IN81N
				"IN403N.netmhcpan.xls", // R1_IN403N
				"IN407N.netmhcpan.xls", // IN407N
				"IN506N.netmhcpan.xls", // IN506N
				"IN524N.netmhcpan.xls", // IN524N
				"IN525N.netmhcpan.xls", // IN525N
				"IN526N.netmhcpan.xls", // IN526N
				"IN529N.netmhcpan.xls", // IN529N
				"IN532N.netmhcpan.xls", // R1_IN532N
				"M004N.netmhcpan.xls", // M004N
				"M009N.netmhcpan.xls" // M009N
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
				"THP1-3",
				
				"DI2T", // DI2T
				"DI5T", // DI5T
				"IN19T", // IN19T
				"IN26T", // IN26T
				"IN81T", // IN81T
				"R1_IN403T", // R1_IN403T
				"R2_IN403T", // R2_IN403T
				"IN407T", // IN407T
				"IN506T", // IN506T
				"IN524T", // IN524T
				"IN525T", // IN525T
				"IN526T", // IN526T
				"IN529T", // IN529T
				"R1_IN532T", // R1_IN532T
				"R2_IN532T", // R2_IN532T
				"M004T", // M004T
				"M009T", // M009T
				
				"DI2N", // DI2N
				"DI5N", // DI5N
				"IN19N", // IN19N
				"IN26N", // IN26N
				"IN81N", // IN81N
				"IN403N", // IN403N
				"IN407N", // IN407N
				"IN506N", // IN506N
				"IN524N", // IN524N
				"IN525N", // IN525N
				"IN526N", // IN526N
				"IN529N", // IN529N
				"IN532N", // IN532N
				"M004N", // M004N
				"M009N" // M009N
		};
		
		String[] outList = {
				"_B-LCL1."+feat+".BA",
				"_B-LCL2."+feat+".BA",
				"_B-LCL3."+feat+".BA",
				"_B-LCL4."+feat+".BA",
				"_HBL1."+feat+".BA",
				"_DOHH2."+feat+".BA",
				"_SUDHL4."+feat+".BA",
				"_THP1-1."+feat+".BA",
				"_THP1-2."+feat+".BA",
				"_THP1-3."+feat+".BA",
				
				"DI2T."+feat+".BA", // DI2T
				"DI5T."+feat+".BA", // DI5T
				"IN19T."+feat+".BA", // IN19T
				"IN26T."+feat+".BA", // IN26T
				"IN81T."+feat+".BA", // IN81T
				"R1_IN403T."+feat+".BA", // R1_IN403T
				"R2_IN403T."+feat+".BA", // R2_IN403T
				"IN407T."+feat+".BA", // IN407T
				"IN506T."+feat+".BA", // IN506T
				"IN524T."+feat+".BA", // IN524T
				"IN525T."+feat+".BA", // IN525T
				"IN526T."+feat+".BA", // IN526T
				"IN529T."+feat+".BA", // IN529T
				"R1_IN532T."+feat+".BA", // R1_IN532T
				"R2_IN532T."+feat+".BA", // R2_IN532T
				"M004T."+feat+".BA", // M004T
				"M009T."+feat+".BA", // M009T
				
				
				"DI2N."+feat+".BA", // DI2N
				"DI5N."+feat+".BA", // DI5N
				"IN19N."+feat+".BA", // IN19N
				"IN26N."+feat+".BA", // IN26N
				"IN81N."+feat+".BA", // IN81N
				"IN403N."+feat+".BA", // R1_IN403N
				"IN407N."+feat+".BA", // IN407N
				"IN506N."+feat+".BA", // IN506N
				"IN524N."+feat+".BA", // IN524N
				"IN525N."+feat+".BA", // IN525N
				"IN526N."+feat+".BA", // IN526N
				"IN529N."+feat+".BA", // IN529N
				"IN532N."+feat+".BA", // IN532N
				"M004N."+feat+".BA", // M004N
				"M009N."+feat+".BA" // M009N
		};
		
		BufferedWriter BWAll = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/_all."+feat+".BA"));
		
		int cnt = 0;
		
		for(int i=0; i<sampleListTarget.length; i++) {
			
			if(selectedSet[i] == false) {
				continue;
			}

			BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/"+outList[i]));
			
			BufferedReader BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/"+sampleListTarget[i]));
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
			BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/"+sampleListDecoy[i]));
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String psmId = fields[0];
				String peptide = fields[4];
				
				percolatorRes.put(psmId+"_"+peptide, line);
			}
			
			System.out.println("Percolator result: "+percolatorRes.size());
			BR.close();
			NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/"+BAList[i]);
			
			BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/"+pxgList[i]+".predfeat."+feat+".pin"));
			String[] pinHeader = BR.readLine().split("\t");
			int specIdIdx = -1;
			int SAIdx = -1;
			int log2ReadIdx = -1;
			int bestDeltaRTIdx = -1;
			int log2MeanQScoreIdx = -1;
			int mLog2BestELRankIdx = -1;
			int peptideIdx = -1;
			int proteinsIdx = -1;
			
			Hashtable<String, String> pinFeatureMapper = new Hashtable<String, String>();
			
			for(int j=0; j<pinHeader.length; j++) {
				if(pinHeader[j].equalsIgnoreCase("SpecId")) {
					specIdIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("Log2Reads")) {
					log2ReadIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("SA")) {
					SAIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("BestDeltaRT")) {
					bestDeltaRTIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("Log2MeanQScore")) {
					log2MeanQScoreIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("mLog2BestELRank")) {
					mLog2BestELRankIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("Peptide")) {
					peptideIdx = j;
				} else if(pinHeader[j].equalsIgnoreCase("Proteins")) {
					proteinsIdx = j;
				}
			}
			
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String genomicId = fields[proteinsIdx];
				if(genomicId.contains("_")) {
					genomicId = genomicId.split("_")[1];
				}
				//String key = fields[specIdIdx]+"_"+genomicId+"_"+fields[peptideIdx];
				//pinFeatureMapper.put(key, fields[log2ReadIdx]+"\t"+fields[SAIdx]+"\t"+
				//fields[bestDeltaRTIdx]+"\t"+fields[log2MeanQScoreIdx]);
				
				String key = fields[specIdIdx]+"_"+genomicId+"_"+fields[peptideIdx];
				pinFeatureMapper.put(key, fields[SAIdx]+"\t"+fields[bestDeltaRTIdx]);
			}
			
			BR.close();
			
			
			BR = new BufferedReader(new FileReader("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/"+pxgList[i]));
			//String header = BR.readLine()+"\tLog2Reads\tSA\tBestDeltaRT\tLog2MeanQScore\tmLog2BestELRank\tpercolator_score\tq-value\tpep\tEL_Rank\tBestHLAType\tSample\tIsSelectedTarget";
			String header = BR.readLine()+"\tSA\tBestDeltaRT\tmLog2BestELRank\tpercolator_score\tq-value\tpep\tEL_Rank\tBestHLAType\tSample\tIsSelectedTarget";
			BW.append(header);
			BW.newLine();
			
			if(cnt == 0) {
				BWAll.append(header);
				BWAll.newLine();
			}
			
			Hashtable<String,String> surviveMapper = parseXML("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut_features/"+sampleListTargetXML[i]);
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String psmId = fields[0];
				String peptide = fields[24];
				String genomicId = fields[1];
				String key = psmId+"_"+peptide;
				String pRecord = percolatorRes.get(key);
				String key2 = psmId+"_"+genomicId+"_"+peptide;
				String key3 = psmId+"_"+genomicId;
				String isSurvived = surviveMapper.get(key3);
				if(isSurvived == null) {
					isSurvived = "-1";
				} else {
					isSurvived = "1";
				}
				
				if(pRecord != null) {
					fields = pRecord.split("\t");
					pRecord = fields[1]+"\t"+fields[2]+"\t"+fields[3];
					BW.append(line).append("\t").append(pinFeatureMapper.get(key2)).append("\t").append(pRecord).append("\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()+"\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestHLAType()).append("\t")
					.append(sampleNames[i]).append("\t").append(isSurvived);
					BW.newLine();
					
					BWAll.append(line).append("\t").append(pinFeatureMapper.get(key2)).append("\t").append(pRecord).append("\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()+"\t")
					.append(netMHCpanResult.peptideToRecord.get(peptide).getBestHLAType()).append("\t")
					.append(sampleNames[i]).append("\t").append(isSurvived);
					BWAll.newLine();
					cnt++;
				}
			}
			
			BR.close();
			BW.close();
		}
		
		BWAll.close();
		
	}
	
	public static Hashtable<String, String> parseXML (String fileName) throws IOException {
		// specId_genomicId
		Hashtable<String, String> mapper = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		String line = null;
		
		String specId = null;
		String genomicId = null;
		while((line = BR.readLine()) != null) {
			if(line.contains("p:psm_id=")) {
				specId = line.split("\"")[1];
			}else if(line.contains("<protein_id>")) {
				genomicId = line.replace("<protein_id>","").replace("</protein_id>", "").trim();
				mapper.put(specId+"_"+genomicId,"");
			}
		}
		BR.close();
		
		return mapper;
	}
	
	public static void main(String[] args) throws IOException {
		//analysisTD("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/S4.RAW.PEAKS.nocut.pxg");
		percolatorAnalysis();
	}
}
