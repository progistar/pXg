package progistar.revision;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class Q4_CompareWithPreviousResults {

	
	
	public static void main(String[] args) throws IOException {
		Hashtable<String, Hashtable<String, String>> map = loadPreviousResult(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/7.Compare/IEDB_IEAtlas_HLAligand_Cuevas_Laumont_Scull.txt"));
		Hashtable<String, Hashtable<String, String>> peaks = loadPEAKS(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/7.Compare").listFiles());
		Hashtable<String, Hashtable<String, String>> pXgs = loadpXg(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/7.Compare/PEAKS_All.pXg"), 6, 51, -1);
		Hashtable<String, Hashtable<String, String>> pXgBAs = loadpXg(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/7.Compare/PEAKS_All.pXg"), 6, 51, 49);
		Hashtable<String, Hashtable<String, String>> pXgFDRs = loadpXg(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/7.Compare/PEAKS_All_FDR.pXg"), 4, 47, -1);
		
		// compare between peaks and previous publications
		printStat(peaks, map);
		printStat(pXgs, map);
		printStat(pXgBAs, map);
		printStat(pXgFDRs, map);
		
	}
	
	public static void printStat (Hashtable<String, Hashtable<String, String>> set1, Hashtable<String, Hashtable<String, String>> set2) {
		System.out.println();
		System.out.println("Sample\tPEAKS only\tPrev. only\tCommon (coverage)\tPrev. only\tCommon (coverage)");
		set1.forEach((sample, peptides)->{
			Hashtable<String, String> previousCanonicalMap = set2.get("Canonical@"+sample);
			Hashtable<String, String> previousNoncanonicalMap = set2.get("Noncanonical@"+sample);
			
			int sizeOfPeaks = peptides.size();
			int sizeOfCanonical = previousCanonicalMap.size();
			int sizeOfNoncanonical = previousNoncanonicalMap.size();
			
			
			int[] commons = new int[2]; // 0 for canonical, 1 for noncanonical
			peptides.forEach((p, nil)->{
				if(previousCanonicalMap.get(p) != null) {
					commons[0]++;
				}
				if(previousNoncanonicalMap.get(p) != null) {
					commons[1]++;
				}
			});
			
			int peaksOnly = sizeOfPeaks - commons[0] - commons[1];
			int canonicalOnly = sizeOfCanonical - commons[0];
			int noncanonicalOnly = sizeOfNoncanonical - commons[1];
			double canonicalCoverage = ((double)commons[0]/(double)(commons[0]+canonicalOnly))*100;
			double noncanonicalCoverage = ((double)commons[1]/(double)(commons[1]+noncanonicalOnly))*100;
			
			System.out.println(sample+"\t"+peaksOnly+"\t"+canonicalOnly+"\t"+commons[0]+" ("+canonicalCoverage+")\t"+noncanonicalOnly+"\t"+commons[1]+" ("+noncanonicalCoverage+")");
		});
	}
	
	public static Hashtable<String, Hashtable<String, String>> loadpXg (File file, int peptIdx, int sampleIdx, int baIdx) throws IOException {
		Hashtable<String, Hashtable<String, String>> map = new Hashtable<String, Hashtable<String, String>>();
	
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[peptIdx].replaceAll("[+-0123456789*.]", "").replace("I", "L");
			String sampleName = fields[sampleIdx];
			
			if(sampleName.startsWith("THP")) {
				sampleName="THP1";
			} else if(sampleName.startsWith("B-LCL2")
					|| sampleName.startsWith("B-LCL3")
					|| sampleName.startsWith("B-LCL4")) {
				continue;
			}
			
			if(baIdx != -1) {
				double ba = Double.parseDouble(fields[baIdx]);
				if(ba >= 2) {
					continue;
				}
			}

			Hashtable<String, String> records = map.get(sampleName);
			
			if(records == null) {
				records = new Hashtable<String, String>();
				map.put(sampleName, records);
			}
			
			records.put(peptide, "");
		}
		
		BR.close();
		
		System.out.println();
		System.out.println("pXg");
		map.forEach((res, records)->{
			System.out.println(res+"\t"+records.size());
		});
		System.out.println("----------");
		
		return map;
	}
	
	public static Hashtable<String, Hashtable<String, String>> loadPEAKS (File[] files) throws IOException {
		Hashtable<String, Hashtable<String, String>> map = new Hashtable<String, Hashtable<String, String>>();
		
		for(File file : files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".csv")) continue;
			
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String sampleName = file.getName().split("_")[1].split("\\.")[0];
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split(",");
				String peptide = fields[3].replaceAll("[+-0123456789*.]", "").replace("I", "L");
				

				Hashtable<String, String> records = map.get(sampleName);
				
				if(records == null) {
					records = new Hashtable<String, String>();
					map.put(sampleName, records);
				}
				
				records.put(peptide, "");
			}
			
			BR.close();
		}
		System.out.println();
		System.out.println("PEAKS");
		map.forEach((res, records)->{
			System.out.println(res+"\t"+records.size());
		});
		System.out.println("----------");
		
		return map;
	}
	
	public static Hashtable<String, Hashtable<String, String>> loadPreviousResult (File file) throws IOException {
		Hashtable<String, Hashtable<String, String>> map = new Hashtable<String, Hashtable<String, String>>();
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0].replaceAll("[+-0123456789*.]", "").replace("I", "L");
			
			String tissue = fields[1];
			String resource = fields[2];
			boolean isTarget = false;
			
			if(resource.startsWith("Noncanonical")){
				resource = "Noncanonical";
				isTarget = true;
			}

			if(resource.startsWith("Canonical")) {
				resource = "Canonical";
				isTarget = true;
			}
			
			if(!isTarget) continue;
			String key = resource+"@"+tissue;
			
			Hashtable<String, String> records = map.get(key);
			
			if(records == null) {
				records = new Hashtable<String, String>();
				map.put(key, records);
			}
			
			records.put(peptide, "");
			
			
			records = map.get(tissue);
			
			if(records == null) {
				records = new Hashtable<String, String>();
				map.put(tissue, records);
			}
			
			records.put(peptide, "");
		}
		
		BR.close();
		
		
		System.out.println("Previous publications");
		map.forEach((res, records)->{
			System.out.println(res+"\t"+records.size());
		});
		System.out.println("----------");
		
		return map;
	}
	
}
