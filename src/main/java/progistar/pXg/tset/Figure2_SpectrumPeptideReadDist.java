package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Figure2_SpectrumPeptideReadDist {

	public static void main(String[] args) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/target_decoy_per_spectrum.tsv"));
		BufferedWriter BWMat = new BufferedWriter(new FileWriter("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/target_decoy_per_spectrum.matrix"));
		File[] files = new File("/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut").listFiles();
		
		
		Hashtable<String, Integer> totalPair = new Hashtable<String, Integer>();
		Hashtable<String, Integer> targetPair = new Hashtable<String, Integer>();
		Hashtable<String, Integer> decoyPair = new Hashtable<String, Integer>();
		
		for(File file : files) {
			if(!file.getName().endsWith("predfeat.pxg")) {
				continue;
			}
			
			System.out.println(file.getName());
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			BR.readLine(); // skip header
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String specId = fields[0];
				String label = fields[2];
				
				Integer totalCnt = totalPair.get(specId);
				if(totalCnt == null) {
					totalCnt = 0;
				}
				totalPair.put(specId,totalCnt+1);
				
				if(label.equalsIgnoreCase("1")) {
					Integer targetCnt = targetPair.get(specId);
					if(targetCnt == null) {
						targetCnt = 0;
					}
					targetPair.put(specId, targetCnt+1);
				}else {
					Integer decoyCnt = decoyPair.get(specId);
					if(decoyCnt == null) {
						decoyCnt = 0;
					}
					decoyPair.put(specId, decoyCnt+1);
				}
			}
			
			BR.close();
		}
		
		double[] totalMatrix = new double[1000];
		double[] targetMatrix = new double[1000];
		double[] decoyMatrix = new double[1000];
		int[][] countMatrix = new int[11][11];
		totalPair.forEach((id, cnt)->{
			totalMatrix[cnt]++;
			Integer targetCnt = targetPair.get(id);
			if(targetCnt == null) {
				targetCnt = 0;
			}
			targetMatrix[cnt] += targetCnt;
			Integer decoyCnt = decoyPair.get(id);
			if(decoyCnt == null) {
				decoyCnt = 0;
			}
			decoyMatrix[cnt] += decoyCnt;
			
			if(cnt == 710) {
				System.out.println(id);
			}
			if(targetCnt <= 10 && decoyCnt <= 10) {
				countMatrix[targetCnt][decoyCnt]++;
			}
			try {
				BW.append(id+"\t"+targetCnt+"\t"+decoyCnt);
				BW.newLine();
			}catch(IOException ioe) {
				
			}
		});
		BW.close();
		for(int i=1; i<totalMatrix.length; i++) {
			System.out.println(i+"\t"+totalMatrix[i]+"\t"+targetMatrix[i]+"\t"+decoyMatrix[i]);
		}
		
		for(int i=0; i<11; i++) {
			BWMat.append("\t"+i);
		}
		BWMat.newLine();
		
		for(int i=0; i<11; i++) {
			BWMat.append(i+"\t");
			for(int j=0; j<11; j++) {
				if(j!= 0) {
					BWMat.append("\t");
				}
				BWMat.append(countMatrix[i][j]+"");
			}
			BWMat.newLine();
		}
		BWMat.close();
	}
}
