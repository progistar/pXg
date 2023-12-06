package progistar.revision;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Q1_MaxQualTest {

	public static void main(String[] args) throws IOException {
		Hashtable<String, String> avgQualMap = loadQualTable(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/2.PhredMax/B-LCL1.pxg"));
		Hashtable<String, String> maxQualMap = loadQualTable(new File("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/2.PhredMax/B-LCL1.phredMAX.pxg.pXg"));
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("B-LCL1.mean.max.tsv"));
		
		BW.append("id\tReads\tMean\tMax");
		BW.newLine();
		avgQualMap.forEach((key, meanQual)->{
			try {
				String maxQual = maxQualMap.get(key);
				if(maxQual != null) {
					BW.append(key+"\t"+meanQual+"\t"+maxQual);
					BW.newLine();
				}
			}catch(IOException ioe) {
				
			}
		});
		
		BW.close();
	}
	
	public static Hashtable<String, String> loadQualTable(File pXgFile) throws IOException {
		Hashtable<String, String> qualMap = new Hashtable<String, String>();
		BufferedReader BR = new BufferedReader(new FileReader(pXgFile));
		String line = null;
		
		BR.readLine();// skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[0]+"_"+fields[24]+"_"+fields[25]+"_"+fields[26]+"_"+fields[27]+"_"+fields[28]+"\t"+fields[39];
			String qual = fields[40];
			
			if(qualMap.get(key) == null) {
				qualMap.put(key, qual);
			} else {
				System.out.println(key+" is duplicated!");
			}
		}
		
		BR.close();
		
		return qualMap;
	}
}
