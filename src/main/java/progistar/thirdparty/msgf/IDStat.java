package progistar.thirdparty.msgf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class IDStat {

	public static void main(String[] args) throws IOException {
		String fileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/MSGF/fdr_at5/S3.fdr";
		
		// MS-GF+ index setting
		int titleIdx = 3;
		int peptideIdx = 9;
		
		
		File file = new File(fileName);
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		BR.readLine();// skip header
		
		
		// read
		Hashtable<String, String> modifiedPeptides = new Hashtable<String, String>();
		Hashtable<String, String> unmodifiedPeptides = new Hashtable<String, String>();
		Hashtable<String, String> modifiedPSMs = new Hashtable<String, String>();
		Hashtable<String, String> unmodifiedPSMs = new Hashtable<String, String>();
		
		// PTM pattern
		Pattern PTM = Pattern.compile("([+|-]+[0-9\\.]+)");
		
		// to enumerate PTMs
		Hashtable<String, String> PTMTypeChecker = new Hashtable<String, String>();
		ArrayList<String> PTMTypeList = new ArrayList<String>();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String title = fields[titleIdx];
			String peptide = fields[peptideIdx];
			
			// modified
			if(peptide.contains("+") || peptide.contains("-")) {
				modifiedPeptides.put(peptide, "");
				modifiedPSMs.put(title, peptide);
				
				Matcher matcher = PTM.matcher(peptide);
				while(matcher.find()) {
					int pos = matcher.start()-1;
					char aa = peptide.charAt(pos);
					String aaPTM = aa+matcher.group();
					
					// store PTM type at once
					if(PTMTypeChecker.get(aaPTM) == null) {
						PTMTypeList.add(aaPTM);
						PTMTypeChecker.put(aaPTM, "");
					}
				}
			}
			// unmodified
			else {
				unmodifiedPeptides.put(peptide, "");
				unmodifiedPSMs.put(title, peptide);
			}
		}
		
		BR.close();
		// PTM category
		Hashtable<String, Integer> PTMCategoryCount = new Hashtable<String, Integer>();
		modifiedPeptides.forEach((peptide, nil) -> {
			int[] counts = new int[PTMTypeList.size()];
			for(int i=0; i<PTMTypeList.size(); i++) {
				// count each PTM type
				while(peptide.contains(PTMTypeList.get(i))) {
					counts[i]++;
					int startIdx = peptide.indexOf(PTMTypeList.get(i));
					int endIdx = startIdx + PTMTypeList.get(i).length();
					String leftSub = "";
					String rightSub = "";
					
					if(startIdx != 0) {
						leftSub = peptide.substring(0, startIdx);
					}
					if(endIdx < peptide.length()) {
						rightSub = peptide.substring(endIdx);
					}
					
					peptide = leftSub+rightSub;
				}
			}
			
			String ptmList = "";
			for(int i=0; i<PTMTypeList.size(); i++) {
				if(counts[i] != 0) {
					if(ptmList.length() != 0) {
						ptmList += "&";
					}
					ptmList += counts[i]+"x"+PTMTypeList.get(i);
				}
			}
			
			Integer categoryCount = PTMCategoryCount.get(ptmList);
			if(categoryCount == null) {
				categoryCount = 0;
			}
			categoryCount++;
			PTMCategoryCount.put(ptmList, categoryCount);
		});
		
		System.out.println("UnmodifiedPSM\tModifiedPSM\tUnmodifiedPeptide\tModifiedPeptide");
		System.out.println(unmodifiedPSMs.size()+"\t"+modifiedPSMs.size()+"\t"+unmodifiedPeptides.size()+"\t"+modifiedPeptides.size());
		
		// list up
		System.out.println("Category\tDetail\tCount");
		PTMCategoryCount.forEach((category, count) -> {
			String broad = "";
			String[] fields = category.split("\\&");
			for(int i=0; i<fields.length; i++) {
				if(i!= 0) {
					broad += "&";
				}
				broad += fields[i].split("x")[1];
			}
			System.out.println(broad+"\t"+category+"\t"+count);
		});
	}
}
