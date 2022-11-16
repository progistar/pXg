package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Laumont {

	public static void main(String[] args) throws IOException {
		String pXgFile = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/pXg_Subject1.tsv";
		String laumontNovelFile = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/Laumont_Cryptic.tsv";
		String laumontPCFile = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/Laumont_Conventional.tsv";
		BufferedReader BR = new BufferedReader(new FileReader(laumontNovelFile));
		String line = null;
		
		Hashtable<String, String> ncLaumont = new Hashtable<String, String>();
		Hashtable<String, String> ncpXg = new Hashtable<String, String>();
		Hashtable<String, String> cpXg = new Hashtable<String, String>();
		Hashtable<String, String> cLaumont = new Hashtable<String, String>();
		
		BR.readLine(); // skip header
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			ncLaumont.put(peptide, line);
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(laumontPCFile));
		
		BR.readLine(); // skip header
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			cLaumont.put(peptide, line);
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(pXgFile));
		
		BR.readLine(); // skip header
		Pattern ensgPattern = Pattern.compile("ENSG+[0-9]*.[0-9]*");
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[20];
			
			// filter
			String genomicLociCountStr = fields[19];
			String geneNameCountStr = fields[30];
			String eventCountStr = fields[32];
			
			Integer genomicLociCount = Integer.parseInt(genomicLociCountStr);
			Integer geneNameCount = Integer.parseInt(geneNameCountStr);
			Integer eventCount = Integer.parseInt(eventCountStr);
			
			if(genomicLociCount > 1) continue;
			if(geneNameCount > 1) continue;
			if(eventCount > 1) continue;
			
			Matcher matcher = ensgPattern.matcher(line);
			String ensg = "ENSG@@";
			String repEnsg = "ENSG@@";
			if(matcher.find()) {
				ensg = matcher.group();
				repEnsg = ensg.split("\\.")[0];
			}
			line = line.replace(ensg, repEnsg);
			
			if(fields[36].equalsIgnoreCase("TRUE")) {
				if(cpXg.get(peptide) == null) {
					cpXg.put(peptide, line);
				}
			} else {
				if(ncpXg.get(peptide) == null) {
					ncpXg.put(peptide, line);
				}
			}
		}
		BR.close();
		// overlap
		Hashtable<String, String> cOverlap = new Hashtable<String, String>();
		Hashtable<String, String> ncOverlap = new Hashtable<String, String>();
		cpXg.forEach((peptide, l)->{
			String laumontL = cLaumont.get(peptide);
			if(laumontL != null) {
				cOverlap.put(peptide, l+"\t"+laumontL);
				System.out.println(l+"\t"+laumontL);
			}
		});
		
		ncpXg.forEach((peptide, l)->{
			String laumontL = ncLaumont.get(peptide);
			if(laumontL != null) {
				ncOverlap.put(peptide, l+"\t"+laumontL);
			}
		});
		
		System.out.println("unique cMAPs");
		System.out.println("Laumont\tpXg\tOverlap");
		System.out.println(cLaumont.size()+"\t"+cpXg.size()+"\t"+cOverlap.size());
		System.out.println("unique ncMAPs");
		System.out.println("Laumont\tpXg\tOverlap");
		System.out.println(ncLaumont.size()+"\t"+ncpXg.size()+"\t"+ncOverlap.size());
		
		
		
	}
}
