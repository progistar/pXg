package progistar.thirdparty.examples;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class pXgToBED {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/tmp.tsv");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter("ncMAP.bed"));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			String[] loci = fields[1].split("\\|");
			String event = fields[2];
			
			if(event.equalsIgnoreCase("unknown")) continue;
			
			for(String locus : loci) {
				String[] split = locus.split("\\:");
				String chr = split[0];
				String leftLocus = (Integer.parseInt(split[1].split("\\-")[0])-1)+"";
				String rightLocus = split[1].split("\\-")[1];
				
				BW.append(chr+"\t"+leftLocus+"\t"+rightLocus+"\t"+peptide+"_"+event);
				BW.newLine();
			}
		}
		
		
		BR.close();
		BW.close();
	}
}
