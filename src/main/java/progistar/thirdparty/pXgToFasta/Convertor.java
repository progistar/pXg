package progistar.thirdparty.pXgToFasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Convertor {

	public static void main(String[] args) throws IOException {
		// Peptide Event
		String laumontPeptideList = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/LaumontResults/tmp.tsv";
		BufferedReader BR = new BufferedReader(new FileReader(laumontPeptideList));
		BufferedWriter BW = new BufferedWriter(new FileWriter("Laumont.ided.MAPs.fasta"));
		String line = null;
		
		BR.readLine(); // skip header
		
		int index = 1;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			BW.append(">"+fields[1]+"-"+index);
			BW.newLine();
			BW.append(fields[0]);
			BW.newLine();
			index++;
		}
		BW.close();
		BR.close();
	}
}
