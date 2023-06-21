package progistar.pXg.thirdparty.iedoimmuno;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class OutputParser {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/tools/iedb_immunogenecity3.0/pred_noncanonical.txt");
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".txt", ".format.txt")));
		String line = null;
		
		String hla = "";
		int stepIdx = 0;
		
		BW.append("HLA\tPeptide\tScore");
		BW.newLine();
		System.out.println("Discarded list below:");
		while((line = BR.readLine()) != null) {
			if(line.startsWith("Allele")) {
				System.out.println(line.split("\\s")[1]);
				stepIdx = 0;
			}
			if(line.startsWith("allele")) {
				hla = line.split("\\s")[1];
				stepIdx = 1;
			} else if(line.startsWith("peptide") && stepIdx == 1) {
				stepIdx ++;
			} else if(stepIdx == 2) {
				String[] fields = line.split("\\,");
				String peptide = fields[0];
				String score = fields[2];
				BW.append(hla+"\t"+peptide+"\t"+score);
				BW.newLine();
			}
		}
		BR.close();
		BW.close();
	}
}
