package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class Laumont {

	public static void main(String[] args) throws IOException {
		String pXgFile = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\subjectM.5ppm.002.rep1.netMHCpan.pXg";
		String laumontNovelFile = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\crypticPeptides.rmUniProt.txt";
		String laumontPCFile = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\convenPeptides.txt";
		BufferedReader BR = new BufferedReader(new FileReader(laumontNovelFile));
		String line = null;
		
		Hashtable<String, String> lMapper = new Hashtable<String, String>();
		
		BR.readLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			lMapper.put(peptide, "Novel");
		}
		
		BR.close();
		
		BR = new BufferedReader(new FileReader(laumontPCFile));
		BR.readLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0];
			String str = lMapper.get(peptide);
			if(str == null) lMapper.put(peptide, "PC");
			else lMapper.put(peptide, "Novel;PC");
		}
		BR.close();
		
		BR = new BufferedReader(new FileReader(pXgFile));
		System.out.println(BR.readLine()+"\tLaumont");
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[19];
			if(lMapper.get(peptide) == null) {
				System.out.println(line+"\tNO");
			} else {
				System.out.println(line+"\t"+lMapper.get(peptide));
			}
		}
		BR.close();
	}
}
