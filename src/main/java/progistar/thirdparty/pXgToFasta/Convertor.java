package progistar.thirdparty.pXgToFasta;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class Convertor {

	public static void main(String[] args) throws IOException {
		String pXgResult = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\convenPeptides.LaumontOnly.txt";
		BufferedReader BR = new BufferedReader(new FileReader(pXgResult));
		String line = null;
		
		BR.readLine(); // skip header
		
		int index = 1;
		int length = 9;
		String selectHLA = "HLA-B08:01";
		
		Hashtable<String, String> duplicatations = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String pept = fields[0];
			String hla = fields[7];
			if(pept.length() != length) continue;
			if(!hla.equalsIgnoreCase(selectHLA)) continue;
			if(duplicatations.get(pept) == null) {
//				System.out.println(">"+index);
				System.out.println(pept);
				index++;
				
				duplicatations.put(pept, "");
			}
		}
		
		BR.close();
	}
}
