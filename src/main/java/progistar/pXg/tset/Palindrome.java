package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;

public class Palindrome {

	public static void main(String[] args) throws IOException {
		String fileName = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Cuevas_CellRep2021\\PEAKS\\HBL1.10ppm.0025.tsv";
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		String line = null;
		
		Hashtable<String, Boolean> isAlready = new Hashtable<String, Boolean>();
		Hashtable<String, String> peptides = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String key = fields[1]+"_"+fields[4];
			
			if(isAlready.get(key) == null) {
				isAlready.put(key, true);
			} else continue;
			
			String pept = fields[3].replaceAll("[+-0123456789.\\(\\)]", "");
			pept = pept.replace("I", "L");
			
			if(pept.length() < 8 || pept.length() > 15) continue;
			peptides.put(pept, pept);
		}
		
		BR.close();
		
		Iterator<String> pepts = (Iterator<String>) peptides.keys();
		int palCount = 0;
		while(pepts.hasNext()) {
			String pept = pepts.next();
			boolean isPalindrome = true;
			System.out.println(pept);
			for(int i=0; i<pept.length(); i++) {
				if(pept.charAt(i) != pept.charAt(pept.length()-i-1)) {
					isPalindrome = false;
					break;
				}
			}
			
			if(isPalindrome) palCount++;
		}
		
		System.out.println(palCount +"/"+peptides.size());
	}
}
