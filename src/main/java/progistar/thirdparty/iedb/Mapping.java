package progistar.thirdparty.iedb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Mapping {

	public static void main(String[] args) throws IOException {
		File iedbFile = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\IEDB\\IEDB.host");
		File pXgFile = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\IEDB\\noncanonical.txt");
		
		Hashtable<String, String> iedbMapper = loadIEDB(iedbFile);
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgFile.getAbsolutePath()+".iedb"));
		String line = null;
		
		// skip header
		BW.append(BR.readLine());
		BW.newLine();
		BW.append(BR.readLine());
		BW.append("\tIEDB");
		BW.newLine();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[20];
			String iedbInfo = iedbMapper.get(peptide);
			if(iedbInfo == null) {
				iedbInfo = "NA";
			}
			BW.append(line+"\t"+iedbInfo);
			BW.newLine();
		}
		BW.close();
		BR.close();
	}
	
	public static Hashtable<String, String> loadIEDB (File file) throws IOException {
		Hashtable<String, String> iedb = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[0].split("\\s")[0];
			iedb.put(peptide, line);
		}
		
		BR.close();
		
		return iedb;
	}
}
