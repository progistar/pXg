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
		File iedbFile = new File("/Users/gistar/projects/pXg/IEDB/EpitopeTable_20220925_MHC.csv");
		File pXgFile = new File("/Users/gistar/projects/pXg/IEDB/198UniquencMAPs.tsv");
		
		Hashtable<String, String> iedbMapper = loadIEDB(iedbFile);
		
		BufferedReader BR = new BufferedReader(new FileReader(pXgFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pXgFile.getAbsolutePath()+".iedb"));
		String line = null;
		
		// skip header
		BW.append(BR.readLine());
		BW.append("\tAtigenName\tAntigenID\tParentProtein\tParentID\tIEDBComment");
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String peptide = fields[20];
			
			String iedbInfo = iedbMapper.get(peptide);
			if(iedbInfo == null) {
				iedbInfo = "NotMatched";
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
		BR.readLine();	
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\\,");
			fields[13] = fields[13].replace("\"", "");
			fields[2] = fields[2].replace("\"", "");
			String peptide = fields[2];
			String antigenName = fields[9].replace("\"", "");
			String antigenID = fields[10].replace("\"", "");
			String parentProtein = fields[11].replace("\"", "");
			String parentID = fields[12].replace("\"", "");
			String comment = fields[16].replace("\"", "");
			
			if(antigenName.length() == 0) {
				antigenName = "NA";
			}
			if(antigenID.length() == 0) {
				antigenID = "NA";
			}
			if(parentProtein.length() == 0) {
				parentProtein = "NA";
			}
			if(parentID.length() == 0) {
				parentID = "NA";
			}
			if(comment.length() == 0) {
				comment = "NA";
			}
			
			iedb.put(peptide, antigenName+"\t"+antigenID+"\t"+parentProtein+"\t"+parentID+"\t"+comment);
		}
		
		BR.close();
		
		return iedb;
	}
}
